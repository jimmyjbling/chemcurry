"""curation workflows"""

import datetime
import os
import pickle
import warnings
from collections.abc import Iterable
from typing import Any, Dict, List, Literal, Optional, Union, overload

import numpy as np
import numpy.typing as npt
import pandas as pd
from rdkit.Chem import Mol

from .molecule import Molecule
from .steps import Filter, Update


class CurationWorkflowError(Exception):
    """
    Exception to raise when a curation workflow has an error

    Primary use case is for when the curation workflow has unmet dependencies
    """

    pass


class CurationWorkflow:
    """
    A curation workflow to curate molecules

    Assembled as a list of steps passed in the order they should be run in.
    In theory, curation steps have specified dependencies, but this is not enforced
    by the workflow (yet).
    The only check the workflow will make is to see if all `Filter` steps come
    after all `Update` steps. If not, a warning will be raised as this is not
    best practice in most cases. If you are sure you want to mix run updates after filters,
    that warning can be suppressed by setting `suppress_warnings=True` when instantiating
    the workflow.
    """

    def __init__(
        self,
        steps: List[Union[Filter, Update]],
        track_history: bool = False,
        suppress_warnings: bool = False,
    ):
        """
        Initialize a curation workflow

        Parameters
        ----------
        steps: list[CurationStep]
            a list of curation steps that should be taken
        track_history: bool, default=False
            enable history tracking of the molecules
        suppress_warnings: bool, default=False
            if True, will suppress any warnings about the workflow
        """
        self.steps: List[Union[Filter, Update]] = steps
        self.track_history = track_history

        _seen_filter_step: bool = False
        for step in self.steps:
            if isinstance(step, Filter):
                _seen_filter_step = True
            if _seen_filter_step and isinstance(step, Update):
                if not suppress_warnings:
                    warnings.warn(
                        f"update step '{step.__class__.__name__}' "
                        f"comes after a filter step; updating molecules after "
                        f"filtering could cause end result to violate the"
                        f"filter;",
                        stacklevel=2,
                    )

    def _run_workflow(self, mols: List[Molecule], from_: str = "Unknown") -> "CuratedMoleculeSet":
        """
        main point of entry into the workflow

        will run mols through all the steps sequentially

        Parameters
        ----------
        mols: list[Molecule, ...]
            molecules to curate
        from_: str, default="Unknown"
            a string to describe where the molecules are coming from
            should be the file path if loaded from a file
            should be "List of Mols" if from `curate_mols`
            should be "List of SMILES" if from `curate_smiles`

        Returns
        -------
        CuratedMoleculeSet
            set of curated molecules
        """
        # how long the pipeline steps took; init with 0 for loading mols
        _timings: List[datetime.timedelta] = [datetime.timedelta(seconds=0.0)]
        # the number failed at each step; init with number of mols that failed
        _issue_counts: List[int] = [sum([mol.failed_curation for mol in mols])]
        _note_counts: List[int] = [0]

        for step in self.steps:
            _t0 = datetime.datetime.now()
            _num_notes, _num_issues = step(mols)
            _timings.append(datetime.datetime.now() - _t0)
            _issue_counts.append(_num_issues)
            _note_counts.append(_num_notes)

        return CuratedMoleculeSet(mols, self, _issue_counts, _note_counts, _timings, from_)

    def curate_smiles(
        self, smis: Iterable[str], ids: Optional[Iterable[Union[int, str]]] = None
    ) -> "CuratedMoleculeSet":
        """
        Given a list of SMILES, run the workflow on them

        Notes
        -----
        If a SMILES cannot be rendered by RDKit,
        it will loaded as an empty molecule
        (with a SMILES of "")

        Parameters
        ----------
        smis: Iterable[str]
            smiles to curate
        ids: Optional[Iterable[Union[int, str]]], default=None
            an optional list of ids to associate with the molecules
            if left as `None` will use the index of the molecule in the list

        Returns
        -------
        CuratedMoleculeSet
            the curated molecules
        """
        if ids is not None:
            mols = [
                Molecule.from_smiles(id_=id_, smiles=smi, track_history=self.track_history)
                for id_, smi in zip(ids, smis)
            ]
        else:
            mols = [
                Molecule.from_smiles(id_=i, smiles=smi, track_history=self.track_history)
                for i, smi in enumerate(smis)
            ]

        return self._run_workflow(mols, from_="List of SMILES")

    def curate_mols(
        self, mols: Iterable[Optional[Mol]], ids: Optional[Iterable[Union[int, str]]] = None
    ) -> "CuratedMoleculeSet":
        """
        Given a list of mols, run the workflow on them

        Notes
        -----
        Can handle None inplace of Mols (as rdkit will do this)
        when encountering a None will redner the molecule as empty
        and automatically flagg it as an issue

        Parameters
        ----------
        mols: Iterable[Optional[Mol]]
            smiles to curate
        ids: Optional[Iterable[Union[int, str]]], default=None
            an optional list of ids to associate with the molecules
            if left as `None` will use the index of the molecule in the list

        Returns
        -------
        CuratedMoleculeSet
            the curated molecules
        """
        _mols: List[Molecule]
        if ids is not None:
            _mols = [
                Molecule(id_=id_, mol=mol, track_history=self.track_history)
                for id_, mol in zip(ids, mols)
            ]
        else:
            _mols = [
                Molecule(id_=i, mol=mol, track_history=self.track_history)
                for i, mol in enumerate(mols)
            ]

        return self._run_workflow(_mols, from_="List of SMILES")

    def to_string(self) -> str:
        """
        Render the workflow as a string

        Links steps by the ' -> ' characters
        Will always start with 'RDKitLoading'

        Returns
        -------
        str
        """
        return "RDKitLoading -> " + " -> ".join(step.__class__.__name__ for step in self.steps)


class CuratedMoleculeSet:
    """
    Hold the resulting molecules from a curation workflow

    Note: Users should not directly initialize this class.
          it should only be initialized by the `CurationWorkflow` class

    Attributes
    ----------
    remaining: List[int]
        the number of molecules remaining after each step
        will have length of number of curation steps + 2
        element at index 0 it raw number of molecules loaded
        element at index 1 is the number of molecules that were loaded successfully
    """

    def __init__(
        self,
        molecules: List[Molecule],
        workflow: CurationWorkflow,
        num_issues: List[int],
        num_notes: List[int],
        timings: List[datetime.timedelta],
        from_: str = "Unknown",
    ):
        """
        Initialize a CuratedMoleculeSet object

        Parameters
        ----------
        molecules: List[Molecule]
            the curated molecules
            should be the same as the ones passed to the workflow
        workflow: CurationWorkflow
            the workflow used to curate the molecules
        num_issues: List[int]
            the number of molecules with issues caused by each step
            should have length of number of curation steps + 1
            value at index 0 is the number of molecules rdkit failed to load
        num_notes: List[int]
            the number of molecules with notes caused by each step
            should have length of number of curation steps
            value at index 0 is should always be 0 (as loading molecules cannot have notes)
            for Filter steps a value of 0 is expect for num_notes
        timings: List[datetime.timedelta]
            the time it took for each step to run
            should have length of number of curation steps + 1
            value at index 0 will be a time delta of 0 (loading time is not tracked)
        from_: str, default="Unknown"
            the source of the molecules
            should be identical to the value passed the curation call
        """
        self.molecules = molecules
        self.workflow = workflow
        self.num_issues = num_issues
        self.num_notes = num_notes
        self.timings = timings

        self.from_ = from_

        self.remaining = [len(molecules)]
        for num_failures in self.num_issues:
            self.remaining.append(self.remaining[-1] - num_failures)

        assert len(self.timings) == len(workflow.steps) + 1
        assert len(self.num_issues) == len(workflow.steps) + 1
        assert len(self.remaining) == len(workflow.steps) + 2

    def to_smiles(self, include_failed: bool = False) -> List[str]:
        """
        Returns curated molecules as SMILES

        Notes
        -----
        if a molecule failed to be loaded by RDKit, its smiles with be and empty string

        Parameters
        ----------
        include_failed : bool, default=False
            Whether to include molecules that failed curation in the output list.

        Returns
        -------
        List[str]
            A list of SMILES strings representing the curated molecules. If `include_failed`
            is set to False, molecules that failed curation are excluded from the list.
        """
        return [
            mol.get_smiles()
            for mol in self.molecules
            if (not include_failed and mol.failed_curation)
        ]

    def to_mols(self, include_failed: bool = False) -> List[Mol]:
        """
        Returns curated molecules as SMILES

        Notes
        -----
        if a molecule failed to be loaded by RDKit, it will be return as
        an empty Mol object, *not* as a `None` value.

        Parameters
        ----------
        include_failed : bool, default=False
            Whether to include molecules that failed curation in the output list.

        Returns
        -------
        List[Mol]
            A list of rdkit Mol objects representing the curated molecules. If `include_failed`
            is set to True, molecules that failed curation are included in the output.
        """
        return [mol.mol for mol in self.molecules if (not include_failed and mol.failed_curation)]

    def to_pandas(
        self,
        include_notes: bool = False,
        include_issues: bool = False,
        include_failed: bool = False,
    ) -> pd.DataFrame:
        """
        Convert to a pandas DataFrame

        will have the following columns:
        - id: the id of the molecule
        - smiles: the smiles of the molecule
        - mol: the rdkit mol object
        - issue: the issue caused by the step ('PASSED' if no issue)
          only if `include_issues` is True
        - notes: the notes for the molecule as a list of strings (empty list if no notes)
          only if `include_notes` is True

        if `include_failed` is True, molecules that failed curation will be included
        otherwise they will be excluded

        Warning: this doesn't not retain all the information collected during curation.
        For example, the workflow used will be lost, as will any tracked molecule history
        Recommend to only use this at the very end after all analysis is complete and a copy
        of the results (in pickle format) is saved.

        Parameters
        ----------
        include_notes: bool, default=False
            include the notes for each molecule as a list of strings
        include_issues: bool, default=False
            include the issue for each molecule as a string ('PASSED' if no issue)
        include_failed: bool, default=False
            include molecules that failed curation in the dataframe

        Returns
        -------
        pd.DataFrame
        """
        _data: Dict[str, Any] = {"id": [], "smiles": [], "mol": [], "issue": [], "notes": []}
        for mol in self.molecules:
            if include_failed or not mol.failed_curation:
                _data["id"].append(mol.id_)
                _data["smiles"].append(mol.get_smiles())
                _data["mol"].append(mol.mol)
                if include_notes:
                    _data["issue"].append(mol.issue if mol.failed_curation else "PASSED")
                if include_issues:
                    _data["notes"].append(mol.notes)

        if not include_notes:
            del _data["notes"]

        if not include_issues:
            del _data["issue"]

        return pd.DataFrame(_data)

    @overload
    def get_passing_mask(self, as_numpy: Literal[True]) -> npt.NDArray[bool]: ...

    @overload
    def get_passing_mask(self, as_numpy: Literal[False]) -> List[bool]: ...

    def get_passing_mask(self, as_numpy: bool = False) -> Union[npt.NDArray[bool], List[bool]]:
        """
        Get a boolean mask of which molecules passed curation

        The order of the mask will match the order of the molecules passed the
        the curation workflow (and the molecules in the `CuratedMoleculeSet` object).
        Mask will have True for molecules that passed curation
        and False for molecules that failed curation.

        Parameters
        ----------
        as_numpy: bool, default=False
            if True, return as a numpy array if True with dtype np.bool_

        Returns
        -------
        Union[npt.NDArray[bool], List[bool]]
            boolean mask
        """
        _mask = [not mol.failed_curation for mol in self.molecules]
        if as_numpy:
            return np.array(_mask, dtype=np.bool_)
        else:
            return _mask

    @overload
    def get_num_issues_at_step(self, idx: int) -> int: ...

    @overload
    def get_num_issues_at_step(self, idx: str) -> List[int]: ...

    def get_num_issues_at_step(self, idx: Union[str, int]) -> Union[List[int], int]:
        """
        the number of molecules with issues caused by a given step

        Can use the index of the step in the workflow or the name of the step
        If using the name, will return the issue count of all steps with
        that name in the order the steps were run

        Parameters
        ----------
        idx: Union[str, int]
            index of the step in the workflow (int or name)

        Returns
        -------
        Union[List[int], int]
            the number(s) of molecules with an issue caused by a given step
            List[int] if using name, int if using index
        """
        if isinstance(idx, int):
            return self.num_issues[idx]
        if isinstance(idx, str):
            return [
                self.num_issues[i + 1]
                for i, step in enumerate(self.workflow.steps)
                if step.__class__.__name__ == idx
            ]

    @overload
    def get_num_notes_at_step(self, idx: int) -> int: ...

    @overload
    def get_num_notes_at_step(self, idx: str) -> List[int]: ...

    def get_num_notes_at_step(self, idx: Union[str, int]) -> Union[List[int], int]:
        """
        the number of molecules with notes caused by a given step

        Can use the index of the step in the workflow or the name of the step
        If using the name, will return the note count of all steps with
        that name in the order the steps were run

        Parameters
        ----------
        idx: Union[str, int]
            index of the step in the workflow (int or name)

        Returns
        -------
        Union[List[int], int]
            the number(s) of molecules with a note caused by a given step
            List[int] if using name, int if using index
        """
        if isinstance(idx, int):
            return self.num_notes[idx]
        if isinstance(idx, str):
            return [
                self.num_notes[i + 1]
                for i, step in enumerate(self.workflow.steps)
                if step.__class__.__name__ == idx
            ]

    @overload
    def get_remaining_molecules_after_step(self, idx: str) -> List[int]: ...

    @overload
    def get_remaining_molecules_after_step(self, idx: int) -> int: ...

    def get_remaining_molecules_after_step(self, idx: Union[str, int]) -> Union[List[int], int]:
        """
        the number of molecules passing curation after a given step

        Can use the index of the step in the workflow or the name of the step
        If using the name, will return the remaining count after all steps with
        that name in the order the steps were run

        Parameters
        ----------
        idx: Union[str, int]
            index of the step in the workflow (int or name)

        Returns
        -------
        Union[List[int], int]
            the number(s) of molecules passing curation after a given step
            List[int] if using name, int if using index
        """
        if isinstance(idx, int):
            return self.remaining[idx + 1]
        if isinstance(idx, str):
            return [
                self.remaining[i + 2]
                for i, step in enumerate(self.workflow.steps)
                if step.__class__.__name__ == idx
            ]

    def get_report_string(self) -> str:
        """
        Generate a report string for the curation workflow results

        The report will contain:
        - the version of ChemCurry AND rdkit used
        - the date and time the curation was run
        - the workflow used (both the string definition and the hash)
        - information about compounds removed, updated and remaining after each step
        - final count summary for whole workflow

        The report will be human readable and is the main way to analysis and
        share the results of your ChemCurry curation workflow run.

        Returns
        -------
        str
        """
        _report_str = (
            f"ChemCurry curation report\n"
            f"Report generated on "
            f"{datetime.datetime.strftime(datetime.datetime.now(), '%H:%M:%S, %B %d, %Y')}\n\n"
            f"Using workflow {self.workflow.to_string()}\nWorkflow hash: {hash(self.workflow)}\n\n"
            f"Loaded {self.remaining[0]} chemicals from '{self.from_}' for curation\n"
            + "\n".join(
                [
                    f"Curation Step {i}: {step.__class__.__name__}; "
                    f"Flagged {num_issues} compounds with issues; "
                    f"{num_remaining} compounds remaining"
                    if isinstance(step, Filter)
                    else f"Curation Step {i}: {step.__class__.__name__}; "
                    f"Updated/altered {num_notes} compounds with issues; "
                    f"Flagged {num_issues} compounds with issues; "
                    f"{num_remaining} compounds remaining"
                    for i, (step, num_issues, num_notes, num_remaining) in enumerate(
                        zip(
                            self.workflow.steps,
                            self.num_issues[1:],
                            self.num_notes[1:],
                            self.remaining[1:-1],
                        )
                    )
                ]
            )
            + f"\n complete workflow in {sum(self.timings, datetime.timedelta())} seconds"
            f"\n REMOVED {sum(self.num_issues)}\n" + f"ALTERED {sum(self.num_notes)}\n"
            if sum(self.num_notes) > 0
            else "" + f"FINAL COMPOUND COUNT: {self.remaining[-1]}\n"
        )
        return _report_str

    def write_report(self, path: os.PathLike):
        """
        Save the curation report to a file

        Parameters
        ----------
        path: os.PathLike
            path to save the report to
        """
        open(path, "w").write(self.get_report_string())

    def save(self, path: os.PathLike):
        """
        save a pickle of this object and all data in it

        this is the best way to save the curation results such that all
        information is available to load the results later.
        However, this is the most memory intensive way to save the results

        Parameters
        ----------
        path: str
            path to save file
        """
        pickle.dump(self, open(path, "wb"))

    def save_as_txt(self, path: os.PathLike):
        r"""
        save a txt file containing all note and issues for each molecule

        Format for each line will be:
        `<molecule_id>\t<SMILES>\t<issue>\t<notes joined by \t>`
        `issue` will be either the issue text (i.e. "contained mixture") or "PASSED"
        if there were no issues for that molecule

        Parameters
        ----------
        path: str
            path to save file
        """
        with open(path, "w") as f:
            for mol in self.molecules:
                _row = (
                    f"{mol.id_}\t"
                    f"{mol.get_smiles()}\t"
                    f"{'PASSED' if mol.failed_curation else mol.issue}\t"
                    + "\t".join(mol.notes)
                    + "\n"
                )
                f.write(_row)

    def save_as_json(self, path: os.PathLike):
        """
        save a json file containing all note and issues for each molecule

        will be a list of dicts with keys:
            id: molecule id
            smiles: molecule smiles
            issue: PASSED if no issue else the issue text
            notes: list of notes (empty list of no notes)

        Parameters
        ----------
        path: str
            path to save file
        """
        import json

        _json = []
        for mol in self.molecules:
            _data: Dict[str, Any] = {
                "id": mol.id_,
                "smiles": mol.get_smiles(),
                "issue": mol.issue if mol.failed_curation else "PASSED",
                "notes": mol.notes,
            }
            _json.append(_data)

        with open(path, "w") as f:
            json.dump(_json, f, indent=4)

    def save_as_csv(self, path: os.PathLike):
        """
        save a csv file containing the curation results with all note and issues for each molecule

        will save the note as a list inside a single element in the csv
        this can make it hard to parse back out into python objects,
        if you need to save and reload into python the `save`, `save_as_json` or
        `save_as_pandas` methods are recommended

        will have the following columns::
            id: molecule id
            smiles: molecule smiles
            issue: PASSED if no issue else the issue text
            notes: list of notes (empty list of no notes)

        Parameters
        ----------
        path: str
            path to save file
        """
        df = self.to_pandas(include_notes=True, include_issues=True, include_failed=True)
        del df["mol"]
        df.to_csv(path, index=False)

    def save_as_pandas(self, path: os.PathLike):
        """
        save a pickle file containing a pandas dataframe with all note and issues for each molecule

        unlike text passed formats, this will also save the rdkit Mol objects

        will have the following columns::
            id: molecule id
            smiles: molecule smiles
            mol: the rdkit Mol object
            issue: PASSED if no issue else the issue text
            notes: list of notes (empty list of no notes)

        Parameters
        ----------
        path: str
            path to save file
        """
        df = self.to_pandas(include_notes=True, include_issues=True, include_failed=True)
        pickle.dump(df, open(path, "wb"))
