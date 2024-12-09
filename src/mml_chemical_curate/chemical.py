"""a wrapper class for RDKit chemicals"""

import importlib
from copy import deepcopy
from typing import List, Optional, Union

from rdkit.Chem import Mol

from .flags import CurationIssue, CurationNote


class SmilesHashingMixin:
    """enables class to generate a Smiles hash"""

    _hash_function = importlib.import_module("rdkit.Chem.rdMolHash").HashFunction.CanonicalSmiles
    _hasher = importlib.import_module("rdkit.Chem.rdMolHash").MolHash
    _smiles_hash_cache: str = ""

    def _get_smiles_hash(self, mol: Mol) -> str:
        """Gets SMILES hash of a Mol object"""
        return self._hasher(mol, self._hash_function)

    def _cache_smiles_hash(self, mol: Mol):
        """Gets SMILES hash of a Mol object"""
        self._smiles_hash_cache = self._get_smiles_hash(mol)


class MolHashingMixin:
    """enables class to generate a mol hash"""

    _mol_hash_cache: int = hash(None)

    @staticmethod
    def _get_mol_hash(mol: Mol) -> int:
        """Gets mol hash of a Mol object"""
        return hash(mol.ToBinary())

    def _cache_mol_hash(self, mol: Mol):
        """Gets mol hash of a Mol object"""
        self._mol_hash_cache = self._get_mol_hash(mol)


class _BaseChemical:
    """abstract base class for Chemicals"""

    def __init__(self, idx: int):
        """
        Initialization

        Parameters
        ----------
        idx: int
            the index in whatever file or iterable object the chemical came from
        """
        self.idx = idx


class InvalidChemical(_BaseChemical):
    """A dummy chemical for when rdkit cannot load in a specific datapoint"""

    def __init__(self, chemical_rep: str, idx: int):
        """
        Initialize a InvalidChemical object

        Parameters
        ----------
        chemical_rep: str
            the str of the chemical representation RDKit failed to load
        """
        super().__init__(idx)
        self.chemical_rep = chemical_rep


class Chemical(_BaseChemical, MolHashingMixin, SmilesHashingMixin):
    """
    Wrapper class for RDKit Mol objects

    RDKit, for reason I don't understand, has some functions that will return
    new mols rather than edit them in place, and other that will just edit them
    inplace and not create a new Mol object. This behavior is not predictable.

    Further, for our curator we want all of our updates to be inplace (which will
    help with tracking history). This class simple enables that functionality
    """

    def __init__(
        self,
        mol: Mol,
        idx: int,
        label: Optional[Union[str, int, float]] = None,
        track_history: bool = False,
    ):
        """
        Initialize a Chemical object

        Parameters
        ----------
        mol: Rdkit.Chem.Mol
            the rdkit mol for the object
        idx: int
            the index in whatever file or iterable object the chemical came from
        label: Optional[Union[str, int, float]]
            the label associated with the chemical (if there is one)
        track_history:
            track the history of molecule and label updates
        """
        super().__init__(idx)
        self.mol = deepcopy(mol)
        self.label = deepcopy(label)

        self.issue: Optional[CurationIssue] = None
        self.notes: List[CurationNote] = []

        self.failed_curation: bool = False

        self._track_history = track_history
        self.mol_history: List[Mol] = []
        self.label_history: List[Optional[Union[str, int, float]]] = []

        self._cache_mol_hash(self.mol)
        self._cache_smiles_hash(self.mol)

    @property
    def track_history(self):
        """track_history parameter"""
        return self._track_history

    @track_history.setter
    def track_history(self, value):
        """Prevent track history from being changed after initialization of obj"""
        raise RuntimeError("'track_history' cannot be change after object initialization")

    def update_mol(self, new_mol: Optional[Mol], note: CurationNote):
        """
        Update the mol to a new mol and take the associate update note

        Will only make an update and attached a note if it is detected that
        the mol has actually changed. If not, no update will occur
        Will also track the history of the mol if track_history is True

        Parameters
        ----------
        new_mol: rdkit.Chem.Mol
            the mol object to be update the current mol to
        note: CurationNote
            the note associated with this update
        """
        _hash = self._get_mol_hash(new_mol)
        if _hash != self._mol_hash_cache:
            self.notes.append(note)
            if self.track_history:
                self.mol_history.append(deepcopy(self.mol))
            self._mol_hash_cache = _hash
            self.mol = new_mol

    def update_label(
        self, new_label: Optional[Union[str, int, float]], note: CurationNote, force: bool = False
    ):
        """
        Update the label to a new mol label take the associate update note

        Will only make an update and attached a note if it is detected that
        the label has actually changed (either value OR type). If not, no update will occur
        Will also track the history of the label if track_history is True

        Parameters
        ----------
        new_label: str, int, float or None
            the label value to be update the current label to
        note: CurationNote
            the note associated with this update
        force: bool, default=False
            force the attachment of the note and history tracking
            even if no change in label detected
        """
        if (new_label is not self.label) or (new_label != self.label) or force:
            self.notes.append(note)
            if self.track_history:
                self.label_history.append(deepcopy(self.label))
            self.label = new_label

    def flag_issue(self, issue: CurationIssue):
        """
        Flag this chemical with an issue

        Notes
        -----
        Every chemical should have only 1 issue (the first one attached to it)
        Once a chemical is flagged with an issue, the 'failed_curation' flag is set
        to True and no more curation will occur

        Parameters
        ----------
        issue: CurationIssue
            the issue to flag for the chemical
        """
        if not self.failed_curation:
            self.issue = issue
            self.failed_curation = True


class ChemicalGroup(MolHashingMixin):
    """A group of chemical objects"""

    def __init__(self, chemicals: List[Chemical]):
        """
        Initialize a ChemicalGroup object

        Parameters
        ----------
        chemicals: List[Chemical]
            a list of chemical objects to be in the group
        """
        self.chemicals = chemicals

    @property
    def labels(self) -> List[Optional[Union[str, int, float]]]:
        """The labels of the chemicals in this group as a list"""
        return [chemical.label for chemical in self.chemicals]

    def update_mol(self, new_mol: Optional[Mol], note: CurationNote):
        """
        Update the mol to a new mol and take the associate update note

        Will only make an update and attached a note if it is detected that
        the mol has actually changed. If not, no update will occur
        Will also track the history of the mol if track_history is True

        Parameters
        ----------
        new_mol: rdkit.Chem.Mol
            the mol object to be update the current mol to
        note: CurationNote
            the note associated with this update
        """
        for chemical in self.chemicals:
            chemical.update_mol(new_mol, note)

    def update_label(
        self, new_label: Optional[Union[str, int, float]], note: CurationNote, force: bool = False
    ):
        """
        Update the label to a new mol label take the associate update note

        Will only make an update and attached a note if it is detected that
        the label has actually changed (either value OR type). If not, no update will occur
        Will also track the history of the label if track_history is True

        Parameters
        ----------
        new_label: str, int, float or None
            the label value to be update the current label to
        note: CurationNote
            the note associated with this update
        force: bool, default=False
            force the attachment of the note and history tracking
            even if no change in label detected
        """
        for chemical in self.chemicals:
            chemical.update_label(new_label, note, force)

    def flag_issue(self, issue: CurationIssue):
        """
        Flag this chemical with an issue

        Notes
        -----
        Every chemical should have only 1 issue (the first one attached to it)
        Once a chemical is flagged with an issue, the 'failed_curation' flag is set
        to True and no more curation will occur

        Parameters
        ----------
        issue: CurationIssue
            the issue to flag for the chemical
        """
        for chemical in self.chemicals:
            chemical.flag_issue(issue)


class ChemicalSmilesGroup(ChemicalGroup, SmilesHashingMixin):
    """A group of chemical objects that all share the same SMILES hash"""

    def __init__(self, chemicals: List[Chemical]):
        """
        Initialize a ChemicalSmilesGroup object

        Parameters
        ----------
        chemicals: List[Chemical]
            a list of chemical objects to be in the group
        """
        super().__init__(chemicals)

        self._cache_smiles_hash(self.chemicals[0].mol)

        # check that group is valid
        for chemical in self.chemicals:
            chemical._cache_smiles_hash(chemical.mol)
            if chemical._smiles_hash_cache != self._smiles_hash_cache:
                raise ValueError(
                    f"ChemicalSmilesGroup must have all chemical share the same smiles hash;"
                    f"found '{chemical._smiles_hash_cache}' and '{self._smiles_hash_cache}'"
                )

    def get_group_smiles(self) -> str:
        """Get the SMILES shared by the whole group"""
        return self._smiles_hash_cache
