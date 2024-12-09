"""a wrapper class for RDKit chemicals"""

import importlib
from copy import deepcopy
from typing import List, Optional, Union

from rdkit.Chem import Mol


class SmilesHashingMixin:
    """enables class to generate a Smiles hash"""

    _hash_function = importlib.import_module("rdkit.Chem.rdMolHash").HashFunction.CanonicalSmiles
    _hasher = importlib.import_module("rdkit.Chem.rdMolHash").MolHash
    smiles_hash_cache: str = ""
    mol: Mol

    def get_smiles_hash(self, mol: Mol) -> str:
        """Gets SMILES hash of a Mol object"""
        return self._hasher(mol, self._hash_function)

    def cache_smiles_hash(self):
        """Gets SMILES hash of a Mol object"""
        self.smiles_hash_cache = self.get_smiles_hash(self.mol)


class MolHashingMixin:
    """enables class to generate a mol hash"""

    mol_hash_cache: int = hash(None)
    mol: Mol

    @staticmethod
    def get_mol_hash(mol: Mol) -> int:
        """Gets mol hash of a Mol object"""
        return hash(mol.ToBinary())

    def cache_mol_hash(self):
        """Gets mol hash of a Mol object"""
        self.mol_hash_cache = self.get_mol_hash(self.mol)


class Chemical(MolHashingMixin, SmilesHashingMixin):
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
        label: Optional[Union[str, int, float]] = None,
        track_history: bool = False,
    ):
        """
        Initialize a Chemical object

        Parameters
        ----------
        mol: Rdkit.Chem.Mol
            the rdkit mol for the object
        label: Optional[Union[str, int, float]]
            the label associated with the chemical (if there is one)
        track_history:
            track the history of molecule and label updates
        """
        self.mol = deepcopy(mol)
        self.label = deepcopy(label)

        self.issue: Optional[str] = None
        self.notes: List[str] = []

        self.failed_curation: bool = False

        self._track_history = track_history
        self.mol_history: List[Mol] = []
        self.label_history: List[Optional[Union[str, int, float]]] = []

        self.cache_mol_hash()
        self.cache_smiles_hash()

    @property
    def track_history(self):
        """track_history parameter"""
        return self._track_history

    @track_history.setter
    def track_history(self, value):
        """Prevent track history from being changed after initialization of obj"""
        raise RuntimeError("'track_history' cannot be change after object initialization")

    def update_mol(self, new_mol: Optional[Mol], note: str):
        """
        Update the mol to a new mol and take the associate update note

        Will only make an update and attached a note if it is detected that
        the mol has actually changed. If not, no update will occur
        Will also track the history of the mol if track_history is True

        Parameters
        ----------
        new_mol: rdkit.Chem.Mol
            the mol object to be update the current mol to
        note: str
            the note associated with this update
        """
        _hash = self.get_mol_hash(new_mol)
        if _hash != self.mol_hash_cache:
            self.notes.append(note)
            if self.track_history:
                self.mol_history.append(deepcopy(self.mol))
            self.mol_hash_cache = _hash
            self.mol = new_mol

    def update_label(
        self, new_label: Optional[Union[str, int, float]], note: str, force: bool = False
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
        note: str
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

    def flag_issue(self, issue: str):
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

    def tag_note(self, note: str):
        """
        Tag a new note to the chemical

        Primarily used for tagging notes that are not
        associated with updates to chemical

        Parameters
        ----------
        note: str
            the note to tag
        """
        self.notes.append(note)


class BaseChemicalGroup:
    """A group of chemical objects"""

    def __init__(self, group_id: Union[int, str], chemicals: List[Chemical]):
        """
        Initialize a BaseChemicalGroup object

        Parameters
        ----------
        group_id: Union[int, str]
            a unique id for the group
            uniqueness is not enforced but if not unique could
            cause unexpected and uncaught issues
        chemicals: List[Chemical]
            a list of chemical objects to be in the group
        """
        self.group_id = group_id
        self.chemicals = chemicals

    @property
    def labels(self) -> List[Optional[Union[str, int, float]]]:
        """The labels of the chemicals in this group as a list"""
        return [chemical.label for chemical in self.chemicals]

    @property
    def passing_chemicals(self):
        """Return a list of all the currently passing chemicals in the groups"""
        return [chemical for chemical in self.chemicals if not chemical.failed_curation]

    @property
    def flagged_chemicals(self):
        """Return a list of all the currently flagged chemicals in the group"""
        return [chemical for chemical in self.chemicals if chemical.failed_curation]

    def update_mol(self, new_mol: Optional[Mol], note: str):
        """
        Update the mol to a new mol and take the associate update note

        Will only make an update and attached a note if it is detected that
        the mol has actually changed. If not, no update will occur
        Will also track the history of the mol if track_history is True

        Parameters
        ----------
        new_mol: rdkit.Chem.Mol
            the mol object to be update the current mol to
        note: str
            the note associated with this update
        """
        for chemical in self.chemicals:
            chemical.update_mol(new_mol, note)

    def update_label(
        self, new_label: Optional[Union[str, int, float]], note: str, force: bool = False
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
        note: str
            the note associated with this update
        force: bool, default=False
            force the attachment of the note and history tracking
            even if no change in label detected
        """
        for chemical in self.chemicals:
            chemical.update_label(new_label, note, force)

    def flag_issue(self, issue: str):
        """
        Flag all the currently passing chemical with an issue

        Notes
        -----
        Every chemical should have only 1 issue (the first one attached to it)
        Once a chemical is flagged with an issue, the 'failed_curation' flag is set
        to True and no more curation will occur

        Parameters
        ----------
        issue: str
            the issue to flag for the chemical
        """
        for chemical in self.passing_chemicals:
            chemical.flag_issue(issue)

    def tag_note(self, note: str):
        """
        Tag a note to all the currently passing chemicals in the group chemical

        Primarily used for tagging notes that are not
        associated with updates to chemical

        Parameters
        ----------
        note: str
            the note to tag
        """
        for chemical in self.passing_chemicals:
            chemical.tag_note(note)


class ChemicalMoleculeGroup(BaseChemicalGroup, MolHashingMixin):
    """A group of chemical objects that share the same mol hash"""

    def __init__(self, group_id: Union[str, int], chemicals: List[Chemical]):
        """
        Initialize a ChemicalMoleculeGroup object

        Parameters
        ----------
        chemicals: List[Chemical]
            a list of chemical objects to be in the group
        """
        super().__init__(group_id, chemicals)
        self.mol = chemicals[0].mol
        self.cache_mol_hash()

        # check that group is valid
        for chemical in self.chemicals:
            chemical.get_mol_hash(chemical.mol)
            if chemical.mol_hash_cache != self.mol_hash_cache:
                raise ValueError(
                    f"ChemicalSmilesGroup must have all chemical share the same smiles hash;"
                    f"found '{chemical.mol_hash_cache}' and '{self.mol_hash_cache}'"
                )

    def get_group_mol_hash(self) -> int:
        """Get the SMILES shared by the whole group"""
        return self.mol_hash_cache


class ChemicalSmilesGroup(BaseChemicalGroup, SmilesHashingMixin):
    """A group of chemical objects that all share the same SMILES hash"""

    def __init__(self, group_id: Union[str, int], chemicals: List[Chemical]):
        """
        Initialize a ChemicalSmilesGroup object

        Parameters
        ----------
        chemicals: List[Chemical]
            a list of chemical objects to be in the group
        """
        super().__init__(group_id, chemicals)
        self.mol = chemicals[0].mol
        self.cache_smiles_hash()

        # check that group is valid
        for chemical in self.chemicals:
            chemical.get_smiles_hash(chemical.mol)
            if chemical.smiles_hash_cache != self.smiles_hash_cache:
                raise ValueError(
                    f"ChemicalSmilesGroup must have all chemical share the same smiles hash;"
                    f"found '{chemical.smiles_hash_cache}' and '{self.smiles_hash_cache}'"
                )

    def get_group_smiles(self) -> str:
        """Get the SMILES shared by the whole group"""
        return self.smiles_hash_cache
