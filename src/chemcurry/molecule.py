"""a wrapper class for RDKit molecules"""

import abc
import importlib
from copy import deepcopy
from typing import List, Optional, Self, Union

from rdkit.Chem import Mol, MolFromSmiles


class SmilesMixin:
    """
    Enables class with mol attribute to generate or instantiated from a canonical SMILES

    This mixin is pickle-able even though it uses RDKit class
    """

    @abc.abstractmethod
    def __init__(self, mol: Mol, **kwargs):
        raise NotImplementedError

    # this is the smiles hashing function settings
    _smiles_hash_function = importlib.import_module(
        "rdkit.Chem.rdMolHash"
    ).HashFunction.CanonicalSmiles
    # this is the hash function that operates on the mol object
    _smiles_hasher = importlib.import_module("rdkit.Chem.rdMolHash").MolHash
    mol: Mol

    @classmethod
    def from_smiles(cls, id_: Union[int, str], smiles: str, **kwargs) -> Self:
        """
        Create an instance of the class from a SMILES string.

        This method parses the given SMILES string to create an RDKit Mol object
        and initializes the class with it.

        Notes
        -----
        Any object that implements this mixin should be able to be
        able to handle cases where a smiles gets rendered as `None`
        by RDKit

        Parameters
        ----------
        id_: Union[int, str]
            identifier for the molecule
        smiles: str
            The SMILES string representing the molecule.
        kwargs: dict
            Additional keyword arguments to pass during initialization.

        Returns
        -------
        Self
            An instance of the class if the SMILES
        """
        _mol = MolFromSmiles(smiles)
        return cls(id_=id_, mol=_mol, **kwargs)

    def get_smiles(self) -> str:
        """
        Generate the canonical SMILES string

        Uses the Mol object in the `mol` attribute to generate the canonical SMILES.

        Returns
        -------
        str
            The canonical SMILES string.
        """
        return self._smiles_hasher(self.mol, self._smiles_hash_function)

    def has_same_smiles(self, other_mol: Mol) -> bool:
        """
        Check if the a different Mol object has the same SMILES as this one

        compares against the Mol object in the `mol` attribute

        Parameters
        ----------
        other_mol: Mol
            the other mol object to compare to

        Returns
        -------
        bool
        """
        return self.get_smiles() == self._smiles_hasher(other_mol, self._smiles_hash_function)


class Molecule(SmilesMixin):
    """
    Wrapper class for RDKit Mol objects

    Notes
    -----
    The main justification for this wrapper class is to enable easy tracking
    the molecule history. This both enables an explict list showing the evolution of the
    molecule over time, but also help with only attaching note to molecule that actually
    were altered as the result of an Update curation step
    """

    def __init__(
        self,
        id_: Union[int, str],
        mol: Optional[Mol],
        track_history: bool = False,
    ):
        """
        Initialize a Molecule object

        Notes
        -----
        `id_` uniqueness is *NOT* enforced
        molecules can share the same id
        this is used to help with tracking compounds
        on the user side

        If mol is None, will flag with issue
        'rdkit failed to render Mol object'

        Parameters
        ----------
        id_: Union[int, str]
            identifier for the molecule
        mol: Optional[Rdkit.Chem.Mol]
            the rdkit mol for the object
            if a None, will create a dummy mol and
            automatically flag this molecule with an issue
        track_history:
            track the history of molecule and label updates
        """
        self.id_: Union[int, str] = id_

        self.issue: str = ""
        self.notes: List[str] = []

        self._track_history = track_history
        self.mol_history: List[Mol] = []

        self.failed_curation: bool = False

        self.mol: Mol
        if mol is None:
            self.mol = MolFromSmiles("")
            self.failed_curation = True
            self.issue = "rdkit failed to render Mol object"
        else:
            self.mol = mol

        self._current_hash = self._generate_mol_hash(self.mol)

    @staticmethod
    def _generate_mol_hash(mol: Mol) -> int:
        return hash(mol.ToBinary())

    @property
    def track_history(self) -> bool:
        """`track_history` parameter; True means history track in enabled"""
        return self._track_history

    @track_history.setter
    def track_history(self, value):
        """Prevent track history from being changed after initialization of obj"""
        raise RuntimeError("'track_history' cannot be change after object initialization")

    def update_mol(self, new_mol: Mol, note: str) -> bool:
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

        Returns
        -------
        bool
            True if molecule updated, False if no update occurred
        """
        if new_mol is None or len(new_mol.GetAtoms()) == 0:
            if new_mol is None:
                raise ValueError(
                    "if molecule becomes invalid, should be caught "
                    "and flagged with issue by curation step; 'None'"
                )
            if len(new_mol.GetAtoms()) == 0:
                raise ValueError(
                    "if molecule becomes invalid, should be caught "
                    "and flagged with issue by curation step; 'Empty Mol'"
                )

        _hash = self._generate_mol_hash(new_mol)
        if _hash != self._current_hash:
            self.notes.append(note)
            if self.track_history:
                self.mol_history.append(deepcopy(self.mol))
            self._current_hash = _hash
            self.mol = new_mol
            return True
        else:
            return False

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
