"""a wrapper class for RDKit chemicals"""

from copy import deepcopy
from typing import List, Optional, Union

from rdkit.Chem import Mol

from .flags import CurationIssue, CurationNote


class Chemical:
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
        mol: Optional[Mol],
        label: Optional[Union[str, int, float]] = None,
        track_history: bool = False,
    ):
        self.mol = deepcopy(mol)
        self.label = deepcopy(label)

        self.issue: Optional[CurationIssue] = None
        self.notes: List[CurationNote] = []

        self.failed_curation: bool = False

        self._track_history = track_history
        self.mol_history: List[Mol] = []
        self._mol_hash = self._get_mol_hash(mol)
        self.label_history: List[Optional[Union[str, int, float]]] = []
        self._label_hash = hash(self.label)

    @property
    def track_history(self):
        """track_history parameter"""
        return self._track_history

    @track_history.setter
    def track_history(self, value):
        """Prevent track history from being changed after initialization of obj"""
        raise RuntimeError("'track_history' cannot be change after object initialization")

    @staticmethod
    def _get_mol_hash(mol: Mol) -> int:
        """Convert mol into a hash for easy equality"""
        if Mol is not None:
            return hash(mol.ToBinary())
        else:
            return hash(mol)

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
        if self._get_mol_hash(new_mol) != self._mol_hash:
            self.notes.append(note)
            if self.track_history:
                self.mol_history.append(deepcopy(self.mol))
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
        if issue is None:
            self.issue = issue
            self.failed_curation = True

    def add_note(self, note: CurationNote):
        """
        Add a note to the current chemical

        Warning
        -------
        This should probably not be used, as notes without update
        are kind of pointless (with how the package is built)

        Notes
        -----
        This is to enable the addition of a note not associated
        with a change to the chemical.

        Parameters
        ----------
        note: CurationNote
            the note to add
        """
        self.notes.append(note)
