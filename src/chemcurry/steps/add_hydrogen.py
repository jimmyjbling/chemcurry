"""adding hydrogen curation functions"""

from typing import Optional

from rdkit.Chem import AddHs, Mol

from .base import Update, check_for_boost_rdkit_error


class AddH(Update):
    """
    Curation step to add explict hydrogen atoms to molecules

    Notes
    -----
    Will add all possible hydrogen atoms based on valance rules.
    If adding hydrogen fails due to an internal RDKit error,
    the molecule is flagged with an issue.

    Attributes
    ----------
    issue : str
        Description of issue related to the curation step
    note : str
        Description of what changes were made when a molecule was updated
    dependency : set
        Set containing the names of preceding required steps.
        If no dependency, will be an empty set.
    """

    def __init__(self):
        """Initialize the curation step"""
        self.issue = "failed to add explicit hydrogen atoms"
        self.note = "added explicit hydrogen atoms"

    def _update(self, mol: Mol) -> Optional[Mol]:
        """Attempts to add explicit hydrogen atoms to molecule; returns None if fails"""
        try:
            return AddHs(mol)
        except TypeError as e:
            if check_for_boost_rdkit_error(str(e)):
                return None
            else:
                raise e
