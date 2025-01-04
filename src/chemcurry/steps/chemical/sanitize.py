"""sanitization curation steps"""

from typing import Optional

from rdkit.Chem import Mol
from rdkit.Chem.rdmolops import SANITIZE_NONE, SanitizeMol

from ..base import Update


class SanitizeMolecule(Update):
    """
    Uses rdkit to sanitize molecules

    Notes
    -----
    If anything other than a SANITIZE_NONE flag is returned,
    will flag with an issue

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
        super().__init__()
        self.issue = "failed to sanitize chemical"
        self.note = "chemical sanitized"
        self.rank = 3

    def _update(self, mol: Mol) -> Optional[Mol]:
        _flags = SanitizeMol(mol)
        if _flags != SANITIZE_NONE:
            return None
        else:
            return mol
