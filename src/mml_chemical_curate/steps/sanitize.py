"""sanitization curation steps"""

from rdkit.Chem.rdmolops import SANITIZE_NONE, SanitizeMol

from ..flags import CurationIssue, CurationNote
from .base import CurationStep


class CurateSanitize(CurationStep):
    """
    Uses rdkit to sanitize molecules

    Notes
    -----
    Will flag with issue if any sanitization flags are return
    """

    def __init__(self):
        super().__init__()
        self.issue = CurationIssue.sanitize_failed
        self.note = CurationNote.sanitized
        self.rank = 3

    def _func(self, molecules):
        for mol in molecules:
            if mol.failed_curation:
                continue
            _flags = SanitizeMol(mol.mol)
            if _flags != SANITIZE_NONE:
                mol.flag_issue(self.issue)
