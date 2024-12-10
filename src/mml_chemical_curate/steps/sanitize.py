"""sanitization curation steps"""

from rdkit.Chem.rdmolops import SANITIZE_NONE, SanitizeMol

from .base import SingleCurationStep


class CurateSanitize(SingleCurationStep):
    """
    Uses rdkit to sanitize molecules

    Notes
    -----
    Will flag with issue if any sanitization flags are return
    """

    def __init__(self):
        super().__init__()
        self.issue = "failed to sanitize chemical"
        self.note = "chemical sanitized"
        self.rank = 3

    def _func(self, chemical):
        _flags = SanitizeMol(chemical.mol)
        if _flags != SANITIZE_NONE:
            chemical.flag_issue(self.get_issue_text())
