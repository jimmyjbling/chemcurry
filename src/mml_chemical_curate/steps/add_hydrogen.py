"""adding hydrogen curation functions"""

from rdkit.Chem import AddHs, Mol

from .base import SingleCurationStep, check_for_boost_rdkit_error


class CurateAddH(SingleCurationStep):
    """
    Curation function that attempts to add explict hydrogen atoms to molecules
    """

    def __init__(self):
        self.issue = "failed to add explicit hydrogen atoms"
        self.note = "added explicit hydrogen atoms"
        self.rank = 3

    def _func(self, chemical):
        try:
            chemical.update_mol(AddHs(Mol), self.get_note_text())
        except TypeError as e:
            if check_for_boost_rdkit_error(str(e)):
                chemical.flag_issue(self.get_issue_text())
            else:
                raise e
