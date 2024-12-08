"""adding hydrogen curation functions"""

from rdkit.Chem import AddHs, Mol

from ..flags import CurationIssue, CurationNote
from .base import CurationStep, check_for_boost_rdkit_error


class CurateAddH(CurationStep):
    """
    Curation function that attempts to add explict hydrogen atoms to molecules
    """

    def __init__(self):
        self.issue = CurationIssue.failed_adding_Hs
        self.note = CurationNote.added_hs
        self.rank = 3

    def _func(self, molecules):
        for mol in molecules:
            if mol.failed_curation:
                continue
            try:
                mol.update_mol(AddHs(Mol), self.note)
            except TypeError as e:
                if check_for_boost_rdkit_error(str(e)):
                    mol.flag_issue(self.issue)
                else:
                    raise e
