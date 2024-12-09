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
