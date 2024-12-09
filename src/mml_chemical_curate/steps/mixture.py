"""mixture based curation functions"""

from rdkit.Chem import GetMolFrags
from rdkit.Chem.MolStandardize.rdMolStandardize import LargestFragmentChooser

from .base import SingleCurationStep


class CurateMixtures(SingleCurationStep):
    """Flag compounds that have a mixture"""

    def __init__(self):
        super().__init__()
        self.issue = "chemical contains a mixture"
        self.dependency = {"CurateRemoveH|CurateRemoveAllH"}
        self.rank = 2

    def _func(self, molecules):
        for mol in molecules:
            if mol.failed_curation:
                continue
            if len(GetMolFrags(mol)) > 1:
                mol.flag_issue(self.issue)


class CurateDemix(SingleCurationStep):
    """Extracts the largest component of a mixture as the new molecule"""

    def __init__(self):
        super().__init__()
        self.issue = "failed to de-mix the chemical"
        self.note = "de-mixed by picking largest chemical compound"
        self.rank = 2

    def _func(self, molecules):
        _chooser = LargestFragmentChooser()
        for mol in molecules:
            if mol.failed_curation:
                continue
            mol.update_mol(_chooser.choose(mol), self.note)
