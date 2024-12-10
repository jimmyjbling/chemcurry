"""mixture based curation functions"""

import importlib

from rdkit.Chem import GetMolFrags

from .base import SingleCurationStep


class CurateMixtures(SingleCurationStep):
    """Flag compounds that have a mixture"""

    def __init__(self):
        super().__init__()
        self.issue = "chemical contains a mixture"
        self.dependency = {"CurateRemoveH|CurateRemoveAllH"}
        self.rank = 2

    def _func(self, chemical):
        if len(GetMolFrags(chemical.mol)) > 1:
            chemical.flag_issue(self.get_issue_text())


class CurateDemix(SingleCurationStep):
    """Extracts the largest component of a mixture as the new molecule"""

    def __init__(self):
        super().__init__()
        self.issue = "failed to de-mix the chemical"
        self.note = "de-mixed by picking largest chemical compound"
        self.rank = 2

        self._chooser = importlib.import_module(
            "rdkit.Chem.MolStandardize.rdMolStandardize"
        ).LargestFragmentChooser()

    def _func(self, chemical):
        chemical.update_mol(self._chooser.choose(chemical.mol), self.get_note_text())
