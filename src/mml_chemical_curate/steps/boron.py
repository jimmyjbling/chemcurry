"""boron curation functions"""

from rdkit.Chem import MolFromSmarts

from .base import SingleCurationStep


class CurateBoron(SingleCurationStep):
    """Flags compounds that have Boron in them"""

    def __init__(self):
        super().__init__()
        self.issue = "contained a Boron atom"
        self.rank = 4

    def _func(self, chemical):
        if chemical.mol.HasSubstructMatch(MolFromSmarts("[#5]")):
            chemical.flag_issue(self.get_issue_text())
