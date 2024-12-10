"""inorganic curation functions"""

from rdkit.Chem import MolFromSmarts

from .base import SingleCurationStep


NON_ORGANIC = MolFromSmarts("[!#6;!#5;!#8;!#7;!#16;!#15;!F;!Cl;!Br;!I;!Na;!K;!Mg;!Ca;!Li;!#1]")


class CurateInorganic(SingleCurationStep):
    """Flags compounds that have inorganic atoms"""

    def __init__(self):
        self.issue = "chemical contained inorganic atoms"
        self.rank = 4

    def _func(self, chemical):
        if chemical.mol.HasSubstructMatch(NON_ORGANIC):
            chemical.flag_issue(self.get_issue_text())
