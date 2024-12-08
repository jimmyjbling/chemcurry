"""inorganic curation functions"""

from rdkit.Chem import MolFromSmarts

from ..flags import CurationIssue
from .base import CurationStep


NON_ORGANIC = MolFromSmarts("[!#6;!#5;!#8;!#7;!#16;!#15;!F;!Cl;!Br;!I;!Na;!K;!Mg;!Ca;!Li;!#1]")


class CurateInorganic(CurationStep):
    """Flags compounds that have inorganic atoms"""

    def __init__(self):
        self.issue = CurationIssue.inorganic
        self.rank = 4

    def _func(self, molecules):
        for mol in molecules:
            if mol.failed_curation:
                continue
            if mol.mol.HasSubstructMatch(NON_ORGANIC):
                mol.flag_issue(self.issue)
