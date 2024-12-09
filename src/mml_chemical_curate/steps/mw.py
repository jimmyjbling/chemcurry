"""molecular weight based curation functions"""

from rdkit.Chem.rdMolDescriptors import CalcExactMolWt

from ..flags import CurationIssue
from .base import SingleCurationStep


class CurateMW(SingleCurationStep):
    """remove compounds with molecular weight above or below some cutoff"""

    def __init__(self, min_mw: float = 1, max_mw: float = float("inf")):
        """
        Initialize a curation step

        Notes
        -----
        Bounds are inclusive

        Parameters
        ----------
        min_mw: float, default=1
            the minimum molecular weight to be considered
        max_mw: float, default=inf
            the maximum molecular weight to be considered
        """
        super().__init__()
        self.issue = CurationIssue.wrong_mw
        self.rank = 4
        self.min_mw = min_mw
        self.max_mw = max_mw

        if self.min_mw <= 0:
            raise ValueError(f"min_mw must be greater than 0; got {self.min_mw}")
        if self.min_mw > self.max_mw:
            raise ValueError(
                f"min_mw cannot be larger than man_mw; "
                f"`min_mw`: {self.min_mw} `max_mw`: {self.max_mw}"
            )

    def _func(self, molecules):
        for mol in molecules:
            if mol.failed_curation:
                continue
            if not (self.min_mw <= CalcExactMolWt(mol.mol) <= self.max_mw):
                mol.flag_issue(self.issue)
