"""molecular weight based curation functions"""

from rdkit.Chem import Mol
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt

from .base import Filter


class FilterMW(Filter):
    """
    Flag compounds with molecular weight above or below some cutoff

    Attributes
    ----------
    issue : str
        Description of issue related to the curation step
    dependency : set
        Set containing the names of preceding required steps.
        If no dependency, will be an empty set.
    """

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
        self.min_mw = min_mw
        self.max_mw = max_mw

        if self.min_mw <= 0:
            raise ValueError(f"min_mw must be greater than 0; got {self.min_mw}")
        if self.min_mw > self.max_mw:
            raise ValueError(
                f"min_mw cannot be larger than man_mw; "
                f"`min_mw`: {self.min_mw} `max_mw`: {self.max_mw}"
            )

        self.issue = (
            f"chemical had a molecular weight below " f"{self.min_mw} or above {self.max_mw}"
        )

    def _filter(self, mol: Mol) -> bool:
        """Returns True if molecular weight is within bounds"""
        return self.min_mw <= CalcExactMolWt(mol) <= self.max_mw
