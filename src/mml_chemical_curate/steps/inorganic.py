"""inorganic curation functions"""

from rdkit.Chem import Mol, MolFromSmarts

from .base import Filter


NON_ORGANIC = MolFromSmarts("[!#6;!#5;!#8;!#7;!#16;!#15;!F;!Cl;!Br;!I;!Na;!K;!Mg;!Ca;!Li;!#1]")


class FlagInorganic(Filter):
    """
    Curation step to filter out molecule containing inorganic atoms

    Notes
    -----
    Inorganic atoms are defined as any atom not in the following list:
    B, C, N, O, H, P, S, F, Cl, Br, I, Na, K, Mg, Ca, Li

    Attributes
    ----------
    issue : str
        Description of issue related to the curation step
    dependency : set
        Set containing the names of preceding required steps.
        If no dependency, will be an empty set.
    """

    def __init__(self):
        self.issue = "chemical contained inorganic atoms"

    def _filter(self, mol: Mol) -> bool:
        """Returns True if molecule contains inorganic atom"""
        return not mol.HasSubstructMatch(NON_ORGANIC)
