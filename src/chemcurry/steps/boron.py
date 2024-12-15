"""boron curation functions"""

from rdkit.Chem import Mol, MolFromSmarts

from .base import Filter


class FlagBoron(Filter):
    """
    Curation step to filter out compounds containing boron atoms.

    Notes
    -----
    Searches for boron atoms based on element number.
    This will detect any isotope of Boron.

    Attributes
    ----------
    issue : str
        Description of issue related to the curation step
    dependency : set
        Set containing the names of preceding required steps.
        If no dependency, will be an empty set.
    """

    def __init__(self):
        """Initialize the curation step"""
        super().__init__()
        self.issue = "contained a Boron atom"

    def _filter(self, mol: Mol) -> bool:
        """Returns True if molecule contains boron atom"""
        return not mol.HasSubstructMatch(MolFromSmarts("[#5]"))
