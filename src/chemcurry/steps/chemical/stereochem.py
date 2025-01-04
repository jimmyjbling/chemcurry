"""stereochemistry curation functions"""

from rdkit.Chem import Mol
from rdkit.Chem.rdmolops import RemoveStereochemistry

from ..base import Update


class RemoveStereochem(Update):
    """
    Removes all specified stereochemistry from molecules

    Notes
    -----
    This converts compounds into racemic mixtures in the eyes of RDKIT

    Attributes
    ----------
    issue : str
        Description of issue related to the curation step
    note : str
        Description of what changes were made when a molecule was updated
    dependency : set
        Set containing the names of preceding required steps.
        If no dependency, will be an empty set.
    """

    def __init__(self):
        """Initialize the curation step"""
        self.issue = "failed to removed stereochemistry from chemical"
        self.note = "all stereochemistry are removed from chemical"

    def _update(self, mol: Mol) -> Mol:
        # this rdkit function is inplace for some reason
        RemoveStereochemistry(mol)
        return mol
