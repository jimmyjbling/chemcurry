"""mixture based curation functions"""

import importlib

from rdkit.Chem import GetMolFrags, Mol

from .base import Filter, Update


class FlagMixtures(Filter):
    """
    Flag compounds that have a mixture

    Mixtures are defined as molecules with more than one disconnected fragment.
    In SMILES is this is usually represented by '.' characters.

    Attributes
    ----------
    issue : str
        Description of issue related to the curation step
    dependency : set
        Set containing the names of preceding required steps.
        If no dependency, will be an empty set.
    """

    def __init__(self):
        super().__init__()
        self.issue = "chemical contains a mixture"
        self.dependency = {"RemoveH|RemoveAllH"}

    def _filter(self, mol: Mol) -> bool:
        """Returns True if molecule contains a mixture"""
        return len(GetMolFrags(mol)) > 1


class DemixLargestFragment(Update):
    """
    Update molecules as the largest component of a mixture

    If a compound is not a mixture will not make any changes

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
        super().__init__()
        self.issue = "failed to de-mix the chemical"
        self.note = "de-mixed by picking largest chemical compound"

        # makes this object pickle-able
        self._chooser = importlib.import_module(
            "rdkit.Chem.MolStandardize.rdMolStandardize"
        ).LargestFragmentChooser()

    def _update(self, mol: Mol) -> Mol:
        """Returns the largest component of a mixture"""
        return self._chooser.choose(mol)
