"""stereochemistry curation functions"""

from copy import deepcopy

from rdkit.Chem.rdmolops import RemoveStereochemistry

from .base import SingleCurationStep


class CurateRemoveStereochem(SingleCurationStep):
    """Removes stereochemistry from chemicals"""

    def __init__(self):
        self.issue = "failed to removed stereochemistry from chemical"
        self.note = "all stereochemistry are removed from chemical"
        self.rank = 3

    def _func(self, molecules):
        for mol in molecules:
            if mol.failed_curation:
                continue
            _tmp = deepcopy(mol.mol)
            RemoveStereochemistry(_tmp)  # this function is inplace for some reason
            mol.update_mol(_tmp, self.note)
