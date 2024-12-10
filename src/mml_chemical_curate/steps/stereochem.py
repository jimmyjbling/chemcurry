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

    def _func(self, chemical):
        # this rdkit function is inplace for some reason
        _tmp = deepcopy(chemical.mol)
        RemoveStereochemistry(_tmp)
        chemical.update_mol(_tmp, self.get_note_text())
