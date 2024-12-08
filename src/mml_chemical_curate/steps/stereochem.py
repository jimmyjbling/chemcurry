"""stereochemistry curation functions"""

from copy import deepcopy

from rdkit.Chem.rdmolops import RemoveStereochemistry

from ..flags import CurationIssue, CurationNote
from .base import CurationStep


class CurateRemoveStereochem(CurationStep):
    """Removes stereochemistry from compounds"""

    def __init__(self):
        self.issue = CurationIssue.flatten_failed
        self.note = CurationNote.stereochem_removed
        self.rank = 3

    def _func(self, molecules):
        for mol in molecules:
            if mol.failed_curation:
                continue
            _tmp = deepcopy(mol.mol)
            RemoveStereochemistry(_tmp)  # this function is inplace for some reason
            mol.update_mol(_tmp, self.note)
