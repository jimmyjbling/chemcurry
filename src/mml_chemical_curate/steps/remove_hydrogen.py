"""removal of hydrogen curation functions"""

import importlib
from typing import Dict, Optional

from rdkit.Chem import Mol
from rdkit.Chem.rdmolops import RemoveAllHs, RemoveHs

from .base import SingleCurationStep, check_for_boost_rdkit_error


DEFAULT_REMOVE_HS_PARAMETERS = {
    "removeAndTrackIsotopes",
    "removeDefiningBondStereo",
    "removeDegreeZero",
    "removeDummyNeighbors",
    "removeHigherDegrees",
    "removeHydrides",
    "removeInSGroups",
    "removeIsotopes",
    "removeMapped",
    "removeNonimplicit",
    "removeNontetrahedralNeighbors",
    "removeOnlyHNeighbors",
    "removeWithQuery",
    "removeWithWedgedBond",
    "showWarnings",
    "updateExplicitCount",
}


class CurateRemoveH(SingleCurationStep):
    """Remove only non-essential explicit hydrogen atoms"""

    def __init__(self, remove_hs_params: Optional[Dict[str, bool]] = None):
        """
        Initialize the curation step

        Parameters
        ----------
        remove_hs_params: dict[str, bool], optional
            a dictionary of RemoveHsParameters keys and boolean values to pass to
            the RemoveHs function in rdkit
            If not passed, the default parameters are used.

        References
        ----------
        https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.RemoveHsParameters
        """
        super().__init__()
        self.issue = "failed to remove non-essential explict hydrogen atoms"
        self.note = "removed non-essential explicit hydrogen atoms"
        self.rank = 3
        self._load_remove_hs_params(remove_hs_params if remove_hs_params else {})

    def _load_remove_hs_params(self, params: Dict[str, bool]):
        """Make the object pickle-able"""
        self.remove_hs_params = importlib.import_module("rdkit.Chem.rdmolops").RemoveHsParameters
        for key, value in params.items():
            setattr(self.remove_hs_params, key, value)

    def _func(self, molecules):
        for mol in molecules:
            if mol.failed_curation:
                continue
            try:
                mol.update_mol(RemoveHs(Mol, self.remove_hs_params), self.note)
            except TypeError as e:
                if check_for_boost_rdkit_error(str(e)):
                    mol.flag_issue(self.issue)
                else:
                    raise e


class CurateRemoveAllH(SingleCurationStep):
    """
    Remove all explicit hydrogen atoms, even if they are required

    DANGER this can create new molecules if you are not careful, use sparingly
    """

    def __init__(self):
        self.issue = "failed to remove all explict hydrogen atoms"
        self.note = "removed all explict hydrogen atoms"
        self.rank = 3

    def _func(self, molecules):
        for mol in molecules:
            if mol.failed_curation:
                continue
            try:
                mol.update_mol(RemoveAllHs(Mol), self.note)
            except TypeError as e:
                if check_for_boost_rdkit_error(str(e)):
                    mol.flag_issue(self.issue)
                else:
                    raise e
