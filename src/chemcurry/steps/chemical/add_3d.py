"""3d curation steps"""

from typing import Optional

from func_timeout import FunctionTimedOut, func_timeout
from rdkit.Chem import Mol
from rdkit.Chem.rdDistGeom import EmbedMolecule, ETKDGv3

from ..base import Update, check_for_boost_rdkit_error


class Add3D(Update):
    """
    Curation step to add 3D conformer to molecules

    Notes
    -----
    Uses the ETKDGv3 torsional angle potentials to add a 3D conformer
    to molecules. If the conformer generation exceeds the timeout duration
    or fails due to an internal RDKit error, the molecule is flagged with an issue.

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

    def __init__(self, timeout: int = 10):
        """
        Initialize the curation step

        Parameters
        ----------
        timeout: int, default=10
            time to wait before conformer generation fails
        """
        self.issue = "failed to generate a 3D conformer"
        self.note = "generated a 3D conformer using ETKDGv3"
        self.timeout = timeout
        self.dependency = {"CurateAddH"}

    def _update(self, mol: Mol) -> Optional[Mol]:
        """Attempts to add 3D conformer to molecule; returns None if fails"""
        try:
            ps = ETKDGv3()
            ps.useRandomCoords = True
            func_timeout(self.timeout, EmbedMolecule, (mol, ps))
            EmbedMolecule(mol, ps)
            if len(mol.GetConformers()) == 0:
                return None
            else:
                return mol
        except TypeError as e:
            if check_for_boost_rdkit_error(str(e)):
                return None
            else:
                raise e
        except FunctionTimedOut:
            return None
