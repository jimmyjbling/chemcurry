"""3d curation steps"""

from func_timeout import FunctionTimedOut, func_timeout
from rdkit.Chem import Mol
from rdkit.Chem.rdDistGeom import EmbedMolecule, ETKDGv3

from ..flags import CurationIssue, CurationNote
from .base import CurationStep, check_for_boost_rdkit_error


def _add_3d(mol: Mol, timeout: int = 10) -> Mol:
    """
    Given a rdkit Mol, generate a random energy minimized 3D conformer, inplace

    Notes
    -----
    Will timeout after 10 seconds and fail to generate a 3D conformer.
    This is a limitation of RDKit sometimes hanging on this function call
    Will return None if RDKit cannot make the 3D pose in the given time limit

    Parameters
    ----------
    mol: Mol
        Mol to add 3D conformer to
    timeout: int
        time to wait before conformer generation fails

    Returns
    -------
    new_mol: Mol
    """
    ps = ETKDGv3()
    ps.useRandomCoords = True
    func_timeout(timeout, EmbedMolecule, (mol, ps))
    EmbedMolecule(mol, ps)
    mol.GetConformer()
    return mol


class CurateAdd3D(CurationStep):
    """
    Curation function to add 3D conformer to molecules

    Notes
    -----
    Uses the ETKDGv3 torsional angle potentials to do this
    """

    def __init__(self, timeout: int = 10):
        """
        Initialize the curation step

        Parameters
        ----------
        timeout: int, default=10
            time to wait before conformer generation fails
        """
        self.issue = CurationIssue.failed_gen_3d_conformer
        self.note = CurationNote.generated_3d_conformer
        self.rank = 3
        self.timeout = timeout
        self.dependency = {"CurateAddH"}

    def _func(self, molecules):
        for mol in molecules:
            if mol.failed_curation:
                continue
            try:
                mol.update_mol(_add_3d(mol.mol, self.timeout), self.note)
            except TypeError as e:
                if check_for_boost_rdkit_error(str(e)):
                    mol.flag_issue(CurationIssue.failed_gen_3d_conformer)
                else:
                    raise e
            except FunctionTimedOut:
                mol.flag_issue(CurationIssue.failed_gen_3d_conformer)
