"""curation functions that involve formal charges"""

from copy import deepcopy

from rdkit.Chem import Mol, MolFromSmarts

from .base import SingleCurationStep


def neutralize_mol(mol: Mol) -> Mol:
    """
    Removes any neutralize-able charge on the molecule in place.

    References
    ----------
    https://www.rdkit.org/docs/Cookbook.html#neutralizing-molecules

    Notes
    -----
    Will remove charge regardless of pka and pH

    Parameters
    ----------
    mol: Mol
        Mol to neutralize

    Returns
    -------
    neutralized_mol: Mol
    """
    _mol = deepcopy(mol)
    pattern = MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),$([!B&-1])!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            h_count = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(h_count - chg)
            atom.UpdatePropertyCache()
    return _mol


class CurateNeutralize(SingleCurationStep):
    """make a charged molecule neutral"""

    def __init__(self):
        super().__init__()
        self.issue = "failed to neutralize chemical"
        self.note = "chemical neutralized"
        self.rank = 3

    def _func(self, chemical):
        chemical.update_mol(neutralize_mol(chemical.mol), self.get_note_text())
