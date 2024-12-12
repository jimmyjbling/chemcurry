"""curation functions that involve formal charges"""

from rdkit.Chem import Mol, MolFromSmarts

from .base import Update


class Neutralize(Update):
    """
    Curation step to neutralize charged atoms in a molecule

    Will only neutralize atoms that can be neutralized
    by adding implicit hydrogen atoms

    Attributes
    ----------
    issue : str
        Description of issue related to the curation step
    note : str
        Description of what changes were made when a molecule was updated
    dependency : set
        Set containing the names of preceding required steps.
        If no dependency, will be an empty set.

    References
    ----------
    https://www.rdkit.org/docs/Cookbook.html#neutralizing-molecules
    """

    def __init__(self):
        super().__init__()
        self.note = "chemical neutralized"
        self.rank = 3

    def _update(self, mol: Mol) -> Mol:
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
        return mol
