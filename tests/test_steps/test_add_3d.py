"""test add 3d curation steps"""

import pytest
from rdkit.Chem.rdmolops import AddHs

from chemcurry.molecule import Molecule
from chemcurry.steps import Add3D


@pytest.fixture
def molecules():
    """Fixture for molecules to test Add3D on"""
    _mols = [
        Molecule("mol1", None),
        Molecule.from_smiles("mol2", "CCCCCC"),
        Molecule.from_smiles("mol3", "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"),
    ]

    # need to add H to stop rdkit from complaining
    for _mol in _mols:
        _mol.mol = AddHs(_mol.mol)

    return _mols


@pytest.mark.unit
class TestAdd3D:
    """Test Add3D curation step"""

    def test_add_3d(self, molecules):
        """Test that Add3D works as expected"""
        step = Add3D(timeout=30)
        num_notes, num_issues = step(molecules)
        assert num_notes == 2
        assert num_issues == 0

    def test_add_3d_timeout(self, molecules):
        """Test that Add3D times out with issues"""
        step = Add3D(timeout=0)
        num_notes, num_issues = step(molecules)
        assert num_notes == 0
        assert num_issues == 2
