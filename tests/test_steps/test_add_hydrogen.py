"""test add 3d curation step"""

import pytest
from rdkit.Chem.rdmolfiles import MolFromSmiles
from rdkit.Chem.rdmolops import AddHs

from chemcurry.molecule import Molecule
from chemcurry.steps import AddH


@pytest.fixture
def molecules():
    """Fixture for molecules to test Add3D on"""
    return [
        Molecule("mol1", None),
        Molecule("mol2", AddHs(MolFromSmiles("CCCC"))),
        Molecule.from_smiles("mol3", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ]


class TestAddH:
    """Test Add curation step"""

    def test_add_3d(self, molecules):
        """Test that Add3D works as expected"""
        step = AddH()
        num_notes, num_issues = step(molecules)
        assert num_notes == 1
        assert num_issues == 0
