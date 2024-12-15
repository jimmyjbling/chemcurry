"""test add hydrogen curation steps"""

import pytest
from rdkit.Chem.rdmolfiles import MolFromSmiles
from rdkit.Chem.rdmolops import AddHs

from chemcurry.molecule import Molecule
from chemcurry.steps import AddH


@pytest.fixture
def molecules():
    """Fixture for molecules"""
    return [
        Molecule("mol1", None),
        Molecule("mol2", AddHs(MolFromSmiles("CCCC"))),
        Molecule.from_smiles("mol3", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ]


@pytest.mark.unit
class TestAddH:
    """Test Add curation step"""

    def test_add_h(self, molecules):
        """Test that AddH works as expected"""
        step = AddH()
        num_notes, num_issues = step(molecules)
        assert num_notes == 1
        assert num_issues == 0

        assert len(molecules[2].notes) == 1
