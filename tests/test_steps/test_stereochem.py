"""test add stereochem steps"""

import pytest

from chemcurry.molecule import Molecule
from chemcurry.steps import RemoveStereochem


@pytest.fixture
def molecules():
    """Fixture for molecules"""
    return [
        Molecule("mol1", None),
        Molecule.from_smiles("mol2", "C[C@H](N)C(=O)O"),
        Molecule.from_smiles("mol3", "C[C@@H](N)C(=O)O"),
        Molecule.from_smiles("mol4", "CC(N)C(=O)O"),
    ]


@pytest.mark.unit
class TestRemoveStereochem:
    """Test RemoveStereochem curation step"""

    def test_remove_stereochem(self, molecules):
        """Test that RemoveStereochem works as expected"""
        step = RemoveStereochem()
        num_notes, num_issues = step(molecules)
        assert num_notes == 2
        assert num_issues == 0

        assert len(molecules[1].notes) == 1
        assert len(molecules[2].notes) == 1
