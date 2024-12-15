"""test add charge steps"""

import pytest

from chemcurry.molecule import Molecule
from chemcurry.steps import Neutralize


@pytest.fixture
def molecules():
    """Fixture for molecules"""
    return [
        Molecule("mol1", None),
        Molecule.from_smiles("mol2", "CCCC(=O)[O-]"),
        Molecule.from_smiles("mol3", "CCCC(=O)O"),
    ]


@pytest.mark.unit
class TestNeutralize:
    """Test Neutralize curation step"""

    def test_neutralize(self, molecules):
        """Test that FlagBoron works as expected"""
        step = Neutralize()
        num_notes, num_issues = step(molecules)
        assert num_notes == 1
        assert num_issues == 0

        assert len(molecules[1].notes) == 1
