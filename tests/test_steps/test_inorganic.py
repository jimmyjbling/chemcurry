"""test add inorganic steps"""

import pytest

from chemcurry.molecule import Molecule
from chemcurry.steps import FlagInorganic


@pytest.fixture
def molecules():
    """Fixture for molecules"""
    return [
        Molecule("mol1", None),
        Molecule.from_smiles("mol2", "[Ni+2].[Cl-].[Cl-]"),
        Molecule.from_smiles("mol3", "CCCC(=O)O"),
    ]


@pytest.mark.unit
class TestFlagInorganic:
    """Test FlagInorganic curation step"""

    def test_flag_inorganic(self, molecules):
        """Test that FlagBoron works as expected"""
        step = FlagInorganic()
        num_notes, num_issues = step(molecules)
        assert num_notes == 0
        assert num_issues == 1

        assert molecules[1].issue != ""
