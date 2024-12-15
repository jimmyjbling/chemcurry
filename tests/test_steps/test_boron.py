"""test add boron steps"""

import pytest

from chemcurry.molecule import Molecule
from chemcurry.steps import FlagBoron


@pytest.fixture
def molecules():
    """Fixture for molecules"""
    return [
        Molecule("mol1", None),
        Molecule.from_smiles("mol2", "O=C(N[C@H](C(=O)N[C@H](B(O)O)CC(C)C)Cc1ccccc1)c2nccnc2"),
        Molecule.from_smiles("mol3", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
        Molecule.from_smiles("mol4", "CCCCBr"),
    ]


@pytest.mark.unit
class TestFlagBoron:
    """Test Add curation step"""

    def test_flag_boron(self, molecules):
        """Test that FlagBoron works as expected"""
        step = FlagBoron()
        num_notes, num_issues = step(molecules)
        assert num_notes == 0
        assert num_issues == 1

        assert molecules[1].issue != ""
