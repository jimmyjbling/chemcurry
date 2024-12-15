"""test sanitize steps"""

import pytest

from chemcurry.molecule import Molecule
from chemcurry.steps import SanitizeMolecule


@pytest.fixture
def molecules():
    """Fixture for molecules"""
    return [
        Molecule("mol1", None),
        Molecule.from_smiles("mol2", "CCCCO"),
        Molecule.from_smiles("mol3", "CCCC(=O)O"),
    ]


@pytest.mark.unit
class TestSanitize:
    """Test Sanitize Molecule curation step"""

    def test_sanitize(self, molecules):
        """Test that Sanitize works as expected"""
        step = SanitizeMolecule()
        num_notes, num_issues = step(molecules)
        assert num_notes == 0
        assert num_issues == 0
