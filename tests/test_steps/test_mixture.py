"""test add mixture steps"""

import pytest

from chemcurry.molecule import Molecule
from chemcurry.steps import DemixLargestFragment, FlagMixtures


@pytest.fixture
def molecules():
    """Fixture for molecules"""
    return [
        Molecule("mol1", None),
        Molecule.from_smiles("mol2", "CCCCCC"),
        Molecule.from_smiles("mol3", "CCCCCC.[H]"),
        Molecule.from_smiles("mol4", "CCCCCC.CCO"),
        Molecule.from_smiles("mol5", "CCO.CCCCCNCCCC"),
    ]


@pytest.mark.unit
class TestFlagMixture:
    """Test FlagMixture curation step"""

    def test_flag_mixture(self, molecules):
        """Test that FlagBoron works as expected"""
        step = FlagMixtures()
        num_notes, num_issues = step(molecules)
        assert num_notes == 0
        assert num_issues == 3

        assert molecules[2].issue != ""
        assert molecules[3].issue != ""
        assert molecules[4].issue != ""


@pytest.mark.unit
class TestDemixLargestFragment:
    """Test DemixLargestFragment curation step"""

    def test_demix_largest_fragment(self, molecules):
        """Test that DemixLargestFragment works as expected"""
        step = DemixLargestFragment()
        num_notes, num_issues = step(molecules)
        assert num_notes == 3
        assert num_issues == 0

        assert len(molecules[2].notes) == 1
        assert len(molecules[3].notes) == 1
        assert len(molecules[4].notes) == 1

        assert molecules[2].get_smiles() == "CCCCCC"
        assert molecules[3].get_smiles() == "CCCCCC"
        assert molecules[4].get_smiles() == "CCCCCNCCCC"
