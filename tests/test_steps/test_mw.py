"""test add mw steps"""

import pytest

from chemcurry.molecule import Molecule
from chemcurry.steps import FilterMW


@pytest.fixture
def molecules():
    """Fixture for molecules"""
    return [
        Molecule("mol1", None),
        Molecule.from_smiles(
            "mol2",
            "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H]"
            "(OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)"
            "C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H]"
            "(NC(=O)c6ccccc6)c7ccccc7)O)(C)C",
        ),
        Molecule.from_smiles("mol3", "CCCC"),
    ]


@pytest.mark.unit
class TestFilterMW:
    """Test FilterMW curation step"""

    def test_filter_mw_0_lower_bound(self):
        """Test that FilterMW raises an error with min_mw=0"""
        with pytest.raises(ValueError, match="min_mw must be greater than 0"):
            FilterMW(min_mw=0)

    def test_filter_mw_bound_invalid(self):
        """Test that FilterMW raises an error with min_mw > man_mw"""
        with pytest.raises(ValueError, match="min_mw cannot be larger than man_mw"):
            FilterMW(min_mw=10, max_mw=5)

    def test_filter_mw(self, molecules):
        """Test that FlagBoron works as expected"""
        step = FilterMW()
        num_notes, num_issues = step(molecules)
        assert num_notes == 0
        assert num_issues == 0

    def test_filter_mw_small_upper(self, molecules):
        """Test that FlagBoron works as expected with small upper bound"""
        step = FilterMW(max_mw=10)
        num_notes, num_issues = step(molecules)
        assert num_notes == 0
        assert num_issues == 2

        assert molecules[1].issue != ""
        assert molecules[2].issue != ""

    def test_filter_mw_large_upper(self, molecules):
        """Test that FlagBoron works as expected with large upper bound"""
        step = FilterMW(max_mw=10000)
        num_notes, num_issues = step(molecules)
        assert num_notes == 0
        assert num_issues == 0

    def test_filter_mw_small_lower(self, molecules):
        """Test that FlagBoron works as expected small lower bound"""
        step = FilterMW(min_mw=100)
        num_notes, num_issues = step(molecules)
        assert num_notes == 0
        assert num_issues == 1

        assert molecules[2].issue != ""

    def test_filter_mw_large_lower(self, molecules):
        """Test that FlagBoron works as expected large lower bound"""
        step = FilterMW(min_mw=10000)
        num_notes, num_issues = step(molecules)
        assert num_notes == 0
        assert num_issues == 2

        assert molecules[1].issue != ""
        assert molecules[2].issue != ""
