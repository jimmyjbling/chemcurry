"""test add hydrogen curation steps"""

import pytest
from rdkit.Chem.rdmolfiles import MolFromSmiles
from rdkit.Chem.rdmolops import AddHs

from chemcurry.molecule import Molecule
from chemcurry.steps import RemoveAllHs, RemoveHs
from chemcurry.steps.remove_hydrogen import DEFAULT_REMOVE_HS_PARAMETERS


@pytest.fixture
def molecules():
    """Fixture for molecules"""
    return [
        Molecule("mol1", None),
        Molecule("mol2", AddHs(MolFromSmiles("CCCC"))),
        Molecule.from_smiles("mol3", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
        Molecule.from_smiles("mol4", "CC(C)Cc1ccc(cc1)C(C)C(=O)O.[H]"),
        Molecule.from_smiles("mol5", "CCC(=O)[NH2]"),
    ]


@pytest.mark.unit
class TestRemoveAllH:
    """Test remove all unnecessary hydrogen step"""

    def test_remove_all_h(self, molecules):
        """Test that AddH works as expected"""
        step = RemoveAllHs()
        num_notes, num_issues = step(molecules)
        assert num_notes == 2
        assert num_issues == 0

        assert len(molecules[1].notes) == 1
        assert len(molecules[3].notes) == 1


@pytest.mark.unit
class TestRemoveH:
    """Test remove custom hydrogen step"""

    def test_remove_h(self, molecules):
        """Test that remove H works as expected"""
        step = RemoveHs()
        num_notes, num_issues = step(molecules)
        assert num_notes == 1
        assert num_issues == 0

        assert len(molecules[1].notes) == 1

    def test_remove_h_with_custom_param(self, molecules):
        """Test that remove H with custom settings works as expected"""
        step = RemoveHs({"removeDegreeZero": True})
        num_notes, num_issues = step(molecules)
        assert num_notes == 2
        assert num_issues == 0

        assert len(molecules[1].notes) == 1
        assert len(molecules[3].notes) == 1

    def test_get_h_params(self):
        """Test that get_remove_h_parameters works as expected"""
        step = RemoveHs()
        step.get_remove_h_parameters()
        for key in step.get_remove_h_parameters().keys():
            assert key in DEFAULT_REMOVE_HS_PARAMETERS
