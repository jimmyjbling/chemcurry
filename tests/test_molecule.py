"""test the molecule.py module"""

import pytest
from rdkit import Chem

from chemcurry.molecule import Molecule, SmilesMixin


@pytest.fixture
def valid_smiles():
    """Get a valid smiles string for testing"""
    return "CCC"


@pytest.fixture
def valid_mol(valid_smiles):
    """Covert the valid_smiles to a mol object for testing"""
    return Chem.MolFromSmiles(valid_smiles)


@pytest.fixture
def invalid_mol():
    """Get an invalid mol object for testing"""
    return None


@pytest.fixture
def valid_mol2():
    """Get a second but unique valid mol object for testing"""
    return Chem.MolFromSmiles("CCOCC")


@pytest.mark.unit
class TestMolecule:
    """test the Molecule class"""

    def test_molecule_init(self, valid_mol, invalid_mol, valid_mol2):
        """Test the __init__ method of the class"""
        # check that empty call fails
        with pytest.raises(TypeError):
            _ = Molecule()

        # check that valid mol is loaded
        molecule1 = Molecule("valid", valid_mol)

        # check parameters set currently
        assert molecule1.mol is not None
        assert molecule1.mol == valid_mol
        assert molecule1.id_ == "valid"
        assert len(molecule1.notes) == 0
        assert molecule1.issue == ""
        assert isinstance(molecule1._current_hash, int)
        assert len(molecule1.mol_history) == 0
        assert molecule1.failed_curation is False
        assert molecule1.track_history is False

        molecule2 = Molecule("invalid", invalid_mol)
        assert molecule2.mol is not None
        assert len(molecule2.mol.GetAtoms()) == 0
        assert molecule2.failed_curation is True

    def test_track_history(self, valid_mol):
        """Test the track_history property of the class"""
        molecule = Molecule("valid", valid_mol, track_history=True)
        assert molecule.track_history is True
        assert len(molecule.mol_history) == 0
        with pytest.raises(
            RuntimeError, match="'track_history' cannot be change after object initialization"
        ):
            molecule.track_history = False

    def test_update_mol(self, valid_mol, valid_mol2, invalid_mol):
        """Test the update_mol method of the class"""
        # update to new valid mol
        molecule = Molecule("valid", valid_mol)
        _first_hash = molecule._current_hash
        _update = molecule.update_mol(valid_mol2, "update to mol")
        assert _update is True
        assert molecule.mol is not None
        assert molecule.id_ == "valid"
        assert molecule.mol == valid_mol2
        assert molecule._current_hash != _first_hash
        assert len(molecule.notes) == 1
        assert molecule.notes[-1] == "update to mol"
        assert molecule.failed_curation is False
        # check that no update return false:
        _update = molecule.update_mol(valid_mol2, "update to mol")
        assert _update is False

        # update to new valid mol with history tracking
        molecule = Molecule("valid", valid_mol, track_history=True)
        _first_hash = molecule._current_hash
        molecule.update_mol(valid_mol2, "update to mol")
        assert len(molecule.mol_history) == 1
        # these should be different object because of deepcopy, but same contents
        assert molecule.mol_history[0] != valid_mol
        assert hash(molecule.mol_history[0].ToBinary()) == _first_hash

        # update with invalid mol
        with pytest.raises(
            ValueError,
            match="if molecule becomes invalid, should be caught "
            "and flagged with issue by curation step; 'None'",
        ):
            molecule.update_mol(invalid_mol, "update to mol None")
        with pytest.raises(
            ValueError,
            match="if molecule becomes invalid, should be caught "
            "and flagged with issue by curation step; 'Empty Mol'",
        ):
            molecule.update_mol(Chem.MolFromSmiles(""), "update to mol empty")

    def test_flag_issue(self, valid_mol):
        """Test the flag_issue method of the class"""
        molecule = Molecule("valid", valid_mol)
        molecule.flag_issue("issue")
        assert molecule.issue == "issue"
        assert molecule.failed_curation is True


@pytest.fixture
def dummy_smiles_mixin():
    """Provide a concrete implementation of SmilesMixin for testing"""

    class DummySmilesMixin(SmilesMixin):
        """dummy smiles mixin class for testing"""

        def __init__(self, id_, mol):
            """Initialize the class"""
            self.id_ = id_
            self.mol = mol

    return DummySmilesMixin


@pytest.fixture
def dummy_extra_args_smiles_mixin():
    """Provide a concrete implementation of SmilesMixin for testing with kwargs"""

    class DummySmilesMixinExtraArg(SmilesMixin):
        """dummy smiles mixin class for testing with and extra argument"""

        def __init__(self, id_, mol, extra_arg=None):
            """Initialize the class"""
            self.id_ = id_
            self.mol = mol
            self.extra_arg = extra_arg

    return DummySmilesMixinExtraArg


class TestSmilesMixin:
    """test the SmilesMixin class"""

    def test_is_abstract(self):
        """Test that the SmilesMixin class is abstract"""
        with pytest.raises(TypeError):
            _ = SmilesMixin()

    def test_from_smiles(
        self, dummy_smiles_mixin, dummy_extra_args_smiles_mixin, valid_smiles, valid_mol
    ):
        """Test the from_smiles class method"""
        mol = dummy_smiles_mixin.from_smiles(0, valid_smiles)
        assert mol.id_ == 0
        assert hash(mol.mol.ToBinary()) == hash(valid_mol.ToBinary())
        assert Chem.MolToSmiles(mol.mol) == valid_smiles

        mol = dummy_extra_args_smiles_mixin.from_smiles(0, valid_smiles, extra_arg=10)
        assert mol.extra_arg == 10

    def test_get_smiles(self, dummy_smiles_mixin, valid_smiles, valid_mol):
        """Test the get_smiles method"""
        mol = dummy_smiles_mixin(0, valid_mol)
        assert mol.get_smiles() == valid_smiles

    def test_has_same_smiles(self, dummy_smiles_mixin, valid_smiles, valid_mol):
        """Test the has_same_smiles method"""
        mol1 = dummy_smiles_mixin(0, valid_mol)
        mol2 = dummy_smiles_mixin(1, valid_mol)
        assert mol1.has_same_smiles(mol2.mol) is True
