"""test base classes for curation steps"""

import pytest
from rdkit import Chem
from rdkit.Chem import Mol

from chemcurry.molecule import Molecule
from chemcurry.steps.base import (
    DEFAULT_ISSUE,
    DEFAULT_NOTE,
    BaseCurationStep,
    CurationStepError,
    Filter,
    IssueMixin,
    NoteMixin,
    Update,
    check_for_boost_rdkit_error,
)


@pytest.fixture
def molecules():
    """Fixture for example Molecule"""
    return [
        Molecule("mol1", Mol()),
        Molecule("mol2", Chem.MolFromSmiles("CC")),
        Molecule("mol2", Chem.MolFromSmiles("CCCCO")),
        Molecule("mol2", Chem.MolFromSmiles("CCCCN")),
    ]


@pytest.mark.unit
@pytest.mark.parametrize(
    "error_message,expected",
    [
        ("Boost error detected in RDKit", True),
        ("random boost message unrelated to RDKit", True),
        ("RDKIT failed", False),
        ("Boost error", False),
    ],
)
def test_check_for_boost_rdkit_error(error_message, expected):
    """Test check_for_boost_rdkit_error"""
    assert check_for_boost_rdkit_error(error_message) == expected


@pytest.mark.unit
def test_curation_step_error():
    """Test CurationStepError"""
    with pytest.raises(CurationStepError):
        raise CurationStepError("Test error")


@pytest.mark.unit
def test_note_mixin():
    """Test NoteMixin"""
    mixin = NoteMixin()
    assert isinstance(mixin.note, str)
    mixin.note = "Note for {} step"
    assert mixin.get_note_text("TestCuration") == "Note for TestCuration step"


@pytest.mark.unit
def test_issue_mixin():
    """Test IssueMixin"""
    mixin = IssueMixin()
    assert isinstance(mixin.issue, str)
    mixin.issue = "Issue in {} curation"
    assert mixin.get_issue_text("TestStep") == "Issue in TestStep curation"


@pytest.fixture
def mock_base_step():
    """Fixture for mocking BaseCurationStep"""

    class MockCurationStep(BaseCurationStep):
        def __call__(self, chemicals):
            return 0, 0

    return MockCurationStep


@pytest.mark.unit
class TestBaseCurationStep:
    """unit tests for BaseCurationStep class"""

    @pytest.mark.filterwarnings("error")
    def test_base_curation_step(self, mock_base_step):
        """Test that BaseCurationStep can be instantiated with no warnings"""
        mock_base_step.issue = "Mock issue"
        mock_base_step.note = "Mock note"
        step = mock_base_step()
        assert str(step) == "MockCurationStep"
        assert repr(step) == "MockCurationStep"
        assert len(step.dependency) == 0


@pytest.fixture
def mock_filter():
    """Fixture for mocking Filter"""

    class MockFilter(Filter):
        issue = "Test issue"

        def _filter(self, mol):
            # simulate a filter that fails if 'O' in smiles
            return "O" not in Chem.MolToSmiles(mol)

    return MockFilter


class TestFilter:
    """test Filter class"""

    def test_filter_abstract(self):
        """Test that Filter is abstract"""
        with pytest.raises(TypeError):
            _ = Filter()

    def test_filter_is_callable(self, mock_filter):
        """Test that Filter is callable"""
        filter_ = mock_filter()
        assert callable(filter_)

    def test_broken_filter(self, mock_filter):
        """Test that Filter can't be instantiated with a note'"""
        mock_filter.note = "Test note"
        with pytest.raises(
            CurationStepError,
            match="Filter curation steps should not implement the 'note' attribute",
        ):
            mock_filter()

    def test_run_filter(self, mock_filter, molecules):
        """Test that Filter can be instantiated with no warnings"""
        filter_ = mock_filter()
        num_notes, num_issues = filter_(molecules)
        assert num_notes == 0
        assert num_issues == 1

        # first mol is invalid, second should pass, third should fail, forth passes
        assert molecules[0].failed_curation is True
        assert molecules[1].failed_curation is False
        assert molecules[2].failed_curation is True
        assert molecules[3].failed_curation is False

        # the issue in the mock is "Test issue"; only molecule2 should get this issue
        assert molecules[0].issue != "Test issue"
        assert molecules[2].issue == "Test issue"

    @pytest.mark.filterwarnings("error")
    def test_default_issue_warning(self, mock_filter):
        """Test that Filter will raise warning if issue description is default"""
        mock_filter.issue = DEFAULT_ISSUE.format(mock_filter.__name__)
        with pytest.warns(UserWarning, match="using default issue description"):
            _ = mock_filter()

    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_invalid_issue(self, mock_filter):
        """Test that Filter will raise error if issue description is invalid"""
        mock_filter.issue = 123
        with pytest.raises(
            CurationStepError, match=r"CurationSteps require that the `issue` attribute is a str;"
        ):
            _ = mock_filter()


@pytest.fixture
def mock_update():
    """Fixture for mocking Update"""

    class MockUpdate(Update):
        issue = "Test issue"
        note = "Test note"

        def _update(self, mol):
            # simulate a filter that changes 'O' to 'N' and fails when 'N' is present
            _smi = Chem.MolToSmiles(mol)
            if "O" in _smi:
                return Chem.MolFromSmiles(_smi.replace("O", "N"))
            elif "N" in _smi:
                return None
            else:
                return mol

    return MockUpdate


class TestUpdate:
    """test Filter class"""

    def test_update_abstract(self):
        """Test that Update is abstract"""
        with pytest.raises(TypeError):
            _ = Update()

    def test_update_is_callable(self, mock_update):
        """Test that Update is callable"""
        update = mock_update()
        assert callable(update)

    def test_run_update(self, mock_update, molecules):
        """Test that Filter can be instantiated with no warnings"""
        update = mock_update()
        num_notes, num_issues = update(molecules)
        assert num_notes == 1
        assert num_issues == 1

        # first mol is invalid, second should pass, third should pass, fourth should fail
        assert molecules[0].failed_curation is True
        assert molecules[1].failed_curation is False
        assert molecules[2].failed_curation is False
        assert molecules[3].failed_curation is True

        # check that issues are set right
        assert molecules[0].issue != "Test issue"
        assert molecules[3].issue == "Test issue"

        # check that notes are set right
        assert len(molecules[1].notes) == 0
        assert molecules[2].notes[-1] == "Test note"

    @pytest.mark.filterwarnings("error")
    def test_default_issue_warning(self, mock_update):
        """Test that Update will raise warning if issue description is default"""
        mock_update.issue = DEFAULT_ISSUE.format(mock_update.__name__)
        with pytest.warns(UserWarning, match="using default issue description"):
            _ = mock_update()

    @pytest.mark.filterwarnings("error")
    def test_default_note_warning(self, mock_update):
        """Test that Update will raise warning if note description is default"""
        mock_update.note = DEFAULT_NOTE.format(mock_update.__name__)
        with pytest.warns(UserWarning, match="using default note description"):
            _ = mock_update()

    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_invalid_issue(self, mock_update):
        """Test that Update will raise error if issue description is invalid"""
        mock_update.issue = 123
        with pytest.raises(
            CurationStepError, match=r"CurationSteps require that the `issue` attribute is a str;"
        ):
            _ = mock_update()

    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_invalid_note(self, mock_update):
        """Test that Update will raise error if note description is invalid"""
        mock_update.note = 123
        with pytest.raises(
            CurationStepError, match=r"CurationSteps require that the `note` attribute is a str;"
        ):
            _ = mock_update()
