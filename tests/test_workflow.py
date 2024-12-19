"""Test workflow module"""

import json

import numpy as np
import pytest
from rdkit.Chem import Mol
from rdkit.Chem.rdmolfiles import MolFromSmiles

from chemcurry.steps import FlagInorganic, FlagMixtures, RemoveStereochem
from chemcurry.workflow import CurationWorkflow, CurationWorkflowError


@pytest.fixture
def workflow():
    """Workflow fixture"""
    step1 = RemoveStereochem()
    step2 = FlagInorganic()
    step3 = FlagMixtures()
    return CurationWorkflow([step1, step2, step3])


@pytest.fixture
def smiles():
    """Smiles fixture"""
    return [
        "None",
        "[Ni+2].[Cl-].[Cl-]",
        "CCCC(=O)O",
        "CCCCCC",
        "CCCCCC.[H]",
        "CCCCCC.CCO",
        "CCO.CCCCCNCCCC",
        "C[C@H](N)C(=O)O",
        "C[C@@H](N)C(=O)O",
        "CC(N)C(=O)O",
    ]


@pytest.fixture
def molecules():
    """Molecules fixture"""
    return [
        MolFromSmiles(_)
        for _ in [
            "None",
            "[Ni+2].[Cl-].[Cl-]",
            "CCCC(=O)O",
            "CCCCCC",
            "CCCCCC.[H]",
            "CCCCCC.CCO",
            "CCO.CCCCCNCCCC",
            "C[C@H](N)C(=O)O",
            "C[C@@H](N)C(=O)O",
            "CC(N)C(=O)O",
        ]
    ]


@pytest.mark.unit
@pytest.fixture
def curated_molecule_set(workflow, smiles):
    """Curated molecule set fixture"""
    return workflow.curate_smiles(smiles)


class TestWorkflow:
    """Test CurationWorkflow class"""

    def test_init(self, workflow):
        """Test that CurationWorkflow is initialized correctly"""
        assert workflow.name == "NA"
        assert workflow.description == "NA"
        assert workflow.repo_url == "NA"

    def test_order_warning(self, workflow):
        """Test that warning is raised when steps are out of order"""
        with pytest.warns(UserWarning):
            CurationWorkflow([FlagInorganic(), RemoveStereochem()])

    @pytest.mark.filterwarnings("error")
    def test_supress_order_warning(self, workflow):
        """Test that a warning is not raised when it is suppressed"""
        CurationWorkflow([FlagInorganic(), RemoveStereochem()], suppress_warnings=True)

    def test_name_property(self, workflow):
        """Test that name property works as expected"""
        # Test set and get
        workflow.name = "Test Name"
        assert workflow.name == "Test Name"
        # Test delete
        del workflow.name
        assert workflow.name == "NA"

    def test_description_property(self, workflow):
        """Test that description property works as expected"""
        # Test set and get
        workflow.description = "Test Description"
        assert workflow.description == "Test Description"
        # Test delete
        del workflow.description
        assert workflow.description == "NA"

    def test_repo_url_property(self, workflow):
        """Test that repo_url property works as expected"""
        # Test set and get
        workflow.repo_url = "https://example.com/repo"
        assert workflow.repo_url == "https://example.com/repo"
        # Test delete
        del workflow.repo_url

    def test_save_and_load_workflow_file(self, workflow, tmpdir):
        """Test that saving and loading a workflow file works as expected"""
        # Save the workflow to a temporary file
        file_path = tmpdir.join("workflow.json")
        workflow.save_workflow_file(file_path)

        # Check that the file exists and contains expected data (mocked step here)
        assert file_path.exists()

        loaded_workflow = CurationWorkflow.load(file_path)

        assert len(loaded_workflow.steps) == 3
        assert isinstance(loaded_workflow.steps[0], RemoveStereochem)
        assert isinstance(loaded_workflow.steps[1], FlagInorganic)
        assert isinstance(loaded_workflow.steps[2], FlagMixtures)

        assert loaded_workflow.name == workflow.name
        assert loaded_workflow.description == workflow.description
        assert loaded_workflow.repo_url == workflow.repo_url

    def test_load_workflow_with_unknown_step(self, workflow, tmpdir):
        """Test that loading a workflow file with an unknown step raises an error"""
        file_path = tmpdir.join("workflow.json")
        workflow.save_workflow_file(file_path)

        # create a fake step
        data = json.load(open(file_path, "r"))
        data["steps"]["0"]["name"] = "MadeUpCurationStep"
        json.dump(data, open(file_path, "w"))

        with pytest.raises(CurationWorkflowError, match="could not find curation step"):
            CurationWorkflow.load(file_path)

    def test_load_workflow_with_wrong_position(self, workflow, tmpdir):
        """Test that loading a workflow file with a step at the wrong position raises an error"""
        file_path = tmpdir.join("workflow.json")
        workflow.save_workflow_file(file_path)

        # create a fake step
        data = json.load(open(file_path, "r"))
        data["steps"]["5"] = data["steps"]["0"]
        json.dump(data, open(file_path, "w"))

        with pytest.raises(CurationWorkflowError, match="steps and positions"):
            CurationWorkflow.load(file_path)

    def test_curate_smiles(self, workflow, smiles):
        """Test that curate_smiles works as expected"""
        workflow.curate_smiles(smiles)

    def test_curate_mols(self, workflow, molecules):
        """Test that curate_mols works as expected"""
        workflow.curate_mols(molecules)


class TestCuratedMoleculeSet:
    """Test CuratedMoleculeSet class"""

    def test_to_smiles(self, curated_molecule_set):
        """Test that to_smiles works as expected"""
        assert curated_molecule_set.to_smiles() == [
            "CCCC(=O)O",
            "CCCCCC",
            "CC(N)C(=O)O",
            "CC(N)C(=O)O",
            "CC(N)C(=O)O",
        ]

    def test_to_smiles_with_failed(self, curated_molecule_set):
        """Test that to_smiles works as expected"""
        assert curated_molecule_set.to_smiles(include_failed=True) == [
            "",
            "[Cl-].[Cl-].[Ni+2]",
            "CCCC(=O)O",
            "CCCCCC",
            "CCCCCC.[H]",
            "CCCCCC.CCO",
            "CCCCCNCCCC.CCO",
            "CC(N)C(=O)O",
            "CC(N)C(=O)O",
            "CC(N)C(=O)O",
        ]

    def test_to_mols(self, curated_molecule_set):
        """Test that to_mols works as expected"""
        mols = curated_molecule_set.to_mols()

        assert len(mols) == 5
        assert isinstance(mols[0], Mol)

    def test_to_pandas(self, curated_molecule_set):
        """Test that to_pandas works as expected"""
        df = curated_molecule_set.to_pandas()

        assert "id" in df.columns
        assert "smiles" in df.columns
        assert "mol" in df.columns
        assert "issue" not in df.columns
        assert "notes" not in df.columns
        assert "passed" not in df.columns

        df = curated_molecule_set.to_pandas(
            include_failed=True, include_issues=True, include_notes=True
        )

        assert "id" in df.columns
        assert "smiles" in df.columns
        assert "mol" in df.columns
        assert "issue" in df.columns
        assert "notes" in df.columns
        assert "passed" in df.columns

    def test_get_passing_mask(self, curated_molecule_set):
        """Test that get_passing_mask works as expected"""
        mask = curated_molecule_set.get_passing_mask()
        assert mask == [False, False, True, True, False, False, False, True, True, True]

        mask = curated_molecule_set.get_passing_mask(as_numpy=True)
        assert isinstance(mask, np.ndarray)

    def test_get_num_issues_at_step(self, curated_molecule_set):
        """Test that get_num_issues_at_step works as expected"""
        assert curated_molecule_set.get_num_issues_at_step(0) == 1
        assert curated_molecule_set.get_num_issues_at_step(3) == 3
        assert curated_molecule_set.get_num_issues_at_step("FlagInorganic")[0] == 1
        assert curated_molecule_set.get_num_issues_at_step("FlagMixtures")[0] == 3

    def test_get_num_notes_at_step(self, curated_molecule_set):
        """Test that get_num_notes_at_step works as expected"""
        assert curated_molecule_set.get_num_notes_at_step(0) == 0
        assert curated_molecule_set.get_num_notes_at_step(1) == 2
        assert curated_molecule_set.get_num_notes_at_step("RemoveStereochem")[0] == 2

    def test_get_num_remaining_molecules_after_step(self, curated_molecule_set):
        """Test that get_num_remaining_molecules_after_step works as expected"""
        assert curated_molecule_set.get_num_remaining_molecules_after_step(0) == 9
        assert curated_molecule_set.get_num_remaining_molecules_after_step(1) == 9
        assert curated_molecule_set.get_num_remaining_molecules_after_step(2) == 8
        assert curated_molecule_set.get_num_remaining_molecules_after_step(3) == 5

    def test_write_report(self, curated_molecule_set, tmpdir):
        """Test that write_report works as expected"""
        file_path = tmpdir.join("curation_report.txt")
        curated_molecule_set.write_report(file_path)
        assert file_path.exists()

        lines = open(file_path, "r").read()
        assert "ChemCurry curation report" in lines

    def test_save(self, curated_molecule_set, tmpdir):
        """Test that save works as expected"""
        file_path = tmpdir.join("curated_molecules.pkl")
        curated_molecule_set.save(file_path)
        assert file_path.exists()

    def test_save_as_txt(self, curated_molecule_set, tmpdir):
        """Test that save_as_txt works as expected"""
        file_path = tmpdir.join("curated_molecules.txt")
        curated_molecule_set.save_as_txt(file_path)
        assert file_path.exists()

    def test_save_as_json(self, curated_molecule_set, tmpdir):
        """Test that save_as_json works as expected"""
        file_path = tmpdir.join("curated_molecules.json")
        curated_molecule_set.save_as_json(file_path)
        assert file_path.exists()

    def test_save_as_csv(self, curated_molecule_set, tmpdir):
        """Test that save_as_csv works as expected"""
        file_path = tmpdir.join("curated_molecules.csv")
        curated_molecule_set.save_as_csv(file_path)
        assert file_path.exists()

    def test_save_as_pandas(self, curated_molecule_set, tmpdir):
        """Test that save_as_pandas works as expected"""
        file_path = tmpdir.join("curated_molecules.pd.pkl")
        curated_molecule_set.save_as_pandas(file_path)
        assert file_path.exists()
