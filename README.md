<picture>
  <source media="(prefers-color-scheme: dark)" srcset="./ChemCurry-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="./ChemCurry-light.svg">
  <img src="./ChemCurry-light.svg?raw=true" alt="ChemCurry logo">
</picture>

# ChemCurry

`chemcurry` is a chemical curation workflow package meant to both streamline building
curation workflows and producing detailed reports about which chemicals where flagged
and when *while* doing so in a manner to enforce reproducibility and easy sharing.
The Molecular Modeling Lab @ UNC often finds itself needing
to generate these reports to show to our PI and share our workflows with new members,
so this package was developed as a way to standardize that process.

While most chemical curation workflows for any project can be built in under 100 lines of code,
the core idea behind `chemcurry` is to assert reproducibility and easy building/sharing
among chemist with any level of coding background. Most cheminformatics projects and
publications will need to do some type of curation, and, frankly, the methods on how
this is done is often not up to par with scientific reproducibility standards.
We believe that lack of reproducibility hurts our filed and `chemcurry` aims to
fix that (for at least on part of it).

Closely related to the philosophy of reproducibility, `chemcurry` was also designed to be
easy to add new curation functions too. There is a simple API that really only requires
you to write the same code you might if you were doing in manually in a notebook or script.

### What about curation with labels or non-chemical properties?

`chemcurry` is designed to operate on explict chemical properties, meaning if the property cannot
be calculated using just the chemical, it will not fit into the workflow.
If you find yourself needing a curation workflow that can use external properties
(say to curated a data set with IC50 values for a machine learning/QSAR model)
look into `chemcurry-learn` which extends `chemcurry` to support this.

# Installing
You can install the `chemcurry` package using pip:
```shell
pip install chemcurry
```

This package was built using poetry, so you can also install it by cloning the
repository and create a poetry environment
(though this is not recommend outside of development).
```shell
git clone https://github.com/jimmyjbling/ChemCurry
cd chemcurry
poetry install
```

# Building and running a workflow

Building a chemical curation workflow with `chemcurry` requires only a few lines of code

```python
smiles = ["CCCC", "CCCO", "CCCCN"]

from chemcurry.workflow import CurationWorkflow
from chemcurry.steps import AddH, Add3D, FilterMW, RemoveStereochem

steps = [
    AddH(),
    Add3D(timeout=30),
    FilterMW(max_mw=100, min_mw=10),
    RemoveStereochem()
]

my_workflow = CurationWorkflow(steps=steps)
curated_chemicals = my_workflow.curate_smiles(smiles)
```

The result of the workflow run, a `CuratedChemicalSet` contains all the info about which
compounds failed curation, which compound were altered and why/how all of it happened.
You can save save that info in a human readable report by simply running
```python
curated_chemicals.write_report("path/to/my/report.txt")
```

You can also extract the curated smiles, either as canonical smiles or rdkit Mols
```python
curated_mols = curated_chemicals.to_mols()
curated_smiles = curated_chemicals.to_smiles()
```

### History tracking
You can optionally turn on history tracking mode if you want extremely detailed
information about the evolution of chemical as they progress through curation.
This comes at the expense of extra memory. All you need to do is set `history_tracking=True`
when initializing your workflow. This will save copies of the molecules after each
update is made to them so you can render the full history of the molecule. This
can be done by looping through the `Molecule` objects attached to the curation output
in the `molecules` attribute.

> **Note:** Right now there is not alot you can do with history. In the future, extra features like
> viewing the history of the molecule as an image might be added.

# Saving, loading and sharing workflows
After making and using a workflow, there is a good chance you will want to save it,
either so you can using it again later without having to redefine it, or so you can
share it as part of a publication or project. You can do this by creating a workflow file
(see [here](./docs/workflow_files.rst) for more info on these files)

All you need to do is
```python
my_workflow.save_workflow_file("path/to/my/workflow.json")
```

To load in an existing one you can use
```python
my_workflow = CurationWorkflow.load("path/to/my/workflow.json")
```

Simplate as that. There are some checks and other things happening under the
hood to help prioritize reproducibility and prevent unexpected behavior. You
can read more about how all that work [here](./docs/workflow_files.rst)

# Creating a custom curation steps
The curation functions that already exist in `chemcurry` are unlikely to
always have everything you need.
chemcurry defines very simple APIs that allow you to easily write your own curation steps
You can read more about how that work [here](./docs/define_curation_step.rst)

If you do make your own, we humbly request you submit them to chemcurry so that
the community can benefit from them. Simply make a fork, push your new function
(and its unit test) and then make a pull request. You can read more about
contributing to chemcurry [here]
