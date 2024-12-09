# MML Chemical Curation Pipeline
This package is meant to both streamline chemical curation of datasets *and*
produce detailed reports about which chemicals were removed by which curation
steps. The Molecular Modeling Lab often finds itself needed to generate these
reports, so this package was developed as a way to standardize our curation
pipelines among members and generated the detailed reports needed for any future
publication and open sharing of our work.

Most curation functions were primarily designed for Quantitative Structure
Activity Relationship (QSAR) modeling. However, we made it is easy to add your own
custom functions that can be more focus for other tasks (and if you do make them,
feel free to contribute them)

Internally, the MML names all their software after food (and their computers after frogs :frog:).
In house, we (or maybe just me) call this package CUPCAKE:
**C**hemical **U**nification and **P**rocessing for **C**uration,
**A**nalysis, and **K**nowledge **E**ngineering
(which also follows our longstanding tradition of forced acronyms).

# Installing
You can install the `mml-chemical-curation` package using pip:
```shell
pip install mml-chemical-curation
```

This package was built using poetry, so you can also install it by cloning the
repository and create a poetry environment
(though this is not recommend outside of development).
```shell
git clone [TODO]
cd mml-chemical-curation
poetry create
poetry install
```

# Creating a Curation Workflow
Using `mml-chemical-curate` to make a curation workflow is simple.
You just need to create a `CurationWorkflow` object and give it the
list of curation steps that you want to include
```python
from mml_chemical_curate import CurationWorkflow

```
