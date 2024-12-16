
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


# Saving, loading and sharing workflows


# Creating a custom curation steps
