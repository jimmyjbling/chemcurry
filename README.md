
![# ChemCurry](ChemCurry.svg)

# ChemCurry

`chemcurry` is a chemical curation workflow package meant to both streamline
building curation workflows *and* producing detailed reports about which chemicals
where flagged and when. The Molecular Modeling Lab @ UNC often finds itself needed
to generate these reports to show to our PI and share our workflows with new members,
so this package was developed as a way to standardize that process.

While most of the chemical curation workflows can be built in under 100 lines of code,
the core idea behind `chemcurry` is to assert reproducibility and easy building/sharing
among chemist will any level of coding background. Most cheminformatics project and
publications will need to do some type of curation, and, frankly, the methods on how
this is done is often not up to par to make it easy to reproduce. `chemcurry` aims to
fix that.

Closely related to the philosophy of reproducibility, `chemcurry` was also designed to be
easy to add new curation functions too. There is a simple API that really only requires
you to write the same code you might if you were doing in manually in a notebook.

`chemcurry` is design to operate on explict chemical properties, meaning if the property cannot
be calculated by using just the chemical, it will not fit into the workflow.
If you find yourself needing a curation workflow that can use external properties
(say to curated a data set for a machine learning or QSAR model) look into `chemcurry-learn`


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
