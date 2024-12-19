## Workflow files
Since ChemCurry places emphasis on reproducibility and sharing of workflows,
a custom file format to store information about workflows was created to help
achieve that goal.

It is not very complicated, its is just a json that records the initialization parameters,
order and name of all the steps in the workflow.
It will additionally hash the workflow object to give it a ID
and record the version number for `rdkit` and `chemcurry` used to build the workflow.
Lastly it will hash all the source code to give the workflow a source code id.

### Source Code ID
The point of the source code ID is to assert that the code between the workflows has not changed,
even if the workflow itself has the same settings/setup.
This is to prevent user that get a workflow from someone else that might have
tweak the code or used a custom fork of ChemCurry from loading in a workflow that
could behave differently from the workflow being loaded by the user.

It would be best practice to never modify ChemCurry source code, or if you create a fork,
to make that fork public and set the `repo_location` for the workflow to be this custom repo
rather than the main chemcurry github so others can access it and update their code if need be.

### Disabling workflow checks
The checking of workflow versions and source code
can be turned off by using the `unsafe=True` parameter in the
`CurationWorkflow.load` function. If this is the case the output workflow files
from any curation will contain an annotation note that the workflow checks were
not run when loading this workflow.
