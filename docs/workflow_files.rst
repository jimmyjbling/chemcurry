Workflow files
==============

Workflows are saved as as json files, but with some extra keys to
help with asserting reproducibility when sharing workflows.

Workflow file keys
------------------
There are several keys in the workflow file format:
 - workflow_name: this is the workflow name (if it is set)
 - workflow_description: this is the workflow description
 - workflow_hash: a hash of the workflow object
 - chemcurry_repo_url: the url (or uri) to the chemcurry repo
   used to build the workflow
 - workflow_source_code_hash: the hash of the source code for the
   workflow object
 - versions: a dict of package name and versions numbers for dependent packages
 - steps: a ordered dict of the step
 - num_steps: the number of steps in the workflow

You should never make or alter the workflow files after they are created.
While human readable, these files are meant for saving, loading and sharing
workflows, not creating them. You should create workflows using the ``chemcurry``
package instead.

Safety checks
-------------

To help with keeping workflows reproducible and prevent unexpected behavior
there are a several checks that are made to see if the workflow will behave
the same way it did when it was first created. The primary way this is done
is through checking the hash of the workflow object itself *and* the hash of
the source code from both the workflow and all its steps.
These hashes all must match for a workflow to be safely loaded.
Further checks are made to make sure version numbers between chemcurry and rdkit
match the workflows required versions.

If you are only ever using the main chemcurry repo, it should rare that
there is a source code conflict (unless you altered the site-package source
yourself). You can generally fix any by making sure the version numbers match.
If you got the workflow from someone else, reach out and see if they are using
a fork or otherwise separate version of chemcurry not attached to the main repo.

If you must run a workflow and cannot figure out the correct settings to make it
load safely, you can load the workflow 'unsafely' instead by using
``CurationWorkflow.load(path, safe=False)``. However, this will flag the workflow
as unsafe, and any curation reports generated will be flagged as unsafe.
