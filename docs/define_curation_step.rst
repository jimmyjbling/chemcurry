Defining Curation Steps
========================

It is very likely that at some point a new curation function will be needed
that does not currently exist in the package. In this case you will need to
add it yourself. This doc outlines how to go about the process of create new
curation step classes to add to your curation workflows

There are two types of curation steps:

 * ``MoleculeCurationStep``
 * ``LabelCurationStep``

These both implement the base ``CurationStep`` class.
A ``MoleculeCurationStep`` is used when the curation only takes into
account the molecules themselves (not label/external information).
``LabelCurationStep`` means that lables are taken into account (and might be altered).

The main entrypoint for curation is the ``_func`` function associated with the class.
For a ``CurationStep`` this function should always have the following signature:
_func(self, molecules, **kwargs) -> Tuple(issue_mask: NDArray, note_mask: NDarray)

Calling a ``CurationStep`` should always return:
    - 1: a boolean mask for which chemicals had issues
    - 2: a boolean mask for which chemicals had notes

CurationSteps must implement 2 methods:

    1. _func(self, molecules, **kwargs) -> issue_mask: NDArray, note_mask: NDarray
        - This function is what will be called for the main curation
          It should return a boolean mask defining which inputs pass (True)
          or failed (False) the curation step
          All CurationSteps operate on rdkit.Mol objects
          so molecules should be some type of indexable iterable of rdkit.Mol objects
          It should be able to handle None objects in the molecules list
          The shape of the mask return must match the shape of molecules (input)
          Its important that _func handles **kwargs, and this can be used to
          add any additional settings to the functions
    2. get_rank() -> int
        - This function is used to determine when this CurationStep should occur with respect to other steps.
          Rank 0 is reserved for rdkit molecule generation
          The order for ranking is as follows:
            1. Separate a multi-component chemical (i.e. mixture handling)
            2. Changes/Standardizes a molecule (e.g. neutralize)
            3. Excludes a molecule (i.e. remove inorganic)
            4. Generate 3D information
            5. Standardize labels (i.e. make numeric/categorical)
            6. Modify labels (i.e. binarize continuous labels)
            7. Duplicate handling
          excludes means the function will not change anything, just flag
          changes means inputs are both flagged and possibly changed

You must also set the corresponding CurationIssue and CurationNote for the function. These will be attached to any
mol that failed the curation step (issue) or is somehow altered by the CurationStep (note). One of the two (or both)
must be set. Each function can only on one of each though. If you find yourself needing two, you probably have two
curation functions mixed into one

Label Curation
------------------------

A ``LabelCurationStep`` is used when the curation takes into
account the molecules `and` some type of label associated with it

The only difference between this class of curation function and the
``MoleculeCurationStep`` is that it requires an additional ``label``
parameter in the ``_func`` signature (this is why **kwargs is required).
This is what will get used to curate any of the labels associated with the molecules

Just like the molecules, the labels should be changed inplace.
When writing the functions make sure this is the case,
otherwise any curation functions that alter labels will not work

Dependency
----------
Lastly, some CurationSteps require an explict dependency before it can be used. A notable example is
``CurateRemoveDisagreeingDuplicatesMinMax``, which requires that labels are numeric. Thus, it requires that the
CurationStep ``CurateStandardizeLabels`` is run before it. To help a user understand these dependencies all
CurationSteps have a ``dependency`` attribute that includes the ``__name__`` attribute of any CurationStep class
that it is dependent on.

Some curation steps might need extra arguments (for example, generating a binary threshold would require a threshold
 parameter is passed). In this case, you can pass that as kwarg during construction, and it will be saved in a kwarg
 dict for later use by _func. This dict will automatically be pass as kwargs to _func, so when declaring your custom
 _func, you can add it to the signature. You can also add it to the constructor

