Defining Curation Steps
========================

It is very likely that at some point you will need a new curation function will
be needed that does not currently exist in the package. In this case you will
need to add it yourself. This doc outlines how to go about the process of create
new curation step classes to add to your curation workflows.

Curation Step Types
-------------------

All curation steps are children of ``BaseCurationStep``.
If you are just adding new steps, you will not have to deal with this class.
Instead ``chemcurry`` provides two abstract types of curation steps your new one
can inherit from:
 * ``Update``
 * ``Filter``

Determining which to use is pretty easy.
If your curation step modifies chemicals at any point, it is an ``Update`` step.
If it just 'flags' chemicals as failing/passing curation it is a ``filter`` step.
There are some slight different when defining one versus the other,
but the main difference is that ``Filter`` steps cannot assign notes
(which means it cannot change the chemicals).
Naming them is just as simple; name them directly after what they do.
For example, ``FlagExplicitHydrogen`` would flag any chemicals that have
explict hydrogen atoms.

Writing ``Filter`` steps
^^^^^^^^^^^^^^^^^^^^^^^^
Custom filter steps inherit from the ``Filter`` class.
When it comes to creating the class,
all you *must* do is define one attributes and one function.
First you must define the ``issue`` attribute, which is the text
you want to render alongside any chemical that fails curation due to this step.
It should reflect the issue the Filter step is meant to filter.
For example, "chemical has explict hydrogen atoms".
Then you need to define the ``filter`` function.
This function has a specific signature:
it takes in a single rdkit Mol object and returns a boolean.
At no point in this function alter the passed molecule.
To help prevent this, the a deepcopy of the mol is passed to this
function, so you don't have to worry about accidentally altering the molecule.

Here is an example of a curation step that will flag any chemicals with oxygen atoms:
::
    from chemcurry.steps import Filter
    from rdkit import Chem

    class FlagOxygen(Filter):
        def __init__(self):
            self.issue = "molecule has oxygen"

        def filter(mol) -> bool:
            return mol.HasSubstructureMatch(Chem.MolFromSmarts("[#8]"))

Writing ``Update`` steps
^^^^^^^^^^^^^^^^^^^^^^^^
Custom update steps inherit from the ``Update`` class.
Unlike ``Filter`` steps, which cannot have a must have an ``issue``
but not a ``note``, ``Update`` steps require an ``note`` and can
*optionally* have a ``issue`` as well.
This is because ``Update`` steps can update molecules (and attached notes
when they do) but also flag chemicals with issue if they cannot make the
required update.
For example, ``Add3D`` has a ``note`` for "added a 3d conformer" and an
issue for when "adding 3d conformer timed out".
Not all ``Update`` steps require issues, as some might not be expected to
fail to make the update. An issue should only be defined if you expect a
the update to fail for a specific reason.
The only other thing that must be defined is the ``update`` function.
Like with ``filter``, this must take in a rdkit Mol object.
However, instead of returning a boolean, it must return either a ``None`` or
an rdkit Mol. The mol passed to this function will be a deepcopy, so no need
to copy it yourself in the function; you cannot modify in place you must return it
The function should only return ``None`` if the update step failed. This is what
tells ``chemcurry`` to flag this chemical with an issue during the ``Update`` step.

Here is an example of a curation step that will convert oxygen atoms into dummy atoms
but only when it is Sunday, otherwise it fails (because why not):
::
    from chemcurry.steps import Update
    from rdkit import Chem
    from datetime import datetime

    class MakeOxygenDummyOnSunday(Update):
        def __init__(self):
            self.issue = "its not sunday, cannot change oxygen atoms to dummy atoms"
            self.note = "oxygen atoms changed to dummy atoms"

        def update(mol):
            editable_mol = Chem.RWMol(mol)

            for atom in editable_mol.GetAtoms():
                if atom.GetAtomicNum() == 8:
                atom.SetAtomicNum(0)
                atom.SetSymbol("*")

            if datetime.note().weekday() == 6:
                return editable_mol.GetMol()
            else
                return None

And that's it. This is all you will need to add new functions.
In fact, most of the curation functions are only 10-20 lines of code


Adding parameters to curation steps
-----------------------------------

Some curation steps might require parameters to initialize them.
For example, the ``Add3D`` step has a ``timeout`` parameter that
sets the timelimit until the ``Add3D`` step will fail.
If your custom curation step needs these, they should be defined
in the ``__init__`` call for your curation function.
They should have a default value if it makes sense, but they can be
required. Neither ``Filter`` nor ``Update`` define an ``__init__``
function, so you do not need to make a call to super (but you still
can if you want).

For example, here is a filter curation function that remove chemicals that have
more than ``N`` atoms::
    from chemcurry.steps import Filter
    from rdkit import Chem

    class FlagOxygen(Filter):
        def __init__(self, max_atoms: int = 10):
            self.max_atoms = max_atoms
            self.issue = f"molecule has more than {max_atoms} atoms"

        def filter(mol) -> bool:
            return mol.GetNumAtoms() <= self.max_atoms

You don't need to assign the parameter to a class attribute of the same name.
``chemcurry`` will remember the values passed to the initialization call and
save them. That way the function can be save and re-initialized the same way
at a later time (this is in the ``_init_params`` private attribute).


Optional description
--------------------
All curation steps define a ``description`` property that can contain
information about what the curation function does.
It defaults to "NA" if not set when defining the function.
This should be just a short 1-2 sentence description of the functions
purpose/use-case.

.. warning ::
    This may be removed in the future as it includes info that should be
    in the docstring of the curation function.


Dependency
----------
You might have noticed that there is support for specifying
dependencies in ``chemcurry``. Right now, these are not enforced
nor will they raise errors if dependencies are missing in a workflow.
Instead it is just a nice annotation to make when defining functions if
there is an expectation of another function being ran first.
