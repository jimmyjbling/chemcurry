 - [*] Need to update all the curation functions to handle the new text based
issue and note system
 - [*] remove the legacy flags.py as we swapped to attaching note and issues
directly to curation functions
- [ ] rework ranking system based on the type of function
- [ ] need to implement some type of description text on to each curation-step so
the report can show a detailed summary of what happened
- [*] should make all the calls to issue or notes be "get_issue_text" or "get_note_text" calls
- [ ] implement some the conditional logic into the dependency checking
- [*] make a type alias for label
- [ ] make the docstring better


### Stretch goals
I want to add some functionality for multilabel stuff
I think the best way to do this is to make the chemicals
hold a "LabelDict" which is a mapping of label to val
Should be easy to update the object for this but would break backwards
compatability and need to get bumped up a major version


### SCRAP LABELS
I think this is getting to be too much when trying to handle
labels. I might rework into just a molecule filtering/standardize package
(think something like stoplight but more stuff) which could be help for
nomination or chemical filter of generated compounds and the like

Instead for the ML side of things I think it would be better if the core
class was essentially just a wrapper of pandas/polars dataframes that just implements functions that containerize simple and common commands will also
implementing the same issue and note tracking as in this package.

Could be part of this package but it could also be part of this one, thought
might be too distinct. Could just finish this then make the machine learning one
and have it require this as a dependency
