- [ ] need to implement some type of description text on to each curation-step so
there is a record of what it is supposed to do
- [ ] add functionality to load and save workflows with version checking
- [ ] add the workflow unit test
- [ ] add the cli
- [ ] add the docs on how to build custom curation functions
- [ ] finish the readme docs

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
