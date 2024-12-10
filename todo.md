 - [*] Need to update all the curation functions to handle the new text based
issue and note system
 - [*] remove the legacy flags.py as we swapped to attaching note and issues
directly to curation functions
- [ ] rework the ranking system it is all broken now that we have group by and
agg functions (might have to forgo it entirely, might remove it and add it back
later as a stretch goal)
- [ ] need to implement some type of description text on to each curation-step so
the report can show a detailed summary of what happened
- [*] should make all the calls to issue or notes be "get_issue_text" or "get_note_text" calls
- [ ] implement some the conditional logic into the dependency checking
