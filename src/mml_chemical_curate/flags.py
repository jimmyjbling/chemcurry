"""curation flag enums"""

import enum


class CurationNote(enum.StrEnum):
    """
    All possible curation notes that can occur

    Notes are used to notify that a compound has be
    altered by the curation pipeline, BUT NOT removed.
    The text associated with a note is what will be
    rendered in the curation report.

    Whenever a new curation function is added that requires
    a note, it should be taken from the Enum. If it is
    absolutely necessary for a new type of note to be added, it
    should be added here
    """

    failed_curate = "compound failed curation"
    flattened = "compound flattened"
    sanitized = "compound sanitized"
    neutralized = "compound neutralized"
    canonical = "compound made canonical"
    added_hs = "explict H added to compound"
    demixed = "seperated out a mixture component"
    generated_3d_conformer = "3d conformer generated for compound"
    label_made_numeric = "label made numeric"
    digitized_label = "digitized label"
    binarized_label = "binarized label"
    removed_unnecessary_hs = "removed unnecessary H"
    removed_all_Hs = "removed all H"
    filled_missing_label = "filled missing label value"

    @classmethod
    def show_all_notes(cls):
        """Prints all note member names and corresponding note text"""
        for var_name, var_text in cls.__members__.items():
            print(f"{var_name}: {var_text}")


class CurationIssue(enum.StrEnum):
    """
    All possible curation issues that can occur

    issues are used to flag a compound that has
    failed the curation pipeline.

    Whenever a new curation function is added that requires
    a issue, it should be taken from the Enum. If it is
    absolutely necessary for a new type of issue to be added, it
    should be added here
    """

    mixture = "compound is a mixture"
    inorganic = "compound is inorganic"
    boron = "compound has boron"
    rdkit_failed = "rdkit failed to read molecule"
    failed_adding_Hs = "rdkit failed to add H to molecule"
    failed_gen_3d_conformer = "rdkit failed to generated a 3d pose for molecule"
    flatten_failed = "compound failed to be flattened"
    sanitize_failed = "compound failed to be sanitized"
    neutralize_failed = "compound failed to be neutralized"
    canonical_failed = "compound failed to be canonicalized"
    duplicate = "compound is duplicate"
    disagreeing_duplicate = "compound has duplicate with disagreeing label"
    missing_label = "label value is missing"
    non_numeric_label = "label value is not numeric"
    failed_custom_label_filter = "label failed custom label filter"
    wrong_mw = "molecule weight too big or small"
    failed_remove_hs = "rdkit failed to remove H from molecule"

    @classmethod
    def show_all_issues(cls):
        """Prints all issue member names and corresponding issue text"""
        for var_name, var_text in cls.__members__.items():
            print(f"{var_name}: {var_text}")
