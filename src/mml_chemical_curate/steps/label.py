"""Generic label curation functions"""

from typing import Any, Callable

import pandas as pd

from .base import CurationStepError, SingleCurationStep


class CurateMissingLabel(SingleCurationStep):
    """flag compound with missing labels"""

    def __init__(self):
        self.issue = "chemical has missing label"
        self.rank = 0

    def _func(self, chemical):
        if pd.isna(chemical.label):
            chemical.flag_issue(self.get_issue_text())


class CurateFillMissingLabel(SingleCurationStep):
    """replace missing labels with a given value"""

    def __init__(self, fill_value: Any):
        self.rank = 0
        self.fill_value = fill_value

        if pd.isna(self.fill_value):
            raise CurationStepError(
                f"missing fill value cannot also be a missing value: '{self.fill_value}'"
            )
        self.note = f"filled missing label with {fill_value}"

    def _func(self, chemical):
        if pd.isna(chemical.label):
            chemical.update_label(self.fill_value, self.get_note_text())
            chemical.flag_issue(self.get_issue_text())


class CurateNumericalLabel(SingleCurationStep):
    """make labels numeric and flag compounds with non-numeric labels"""

    def __init__(self):
        self.issue = "chemical label could not be cast to numeric"
        self.note = "chemical label made numeric"
        self.rank = 5

        self.dependency = {"CurateMissingLabel|CurateFillMissingLabel"}

    def _func(self, chemical):
        try:
            chemical.update_label(float(chemical.label), self.get_note_text())
        except ValueError:
            chemical.flag_issue(self.get_issue_text())


class CurateFilterLabel(SingleCurationStep):
    """filter compound based on labels and some custom filter function"""

    def __init__(self, filter_func: Callable):
        """
        Initialize curation step

        Parameters
        ----------
        filter_func: Callable
            should be a function that takes in a single float, str, int or None value
            and return True or False based on some static rule
        """
        self.filter_func = filter_func
        self.issue = "chemical label filtered out by custom filter"
        self.rank = 6

        self.dependency = {"CurateMissingLabel|CurateFillMissingLabel"}

    def _func(self, chemical):
        if not self.filter_func(chemical.label):
            chemical.flag_issue(self.get_issue_text())


class CurateBinarizeLabel(SingleCurationStep):
    """convert numerical labels in binary (0, 1) class labels"""

    def __init__(self, threshold: float, greater: bool = True):
        """
        Initialize curation step

        Notes
        -----
        All comparisons are exclusion (i.e. < not <=)

        Parameters
        ----------
        threshold: float
            threshold to use for binary cutoff
        greater: bool, default=True
            if True, greater will be class '1' and less than class '0'
            if False, greater will be class '0' and less than class '1'
        """
        super().__init__()
        self.threshold = threshold
        self.greater = greater

        self.note = f"label was barbarized with threshold {threshold}"
        self.rank = 6
        self.dependency = {"CurateNumericalLabel"}

    def _func(self, chemical):
        chemical.update_label(
            int(
                chemical.label > self.threshold
                if self.greater
                else chemical.label < self.threshold
            ),
            self.note,
            force=True,
        )
