"""Generic label curation functions"""

from typing import Any, Callable

import pandas as pd

from ..flags import CurationIssue, CurationNote
from .base import CurationStep, CurationStepError


class CurateMissingLabel(CurationStep):
    """flag compound with missing labels"""

    def __init__(self):
        self.issue = CurationIssue.missing_label
        self.rank = 0

    def _func(self, molecules):
        for mol in molecules:
            if mol.failed_curation:
                continue
            if pd.isna(mol.label):
                mol.flag_issue(self.issue)


class CurateFillMissingLabel(CurationStep):
    """replace missing labels with a given value"""

    def __init__(self, fill_value: Any):
        self.note = CurationNote.filled_missing_label
        self.rank = 0
        self.fill_value = fill_value

        if pd.isna(self.fill_value):
            raise CurationStepError(
                f"missing fill value cannot also be a missing value: '{self.fill_value}'"
            )

    def _func(self, molecules):
        for mol in molecules:
            if mol.failed_curation:
                continue
            if pd.isna(mol.label):
                mol.update_label(self.fill_value, self.note)
                mol.flag_issue(self.issue)


class CurateNumericalLabel(CurationStep):
    """make labels numeric and flag compounds with non-numeric labels"""

    def __init__(self):
        self.issue = CurationIssue.non_numeric_label
        self.note = CurationNote.label_made_numeric
        self.rank = 5

    def _func(self, molecules):
        for mol in molecules:
            if mol.failed_curation:
                continue
            try:
                mol.update_label(float(mol.label), self.note)
            except ValueError:
                mol.flag_issue(self.issue)


class CurateFilterLabel(CurationStep):
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
        self.issue = CurationIssue.failed_custom_label_filter
        self.rank = 6

    def _func(self, molecules):
        for mol in molecules:
            if mol.failed_curation:
                continue
            if not self.filter_func(mol.label):
                mol.flag_issue(self.issue)


class CurateBinarizeLabel(CurationStep):
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

        self.note = CurationNote.binarized_label
        self.rank = 6
        self.dependency = {"CurateNumericalLabel"}

    def _func(self, molecules):
        for mol in molecules:
            if mol.failed_curation:
                continue
            mol.update_label(
                int(mol.label > self.threshold if self.greater else mol.label < self.threshold),
                self.note,
                force=True,
            )
