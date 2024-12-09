"""base classes and functions for curation steps"""

import abc
from collections import defaultdict
from typing import DefaultDict, List, Union

from ..chemical import BaseChemicalGroup, Chemical
from ..flags import CurationIssue, CurationNote


def check_for_boost_rdkit_error(error_message: str) -> bool:
    """
    Return True if TypeError message is a RDKit Boost Error

    This is because you cannot directly catch RDKit errors, but
    catching all errors could make resolving bugs harder.
    This function should help enable the identification of
    RDKit errors

    Notes
    -----
    This is not perfect, there is a possibility that a non-rdkit error
    can get caught since it just looks for text 'rdkit' and 'boost' in
    the error message. If you have a weird bug you cannot sort out it is
    worth taking a look and seeing if this is catching a error it should
    not be

    Parameters
    ----------
    error_message: str
        the error message as a string

    Returns
    -------
    bool
    """
    error_message = error_message.lower()
    return ("boost" in error_message) and ("rdkit" in error_message)


class CurationStepError(Exception):
    """
    Default exception to throw if there is an error raised with a CurationStep

    This should only be raised if there is an error that is caused by CurationStep itself.
    For example, the curation step lacks both a note and and issue message during instantiation.
    """

    pass


class PostInitMeta(abc.ABCMeta, type):
    """
    Enables a '__post_init__' hook that is called after '__init__'

    To help with error handling when defining new CurationSteps it would be nice
    to enable some auto check that specific class attributes have been declared properly

    But this requires that a "check" function is automatically called after initialization.
    Generic python objects don't have this hook, so we can add it with a MetaClass

    I generally tend to avoid meta classes because they are python black magic,
    but this one is pretty simple and gives us the utility we want
    """

    def __call__(cls, *args, **kwargs):
        """Add post-init hook"""
        instance = super().__call__(*args, **kwargs)
        if post := getattr(cls, "__post_init__", None):
            post(instance, *args, **kwargs)
        return instance


class CurationStep(abc.ABC, metaclass=PostInitMeta):
    """
    The base abstract class for all CurationSteps.

    On a high level, a CurationStep is a callable object that
    wraps/implements some curation function.

    All curation functions can flag chemicals with either an 'issue' or a 'note'.
    Issue flags means the chemical 'failed' that curation
    step and should be removed from the final dataset.
    Note flags mean that the curation step has somehow altered the chemical
    (or its representation).
    The curation step will not remove any compounds, just flag them.

    To avoid attaching flags to objects (which allows for easy compatability with RDKit),
    a CurationStep instead calculates a boolean mask to identify which compounds get a flag.
    A operate mask is also calculated for issues.

    If you want to make your own curationStep, see "Defining Curation Steps" in the docs

    Attributes
    ----------
    issue: CurationIssue
        the associated curation issue to attached to mol that gets flagged
        (None if no flagging occurs)
    note: CurationNote
        the associated CurationNote to attach to a mol that gets changed
        (None if no change occurs)
    dependency: set[str], default=set()
        the set of __name__ attributes for the CurationSteps this
        CurationStep is dependent on
    rank: int
        the rank of this CurationStep (lower is higher priority)
        must be a positive non-zero integer
    """

    issue: CurationIssue
    note: CurationNote
    dependency: set[str] = set()
    rank: int

    def __post_init__(self):
        """
        Called after __init__ finishes for object

        This is primary to check that user defined CurationSteps are
        valid and compatible within the CurationWorkflow
        """
        if not hasattr(self, "rank"):
            raise CurationStepError(
                f"CurationSteps must declare a `self.rank` attribute; "
                f"CurationStep '{self.__class__.__name__}' "
                f"does not have a `self.rank` attribute"
            )
        if not isinstance(self.rank, int) or self.rank < 0:
            raise CurationStepError(
                f"CurationSteps require that the `self.rank` parameter is "
                f"declared as a non-zero positive integer in the `__init__` method; "
                f"CurationStep '{self.__class__.__name__}' invalidates this"
            )

        if hasattr(self, "issue"):
            _issue_is_none = self.issue is None
            if not _issue_is_none and not isinstance(self.issue, CurationIssue):
                if isinstance(self.issue, (list, tuple)):
                    raise CurationStepError(
                        "CurationSteps can only handle a single "
                        "`CurationIssue` instance, not multiple"
                    )
                raise CurationStepError(
                    f"CurationSteps require that the `self.issue`"
                    f" parameter is a CurationIssue enum; "
                    f"not a {type(self.issue)}"
                )
        else:
            _issue_is_none = True

        if hasattr(self, "note"):
            _note_is_none = self.note is None
            if not _note_is_none and not isinstance(self.note, CurationNote):
                if isinstance(self.note, (list, tuple)):
                    raise CurationStepError(
                        "CurationSteps can only handle a single "
                        "`CurationNote` instance, not multiple"
                    )
                raise CurationStepError(
                    f"CurationSteps require that the `self.note`"
                    f" parameter is a CurationNote enum; "
                    f"not a {type(self.note)}"
                )
        else:
            _note_is_none = True

        if _issue_is_none and _note_is_none:
            raise CurationStepError(
                f"CurationSteps require that either `self.issue` or `"
                f"self.note` attributes is not None; "
                f"CurationStep '{self.__class__.__name__}' have both "
                f"`self.issue` and `self.note` declared as None"
            )

    @abc.abstractmethod
    def _func(self, chemicals: Union[Chemical, BaseChemicalGroup]) -> None:
        raise NotImplementedError

    @abc.abstractmethod
    def __call__(self, chemicals: Union[Chemical, BaseChemicalGroup]) -> None:
        """Makes CurationStep callable; calls the `_func` function"""
        raise NotImplementedError

    def __str__(self) -> str:
        """Return the name of the CurationStep class as a str"""
        return str(self.__class__.__name__)

    def __repr__(self) -> str:
        """Return the str representation of the CurationStep class"""
        return self.__str__()

    # TODO add support for '|' (or) based dependencies and '&' (and) based
    #  where you can specify specific sets of dependencies together
    def missing_dependency(self, steps: list) -> set[str]:
        """
        Finds all the missing dependency from a given list of steps for this CurationsStep

        Parameters
        ----------
        steps: list[CurationStep]
            the steps that you want to check for missing dependencies

        Returns
        -------
        missing_dependencies: set[str]
            the set of missing_dependencies (will be empty set is none are missing)

        """
        return self.dependency - set([str(step) for step in steps])


class SingleCurationStep(CurationStep, abc.ABC):
    """
    The base abstract class for all CurationSteps that operate on individual chemicals.

    This means that the curation function is independent of the information present in
    other chemicals
    """

    @abc.abstractmethod
    def _func(self, chemical: Chemical):  # type: ignore[override]
        raise NotImplementedError

    def __call__(self, chemicals: List[Chemical]):  # type: ignore[override]
        """Makes CurationStep callable; calls the `_func` function"""
        for chemical in chemicals:
            self._func(chemical)


class GroupCurationStep(CurationStep, abc.ABC):
    """
    The base abstract class for all CurationSteps that operate on groups of chemicals.

    This means that the curation function is dependent of the information present in
    other chemicals of the group
    """

    @abc.abstractmethod
    def _func(self, chemical_group: BaseChemicalGroup):  # type: ignore[override]
        raise NotImplementedError

    def __call__(self, chemical_groups: List[BaseChemicalGroup]):  # type: ignore[override]
        """Makes CurationStep callable; calls the `_func` function"""
        for chemical_group in chemical_groups:
            self._func(chemical_group)


class GroupBy(abc.ABC):
    """abstract curation function for steps that combine Chemicals into groups"""

    group_class: type

    @abc.abstractmethod
    def _get_group_attribute(self, chemical: Chemical) -> Union[int, str, float]:
        """
        Convert Chemical obj into its group key

        Parameters
        ----------
        chemical: Chemical
            the chemical to get the group key for

        Returns
        -------
        Union[int, str, float]
        """
        raise NotImplementedError

    def __call__(self, chemicals: List[Chemical]) -> List[BaseChemicalGroup]:
        """Call group-by on list of chemical groups"""
        _hash_map: DefaultDict[Union[int, str, float], List[Chemical]] = defaultdict(list)

        for chemical in chemicals:
            _hash_map[self._get_group_attribute(chemical)].append(chemical)

        return [self.group_class(group) for group in _hash_map.values()]


class Aggregate:
    """abstract class for steps that convert chemical groups to single chemicals"""

    @abc.abstractmethod
    def _get_keep_chemical(self, chemical_group: BaseChemicalGroup) -> Chemical:
        """
        take a chemical group and select a single chemical

        Parameters
        ----------
        chemical_group: BaseChemicalGroup
            the chemical group to aggregate

        Returns
        -------
        Chemical
            the single chemical to represent this group
        """
        raise NotImplementedError

    def __call__(self, chemical_group: BaseChemicalGroup) -> Chemical:
        """Aggregate chemical group into a single chemical"""
        return self._get_keep_chemical(chemical_group)
