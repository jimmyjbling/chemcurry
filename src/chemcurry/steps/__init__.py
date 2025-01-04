"""curation step functions and classes"""

from .base import BaseCurationStep, Filter, Update

from .chemical.add_3d import Add3D
from .chemical.add_hydrogen import AddH
from .chemical.boron import FlagBoron
from .chemical.charge import Neutralize
from .chemical.inorganic import FlagInorganic
from .chemical.mixture import DemixLargestFragment, FlagMixtures
from .chemical.mw import FilterMW
from .chemical.remove_hydrogen import RemoveAllHs, RemoveHs
from .chemical.sanitize import SanitizeMolecule
from .chemical.stereochem import RemoveStereochem


def get_step(name, *args, **kwargs):
    """Get a curation step by name"""
    try:
        return globals()[name](*args, **kwargs)
    except KeyError as e:
        raise ValueError(f"Unknown curation step: {name}") from e


__all__ = [
    "Filter",
    "Update",
    "BaseCurationStep",
    "Add3D",
    "AddH",
    "FlagBoron",
    "Neutralize",
    "FlagInorganic",
    "FlagMixtures",
    "DemixLargestFragment",
    "FilterMW",
    "RemoveHs",
    "RemoveAllHs",
    "SanitizeMolecule",
    "RemoveStereochem",
    "get_step",
]
