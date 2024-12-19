"""curation step functions and classes"""

from .add_3d import Add3D
from .add_hydrogen import AddH
from .base import BaseCurationStep, Filter, Update
from .boron import FlagBoron
from .charge import Neutralize
from .inorganic import FlagInorganic
from .mixture import DemixLargestFragment, FlagMixtures
from .mw import FilterMW
from .remove_hydrogen import RemoveAllHs, RemoveHs
from .sanitize import SanitizeMolecule
from .stereochem import RemoveStereochem


def get_step(name) -> BaseCurationStep:
    """Get a curation step by name"""
    try:
        return globals()[name]
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
