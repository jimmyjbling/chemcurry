"""group by curation steps"""

from typing import Union

from ..chemical import Chemical
from .base import GroupBy


class GroupBySmiles(GroupBy):
    """
    Group chemicals by their canonical SMILES"""

    def __init__(self):
        self.note = "grouped chemical into group {}"
        self.rank = 9

    def _get_group_attribute(self, chemical: Chemical) -> Union[int, str, float]:
        chemical.cache_smiles_hash()
        return chemical.smiles_hash_cache


class GroupByMolHash(GroupBy):
    """
    Group chemicals by their RDKit Mol object hash

    This is a more rigours grouping than SMILES, as it takes into account
    3D conformers, advanced stereochem and other properties not encoded in SMILES
    """

    def __init__(self):
        """Initialize GroupByMolHash"""
        self.note = "grouped chemical into group {}"
        self.rank = 9

    def _get_group_attribute(self, chemical: Chemical) -> int:
        chemical.cache_smiles_hash()
        return chemical.mol_hash_cache
