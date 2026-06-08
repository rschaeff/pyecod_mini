"""
pyECOD Mini - Clean Domain Partitioning Tool

A minimal, validated domain partitioning tool for ECOD protein classification.

Library API:
    partition_protein() - Main partitioning function (simple, loads refs each call)
    Partitioner         - Batch-optimized partitioner with reference caching
    PartitionResult     - Result dataclass
    PartitionError      - Exception for partition failures
    Domain              - Domain result dataclass

Batch Processing (recommended for multiple proteins):
    from pyecod_mini import Partitioner

    with Partitioner() as p:
        p.load_references_from_config()
        for chain in chains:
            result = p.partition(summary_xml, output_xml, pdb_id, chain_id)
"""

__version__ = "2.1.0"
__author__ = "pyECOD Mini Development Team"

# Export library API
from pyecod_mini.api import (
    Domain,
    PartitionError,
    Partitioner,
    PartitionResult,
    partition_protein,
)
from pyecod_mini.core.exclusions import ExclusionPolicy

__all__ = [
    "partition_protein",
    "Partitioner",
    "PartitionResult",
    "PartitionError",
    "Domain",
    "ExclusionPolicy",
    "__version__",
    "__author__",
]
