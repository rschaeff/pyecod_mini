"""
pyECOD Mini - Clean Domain Partitioning Tool

A minimal, validated domain partitioning tool for ECOD protein classification.

Library API:
    partition_protein() - Main partitioning function
    PartitionResult     - Result dataclass
    PartitionError      - Exception for partition failures
    Domain              - Domain result dataclass
"""

__version__ = "2.0.1"
__author__ = "pyECOD Mini Development Team"

# Export library API
from pyecod_mini.api import Domain, PartitionError, PartitionResult, partition_protein

__all__ = [
    "partition_protein",
    "PartitionResult",
    "PartitionError",
    "Domain",
    "__version__",
    "__author__",
]
