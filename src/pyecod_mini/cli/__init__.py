"""
Command-line interface for pyecod_mini

Provides the pyecod-mini command for domain partitioning.
"""

from .config import BatchFinder, PyEcodMiniConfig
from .main import main
from .partition import analyze_protein_batches, partition_protein
from .utils import run_test_suite, setup_references

__all__ = [
    "main",
    "PyEcodMiniConfig",
    "BatchFinder",
    "partition_protein",
    "analyze_protein_batches",
    "setup_references",
    "run_test_suite",
]
