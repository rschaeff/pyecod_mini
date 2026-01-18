#!/usr/bin/env python3
"""
Reference data caching for batch processing optimization.

This module provides caching of reference data (domain definitions, reference lengths,
protein lengths) to avoid redundant I/O when processing multiple proteins in a batch.

Loading reference data (~50MB) takes 60-120 seconds per call. By caching this data,
batch processing performance improves from ~90 seconds/chain to ~2-5 seconds/chain.
"""

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .decomposer import DomainReference, load_domain_definitions
from .parser import load_protein_lengths, load_reference_lengths


@dataclass
class ReferenceData:
    """Container for cached reference data."""

    domain_definitions: dict[tuple[str, str], list[DomainReference]] = field(
        default_factory=dict
    )
    reference_lengths: dict[str, int] = field(default_factory=dict)
    protein_lengths: dict[tuple[str, str], int] = field(default_factory=dict)

    # Loading metadata
    load_time_seconds: float = 0.0
    domain_definitions_count: int = 0
    reference_lengths_count: int = 0
    protein_lengths_count: int = 0

    def is_loaded(self) -> bool:
        """Check if reference data has been loaded."""
        return (
            bool(self.domain_definitions)
            or bool(self.reference_lengths)
            or bool(self.protein_lengths)
        )

    def summary(self) -> str:
        """Return a summary string of loaded data."""
        return (
            f"ReferenceData: {self.domain_definitions_count} domain definitions, "
            f"{self.reference_lengths_count} reference lengths, "
            f"{self.protein_lengths_count} protein lengths "
            f"(loaded in {self.load_time_seconds:.1f}s)"
        )


class ReferenceCache:
    """
    Cache for reference data used in domain partitioning.

    This class loads reference data once and reuses it across multiple
    partition operations, significantly improving batch processing performance.

    Usage:
        # Create cache and load data
        cache = ReferenceCache()
        cache.load(
            domain_definitions_file="path/to/domain_definitions.csv",
            reference_lengths_file="path/to/domain_lengths.csv",
            protein_lengths_file="path/to/protein_lengths.csv",
        )

        # Use cached data for multiple proteins
        for protein in proteins:
            result = partition_with_cache(protein, cache)

        # Clear cache when done (optional - helps with memory)
        cache.clear()

    Context manager usage:
        with ReferenceCache() as cache:
            cache.load(...)
            for protein in proteins:
                result = partition_with_cache(protein, cache)
        # Cache automatically cleared on exit
    """

    def __init__(self):
        self._data: Optional[ReferenceData] = None
        self._source_files: dict[str, str] = {}

    def load(
        self,
        domain_definitions_file: Optional[str] = None,
        reference_lengths_file: Optional[str] = None,
        protein_lengths_file: Optional[str] = None,
        blacklist_file: Optional[str] = None,
        verbose: bool = False,
    ) -> ReferenceData:
        """
        Load reference data from files.

        Args:
            domain_definitions_file: Path to domain definitions CSV
            reference_lengths_file: Path to reference lengths CSV
            protein_lengths_file: Path to protein lengths CSV
            blacklist_file: Path to reference blacklist CSV (optional)
            verbose: Print loading progress

        Returns:
            ReferenceData containing all loaded data
        """
        start_time = time.time()

        if verbose:
            print("Loading reference data (this may take a minute)...")

        self._data = ReferenceData()

        # Load domain definitions
        if domain_definitions_file and Path(domain_definitions_file).exists():
            if verbose:
                print(f"  Loading domain definitions from {domain_definitions_file}...")
            self._data.domain_definitions = load_domain_definitions(
                domain_definitions_file,
                verbose=verbose,
                blacklist_path=blacklist_file,
            )
            self._data.domain_definitions_count = len(self._data.domain_definitions)
            self._source_files["domain_definitions"] = domain_definitions_file
            if verbose:
                print(
                    f"    Loaded {self._data.domain_definitions_count} chain definitions"
                )

        # Load reference lengths
        if reference_lengths_file and Path(reference_lengths_file).exists():
            if verbose:
                print(f"  Loading reference lengths from {reference_lengths_file}...")
            self._data.reference_lengths = load_reference_lengths(reference_lengths_file)
            self._data.reference_lengths_count = len(self._data.reference_lengths)
            self._source_files["reference_lengths"] = reference_lengths_file
            if verbose:
                print(f"    Loaded {self._data.reference_lengths_count} reference lengths")

        # Load protein lengths
        if protein_lengths_file and Path(protein_lengths_file).exists():
            if verbose:
                print(f"  Loading protein lengths from {protein_lengths_file}...")
            self._data.protein_lengths = load_protein_lengths(protein_lengths_file)
            self._data.protein_lengths_count = len(self._data.protein_lengths)
            self._source_files["protein_lengths"] = protein_lengths_file
            if verbose:
                print(f"    Loaded {self._data.protein_lengths_count} protein lengths")

        self._data.load_time_seconds = time.time() - start_time

        if verbose:
            print(f"  Reference data loaded in {self._data.load_time_seconds:.1f}s")

        return self._data

    def is_loaded(self) -> bool:
        """Check if reference data has been loaded."""
        return self._data is not None and self._data.is_loaded()

    def get_data(self) -> ReferenceData:
        """Get the cached reference data."""
        if self._data is None:
            raise RuntimeError(
                "Reference data not loaded. Call load() first or use load_from_config()."
            )
        return self._data

    @property
    def domain_definitions(self) -> dict[tuple[str, str], list[DomainReference]]:
        """Get cached domain definitions."""
        return self.get_data().domain_definitions

    @property
    def reference_lengths(self) -> dict[str, int]:
        """Get cached reference lengths."""
        return self.get_data().reference_lengths

    @property
    def protein_lengths(self) -> dict[tuple[str, str], int]:
        """Get cached protein lengths."""
        return self.get_data().protein_lengths

    def clear(self) -> None:
        """Clear all cached data to free memory."""
        self._data = None
        self._source_files.clear()

    def summary(self) -> str:
        """Return a summary of cached data."""
        if self._data is None:
            return "ReferenceCache: No data loaded"
        return self._data.summary()

    def __enter__(self) -> "ReferenceCache":
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Context manager exit - clear cache."""
        self.clear()

    def __repr__(self) -> str:
        if self._data is None:
            return "ReferenceCache(not loaded)"
        return f"ReferenceCache({self._data.domain_definitions_count} defs, {self._data.reference_lengths_count} refs)"
