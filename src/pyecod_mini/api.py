#!/usr/bin/env python3
"""
Library API for pyecod_mini - Clean interface for integration with pyecod_prod.

Per PYECOD_MINI_API_SPEC.md.

This module provides a stable API for programmatic access to pyecod_mini's
domain partitioning algorithm, separate from the CLI interface.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional
import xml.etree.ElementTree as ET

import pyecod_mini
from pyecod_mini.core.reference_cache import ReferenceCache, ReferenceData


class PartitionError(Exception):
    """Raised when partitioning fails"""
    pass


@dataclass
class Domain:
    """A single partitioned domain (API result format)"""
    domain_id: str
    range_string: str  # e.g., "10-110" or "10-50,60-110"
    residue_count: int
    source: str  # 'chain_blast', 'domain_blast', 'hhsearch', 'chain_blast_decomposed'
    family_name: str
    confidence: Optional[float] = None


@dataclass
class PartitionResult:
    """Result from domain partitioning (API result format)"""
    success: bool
    pdb_id: str
    chain_id: str
    sequence_length: int
    domains: List[Domain]
    coverage: float  # 0.0-1.0
    partition_xml_path: str
    algorithm_version: str  # e.g., "2.0.0"
    error_message: Optional[str] = None


def partition_protein(
    summary_xml: str,
    output_xml: str,
    pdb_id: str,
    chain_id: str,
    batch_id: Optional[str] = None,
    blast_dir: Optional[str] = None,
) -> PartitionResult:
    """
    Partition a protein into domains using evidence from domain_summary.xml.

    This is the stable library API for pyecod_mini. It wraps the internal
    partitioning logic and provides a clean interface for programmatic use.

    Args:
        summary_xml: Path to domain_summary.xml (input)
        output_xml: Path to partition.xml (output)
        pdb_id: PDB ID
        chain_id: Chain ID
        batch_id: Optional batch ID for tracking
        blast_dir: Optional path to directory containing BLAST XML files
                   (enables chain BLAST decomposition with alignment data)

    Returns:
        PartitionResult with domains, coverage, and metadata

    Raises:
        PartitionError: If partitioning fails
        FileNotFoundError: If summary_xml doesn't exist

    Example:
        >>> result = partition_protein(
        ...     summary_xml="/path/to/8abc_A.summary.xml",
        ...     output_xml="/path/to/8abc_A.partition.xml",
        ...     pdb_id="8abc",
        ...     chain_id="A",
        ...     blast_dir="/path/to/blast",  # Optional: enables chain BLAST
        ... )
        >>> print(f"Found {len(result.domains)} domains, {result.coverage:.1%} coverage")
    """

    # Validate inputs
    summary_path = Path(summary_xml)
    if not summary_path.exists():
        raise FileNotFoundError(f"Summary XML not found: {summary_xml}")

    output_path = Path(output_xml)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Import internal partition logic
    try:
        from pyecod_mini.cli.partition import partition_protein as cli_partition
        from pyecod_mini.cli.config import PyEcodMiniConfig
        from pyecod_mini.core.models import DomainLayout
    except ImportError as e:
        raise PartitionError(f"Failed to import pyecod_mini internals: {e}") from e

    # Create minimal config for library API
    # Note: This uses default paths for reference data
    try:
        config = PyEcodMiniConfig()
    except Exception as e:
        raise PartitionError(f"Failed to initialize pyecod_mini config: {e}") from e

    # Call internal partition logic
    try:
        # Build protein ID
        protein_id = f"{pdb_id}_{chain_id}"

        # Call CLI partition function with custom paths
        domains = cli_partition(
            protein_id=protein_id,
            config=config,
            batch_id=batch_id,
            verbose=False,  # Library API is silent by default
            visualize=False,
            summary_xml=summary_xml,
            output_path=output_xml,
            blast_dir=blast_dir,  # Pass BLAST directory for alignment data
        )

        # Check if partitioning succeeded
        if domains is None:
            raise PartitionError("Partitioning returned None - likely parsing error")

        # Parse the output XML to get complete results
        # This ensures we return exactly what was written to disk
        import xml.etree.ElementTree as ET

        if not output_path.exists():
            raise PartitionError(f"Partition XML was not created: {output_xml}")

        tree = ET.parse(output_xml)
        root = tree.getroot()

        # Extract metadata
        metadata_elem = root.find("metadata")
        version_elem = metadata_elem.find("version") if metadata_elem is not None else None
        algorithm_version = (
            version_elem.get("algorithm") if version_elem is not None
            else pyecod_mini.__version__
        )

        # Extract sequence length and coverage
        stats_elem = metadata_elem.find("statistics") if metadata_elem is not None else None
        if stats_elem is not None:
            sequence_length = int(stats_elem.get("sequence_length", "0"))
            coverage = float(stats_elem.get("total_coverage", "0.0"))
        else:
            # Fallback: estimate from domains
            if domains:
                max_pos = max(d.range.segments[-1].end for d in domains)
                sequence_length = int(max_pos * 1.1)
                total_assigned = sum(d.length for d in domains)
                coverage = total_assigned / sequence_length if sequence_length > 0 else 0.0
            else:
                sequence_length = 0
                coverage = 0.0

        # Convert internal Domain objects to API Domain format
        api_domains = []
        for domain in domains:
            api_domain = Domain(
                domain_id=domain.id,
                range_string=str(domain.range),
                residue_count=domain.length,
                source=domain.source,
                family_name=domain.family,
                confidence=domain.confidence_score,
            )
            api_domains.append(api_domain)

        # Return successful result
        return PartitionResult(
            success=True,
            pdb_id=pdb_id,
            chain_id=chain_id,
            sequence_length=sequence_length,
            domains=api_domains,
            coverage=coverage,
            partition_xml_path=str(output_path),
            algorithm_version=algorithm_version,
            error_message=None,
        )

    except FileNotFoundError as e:
        # Re-raise as-is
        raise

    except Exception as e:
        # Wrap all other errors as PartitionError
        error_msg = f"Partitioning failed: {e}"

        # Try to return partial result if output was created
        if output_path.exists():
            return PartitionResult(
                success=False,
                pdb_id=pdb_id,
                chain_id=chain_id,
                sequence_length=0,
                domains=[],
                coverage=0.0,
                partition_xml_path=str(output_path),
                algorithm_version=pyecod_mini.__version__,
                error_message=error_msg,
            )
        else:
            raise PartitionError(error_msg) from e


class Partitioner:
    """
    Batch-optimized domain partitioner with reference data caching.

    This class loads reference data once and reuses it across multiple
    partition operations, providing ~50x speedup for batch processing.

    Usage:
        # Initialize and load references once
        partitioner = Partitioner()
        partitioner.load_references(
            domain_definitions_file="path/to/domain_definitions.csv",
            reference_lengths_file="path/to/domain_lengths.csv",
            protein_lengths_file="path/to/protein_lengths.csv",
        )

        # Process many proteins efficiently
        for summary_xml, output_xml in batch:
            result = partitioner.partition(
                summary_xml=summary_xml,
                output_xml=output_xml,
                pdb_id=pdb_id,
                chain_id=chain_id,
            )

        # Clean up when done
        partitioner.close()

    Context manager usage:
        with Partitioner() as p:
            p.load_references(...)
            for chain in chains:
                result = p.partition(...)
        # References automatically released

    Performance:
        - Without caching: ~60-120 seconds per chain (reference loading dominates)
        - With caching: ~2-5 seconds per chain (algorithm time only)
        - For 3,677 chains: 92 hours -> ~30 minutes
    """

    def __init__(self, verbose: bool = False):
        """
        Initialize the Partitioner.

        Args:
            verbose: Whether to print progress information
        """
        self._cache = ReferenceCache()
        self._verbose = verbose
        self._partition_count = 0

    def load_references(
        self,
        domain_definitions_file: Optional[str] = None,
        reference_lengths_file: Optional[str] = None,
        protein_lengths_file: Optional[str] = None,
        blacklist_file: Optional[str] = None,
    ) -> "Partitioner":
        """
        Load reference data for partitioning.

        This should be called once before processing multiple proteins.
        The data is cached and reused for all subsequent partition() calls.

        Args:
            domain_definitions_file: Path to domain definitions CSV
            reference_lengths_file: Path to reference lengths CSV
            protein_lengths_file: Path to protein lengths CSV
            blacklist_file: Path to reference blacklist CSV (optional)

        Returns:
            self (for method chaining)

        Example:
            partitioner = Partitioner().load_references(
                domain_definitions_file="domain_definitions.csv",
                reference_lengths_file="domain_lengths.csv",
                protein_lengths_file="protein_lengths.csv",
            )
        """
        self._cache.load(
            domain_definitions_file=domain_definitions_file,
            reference_lengths_file=reference_lengths_file,
            protein_lengths_file=protein_lengths_file,
            blacklist_file=blacklist_file,
            verbose=self._verbose,
        )
        return self

    def load_references_from_config(self) -> "Partitioner":
        """
        Load reference data using default configuration paths.

        This is a convenience method that uses the standard pyecod_mini
        configuration to locate reference files.

        Returns:
            self (for method chaining)
        """
        from pyecod_mini.cli.config import PyEcodMiniConfig

        config = PyEcodMiniConfig()
        return self.load_references(
            domain_definitions_file=str(config.domain_definitions_file),
            reference_lengths_file=str(config.domain_lengths_file),
            protein_lengths_file=str(config.protein_lengths_file),
            blacklist_file=(
                str(config.reference_blacklist_file)
                if config.reference_blacklist_file.exists()
                else None
            ),
        )

    def is_loaded(self) -> bool:
        """Check if reference data has been loaded."""
        return self._cache.is_loaded()

    def partition(
        self,
        summary_xml: str,
        output_xml: str,
        pdb_id: str,
        chain_id: str,
        batch_id: Optional[str] = None,
        blast_dir: Optional[str] = None,
    ) -> PartitionResult:
        """
        Partition a protein into domains using cached reference data.

        This method uses pre-loaded reference data for efficient batch processing.
        Call load_references() or load_references_from_config() before using this.

        Args:
            summary_xml: Path to domain_summary.xml (input)
            output_xml: Path to partition.xml (output)
            pdb_id: PDB ID
            chain_id: Chain ID
            batch_id: Optional batch ID for tracking
            blast_dir: Optional path to directory containing BLAST XML files

        Returns:
            PartitionResult with domains, coverage, and metadata

        Raises:
            PartitionError: If partitioning fails
            RuntimeError: If reference data not loaded
            FileNotFoundError: If summary_xml doesn't exist
        """
        if not self._cache.is_loaded():
            raise RuntimeError(
                "Reference data not loaded. Call load_references() or "
                "load_references_from_config() before partition()."
            )

        # Validate inputs
        summary_path = Path(summary_xml)
        if not summary_path.exists():
            raise FileNotFoundError(f"Summary XML not found: {summary_xml}")

        output_path = Path(output_xml)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # Import internal modules
            from pyecod_mini.core.blast_parser import load_chain_blast_alignments
            from pyecod_mini.core.boundary_optimizer import BoundaryOptimizer
            from pyecod_mini.core.models import DomainLayout
            from pyecod_mini.core.parser import parse_domain_summary
            from pyecod_mini.core.partitioner import partition_domains
            from pyecod_mini.core.writer import (
                create_metadata_from_batch,
                write_domain_partition,
                write_domain_partition_from_layout,
            )

            # Get cached reference data
            ref_data = self._cache.get_data()
            protein_id = f"{pdb_id}_{chain_id}"

            # Determine BLAST directory
            if blast_dir:
                blast_path = Path(blast_dir)
            else:
                blast_path = summary_path.parent.parent / "blast"

            # Load BLAST alignments (not cached - small per protein)
            blast_alignments = {}
            if blast_path.exists():
                blast_alignments = load_chain_blast_alignments(
                    str(blast_path), pdb_id, chain_id, verbose=self._verbose
                )

            # Parse evidence using cached reference data
            evidence = parse_domain_summary(
                str(summary_path),
                reference_lengths=ref_data.reference_lengths,
                protein_lengths=ref_data.protein_lengths,
                blast_alignments=blast_alignments,
                require_reference_lengths=True,
                verbose=self._verbose,
            )

            # Read sequence length from summary XML
            tree = ET.parse(str(summary_path))
            root = tree.getroot()
            protein_elem = root.find("protein")

            if protein_elem is not None and protein_elem.get("length"):
                sequence_length = int(protein_elem.get("length"))
            elif evidence:
                max_pos = max(ev.query_range.segments[-1].end for ev in evidence)
                sequence_length = int(max_pos * 1.1)
            else:
                sequence_length = 0

            # Determine batch info
            batch_dir = summary_path.parent.parent
            batch_name = batch_id or batch_dir.name

            # Handle no evidence case
            if not evidence:
                metadata = create_metadata_from_batch(
                    pdb_id, chain_id, str(batch_dir), batch_name
                )
                metadata.sequence_length = sequence_length
                metadata.process_parameters.update({
                    "evidence_items_found": 0,
                    "domains_assigned": 0,
                    "boundary_optimization_enabled": False,
                    "cached_references_used": True,
                })
                write_domain_partition([], metadata, str(output_path))

                self._partition_count += 1
                return PartitionResult(
                    success=True,
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    sequence_length=sequence_length,
                    domains=[],
                    coverage=0.0,
                    partition_xml_path=str(output_path),
                    algorithm_version=pyecod_mini.__version__,
                    error_message=None,
                )

            # Partition domains using cached domain definitions
            domains = partition_domains(
                evidence,
                sequence_length=sequence_length,
                domain_definitions=ref_data.domain_definitions or None,
                verbose=self._verbose,
            )

            # Handle no domains after filtering
            if not domains:
                metadata = create_metadata_from_batch(
                    pdb_id, chain_id, str(batch_dir), batch_name
                )
                metadata.sequence_length = sequence_length
                metadata.process_parameters.update({
                    "evidence_items_processed": len(evidence),
                    "domains_assigned": 0,
                    "boundary_optimization_enabled": False,
                    "quality_filtering_rejected_all_evidence": True,
                    "cached_references_used": True,
                })
                write_domain_partition([], metadata, str(output_path))

                self._partition_count += 1
                return PartitionResult(
                    success=True,
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    sequence_length=sequence_length,
                    domains=[],
                    coverage=0.0,
                    partition_xml_path=str(output_path),
                    algorithm_version=pyecod_mini.__version__,
                    error_message=None,
                )

            # Apply boundary optimization
            layout = DomainLayout.from_domains(domains, sequence_length)
            optimizer = BoundaryOptimizer()
            optimized_layout = optimizer.optimize_boundaries(
                layout, min_domain_size=25, neighbor_tolerance=5, verbose=self._verbose
            )

            final_domains = optimized_layout.domains
            final_stats = optimized_layout.get_coverage_stats()

            # Create metadata
            metadata = create_metadata_from_batch(
                pdb_id, chain_id, str(batch_dir), batch_name
            )
            metadata.sequence_length = sequence_length
            metadata.process_parameters.update({
                "evidence_items_processed": len(evidence),
                "blast_alignments_loaded": len(blast_alignments),
                "domain_definitions_available": ref_data.domain_definitions_count,
                "reference_lengths_available": ref_data.reference_lengths_count,
                "boundary_optimization_enabled": True,
                "min_domain_size": 25,
                "neighbor_tolerance": 5,
                "domains_before_optimization": len(domains),
                "domains_after_optimization": len(final_domains),
                "cached_references_used": True,
            })

            # Write output
            write_domain_partition_from_layout(
                layout=optimized_layout, metadata=metadata, output_path=str(output_path)
            )

            # Convert to API domain format
            api_domains = []
            for domain in final_domains:
                api_domain = Domain(
                    domain_id=domain.id,
                    range_string=str(domain.range),
                    residue_count=domain.length,
                    source=domain.source,
                    family_name=domain.family,
                    confidence=domain.confidence_score,
                )
                api_domains.append(api_domain)

            self._partition_count += 1

            return PartitionResult(
                success=True,
                pdb_id=pdb_id,
                chain_id=chain_id,
                sequence_length=sequence_length,
                domains=api_domains,
                coverage=final_stats["coverage_percent"] / 100.0,
                partition_xml_path=str(output_path),
                algorithm_version=pyecod_mini.__version__,
                error_message=None,
            )

        except FileNotFoundError:
            raise
        except Exception as e:
            error_msg = f"Partitioning failed: {e}"
            if output_path.exists():
                return PartitionResult(
                    success=False,
                    pdb_id=pdb_id,
                    chain_id=chain_id,
                    sequence_length=0,
                    domains=[],
                    coverage=0.0,
                    partition_xml_path=str(output_path),
                    algorithm_version=pyecod_mini.__version__,
                    error_message=error_msg,
                )
            raise PartitionError(error_msg) from e

    def close(self) -> None:
        """
        Release cached reference data to free memory.

        Call this when done processing a batch to release memory.
        The Partitioner can be reused by calling load_references() again.
        """
        self._cache.clear()

    @property
    def partition_count(self) -> int:
        """Number of proteins partitioned since loading references."""
        return self._partition_count

    def summary(self) -> str:
        """Return a summary of the partitioner state."""
        cache_summary = self._cache.summary()
        return f"Partitioner: {self._partition_count} partitions completed. {cache_summary}"

    def __enter__(self) -> "Partitioner":
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Context manager exit - release cached data."""
        self.close()

    def __repr__(self) -> str:
        status = "loaded" if self._cache.is_loaded() else "not loaded"
        return f"Partitioner({status}, {self._partition_count} partitions)"
