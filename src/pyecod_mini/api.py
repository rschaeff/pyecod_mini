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

import pyecod_mini


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
