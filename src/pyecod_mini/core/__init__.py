"""
Core domain partitioning algorithm

This module contains the proven mini algorithm for domain partitioning,
including evidence parsing, domain boundary detection, and optimization.
"""

# Data models
from .models import (
    Evidence,
    Domain,
    DomainLayout,
    UnassignedSegment,
    PartitionMetadata,
    AlignmentData,
    SegmentType,
    FragmentSize,
)
from .sequence_range import SequenceRange, SequenceSegment, parse_range, positions_to_range

# Core functionality
from .parser import parse_domain_summary, load_reference_lengths, load_protein_lengths
from .partitioner import partition_domains
from .decomposer import load_domain_definitions
from .boundary_optimizer import BoundaryOptimizer
from .writer import write_domain_partition, write_domain_partition_from_layout
from .blast_parser import load_chain_blast_alignments
from .evidence_utils import calculate_evidence_confidence
from .ecod_domains_parser import load_ecod_classifications

__all__ = [
    # Data models
    "Evidence",
    "Domain",
    "DomainLayout",
    "UnassignedSegment",
    "PartitionMetadata",
    "AlignmentData",
    "SegmentType",
    "FragmentSize",
    "SequenceRange",
    "SequenceSegment",
    # Parsers
    "parse_domain_summary",
    "load_reference_lengths",
    "load_protein_lengths",
    "load_domain_definitions",
    "load_chain_blast_alignments",
    "load_ecod_classifications",
    # Core algorithm
    "partition_domains",
    "BoundaryOptimizer",
    # Writers
    "write_domain_partition",
    "write_domain_partition_from_layout",
    # Utilities
    "parse_range",
    "positions_to_range",
    "calculate_evidence_confidence",
]
