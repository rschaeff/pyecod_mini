#!/usr/bin/env python3
"""
Domain partitioning functions for pyecod_mini CLI

Main processing logic for partitioning proteins with provenance tracking.
"""

from typing import Optional

from pyecod_mini.core.blast_parser import load_chain_blast_alignments
from pyecod_mini.core.boundary_optimizer import BoundaryOptimizer
from pyecod_mini.core.decomposer import load_domain_definitions
from pyecod_mini.core.models import DomainLayout
from pyecod_mini.core.parser import (
    load_protein_lengths,
    load_reference_lengths,
    parse_domain_summary,
)
from pyecod_mini.core.partitioner import partition_domains
from pyecod_mini.core.writer import (
    create_metadata_from_batch,
    write_domain_partition,
    write_domain_partition_from_layout,
)

from .config import PyEcodMiniConfig


def partition_protein(
    protein_id: str,
    config: PyEcodMiniConfig,
    batch_id: Optional[str] = None,
    verbose: bool = False,
    visualize: bool = False,
    summary_xml: Optional[str] = None,
    output_path: Optional[str] = None,
    blast_dir: Optional[str] = None,
) -> Optional[list]:
    """Partition domains for a single protein with enhanced provenance tracking

    Args:
        protein_id: Protein ID to process (e.g., '8ovp_A')
        config: Configuration object
        batch_id: Optional batch ID for batch detection
        verbose: Show detailed output
        visualize: Generate PyMOL visualization
        summary_xml: Optional path to custom summary XML (overrides batch detection)
        output_path: Optional path to custom output XML (overrides batch detection)
        blast_dir: Optional path to directory containing BLAST XML files
                    (enables chain BLAST decomposition with alignment data)
    """

    try:
        from pathlib import Path

        # When custom summary_xml is provided, we can skip batch detection for input files
        # but still need reference files from config
        if summary_xml:
            # Use custom paths without batch detection
            # Determine BLAST directory
            if blast_dir:
                # Explicitly provided
                blast_path = Path(blast_dir)
            else:
                # Infer from summary_xml path
                blast_path = Path(summary_xml).parent.parent / "blast"

            paths = {
                'domain_summary': Path(summary_xml),
                'output': Path(output_path) if output_path else Path(f"/tmp/{protein_id}.domains.xml"),
                'batch_name': batch_id or 'custom',
                'batch_dir': Path(summary_xml).parent.parent if summary_xml else Path("/tmp"),
                # Reference files from config
                'domain_lengths': config.domain_lengths_file,
                'protein_lengths': config.protein_lengths_file,
                'domain_definitions': config.domain_definitions_file,
                # Optional BLAST files (may not exist)
                'blast_xml': blast_path / f"{protein_id}.chain.blast.xml",
                'blast_dir': blast_path,
            }
        else:
            # Get all paths with smart batch detection
            paths = config.get_paths_for_protein(protein_id, batch_id, verbose)

            # Override output path if custom one is provided
            if output_path:
                paths['output'] = Path(output_path)

        if verbose:
            print(f"Processing: {protein_id}")
            print(f"Batch: {paths['batch_name']}")
            print(f"Domain summary: {paths['domain_summary']}")
            print(f"BLAST XML: {paths['blast_xml']}")

        # Validate required files exist
        required_files = [
            ("domain_summary", "Domain summary"),
            ("domain_lengths", "Domain lengths"),
            ("protein_lengths", "Protein lengths"),
        ]

        for key, description in required_files:
            if not paths[key].exists():
                print(f"ERROR: {description} file not found: {paths[key]}")
                return None

        # Load reference data
        if verbose:
            print("\nLoading reference data...")

        domain_lengths = load_reference_lengths(str(paths["domain_lengths"]))
        protein_lengths = load_protein_lengths(str(paths["protein_lengths"]))

        # Load domain definitions (optional but recommended)
        domain_definitions = {}
        if paths["domain_definitions"].exists():
            blacklist_path = (
                str(config.reference_blacklist_file)
                if config.reference_blacklist_file.exists()
                else None
            )
            domain_definitions = load_domain_definitions(
                str(paths["domain_definitions"]), verbose=verbose, blacklist_path=blacklist_path
            )
            if verbose:
                print(f"Loaded domain definitions for {len(domain_definitions)} protein chains")
        else:
            print("WARNING: Domain definitions not found - chain BLAST decomposition disabled")

        # Load BLAST alignments (optional but needed for accurate decomposition)
        blast_alignments = {}
        if paths["blast_dir"].exists():
            # Parse protein ID
            parts = protein_id.split("_")
            pdb_id, chain_id = parts[0], parts[1] if len(parts) > 1 else "A"

            blast_alignments = load_chain_blast_alignments(
                str(paths["blast_dir"]), pdb_id, chain_id, verbose=verbose
            )
            if verbose:
                print(f"Loaded {len(blast_alignments)} BLAST alignments")
        else:
            if verbose:
                print("WARNING: BLAST directory not found - alignment-based decomposition disabled")

        # Parse evidence
        if verbose:
            print(f"\nParsing evidence from {paths['domain_summary'].name}...")

        evidence = parse_domain_summary(
            str(paths["domain_summary"]),
            reference_lengths=domain_lengths,
            protein_lengths=protein_lengths,
            blast_alignments=blast_alignments,
            require_reference_lengths=True,
            verbose=verbose,
        )

        if not evidence:
            print("No homology evidence found - protein has 0 domains")

            # Create metadata even for empty results
            parts = protein_id.split("_")
            pdb_id, chain_id = parts[0], parts[1] if len(parts) > 1 else "A"

            metadata = create_metadata_from_batch(
                pdb_id, chain_id, str(paths["batch_dir"]), paths["batch_name"]
            )
            metadata.process_parameters.update(
                {
                    "evidence_items_found": 0,
                    "domains_assigned": 0,
                    "boundary_optimization_enabled": False,
                }
            )

            # Use enhanced writer even for empty results
            write_domain_partition([], metadata, str(paths["output"]))
            print(f"✅ Successfully processed {protein_id} (0 domains)")
            return []  # Success with empty domain list

        # Show evidence summary
        evidence_by_type = {}
        for ev in evidence:
            evidence_by_type[ev.type] = evidence_by_type.get(ev.type, 0) + 1

        print(f"\nFound {len(evidence)} evidence items:")
        for etype, count in sorted(evidence_by_type.items()):
            print(f"  {etype}: {count}")

        # Read sequence length from summary XML (NOT estimated from evidence!)
        import xml.etree.ElementTree as ET
        try:
            tree = ET.parse(str(paths["domain_summary"]))
            root = tree.getroot()
            protein_elem = root.find("protein")

            if protein_elem is not None and protein_elem.get("length"):
                sequence_length = int(protein_elem.get("length"))
                if verbose:
                    print(f"Sequence length from summary XML: {sequence_length}")
            else:
                # Fallback: estimate from evidence (old behavior)
                max_pos = max(ev.query_range.segments[-1].end for ev in evidence)
                sequence_length = int(max_pos * 1.1)
                print(f"WARNING: Could not read sequence length from summary XML, estimating: {sequence_length}")
        except Exception as e:
            # Fallback: estimate from evidence (old behavior)
            max_pos = max(ev.query_range.segments[-1].end for ev in evidence)
            sequence_length = int(max_pos * 1.1)
            print(f"WARNING: Error reading summary XML length ({e}), estimating: {sequence_length}")

        # Show decomposition readiness
        chain_blast_count = evidence_by_type.get("chain_blast", 0)
        with_alignments = sum(1 for e in evidence if e.type == "chain_blast" and e.alignment)

        print("\nDecomposition status:")
        print(f"  Chain BLAST evidence: {chain_blast_count}")
        print(f"  With alignment data: {with_alignments}")
        print(f"  Domain definitions: {'✓' if domain_definitions else '✗'}")

        # Partition domains
        print("\nPartitioning domains...")

        domains = partition_domains(
            evidence,
            sequence_length=sequence_length,
            domain_definitions=domain_definitions if domain_definitions else None,
            verbose=verbose,
        )

        if not domains:
            print("No domains found")
            return None

        # Apply boundary optimization with provenance tracking
        if verbose:
            print("\nApplying boundary optimization...")

        layout = DomainLayout.from_domains(domains, sequence_length)
        optimizer = BoundaryOptimizer()
        optimized_layout = optimizer.optimize_boundaries(
            layout, min_domain_size=25, neighbor_tolerance=5, verbose=verbose
        )

        final_domains = optimized_layout.domains

        # Show results
        print(f"\n{'='*50}")
        print(f"RESULTS: {len(final_domains)} domains found")
        print(f"{'='*50}")

        # Get final statistics from layout
        final_stats = optimized_layout.get_coverage_stats()
        total_coverage = final_stats["assigned_residues"]
        coverage_pct = final_stats["coverage_percent"]

        for i, domain in enumerate(final_domains, 1):
            optimization_info = " (optimized)" if domain.was_optimized() else ""
            print(f"\n{i}. Domain {domain.id}:")
            print(f"   Family: {domain.family}")
            print(f"   Range: {domain.range}")
            print(f"   Size: {domain.range.total_length} residues")
            print(f"   Source: {domain.source}{optimization_info}")
            if domain.range.is_discontinuous:
                print(f"   Discontinuous: {len(domain.range.segments)} segments")

            # Show optimization details if optimized
            if domain.was_optimized() and verbose:
                print(f"   Original range: {domain.original_range}")
                if domain.optimization_actions:
                    print(f"   Optimization actions: {', '.join(domain.optimization_actions)}")

        print(
            f"\nTotal coverage: {total_coverage}/{sequence_length} residues ({coverage_pct:.1f}%)"
        )

        # Create comprehensive metadata and use enhanced writer
        parts = protein_id.split("_")
        pdb_id, chain_id = parts[0], parts[1] if len(parts) > 1 else "A"

        metadata = create_metadata_from_batch(
            pdb_id, chain_id, str(paths["batch_dir"]), paths["batch_name"]
        )
        metadata.sequence_length = sequence_length
        metadata.process_parameters.update(
            {
                "evidence_items_processed": len(evidence),
                "blast_alignments_loaded": len(blast_alignments),
                "domain_definitions_available": len(domain_definitions),
                "reference_lengths_available": len(domain_lengths),
                "boundary_optimization_enabled": True,
                "min_domain_size": 25,
                "neighbor_tolerance": 5,
                "domains_before_optimization": len(domains),
                "domains_after_optimization": len(final_domains),
                "optimization_actions_taken": sum(
                    len(d.optimization_actions) for d in final_domains
                ),
            }
        )

        # Use enhanced writer with comprehensive provenance
        write_domain_partition_from_layout(
            layout=optimized_layout, metadata=metadata, output_path=str(paths["output"])
        )

        print(f"\n✓ Output written to: {paths['output']}")

        # Generate visualization if requested
        if visualize:
            print("\nGenerating PyMOL comparison...")
            try:
                from pyecod_mini.core.visualization import quick_comparison

                script_path = quick_comparison(protein_id, str(paths["batch_dir"]))
                print(f"✓ PyMOL script: {script_path}")
                print(f"  Run: pymol {script_path}")
            except Exception as e:
                print(f"⚠️  Visualization failed: {e}")

        return final_domains

    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        return None
    except Exception as e:
        print(f"ERROR processing {protein_id}: {e}")
        if verbose:
            import traceback

            traceback.print_exc()
        return None


def analyze_protein_batches(protein_id: str, config: PyEcodMiniConfig) -> bool:
    """Analyze a protein across multiple batches"""

    analysis = config.batch_finder.analyze_protein_batches(protein_id)

    if not analysis["multi_batch"]:
        if analysis["batches"]:
            print(f"{protein_id} exists in only one batch: {analysis['batches'][0]}")
        else:
            print(f"{protein_id} not found in any batch")
        return True

    print(f"{protein_id} exists in {len(analysis['batches'])} batches:")
    print("=" * 60)

    for batch_name in analysis["batches"]:
        print(f"  {batch_name}")

    print("\n💡 Recommendations:")
    print(f"  • Test different batches: pyecod-mini {protein_id} --batch-id BATCH_NAME")
    print("  • Compare results to choose the best batch")
    print("  • For test cases, we use: ecod_batch_036_20250406_1424 (known stable)")

    return True
