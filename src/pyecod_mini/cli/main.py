#!/usr/bin/env python3
"""
pyecod_mini - Clean Minimal Domain Partitioning Tool

A standalone domain partitioning tool with smart batch detection,
visualization, and testing capabilities.

Usage:
    pyecod-mini 8ovp_A                    # Basic domain partitioning
    pyecod-mini 8ovp_A --visualize        # With PyMOL comparison
    pyecod-mini --test-suite               # Run formal test cases
    pyecod-mini --setup-references         # Generate reference files
"""

import argparse
import sys

from .config import PyEcodMiniConfig
from .partition import analyze_protein_batches, partition_protein
from .utils import run_test_suite, setup_references


def main():
    """Main entry point"""

    parser = argparse.ArgumentParser(
        description="pyECOD Mini - Clean Domain Partitioning Tool with Enhanced Provenance Tracking",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  pyecod-mini 8ovp_A                    # Basic domain partitioning
  pyecod-mini 8ovp_A --visualize        # With PyMOL comparison
  pyecod-mini 8ovp_A --batch-id 036     # Use specific batch
  pyecod-mini 8ovp_A --analyze-batches  # Check if protein exists in multiple batches
  pyecod-mini --test-suite               # Run formal test cases
  pyecod-mini --setup-references         # Generate reference files
  pyecod-mini --list-batches             # Show available batches

Enhanced Features:
  • Comprehensive provenance tracking with git versioning
  • Evidence chain documentation and boundary optimization audit trails
  • File integrity verification with SHA256 hashes
  • Algorithm parameter and performance metrics capture
        """,
    )

    # Main action
    parser.add_argument("protein_id", nargs="?", help="Protein ID to process (e.g., 8ovp_A)")

    # Options
    parser.add_argument("--batch-id", help='Batch ID (number like "036" or full name)')
    parser.add_argument("-v", "--verbose", action="store_true", help="Show detailed output")
    parser.add_argument(
        "--visualize", action="store_true", help="Generate PyMOL comparison visualization"
    )

    # Custom input/output paths (for integration with pyecod_prod)
    parser.add_argument(
        "--summary-xml", help="Path to input domain summary XML file (overrides batch detection)"
    )
    parser.add_argument(
        "--output", help="Path to output partition XML file (overrides batch detection)"
    )

    # Utility commands
    parser.add_argument(
        "--list-batches", action="store_true", help="List available batch directories"
    )
    parser.add_argument(
        "--validate", action="store_true", help="Validate configuration and show status"
    )
    parser.add_argument("--test-suite", action="store_true", help="Run formal test suite")
    parser.add_argument(
        "--setup-references", action="store_true", help="Generate reference files from ECOD cache"
    )
    parser.add_argument(
        "--analyze-batches", action="store_true", help="Analyze protein across multiple batches"
    )
    parser.add_argument(
        "--cache-file",
        default="/data/ecod/database_versions/v291/ecod.develop291.range_cache.txt",
        help="ECOD range cache file (for --setup-references)",
    )

    args = parser.parse_args()

    # Initialize configuration
    config = PyEcodMiniConfig()

    # Handle utility commands
    if args.list_batches:
        batch_info = config.list_available_batches()
        if batch_info:
            print("Available batches:")
            for batch_name, protein_count in batch_info:
                print(f"  {batch_name} ({protein_count:,} proteins)")
        else:
            print("No batch directories found")
        return

    if args.validate:
        print("Validating pyECOD Mini configuration...")
        print(f"Base directory: {config.base_dir}")
        print(f"Test data directory: {config.test_data_dir}")

        # Show some available batches
        batch_info = config.list_available_batches()
        if batch_info:
            print(f"Available batches: {len(batch_info)} found")
            # Show most recent few
            for batch_name, protein_count in batch_info[-3:]:
                print(f"  {batch_name} ({protein_count:,} proteins)")

        is_valid, issues = config.validate_setup()
        if is_valid:
            print("✓ Configuration is valid")
        else:
            print("✗ Configuration issues found:")
            for issue in issues:
                print(f"  - {issue}")
        return

    if args.setup_references:
        # Find output directory
        from pathlib import Path

        current_dir = Path(__file__).parent
        project_root = None
        for parent in [current_dir] + list(current_dir.parents):
            if (parent / "pyproject.toml").exists():
                project_root = parent
                break
        output_dir = str(project_root / "test_data") if project_root else None

        success = setup_references(args.cache_file, output_dir)
        if not success:
            sys.exit(1)
        return

    if args.test_suite:
        success = run_test_suite(config, args.verbose)
        if not success:
            sys.exit(1)
        return

    if args.analyze_batches:
        if not args.protein_id:
            print("ERROR: --analyze-batches requires a protein_id")
            print("Usage: pyecod-mini PROTEIN_ID --analyze-batches")
            sys.exit(1)

        success = analyze_protein_batches(args.protein_id, config)
        if not success:
            sys.exit(1)
        return

    # Main processing
    if not args.protein_id:
        parser.print_help()
        print("\nERROR: protein_id is required for domain partitioning")
        print("Use --test-suite to run tests or --setup-references to prepare data")
        sys.exit(1)

    # Validate configuration
    is_valid, issues = config.validate_setup()
    if not is_valid:
        print("Configuration issues found:")
        for issue in issues:
            print(f"  - {issue}")
        print("\nUse --validate for more details or --setup-references to fix")
        sys.exit(1)

    # Process protein
    result = partition_protein(
        args.protein_id,
        config,
        args.batch_id,
        args.verbose,
        args.visualize,
        summary_xml=args.summary_xml,
        output_path=args.output,
    )

    if result is None:
        sys.exit(1)
    else:
        print(f"\n✅ Successfully processed {args.protein_id}")
        if args.visualize:
            print("✅ Visualization generated")


if __name__ == "__main__":
    main()
