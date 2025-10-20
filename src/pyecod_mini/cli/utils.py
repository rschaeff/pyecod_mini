#!/usr/bin/env python3
"""
Utility functions for pyecod_mini CLI

Setup, testing, and helper functions.
"""

import os
from pathlib import Path

from .config import PyEcodMiniConfig


def setup_references(cache_file: str = None, output_dir: str = None) -> bool:
    """Setup reference files from ECOD range cache"""

    if cache_file is None:
        cache_file = "/data/ecod/database_versions/v291/ecod.develop291.range_cache.txt"

    if output_dir is None:
        # Use project test_data directory
        current_dir = Path(__file__).parent
        for parent in [current_dir] + list(current_dir.parents):
            if (parent / "pyproject.toml").exists():
                output_dir = str(parent / "test_data")
                break
        if output_dir is None:
            output_dir = str(current_dir.parent.parent.parent / "test_data")

    print("Setting up reference files...")
    print(f"  Cache file: {cache_file}")
    print(f"  Output directory: {output_dir}")

    if not os.path.exists(cache_file):
        print(f"ERROR: Cache file not found: {cache_file}")
        return False

    try:
        from pyecod_mini.core.range_cache_parser import (
            create_domain_definitions_from_cache,
            create_domain_lengths_from_cache,
            extract_protein_lengths_from_cache,
        )

        os.makedirs(output_dir, exist_ok=True)

        print("Generating domain lengths...")
        create_domain_lengths_from_cache(cache_file, f"{output_dir}/domain_lengths.csv")

        print("Generating domain definitions...")
        create_domain_definitions_from_cache(cache_file, f"{output_dir}/domain_definitions.csv")

        print("Generating protein lengths...")
        extract_protein_lengths_from_cache(cache_file, f"{output_dir}/protein_lengths.csv")

        print("✅ Reference files generated successfully")
        return True

    except Exception as e:
        print(f"ERROR: Failed to generate reference files: {e}")
        return False


def run_test_suite(config: PyEcodMiniConfig, verbose: bool = False) -> bool:
    """Run the formal test suite"""

    print("Running formal test suite...")

    try:
        # Use pytest to run the tests
        import subprocess

        # Find project root
        current_dir = Path(__file__).parent
        project_root = None
        for parent in [current_dir] + list(current_dir.parents):
            if (parent / "pyproject.toml").exists():
                project_root = parent
                break

        if project_root is None:
            print("ERROR: Could not find project root")
            return False

        tests_dir = project_root / "tests"
        if not tests_dir.exists():
            print(f"ERROR: Tests directory not found: {tests_dir}")
            return False

        # Run pytest
        cmd = ["pytest", str(tests_dir), "-v"]
        if verbose:
            cmd.append("-vv")

        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, cwd=str(project_root))

        if result.returncode == 0:
            print("\n✅ All tests passed")
            return True
        print("\n❌ Some tests failed")
        return False

    except Exception as e:
        print(f"ERROR: Failed to run test suite: {e}")
        return False
