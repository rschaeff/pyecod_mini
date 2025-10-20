#!/usr/bin/env python3
"""
CLI tests for mini_pyecod

Basic smoke tests for the command-line interface.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


class TestCLIBasics:
    """Test basic CLI functionality"""

    @pytest.mark.unit
    def test_cli_help(self):
        """Test that help works"""
        result = subprocess.run(
            [sys.executable, "-m", "pyecod_mini.cli.main", "--help"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
        )

        assert result.returncode == 0
        assert "pyECOD Mini" in result.stdout or "pyecod-mini" in result.stdout.lower()
        assert "--verbose" in result.stdout
        assert "--batch-id" in result.stdout

    @pytest.mark.unit
    def test_cli_validate(self):
        """Test configuration validation"""
        result = subprocess.run(
            [sys.executable, "-m", "pyecod_mini.cli.main", "--validate"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
        )

        # Should run without crashing
        assert result.returncode == 0
        assert "Configuration" in result.stdout or "configuration" in result.stdout.lower()

    @pytest.mark.unit
    def test_cli_list_batches(self):
        """Test listing batches"""
        result = subprocess.run(
            [sys.executable, "-m", "pyecod_mini.cli.main", "--list-batches"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
        )

        # Should run without crashing
        assert result.returncode == 0
        # Either shows batches or says none found
        assert "batch" in result.stdout.lower()

    @pytest.mark.unit
    def test_cli_missing_protein_id(self):
        """Test error when protein ID is missing"""
        result = subprocess.run(
            [sys.executable, "-m", "pyecod_mini.cli.main"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
        )

        # Should show error
        assert result.returncode != 0
        output = result.stderr + result.stdout
        assert "protein_id" in output.lower() or "protein ID" in output or "required" in output

    @pytest.mark.unit
    def test_cli_version_flag(self):
        """Test that --version displays package version"""
        result = subprocess.run(
            [sys.executable, "-m", "pyecod_mini.cli.main", "--version"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
        )

        # Should succeed and display version
        assert result.returncode == 0
        output = result.stdout + result.stderr
        assert "pyecod" in output.lower() or "2.0.0" in output
        assert "2.0.0" in output  # Current version

    @pytest.mark.integration
    @pytest.mark.slow
    def test_cli_with_protein(self, stable_batch_dir):
        """Test running with a real protein ID"""
        # Only run if we have batch data
        if not Path(stable_batch_dir).exists():
            pytest.skip("No batch directory available")

        result = subprocess.run(
            [sys.executable, "pyecod_mini.py", "8ovp_A", "--batch-id", "036"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
            timeout=60,  # Timeout after 60 seconds
        )

        # Check that it ran (might fail if data missing, but shouldn't crash)
        output = result.stdout + result.stderr

        # Should either succeed or give meaningful error
        assert "RESULTS:" in output or "ERROR:" in output or "not found" in output


class TestBatchFinder:
    """Test the BatchFinder functionality"""

    @pytest.mark.unit
    def test_batch_finder_import(self):
        """Test that BatchFinder can be imported and instantiated"""
        from mini.pyecod_mini import BatchFinder

        finder = BatchFinder("/tmp/test")
        assert finder.base_dir == Path("/tmp/test")
        assert isinstance(finder.stable_batches, dict)

    @pytest.mark.unit
    def test_config_import(self):
        """Test that PyEcodMiniConfig can be imported"""
        from mini.pyecod_mini import PyEcodMiniConfig

        config = PyEcodMiniConfig()
        assert config.base_dir == Path("/data/ecod/pdb_updates/batches")
        assert config.test_data_dir.name == "test_data"


class TestMainFunctions:
    """Test main module functions"""

    @pytest.mark.unit
    def test_partition_protein_import(self):
        """Test that main functions can be imported"""
        from mini.pyecod_mini import analyze_protein_batches, partition_protein

        # Functions should exist
        assert callable(partition_protein)
        assert callable(analyze_protein_batches)


class TestCLICustomPaths:
    """Test CLI with custom input/output paths for pyecod_prod integration"""

    @pytest.mark.integration
    def test_cli_custom_summary_xml_path(self, domain_summary_path, temp_output_dir):
        """Test --summary-xml custom input path"""
        output_path = os.path.join(temp_output_dir, "custom_summary_test.xml")

        # Use custom summary path instead of batch detection
        result = subprocess.run(
            [
                sys.executable,
                "pyecod_mini.py",
                "8ovp_A",
                "--summary-xml",
                domain_summary_path,
                "--output",
                output_path,
            ],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
            timeout=60,
        )

        output = result.stdout + result.stderr

        # Should succeed or give meaningful error
        if result.returncode == 0:
            # If successful, output file should exist
            assert os.path.exists(output_path)
        else:
            # If failed, should have error message
            assert "ERROR" in output or "error" in output

    @pytest.mark.integration
    def test_cli_custom_output_path(self, domain_summary_path, temp_output_dir):
        """Test --output custom output path"""
        custom_output = os.path.join(temp_output_dir, "nested", "custom", "result.xml")

        result = subprocess.run(
            [
                sys.executable,
                "pyecod_mini.py",
                "8ovp_A",
                "--summary-xml",
                domain_summary_path,
                "--output",
                custom_output,
            ],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
            timeout=60,
        )

        # If successful, custom output path should be created
        if result.returncode == 0:
            assert os.path.exists(custom_output)
            # Directory should have been created
            assert os.path.exists(os.path.dirname(custom_output))

    @pytest.mark.integration
    def test_cli_custom_paths_with_batch_id(self, domain_summary_path, temp_output_dir):
        """Test custom paths work with --batch-id"""
        output_path = os.path.join(temp_output_dir, "batch_id_custom_test.xml")

        result = subprocess.run(
            [
                sys.executable,
                "pyecod_mini.py",
                "8ovp_A",
                "--batch-id",
                "036",
                "--summary-xml",
                domain_summary_path,
                "--output",
                output_path,
            ],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
            timeout=60,
        )

        # Should not conflict - custom paths override batch detection
        if result.returncode == 0:
            assert os.path.exists(output_path)

    @pytest.mark.unit
    def test_cli_analyze_batches_flag(self):
        """Test --analyze-batches flag requires protein_id"""
        result = subprocess.run(
            [sys.executable, "-m", "pyecod_mini.cli.main", "--analyze-batches"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
        )

        # Should fail without protein_id
        assert result.returncode != 0
        output = result.stdout + result.stderr
        assert "protein" in output.lower() or "required" in output.lower()


class TestCLIIntegrationWorkflows:
    """Test realistic CLI integration workflows"""

    @pytest.mark.integration
    def test_cli_pyecod_prod_workflow(self, domain_summary_path, temp_output_dir):
        """Test typical pyecod_prod integration workflow"""
        # Simulate pyecod_prod calling pyecod_mini with custom paths
        custom_input = domain_summary_path
        custom_output = os.path.join(temp_output_dir, "pyecod_prod_output.xml")

        result = subprocess.run(
            [
                sys.executable,
                "pyecod_mini.py",
                "8ovp_A",
                "--summary-xml",
                custom_input,
                "--output",
                custom_output,
            ],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
            timeout=60,
        )

        # This should work without batch detection
        if result.returncode == 0:
            # Output should be created at custom path
            assert os.path.exists(custom_output)

            # Output should be valid XML
            import xml.etree.ElementTree as ET

            tree = ET.parse(custom_output)
            root = tree.getroot()
            assert root.tag == "domain_partition"


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])
