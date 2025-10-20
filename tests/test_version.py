#!/usr/bin/env python3
"""
Comprehensive tests for version tracking across pyecod_mini.

Tests version consistency between package metadata, XML output, and API results.
"""

import os
import tempfile
from pathlib import Path

import pytest

import pyecod_mini
from pyecod_mini.core.writer import get_git_version


@pytest.mark.unit
class TestPackageVersion:
    """Test package version metadata"""

    def test_version_exists(self):
        """Test that __version__ is defined"""
        assert hasattr(pyecod_mini, "__version__")
        assert pyecod_mini.__version__ is not None

    def test_version_format(self):
        """Test that version follows semantic versioning"""
        version = pyecod_mini.__version__
        assert isinstance(version, str)

        # Should be in format X.Y.Z
        parts = version.split(".")
        assert len(parts) >= 2  # At least major.minor

        # First two parts should be numbers
        assert parts[0].isdigit()
        assert parts[1].isdigit()

    def test_version_is_2_0_0(self):
        """Test that current version is 2.0.0"""
        assert pyecod_mini.__version__ == "2.0.0"

    def test_version_accessible_from_init(self):
        """Test that version is accessible from package __init__"""
        from pyecod_mini import __version__

        assert __version__ == "2.0.0"


@pytest.mark.unit
class TestWriterVersionFunctions:
    """Test version tracking functions in writer.py"""

    def test_get_git_version_returns_package_version(self):
        """Test that get_git_version() prioritizes package version"""
        version = get_git_version()

        # Should return package version (not git version)
        assert version == pyecod_mini.__version__
        assert version == "2.0.0"

    def test_get_git_version_returns_string(self):
        """Test that get_git_version() always returns a string"""
        version = get_git_version()
        assert isinstance(version, str)
        assert len(version) > 0


@pytest.mark.integration
class TestVersionInXMLOutput:
    """Test that version is correctly written to XML output"""

    def test_version_in_partition_xml(self, domain_summary_path, temp_output_dir):
        """Test that algorithm version is written to partition.xml"""
        from pyecod_mini import partition_protein

        output_path = os.path.join(temp_output_dir, "version_test.xml")

        result = partition_protein(
            summary_xml=domain_summary_path,
            output_xml=output_path,
            pdb_id="8ovp",
            chain_id="A",
        )

        # Read the XML file
        import xml.etree.ElementTree as ET

        tree = ET.parse(output_path)
        root = tree.getroot()

        # Find version element
        metadata = root.find("metadata")
        assert metadata is not None, "Metadata element not found in XML"

        version_elem = metadata.find("version")
        assert version_elem is not None, "Version element not found in metadata"

        # Check algorithm version attribute
        algorithm_version = version_elem.get("algorithm")
        assert algorithm_version is not None, "Algorithm version attribute not found"
        assert algorithm_version == "2.0.0"
        assert algorithm_version == pyecod_mini.__version__

    def test_version_consistency_in_xml(self, domain_summary_path, temp_output_dir):
        """Test that version in XML matches package version"""
        from pyecod_mini import partition_protein

        output_path = os.path.join(temp_output_dir, "version_consistency_test.xml")

        result = partition_protein(
            summary_xml=domain_summary_path,
            output_xml=output_path,
            pdb_id="8ovp",
            chain_id="A",
        )

        import xml.etree.ElementTree as ET

        tree = ET.parse(output_path)
        root = tree.getroot()

        version_elem = root.find("metadata/version")
        xml_version = version_elem.get("algorithm")

        # XML version should match package version
        assert xml_version == pyecod_mini.__version__

        # XML version should also match result version
        assert xml_version == result.algorithm_version


@pytest.mark.integration
class TestVersionConsistency:
    """Test version consistency across all components"""

    def test_version_consistency_across_components(self, domain_summary_path, temp_output_dir):
        """Test that version is consistent across package, API, and XML"""
        from pyecod_mini import __version__, partition_protein

        output_path = os.path.join(temp_output_dir, "consistency_test.xml")

        # Get package version
        package_version = __version__

        # Get version from writer
        writer_version = get_git_version()

        # Run partition and get API result
        result = partition_protein(
            summary_xml=domain_summary_path,
            output_xml=output_path,
            pdb_id="8ovp",
            chain_id="A",
        )
        api_version = result.algorithm_version

        # Read XML
        import xml.etree.ElementTree as ET

        tree = ET.parse(output_path)
        root = tree.getroot()
        xml_version = root.find("metadata/version").get("algorithm")

        # All versions should match
        assert package_version == writer_version == api_version == xml_version == "2.0.0"

    def test_version_in_all_outputs(self, domain_summary_path, temp_output_dir):
        """Test that version appears in all expected outputs"""
        from pyecod_mini import partition_protein

        output_path = os.path.join(temp_output_dir, "all_outputs_test.xml")

        result = partition_protein(
            summary_xml=domain_summary_path,
            output_xml=output_path,
            pdb_id="8ovp",
            chain_id="A",
        )

        # Check result object
        assert result.algorithm_version == "2.0.0"

        # Check XML file exists and has version
        assert os.path.exists(output_path)

        with open(output_path, "r") as f:
            xml_content = f.read()
            assert '2.0.0' in xml_content
            assert 'algorithm="2.0.0"' in xml_content


@pytest.mark.unit
class TestVersionInMetadata:
    """Test version tracking in partition metadata"""

    def test_partition_metadata_includes_version(self):
        """Test that PartitionMetadata stores version information"""
        from pyecod_mini.core.models import PartitionMetadata

        metadata = PartitionMetadata(pdb_id="8abc", chain_id="A")

        # Initially None (set by writer)
        assert metadata.algorithm_version is None

        # Can be set
        metadata.algorithm_version = "2.0.0"
        assert metadata.algorithm_version == "2.0.0"

    def test_writer_sets_version_in_metadata(self, temp_output_dir):
        """Test that writer automatically sets version in metadata"""
        from pyecod_mini.core.models import Domain, Evidence, PartitionMetadata
        from pyecod_mini.core.sequence_range import SequenceRange
        from pyecod_mini.core.writer import write_domain_partition

        # Create minimal evidence
        evidence = Evidence(
            type="chain_blast",
            source_pdb="6dgv",
            query_range=SequenceRange.parse("1-100"),
            confidence=0.85,
        )

        # Create minimal domain
        domain = Domain(
            id="e8ovpA1",
            range=SequenceRange.parse("1-100"),
            family="e6dgvA1",
            source="chain_blast",
            evidence_count=1,
            evidence_items=[evidence],
        )

        # Create metadata without version
        metadata = PartitionMetadata(pdb_id="8ovp", chain_id="A", sequence_length=569)
        assert metadata.algorithm_version is None

        # Write partition
        output_path = os.path.join(temp_output_dir, "metadata_version_test.xml")
        write_domain_partition([domain], metadata, output_path)

        # Metadata should now have version set
        assert metadata.algorithm_version is not None
        assert metadata.algorithm_version == "2.0.0"


@pytest.mark.integration
class TestVersionInCLI:
    """Test version display in CLI"""

    def test_cli_displays_version(self):
        """Test that CLI --version displays correct version"""
        import subprocess
        import sys

        result = subprocess.run(
            [sys.executable, "-m", "pyecod_mini.cli.main", "--version"],
            capture_output=True,
            text=True,
        )

        output = result.stdout + result.stderr
        assert "2.0.0" in output
        assert "pyecod-mini" in output.lower()


@pytest.mark.unit
class TestVersionDocumentation:
    """Test that version is properly documented"""

    def test_version_in_package_docstring(self):
        """Test that package has version info in __init__"""
        import pyecod_mini

        assert hasattr(pyecod_mini, "__version__")
        assert hasattr(pyecod_mini, "__author__")

    def test_version_accessible_from_all_exports(self):
        """Test that version is in __all__ exports"""
        from pyecod_mini import __all__

        assert "__version__" in __all__
