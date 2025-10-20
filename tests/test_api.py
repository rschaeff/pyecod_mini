#!/usr/bin/env python3
"""
Comprehensive tests for the library API (api.py).

Tests the public-facing interface for pyecod_mini library integration.
"""

import os
import tempfile
from pathlib import Path

import pytest

from pyecod_mini import Domain, PartitionError, PartitionResult, __version__, partition_protein


@pytest.mark.unit
class TestAPIDataclasses:
    """Test API dataclass structures"""

    def test_domain_dataclass_structure(self):
        """Verify Domain dataclass has correct fields"""
        domain = Domain(
            domain_id="e8ovpA1",
            range_string="10-110",
            residue_count=101,
            source="chain_blast",
            family_name="e6dgvA1",
            confidence=0.95,
        )

        assert domain.domain_id == "e8ovpA1"
        assert domain.range_string == "10-110"
        assert domain.residue_count == 101
        assert domain.source == "chain_blast"
        assert domain.family_name == "e6dgvA1"
        assert domain.confidence == 0.95

    def test_domain_dataclass_optional_confidence(self):
        """Verify Domain confidence is optional"""
        domain = Domain(
            domain_id="e8ovpA1",
            range_string="10-110",
            residue_count=101,
            source="chain_blast",
            family_name="e6dgvA1",
        )

        assert domain.confidence is None

    def test_partition_result_dataclass_structure(self):
        """Verify PartitionResult dataclass has correct fields"""
        result = PartitionResult(
            success=True,
            pdb_id="8ovp",
            chain_id="A",
            sequence_length=569,
            domains=[],
            coverage=0.85,
            partition_xml_path="/path/to/output.xml",
            algorithm_version="2.0.0",
            error_message=None,
        )

        assert result.success is True
        assert result.pdb_id == "8ovp"
        assert result.chain_id == "A"
        assert result.sequence_length == 569
        assert result.domains == []
        assert result.coverage == 0.85
        assert result.partition_xml_path == "/path/to/output.xml"
        assert result.algorithm_version == "2.0.0"
        assert result.error_message is None

    def test_partition_result_with_error(self):
        """Verify PartitionResult can represent errors"""
        result = PartitionResult(
            success=False,
            pdb_id="8abc",
            chain_id="A",
            sequence_length=0,
            domains=[],
            coverage=0.0,
            partition_xml_path="/path/to/output.xml",
            algorithm_version="2.0.0",
            error_message="Partitioning failed: test error",
        )

        assert result.success is False
        assert result.error_message is not None
        assert "test error" in result.error_message


@pytest.mark.unit
class TestAPIExceptions:
    """Test API exception classes"""

    def test_partition_error_exception(self):
        """Verify PartitionError can be raised and caught"""
        with pytest.raises(PartitionError) as exc_info:
            raise PartitionError("Test error message")

        assert "Test error message" in str(exc_info.value)

    def test_partition_error_with_cause(self):
        """Verify PartitionError can wrap other exceptions"""
        original_error = ValueError("Original error")

        with pytest.raises(PartitionError) as exc_info:
            try:
                raise original_error
            except ValueError as e:
                raise PartitionError("Wrapped error") from e

        assert exc_info.value.__cause__ is original_error


@pytest.mark.integration
class TestAPIPartitionProtein:
    """Test partition_protein() function with real data"""

    def test_partition_protein_success(self, domain_summary_path, temp_output_dir):
        """Test successful partition_protein() call"""
        output_path = os.path.join(temp_output_dir, "8ovp_A_api_test.partition.xml")

        result = partition_protein(
            summary_xml=domain_summary_path,
            output_xml=output_path,
            pdb_id="8ovp",
            chain_id="A",
        )

        # Verify result structure
        assert isinstance(result, PartitionResult)
        assert result.success is True
        assert result.pdb_id == "8ovp"
        assert result.chain_id == "A"
        assert result.error_message is None

        # Verify output file was created
        assert os.path.exists(output_path)
        assert result.partition_xml_path == output_path

        # Verify domains were found
        assert len(result.domains) > 0
        assert all(isinstance(d, Domain) for d in result.domains)

        # Verify coverage is reasonable
        assert 0.0 <= result.coverage <= 1.0

        # Verify sequence length is positive
        assert result.sequence_length > 0

        # Verify version is set
        assert result.algorithm_version == __version__

    def test_partition_protein_with_batch_id(self, domain_summary_path, temp_output_dir):
        """Test partition_protein() with batch_id tracking"""
        output_path = os.path.join(temp_output_dir, "8ovp_A_batch_test.partition.xml")

        result = partition_protein(
            summary_xml=domain_summary_path,
            output_xml=output_path,
            pdb_id="8ovp",
            chain_id="A",
            batch_id="test_batch_001",
        )

        assert result.success is True
        assert os.path.exists(output_path)

        # Verify batch_id was written to XML
        with open(output_path, "r") as f:
            xml_content = f.read()
            assert "test_batch_001" in xml_content

    def test_partition_protein_domain_structure(self, domain_summary_path, temp_output_dir):
        """Test that returned Domain objects have correct structure"""
        output_path = os.path.join(temp_output_dir, "8ovp_A_domain_test.partition.xml")

        result = partition_protein(
            summary_xml=domain_summary_path,
            output_xml=output_path,
            pdb_id="8ovp",
            chain_id="A",
        )

        assert len(result.domains) > 0

        # Check first domain structure
        domain = result.domains[0]
        assert isinstance(domain.domain_id, str)
        assert len(domain.domain_id) > 0

        assert isinstance(domain.range_string, str)
        assert len(domain.range_string) > 0

        assert isinstance(domain.residue_count, int)
        assert domain.residue_count > 0

        assert isinstance(domain.source, str)
        assert domain.source in [
            "chain_blast",
            "domain_blast",
            "hhsearch",
            "chain_blast_decomposed",
        ]

        assert isinstance(domain.family_name, str)

        # Confidence can be None or float
        assert domain.confidence is None or isinstance(domain.confidence, float)

    def test_partition_protein_coverage_calculation(self, domain_summary_path, temp_output_dir):
        """Test that coverage is calculated correctly"""
        output_path = os.path.join(temp_output_dir, "8ovp_A_coverage_test.partition.xml")

        result = partition_protein(
            summary_xml=domain_summary_path,
            output_xml=output_path,
            pdb_id="8ovp",
            chain_id="A",
        )

        # Coverage should be between 0 and 1
        assert 0.0 <= result.coverage <= 1.0

        # If domains exist, coverage should be > 0
        if result.domains:
            assert result.coverage > 0.0

        # Calculate coverage manually from domains
        if result.domains and result.sequence_length > 0:
            total_residues = sum(d.residue_count for d in result.domains)
            expected_coverage = total_residues / result.sequence_length
            # Allow small floating point differences
            assert abs(result.coverage - expected_coverage) < 0.01

    def test_partition_protein_version_in_result(self, domain_summary_path, temp_output_dir):
        """Test that algorithm version is included in result"""
        output_path = os.path.join(temp_output_dir, "8ovp_A_version_test.partition.xml")

        result = partition_protein(
            summary_xml=domain_summary_path,
            output_xml=output_path,
            pdb_id="8ovp",
            chain_id="A",
        )

        # Verify version matches package version
        assert result.algorithm_version == __version__
        assert result.algorithm_version == "2.0.0"

    def test_partition_protein_version_in_xml(self, domain_summary_path, temp_output_dir):
        """Test that algorithm version is written to XML output"""
        output_path = os.path.join(temp_output_dir, "8ovp_A_xml_version_test.partition.xml")

        result = partition_protein(
            summary_xml=domain_summary_path,
            output_xml=output_path,
            pdb_id="8ovp",
            chain_id="A",
        )

        # Read XML and verify version is present
        import xml.etree.ElementTree as ET

        tree = ET.parse(output_path)
        root = tree.getroot()

        metadata = root.find("metadata")
        assert metadata is not None

        version_elem = metadata.find("version")
        assert version_elem is not None

        algorithm_version = version_elem.get("algorithm")
        assert algorithm_version == "2.0.0"
        assert algorithm_version == result.algorithm_version


@pytest.mark.unit
class TestAPIErrorHandling:
    """Test error handling in partition_protein()"""

    def test_partition_protein_missing_summary_xml(self, temp_output_dir):
        """Test FileNotFoundError when summary_xml doesn't exist"""
        nonexistent_summary = "/nonexistent/path/summary.xml"
        output_path = os.path.join(temp_output_dir, "output.xml")

        with pytest.raises(FileNotFoundError) as exc_info:
            partition_protein(
                summary_xml=nonexistent_summary,
                output_xml=output_path,
                pdb_id="8abc",
                chain_id="A",
            )

        assert "Summary XML not found" in str(exc_info.value)
        assert nonexistent_summary in str(exc_info.value)

    def test_partition_protein_creates_output_directory(self, domain_summary_path):
        """Test that partition_protein creates output directory if needed"""
        with tempfile.TemporaryDirectory() as tmpdir:
            nested_output = os.path.join(tmpdir, "nested", "subdir", "output.xml")

            result = partition_protein(
                summary_xml=domain_summary_path,
                output_xml=nested_output,
                pdb_id="8ovp",
                chain_id="A",
            )

            # Verify nested directory was created
            assert os.path.exists(os.path.dirname(nested_output))
            assert os.path.exists(nested_output)
            assert result.success is True

    def test_partition_protein_invalid_xml(self, temp_output_dir):
        """Test handling of invalid XML input"""
        # Create a file with invalid XML content
        invalid_xml = os.path.join(temp_output_dir, "invalid.xml")
        with open(invalid_xml, "w") as f:
            f.write("This is not valid XML!")

        output_path = os.path.join(temp_output_dir, "output.xml")

        # API handles parse errors gracefully - returns success but with 0 domains
        result = partition_protein(
            summary_xml=invalid_xml,
            output_xml=output_path,
            pdb_id="8abc",
            chain_id="A",
        )

        # Should succeed but with no domains found
        assert result.success is True
        assert len(result.domains) == 0
        assert os.path.exists(output_path)


@pytest.mark.unit
class TestAPIVersionExport:
    """Test that version is properly exported from package"""

    def test_version_exported(self):
        """Test that __version__ is accessible from package"""
        from pyecod_mini import __version__

        assert __version__ is not None
        assert isinstance(__version__, str)
        assert __version__ == "2.0.0"

    def test_api_exports(self):
        """Test that API functions are exported from package __init__"""
        import pyecod_mini

        # Check all expected exports
        assert hasattr(pyecod_mini, "partition_protein")
        assert hasattr(pyecod_mini, "PartitionResult")
        assert hasattr(pyecod_mini, "PartitionError")
        assert hasattr(pyecod_mini, "Domain")
        assert hasattr(pyecod_mini, "__version__")

        # Verify they're the correct types
        assert callable(pyecod_mini.partition_protein)
        assert isinstance(pyecod_mini.__version__, str)


@pytest.mark.integration
class TestAPICustomPaths:
    """Test partition_protein() with custom input/output paths"""

    def test_partition_protein_custom_paths(self, domain_summary_path, temp_output_dir):
        """Test that custom paths are correctly used"""
        custom_output = os.path.join(temp_output_dir, "custom", "path", "result.xml")

        result = partition_protein(
            summary_xml=domain_summary_path,
            output_xml=custom_output,
            pdb_id="8ovp",
            chain_id="A",
        )

        # Verify custom output path was used
        assert result.partition_xml_path == custom_output
        assert os.path.exists(custom_output)
        assert result.success is True

    def test_partition_protein_preserves_paths(self, domain_summary_path, temp_output_dir):
        """Test that exact paths are preserved in results"""
        # Use absolute path
        abs_output = os.path.abspath(os.path.join(temp_output_dir, "absolute_path.xml"))

        result = partition_protein(
            summary_xml=domain_summary_path,
            output_xml=abs_output,
            pdb_id="8ovp",
            chain_id="A",
        )

        assert result.partition_xml_path == abs_output
        assert os.path.isabs(result.partition_xml_path)


@pytest.mark.integration
class TestAPIIntegrationScenarios:
    """Test realistic integration scenarios"""

    def test_api_basic_workflow(self, domain_summary_path, temp_output_dir):
        """Test basic workflow: read summary, partition, check results"""
        output_path = os.path.join(temp_output_dir, "workflow_test.xml")

        # Step 1: Partition protein
        result = partition_protein(
            summary_xml=domain_summary_path,
            output_xml=output_path,
            pdb_id="8ovp",
            chain_id="A",
        )

        # Step 2: Verify results
        assert result.success is True
        assert len(result.domains) > 0

        # Step 3: Check output file
        assert os.path.exists(result.partition_xml_path)
        assert os.path.getsize(result.partition_xml_path) > 0

        # Step 4: Verify we can access domain details
        for domain in result.domains:
            assert domain.domain_id is not None
            assert domain.residue_count > 0
            assert domain.range_string is not None

    def test_api_result_serialization(self, domain_summary_path, temp_output_dir):
        """Test that PartitionResult can be easily serialized"""
        output_path = os.path.join(temp_output_dir, "serialization_test.xml")

        result = partition_protein(
            summary_xml=domain_summary_path,
            output_xml=output_path,
            pdb_id="8ovp",
            chain_id="A",
        )

        # Verify result can be converted to dict (for JSON serialization)
        from dataclasses import asdict

        result_dict = asdict(result)

        assert result_dict["success"] is True
        assert result_dict["pdb_id"] == "8ovp"
        assert result_dict["chain_id"] == "A"
        assert isinstance(result_dict["domains"], list)

        # Verify domains are serializable too
        if result_dict["domains"]:
            domain_dict = result_dict["domains"][0]
            assert "domain_id" in domain_dict
            assert "range_string" in domain_dict
            assert "residue_count" in domain_dict


@pytest.mark.unit
class TestAPIDocumentation:
    """Test that API has proper documentation"""

    def test_partition_protein_has_docstring(self):
        """Verify partition_protein has comprehensive docstring"""
        assert partition_protein.__doc__ is not None
        doc = partition_protein.__doc__

        # Check for key documentation elements
        assert "Args:" in doc or "Parameters:" in doc
        assert "Returns:" in doc
        assert "Raises:" in doc
        assert "Example:" in doc or "Examples:" in doc

        # Check for parameter documentation
        assert "summary_xml" in doc
        assert "output_xml" in doc
        assert "pdb_id" in doc
        assert "chain_id" in doc

    def test_exception_classes_have_docstrings(self):
        """Verify exception classes have docstrings"""
        assert PartitionError.__doc__ is not None

    def test_dataclasses_have_docstrings(self):
        """Verify dataclasses have docstrings"""
        assert Domain.__doc__ is not None
        assert PartitionResult.__doc__ is not None
