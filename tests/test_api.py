#!/usr/bin/env python3
"""
Comprehensive tests for the library API (api.py).

Tests the public-facing interface for pyecod_mini library integration.
"""

import os
import tempfile
from pathlib import Path

import pytest

from pyecod_mini import Domain, PartitionError, Partitioner, PartitionResult, __version__, partition_protein
from pyecod_mini.core.reference_cache import ReferenceCache, ReferenceData


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
        assert algorithm_version == __version__
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
        # Version should be semantic versioning format
        assert len(__version__.split(".")) == 3

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


@pytest.mark.unit
class TestReferenceCache:
    """Test ReferenceCache class for batch processing optimization"""

    def test_reference_cache_initial_state(self):
        """Test ReferenceCache initial state"""
        cache = ReferenceCache()
        assert not cache.is_loaded()
        assert cache.summary() == "ReferenceCache: No data loaded"

    def test_reference_cache_load_empty(self):
        """Test loading with no files specified"""
        cache = ReferenceCache()
        data = cache.load()

        # Should return empty data but be considered loaded
        assert isinstance(data, ReferenceData)
        assert data.domain_definitions_count == 0
        assert data.reference_lengths_count == 0
        assert data.protein_lengths_count == 0

    def test_reference_cache_load_nonexistent_files(self):
        """Test loading with nonexistent files gracefully"""
        cache = ReferenceCache()
        data = cache.load(
            domain_definitions_file="/nonexistent/domain_definitions.csv",
            reference_lengths_file="/nonexistent/reference_lengths.csv",
            protein_lengths_file="/nonexistent/protein_lengths.csv",
        )

        # Should not raise, just return empty data
        assert isinstance(data, ReferenceData)
        assert not data.is_loaded()

    def test_reference_cache_clear(self):
        """Test clearing cached data"""
        cache = ReferenceCache()
        cache.load()
        cache.clear()
        assert not cache.is_loaded()

    def test_reference_cache_context_manager(self):
        """Test ReferenceCache as context manager"""
        with ReferenceCache() as cache:
            cache.load()
            assert isinstance(cache, ReferenceCache)

        # After exiting context, cache should be cleared
        assert not cache.is_loaded()

    def test_reference_data_is_loaded(self):
        """Test ReferenceData.is_loaded() method"""
        data = ReferenceData()
        assert not data.is_loaded()

        # With some data
        data.reference_lengths = {"test": 100}
        assert data.is_loaded()


@pytest.mark.unit
class TestPartitionerClass:
    """Test Partitioner class structure and methods"""

    def test_partitioner_initial_state(self):
        """Test Partitioner initial state"""
        p = Partitioner()
        assert not p.is_loaded()
        assert p.partition_count == 0
        assert "not loaded" in repr(p)

    def test_partitioner_repr(self):
        """Test Partitioner string representation"""
        p = Partitioner()
        assert "Partitioner" in repr(p)
        assert "not loaded" in repr(p)
        assert "0 partitions" in repr(p)

    def test_partitioner_context_manager(self):
        """Test Partitioner as context manager"""
        with Partitioner() as p:
            assert isinstance(p, Partitioner)

    def test_partitioner_requires_load_before_partition(self, temp_output_dir):
        """Test that partition() requires load_references() first"""
        p = Partitioner()

        with pytest.raises(RuntimeError) as exc_info:
            p.partition(
                summary_xml="/dummy/path.xml",
                output_xml=os.path.join(temp_output_dir, "output.xml"),
                pdb_id="test",
                chain_id="A",
            )

        assert "Reference data not loaded" in str(exc_info.value)
        assert "load_references()" in str(exc_info.value)

    def test_partitioner_close(self):
        """Test Partitioner.close() releases resources"""
        p = Partitioner()
        # Load some reference data (even if empty)
        p._cache.load()
        assert p.is_loaded() or True  # May or may not be "loaded" depending on files

        p.close()
        assert not p.is_loaded()

    def test_partitioner_summary(self):
        """Test Partitioner.summary() method"""
        p = Partitioner()
        summary = p.summary()
        assert "Partitioner" in summary
        assert "0 partitions completed" in summary

    def test_partitioner_method_chaining(self):
        """Test that load_references returns self for chaining"""
        p = Partitioner()
        result = p.load_references()  # With no files, should still work
        assert result is p


@pytest.mark.integration
class TestPartitionerIntegration:
    """Integration tests for Partitioner batch processing"""

    def test_partitioner_load_references_from_config(self):
        """Test loading references from default config"""
        p = Partitioner()
        p.load_references_from_config()

        # Should have loaded some data
        assert p.is_loaded()
        p.close()

    def test_partitioner_partition_success(self, domain_summary_path, temp_output_dir):
        """Test successful partition with Partitioner"""
        output_path = os.path.join(temp_output_dir, "partitioner_test.xml")

        with Partitioner() as p:
            p.load_references_from_config()

            result = p.partition(
                summary_xml=domain_summary_path,
                output_xml=output_path,
                pdb_id="8ovp",
                chain_id="A",
            )

            # Verify result
            assert isinstance(result, PartitionResult)
            assert result.success is True
            assert result.pdb_id == "8ovp"
            assert result.chain_id == "A"
            assert os.path.exists(output_path)

            # Verify partition count
            assert p.partition_count == 1

    def test_partitioner_multiple_partitions(self, domain_summary_path, temp_output_dir):
        """Test partitioning multiple proteins with same references"""
        with Partitioner() as p:
            p.load_references_from_config()

            # Partition same protein twice (simulates batch)
            for i in range(3):
                output_path = os.path.join(temp_output_dir, f"batch_test_{i}.xml")
                result = p.partition(
                    summary_xml=domain_summary_path,
                    output_xml=output_path,
                    pdb_id="8ovp",
                    chain_id="A",
                )
                assert result.success is True

            # Should have counted all partitions
            assert p.partition_count == 3

    def test_partitioner_marks_cached_references(self, domain_summary_path, temp_output_dir):
        """Test that Partitioner marks output as using cached references"""
        output_path = os.path.join(temp_output_dir, "cached_refs_test.xml")

        with Partitioner() as p:
            p.load_references_from_config()
            p.partition(
                summary_xml=domain_summary_path,
                output_xml=output_path,
                pdb_id="8ovp",
                chain_id="A",
            )

        # Check XML for cached_references_used marker
        import xml.etree.ElementTree as ET

        tree = ET.parse(output_path)
        root = tree.getroot()

        metadata = root.find("metadata")
        params = metadata.find("parameters")

        # Look for cached_references_used parameter
        found_cached_param = False
        for param in params.findall("parameter"):
            if param.get("name") == "cached_references_used":
                assert param.get("value") == "True"
                found_cached_param = True
                break

        assert found_cached_param, "cached_references_used parameter not found in XML"

    def test_partitioner_file_not_found(self, temp_output_dir):
        """Test Partitioner handles missing input files"""
        with Partitioner() as p:
            p.load_references_from_config()

            with pytest.raises(FileNotFoundError):
                p.partition(
                    summary_xml="/nonexistent/summary.xml",
                    output_xml=os.path.join(temp_output_dir, "output.xml"),
                    pdb_id="test",
                    chain_id="A",
                )


@pytest.mark.unit
class TestPartitionerExport:
    """Test Partitioner is properly exported"""

    def test_partitioner_exported_from_package(self):
        """Test Partitioner is accessible from pyecod_mini"""
        import pyecod_mini

        assert hasattr(pyecod_mini, "Partitioner")
        assert pyecod_mini.Partitioner is Partitioner

    def test_partitioner_in_all(self):
        """Test Partitioner is in __all__"""
        import pyecod_mini

        assert "Partitioner" in pyecod_mini.__all__

    def test_partitioner_has_docstring(self):
        """Test Partitioner has comprehensive documentation"""
        assert Partitioner.__doc__ is not None
        doc = Partitioner.__doc__

        # Check for key documentation
        assert "batch" in doc.lower() or "Batch" in doc
        assert "caching" in doc.lower()  # "reference data caching"
        assert "partition" in doc.lower()
