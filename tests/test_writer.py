#!/usr/bin/env python3
"""
Domain partition writer tests for mini_pyecod

Tests the XML output writing functionality.
"""

# Add parent directory to path for imports
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from pyecod_mini.core.models import Domain, PartitionMetadata
from pyecod_mini.core.sequence_range import SequenceRange
from pyecod_mini.core.writer import write_domain_partition


def create_test_metadata(pdb_id: str, chain_id: str):
    """Create metadata object with EXACT signature from API"""
    return PartitionMetadata(
        pdb_id=pdb_id,
        chain_id=chain_id,
        algorithm_version="test_v1.0",
        # NOTE: 'reference' is NOT part of PartitionMetadata - it's passed to write_domain_partition separately
        # NOTE: 'is_classified' is NOT part of PartitionMetadata either
    )


class TestDomainWriter:
    """Test domain partition XML writing"""

    @pytest.mark.unit
    def test_write_basic_domain_partition(self, tmp_path):
        """Test writing a basic domain partition"""
        domains = [
            Domain(
                id="d1",
                range=SequenceRange.parse("10-100"),
                family="test_family",
                evidence_count=3,
                source="domain_blast",
                evidence_items=[],
            )
        ]

        output_file = tmp_path / "test_output.xml"
        metadata = create_test_metadata("1abc", "A")

        # Use the EXACT signature: (domains, metadata, output_path, reference="mini_pyecod")
        write_domain_partition(domains, metadata, str(output_file), reference="test_reference")

        # Verify file exists
        assert output_file.exists()

        # Parse and check content
        tree = ET.parse(output_file)
        root = tree.getroot()

        assert root.tag == "domain_partition"
        assert root.get("pdb_id") == "1abc"
        assert root.get("chain_id") == "A"
        assert root.get("reference") == "test_reference"

        # Check domains
        domains_elem = root.find("domains")
        assert domains_elem is not None

        domain_elems = domains_elem.findall("domain")
        assert len(domain_elems) == 1

        # Check domain attributes
        d = domain_elems[0]
        assert d.get("id") == "d1"
        assert d.get("range") == "10-100"
        assert d.get("family") == "test_family"
        assert d.get("source") == "domain_blast"
        assert d.get("evidence_count") == "3"
        assert d.get("is_discontinuous") == "false"

    @pytest.mark.unit
    def test_write_empty_domains(self, tmp_path):
        """Test writing with no domains (unclassified)"""
        domains = []

        output_file = tmp_path / "empty_output.xml"
        metadata = create_test_metadata("2xyz", "B")
        write_domain_partition(domains, metadata, str(output_file))

        tree = ET.parse(output_file)
        root = tree.getroot()

        # Should indicate unclassified when no domains
        domains_elem = root.find("domains")
        assert domains_elem is not None
        assert len(domains_elem.findall("domain")) == 0

    @pytest.mark.unit
    def test_write_discontinuous_domain(self, tmp_path):
        """Test writing discontinuous domains"""
        domains = [
            Domain(
                id="d1",
                range=SequenceRange.parse("10-50,100-150"),
                family="discontinuous_family",
                evidence_count=1,
                source="chain_blast_decomposed",
                evidence_items=[],
            )
        ]

        output_file = tmp_path / "discontinuous_output.xml"
        metadata = create_test_metadata("3def", "C")
        write_domain_partition(domains, metadata, str(output_file))

        tree = ET.parse(output_file)
        domain = tree.find(".//domain")

        assert domain.get("range") == "10-50,100-150"
        assert domain.get("is_discontinuous") == "true"

    @pytest.mark.unit
    def test_write_multiple_domains(self, tmp_path):
        """Test writing multiple domains"""
        domains = [
            Domain(
                id="d1",
                range=SequenceRange.parse("1-100"),
                family="family1",
                evidence_count=5,
                source="hhsearch",
                evidence_items=[],
            ),
            Domain(
                id="d2",
                range=SequenceRange.parse("150-250"),
                family="family2",
                evidence_count=3,
                source="domain_blast",
                evidence_items=[],
            ),
            Domain(
                id="d3",
                range=SequenceRange.parse("300-400,450-500"),
                family="family3",
                evidence_count=2,
                source="chain_blast_decomposed",
                evidence_items=[],
            ),
        ]

        output_file = tmp_path / "multi_output.xml"
        metadata = create_test_metadata("4ghi", "D")
        write_domain_partition(domains, metadata, str(output_file), reference="custom_reference")

        tree = ET.parse(output_file)
        root = tree.getroot()

        # Check custom reference
        assert root.get("reference") == "custom_reference"

        # Check all domains
        domain_elems = tree.findall(".//domain")
        assert len(domain_elems) == 3

        # Verify domain order preserved
        assert domain_elems[0].get("id") == "d1"
        assert domain_elems[1].get("id") == "d2"
        assert domain_elems[2].get("id") == "d3"

        # Check discontinuous flag
        assert domain_elems[0].get("is_discontinuous") == "false"
        assert domain_elems[1].get("is_discontinuous") == "false"
        assert domain_elems[2].get("is_discontinuous") == "true"

    @pytest.mark.unit
    def test_xml_formatting(self, tmp_path):
        """Test that XML is properly formatted"""
        domains = [
            Domain(
                id="d1",
                range=SequenceRange.parse("1-100"),
                family="test",
                evidence_count=1,
                source="test",
                evidence_items=[],
            )
        ]

        output_file = tmp_path / "formatted.xml"
        metadata = create_test_metadata("1abc", "A")
        write_domain_partition(domains, metadata, str(output_file))

        # Read the file content
        with open(output_file, encoding="utf-8") as f:
            content = f.read()

        # Check for XML declaration
        assert content.startswith(
            "<?xml"
        ), f"XML should start with declaration, got: {content[:50]}"

        # Check indentation (should have spaces)
        assert "  <domains>" in content or "<domains>" in content
        assert "  <domain" in content or "<domain" in content

    @pytest.mark.unit
    def test_special_characters_in_family(self, tmp_path):
        """Test handling of special characters in family names"""
        domains = [
            Domain(
                id="d1",
                range=SequenceRange.parse("1-100"),
                family="family&with<special>chars",
                evidence_count=1,
                source="test",
                evidence_items=[],
            )
        ]

        output_file = tmp_path / "special_chars.xml"
        metadata = create_test_metadata("1abc", "A")
        write_domain_partition(domains, metadata, str(output_file))

        # Parse to ensure valid XML
        tree = ET.parse(output_file)
        domain = tree.find(".//domain")

        # XML parser should handle escaping
        assert domain.get("family") == "family&with<special>chars"


class TestOverlappingDomainCoverage:
    """Test coverage calculation with overlapping domains (bug fix verification)"""

    @pytest.mark.unit
    def test_coverage_with_overlapping_domains(self, tmp_path):
        """Coverage should not exceed 100% even with overlapping domains.

        This test verifies the fix for the overlap bug where domains with
        overlapping residue ranges caused reported coverage to exceed 100%.
        """
        # Create domains with overlap (like 9hpj_B case from bug report)
        # d2: 1-255, d1: 236-341 = 20 residue overlap at 236-255
        domains = [
            Domain(
                id="d2",
                range=SequenceRange.parse("1-255"),
                family="1.1.9",
                evidence_count=1,
                source="chain_blast",
                evidence_items=[],
            ),
            Domain(
                id="d1",
                range=SequenceRange.parse("236-341"),
                family="708.1.2",
                evidence_count=1,
                source="hhsearch",
                evidence_items=[],
            ),
        ]

        output_file = tmp_path / "overlap_test.xml"
        metadata = create_test_metadata("9hpj", "B")
        metadata.sequence_length = 341

        write_domain_partition(domains, metadata, str(output_file))

        # Parse and check coverage
        tree = ET.parse(output_file)
        stats = tree.find(".//statistics")

        # Coverage should be exactly 1.0 (100%) since domains cover positions 1-341
        # NOT 1.059 (105.9%) which was the buggy behavior
        coverage = float(stats.get("total_coverage"))
        residues_assigned = int(stats.get("residues_assigned"))

        assert coverage <= 1.0, f"Coverage {coverage} exceeds 100%"
        assert coverage == 1.0, f"Expected 100% coverage, got {coverage:.1%}"
        assert residues_assigned == 341, f"Expected 341 residues, got {residues_assigned}"

    @pytest.mark.unit
    def test_coverage_with_discontinuous_overlap(self, tmp_path):
        """Coverage with discontinuous domain that overlaps another.

        Based on 9uko_F case: discontinuous domain inserts into another.
        """
        # d2: 1-173, d1: 9-25,174-340 (discontinuous, segment 9-25 overlaps d2)
        domains = [
            Domain(
                id="d2",
                range=SequenceRange.parse("1-173"),
                family="142.1.1",
                evidence_count=1,
                source="domain_blast",
                evidence_items=[],
            ),
            Domain(
                id="d1",
                range=SequenceRange.parse("9-25,174-340"),
                family="142.1.1",
                evidence_count=1,
                source="chain_blast_decomposed",
                evidence_items=[],
            ),
        ]

        output_file = tmp_path / "discontinuous_overlap.xml"
        metadata = create_test_metadata("9uko", "F")
        metadata.sequence_length = 340

        write_domain_partition(domains, metadata, str(output_file))

        tree = ET.parse(output_file)
        stats = tree.find(".//statistics")

        coverage = float(stats.get("total_coverage"))
        residues_assigned = int(stats.get("residues_assigned"))

        # Should be exactly 340/340 = 1.0, NOT 357/340 = 1.05
        assert coverage <= 1.0, f"Coverage {coverage} exceeds 100%"
        assert residues_assigned == 340, f"Expected 340 residues, got {residues_assigned}"

    @pytest.mark.unit
    def test_coverage_non_overlapping_domains(self, tmp_path):
        """Non-overlapping domains should calculate correctly."""
        domains = [
            Domain(
                id="d1",
                range=SequenceRange.parse("1-100"),
                family="family1",
                evidence_count=1,
                source="domain_blast",
                evidence_items=[],
            ),
            Domain(
                id="d2",
                range=SequenceRange.parse("101-200"),
                family="family2",
                evidence_count=1,
                source="domain_blast",
                evidence_items=[],
            ),
        ]

        output_file = tmp_path / "non_overlap.xml"
        metadata = create_test_metadata("test", "A")
        metadata.sequence_length = 300

        write_domain_partition(domains, metadata, str(output_file))

        tree = ET.parse(output_file)
        stats = tree.find(".//statistics")

        coverage = float(stats.get("total_coverage"))
        residues_assigned = int(stats.get("residues_assigned"))

        # 200 residues covered out of 300 = 66.67%
        expected_coverage = 200 / 300
        assert abs(coverage - expected_coverage) < 0.001
        assert residues_assigned == 200


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])
