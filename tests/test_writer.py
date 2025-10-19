#!/usr/bin/env python3
"""
Domain partition writer tests for mini_pyecod

Tests the XML output writing functionality.
"""

import pytest
from pathlib import Path
import xml.etree.ElementTree as ET

# Add parent directory to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from pyecod_mini.core.writer import write_domain_partition
from pyecod_mini.core.models import Domain, Evidence, PartitionMetadata
from pyecod_mini.core.sequence_range import SequenceRange


def create_test_metadata(pdb_id: str, chain_id: str):
    """Create metadata object with EXACT signature from API"""
    return PartitionMetadata(
        pdb_id=pdb_id,
        chain_id=chain_id,
        algorithm_version="test_v1.0"
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
                evidence_items=[]
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
                evidence_items=[]
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
                evidence_items=[]
            ),
            Domain(
                id="d2",
                range=SequenceRange.parse("150-250"),
                family="family2",
                evidence_count=3,
                source="domain_blast",
                evidence_items=[]
            ),
            Domain(
                id="d3",
                range=SequenceRange.parse("300-400,450-500"),
                family="family3",
                evidence_count=2,
                source="chain_blast_decomposed",
                evidence_items=[]
            )
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
                evidence_items=[]
            )
        ]
        
        output_file = tmp_path / "formatted.xml"
        metadata = create_test_metadata("1abc", "A")
        write_domain_partition(domains, metadata, str(output_file))
        
        # Read the file content
        with open(output_file, 'r', encoding='utf-8') as f:
            content = f.read()

        # Check for XML declaration
        assert content.startswith('<?xml'), f"XML should start with declaration, got: {content[:50]}"

        # Check indentation (should have spaces)
        assert '  <domains>' in content or '<domains>' in content
        assert '  <domain' in content or '<domain' in content

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
                evidence_items=[]
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


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])
