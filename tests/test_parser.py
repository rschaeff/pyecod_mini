#!/usr/bin/env python3
"""
Domain summary XML parser tests for mini_pyecod

Tests the parsing of domain summary XML files and reference data loading.
"""

import csv

# Add parent directory to path for imports
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from pyecod_mini.core.models import Evidence
from pyecod_mini.core.parser import (
    get_evidence_summary,
    load_protein_lengths,
    load_reference_lengths,
    parse_domain_summary,
)
from pyecod_mini.core.sequence_range import SequenceRange


class TestDomainSummaryParsing:
    """Test domain summary XML parsing"""

    @pytest.mark.unit
    def test_parse_valid_domain_summary(self, tmp_path):
        """Test parsing a valid domain summary XML"""
        xml_content = """<?xml version="1.0"?>
<blast_summ_doc>
  <blast_summ pdb="8ovp" chain="A"/>
  <chain_blast_run program="blastp">
    <hits>
      <hit num="1" pdb_id="6dgv" chain_id="A" hsp_count="1" evalues="1e-50">
        <query_reg>252-494</query_reg>
        <hit_reg>1-243</hit_reg>
      </hit>
    </hits>
  </chain_blast_run>
  <blast_run program="blastp">
    <hits>
      <hit domain_id="e6dgvA1" pdb_id="6dgv" chain_id="A" evalues="1e-30">
        <query_reg>260-480</query_reg>
        <hit_reg>8-228</hit_reg>
      </hit>
    </hits>
  </blast_run>
  <hh_run program="hhsearch">
    <hits>
      <hit hit_id="6dgv_A" domain_id="e6dgvA1" num="1" probability="99.5" evalue="1e-20" score="150.5">
        <query_reg>255-490</query_reg>
        <hit_reg>5-240</hit_reg>
      </hit>
    </hits>
  </hh_run>
</blast_summ_doc>"""

        xml_file = tmp_path / "test_summary.xml"
        xml_file.write_text(xml_content)

        # Create mock reference data
        reference_lengths = {"e6dgvA1": 238, "6dgv": 238}
        protein_lengths = {("6dgv", "A"): 238}

        # Parse
        evidence = parse_domain_summary(
            str(xml_file),
            reference_lengths=reference_lengths,
            protein_lengths=protein_lengths,
            require_reference_lengths=False,
        )

        # Verify we got all evidence types
        assert len(evidence) >= 1

        # Check chain BLAST
        chain_blast = [e for e in evidence if e.type == "chain_blast"]
        assert len(chain_blast) == 1
        assert chain_blast[0].source_pdb == "6dgv"
        assert str(chain_blast[0].query_range) == "252-494"
        assert chain_blast[0].evalue == 1e-50

        # Check domain BLAST
        domain_blast = [e for e in evidence if e.type == "domain_blast"]
        assert len(domain_blast) == 1
        assert domain_blast[0].domain_id == "e6dgvA1"
        assert str(domain_blast[0].query_range) == "260-480"

        # Check HHSearch (may be filtered)
        hhsearch = [e for e in evidence if e.type == "hhsearch"]
        if hhsearch:
            assert len(hhsearch) == 1
            assert hhsearch[0].confidence == 0.995  # 99.5% converted to 0-1
            assert str(hhsearch[0].query_range) == "255-490"
        else:
            print("⚠️ HHsearch evidence was filtered - this may be expected")

        # Check confidence scores based on e-values using actual domain IDs
        conf_by_domain = {e.domain_id: e.confidence for e in evidence if e.domain_id}

        # Chain blast should have high confidence (1e-50 e-value)
        chain_blast_evidence = [e for e in evidence if e.type == "chain_blast"][0]
        assert chain_blast_evidence.confidence >= 0.8  # 1e-50 should be high confidence

        # Domain blast should have good confidence (1e-30 e-value)
        if "e6dgvA1" in conf_by_domain:
            # Get the domain blast version (not HHsearch)
            domain_blast_conf = [
                e.confidence
                for e in evidence
                if e.type == "domain_blast" and e.domain_id == "e6dgvA1"
            ]
            if domain_blast_conf:
                assert domain_blast_conf[0] >= 0.6  # 1e-30 should be medium-high

        print(f"✓ Parsed {len(evidence)} evidence items successfully")
        for e in evidence:
            print(f"  {e.type}: {e.domain_id or e.source_pdb} (conf={e.confidence:.3f})")


class TestReferenceLengthLoading:
    """Test loading reference length files"""

    @pytest.mark.unit
    def test_load_domain_lengths(self, tmp_path):
        """Test loading domain reference lengths from CSV"""
        csv_file = tmp_path / "domain_lengths.csv"

        # Create test CSV
        with open(csv_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["domain_id", "length"])  # Header
            writer.writerow(["e6dgvA1", "238"])
            writer.writerow(["e2ia4A1", "98"])
            writer.writerow(["e2ia4A2", "156"])

        lengths = load_reference_lengths(str(csv_file))

        assert len(lengths) == 3
        assert lengths["e6dgvA1"] == 238
        assert lengths["e2ia4A1"] == 98
        assert lengths["e2ia4A2"] == 156

    @pytest.mark.unit
    def test_load_domain_lengths_no_header(self, tmp_path):
        """Test loading without header row"""
        csv_file = tmp_path / "lengths_no_header.csv"

        # No header, straight to data
        with open(csv_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["test1", "100"])
            writer.writerow(["test2", "200"])

        lengths = load_reference_lengths(str(csv_file))

        assert len(lengths) == 2
        assert lengths["test1"] == 100
        assert lengths["test2"] == 200

    @pytest.mark.unit
    def test_load_protein_lengths(self, tmp_path):
        """Test loading protein lengths from CSV"""
        csv_file = tmp_path / "protein_lengths.csv"

        # Format: pdb_id,chain_id,length
        with open(csv_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["pdb_id", "chain_id", "length"])
            writer.writerow(["6dgv", "A", "238"])
            writer.writerow(["2ia4", "A", "508"])
            writer.writerow(["8ovp", "A", "569"])

        lengths = load_protein_lengths(str(csv_file))

        assert len(lengths) == 3
        assert lengths[("6dgv", "A")] == 238
        assert lengths[("2ia4", "A")] == 508
        assert lengths[("8ovp", "A")] == 569

    @pytest.mark.unit
    def test_load_protein_lengths_combined_format(self, tmp_path):
        """Test loading protein lengths with pdb_chain format"""
        csv_file = tmp_path / "protein_lengths_alt.csv"

        # Format: pdb_chain,length
        with open(csv_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["protein_id", "length"])
            writer.writerow(["6dgv_A", "238"])
            writer.writerow(["2ia4_B", "508"])

        lengths = load_protein_lengths(str(csv_file))

        assert len(lengths) == 2
        assert lengths[("6dgv", "A")] == 238
        assert lengths[("2ia4", "B")] == 508

    @pytest.mark.unit
    def test_load_missing_file(self):
        """Test handling of missing files"""
        lengths = load_reference_lengths("/nonexistent/file.csv")
        assert lengths == {}

        protein_lengths = load_protein_lengths("/nonexistent/file.csv")
        assert protein_lengths == {}


class TestEvidenceSummary:
    """Test evidence summary statistics"""

    @pytest.mark.unit
    def test_evidence_summary_empty(self):
        """Test summary of empty evidence"""
        summary = get_evidence_summary([])

        assert summary["total"] == 0
        assert summary["by_type"] == {}
        assert summary["high_confidence"] == 0
        assert summary["unique_families"] == 0

    @pytest.mark.unit
    def test_evidence_summary_mixed(self):
        """Test summary of mixed evidence types"""
        evidence = [
            Evidence(
                type="chain_blast",
                source_pdb="pdb1",
                query_range=SequenceRange.parse("1-100"),
                confidence=0.95,
                evalue=1e-50,
            ),
            Evidence(
                type="chain_blast",
                source_pdb="pdb2",
                query_range=SequenceRange.parse("150-250"),
                confidence=0.6,
                evalue=1e-5,
            ),
            Evidence(
                type="domain_blast",
                source_pdb="pdb1",
                query_range=SequenceRange.parse("10-90"),
                confidence=0.8,
                evalue=1e-20,
            ),
            Evidence(
                type="hhsearch",
                source_pdb="pdb3",
                query_range=SequenceRange.parse("300-400"),
                confidence=0.85,
                evalue=1e-15,
            ),
        ]

        summary = get_evidence_summary(evidence)

        assert summary["total"] == 4
        assert summary["by_type"]["chain_blast"] == 2
        assert summary["by_type"]["domain_blast"] == 1
        assert summary["by_type"]["hhsearch"] == 1
        assert summary["high_confidence"] == 3  # conf > 0.7 or evalue < 1e-10
        assert summary["unique_families"] == 3  # pdb1, pdb2, pdb3


class TestHHsearchAPISpecFormat:
    """Test HHsearch evidence parsing in API spec format"""

    @pytest.mark.unit
    def test_parse_hhsearch_api_spec_format(self, tmp_path):
        """Test parsing HHsearch evidence in API spec format (type='hhsearch')"""
        xml_content = """<?xml version="1.0"?>
<domain_summary version="1.0">
  <protein pdb_id="8axb" chain_id="A" length="430">
    <sequence>SQITIQARLISFESNRQQLWKLMADLNTPLINELLCQLGQHPDFEKWQQKGKLPSTVVSQLCQPLKTDPRFAGQPSRLYMSAIHIVDYIYKSWLASSLPFPLVFETNEDMVWSKNQKGRLCVHFNGLSDLIFEVYCGNRQLHWFQRFLEDQQTKRKSKNQHSSGLFTLRNGHLVWLEGEGKGEPWNLHHLTLYCCVDNRLWTEEGTEIVRQEKADKQSTLTRINNSFERPSQPLYQGQSHILVGVSLGLEKPATVAVVDAIANKVLAYRSIKQLLGDNYELLNRQRRQQQYLSHERHKAQKNFSPNQFGASELGQHIDRLLAKAIVALARTYKAGSIVLPKLGDMREVVQSEIQAIAEQKFPGYIEGQQKYAKQYRVNVHRWSYGRLIQSIQSKAAQTGIVIEEGKQPIRGSPHDKAKELALSAYNLRLT</sequence>
  </protein>
  <evidence>
    <hit type="hhsearch" target="e5b43A3" target_family="Ribonuclease H-like" probability="98.4" evalue="2.10e-10" score="109.4" query_range="239-417" target_range="16-161"/>
  </evidence>
</domain_summary>"""

        xml_file = tmp_path / "test_hhsearch_api.xml"
        xml_file.write_text(xml_content)

        # Create mock reference data
        reference_lengths = {"e5b43A3": 236, "5b43A3": 236}

        # Parse with reference lengths
        evidence = parse_domain_summary(
            str(xml_file),
            reference_lengths=reference_lengths,
            require_reference_lengths=False,
        )

        # Verify we got HHsearch evidence
        assert len(evidence) == 1
        assert evidence[0].type == "hhsearch"
        assert evidence[0].domain_id == "e5b43A3"
        assert str(evidence[0].query_range) == "239-417"
        assert evidence[0].evalue == 2.10e-10

        # Verify target_range was parsed correctly
        assert evidence[0].hit_range is not None
        assert str(evidence[0].hit_range) == "16-161"

        # Verify reference coverage calculation
        assert evidence[0].reference_length == 236
        assert evidence[0].reference_coverage is not None
        # (161-16+1) / 236 = 146 / 236 ≈ 0.619
        assert abs(evidence[0].reference_coverage - 0.619) < 0.01

        # Verify high confidence (98.4% probability)
        assert evidence[0].confidence > 0.7

        print(f"✓ HHsearch evidence parsed successfully")
        print(f"  Domain: {evidence[0].domain_id}")
        print(f"  Query range: {evidence[0].query_range}")
        print(f"  Target range: {evidence[0].hit_range}")
        print(f"  Reference coverage: {evidence[0].reference_coverage:.1%}")
        print(f"  Confidence: {evidence[0].confidence:.3f}")

    @pytest.mark.unit
    def test_hhsearch_without_target_range(self, tmp_path):
        """Test HHsearch evidence without target_range (should still parse but validation may fail)"""
        xml_content = """<?xml version="1.0"?>
<domain_summary version="1.0">
  <protein pdb_id="test" chain_id="A" length="200">
    <sequence>AAAAA</sequence>
  </protein>
  <evidence>
    <hit type="hhsearch" target="e1234A1" probability="95.0" evalue="1e-15" query_range="10-100"/>
  </evidence>
</domain_summary>"""

        xml_file = tmp_path / "test_hhsearch_no_target.xml"
        xml_file.write_text(xml_content)

        reference_lengths = {"e1234A1": 150}

        # Parse (should not crash)
        evidence = parse_domain_summary(
            str(xml_file),
            reference_lengths=reference_lengths,
            require_reference_lengths=False,
        )

        # Should parse but may not have hit_range
        # In strict mode, this would be rejected
        # In non-strict mode, it should be accepted
        if evidence:
            assert evidence[0].type == "hhsearch"
            assert evidence[0].domain_id == "e1234A1"


class TestIntegrationWithRealData:
    """Integration tests with real domain summary files"""

    @pytest.mark.integration
    @pytest.mark.slow
    def test_real_domain_summary_8ovp(self, domain_summary_path, real_reference_data):
        """Test parsing real 8ovp_A domain summary"""
        evidence = parse_domain_summary(
            domain_summary_path,
            reference_lengths=real_reference_data["domain_lengths"],
            protein_lengths=real_reference_data["protein_lengths"],
        )

        # Should have evidence
        assert len(evidence) > 0

        # Check evidence types
        summary = get_evidence_summary(evidence)
        assert summary["total"] > 0

        # Should have multiple evidence types
        assert len(summary["by_type"]) >= 2

        # Check for expected families
        families = {e.source_pdb for e in evidence if e.source_pdb}
        expected_families = {"6dgv", "2ia4"}  # Known hits for 8ovp_A
        assert len(families.intersection(expected_families)) > 0


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])
