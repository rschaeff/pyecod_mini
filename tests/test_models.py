#!/usr/bin/env python3
"""
Model tests for mini_pyecod

Tests the data model classes.
"""

# Add parent directory to path for imports
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from pyecod_mini.core.models import AlignmentData, Domain, DomainLayout, Evidence
from pyecod_mini.core.sequence_range import SequenceRange
from pyecod_mini.core.domain_utils import get_domain_coverage_stats


class TestEvidenceModel:
    """Test the Evidence data model"""

    @pytest.mark.unit
    def test_evidence_creation_minimal(self):
        """Test creating evidence with minimal fields"""
        evidence = Evidence(
            type="domain_blast",
            source_pdb="1abc",
            query_range=SequenceRange.parse("10-100"),
        )

        assert evidence.type == "domain_blast"
        assert evidence.source_pdb == "1abc"
        assert str(evidence.query_range) == "10-100"
        assert evidence.confidence == 0.0  # Default
        assert evidence.evalue is None
        assert evidence.domain_id is None  # No auto-generation for minimal evidence

    @pytest.mark.unit
    def test_evidence_creation_full(self):
        """Test creating evidence with all fields"""
        evidence = Evidence(
            type="hhsearch",
            source_pdb="2def",
            query_range=SequenceRange.parse("1-50,100-150"),
            confidence=0.95,
            evalue=1e-50,
            domain_id="e2defA1",
            t_group="1234.5.6",
            h_group="1234.5",
            reference_length=101,
            alignment_coverage=0.85,
            alignment=AlignmentData(
                query_seq="ACDEF",
                hit_seq="ACDEF",
                query_start=1,
                query_end=50,
                hit_start=10,
                hit_end=60,
            ),
        )

        assert evidence.type == "hhsearch"
        assert evidence.confidence == 0.95
        assert evidence.evalue == 1e-50
        assert evidence.t_group == "1234.5.6"
        assert evidence.reference_length == 101
        assert evidence.alignment_coverage == 0.85
        assert evidence.alignment is not None
        assert evidence.alignment.query_seq == "ACDEF"

    @pytest.mark.unit
    def test_evidence_with_discontinuous_range(self):
        """Test evidence with discontinuous range - FIXED calculation"""
        evidence = Evidence(
            type="chain_blast",
            source_pdb="3ghi",
            query_range=SequenceRange.parse("10-50,60-100,150-200"),
            domain_id="e3ghiA1",  # Proper domain_id example
        )

        assert evidence.query_range.is_discontinuous
        assert len(evidence.query_range.segments) == 3
        # FIXED: (50-10+1) + (100-60+1) + (200-150+1) = 41 + 41 + 51 = 133
        assert evidence.query_range.total_length == 133


class TestDomainModel:
    """Test the Domain data model"""

    @pytest.mark.unit
    def test_domain_creation(self):
        """Test creating a domain"""
        evidence_items = [
            Evidence(
                type="domain_blast",
                source_pdb="test",
                query_range=SequenceRange.parse("1-100"),
                confidence=0.9,
                domain_id="eTestA1",  # Proper domain_id
            )
        ]

        domain = Domain(
            id="d1",
            range=SequenceRange.parse("1-100"),
            family="test_family",
            evidence_count=1,
            source="domain_blast",
            evidence_items=evidence_items,
            t_group=None,
            h_group=None,
            x_group=None,
        )

        assert domain.id == "d1"
        assert str(domain.range) == "1-100"
        assert domain.family == "test_family"
        assert domain.evidence_count == 1
        assert domain.source == "domain_blast"
        assert len(domain.evidence_items) == 1

    @pytest.mark.unit
    def test_domain_with_empty_evidence(self):
        """Test domain with no evidence items"""
        domain = Domain(
            id="d2",
            range=SequenceRange.parse("150-250"),
            family="another_family",
            evidence_count=0,
            source="unknown",
            evidence_items=[],
        )

        assert len(domain.evidence_items) == 0
        assert domain.evidence_count == 0

    @pytest.mark.unit
    def test_domain_discontinuous(self):
        """Test discontinuous domain"""
        domain = Domain(
            id="d3",
            range=SequenceRange.parse("1-50,100-150"),
            family="disc_family",
            evidence_count=2,
            source="chain_blast_decomposed",
            evidence_items=[],
        )

        assert domain.range.is_discontinuous
        assert domain.range.total_length == 101


class TestAlignmentDataModel:
    """Test the AlignmentData model"""

    @pytest.mark.unit
    def test_alignment_data_creation(self):
        """Test creating alignment data"""
        alignment = AlignmentData(
            query_seq="ACDEFGHIKLMNPQRSTVWY",
            hit_seq="ACDEFGHIKLMNPQRSTVWY",
            query_start=10,
            query_end=29,
            hit_start=1,
            hit_end=20,
        )

        assert alignment.query_seq == "ACDEFGHIKLMNPQRSTVWY"
        assert alignment.hit_seq == "ACDEFGHIKLMNPQRSTVWY"
        assert alignment.query_start == 10
        assert alignment.query_end == 29
        assert alignment.hit_start == 1
        assert alignment.hit_end == 20

    @pytest.mark.unit
    def test_alignment_with_gaps(self):
        """Test alignment with gaps"""
        alignment = AlignmentData(
            query_seq="ACD-FGH",
            hit_seq="ACDEFGH",
            query_start=1,
            query_end=6,
            hit_start=1,
            hit_end=7,
        )

        assert "-" in alignment.query_seq
        assert "-" not in alignment.hit_seq
        assert len(alignment.query_seq) == len(alignment.hit_seq)


class TestModelRelationships:
    """Test relationships between models"""

    @pytest.mark.unit
    def test_evidence_to_domain_relationship(self):
        """Test that evidence properly relates to domains"""
        # Create evidence
        evidence1 = Evidence(
            type="domain_blast",
            source_pdb="1abc",
            query_range=SequenceRange.parse("10-60"),
            confidence=0.8,
            t_group="1111.1.1",
            domain_id="e1abcA1",
        )

        evidence2 = Evidence(
            type="hhsearch",
            source_pdb="1abc",
            query_range=SequenceRange.parse("15-55"),
            confidence=0.9,
            t_group="1111.1.1",
            domain_id="e1abcA1",
        )

        # Create domain from evidence
        domain = Domain(
            id="d1",
            range=SequenceRange.parse("10-60"),  # Uses evidence1's range
            family=evidence1.t_group or evidence1.source_pdb,
            evidence_count=2,
            source=evidence1.type,
            evidence_items=[evidence1, evidence2],
        )

        # Verify relationships
        assert domain.family == "1111.1.1"
        assert domain.evidence_count == len(domain.evidence_items)
        assert domain.source == "domain_blast"
        assert all(e.t_group == "1111.1.1" for e in domain.evidence_items)

    @pytest.mark.unit
    def test_alignment_attached_to_evidence(self):
        """Test alignment data attached to evidence"""
        alignment = AlignmentData(
            query_seq="ABCDEFG",
            hit_seq="ABCDEFG",
            query_start=1,
            query_end=7,
            hit_start=10,
            hit_end=16,
        )

        evidence = Evidence(
            type="chain_blast",
            source_pdb="2def",
            query_range=SequenceRange.parse("1-7"),
            alignment=alignment,
            domain_id="2def_A",  # Chain-level identifier for chain BLAST
        )

        assert evidence.alignment is not None
        assert evidence.alignment.query_start == 1
        assert evidence.alignment.hit_start == 10


class TestModelDefaults:
    """Test default values and optional fields"""

    @pytest.mark.unit
    def test_evidence_optional_fields(self):
        """Test that optional fields have proper defaults"""
        evidence = Evidence(type="test", source_pdb="test", query_range=SequenceRange.parse("1-10"))

        # Check defaults
        assert evidence.confidence == 0.0
        assert evidence.evalue is None
        assert evidence.domain_id is None
        assert evidence.t_group is None
        assert evidence.h_group is None
        assert evidence.reference_length is None
        assert evidence.alignment_coverage is None
        assert evidence.alignment is None

    @pytest.mark.unit
    def test_domain_evidence_list_initialization(self):
        """Test that evidence_items list is properly initialized"""
        # Domain with explicit empty list
        domain1 = Domain(
            id="d1",
            range=SequenceRange.parse("1-100"),
            family="test",
            evidence_count=0,
            source="test",
            evidence_items=[],
        )

        assert domain1.evidence_items == []
        assert isinstance(domain1.evidence_items, list)

        # Domain with evidence
        evidence = Evidence(
            type="test",
            source_pdb="test",
            query_range=SequenceRange.parse("1-100"),
            domain_id="eTestA1",
        )

        domain2 = Domain(
            id="d2",
            range=SequenceRange.parse("1-100"),
            family="test",
            evidence_count=1,
            source="test",
            evidence_items=[evidence],
        )

        assert len(domain2.evidence_items) == 1
        assert domain2.evidence_items[0] is evidence


class TestOverlappingDomainCoverage:
    """Test coverage calculations with overlapping domains (bug fix verification)

    These tests verify the fix for the overlap bug where domains with
    overlapping residue ranges caused reported coverage to exceed 100%.
    """

    @pytest.mark.unit
    def test_domain_layout_coverage_with_overlapping_domains(self):
        """DomainLayout.get_coverage_stats() should handle overlapping domains correctly."""
        # Create domains with overlap (like 9hpj_B case)
        # d2: 1-255, d1: 236-341 = 20 residue overlap at 236-255
        d2 = Domain(
            id="d2",
            range=SequenceRange.parse("1-255"),
            family="1.1.9",
            evidence_count=1,
            source="chain_blast",
            evidence_items=[],
        )
        d1 = Domain(
            id="d1",
            range=SequenceRange.parse("236-341"),
            family="708.1.2",
            evidence_count=1,
            source="hhsearch",
            evidence_items=[],
        )

        layout = DomainLayout.from_domains([d1, d2], sequence_length=341)
        stats = layout.get_coverage_stats()

        # Coverage should be 100% (341/341), NOT 105.9% (361/341)
        assert stats["coverage_percent"] <= 100.0, f"Coverage {stats['coverage_percent']}% exceeds 100%"
        assert stats["coverage_percent"] == 100.0
        assert stats["assigned_residues"] == 341

    @pytest.mark.unit
    def test_domain_layout_coverage_with_discontinuous_overlap(self):
        """DomainLayout handles discontinuous domains that overlap another domain."""
        # Based on 9uko_F case
        d2 = Domain(
            id="d2",
            range=SequenceRange.parse("1-173"),
            family="142.1.1",
            evidence_count=1,
            source="domain_blast",
            evidence_items=[],
        )
        d1 = Domain(
            id="d1",
            range=SequenceRange.parse("9-25,174-340"),
            family="142.1.1",
            evidence_count=1,
            source="chain_blast_decomposed",
            evidence_items=[],
        )

        layout = DomainLayout.from_domains([d1, d2], sequence_length=340)
        stats = layout.get_coverage_stats()

        # Coverage should be 100% (340/340), NOT 105% (357/340)
        assert stats["coverage_percent"] <= 100.0
        assert stats["coverage_percent"] == 100.0
        assert stats["assigned_residues"] == 340

    @pytest.mark.unit
    def test_get_domain_coverage_stats_with_overlap(self):
        """get_domain_coverage_stats() should handle overlapping domains correctly."""
        # Create overlapping domains
        d1 = Domain(
            id="d1",
            range=SequenceRange.parse("1-100"),
            family="family1",
            evidence_count=1,
            source="domain_blast",
            evidence_items=[],
        )
        d2 = Domain(
            id="d2",
            range=SequenceRange.parse("80-150"),
            family="family2",
            evidence_count=1,
            source="hhsearch",
            evidence_items=[],
        )

        stats = get_domain_coverage_stats([d1, d2], sequence_length=200)

        # Domains cover 1-100 and 80-150, union is 1-150 = 150 residues
        # Coverage should be 150/200 = 75%, NOT (100+71)/200 = 85.5%
        assert stats["total_coverage"] == 150, f"Expected 150 residues, got {stats['total_coverage']}"
        assert stats["coverage_percentage"] == 75.0, f"Expected 75%, got {stats['coverage_percentage']}%"

    @pytest.mark.unit
    def test_get_domain_coverage_stats_non_overlapping(self):
        """get_domain_coverage_stats() works correctly for non-overlapping domains."""
        d1 = Domain(
            id="d1",
            range=SequenceRange.parse("1-100"),
            family="family1",
            evidence_count=1,
            source="domain_blast",
            evidence_items=[],
        )
        d2 = Domain(
            id="d2",
            range=SequenceRange.parse("101-200"),
            family="family2",
            evidence_count=1,
            source="domain_blast",
            evidence_items=[],
        )

        stats = get_domain_coverage_stats([d1, d2], sequence_length=300)

        # 200 residues out of 300
        assert stats["total_coverage"] == 200
        assert abs(stats["coverage_percentage"] - 66.666666) < 0.01

    @pytest.mark.unit
    def test_get_domain_coverage_stats_empty(self):
        """get_domain_coverage_stats() handles empty domain list."""
        stats = get_domain_coverage_stats([], sequence_length=100)

        assert stats["total_domains"] == 0
        assert stats["total_coverage"] == 0
        assert stats["coverage_percentage"] == 0.0


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])
