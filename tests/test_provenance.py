#!/usr/bin/env python3
"""
Comprehensive tests for provenance tracking features.

Tests evidence tracking, domain optimization tracking, and metadata management.
"""

import os
from datetime import datetime

import pytest

from pyecod_mini.core.models import (
    Domain,
    DomainLayout,
    Evidence,
    PartitionMetadata,
    UnassignedSegment,
)
from pyecod_mini.core.sequence_range import SequenceRange


@pytest.mark.unit
class TestEvidenceProvenance:
    """Test Evidence provenance tracking methods"""

    def test_evidence_get_reference_coverage(self):
        """Test get_reference_coverage() calculation"""
        evidence = Evidence(
            type="domain_blast",
            source_pdb="6dgv",
            query_range=SequenceRange.parse("10-110"),
            hit_range=SequenceRange.parse("1-100"),
            reference_length=200,
        )

        coverage = evidence.get_reference_coverage()

        # hit_range length is 100, reference length is 200
        assert coverage == 0.5

    def test_evidence_get_reference_coverage_no_data(self):
        """Test get_reference_coverage() returns None without data"""
        evidence = Evidence(
            type="domain_blast",
            source_pdb="6dgv",
            query_range=SequenceRange.parse("10-110"),
        )

        coverage = evidence.get_reference_coverage()

        assert coverage is None

    def test_evidence_get_reference_coverage_cached(self):
        """Test get_reference_coverage() uses cached value"""
        evidence = Evidence(
            type="domain_blast",
            source_pdb="6dgv",
            query_range=SequenceRange.parse("10-110"),
            hit_range=SequenceRange.parse("1-100"),
            reference_length=200,
            reference_coverage=0.75,  # Pre-set value
        )

        coverage = evidence.get_reference_coverage()

        # Should use cached value, not calculate
        assert coverage == 0.75

    def test_evidence_get_quality_metrics(self):
        """Test get_quality_metrics() returns comprehensive metrics"""
        evidence = Evidence(
            type="domain_blast",
            source_pdb="6dgv",
            query_range=SequenceRange.parse("10-110"),
            hit_range=SequenceRange.parse("1-100"),
            reference_length=200,
            confidence=0.85,
            evalue=1e-50,
            discontinuous=False,
        )

        metrics = evidence.get_quality_metrics()

        assert metrics["confidence"] == 0.85
        assert metrics["evalue"] == 1e-50
        assert metrics["reference_coverage"] == 0.5
        assert metrics["query_length"] == 101
        assert metrics["hit_length"] == 100
        assert metrics["reference_length"] == 200
        assert metrics["discontinuous"] is False

    def test_evidence_get_quality_metrics_hhsearch(self):
        """Test get_quality_metrics() includes HHsearch probability"""
        evidence = Evidence(
            type="hhsearch",
            source_pdb="6dgv",
            query_range=SequenceRange.parse("10-110"),
            confidence=0.90,  # confidence = original_prob / 100 * 0.95
        )

        metrics = evidence.get_quality_metrics()

        # Should reverse calculate original probability
        assert "original_hhsearch_probability" in metrics
        # 0.90 / 0.95 * 100 ≈ 94.7
        assert 94 <= metrics["original_hhsearch_probability"] <= 95

    def test_evidence_to_provenance_dict(self):
        """Test to_provenance_dict() export"""
        evidence = Evidence(
            type="domain_blast",
            source_pdb="6dgv",
            source_chain_id="A",
            domain_id="e6dgvA1",
            query_range=SequenceRange.parse("10-110"),
            hit_range=SequenceRange.parse("1-100"),
            reference_length=200,
            confidence=0.85,
            evalue=1e-50,
            hsp_count=3,
            discontinuous=True,
        )

        prov_dict = evidence.to_provenance_dict()

        assert prov_dict["source_type"] == "domain_blast"
        assert prov_dict["source_id"] == "6dgv_A"
        assert prov_dict["domain_id"] == "e6dgvA1"
        assert prov_dict["evalue"] == 1e-50
        assert prov_dict["confidence"] == 0.85
        assert prov_dict["query_range"] == "10-110"
        assert prov_dict["hit_range"] == "1-100"
        assert prov_dict["reference_coverage"] == 0.5
        assert prov_dict["hsp_count"] == 3
        assert prov_dict["discontinuous"] is True


@pytest.mark.unit
class TestDomainProvenanceTracking:
    """Test Domain provenance and optimization tracking"""

    def test_domain_record_optimization_action(self):
        """Test record_optimization_action() adds actions"""
        domain = Domain(
            id="e8ovpA1",
            range=SequenceRange.parse("10-110"),
            family="e6dgvA1",
            source="chain_blast",
            evidence_count=1,
            evidence_items=[],
        )

        # Initially no optimization actions
        assert domain.optimization_actions == []

        # Record some actions
        domain.record_optimization_action("merge_segment", "nterm_1-5")
        domain.record_optimization_action("overlap_trim", "lost_3_to_e8ovpA2")

        assert len(domain.optimization_actions) == 2
        assert "merge_segment:nterm_1-5" in domain.optimization_actions
        assert "overlap_trim:lost_3_to_e8ovpA2" in domain.optimization_actions

    def test_domain_record_optimization_action_no_details(self):
        """Test record_optimization_action() without details"""
        domain = Domain(
            id="e8ovpA1",
            range=SequenceRange.parse("10-110"),
            family="e6dgvA1",
            source="chain_blast",
            evidence_count=1,
            evidence_items=[],
        )

        domain.record_optimization_action("boundary_adjusted")

        assert "boundary_adjusted" in domain.optimization_actions

    def test_domain_was_optimized_true(self):
        """Test was_optimized() returns True when range changed"""
        domain = Domain(
            id="e8ovpA1",
            range=SequenceRange.parse("10-110"),
            family="e6dgvA1",
            source="chain_blast",
            evidence_count=1,
            evidence_items=[],
        )

        # Record original range
        original_range = domain.original_range

        # Modify the range
        domain.range = SequenceRange.parse("10-115")

        assert domain.was_optimized() is True

    def test_domain_was_optimized_false(self):
        """Test was_optimized() returns False when range unchanged"""
        domain = Domain(
            id="e8ovpA1",
            range=SequenceRange.parse("10-110"),
            family="e6dgvA1",
            source="chain_blast",
            evidence_count=1,
            evidence_items=[],
        )

        # Range not modified
        assert domain.was_optimized() is False

    def test_domain_get_optimization_summary(self):
        """Test get_optimization_summary() returns complete summary"""
        evidence = Evidence(
            type="chain_blast",
            source_pdb="6dgv",
            query_range=SequenceRange.parse("10-110"),
            confidence=0.85,
        )

        domain = Domain(
            id="e8ovpA1",
            range=SequenceRange.parse("10-110"),
            family="e6dgvA1",
            source="chain_blast",
            evidence_count=1,
            evidence_items=[evidence],
        )

        # Optimize the domain
        domain.record_optimization_action("merge_segment", "nterm_1-5")
        domain.range = SequenceRange.parse("1-110")
        domain.assigned_positions = set(range(1, 111))  # Update assigned positions

        summary = domain.get_optimization_summary()

        assert summary["was_optimized"] is True
        assert summary["original_range"] == "10-110"
        assert summary["final_range"] == "1-110"
        assert len(summary["actions"]) == 1
        assert "merge_segment:nterm_1-5" in summary["actions"]
        # Position change: new range has 110 positions, original had 101
        assert summary["position_change"] == 110 - 101

    def test_domain_get_quality_assessment(self):
        """Test get_quality_assessment() comprehensive quality check"""
        evidence = Evidence(
            type="domain_blast",
            source_pdb="6dgv",
            query_range=SequenceRange.parse("10-110"),
            hit_range=SequenceRange.parse("1-100"),
            reference_length=200,
            confidence=0.85,
            evalue=1e-50,
        )

        domain = Domain(
            id="e8ovpA1",
            range=SequenceRange.parse("10-110"),
            family="e6dgvA1",
            source="chain_blast",
            evidence_count=1,
            evidence_items=[evidence],
        )

        assessment = domain.get_quality_assessment()

        assert assessment["domain_id"] == "e8ovpA1"
        assert assessment["primary_evidence_quality"] is not None
        assert len(assessment["all_evidence_quality"]) == 1
        assert assessment["overall_assessment"] == "good"
        assert assessment["quality_issues"] == []

    def test_domain_get_quality_assessment_poor(self):
        """Test get_quality_assessment() identifies poor quality"""
        evidence = Evidence(
            type="domain_blast",
            source_pdb="6dgv",
            query_range=SequenceRange.parse("10-110"),
            hit_range=SequenceRange.parse("1-50"),
            reference_length=200,
            confidence=0.3,  # Low confidence
            evalue=5.0,  # Poor e-value
        )

        domain = Domain(
            id="e8ovpA1",
            range=SequenceRange.parse("10-110"),
            family="e6dgvA1",
            source="chain_blast",
            evidence_count=1,
            evidence_items=[evidence],
        )

        assessment = domain.get_quality_assessment()

        assert assessment["overall_assessment"] == "poor"
        assert "low_confidence" in assessment["quality_issues"]
        assert "poor_reference_coverage" in assessment["quality_issues"]
        assert "poor_evalue" in assessment["quality_issues"]

    def test_domain_to_provenance_dict(self):
        """Test to_provenance_dict() exports complete provenance"""
        evidence = Evidence(
            type="domain_blast",
            source_pdb="6dgv",
            source_chain_id="A",
            domain_id="e6dgvA1",
            query_range=SequenceRange.parse("10-110"),
            confidence=0.85,
        )

        domain = Domain(
            id="e8ovpA1",
            range=SequenceRange.parse("10-110"),
            family="e6dgvA1",
            source="chain_blast",
            evidence_count=1,
            evidence_items=[evidence],
            t_group="001",
            h_group="001.001",
            x_group="001.001.001",
            reference_ecod_domain_id="e6dgvA1",
        )

        # Add optimization
        domain.record_optimization_action("merge_segment")
        domain.range = SequenceRange.parse("5-110")

        prov_dict = domain.to_provenance_dict()

        assert prov_dict["id"] == "e8ovpA1"
        assert prov_dict["range"] == "5-110"
        assert prov_dict["family"] == "e6dgvA1"
        assert prov_dict["source"] == "chain_blast"
        assert prov_dict["evidence_count"] == 1
        assert prov_dict["confidence"] == 0.85
        assert prov_dict["t_group"] == "001"
        assert prov_dict["h_group"] == "001.001"
        assert prov_dict["x_group"] == "001.001.001"
        assert prov_dict["reference_ecod_domain_id"] == "e6dgvA1"
        assert prov_dict["was_optimized"] is True
        assert prov_dict["original_range"] == "10-110"
        assert len(prov_dict["optimization_actions"]) == 1
        assert "primary_evidence" in prov_dict


@pytest.mark.unit
class TestDomainLayoutProvenanceTracking:
    """Test DomainLayout provenance tracking during optimization"""

    def test_layout_merge_segment_with_domain_tracking(self):
        """Test merge_segment_with_domain() records provenance"""
        evidence = Evidence(
            type="chain_blast",
            source_pdb="6dgv",
            query_range=SequenceRange.parse("10-110"),
            confidence=0.85,
        )

        domain = Domain(
            id="e8ovpA1",
            range=SequenceRange.parse("10-110"),
            family="e6dgvA1",
            source="chain_blast",
            evidence_count=1,
            evidence_items=[evidence],
        )

        layout = DomainLayout.from_domains([domain], sequence_length=150)

        # Find an unassigned segment
        layout.analyze_gaps()
        assert len(layout.unassigned_segments) > 0

        segment = layout.unassigned_segments[0]
        original_actions_count = len(domain.optimization_actions)

        # Merge segment
        layout.merge_segment_with_domain(segment, domain)

        # Should have recorded the optimization action
        assert len(domain.optimization_actions) > original_actions_count
        assert any("merge_segment" in action for action in domain.optimization_actions)

    def test_layout_split_segment_tracking(self):
        """Test split_segment_between_domains() records provenance"""
        evidence1 = Evidence(
            type="chain_blast",
            source_pdb="6dgv",
            query_range=SequenceRange.parse("10-50"),
            confidence=0.85,
        )
        evidence2 = Evidence(
            type="chain_blast",
            source_pdb="2ia4",
            query_range=SequenceRange.parse("70-110"),
            confidence=0.80,
        )

        domain1 = Domain(
            id="e8ovpA1",
            range=SequenceRange.parse("10-50"),
            family="e6dgvA1",
            source="chain_blast",
            evidence_count=1,
            evidence_items=[evidence1],
        )

        domain2 = Domain(
            id="e8ovpA2",
            range=SequenceRange.parse("70-110"),
            family="e2ia4A1",
            source="chain_blast",
            evidence_count=1,
            evidence_items=[evidence2],
        )

        layout = DomainLayout.from_domains([domain1, domain2], sequence_length=150)
        layout.analyze_gaps()

        # Find segment between domains
        inter_segments = [
            s
            for s in layout.unassigned_segments
            if s.start > domain1.end_position and s.end < domain2.start_position
        ]

        if inter_segments:
            segment = inter_segments[0]

            # Split segment between domains
            split_pos = (segment.start + segment.end) // 2
            positions_to_d1 = set(range(segment.start, split_pos + 1))
            positions_to_d2 = set(range(split_pos + 1, segment.end + 1))

            layout.split_segment_between_domains(
                segment, domain1, domain2, (positions_to_d1, positions_to_d2)
            )

            # Both domains should have recorded the split
            assert any("split_merge" in action for action in domain1.optimization_actions)
            assert any("split_merge" in action for action in domain2.optimization_actions)

    def test_layout_resolve_overlaps_tracking(self):
        """Test resolve_small_overlaps() records provenance"""
        evidence1 = Evidence(
            type="chain_blast",
            source_pdb="6dgv",
            query_range=SequenceRange.parse("10-60"),
            confidence=0.85,
        )
        evidence2 = Evidence(
            type="chain_blast",
            source_pdb="2ia4",
            query_range=SequenceRange.parse("58-110"),  # Overlaps 58-60
            confidence=0.80,
        )

        domain1 = Domain(
            id="e8ovpA1",
            range=SequenceRange.parse("10-60"),
            family="e6dgvA1",
            source="chain_blast",
            evidence_count=1,
            evidence_items=[evidence1],
        )

        domain2 = Domain(
            id="e8ovpA2",
            range=SequenceRange.parse("58-110"),
            family="e2ia4A1",
            source="chain_blast",
            evidence_count=1,
            evidence_items=[evidence2],
        )

        layout = DomainLayout.from_domains([domain1, domain2], sequence_length=150)

        # Resolve overlaps
        resolved_count = layout.resolve_small_overlaps(max_overlap_size=5)

        # Should have resolved the overlap
        assert resolved_count > 0

        # The domain with lower confidence should have recorded the trim
        assert any("overlap_trim" in action for action in domain2.optimization_actions)


@pytest.mark.unit
class TestPartitionMetadata:
    """Test PartitionMetadata provenance container"""

    def test_partition_metadata_creation(self):
        """Test PartitionMetadata initialization"""
        metadata = PartitionMetadata(pdb_id="8ovp", chain_id="A")

        assert metadata.pdb_id == "8ovp"
        assert metadata.chain_id == "A"
        assert metadata.algorithm_version is None
        assert metadata.processing_timestamp is not None
        assert isinstance(metadata.processing_timestamp, datetime)

    def test_partition_metadata_with_version(self):
        """Test PartitionMetadata with versioning"""
        metadata = PartitionMetadata(
            pdb_id="8ovp",
            chain_id="A",
            algorithm_version="2.0.0",
            git_commit_hash="abc123def456",
        )

        assert metadata.algorithm_version == "2.0.0"
        assert metadata.git_commit_hash == "abc123def456"

    def test_partition_metadata_with_file_provenance(self):
        """Test PartitionMetadata with file provenance"""
        metadata = PartitionMetadata(
            pdb_id="8ovp",
            chain_id="A",
            source_domain_summary_path="/path/to/summary.xml",
            source_domain_summary_hash="sha256hash",
            output_xml_path="/path/to/output.xml",
            batch_id="ecod_batch_036",
        )

        assert metadata.source_domain_summary_path == "/path/to/summary.xml"
        assert metadata.source_domain_summary_hash == "sha256hash"
        assert metadata.output_xml_path == "/path/to/output.xml"
        assert metadata.batch_id == "ecod_batch_036"

    def test_partition_metadata_with_parameters(self):
        """Test PartitionMetadata with process parameters"""
        metadata = PartitionMetadata(
            pdb_id="8ovp",
            chain_id="A",
            sequence_length=569,
            process_parameters={
                "boundary_optimization_enabled": True,
                "min_domain_size": 25,
                "overlap_threshold": 5,
            },
        )

        assert metadata.sequence_length == 569
        assert metadata.process_parameters["boundary_optimization_enabled"] is True
        assert metadata.process_parameters["min_domain_size"] == 25
        assert metadata.process_parameters["overlap_threshold"] == 5


@pytest.mark.integration
class TestProvenanceInXML:
    """Test that provenance is correctly written to XML"""

    def test_provenance_in_partition_xml(self, domain_summary_path, temp_output_dir):
        """Test that provenance data appears in partition.xml"""
        from pyecod_mini import partition_protein

        output_path = os.path.join(temp_output_dir, "provenance_test.xml")

        result = partition_protein(
            summary_xml=domain_summary_path,
            output_xml=output_path,
            pdb_id="8ovp",
            chain_id="A",
            batch_id="test_batch",
        )

        # Read the XML
        import xml.etree.ElementTree as ET

        tree = ET.parse(output_path)
        root = tree.getroot()

        # Check metadata contains provenance information
        metadata = root.find("metadata")
        assert metadata is not None

        # Version tracking
        version_elem = metadata.find("version")
        assert version_elem is not None
        assert version_elem.get("algorithm") is not None
        assert version_elem.get("git_commit") is not None
        assert version_elem.get("timestamp") is not None

        # Source provenance
        source_elem = metadata.find("source")
        if source_elem is not None:
            # If source element exists, should have batch_id
            assert source_elem.get("batch_id") is not None

        # Statistics
        stats_elem = metadata.find("statistics")
        assert stats_elem is not None
        assert stats_elem.get("sequence_length") is not None
        assert stats_elem.get("domain_count") is not None
        assert stats_elem.get("total_coverage") is not None

    def test_domain_provenance_in_xml(self, domain_summary_path, temp_output_dir):
        """Test that domain-level provenance appears in XML"""
        from pyecod_mini import partition_protein

        output_path = os.path.join(temp_output_dir, "domain_provenance_test.xml")

        result = partition_protein(
            summary_xml=domain_summary_path,
            output_xml=output_path,
            pdb_id="8ovp",
            chain_id="A",
        )

        # Read the XML
        import xml.etree.ElementTree as ET

        tree = ET.parse(output_path)
        root = tree.getroot()

        # Find domains
        domains_elem = root.find("domains")
        assert domains_elem is not None

        domain_elems = domains_elem.findall("domain")
        assert len(domain_elems) > 0

        # Check first domain for provenance
        domain = domain_elems[0]

        # Should have basic attributes
        assert domain.get("id") is not None
        assert domain.get("range") is not None
        assert domain.get("family") is not None
        assert domain.get("source") is not None
        assert domain.get("evidence_count") is not None

        # Check for primary evidence
        primary_evidence = domain.find("primary_evidence")
        if primary_evidence is not None:
            # Should have evidence details
            assert primary_evidence.get("source_type") is not None
            assert primary_evidence.get("confidence") is not None

        # Check for optimization tracking
        optimization = domain.find("boundary_optimization")
        if optimization is not None:
            # Should have optimization details
            assert optimization.get("original_range") is not None
            assert optimization.get("optimized_range") is not None
