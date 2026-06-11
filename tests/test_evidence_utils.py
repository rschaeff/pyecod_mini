#!/usr/bin/env python3
"""
Tests for confidence calculation and evidence standardization.

Focus: confidence must be grounded in a real score (e-value or probability).
There is no silent floor/fallback — missing data fails loudly, and a
caller-supplied confidence is never silently overwritten when it cannot be
recomputed.
"""

import pytest

from pyecod_mini.core.evidence_utils import (
    calculate_evidence_confidence,
    populate_evidence_provenance,
)
from pyecod_mini.core.models import Evidence
from pyecod_mini.core.sequence_range import SequenceRange


class TestCalculateEvidenceConfidence:
    """calculate_evidence_confidence contract"""

    @pytest.mark.unit
    def test_evalue_produces_confidence(self):
        conf = calculate_evidence_confidence(evalue=1e-20, evidence_type="domain_blast")
        assert 0.05 <= conf <= 0.95
        assert conf > 0.5  # strong e-value -> high confidence

    @pytest.mark.unit
    def test_probability_produces_confidence(self):
        conf = calculate_evidence_confidence(probability=99.5, evidence_type="hhsearch")
        assert conf > 0.5

    @pytest.mark.unit
    def test_no_score_basis_raises(self):
        """No e-value and no probability -> fail loudly, do not fabricate a floor."""
        with pytest.raises(ValueError, match="no e-value or probability"):
            calculate_evidence_confidence(evidence_type="domain_blast")

    @pytest.mark.unit
    def test_negative_evalue_without_probability_raises(self):
        with pytest.raises(ValueError):
            calculate_evidence_confidence(evalue=-1.0, evidence_type="domain_blast")


class TestPopulateEvidenceProvenanceConfidence:
    """populate_evidence_provenance must not silently overwrite confidence"""

    def _evidence(self, confidence, evalue):
        return Evidence(
            type="domain_blast",
            source_pdb="test",
            query_range=SequenceRange.parse("1-50"),
            confidence=confidence,
            evalue=evalue,
            domain_id="test_A",
        )

    @pytest.mark.unit
    def test_confidence_preserved_when_no_evalue(self):
        """Without an e-value there is nothing to recompute from -> keep caller value."""
        ev = self._evidence(confidence=0.95, evalue=None)
        ev = populate_evidence_provenance(ev)
        assert ev.confidence == 0.95

    @pytest.mark.unit
    def test_confidence_recalculated_when_evalue_present(self):
        """With an e-value, confidence is (re)derived from it."""
        ev = self._evidence(confidence=0.5, evalue=1e-20)
        ev = populate_evidence_provenance(ev)
        assert ev.confidence > 0.5
