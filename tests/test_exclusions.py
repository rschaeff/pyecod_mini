#!/usr/bin/env python3
"""
Tests for evidence exclusion (self / F-group / T-group / explicit domain id).

These cover the non-circular validation feature used to validate existing ECOD
representatives without the query trivially self-matching its own reference entry.
"""

from pyecod_mini.core.exclusions import (
    ExclusionPolicy,
    apply_exclusions,
    mark_top_evidence_masked,
)
from pyecod_mini.core.models import Domain, Evidence
from pyecod_mini.core.sequence_range import SequenceRange


def _ev(domain_id, source_pdb, query_range="10-110", t_group=None, f_group=None):
    return Evidence(
        type="domain_blast",
        source_pdb=source_pdb,
        query_range=SequenceRange.parse(query_range),
        domain_id=domain_id,
        evalue=0.0,
        t_group=t_group,
        f_group=f_group,
    )


class TestExclusionPolicy:
    def test_inactive_by_default(self):
        assert ExclusionPolicy().is_active is False

    def test_active_flags(self):
        assert ExclusionPolicy(exclude_self=True).is_active
        assert ExclusionPolicy(exclude_domain_ids=frozenset({"e1gcyA2"})).is_active
        assert ExclusionPolicy(exclude_fgroups=frozenset({"1.1.1.1"})).is_active
        assert ExclusionPolicy(exclude_tgroups=frozenset({"1.1.1"})).is_active

    def test_describe(self):
        p = ExclusionPolicy(exclude_self=True, exclude_domain_ids=frozenset({"a", "b"}))
        d = p.describe()
        assert "self" in d and "domains:2" in d


class TestApplyExclusions:
    def test_default_is_noop(self):
        evidence = [_ev("e1gcyA2", "1gcy"), _ev("e2xyzB1", "2xyz")]
        kept, masked = apply_exclusions(evidence, "1gcy", "A", None)
        assert len(kept) == 2
        assert masked == []

    def test_inactive_policy_is_noop(self):
        evidence = [_ev("e1gcyA2", "1gcy")]
        kept, masked = apply_exclusions(evidence, "1gcy", "A", ExclusionPolicy())
        assert len(kept) == 1 and masked == []

    def test_exclude_self_removes_own_structure(self):
        # The 1gcy_A -> e1gcyA2 self-hit (the circularity case from the FR)
        evidence = [_ev("e1gcyA2", "1gcy"), _ev("e2xyzB1", "2xyz")]
        kept, masked = apply_exclusions(
            evidence, "1gcy", "A", ExclusionPolicy(exclude_self=True)
        )
        assert [e.domain_id for e in kept] == ["e2xyzB1"]
        assert [e.domain_id for e in masked] == ["e1gcyA2"]

    def test_exclude_self_is_case_insensitive(self):
        evidence = [_ev("e1GCYA2", "1GCY")]
        kept, masked = apply_exclusions(
            evidence, "1gcy", "A", ExclusionPolicy(exclude_self=True)
        )
        assert kept == [] and len(masked) == 1

    def test_exclude_domain_ids_with_and_without_e_prefix(self):
        evidence = [_ev("e1gcyA2", "1gcy"), _ev("e2xyzB1", "2xyz")]
        # supply id without leading 'e'
        kept, masked = apply_exclusions(
            evidence, "9zzz", "A",
            ExclusionPolicy(exclude_domain_ids=frozenset({"1gcyA2"})),
        )
        assert [e.domain_id for e in masked] == ["e1gcyA2"]
        assert [e.domain_id for e in kept] == ["e2xyzB1"]

    def test_exclude_fgroups(self):
        evidence = [
            _ev("e1gcyA2", "1gcy", f_group="2008.1.1.1"),
            _ev("e2xyzB1", "2xyz", f_group="3000.1.1.1"),
        ]
        kept, masked = apply_exclusions(
            evidence, "9zzz", "A",
            ExclusionPolicy(exclude_fgroups=frozenset({"2008.1.1.1"})),
        )
        assert [e.domain_id for e in masked] == ["e1gcyA2"]

    def test_exclude_tgroups(self):
        evidence = [
            _ev("e1gcyA2", "1gcy", t_group="2008.1.1"),
            _ev("e2xyzB1", "2xyz", t_group="3000.1.1"),
        ]
        kept, masked = apply_exclusions(
            evidence, "9zzz", "A",
            ExclusionPolicy(exclude_tgroups=frozenset({"2008.1.1"})),
        )
        assert [e.domain_id for e in masked] == ["e1gcyA2"]

    def test_fgroup_exclusion_noop_when_field_absent(self):
        # Evidence with no f_group should never be masked by f-group policy
        evidence = [_ev("e1gcyA2", "1gcy", f_group=None)]
        kept, masked = apply_exclusions(
            evidence, "9zzz", "A",
            ExclusionPolicy(exclude_fgroups=frozenset({"2008.1.1.1"})),
        )
        assert len(kept) == 1 and masked == []


class TestMarkTopEvidenceMasked:
    def _domain(self, range_str):
        return Domain(
            id="d1",
            range=SequenceRange.parse(range_str),
            family="test",
            evidence_count=1,
            source="domain_blast",
            evidence_items=[],
        )

    def test_marks_overlapping_domain(self):
        domain = self._domain("10-110")
        masked = [_ev("e1gcyA2", "1gcy", query_range="20-90")]
        flagged = mark_top_evidence_masked([domain], masked)
        assert flagged == 1
        assert domain.top_evidence_masked is True

    def test_does_not_mark_non_overlapping_domain(self):
        domain = self._domain("200-300")
        masked = [_ev("e1gcyA2", "1gcy", query_range="20-90")]
        flagged = mark_top_evidence_masked([domain], masked)
        assert flagged == 0
        assert domain.top_evidence_masked is False

    def test_empty_masked_is_noop(self):
        domain = self._domain("10-110")
        assert mark_top_evidence_masked([domain], []) == 0
        assert domain.top_evidence_masked is False
