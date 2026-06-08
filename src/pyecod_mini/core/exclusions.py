#!/usr/bin/env python3
"""
Evidence exclusion (self / F-group / T-group / explicit domain id) for
non-circular validation of existing ECOD representatives.

Motivation
----------
Every ECOD representative domain is itself in the reference library, so a naive
``partition_protein`` run is *circular*: the query self-matches its own reference
entry and the algorithm trivially recovers the existing label. To validate an
existing rep we must classify it from *independent* evidence only.

This module provides a small, DB-free filter that drops evidence items matching
an :class:`ExclusionPolicy` before partitioning. It operates purely on fields the
evidence already carries (``source_pdb``, ``domain_id``, ``t_group``, ``f_group``),
consistent with PYECOD_MINI_API_SPEC.md (no batch/DB dependency).

The default behavior of the partitioner is unchanged — exclusion is opt-in.
"""

from dataclasses import dataclass
from typing import Iterable, List, Optional, Tuple

from .models import Evidence


@dataclass(frozen=True)
class ExclusionPolicy:
    """Describes which reference evidence to mask before partitioning.

    Attributes:
        exclude_self: Drop hits to the query's own structure (same PDB id). This
            removes the trivial self-hit that makes rep validation circular.
        exclude_domain_ids: Drop hits whose reference ECOD domain id is in this set
            (leading ``e`` optional; matched with and without it).
        exclude_fgroups: Drop hits whose ``f_group`` is in this set. Requires the
            evidence to carry ``f_group`` (emitted by pyecod_prod summaries that
            include classification attributes).
        exclude_tgroups: Drop hits whose ``t_group`` is in this set.
    """

    exclude_self: bool = False
    exclude_domain_ids: frozenset = frozenset()
    exclude_fgroups: frozenset = frozenset()
    exclude_tgroups: frozenset = frozenset()

    @property
    def is_active(self) -> bool:
        """True if any exclusion is configured."""
        return bool(
            self.exclude_self
            or self.exclude_domain_ids
            or self.exclude_fgroups
            or self.exclude_tgroups
        )

    def describe(self) -> str:
        """Compact, human/machine-readable description for output metadata."""
        parts: List[str] = []
        if self.exclude_self:
            parts.append("self")
        if self.exclude_domain_ids:
            parts.append(f"domains:{len(self.exclude_domain_ids)}")
        if self.exclude_fgroups:
            parts.append(f"fgroups:{len(self.exclude_fgroups)}")
        if self.exclude_tgroups:
            parts.append(f"tgroups:{len(self.exclude_tgroups)}")
        return ",".join(parts) if parts else "none"


def _normalize_domain_id(domain_id: Optional[str]) -> Optional[str]:
    """Strip a leading ``e`` so 'e1gcyA2' and '1gcyA2' compare equal."""
    if not domain_id:
        return None
    d = domain_id.strip()
    return d[1:] if d.startswith("e") else d


def _expand_domain_ids(domain_ids: Iterable[str]) -> frozenset:
    """Build a match set containing each id both with and without a leading ``e``."""
    expanded = set()
    for d in domain_ids:
        if not d:
            continue
        d = d.strip()
        if not d:
            continue
        expanded.add(d)
        normalized = _normalize_domain_id(d)
        if normalized:
            expanded.add(normalized)
    return frozenset(expanded)


def exclusion_reason(
    evidence: Evidence,
    query_pdb: str,
    expanded_domain_ids: frozenset,
    policy: ExclusionPolicy,
) -> Optional[str]:
    """Return the reason this evidence item is excluded, or None if kept.

    Reasons (checked in priority order): 'self', 'domain_id', 'f_group', 't_group'.
    """
    if policy.exclude_self and evidence.source_pdb:
        if evidence.source_pdb.lower() == query_pdb:
            return "self"

    if expanded_domain_ids and evidence.domain_id:
        if (
            evidence.domain_id in expanded_domain_ids
            or _normalize_domain_id(evidence.domain_id) in expanded_domain_ids
        ):
            return "domain_id"

    if policy.exclude_fgroups:
        f_group = getattr(evidence, "f_group", None)
        if f_group and f_group in policy.exclude_fgroups:
            return "f_group"

    if policy.exclude_tgroups and evidence.t_group:
        if evidence.t_group in policy.exclude_tgroups:
            return "t_group"

    return None


def apply_exclusions(
    evidence: List[Evidence],
    query_pdb: str,
    query_chain: str,
    policy: Optional[ExclusionPolicy],
) -> Tuple[List[Evidence], List[Evidence]]:
    """Split evidence into (kept, masked) according to ``policy``.

    Args:
        evidence: Parsed evidence items.
        query_pdb: PDB id of the query (used for ``exclude_self``).
        query_chain: Chain id of the query (reserved; self-exclusion is PDB-level
            because all chains of the deposited structure are in the reference).
        policy: The exclusion policy, or None.

    Returns:
        (kept, masked) lists. If ``policy`` is None/inactive, kept == evidence and
        masked is empty.
    """
    if not policy or not policy.is_active:
        return list(evidence), []

    query_pdb_lower = (query_pdb or "").lower()
    expanded_domain_ids = _expand_domain_ids(policy.exclude_domain_ids)

    kept: List[Evidence] = []
    masked: List[Evidence] = []
    for ev in evidence:
        reason = exclusion_reason(ev, query_pdb_lower, expanded_domain_ids, policy)
        if reason:
            masked.append(ev)
        else:
            kept.append(ev)
    return kept, masked


def mark_top_evidence_masked(domains, masked_evidence: List[Evidence]) -> int:
    """Flag domains whose range overlaps any masked evidence.

    Sets ``domain.top_evidence_masked = True`` on each domain that overlaps a
    masked hit. This is the per-domain note requested in the FR: it tells a
    reviewer that the domain's region had independent (e.g. self) evidence removed,
    so the surviving assignment came from other evidence.

    Returns the number of domains flagged.
    """
    if not masked_evidence:
        return 0

    masked_positions = [ev.get_positions() for ev in masked_evidence]
    flagged = 0
    for domain in domains:
        domain_positions = domain.get_positions()
        for mpos in masked_positions:
            if domain_positions & mpos:
                domain.top_evidence_masked = True
                flagged += 1
                break
    return flagged


def load_id_list(path: str) -> frozenset:
    """Load a newline-delimited list of ids from a file.

    Blank lines and lines starting with '#' are ignored. Used by the CLI for
    ``--exclude-domains``, ``--exclude-fgroups`` and ``--exclude-tgroups``.
    """
    ids = set()
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith("#"):
                ids.add(s)
    return frozenset(ids)
