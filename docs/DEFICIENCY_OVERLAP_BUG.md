# Deficiency Report: Domain Boundary Overlap Bug

**Date**: 2026-01-16
**Severity**: Medium
**Affected Version**: pyecod_mini 2.0.2
**Reporter**: Partition quality analysis of Q4 2025 + Q1 2026 batch

## Summary

The boundary optimization algorithm in pyecod_mini can produce overlapping domain ranges, causing the reported coverage to exceed 100%. The coverage calculation sums domain lengths rather than computing the union of assigned residue positions.

## Impact

- **Affected chains**: 31 out of 3,677 (0.84%)
- **Average overlap**: 10.2 residues per affected chain
- **Maximum overlap**: 21 residues (9mni_R)
- **Coverage inflation**: Up to 106.2%

## Root Cause

The `_calculate_coverage()` method sums the lengths of all domain ranges without checking for overlaps:

```python
# Current implementation (buggy)
total_assigned = sum(domain.length for domain in domains)
coverage = total_assigned / sequence_length
```

This fails when boundary optimization extends one domain into territory already claimed by another domain.

## Examples

### Example 1: 9mni_R (106.2% coverage)

```
Sequence length: 338 residues
Reported assigned: 359 residues
Overlap: +21 residues

Domain assignments:
  d2: 1-73, 98-99      (family 913.1.1)
  d1: 74-97, 100-140   (family 913.1.1)
  d3: 120-338          (family 5001.1.1)

Overlap detected:
  d1(100-140) ∩ d3(120-338) = 21 residues overlap at positions 120-140
```

### Example 2: 9hpj_B (105.9% coverage)

```
Sequence length: 341 residues
Reported assigned: 361 residues
Overlap: +20 residues

Domain assignments:
  d2: 1-255    (family 1.1.9)
  d1: 236-341  (family 708.1.2)

Overlap detected:
  d2(1-255) ∩ d1(236-341) = 20 residues overlap at positions 236-255
```

### Example 3: 9uko_F (105.0% coverage)

```
Sequence length: 340 residues
Reported assigned: 357 residues
Overlap: +17 residues

Domain assignments:
  d2: 1-173           (family 142.1.1)
  d1: 9-25, 174-340   (family 142.1.1)  <- discontinuous, overlaps with d2

Overlap detected:
  d2(1-173) ∩ d1(9-25) = discontinuous domain inserts into another
```

## Affected Chains (Complete List)

| Chain | Length | Assigned | Coverage | Domains | Overlap |
|-------|-------:|---------:|---------:|--------:|--------:|
| 9mni_R | 338 | 359 | 106.2% | 3 | +21 |
| 9hpj_B | 341 | 361 | 105.9% | 2 | +20 |
| 9n70_SW | 219 | 230 | 105.0% | 2 | +11 |
| 9uko_F | 340 | 357 | 105.0% | 2 | +17 |
| 9eg8_J | 338 | 354 | 104.7% | 4 | +16 |
| 9efm_J | 345 | 361 | 104.6% | 4 | +16 |
| 9eg1_J | 367 | 383 | 104.4% | 4 | +16 |
| 9qf4_Ab | 195 | 203 | 104.1% | 2 | +8 |
| 9nyz_B | 281 | 292 | 103.9% | 2 | +11 |
| 9iz8_A | 217 | 225 | 103.7% | 2 | +8 |
| 9vg4_A | 432 | 444 | 102.8% | 2 | +12 |
| 9v69_A | 272 | 279 | 102.6% | 2 | +7 |
| 9zb7_A | 238 | 244 | 102.5% | 2 | +6 |
| 9mrh_B | 370 | 379 | 102.4% | 3 | +9 |
| 9efm_E | 306 | 313 | 102.3% | 2 | +7 |
| 9hs1_A | 454 | 464 | 102.2% | 2 | +10 |
| 9mtx_B | 400 | 408 | 102.0% | 3 | +8 |
| 9ujd_C | 323 | 329 | 101.9% | 2 | +6 |
| 9jg1_R | 386 | 393 | 101.8% | 2 | +7 |
| 9krv_E | 752 | 765 | 101.7% | 5 | +13 |
| 9n6v_LN | 663 | 674 | 101.7% | 2 | +11 |
| 9i78_LC | 363 | 369 | 101.6% | 2 | +6 |
| 9s3g_Z | 498 | 506 | 101.6% | 7 | +8 |
| 9mkb_WI | 602 | 611 | 101.5% | 4 | +9 |
| 9oou_C | 426 | 432 | 101.4% | 3 | +6 |
| 9lx0_A | 727 | 737 | 101.4% | 4 | +10 |
| 9iqq_D | 506 | 512 | 101.2% | 3 | +6 |
| 9nir_A | 814 | 823 | 101.1% | 5 | +9 |
| 9nq7_A | 758 | 766 | 101.1% | 5 | +8 |
| 9r78_A | 941 | 950 | 101.0% | 2 | +9 |
| 8zdg_B | 684 | 690 | 100.9% | 8 | +6 |

## Proposed Fix

### Option A: Fix Coverage Calculation (Recommended)

Calculate coverage using union of residue positions:

```python
def _calculate_coverage(self, domains: List[Domain], sequence_length: int) -> float:
    """Calculate coverage as union of all assigned positions."""
    assigned_positions = set()

    for domain in domains:
        for segment in domain.segments:
            for pos in range(segment.start, segment.end + 1):
                assigned_positions.add(pos)

    return len(assigned_positions) / sequence_length
```

**Pros**: Accurate coverage, preserves domain assignments
**Cons**: Allows overlapping domains in output (may confuse downstream tools)

### Option B: Prevent Overlaps During Optimization

Modify boundary optimization to respect existing domain boundaries:

```python
def _optimize_boundary(self, domain: Domain, other_domains: List[Domain]) -> Domain:
    """Extend domain boundaries without overlapping other domains."""
    occupied = set()
    for other in other_domains:
        if other.id != domain.id:
            occupied.update(other.get_positions())

    # Only extend into unoccupied positions
    new_segments = []
    for segment in domain.segments:
        start, end = segment.start, segment.end
        while start > 1 and (start - 1) not in occupied:
            start -= 1
            if not self._should_extend(start):
                break
        # ... similar for end
        new_segments.append(Segment(start, end))

    return Domain(domain.id, new_segments, ...)
```

**Pros**: Clean output, no overlaps
**Cons**: May leave small gaps, more complex logic

### Option C: Post-Process to Resolve Overlaps

After optimization, detect and resolve overlaps by assigning disputed residues to the domain with stronger evidence:

```python
def _resolve_overlaps(self, domains: List[Domain]) -> List[Domain]:
    """Resolve overlapping regions by confidence score."""
    position_owners = {}  # pos -> (domain_id, confidence)

    for domain in domains:
        for pos in domain.get_positions():
            if pos not in position_owners or domain.confidence > position_owners[pos][1]:
                position_owners[pos] = (domain.id, domain.confidence)

    # Rebuild domains from position ownership
    # ...
```

**Pros**: Resolves conflicts intelligently
**Cons**: Most complex, may fragment domains

## Recommended Action

1. **Immediate**: Implement Option A (fix coverage calculation) - this correctly reports coverage without changing domain assignments

2. **Follow-up**: Consider Option B or C to prevent overlaps entirely, as overlapping domain assignments may indicate an algorithmic issue in the partitioning logic

## Test Cases

Add unit tests for overlap detection:

```python
def test_coverage_with_overlapping_domains():
    """Coverage should be ≤100% even with overlapping domains."""
    # Create domains with overlap
    d1 = Domain("d1", [Segment(1, 100)])
    d2 = Domain("d2", [Segment(80, 150)])  # Overlaps d1 by 21 residues

    result = partition_protein(domains=[d1, d2], sequence_length=150)

    assert result.coverage <= 1.0, f"Coverage {result.coverage} exceeds 100%"
    assert result.coverage == 1.0, f"Expected 100% coverage (150/150)"

def test_discontinuous_domain_overlap():
    """Discontinuous domains should not create overlaps."""
    # Based on 9uko_F case
    d1 = Domain("d1", [Segment(9, 25), Segment(174, 340)])  # Discontinuous
    d2 = Domain("d2", [Segment(1, 173)])

    result = partition_protein(domains=[d1, d2], sequence_length=340)

    # Positions 9-25 are claimed by both d1 and d2
    assert result.coverage <= 1.0
```

## References

- Batch: `/data/ecod/pdb_updates/batches/ecod_q4_2025_q1_2026/`
- Analysis: `/data/ecod/pdb_updates/batches/ecod_q4_2025_q1_2026/analysis/partition_categories.json`
- Partition files: `/data/ecod/pdb_updates/batches/ecod_q4_2025_q1_2026/partitions/`
