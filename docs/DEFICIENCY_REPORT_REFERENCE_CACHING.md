# Deficiency Report: Reference Data Caching for Batch Processing

**Date**: 2026-01-18
**Reporter**: Production Pipeline (pyecod_prod)
**Severity**: High (Performance)
**Component**: pyecod_mini core partitioning

## Summary

The `partition_protein()` API function loads reference data (domain definitions, reference lengths) fresh for every call, causing significant performance degradation in batch processing scenarios. Each partition operation takes ~60-120 seconds, with ~90% of that time spent loading reference data that could be cached and reused.

## Impact

### Current Behavior

When processing a batch of N chains:
- Each call to `partition_protein()` loads:
  - 718,436 domain definitions (~30-60 seconds)
  - 1,083,021 reference lengths (~10-20 seconds)
  - ECOD hierarchy data
- Total time: N × (60-120 seconds) = **~2 minutes per chain**

### Production Impact

For the Q4 2025/Q1 2026 batch (3,677 chains):
- Current estimate: 3,677 × 90 sec = **92 hours** (single-threaded)
- With 50 concurrent SLURM jobs: **~2 hours** but wastes cluster resources
- Reference data loaded 3,677 times instead of once

### Desired Behavior

- Load reference data once at batch start
- Reuse for all chains in batch
- Per-chain processing: ~2-5 seconds (actual algorithm time)
- 3,677 chains: **~20-30 minutes** (50× improvement)

## Technical Details

### Current Code Path

```python
# api.py - partition_protein()
def partition_protein(summary_xml, output_xml, pdb_id, chain_id, ...):
    # Every call does this:
    parser = DomainSummaryParser()  # Loads reference data
    result = parser.parse(summary_xml)  # Uses reference data
    domains = partition_domains(...)  # Algorithm
    writer.write(output_xml)  # Output
```

### Reference Data Sources

1. **Domain definitions** (`ecod_domain_definitions.tsv`)
   - 718,436 entries mapping (pdb_id, chain_id) → domain list
   - ~50MB file, ~30-60 sec load time

2. **Reference lengths** (`domain_family_lookup.tsv` or similar)
   - 1,083,021 domain → length mappings
   - Used for coverage calculations

3. **ECOD hierarchy** (optional, for family names)
   - X-group, H-group, T-group, F-group mappings

## Proposed Solutions

### Option 1: Singleton Cache (Minimal Change)

```python
# core/reference_cache.py
class ReferenceCache:
    _instance = None
    _domain_defs = None
    _reference_lengths = None

    @classmethod
    def get_instance(cls):
        if cls._instance is None:
            cls._instance = cls()
            cls._instance._load_data()
        return cls._instance

    def _load_data(self):
        # Load once, cache forever
        self._domain_defs = load_domain_definitions()
        self._reference_lengths = load_reference_lengths()
```

**Pros**: Minimal API change, automatic caching
**Cons**: Memory persists for process lifetime, no explicit control

### Option 2: Context Manager for Batch Processing

```python
# New batch API
from pyecod_mini import PartitionBatch

with PartitionBatch() as batch:
    for chain in chains:
        result = batch.partition(summary_xml, output_xml, ...)
```

**Pros**: Explicit control, clean memory release
**Cons**: New API pattern, requires code changes

### Option 3: Pre-loaded Partitioner Class

```python
# New class-based API
from pyecod_mini import Partitioner

partitioner = Partitioner()  # Loads reference data once
partitioner.load_references()  # Explicit loading

for chain in chains:
    result = partitioner.partition(summary_xml, output_xml, ...)

partitioner.close()  # Release memory
```

**Pros**: OOP pattern familiar to users, explicit lifecycle
**Cons**: Different from current functional API

### Option 4: Lazy Loading with LRU Cache

```python
from functools import lru_cache

@lru_cache(maxsize=1)
def get_domain_definitions():
    return load_domain_definitions()

@lru_cache(maxsize=1)
def get_reference_lengths():
    return load_reference_lengths()
```

**Pros**: Simple implementation, automatic
**Cons**: Less control, memory not explicitly releasable

## Recommendation

**Option 3 (Pre-loaded Partitioner Class)** is recommended because:

1. Matches pyecod_prod's existing `PartitionRunner` pattern
2. Explicit lifecycle control important for long-running batch jobs
3. Can preserve backward compatibility by keeping `partition_protein()` as convenience wrapper
4. Memory management is explicit and predictable

### Proposed API

```python
# Backward compatible - works as before (loads each time)
from pyecod_mini import partition_protein
result = partition_protein(summary_xml, output_xml, ...)

# New batch-optimized API
from pyecod_mini import Partitioner

# Initialize once, process many
partitioner = Partitioner()

for summary_xml, output_xml in batch:
    result = partitioner.partition(
        summary_xml=summary_xml,
        output_xml=output_xml,
        pdb_id=pdb_id,
        chain_id=chain_id,
    )

# Explicit cleanup (optional - also works with context manager)
partitioner.close()

# Or with context manager
with Partitioner() as p:
    for chain in chains:
        result = p.partition(...)
```

## Files to Modify

1. `src/pyecod_mini/core/partitioner.py` - Add `Partitioner` class
2. `src/pyecod_mini/core/parser.py` - Extract reference loading to reusable functions
3. `src/pyecod_mini/api.py` - Add new exports, update `partition_protein()` to optionally use cache
4. `src/pyecod_mini/__init__.py` - Export new `Partitioner` class

## Testing Requirements

1. Verify batch processing performance improvement
2. Ensure memory is properly released after batch
3. Backward compatibility - existing `partition_protein()` calls unchanged
4. Thread safety if applicable

## Timeline

This is blocking production processing of the Q4 2025/Q1 2026 batch (3,677 chains). Current workaround is running many parallel single-chain SLURM jobs, which is inefficient but functional.

## Related Issues

- pyecod_prod batch processing: `/data/ecod/pdb_updates/batches/ecod_q4_2025_q1_2026/`
- Current workaround: `scripts/submit_partition_v293.sh` (single-chain jobs)
