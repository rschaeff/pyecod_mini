# Deficiency Report: Inefficient Reference Data Loading in Batch Processing

**Date**: 2026-01-15
**Reporter**: Production Pipeline (pyecod_prod)
**Severity**: Performance - Critical for batch processing
**Affects**: pyecod_mini v2.x

## Summary

When processing multiple chains in batch mode, pyecod_mini loads the full domain definitions database (~718,436 entries) for every single chain. This creates a significant performance bottleneck, adding approximately 10-15 seconds of overhead per chain just for data loading.

## Observed Behavior

From SLURM job logs during Q4 2025 batch processing (3,677 chains):

```
Processing 8zmz_A (task 100, line 100)
Loaded domain definitions for 718436 chain entries    <-- This line repeats for EVERY chain
Excluded 2 domains from 1 blacklisted chains
...
```

### Impact

| Metric | Single Chain | 3,677 Chain Batch |
|--------|-------------|-------------------|
| Domain loading time | ~10-15 sec | ~10-15 hours overhead |
| Evidence processing | 1-2 min | Variable |
| Total per chain | ~2 min | ~2 min each |
| Batch total | N/A | ~12+ hours |

With the domain definitions loaded once and reused:
- Estimated batch time: ~4-6 hours (3x faster)

## Root Cause

The `partition_protein()` function and `PartitionRunner` class load reference data on every invocation:

```python
# In pyecod_mini partitioner - pseudocode
def partition_protein(summary_xml, output_xml, ...):
    domain_defs = load_domain_definitions()  # 718K entries loaded every time
    evidence = parse_evidence(summary_xml)
    result = run_partitioning(evidence, domain_defs)
    return result
```

## Expected Behavior

For batch processing, reference data should be loaded once and reused:

```python
# Preferred API for batch processing
class BatchPartitioner:
    def __init__(self, reference_version="develop291"):
        self.domain_defs = load_domain_definitions(reference_version)  # Load once

    def partition(self, summary_xml, output_xml, ...):
        evidence = parse_evidence(summary_xml)
        return run_partitioning(evidence, self.domain_defs)  # Reuse
```

## Recommended Fix

### Option 1: Add Batch Mode API (Recommended)

Add a class-based API that maintains state across multiple partitioning calls:

```python
from pyecod_mini import BatchPartitioner

# Initialize once - loads domain definitions
partitioner = BatchPartitioner(reference_version="develop291")

# Process multiple chains efficiently
for summary_file in summary_files:
    result = partitioner.partition(summary_file, output_file)
```

### Option 2: Lazy Loading with Cache

Implement module-level caching of domain definitions:

```python
_DOMAIN_CACHE = {}

def get_domain_definitions(version="develop291"):
    if version not in _DOMAIN_CACHE:
        _DOMAIN_CACHE[version] = load_domain_definitions(version)
    return _DOMAIN_CACHE[version]
```

### Option 3: Memory-Mapped Reference Data

For very large batches, consider memory-mapping the domain definitions file to reduce per-process memory overhead when running parallel SLURM jobs.

## Workaround (Current)

The current workaround is to accept the overhead and run parallel SLURM jobs. Each job pays the loading cost but runs independently. This is suboptimal but functional.

## Affected Code Paths

1. `pyecod_mini.partition_protein()` - Main entry point
2. `pyecod_mini.partitioner.DomainPartitioner.__init__()` - Loads definitions
3. SLURM batch scripts calling pyecod_mini per-chain

## Test Case

```python
import time
from pyecod_mini import partition_protein

summary_files = [f"chain_{i}.summary.xml" for i in range(100)]

# Current (slow) approach
start = time.time()
for f in summary_files:
    partition_protein(f, f.replace('.summary.', '.partition.'))
print(f"Current: {time.time() - start:.1f}s")  # ~200s (100 * 2s)

# With batch API (expected improvement)
# start = time.time()
# partitioner = BatchPartitioner()
# for f in summary_files:
#     partitioner.partition(f, f.replace('.summary.', '.partition.'))
# print(f"Batch: {time.time() - start:.1f}s")  # ~100s (10s load + 100 * 0.9s)
```

## Priority

**High** - This affects all batch processing workflows and is the primary bottleneck for large-scale ECOD updates.

## Related Issues

- Q4 2025 + Q1 2026 batch: 3,677 chains taking 12+ hours instead of ~4 hours
- 2-year backfill: 9,656 representatives affected
