# Performance Issue: Domain Definitions Reloaded on Every partition_protein() Call

**Date**: 2025-10-25
**Version**: pyecod_mini 2.0.0
**Reporter**: Production batch processing (pyecod_prod)

## Issue Summary

When processing large batches of chains (~4,000+), `partition_protein()` reloads 718,436 ECOD domain definitions on **every call**, adding ~13-14 seconds overhead per chain.

## Expected Behavior

Domain definitions should be loaded **once** (either at module import or on first call) and cached for reuse across multiple `partition_protein()` calls.

## Actual Behavior

Each call to `partition_protein()` prints:
```
Loaded domain definitions for 718436 chain entries
Excluded 2 domains from 1 blacklisted chains
```

This happens even when:
- Using the library API (not CLI)
- Reusing the same Python process
- Calling `partition_protein()` multiple times sequentially

## Performance Impact

**Test case**: Process 4,038 chains with HHsearch results

| Approach | Time per chain | Total time | Notes |
|----------|---------------|------------|-------|
| Sequential (no caching) | 14 sec | 15.7 hours | Current behavior |
| Sequential (with caching) | <1 sec | ~1 hour | Expected if cached |
| Parallel 32 workers (no caching) | 14 sec | 29 minutes | Workaround |
| Parallel 32 workers (with caching) | <1 sec | <2 minutes | Ideal |

## Reproduction

```python
from pyecod_mini import partition_protein
import time

# Test with reused runner
for i in range(3):
    start = time.time()
    result = partition_protein(
        summary_xml=f"test_chain_{i}.summary.xml",
        output_xml=f"test_chain_{i}.partition.xml",
        pdb_id="test",
        chain_id=str(i),
        batch_id="test"
    )
    print(f"Call {i+1}: {time.time() - start:.2f}s")
    # Expected: Call 1: 14s, Call 2-3: <1s
    # Actual: Call 1-3: 14s each
```

## Suggested Fixes

### Option 1: Module-level singleton cache (Recommended)
```python
# In core/domain_lookup.py or similar
_DOMAIN_CACHE = None

def get_domain_definitions():
    global _DOMAIN_CACHE
    if _DOMAIN_CACHE is None:
        print("Loading domain definitions (once)...")
        _DOMAIN_CACHE = _load_definitions()
    return _DOMAIN_CACHE
```

### Option 2: Explicit DomainLoader class
```python
from pyecod_mini import DomainLoader, partition_protein

# User creates loader once
loader = DomainLoader()  # Loads 718K definitions once

# Reuse loader across calls
for chain in chains:
    result = partition_protein(
        summary_xml=chain.summary,
        output_xml=chain.output,
        domain_loader=loader  # Pass cached loader
    )
```

### Option 3: Lazy loading with TTL cache
```python
from functools import lru_cache

@lru_cache(maxsize=1)
def _load_domain_definitions():
    # Load once, cache forever
    return load_ecod_domains()
```

## Workaround (Current)

Use multiprocessing to parallelize across chains:
```bash
# 32 workers = ~29 minutes for 4,038 chains
python3 partition_parallel.py --workers 32
```

Each worker still pays the 14-second penalty, but parallelism compensates.

## Impact on Production

- **Batch processing**: 15.7 hours → <1 hour with caching (15x speedup)
- **Interactive use**: Less critical (single chains are fine)
- **Weekly releases**: ~1,700 chains/week → 6.6 hours → <30 minutes

## Priority

**High** - Affects all batch processing workflows in production environment.

## Related Code

- `src/pyecod_mini/core/partition.py` - Main partition logic
- `src/pyecod_mini/production/ecod_loader.py` - Domain loading (probable location)
- `src/pyecod_mini/api.py` - Public API entry point

---

**Next Steps**:
1. Locate domain loading code in pyecod_mini
2. Implement module-level caching (Option 1)
3. Add tests to verify caching works across calls
4. Update API documentation if caching behavior changes
