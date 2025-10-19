# Step 3 Complete: Test Suite Extracted

**Date**: 2025-10-19
**Status**: ✅ Complete
**Commit**: 84a1f07

## What We Accomplished

### Extracted Complete Test Suite

All test files successfully copied from `../pyecod/mini/tests/` with updated imports:

#### Test Configuration
- **conftest.py** (14 KB) - pytest fixtures and test environment setup
  - Updated imports: `mini.core` → `pyecod_mini.core`
  - TestEnvironment class for batch detection
  - ReferenceDataLoader for lazy-loading test data
  - Mock and real data fixtures

#### Critical Regression Tests ⭐
- **test_ecod_regression.py** (28 KB) - **The most important file!**
  - Contains the 6 critical regression tests
  - Validates algorithm correctness
  - Tests 8ovp_A (GFP-PBP fusion) as primary case
  - Compares mini vs current engine
  - Expected: 3 domains, ~80% coverage

#### Core Algorithm Tests (Unit Tests)
- **test_core.py** (21 KB) - Core partitioning logic
- **test_parser.py** (11 KB) - Evidence parsing
- **test_models.py** (9.6 KB) - Data models
- **test_decomposer.py** (11 KB) - Chain decomposition
- **test_discontinuous.py** (13 KB) - Multi-segment domains
- **test_blast_parser.py** (9.4 KB) - BLAST XML parsing
- **test_writer.py** (7.8 KB) - XML output generation

#### Classification & Integration Tests
- **test_ecod_tgroup.py** (8.4 KB) - ECOD T-group classification
- **test_range_cache_parser.py** (9.7 KB) - Reference data parsing
- **test_standalone_regression.py** (18 KB) - Standalone regression tests
- **test_batch_proteins.py** (5.9 KB) - Batch processing
- **test_cli.py** (4.3 KB) - CLI functionality
- **test_cases.py** (18 KB) - Test case definitions

**Total**: 16 test files (15 tests + conftest.py)

### Import Updates

All imports successfully updated using sed:

**Before:**
```python
from mini.core.models import Evidence, Domain
from mini.core.parser import parse_domain_summary
from ecod.core.sequence_range import SequenceRange
```

**After:**
```python
from pyecod_mini.core.models import Evidence, Domain
from pyecod_mini.core.parser import parse_domain_summary
from pyecod_mini.core.sequence_range import SequenceRange
```

### Reference Data Copied

All test reference data copied to `data/` directory:

| File | Size | Purpose |
|------|------|---------|
| domain_lengths.csv | 14 MB | ECOD domain reference lengths |
| protein_lengths.csv | 8.3 MB | Protein chain lengths |
| domain_definitions.csv | 76 MB | Multi-domain decomposition rules |
| ecod_classifications.csv | 96 B | ECOD T-group/H-group classifications |
| reference_blacklist.csv | 117 B | Reference quality blacklist |

Plus additional metadata files.

**Total reference data**: ~99 MB

## Current Status

### ✅ Completed
- [x] All 15 test files extracted
- [x] conftest.py with fixtures
- [x] All imports updated to pyecod_mini.core
- [x] Reference data copied to data/
- [x] Git committed

### ⏳ Pending (Next Steps)
- [ ] Install development dependencies: `pip install -e '.[dev]'`
- [ ] Run full test suite: `pytest tests/ -v`
- [ ] Run regression tests: `pytest tests/test_ecod_regression.py -v`
- [ ] Verify 6 regression tests pass
- [ ] Debug any failing tests
- [ ] Generate test coverage report

## How to Run Tests

Once pytest is installed in a proper Python environment:

### 1. Install Development Dependencies
```bash
pip install -e '.[dev]'
```

This installs:
- pytest (testing framework)
- pytest-cov (coverage reporting)
- mypy (type checking)
- black (code formatting)
- ruff (linting)

### 2. Run All Tests
```bash
pytest tests/ -v
```

### 3. Run Regression Tests Only
```bash
pytest tests/test_ecod_regression.py -v -m regression
```

### 4. Run Unit Tests Only
```bash
pytest tests/ -v -m unit
```

### 5. Run with Coverage
```bash
pytest tests/ --cov=src/pyecod_mini --cov-report=html
```

## Test Markers

From conftest.py, tests can be marked as:
- `unit` - Fast unit tests with mock data
- `integration` - Integration tests requiring real data
- `regression` - Critical regression tests
- `slow` - Slow-running tests (skip with `-m "not slow"`)
- `performance` - Performance benchmarks
- `visualization` - Requires PyMOL

## Expected Regression Test Results

### Primary Test: 8ovp_A

**Expected Output:**
- **Domains found**: 3 (or 2 if discontinuous merged)
- **Coverage**: >= 80% of sequence
- **Families identified**:
  - GFP family (6dgv)
  - PBP family (2vha/2ia4)
- **Discontinuous handling**: Correct (PBP domain split by GFP)

**Domain ranges (approximate):**
1. PBP N-terminal: ~2-248
2. GFP insertion: ~253-499
3. PBP C-terminal: ~500-517

### Full Regression Suite

Once complete, we expect:
- ✅ All 6 regression tests pass
- ✅ Domain count within ±1 of expected
- ✅ Coverage >= 75% average
- ✅ ECOD classification accuracy >= 80%
- ✅ No crashes or errors

## File Structure

```
pyecod_mini/
├── tests/
│   ├── __init__.py
│   ├── conftest.py                    # Test fixtures
│   ├── test_ecod_regression.py        # ⭐ CRITICAL
│   ├── test_core.py                   # Core algorithm
│   ├── test_parser.py                 # Evidence parsing
│   ├── test_models.py                 # Data models
│   ├── test_decomposer.py             # Decomposition
│   ├── test_discontinuous.py          # Multi-segment
│   ├── test_blast_parser.py           # BLAST parsing
│   ├── test_writer.py                 # XML output
│   ├── test_ecod_tgroup.py            # Classification
│   ├── test_range_cache_parser.py     # Reference data
│   ├── test_standalone_regression.py  # Standalone tests
│   ├── test_batch_proteins.py         # Batch processing
│   ├── test_cli.py                    # CLI
│   └── test_cases.py                  # Test cases
│
├── data/                              # Reference data
│   ├── domain_lengths.csv             # 14 MB
│   ├── protein_lengths.csv            # 8.3 MB
│   ├── domain_definitions.csv         # 76 MB
│   ├── ecod_classifications.csv
│   ├── reference_blacklist.csv
│   └── ... (metadata files)
│
└── src/pyecod_mini/core/              # Algorithm (from Step 2)
    └── ... (15 core modules)
```

## Statistics

- **Test files**: 16 (15 tests + conftest)
- **Lines of test code**: ~5,200 lines
- **Core modules**: 15 files
- **Reference data**: ~99 MB (17 CSV files)
- **Total commits**: 4

## Next Steps

According to NEXT_STEPS.md, we have two paths:

### Option A: Validate Now (Recommended)
1. Set up dev environment: `pip install -e '.[dev]'`
2. Run regression tests: `pytest tests/test_ecod_regression.py -v`
3. Debug and fix any failures
4. **Goal**: Prove algorithm works in new structure

### Option B: Polish First
1. Add complete type hints to core modules
2. Run mypy type checking
3. Format with black
4. **Then** validate with tests

**Recommended**: **Option A** - Run the regression tests now to validate the extracted algorithm works!

## Success Criteria

For the regression tests to pass, we need:
- ✅ pytest installed and working
- ✅ Reference data files accessible
- ✅ ECOD batch directories available (or tests skip gracefully)
- ✅ All imports resolve correctly
- ✅ Core algorithm functions correctly

## Troubleshooting

### "No module named 'pytest'"
**Solution**: Install dev dependencies:
```bash
pip install -e '.[dev]'
```

### "No ECOD batch directories available"
**Solution**: Set environment variable or tests will skip:
```bash
export ECOD_BATCH_DIR=/data/ecod/pdb_updates/batches/ecod_batch_036_20250406_1424
pytest tests/ -v
```

### Import errors
**Solution**: Ensure you're running from project root and package is installed:
```bash
cd /home/rschaeff/dev/pyecod_mini
pip install -e .
pytest tests/ -v
```

## Summary

Step 3 is **successfully complete**! We have:
- ✅ All 15 test files extracted with updated imports
- ✅ conftest.py with comprehensive fixtures
- ✅ ~99 MB of reference data
- ✅ Ready to validate algorithm correctness
- ✅ Committed to git (commit 84a1f07)

The proven mini algorithm and its complete test suite are now in the new repository structure!

**Current status**: 🎯 **Tests Extracted - Ready to Validate**

**Next**: Install dev dependencies and run regression tests to prove the algorithm works!

```bash
pip install -e '.[dev]'
pytest tests/test_ecod_regression.py -v
```
