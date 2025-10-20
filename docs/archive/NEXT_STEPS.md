# Next Steps - Implementation Guide

**Date**: 2025-10-19
**Current Status**: 📋 Planning Complete

## Important: Separation of Concerns

**This repository contains the domain partitioning ALGORITHM only.**

Production infrastructure (SLURM, batch orchestration, BLAST/HHsearch execution) lives in the separate **pyecod_prod** repository.

**See**: [PYECOD_MINI_API_SPEC.md](PYECOD_MINI_API_SPEC.md) for the formal API contract between pyecod_mini and pyecod_prod.

## What We've Accomplished

✅ **Complete Planning Documentation**:
1. [EXTRACTION_PLAN.md](EXTRACTION_PLAN.md) - Detailed roadmap for extracting mini from legacy repo
2. [ALGORITHM_VALIDATION.md](ALGORITHM_VALIDATION.md) - Comprehensive validation strategy
3. [PYECOD_MINI_API_SPEC.md](PYECOD_MINI_API_SPEC.md) - Formal API specification for integration with pyecod_prod
4. [README.md](README.md) - Project overview and quick start

## Immediate Next Steps

### Step 1: Create Repository Structure (Day 1)

**Goal**: Set up the modern Python project structure

**Tasks**:
```bash
# 1. Create directory structure (NO production/ - that's in pyecod_prod)
mkdir -p src/pyecod_mini/{core,cli}
mkdir -p tests data scripts docs

# 2. Create __init__.py files
touch src/pyecod_mini/__init__.py
touch src/pyecod_mini/core/__init__.py
touch src/pyecod_mini/cli/__init__.py
touch tests/__init__.py

# 3. Create pyproject.toml
# (Copy from EXTRACTION_PLAN.md section 1.2)

# 4. Initialize git repository (if not already done)
git init
git add .
git commit -m "Initial repository structure"
```

**Deliverable**: Clean repository structure with modern Python packaging

---

### Step 2: Extract Core Algorithm (Days 2-3)

**Goal**: Extract core modules from ../pyecod/mini with type hints

**Priority Order** (from EXTRACTION_PLAN.md Section 2.1):

#### 2A. Core Data Models (Critical Foundation)

```bash
# Start with sequence_range (already well-structured)
cp ../pyecod/mini/core/sequence_range.py src/pyecod_mini/core/
# Minor cleanup, add type hints if missing

# Then models.py (MOST CRITICAL - needs significant refactoring)
# This file needs careful conversion from dict-based to pure dataclasses
```

**For models.py refactoring**:
1. Review current dataclass definitions
2. Identify any remaining dict usage
3. Add comprehensive type hints to ALL attributes and methods
4. Add Google-style docstrings
5. Run `mypy src/pyecod_mini/core/models.py` and fix all errors

#### 2B. Evidence Parsing & Processing

```bash
# Extract in order (dependencies flow downward):
cp ../pyecod/mini/core/blast_parser.py src/pyecod_mini/core/
cp ../pyecod/mini/core/parser.py src/pyecod_mini/core/
cp ../pyecod/mini/core/evidence_utils.py src/pyecod_mini/core/
```

**Refactoring checklist for each**:
- [ ] Add type hints to all function signatures
- [ ] Add type hints to all variables
- [ ] Convert any dict returns to dataclasses
- [ ] Add/update docstrings
- [ ] Run mypy and fix errors
- [ ] Run black for formatting

#### 2C. Domain Partitioning Logic

```bash
cp ../pyecod/mini/core/decomposer.py src/pyecod_mini/core/
cp ../pyecod/mini/core/partitioner.py src/pyecod_mini/core/
cp ../pyecod/mini/core/boundary_optimizer.py src/pyecod_mini/core/
cp ../pyecod/mini/core/domain_utils.py src/pyecod_mini/core/
cp ../pyecod/mini/core/gap_analyzer.py src/pyecod_mini/core/
```

**Same refactoring checklist applies**

#### 2D. Output & Utilities

```bash
cp ../pyecod/mini/core/writer.py src/pyecod_mini/core/
cp ../pyecod/mini/core/ecod_domains_parser.py src/pyecod_mini/core/
cp ../pyecod/mini/core/range_cache_parser.py src/pyecod_mini/core/
cp ../pyecod/mini/core/visualization.py src/pyecod_mini/core/  # Optional
```

**Deliverable**: All core algorithm modules extracted with type hints and dataclasses

**Validation**: After each file extraction:
```bash
# Type check
mypy src/pyecod_mini/core/<filename>.py

# Format
black src/pyecod_mini/core/<filename>.py

# Verify imports work
python -c "from pyecod_mini.core.<modulename> import *"
```

---

### Step 3: Extract and Adapt Tests (Day 3)

**Goal**: Migrate test suite and ensure regression tests pass

#### 3A. Test Infrastructure

```bash
# Copy conftest.py with path updates
cp ../pyecod/mini/tests/conftest.py tests/
# Update paths to reference new src/pyecod_mini structure

# Install test dependencies
pip install pytest pytest-cov
```

#### 3B. Critical Regression Tests (PRIORITY 1)

```bash
# This is THE most important test file
cp ../pyecod/mini/tests/test_ecod_regression.py tests/

# Update imports in test file:
# OLD: from mini.core.parser import ...
# NEW: from pyecod_mini.core.parser import ...
```

**Run regression tests**:
```bash
pytest tests/test_ecod_regression.py -v

# Expected output:
# test_8ovp_A_regression PASSED
# (Plus 5 more regression tests once defined)
```

#### 3C. Unit Tests

```bash
# Copy all test files
cp ../pyecod/mini/tests/test_core.py tests/
cp ../pyecod/mini/tests/test_decomposer.py tests/
cp ../pyecod/mini/tests/test_discontinuous.py tests/
cp ../pyecod/mini/tests/test_parser.py tests/
cp ../pyecod/mini/tests/test_writer.py tests/
cp ../pyecod/mini/tests/test_models.py tests/
cp ../pyecod/mini/tests/test_batch_proteins.py tests/
cp ../pyecod/mini/tests/test_ecod_tgroup.py tests/

# Update imports in ALL test files
# Use find-and-replace: mini.core -> pyecod_mini.core
```

**Run full test suite**:
```bash
pytest tests/ -v

# Goal: All tests pass (may need some import/path fixes)
```

**Deliverable**: Complete test suite passing with new structure

---

### Step 4: Extract Reference Data (Day 3)

**Goal**: Copy reference data files needed for algorithm

```bash
# Copy reference data
cp ../pyecod/mini/test_data/domain_lengths.csv data/
cp ../pyecod/mini/test_data/protein_lengths.csv data/
cp ../pyecod/mini/test_data/domain_definitions.csv data/

# Copy blacklist if it exists
[ -f ../pyecod/mini/test_data/reference_blacklist.csv ] && \
  cp ../pyecod/mini/test_data/reference_blacklist.csv data/

# Verify data files
ls -lh data/
```

**Deliverable**: All reference data in place

---

### Step 5: Extract and Refactor CLI (Day 4)

**Goal**: Extract CLI functionality into modular structure

#### 5A. Split pyecod_mini.py into modules

The current `../pyecod/mini/pyecod_mini.py` is ~800 lines. We need to split it:

```python
# src/pyecod_mini/cli/batch_finder.py
# Extract: BatchFinder class (lines ~39-162)

# src/pyecod_mini/cli/config.py
# Extract: PyEcodMiniConfig class (lines ~164-300)

# src/pyecod_mini/cli/main.py
# Extract: Main CLI logic, partition_protein function, argparse setup
```

**Steps**:
1. Create `batch_finder.py` with BatchFinder class + type hints
2. Create `config.py` with PyEcodMiniConfig class + type hints
3. Create `main.py` with CLI entry point
4. Update `src/pyecod_mini/__main__.py` to call main.py
5. Test CLI works: `python -m pyecod_mini --help`

**Deliverable**: Modular CLI that works like old mini

**Validation**:
```bash
# Should work after refactoring
python -m pyecod_mini 8ovp_A --verbose

# Or via entry point (once installed)
pyecod-mini 8ovp_A --verbose
```

---

### Step 6: Algorithm Validation (Days 5-6)

**Goal**: Comprehensive validation per ALGORITHM_VALIDATION.md

#### 6A. Run Regression Tests

```bash
# Run the critical 6 regression tests
pytest tests/test_ecod_regression.py -v --tb=short

# All 6 must pass
```

**Expected results** (from ALGORITHM_VALIDATION.md):
- Test 1 (8ovp_A): 3 domains, 80%+ coverage, GFP+PBP families ✅
- Tests 2-6: TBD (need to select and validate)

#### 6B. Performance Benchmarks

```bash
# Run performance tests
pytest tests/ -m performance -v

# Check metrics:
# - Single protein: < 10 seconds
# - Memory: < 500 MB
```

#### 6C. Run on 100-Protein Test Set

```bash
# Create 100-protein test set (diverse cases)
python scripts/create_test_set.py --count 100 --output test_proteins_100.txt

# Process all 100
python scripts/batch_process.py --proteins test_proteins_100.txt

# Generate validation report
python scripts/validate_algorithm.py --proteins test_proteins_100.txt --report validation_report.md
```

**Success Criteria** (from ALGORITHM_VALIDATION.md):
- ✅ Success rate >= 95%
- ✅ Average coverage >= 75%
- ✅ ECOD classification accuracy >= 80%
- ✅ No crashes or errors

**Deliverable**: Validation report confirming algorithm quality

---

### Step 7: Implement Library API (Day 7)

**Goal**: Implement public library API per PYECOD_MINI_API_SPEC.md

#### 7A. Create Public API Module

```python
# src/pyecod_mini/__init__.py
"""
pyecod_mini - Domain partitioning algorithm for ECOD.

Public API for integration with pyecod_prod and other tools.
"""

from .api import partition_protein, PartitionResult, Domain, PartitionError

__version__ = "1.0.0"
__all__ = ["partition_protein", "PartitionResult", "Domain", "PartitionError"]
```

#### 7B. Implement partition_protein Function

```python
# src/pyecod_mini/api.py
"""
Public API for pyecod_mini domain partitioning.

See PYECOD_MINI_API_SPEC.md for complete specification.
"""

from pathlib import Path
from typing import Optional
from dataclasses import dataclass
import datetime

def partition_protein(
    summary_xml: str | Path,
    output_xml: str | Path | None = None,
    *,
    pdb_id: str | None = None,
    chain_id: str | None = None,
    batch_id: str | None = None,
    verbose: bool = False,
) -> PartitionResult:
    """
    Partition a protein chain into domains based on evidence.

    See PYECOD_MINI_API_SPEC.md for complete documentation.
    """
    # Implementation delegates to existing core modules
    pass

@dataclass
class PartitionResult:
    """Result from domain partitioning (see API spec)."""
    pdb_id: str
    chain_id: str
    sequence_length: int
    sequence: str
    domains: list[Domain]
    coverage: float
    algorithm_version: str
    success: bool
    error_message: Optional[str] = None
    batch_id: Optional[str] = None
    timestamp: Optional[str] = None

@dataclass
class Domain:
    """A partitioned domain (see API spec)."""
    domain_id: str
    range_string: str
    residue_count: int
    ecod_domain_id: str
    family_name: str
    source: str
    confidence: Optional[float] = None

class PartitionError(Exception):
    """Raised when domain partitioning fails."""
    pass
```

#### 7C. Update XML Writer

Modify `core/writer.py` to output partition.xml format per spec:
- Include `<coverage>` element
- Include `algorithm_version` attribute
- Follow exact schema from PYECOD_MINI_API_SPEC.md

#### 7D. Test Library API

```python
# tests/test_api.py
"""Test public library API."""

from pyecod_mini import partition_protein, PartitionError

def test_partition_protein_success():
    """Test successful partitioning via library API."""
    result = partition_protein(
        summary_xml="tests/data/8ovp_A.summary.xml",
        output_xml="tests/output/8ovp_A.partition.xml"
    )

    assert result.success
    assert result.coverage > 0.0
    assert len(result.domains) > 0
    assert result.algorithm_version == "1.0.0"

def test_partition_protein_file_not_found():
    """Test error handling for missing file."""
    with pytest.raises(FileNotFoundError):
        partition_protein(summary_xml="nonexistent.xml")

def test_partition_protein_invalid_xml():
    """Test error handling for corrupt XML."""
    with pytest.raises(PartitionError):
        partition_protein(summary_xml="tests/data/corrupt.xml")
```

**Deliverable**: Clean library API that matches PYECOD_MINI_API_SPEC.md

---

### Step 8: Integration with pyecod_prod (Day 8)

**Goal**: Ensure pyecod_mini works correctly when called by pyecod_prod

#### 8A. Create Integration Test

```python
# tests/test_pyecod_prod_integration.py
"""
Test integration with pyecod_prod.

Simulates how pyecod_prod will call pyecod_mini.
"""

def test_prod_workflow_simulation():
    """Simulate pyecod_prod calling pyecod_mini."""

    # 1. pyecod_prod generates summary XML (simulate this)
    summary_xml = create_test_summary_xml("tests/output/test_summary.xml")

    # 2. pyecod_prod calls pyecod_mini library
    from pyecod_mini import partition_protein

    result = partition_protein(
        summary_xml=summary_xml,
        output_xml="tests/output/test_partition.xml",
        batch_id="test_batch"
    )

    # 3. pyecod_prod uses results (coverage, domains)
    assert result.coverage >= 0.0
    assert result.coverage <= 1.0

    # 4. pyecod_prod applies quality assessment (its policy, not ours)
    if result.coverage >= 0.80:
        quality = "good"  # This is pyecod_prod's decision
    else:
        quality = "low_coverage"

    print(f"Quality (pyecod_prod policy): {quality}")
```

#### 8B. Coordinate with pyecod_prod

1. Ensure pyecod_prod has updated `partition_runner.py` to use library API
2. Test on small batch through full pyecod_prod pipeline
3. Verify manifest tracking works correctly
4. Validate XML format compatibility

#### 8C. Version Pinning Test

```bash
# In pyecod_prod, test version pinning
pip install pyecod-mini==1.0.0

# Run pyecod_prod workflow
cd ../pyecod_prod
python scripts/run_small_test.py

# Verify partition step uses correct version
grep "algorithm_version" /data/ecod/test_batches/*/partitions/*.xml
```

**Deliverable**: Validated integration with pyecod_prod

---

### Step 9: Documentation & Release (Day 9)

**Goal**: Complete documentation and prepare for release

#### 9A. Generate API Documentation

```bash
# Use pdoc or sphinx
pdoc --html src/pyecod_mini -o docs/api
```

#### 9B. Update Documentation

- [ ] Algorithm description (docs/algorithm.md)
- [ ] Library API usage guide (docs/api_usage.md)
- [ ] Development guide (docs/development.md)
- [ ] Troubleshooting guide (docs/troubleshooting.md)

#### 9C. Create Release Checklist

- [ ] All tests pass (`pytest tests/ -v`)
- [ ] Type checking passes (`mypy src/`)
- [ ] Code formatted (`black src/ tests/`)
- [ ] Linting clean (`ruff check src/`)
- [ ] API spec compliance verified
- [ ] Integration with pyecod_prod tested
- [ ] Documentation complete
- [ ] Version tagged (v1.0.0)

#### 9D. Package and Publish

```bash
# Build package
python -m build

# Test installation
pip install dist/pyecod_mini-1.0.0-py3-none-any.whl

# Verify
python -c "from pyecod_mini import partition_protein; print('OK')"
```

**Deliverable**: Released pyecod_mini v1.0.0 ready for production use

---

## Success Milestones

### Milestone 1: Core Algorithm Extraction ✨
- ✅ All core modules extracted
- ✅ Type hints on all functions
- ✅ mypy passes with no errors
- ✅ All unit tests pass

### Milestone 2: Validation Complete ✨
- ✅ 6 regression tests pass
- ✅ 100-protein test batch successful
- ✅ Performance benchmarks met
- ✅ Validation report generated

### Milestone 3: Library API Complete ✨
- ✅ Public API implemented per PYECOD_MINI_API_SPEC.md
- ✅ Library API tests pass
- ✅ XML output format matches spec
- ✅ Coverage calculation included in output

### Milestone 4: Integration Complete ✨
- ✅ Integration with pyecod_prod validated
- ✅ Small batch test through full pipeline
- ✅ Version pinning tested
- ✅ Documentation complete
- ✅ v1.0.0 released

## Quick Reference Commands

```bash
# Development
pip install -e ".[dev]"               # Install in dev mode
pytest tests/ -v                      # Run all tests
pytest tests/ -m regression           # Run regression tests only
mypy src/                             # Type checking
black src/ tests/                     # Format code
ruff check src/                       # Lint code

# Algorithm Testing (CLI)
pyecod-mini 8ovp_A --verbose          # Test single protein
pyecod-mini --test-suite              # Run formal test suite
python scripts/validate_algorithm.py  # Comprehensive validation

# Library API Testing (Python)
python -c "from pyecod_mini import partition_protein; print('OK')"
python scripts/test_library_api.py    # Test library integration

# Package Build
python -m build                       # Build distribution
pip install dist/pyecod_mini-*.whl    # Install built package
```

## Troubleshooting

### Common Issues

**Issue**: Tests fail with import errors
**Solution**: Ensure you're in project root and package is installed: `pip install -e .`

**Issue**: mypy reports type errors
**Solution**: Add type hints incrementally, use `# type: ignore` for complex cases temporarily

**Issue**: Regression tests fail
**Solution**: Compare with old mini output, may need to adjust paths or reference data

**Issue**: Integration with pyecod_prod fails
**Solution**: Verify API spec compliance, check XML format, ensure version compatibility

## Getting Help

- Check planning documents (EXTRACTION_PLAN.md, ALGORITHM_VALIDATION.md)
- Review PYECOD_MINI_API_SPEC.md for integration questions
- Run with `--verbose` for detailed output
- Check test output for validation errors

---

## Separation of Concerns Reminder

**pyecod_mini** (this repo):
- ✅ Domain partitioning algorithm
- ✅ Evidence parsing (BLAST XML, HHR files)
- ✅ Coverage calculation
- ✅ Library API + simple CLI
- ❌ NO SLURM, NO batch orchestration, NO production infrastructure

**pyecod_prod** (separate repo):
- ✅ PDB data acquisition
- ✅ BLAST/HHsearch execution via SLURM
- ✅ Batch workflow orchestration
- ✅ Quality policy decisions
- ✅ Database integration
- ✅ Calls pyecod_mini as library/CLI

See [PYECOD_MINI_API_SPEC.md](PYECOD_MINI_API_SPEC.md) for the formal contract.

---

**Ready to Begin!** Start with Step 1 and work through systematically.

**Estimated Timeline**: 9 days from start to release
