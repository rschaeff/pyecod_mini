# Next Steps - Implementation Guide

**Date**: 2025-10-19
**Current Status**: 📋 Planning Complete

## What We've Accomplished

✅ **Complete Planning Documentation**:
1. [EXTRACTION_PLAN.md](EXTRACTION_PLAN.md) - Detailed roadmap for extracting mini from legacy repo
2. [ALGORITHM_VALIDATION.md](ALGORITHM_VALIDATION.md) - Comprehensive validation strategy
3. [PRODUCTION_DESIGN.md](PRODUCTION_DESIGN.md) - Production framework architecture
4. [README.md](README.md) - Project overview and quick start

## Immediate Next Steps

### Step 1: Create Repository Structure (Day 1)

**Goal**: Set up the modern Python project structure

**Tasks**:
```bash
# 1. Create directory structure
mkdir -p src/pyecod_mini/{core,cli,production}
mkdir -p tests data scripts docs config

# 2. Create __init__.py files
touch src/pyecod_mini/__init__.py
touch src/pyecod_mini/core/__init__.py
touch src/pyecod_mini/cli/__init__.py
touch src/pyecod_mini/production/__init__.py
touch tests/__init__.py

# 3. Create pyproject.toml
# (Copy from EXTRACTION_PLAN.md section 1.2)

# 4. Create configuration templates
cp config.template.yml config/config.template.yml
# Edit with production settings

# 5. Initialize git repository
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

### Step 7: Production Framework Implementation (Days 7-10)

**Goal**: Build production components per PRODUCTION_DESIGN.md

#### 7A. Scanner Component (Day 7)

```python
# src/pyecod_mini/production/scanner.py
# Implement:
# - ProteinTask dataclass
# - ScanResult dataclass
# - BatchScanner class

# Test with:
pytest tests/test_production.py::test_scanner -v
```

#### 7B. SLURM Manager (Day 8)

```python
# src/pyecod_mini/production/slurm_manager.py
# Implement:
# - SlurmJob dataclass
# - JobBatch dataclass
# - SlurmManager class

# Test with small batch:
pytest tests/test_production.py::test_slurm_manager -v
```

#### 7C. Tracking Database (Day 8)

```python
# src/pyecod_mini/production/tracking_db.py
# Implement:
# - SQLite schema
# - TrackingDatabase class

# Test:
pytest tests/test_production.py::test_tracking_db -v
```

#### 7D. Monitor (Day 9)

```python
# src/pyecod_mini/production/monitor.py
# Implement:
# - ProcessingStats dataclass
# - ProductionMonitor class (with Rich dashboard)

# Manual test:
python -m pyecod_mini.production.monitor --test-mode
```

#### 7E. Database Importer (Day 9)

```python
# src/pyecod_mini/production/database.py
# Implement:
# - ImportResult dataclass
# - DatabaseImporter class

# Test with test database:
pytest tests/test_production.py::test_database_importer -v
```

#### 7F. Production CLI (Day 10)

```bash
# Add production subcommands to CLI
pyecod-mini production scan --help
pyecod-mini production process --help
pyecod-mini production monitor --help
pyecod-mini production import --help
```

**Deliverable**: Complete production framework

---

### Step 8: Integration Testing (Day 11)

**Goal**: End-to-end testing of production workflow

#### 8A. Small Test Batch (10 proteins)

```bash
# Full workflow test
pyecod-mini production scan --max-proteins 10 --output test_tasks.json
pyecod-mini production process --tasks test_tasks.json --max-concurrent 2
pyecod-mini production monitor --stats
pyecod-mini production import --dry-run
```

**Verify**:
- ✅ All 10 proteins processed
- ✅ Output XML files created
- ✅ Tracking database updated
- ✅ Quality metrics reasonable

#### 8B. Medium Test Batch (100 proteins)

```bash
# Scaling test
pyecod-mini production process --max-proteins 100 --max-concurrent 10

# Monitor performance:
# - Throughput: >= 400 proteins/hour
# - Success rate: >= 95%
# - Quality: >= 60% good
```

**Deliverable**: Proven production workflow on test data

---

### Step 9: Production Deployment (Days 12-14)

**Goal**: Process all ~40k representative proteins

#### 9A. Production Configuration

```bash
# Set up production config
cp config/production.template.yml config/production.yml
# Edit with real database credentials, paths

# Validate configuration
pyecod-mini production validate-config
```

#### 9B. Staged Rollout

```bash
# Stage 1: 1000 proteins
pyecod-mini production process --max-proteins 1000 --max-concurrent 50
# Monitor, assess quality, fix any issues

# Stage 2: 5000 proteins
pyecod-mini production process --max-proteins 5000 --max-concurrent 50
# Monitor, assess quality

# Stage 3: All remaining (~34k proteins)
pyecod-mini production process --reps-only --skip-existing --max-concurrent 50
# This will take ~24-48 hours
```

#### 9C. Quality Assessment

```bash
# Generate quality report
pyecod-mini production quality --output production_quality_report.json

# Verify metrics meet targets:
# - Success rate >= 95%
# - Good quality >= 60%
# - Poor quality <= 10%
```

#### 9D. Database Import

```bash
# Check for collisions
pyecod-mini production import --check-collisions --report collisions.txt

# Import with quality filter
pyecod-mini production import --quality-filter --collision-strategy skip

# Verify import
pyecod-mini production verify --count 1000
```

**Deliverable**: All representative proteins processed and imported

---

### Step 10: Documentation & Cleanup (Day 14)

**Goal**: Finalize documentation and prepare for production use

#### 10A. Generate API Documentation

```bash
# Use pdoc or sphinx
pdoc --html src/pyecod_mini -o docs/api
```

#### 10B. Update Documentation

- [ ] Algorithm description (docs/algorithm.md)
- [ ] Production usage guide (docs/production.md)
- [ ] Development guide (docs/development.md)
- [ ] Troubleshooting guide (docs/troubleshooting.md)

#### 10C. Final Validation Report

Create comprehensive report including:
- Algorithm validation results
- Production processing statistics
- Quality metrics distribution
- Performance benchmarks
- Known limitations
- Future improvements

**Deliverable**: Production-ready repository with complete documentation

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

### Milestone 3: Production Framework Ready ✨
- ✅ All production components implemented
- ✅ Integration tests pass
- ✅ 100-protein test batch processed via SLURM

### Milestone 4: Production Deployment ✨
- ✅ 40k proteins processed
- ✅ Quality metrics acceptable
- ✅ Results imported to database
- ✅ Documentation complete

## Quick Reference Commands

```bash
# Development
pip install -e ".[dev]"               # Install in dev mode
pytest tests/ -v                      # Run all tests
pytest tests/ -m regression           # Run regression tests only
mypy src/                             # Type checking
black src/ tests/                     # Format code
ruff check src/                       # Lint code

# Algorithm Testing
pyecod-mini 8ovp_A --verbose          # Test single protein
pyecod-mini --test-suite              # Run formal test suite
python scripts/validate_algorithm.py  # Comprehensive validation

# Production Processing
pyecod-mini production scan           # Scan for proteins
pyecod-mini production process        # Submit SLURM jobs
pyecod-mini production monitor        # Monitor progress
pyecod-mini production import         # Import to database

# Utilities
pyecod-mini --list-batches            # Show available batches
pyecod-mini --validate                # Validate configuration
pyecod-mini production stats          # Show processing stats
```

## Troubleshooting

### Common Issues

**Issue**: Tests fail with import errors
**Solution**: Ensure you're in project root and package is installed: `pip install -e .`

**Issue**: mypy reports type errors
**Solution**: Add type hints incrementally, use `# type: ignore` for complex cases temporarily

**Issue**: Regression tests fail
**Solution**: Compare with old mini output, may need to adjust paths or reference data

**Issue**: SLURM jobs fail
**Solution**: Check logs in `/tmp/pyecod_mini_logs/`, verify environment setup

## Getting Help

- Check planning documents (EXTRACTION_PLAN.md, etc.)
- Review CLAUDE.md for lessons learned
- Run with `--verbose` for detailed output
- Check logs and tracking database

---

**Ready to Begin!** Start with Step 1 and work through systematically.

**Estimated Timeline**: 14 days from start to production deployment
