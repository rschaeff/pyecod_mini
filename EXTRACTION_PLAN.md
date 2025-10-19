# PyECOD Mini Extraction & Setup Plan

**Date**: 2025-10-19
**Source**: `/home/rschaeff/dev/pyecod/mini`
**Target**: `/home/rschaeff/dev/pyecod_mini` (this repository)

## Executive Summary

This plan outlines the extraction of the proven `mini` implementation from the legacy pyecod repository into a new, clean `pyecod_mini` repository with:
1. **Modern Python practices** (dataclasses, type hints, mypy compliance)
2. **Validated algorithm** (6/6 regression tests passing, ~80% accuracy)
3. **New production framework** (SLURM-based processing for ~40k proteins)

## Phase 1: Repository Structure & Setup

### 1.1 Directory Structure

```
pyecod_mini/                          # This repository
├── pyproject.toml                    # Modern Python packaging
├── README.md                         # User-facing documentation
├── EXTRACTION_PLAN.md                # This file
├── ALGORITHM_VALIDATION.md           # Validation report
├── PRODUCTION_DESIGN.md              # Production framework design
│
├── src/
│   └── pyecod_mini/                  # Main package
│       ├── __init__.py
│       ├── __main__.py               # Enable `python -m pyecod_mini`
│       │
│       ├── core/                     # Core partitioning algorithm
│       │   ├── __init__.py
│       │   ├── models.py             # ALL dataclasses with type hints
│       │   ├── sequence_range.py     # Sequence range handling
│       │   ├── parser.py             # Evidence parsing
│       │   ├── partitioner.py        # Main partitioning logic
│       │   ├── decomposer.py         # Chain BLAST decomposition
│       │   ├── boundary_optimizer.py # Boundary optimization
│       │   ├── writer.py             # Output XML generation
│       │   ├── blast_parser.py       # BLAST XML parsing
│       │   ├── evidence_utils.py     # Evidence processing
│       │   ├── domain_utils.py       # Domain utilities
│       │   ├── gap_analyzer.py       # Gap analysis
│       │   ├── ecod_domains_parser.py # ECOD classification parsing
│       │   ├── range_cache_parser.py # Reference data parsing
│       │   └── visualization.py      # PyMOL visualization
│       │
│       ├── cli/                      # Command-line interface
│       │   ├── __init__.py
│       │   ├── main.py               # Main CLI entry point
│       │   ├── batch_finder.py       # Batch detection logic
│       │   └── config.py             # Configuration management
│       │
│       └── production/               # Production processing framework (NEW)
│           ├── __init__.py
│           ├── processor.py          # Main batch processor
│           ├── slurm_manager.py      # SLURM job management
│           ├── database.py           # Database integration
│           ├── monitor.py            # Progress monitoring
│           └── quality_filter.py     # Quality filtering
│
├── tests/                            # Test suite
│   ├── __init__.py
│   ├── conftest.py                   # pytest configuration
│   ├── test_core.py                  # Core algorithm tests
│   ├── test_ecod_regression.py       # Regression tests (6 cases)
│   ├── test_decomposer.py            # Decomposition tests
│   ├── test_discontinuous.py         # Multi-segment domain tests
│   ├── test_parser.py                # Evidence parsing tests
│   ├── test_writer.py                # Output generation tests
│   ├── test_models.py                # Data model tests
│   ├── test_cli.py                   # CLI tests
│   ├── test_batch_proteins.py        # Batch processing tests
│   ├── test_ecod_tgroup.py           # ECOD classification tests
│   └── test_production.py            # Production framework tests (NEW)
│
├── data/                             # Reference data
│   ├── domain_lengths.csv            # ECOD domain lengths
│   ├── protein_lengths.csv           # Protein chain lengths
│   ├── domain_definitions.csv        # Multi-domain decomposition rules
│   ├── reference_blacklist.csv       # Reference quality blacklist
│   └── ecod_classifications.csv      # ECOD T-group/H-group data
│
├── scripts/                          # Utility scripts
│   ├── setup_production.py           # Production setup wizard
│   ├── generate_references.py        # Generate reference data
│   └── validate_algorithm.py         # Algorithm validation
│
├── docs/                             # Documentation
│   ├── algorithm.md                  # Algorithm description
│   ├── production.md                 # Production usage guide
│   ├── api.md                        # API documentation
│   └── development.md                # Development guide
│
└── config/                           # Configuration templates
    ├── config.template.yml           # Main config template
    ├── slurm.template.yml            # SLURM config template
    └── database.template.yml         # Database config template
```

### 1.2 Modern Python Packaging

**pyproject.toml** (PEP 517/518 compliant):
```toml
[build-system]
requires = ["setuptools>=65.0", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "pyecod_mini"
version = "2.0.0"
description = "Minimal domain partitioning tool for ECOD protein classification"
authors = [{name = "Your Name", email = "your.email@example.com"}]
requires-python = ">=3.9"
dependencies = [
    "lxml>=4.9.0",
    "pyyaml>=6.0",
    "psycopg2-binary>=2.9.0",
    "requests>=2.31.0",
]

[project.optional-dependencies]
dev = [
    "pytest>=7.4.0",
    "pytest-cov>=4.1.0",
    "mypy>=1.5.0",
    "black>=23.7.0",
    "ruff>=0.0.284",
]
production = [
    "sqlalchemy>=2.0.0",
    "pandas>=2.0.0",
]

[project.scripts]
pyecod-mini = "pyecod_mini.cli.main:main"

[tool.setuptools.packages.find]
where = ["src"]

[tool.mypy]
python_version = "3.9"
warn_return_any = true
warn_unused_configs = true
disallow_untyped_defs = true

[tool.black]
line-length = 100
target-version = ['py39']

[tool.ruff]
line-length = 100
target-version = "py39"

[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = ["test_*.py"]
python_classes = ["Test*"]
python_functions = ["test_*"]
markers = [
    "unit: Unit tests",
    "integration: Integration tests",
    "regression: Regression tests",
    "slow: Slow-running tests",
    "production: Production framework tests",
]
```

### 1.3 Code Quality Standards

**Mandatory for all new/refactored code:**
1. ✅ **Type hints** on all function signatures and class attributes
2. ✅ **Dataclasses** for all structured data (no dicts)
3. ✅ **Docstrings** (Google style) for all public functions/classes
4. ✅ **mypy** type checking passes with no errors
5. ✅ **black** code formatting
6. ✅ **ruff** linting passes
7. ✅ **pytest** coverage >= 80%

## Phase 2: Core Algorithm Extraction

### 2.1 Files to Extract (with improvements)

**Priority 1: Core Algorithm (Proven & Working)**

| Source File | Target File | Refactoring Required |
|------------|-------------|---------------------|
| `mini/core/models.py` | `src/pyecod_mini/core/models.py` | ⚠️ High - convert dicts to dataclasses |
| `mini/core/sequence_range.py` | `src/pyecod_mini/core/sequence_range.py` | ✅ Low - already well-structured |
| `mini/core/parser.py` | `src/pyecod_mini/core/parser.py` | ⚠️ Medium - add type hints |
| `mini/core/partitioner.py` | `src/pyecod_mini/core/partitioner.py` | ⚠️ Medium - add type hints |
| `mini/core/decomposer.py` | `src/pyecod_mini/core/decomposer.py` | ⚠️ Medium - add type hints |
| `mini/core/boundary_optimizer.py` | `src/pyecod_mini/core/boundary_optimizer.py` | ⚠️ Medium - add type hints |
| `mini/core/writer.py` | `src/pyecod_mini/core/writer.py` | ⚠️ Medium - add type hints |
| `mini/core/blast_parser.py` | `src/pyecod_mini/core/blast_parser.py` | ⚠️ Medium - add type hints |
| `mini/core/evidence_utils.py` | `src/pyecod_mini/core/evidence_utils.py` | ⚠️ Medium - add type hints |
| `mini/core/domain_utils.py` | `src/pyecod_mini/core/domain_utils.py` | ⚠️ Low - minor cleanup |
| `mini/core/gap_analyzer.py` | `src/pyecod_mini/core/gap_analyzer.py` | ⚠️ Low - minor cleanup |
| `mini/core/ecod_domains_parser.py` | `src/pyecod_mini/core/ecod_domains_parser.py` | ⚠️ Medium - add type hints |
| `mini/core/range_cache_parser.py` | `src/pyecod_mini/core/range_cache_parser.py` | ⚠️ Medium - add type hints |
| `mini/core/visualization.py` | `src/pyecod_mini/core/visualization.py` | ⚠️ Low - optional feature |

**Priority 2: CLI & Configuration**

| Source File | Target File | Refactoring Required |
|------------|-------------|---------------------|
| `mini/pyecod_mini.py` | `src/pyecod_mini/cli/main.py` | ⚠️ High - split into modules |
| (extract) | `src/pyecod_mini/cli/batch_finder.py` | ⚠️ High - extract from pyecod_mini.py |
| (extract) | `src/pyecod_mini/cli/config.py` | ⚠️ High - extract from pyecod_mini.py |

**Priority 3: Test Suite**

| Source File | Target File | Status |
|------------|-------------|--------|
| `mini/tests/test_core.py` | `tests/test_core.py` | ✅ Copy with minor updates |
| `mini/tests/test_ecod_regression.py` | `tests/test_ecod_regression.py` | ✅ **Critical - 6 regression tests** |
| `mini/tests/test_decomposer.py` | `tests/test_decomposer.py` | ✅ Copy |
| `mini/tests/test_discontinuous.py` | `tests/test_discontinuous.py` | ✅ Copy |
| `mini/tests/test_parser.py` | `tests/test_parser.py` | ✅ Copy |
| `mini/tests/test_writer.py` | `tests/test_writer.py` | ✅ Copy |
| `mini/tests/test_models.py` | `tests/test_models.py` | ✅ Copy |
| `mini/tests/test_cli.py` | `tests/test_cli.py` | ⚠️ Update paths |
| `mini/tests/test_batch_proteins.py` | `tests/test_batch_proteins.py` | ✅ Copy |
| `mini/tests/test_ecod_tgroup.py` | `tests/test_ecod_tgroup.py` | ✅ Copy |
| `mini/tests/conftest.py` | `tests/conftest.py` | ⚠️ Update for new structure |

**Priority 4: Reference Data**

| Source File | Target File | Status |
|------------|-------------|--------|
| `mini/test_data/domain_lengths.csv` | `data/domain_lengths.csv` | ✅ Copy |
| `mini/test_data/protein_lengths.csv` | `data/protein_lengths.csv` | ✅ Copy |
| `mini/test_data/domain_definitions.csv` | `data/domain_definitions.csv` | ✅ Copy |
| `mini/test_data/reference_blacklist.csv` | `data/reference_blacklist.csv` | ✅ Copy (if exists) |
| (generate) | `data/ecod_classifications.csv` | ⚠️ Generate from ECOD database |

### 2.2 Files to Leave Behind

**DO NOT extract from the old repository:**
- ❌ Everything in `ecod/` (failed pipeline)
- ❌ `mini_production/` old scripts (will be redesigned)
- ❌ Old configuration files
- ❌ Legacy scripts that duplicate functionality
- ❌ Anything not proven to work

### 2.3 Refactoring Priority: Dataclasses & Type Hints

**Critical Insight from CLAUDE.md:**
> The main pyECOD suffered from inconsistent use of dataclasses vs dictionaries. This created technical debt and confusion.

**Refactoring Standard:**

**BEFORE (dict-based, common in old code):**
```python
def process_evidence(evidence: dict) -> dict:
    return {
        "type": evidence["type"],
        "range": evidence["range"],
        "confidence": evidence.get("confidence", 0.0)
    }
```

**AFTER (dataclass-based, mandatory for new repo):**
```python
from dataclasses import dataclass
from typing import Optional

@dataclass
class Evidence:
    """Evidence for domain assignment"""
    type: str
    range: SequenceRange
    confidence: float = 0.0
    evalue: Optional[float] = None

    def get_quality_score(self) -> float:
        """Calculate quality score from evidence metrics"""
        # Type-safe, autocomplete works, clear contract
        return self.confidence * (1.0 / max(self.evalue or 1.0, 1e-10))

def process_evidence(evidence: Evidence) -> Evidence:
    """Process evidence item with validation"""
    # Type hints enable IDE support and catch errors early
    return Evidence(
        type=evidence.type,
        range=evidence.range,
        confidence=min(evidence.confidence, 1.0)
    )
```

**Refactoring Checklist for Each Module:**
- [ ] Replace all dict-based data structures with dataclasses
- [ ] Add type hints to all functions (params + return types)
- [ ] Add docstrings (Google style)
- [ ] Run mypy and fix all type errors
- [ ] Run black to format code
- [ ] Update tests to match new signatures

## Phase 3: Algorithm Validation

### 3.1 Regression Test Suite

**The Critical 6 Tests** (from `test_ecod_regression.py`):

| Test Case | Protein | Expected Domains | Coverage | Classification Accuracy |
|-----------|---------|------------------|----------|------------------------|
| 1 | 8ovp_A | 3 domains | ~80% | GFP + PBP families |
| 2-6 | TBD | TBD | TBD | TBD |

**Validation Criteria:**
- ✅ All 6 regression tests pass
- ✅ Domain count matches expected (±1 acceptable)
- ✅ Coverage >= 75%
- ✅ ECOD T-group classification accuracy >= 80%
- ✅ No crashes or errors
- ✅ Performance < 10 seconds per protein

### 3.2 Extended Validation

**Additional Validation Tests:**
1. **Unit tests**: All core modules individually tested
2. **Integration tests**: Full pipeline end-to-end
3. **Edge cases**:
   - Single domain proteins
   - Multi-domain proteins (5+ domains)
   - Discontinuous domains
   - Overlapping evidence
   - Missing reference data
   - Empty evidence

**Performance Benchmarks:**
- Process 100 proteins in < 15 minutes
- Memory usage < 500MB per protein
- No memory leaks during batch processing

### 3.3 Validation Report

Create `ALGORITHM_VALIDATION.md` with:
- Regression test results
- Performance benchmarks
- Edge case handling
- Comparison with ECOD reference data
- Known limitations and future improvements

## Phase 4: Production Framework Design

### 4.1 Goals

**The new production framework should:**
1. ✅ **Simple**: No complex service architecture
2. ✅ **Proven**: Built on working mini algorithm
3. ✅ **Scalable**: SLURM integration for 40k+ proteins
4. ✅ **Reliable**: Robust error handling and recovery
5. ✅ **Observable**: Real-time monitoring and logging
6. ✅ **Maintainable**: Clean code, good documentation

### 4.2 Architecture Overview

```
┌─────────────────────────────────────────────────────────────┐
│                   Production Processor                       │
│                                                              │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐ │
│  │   Scanner    │───>│ Job Manager  │───>│   Monitor    │ │
│  │              │    │              │    │              │ │
│  │ - Find files │    │ - SLURM      │    │ - Progress   │ │
│  │ - Check DB   │    │ - Tracking   │    │ - Quality    │ │
│  │ - Filter     │    │ - Retry      │    │ - Alerts     │ │
│  └──────────────┘    └──────────────┘    └──────────────┘ │
│         │                    │                    │         │
│         v                    v                    v         │
│  ┌──────────────────────────────────────────────────────┐  │
│  │              Core Algorithm (pyecod_mini)           │  │
│  └──────────────────────────────────────────────────────┘  │
│         │                                          │         │
│         v                                          v         │
│  ┌──────────────┐                          ┌──────────────┐ │
│  │  XML Output  │                          │   Database   │ │
│  │              │                          │   Import     │ │
│  └──────────────┘                          └──────────────┘ │
└─────────────────────────────────────────────────────────────┘
```

### 4.3 Component Design

**Scanner Component** (`src/pyecod_mini/production/scanner.py`):
```python
@dataclass
class ScanResult:
    """Result of batch scanning"""
    total_proteins: int
    new_proteins: List[str]
    existing_proteins: List[str]
    failed_proteins: List[str]
    batch_metadata: Dict[str, Any]

class BatchScanner:
    """Scan batch directories for proteins to process"""

    def scan_batches(
        self,
        batch_dirs: List[Path],
        reps_only: bool = False,
        skip_existing: bool = True
    ) -> ScanResult:
        """Scan batch directories and identify proteins to process"""
        ...
```

**SLURM Manager** (`src/pyecod_mini/production/slurm_manager.py`):
```python
@dataclass
class SlurmJob:
    """SLURM job specification"""
    job_id: str
    protein_id: str
    batch_name: str
    script_path: Path
    status: str  # pending, running, completed, failed
    submit_time: datetime
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    exit_code: Optional[int] = None

class SlurmManager:
    """Manage SLURM job submission and tracking"""

    def submit_job(self, protein_id: str, batch_dir: Path) -> SlurmJob:
        """Submit a single protein processing job"""
        ...

    def monitor_jobs(self, jobs: List[SlurmJob]) -> Dict[str, int]:
        """Monitor job status and update tracking"""
        ...

    def retry_failed_jobs(self, jobs: List[SlurmJob], max_retries: int = 2) -> List[SlurmJob]:
        """Retry failed jobs with exponential backoff"""
        ...
```

**Database Integration** (`src/pyecod_mini/production/database.py`):
```python
@dataclass
class PartitionRecord:
    """Database record for partition result"""
    pdb_id: str
    chain_id: str
    process_version: str = "mini_pyecod_2.0"
    is_classified: bool = False
    domain_count: int = 0
    coverage_fraction: float = 0.0
    domains: List[Dict[str, Any]] = field(default_factory=list)
    timestamp: datetime = field(default_factory=datetime.now)

class DatabaseImporter:
    """Import mini results to PostgreSQL database"""

    def import_batch(
        self,
        results: List[PartitionRecord],
        collision_strategy: str = "skip"  # skip, overwrite, version
    ) -> ImportReport:
        """Import batch of results with collision handling"""
        ...
```

**Monitor** (`src/pyecod_mini/production/monitor.py`):
```python
@dataclass
class ProcessingStats:
    """Real-time processing statistics"""
    total_proteins: int
    completed: int
    failed: int
    pending: int
    running: int
    avg_time_per_protein: float
    estimated_completion: datetime

class ProductionMonitor:
    """Real-time monitoring of production processing"""

    def get_stats(self) -> ProcessingStats:
        """Get current processing statistics"""
        ...

    def watch(self, refresh_interval: int = 10) -> None:
        """Live monitoring dashboard in terminal"""
        ...
```

### 4.4 Configuration Management

**config.yml** structure:
```yaml
database:
  host: dione
  port: 45000
  database: ecod_protein
  user: ecod
  password_env: ECOD_DB_PASSWORD  # Read from environment

paths:
  batch_base_dir: /data/ecod/pdb_updates/batches
  output_base_dir: /data/ecod/mini_v2_results
  log_dir: /var/log/pyecod_mini

slurm:
  partition: All
  time: "1:00:00"
  memory: 4G
  cpus_per_task: 1
  max_concurrent_jobs: 50
  job_script_dir: /tmp/pyecod_mini_jobs
  log_dir: /tmp/pyecod_mini_logs

processing:
  reps_only: true
  skip_existing: true
  max_retries: 2
  retry_delay: 300  # seconds

quality:
  min_coverage: 0.5
  min_confidence: 0.6
  require_reference_lengths: true
```

### 4.5 CLI for Production

```bash
# Scan and submit jobs
pyecod-mini production scan --batch-dir /data/ecod/pdb_updates/batches/batch_036

# Submit all representative proteins
pyecod-mini production process --reps-only --max-jobs 50

# Monitor progress
pyecod-mini production monitor --watch

# Import results to database
pyecod-mini production import --check-collisions --dry-run
pyecod-mini production import --collision-strategy skip

# Quality report
pyecod-mini production quality --output report.json

# Retry failed jobs
pyecod-mini production retry --max-retries 2
```

## Phase 5: Migration Checklist

### 5.1 Pre-Migration

- [ ] Create new repository structure
- [ ] Set up pyproject.toml with dependencies
- [ ] Configure mypy, black, ruff, pytest
- [ ] Create initial README and documentation stubs

### 5.2 Core Algorithm Migration

- [ ] Extract and refactor core modules (Priority 1)
- [ ] Add type hints to all extracted modules
- [ ] Convert dicts to dataclasses in models.py
- [ ] Run mypy and fix all type errors
- [ ] Run black to format all code
- [ ] Extract and update tests
- [ ] Run test suite - ensure 6 regression tests pass

### 5.3 CLI Migration

- [ ] Extract CLI components from pyecod_mini.py
- [ ] Create modular CLI structure
- [ ] Update configuration management
- [ ] Test batch finder functionality
- [ ] Update CLI tests

### 5.4 Production Framework

- [ ] Implement Scanner component
- [ ] Implement SLURM Manager
- [ ] Implement Database Importer
- [ ] Implement Monitor
- [ ] Create configuration templates
- [ ] Write production tests
- [ ] Document production usage

### 5.5 Validation & Testing

- [ ] Run full test suite
- [ ] Validate 6 regression tests pass
- [ ] Run on 100 test proteins
- [ ] Performance benchmarks
- [ ] Memory leak testing
- [ ] Create validation report

### 5.6 Documentation

- [ ] Update README with quick start
- [ ] Write algorithm documentation
- [ ] Write production usage guide
- [ ] Write API documentation
- [ ] Write development guide
- [ ] Document known limitations

### 5.7 Production Deployment

- [ ] Set up production configuration
- [ ] Test on small batch (10 proteins)
- [ ] Test on medium batch (100 proteins)
- [ ] Full production run (40k proteins)
- [ ] Monitor and validate results
- [ ] Import to database
- [ ] Generate quality report

## Phase 6: Key Lessons Applied

### From CLAUDE.md:

1. ✅ **Algorithm First, Infrastructure Second**
   - Core algorithm already proven (6/6 tests passing)
   - Production framework built around proven algorithm
   - No complex architecture before validation

2. ✅ **Consistent Dataclass Usage**
   - All structured data as dataclasses
   - No dict/dataclass confusion
   - Type hints everywhere
   - mypy compliance from day 1

3. ✅ **Test-Driven Validation**
   - Regression tests validate correctness
   - Performance tests validate scalability
   - Quality tests validate production readiness

4. ✅ **Simple, Maintainable Code**
   - Clean separation of concerns
   - Modular design
   - Well-documented
   - Easy to extend

5. ✅ **Production Readiness**
   - Robust error handling
   - Monitoring and observability
   - Database integration
   - SLURM scaling

## Success Criteria

**Phase 1-3 (Algorithm Extraction & Validation):**
- ✅ All 6 regression tests pass
- ✅ mypy type checking passes
- ✅ pytest coverage >= 80%
- ✅ Code formatted with black
- ✅ No ruff linting errors
- ✅ Performance benchmarks met

**Phase 4-5 (Production Framework):**
- ✅ Successfully process 100 test proteins
- ✅ SLURM integration working
- ✅ Database import working
- ✅ Monitoring dashboard functional
- ✅ Quality filtering effective

**Phase 6 (Production Deployment):**
- ✅ Successfully process 40k representative proteins
- ✅ Results imported to database
- ✅ Quality metrics within acceptable ranges
- ✅ No data corruption or collisions
- ✅ Production documentation complete

## Timeline Estimate

- **Phase 1-2 (Structure & Extraction)**: 2-3 days
- **Phase 3 (Validation)**: 1-2 days
- **Phase 4 (Production Framework)**: 3-4 days
- **Phase 5 (Testing & Documentation)**: 2-3 days
- **Phase 6 (Production Deployment)**: 1-2 days

**Total**: 9-14 days for complete extraction, validation, and production deployment

## Next Steps

1. **Immediate**: Create repository structure and pyproject.toml
2. **Next**: Extract core algorithm with type hints
3. **Then**: Run and validate regression tests
4. **Finally**: Build production framework

---

**Note**: This plan is a living document. Update as we progress through the migration.
