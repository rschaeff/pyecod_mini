# pyECOD Mini - Clean Domain Partitioning Tool

**Version**: 2.0.0
**Status**: ✅ Production Ready - Integrated with pyecod_prod Framework

A minimal, validated domain partitioning tool for ECOD protein classification with production-ready scaling.

## Overview

pyECOD Mini is a clean extraction of the proven domain partitioning algorithm from the legacy pyECOD project. This repository represents a fresh start with:

- ✅ **Validated Algorithm**: 6/6 regression tests passing, ~80% domain boundary accuracy
- ✅ **Modern Python**: Type hints, dataclasses, mypy compliance
- ✅ **Production Ready**: Integrated with pyecod_prod for large-scale processing
- ✅ **Clean Architecture**: Simple, maintainable, well-documented
- ✅ **External Workflow Support**: CLI arguments for integration with batch processing systems

## Quick Start

### Installation

```bash
# Clone repository
git clone git@github.com:rschaeff/pyecod_mini.git
cd pyecod_mini

# Install in development mode
pip install -e . --user

# Verify installation
pyecod-mini --validate
```

### Basic Usage

**Partition a single protein from a batch:**
```bash
pyecod-mini 8ovp_A --batch-id ecod_weekly_20250905 --verbose
```

**Use custom input/output paths (for integration with external workflows):**
```bash
pyecod-mini 8s72_N \
    --summary-xml /path/to/domains/8s72_N.develop291.domain_summary.xml \
    --output /path/to/partitions/8s72_N.domains.xml \
    --verbose
```

**Generate PyMOL visualization:**
```bash
pyecod-mini 8ovp_A --batch-id ecod_weekly_20250905 --visualize
pymol /data/ecod/pdb_updates/batches/ecod_weekly_20250905/comparison_8ovp_A.pml
```

### Current Status

✅ **Production Ready** - Fully integrated with pyecod_prod framework

The algorithm currently:
- ✅ Partitions domains with ~80% accuracy
- ✅ Passes 6 regression tests
- ✅ Handles discontinuous domains
- ✅ Optimizes domain boundaries
- ✅ Produces ECOD-compliant XML output
- ✅ Integrates with batch processing workflows
- ✅ Processes 15-chain test batch successfully (100%)

## Goals

### Phase 1: Algorithm Extraction ✨ Current Phase

Extract the proven mini algorithm with improvements:

1. **Modern Python Standards**
   - All data structures as dataclasses (no dicts)
   - Comprehensive type hints (mypy compliant)
   - Black code formatting
   - Ruff linting

2. **Validated Correctness**
   - All 6 regression tests passing
   - 80% test coverage
   - Performance benchmarks met

3. **Clean Architecture**
   - Modular design
   - Clear separation of concerns
   - Comprehensive documentation

### Phase 2: Production Framework (Next)

Build scalable production processing:

1. **SLURM Integration**
   - Parallel processing of 40k+ proteins
   - Job tracking and monitoring
   - Automatic retry on failure

2. **Database Integration**
   - Import results to PostgreSQL
   - Quality filtering
   - Collision detection

3. **Observability**
   - Real-time progress monitoring
   - Quality metrics tracking
   - Performance dashboards

## Key Features

### Core Algorithm

- **Evidence Integration**: Combine BLAST and HHsearch evidence
- **Chain Decomposition**: Handle multi-domain chains via BLAST alignment
- **Boundary Optimization**: Refine domain boundaries using alignment quality
- **Discontinuous Domains**: Support multi-segment domains
- **Provenance Tracking**: Complete audit trail of decisions

### CLI Arguments

**Batch Mode (Default):**
```bash
pyecod-mini PROTEIN_ID [--batch-id BATCH_ID] [--verbose] [--visualize]
```

- `PROTEIN_ID`: Protein to partition (e.g., `8ovp_A`, `8s72_N`)
- `--batch-id`: Optional batch ID for batch detection
- `--verbose`: Show detailed processing information
- `--visualize`: Generate PyMOL comparison script

**Integration Mode (Custom Paths):**
```bash
pyecod-mini PROTEIN_ID \
    --summary-xml PATH_TO_DOMAIN_SUMMARY_XML \
    --output PATH_TO_OUTPUT_XML \
    [--verbose]
```

- `--summary-xml`: Path to input domain summary XML (overrides batch detection)
- `--output`: Path to output partition XML (overrides batch detection)

This mode enables integration with external batch processing workflows like `pyecod_prod`.

### Integration with pyecod_prod

pyecod-mini is designed to integrate seamlessly with the `pyecod_prod` batch processing framework:

1. **pyecod_prod** runs BLAST and HHsearch, generates domain summaries
2. **pyecod_prod** calls `pyecod-mini` with `--summary-xml` and `--output` for each chain
3. **pyecod-mini** partitions domains and writes output XML with provenance metadata
4. **pyecod_prod** tracks completion in batch manifest

**Example Integration Call:**
```python
from pyecod_prod.partition.partition_runner import PartitionRunner

runner = PartitionRunner(pyecod_mini_path="/home/user/.local/bin/pyecod-mini")
runner.partition_protein(
    pdb_id="8s72",
    chain_id="N",
    summary_xml="/data/ecod/batches/ecod_weekly_20250905/domains/8s72_N.develop291.domain_summary.xml",
    output_path="/data/ecod/batches/ecod_weekly_20250905/partitions/8s72_N.domains.xml"
)
```

### Production Processing

- **Batch Scanning**: Identify proteins to process
- **External Workflow Integration**: Custom file paths via CLI
- **Progress Monitoring**: Real-time dashboard in pyecod_prod
- **Quality Control**: Automatic quality assessment
- **Database Import**: Safe import with collision detection

## Project Structure

```
pyecod_mini/
├── src/pyecod_mini/          # Main package
│   ├── core/                 # Core algorithm
│   ├── cli/                  # Command-line interface
│   └── production/           # Production framework
├── tests/                    # Test suite
├── data/                     # Reference data
├── scripts/                  # Utility scripts
├── docs/                     # Documentation
└── config/                   # Configuration templates
```

See [EXTRACTION_PLAN.md](EXTRACTION_PLAN.md) for detailed structure.

## Development

### Prerequisites

- Python >= 3.9
- PostgreSQL (for production database)
- SLURM cluster (for production processing)

### Setup

```bash
# Clone repository
git clone git@github.com:rschaeff/pyecod_mini.git
cd pyecod_mini

# Install in development mode
pip install -e . --user

# Validate installation
pyecod-mini --validate

# Run tests
pytest tests/

# Type checking
mypy src/

# Format code
black src/ tests/
```

## Testing

### Regression Tests

The 6 critical regression tests validate algorithm correctness:

```bash
pytest tests/test_ecod_regression.py -v
```

**Test Case 1: 8ovp_A** (GFP-PBP fusion protein)
- Expected: 3 domains (or 2 if discontinuous merged)
- Coverage: >= 80%
- Classification: GFP + PBP families correctly identified

### Performance Tests

```bash
pytest tests/ -m performance
```

Target metrics:
- Single protein: < 10 seconds
- 100 proteins: < 15 minutes
- Memory: < 500 MB per protein

## Documentation

- [EXTRACTION_PLAN.md](EXTRACTION_PLAN.md) - Repository setup and extraction plan
- [ALGORITHM_VALIDATION.md](ALGORITHM_VALIDATION.md) - Validation strategy and test plan
- [PRODUCTION_DESIGN.md](PRODUCTION_DESIGN.md) - Production framework architecture
- [CLAUDE.md](../pyecod/CLAUDE.md) - Lessons learned from original pyECOD

## Key Lessons Applied

From the original pyECOD failure:

1. ✅ **Algorithm First, Infrastructure Second**
   - Validate algorithm before building complex infrastructure
   - Use proven algorithm as foundation

2. ✅ **Consistent Data Structures**
   - All structured data as dataclasses
   - No dict/dataclass confusion
   - Type safety throughout

3. ✅ **Test-Driven Validation**
   - Regression tests prove correctness
   - Performance tests ensure scalability
   - Quality tests validate production readiness

4. ✅ **Simple > Complex**
   - No premature abstraction
   - Clear, maintainable code
   - Avoid over-engineering

## Roadmap

### Phase 1: Extraction & Validation ✅ COMPLETED
- [x] Create repository structure
- [x] Design extraction plan
- [x] Design validation plan
- [x] Design production framework
- [x] Extract core algorithm with type hints
- [x] Migrate test suite
- [x] Run regression tests (6/6 passing)
- [x] Performance benchmarks

### Phase 2: Production Framework ✅ COMPLETED
- [x] Integration with pyecod_prod
- [x] CLI arguments for external workflows
- [x] Custom path support (--summary-xml, --output)
- [x] Batch processing integration
- [x] End-to-end validation (15/15 chains, 100% success)

### Phase 3: Production Deployment 🚧 IN PROGRESS
- [x] Process 15-chain test batch (100% success)
- [ ] Process 100-protein test batch
- [ ] Process 1000-protein staging batch
- [ ] Full production run (40k+ proteins)
- [ ] Quality assessment
- [ ] Database import
- [ ] Documentation finalization

## Contributing

This is an internal research project. Development follows:

- **Code Style**: Black (line length 100)
- **Type Checking**: mypy strict mode
- **Testing**: pytest with >= 80% coverage
- **Documentation**: Google-style docstrings

## License

Internal research project - not for public distribution.

## Authors

Based on the original pyECOD Mini implementation, redesigned and refactored for production use.

## Acknowledgments

- Original mini algorithm developers
- ECOD database team
- Lessons learned from the failed main pyECOD pipeline (see ../pyecod/CLAUDE.md)

---

**Current Status**: ✅ Production Ready - Integrated with pyecod_prod

- Algorithm: 6/6 regression tests passing
- Integration: End-to-end validation complete (15/15 chains, 100%)
- CLI: Enhanced with --summary-xml and --output arguments
- Next: Scale to full production batch (40k+ proteins)
