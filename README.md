# pyECOD Mini - Clean Domain Partitioning Tool

**Version**: 2.0.0
**Status**: 🚧 In Development - Algorithm Extraction Phase

A minimal, validated domain partitioning tool for ECOD protein classification with production-ready scaling.

## Overview

pyECOD Mini is a clean extraction of the proven domain partitioning algorithm from the legacy pyECOD project. This repository represents a fresh start with:

- ✅ **Validated Algorithm**: 6/6 regression tests passing, ~80% domain boundary accuracy
- ✅ **Modern Python**: Type hints, dataclasses, mypy compliance
- ✅ **Production Ready**: SLURM integration for processing 40k+ proteins
- ✅ **Clean Architecture**: Simple, maintainable, well-documented

## Quick Start

**Note**: This repository is currently under development. The algorithm is being extracted from `/home/rschaeff/dev/pyecod/mini`.

### Current Status

We are in the **extraction and setup phase**. See planning documents:

- 📋 [EXTRACTION_PLAN.md](EXTRACTION_PLAN.md) - Complete extraction roadmap
- ✅ [ALGORITHM_VALIDATION.md](ALGORITHM_VALIDATION.md) - Validation strategy
- 🏭 [PRODUCTION_DESIGN.md](PRODUCTION_DESIGN.md) - Production framework design

### What Works (in source repo)

The mini algorithm in `/home/rschaeff/dev/pyecod/mini` currently:
- ✅ Partitions domains with ~80% accuracy
- ✅ Passes 6 regression tests
- ✅ Handles discontinuous domains
- ✅ Optimizes domain boundaries
- ✅ Produces ECOD-compliant XML output

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

### Production Processing

- **Batch Scanning**: Identify proteins to process
- **SLURM Jobs**: Parallel cluster processing
- **Progress Monitoring**: Real-time dashboard
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

### Setup (Coming Soon)

```bash
# Clone repository
git clone <repo-url>
cd pyecod_mini

# Install in development mode
pip install -e ".[dev]"

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

### Phase 1: Extraction & Validation (Weeks 1-2)
- [x] Create repository structure
- [x] Design extraction plan
- [x] Design validation plan
- [x] Design production framework
- [ ] Extract core algorithm with type hints
- [ ] Migrate test suite
- [ ] Run regression tests
- [ ] Performance benchmarks

### Phase 2: Production Framework (Weeks 3-4)
- [ ] Implement Scanner component
- [ ] Implement SLURM Manager
- [ ] Implement Monitoring dashboard
- [ ] Implement Database Importer
- [ ] Integration testing

### Phase 3: Production Deployment (Week 5)
- [ ] Process 100-protein test batch
- [ ] Process 1000-protein staging batch
- [ ] Full production run (40k proteins)
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

**Current Status**: 📋 Planning Complete - Ready for Algorithm Extraction

See [EXTRACTION_PLAN.md](EXTRACTION_PLAN.md) for next steps.
