# Changelog

All notable changes to pyECOD Mini will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planning Phase
- Created comprehensive extraction plan (EXTRACTION_PLAN.md)
- Designed algorithm validation strategy (ALGORITHM_VALIDATION.md)
- Designed production framework architecture (PRODUCTION_DESIGN.md)
- Created repository structure with modern Python packaging
- Set up pyproject.toml with mypy, black, ruff, pytest configurations

## [2.0.0] - 2025-10-19

### Major Changes
- **Library API**: Stable programmatic API for integration with pyecod_prod
- **Version Tracking**: Algorithm version embedded in all outputs
- **API Spec Compliance**: Formal contract via PYECOD_MINI_API_SPEC.md
- Complete separation from production infrastructure (now in pyecod_prod)

### Added
- **Library API** (`api.py`):
  - `partition_protein()` - Main partitioning function
  - `PartitionResult` - Structured result dataclass
  - `PartitionError` - Partition-specific exception
  - `Domain` - API domain dataclass (simpler than internal model)
- **CLI Version Support**: `pyecod-mini --version` displays package version
- **Version Tracking**:
  - Package version (`__version__ = "2.0.0"`)
  - Writer prioritizes package version over git version
  - Algorithm version in partition XML metadata
- **Enhanced Provenance**:
  - Reference coverage tracking
  - Evidence quality thresholds
  - Complete boundary optimization audit trails
- **Quality Assessment**:
  - HHsearch probability reporting
  - Reference coverage calculation
  - Evidence quality flags

### Changed
- **Architecture**: Separated algorithm (this repo) from production infrastructure (pyecod_prod)
- **Writer**: `get_git_version()` now returns package version for reproducibility
- **Exports**: Clean `__all__` in `__init__.py` with library API components
- **Documentation**: Updated for library usage and API specification

### Fixed
- Version tracking now uses package version (not git tags)
- Evidence standardization before processing
- Reference coverage calculation for quality assessment

### Migration Notes
- **Breaking**: Production workflows moved to pyecod_prod repository
- **New**: Use library API for programmatic access:
  ```python
  from pyecod_mini import partition_protein
  result = partition_protein(summary_xml, output_xml, pdb_id, chain_id)
  ```
- **CLI**: Unchanged for backward compatibility

## [1.0.0] - 2024 (Legacy Mini)

### Legacy Version
- Original mini implementation in ../pyecod/mini
- Proven algorithm with 6/6 regression tests
- ~80% domain boundary accuracy
- Basic SLURM integration via mini_production
- Issues: inconsistent dict/dataclass usage, missing type hints

---

**Note**: Version 2.0.0 represents a clean extraction and modernization of the proven 1.0.0 mini algorithm.
