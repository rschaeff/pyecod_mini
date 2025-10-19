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

## [2.0.0] - TBD

### Major Changes
- Complete rewrite and extraction from legacy pyECOD repository
- Modern Python packaging with pyproject.toml
- Consistent dataclass usage throughout (no dict confusion)
- Comprehensive type hints (mypy compliant)
- Production-ready SLURM framework

### Added
- Core domain partitioning algorithm (proven, 6/6 regression tests)
- Boundary optimization with provenance tracking
- Production processing framework (Scanner, SLURM Manager, Monitor)
- Database import with collision detection
- Real-time monitoring dashboard

### Changed
- Migrated from dict-based to dataclass-based data structures
- Added type hints to all functions
- Modernized CLI interface
- Improved error handling and logging

### Fixed
- Inconsistent data structure usage
- Missing type annotations
- Unclear provenance tracking

## [1.0.0] - 2024 (Legacy Mini)

### Legacy Version
- Original mini implementation in ../pyecod/mini
- Proven algorithm with 6/6 regression tests
- ~80% domain boundary accuracy
- Basic SLURM integration via mini_production
- Issues: inconsistent dict/dataclass usage, missing type hints

---

**Note**: Version 2.0.0 represents a clean extraction and modernization of the proven 1.0.0 mini algorithm.
