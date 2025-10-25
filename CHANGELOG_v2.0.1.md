# pyecod_mini v2.0.1 - HHsearch Evidence Bug Fix

**Release Date**: 2025-10-25  
**Priority**: CRITICAL

## Overview

This release fixes a critical bug where HHsearch evidence was being parsed but immediately rejected during validation, resulting in 0% utilization of HHsearch results in domain partitioning.

## Critical Bug Fix

### Issue #2: HHsearch Evidence Not Used in Partitioning

**Symptoms**:
- HHsearch hits present in summary XML files
- Parser reported "Evidence validation failed: HHsearch evidence missing reference coverage data"
- Result: 0 domains found despite valid HHsearch evidence

**Root Cause**:
- Evidence validation requires `reference_coverage` (calculated from `hit_range` + `reference_length`)
- HHsearch parser failed to extract target range in TWO locations:
  1. API spec format: Missing `target_range` attribute extraction
  2. Legacy format: Missing `<hit_reg>` element extraction

**Fix**:
- **File**: `src/pyecod_mini/core/parser.py`
- **API spec HHsearch** (lines 280-348):
  - Added `target_range_str = hit.get("target_range", "")`
  - Parse to `hit_range` and pass to `populate_evidence_provenance()`
- **Legacy HHsearch** (lines 512-590):
  - Added `hit_reg = hit.find("hit_reg")`
  - Parse to `hit_range` and pass to `populate_evidence_provenance()`

**Impact**:
```
Test Case: 8axb_A (Ribonuclease H-like domain)

BEFORE FIX:
  Evidence parsed: 0
  Domains found: 0
  Coverage: 0%
  is_classified: false

AFTER FIX:
  Evidence parsed: 1 (HHsearch, e5b43A3, prob=98.4%)
  Domains found: 1
  Coverage: 44.7% (192/430 residues)
  Reference coverage: 61.9%
  is_classified: true
```

## Testing

### New Tests Added
- `TestHHsearchAPISpecFormat` class in `tests/test_parser.py`:
  - `test_parse_hhsearch_api_spec_format()`: Validates API spec format parsing
  - `test_hhsearch_without_target_range()`: Edge case testing

### Regression Tests
All formats verified working:
- ✅ Legacy BLAST format (`<blast_run>`)
- ✅ Legacy HHsearch format (`<hh_run>`)
- ✅ API spec BLAST format (`type="domain_blast"`)
- ✅ API spec HHsearch format (`type="hhsearch"`)
- ✅ Mixed evidence types

## Production Impact

### For 2023-2025 PDB Backfill
- **Affected chains**: 4,038 chains with HHsearch evidence
- **Previous coverage**: 0% HHsearch contribution (bug caused rejection)
- **Expected improvement**: Significant increase in domain coverage
- **Recovery time**: ~30 minutes to re-partition 3,920 chains

### Immediate Actions Required
1. ✅ Update pyecod_mini to v2.0.1
2. 🔄 Re-run partitioning for HHsearch-enabled chains
3. 🔄 Re-analyze coverage statistics
4. 🔄 Update production metrics

## Files Changed

### Modified
- `src/pyecod_mini/core/parser.py` (+14 lines)
  - API spec HHsearch: Added target_range parsing
  - Legacy HHsearch: Added hit_reg parsing
- `tests/test_parser.py` (+89 lines)
  - New test class for HHsearch API spec format
  - Comprehensive validation tests

### Added
- `CHANGELOG_v2.0.1.md` (this file)

## Upgrade Instructions

```bash
# Update from v2.0.0 to v2.0.1
cd /home/rschaeff/dev/pyecod_mini
git pull origin main

# Verify version
python3 -c "import pyecod_mini; print(pyecod_mini.__version__)"  # Should show 2.0.1

# Test with known case
python3 -m pyecod_mini.cli.main 8axb_A \
  --summary-xml /data/ecod/pdb_updates/backfill_2023_2025/blast/summaries_with_hhsearch/8axb_A.summary.xml \
  --output /tmp/test_8axb_A.partition.xml

# Should show: "HHsearch results: 1 domains"
```

## Technical Details

### Evidence Validation Flow
1. Parser extracts evidence from XML
2. `populate_evidence_provenance()` calculates `reference_coverage`
3. `validate_evidence_provenance()` checks for complete data
4. **Previous behavior**: Rejected if `reference_coverage = None`
5. **New behavior**: Calculates `reference_coverage` from `hit_range` + `reference_length`

### Reference Coverage Calculation
```python
reference_coverage = hit_range.total_length / reference_length

Example (8axb_A):
  hit_range: 16-161 (146 residues)
  reference_length: 236 residues
  coverage: 146/236 = 61.9% ✓
```

## Notes

- No breaking changes to API or CLI
- Fully backward compatible with v2.0.0
- No migration required for existing code
- Only partitioning results will change (improve) after update

---

**Questions?** Contact: rschaeff  
**Related Issues**: PYECOD_MINI_ISSUES.md (Issue #2)
