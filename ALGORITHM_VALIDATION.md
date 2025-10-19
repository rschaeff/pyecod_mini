# Algorithm Validation Plan

**Date**: 2025-10-19
**Algorithm**: pyECOD Mini Domain Partitioning
**Version**: 2.0.0

## Overview

This document outlines the comprehensive validation strategy for the pyECOD Mini domain partitioning algorithm to ensure production readiness.

## Validation Goals

1. ✅ **Correctness**: Algorithm produces accurate domain boundaries
2. ✅ **Consistency**: Reproducible results across runs
3. ✅ **Performance**: Acceptable speed and memory usage
4. ✅ **Robustness**: Handles edge cases gracefully
5. ✅ **Quality**: ECOD classification accuracy meets thresholds

## Test Pyramid

```
                  ┌─────────────┐
                  │  Production │  (~100 proteins)
                  │    Tests    │
                  └─────────────┘
                 ┌───────────────┐
                 │  Integration  │  (End-to-end workflows)
                 │     Tests     │
                 └───────────────┘
              ┌─────────────────────┐
              │   Regression Tests  │  (6 validated cases)
              └─────────────────────┘
           ┌──────────────────────────┐
           │      Unit Tests          │  (Individual functions)
           └──────────────────────────┘
```

## Level 1: Unit Tests

### Core Components to Test

**1. Sequence Range Handling** (`test_sequence_range.py`)
```python
def test_continuous_range():
    """Test simple continuous range"""
    range = SequenceRange.from_string("10-50")
    assert range.start == 10
    assert range.end == 50
    assert range.total_length == 41
    assert not range.is_discontinuous

def test_discontinuous_range():
    """Test multi-segment range"""
    range = SequenceRange.from_string("10-30,50-70")
    assert len(range.segments) == 2
    assert range.total_length == 42
    assert range.is_discontinuous
```

**2. Evidence Parsing** (`test_parser.py`)
```python
def test_parse_blast_evidence():
    """Test BLAST evidence parsing"""
    # Test parsing valid BLAST XML
    # Test coverage calculation
    # Test reference length lookup

def test_parse_hhsearch_evidence():
    """Test HHsearch evidence parsing"""
    # Test probability conversion
    # Test T-group extraction
```

**3. Domain Partitioning** (`test_partitioner.py`)
```python
def test_simple_partition():
    """Test partitioning with non-overlapping evidence"""

def test_overlapping_resolution():
    """Test overlap resolution by confidence"""

def test_gap_handling():
    """Test handling of unassigned regions"""
```

**4. Boundary Optimization** (`test_boundary_optimizer.py`)
```python
def test_small_gap_merging():
    """Test merging small gaps with adjacent domains"""

def test_interstitial_assignment():
    """Test assignment of interstitial regions"""
```

**Target Coverage**: >= 80% for all core modules

## Level 2: Regression Tests

### The Critical 6 Test Cases

These are real proteins with validated ECOD classifications that the algorithm must handle correctly.

#### Test Case 1: 8ovp_A (GFP-PBP Fusion)

**Expected Outcome:**
```yaml
protein_id: 8ovp_A
expected_domains: 3
domain_ranges:
  - range: "2-248"
    family: "2vha/2ia4"
    t_group: "PBP-like"
    classification: "correct"
  - range: "253-499"
    family: "6dgv"
    t_group: "GFP-like"
    classification: "correct"
  - range: "500-517"
    family: "2vha/2ia4"
    t_group: "PBP-like"
    classification: "correct"
coverage_target: ">= 80%"
discontinuous: true  # Domain 1+3 form discontinuous domain
```

**Validation Criteria:**
- ✅ Finds exactly 3 domains (or 2 if discontinuous merged)
- ✅ Coverage >= 80% of sequence
- ✅ Correctly identifies GFP and PBP families
- ✅ Handles discontinuous PBP domain
- ✅ Boundary accuracy within ±5 residues

**Known Challenges:**
- PBP domain split by GFP insertion
- Requires proper discontinuous domain handling
- Chain BLAST decomposition needed

#### Test Cases 2-6: TBD

**Selection Criteria for Additional Test Cases:**
1. Single domain protein (baseline case)
2. Multi-domain protein (4-5 domains)
3. All continuous domains
4. Weak evidence (low E-values)
5. Remote homology (HHsearch-only)

**Process:**
1. Select candidate proteins from validated batches
2. Manual review of ECOD classifications
3. Run mini algorithm
4. Verify results against ECOD reference
5. Document expected outcomes
6. Add to regression test suite

### Regression Test Execution

**Automated Test:**
```python
@pytest.mark.regression
def test_8ovp_A_regression():
    """Regression test for 8ovp_A"""
    result = run_mini_algorithm("8ovp_A", batch_dir)

    # Critical assertions
    assert result.success, "Algorithm failed"
    assert len(result.domains) >= 2, "Too few domains found"
    assert result.coverage_fraction >= 0.80, "Coverage too low"

    # ECOD classification
    families = [d.family for d in result.domains]
    assert "6dgv" in families or "GFP" in str(families), "GFP domain not found"
    assert "2vha" in families or "2ia4" in families, "PBP domain not found"
```

**Manual Validation:**
- Visual inspection of domain boundaries
- Comparison with PyMOL structure visualization
- Review of evidence quality
- Check for biological plausibility

**Success Criteria:**
- ✅ All 6 regression tests pass
- ✅ Domain count within expected range (±1)
- ✅ Coverage >= target for each case
- ✅ ECOD family classification correct
- ✅ No algorithm crashes or errors

## Level 3: Integration Tests

### End-to-End Workflows

**Test 1: Complete Pipeline**
```python
def test_full_pipeline_8ovp_A():
    """Test complete pipeline from XML to output"""
    # Parse domain summary XML
    # Load reference data
    # Load BLAST alignments
    # Partition domains
    # Optimize boundaries
    # Write output XML
    # Validate output format
```

**Test 2: Batch Processing**
```python
def test_batch_processing():
    """Test processing multiple proteins"""
    proteins = ["8ovp_A", "protein2", "protein3"]
    results = process_batch(proteins)

    assert all(r.success for r in results)
    assert len(results) == len(proteins)
```

**Test 3: Error Recovery**
```python
def test_missing_reference_data():
    """Test graceful handling of missing data"""

def test_corrupt_xml_input():
    """Test handling of malformed input"""

def test_no_evidence_found():
    """Test handling of proteins with no evidence"""
```

## Level 4: Production Tests

### Protein Diversity Testing

**Test on 100 diverse proteins:**
- 20 single-domain proteins
- 40 multi-domain proteins (2-3 domains)
- 20 complex proteins (4+ domains)
- 10 discontinuous domain cases
- 10 edge cases (very long, very short, weak evidence)

**Metrics to Track:**
```yaml
success_rate: ">= 95%"
average_coverage: ">= 75%"
average_domains_per_protein: "2-3"
classification_accuracy: ">= 80%"
```

### Performance Benchmarks

**Speed Targets:**
- Single protein: < 10 seconds
- 100 proteins: < 15 minutes
- Average throughput: >= 400 proteins/hour

**Memory Targets:**
- Peak memory per protein: < 500 MB
- No memory leaks (constant memory for batch)

**Test Execution:**
```python
def test_performance_benchmark():
    """Benchmark processing speed"""
    start = time.time()
    results = process_batch(test_proteins_100)
    duration = time.time() - start

    avg_time = duration / len(test_proteins_100)
    assert avg_time < 10.0, f"Too slow: {avg_time}s per protein"
```

### Scalability Testing

**SLURM Integration:**
- Submit 1000 jobs
- Monitor completion
- Verify no job failures
- Check output quality

## Validation Metrics

### Correctness Metrics

1. **Domain Count Accuracy**
   - Compare with ECOD reference
   - Tolerance: ±1 domain acceptable
   - Target: >= 80% exact match

2. **Boundary Accuracy**
   - Compare domain start/end with ECOD
   - Tolerance: ±5 residues
   - Target: >= 85% within tolerance

3. **Coverage Accuracy**
   - Fraction of sequence assigned to domains
   - Target: >= 75% average coverage
   - Minimum: >= 50% per protein

4. **Classification Accuracy**
   - ECOD T-group/H-group matching
   - Target: >= 80% correct classification
   - Method: Compare with ECOD reference database

### Quality Metrics

1. **Evidence Quality**
   ```python
   def calculate_evidence_quality(evidence: Evidence) -> str:
       """Classify evidence quality"""
       if evidence.confidence > 0.8 and evidence.evalue < 1e-10:
           return "high"
       elif evidence.confidence > 0.5 and evidence.evalue < 1e-5:
           return "medium"
       else:
           return "low"
   ```

2. **Domain Quality**
   ```python
   def assess_domain_quality(domain: Domain) -> Dict[str, Any]:
       """Comprehensive domain quality assessment"""
       return {
           'confidence': domain.confidence,
           'reference_coverage': domain.primary_evidence.get_reference_coverage(),
           'evidence_count': len(domain.evidence_items),
           'size_reasonable': 25 <= domain.length <= 1000,
           'quality_category': 'good' | 'questionable' | 'poor'
       }
   ```

3. **Overall Quality Distribution**
   - Good quality: >= 60%
   - Questionable: <= 30%
   - Poor quality: <= 10%

## Edge Cases & Stress Tests

### Edge Case Coverage

1. **Empty Evidence**
   - Protein with no BLAST/HHsearch hits
   - Expected: Graceful handling, 0 domains

2. **Single Evidence Item**
   - Only one weak hit
   - Expected: Conservative assignment or skip

3. **Highly Overlapping Evidence**
   - 10+ evidence items for same region
   - Expected: Proper conflict resolution

4. **Very Long Proteins**
   - Sequence length > 2000 residues
   - Expected: No performance degradation

5. **Very Short Proteins**
   - Sequence length < 50 residues
   - Expected: Proper handling, reasonable domains

6. **Discontinuous Domains**
   - Multi-segment domains
   - Expected: Correct segment merging

7. **Missing Reference Lengths**
   - Evidence without reference data
   - Expected: Skip or use alternative metrics

### Stress Tests

1. **Large Batch Processing**
   - 10,000 proteins
   - Monitor: Memory, speed, crash rate

2. **Concurrent Processing**
   - 50 parallel SLURM jobs
   - Monitor: Resource contention, failures

3. **Repeated Runs**
   - Same protein 100 times
   - Expected: Identical results (deterministic)

## Validation Workflow

### Automated Validation

```bash
# Run all tests
pytest tests/ --cov=src/pyecod_mini --cov-report=html

# Run only regression tests
pytest tests/ -m regression -v

# Run performance benchmarks
pytest tests/ -m performance --durations=10

# Generate validation report
python scripts/validate_algorithm.py --full-report
```

### Manual Validation

1. **Visual Inspection**
   - Review 10 random results
   - Compare with PyMOL structures
   - Check biological plausibility

2. **Expert Review**
   - Submit results to ECOD expert
   - Review flagged cases
   - Adjust algorithm if needed

3. **Literature Comparison**
   - Compare with published domain definitions
   - Verify against PDB annotations

## Validation Report Template

```markdown
# Algorithm Validation Report
Date: YYYY-MM-DD
Version: 2.0.0

## Test Summary
- Total Tests Run: XXX
- Tests Passed: XXX
- Tests Failed: XXX
- Success Rate: XX%

## Regression Tests (Critical)
- 8ovp_A: ✅ PASS (3 domains, 82% coverage)
- Test 2: ✅ PASS
- Test 3: ✅ PASS
- Test 4: ✅ PASS
- Test 5: ✅ PASS
- Test 6: ✅ PASS

## Production Tests (100 proteins)
- Success Rate: 96%
- Average Coverage: 78%
- Average Domains: 2.3
- Classification Accuracy: 83%

## Performance
- Average Time: 6.2s per protein
- Peak Memory: 320 MB
- Throughput: 580 proteins/hour

## Quality Distribution
- Good: 68%
- Questionable: 24%
- Poor: 8%

## Known Issues
1. Issue description
2. Issue description

## Recommendations
1. Recommendation
2. Recommendation

## Conclusion
[PASS/FAIL with justification]
```

## Success Criteria Summary

**Required for Production Release:**
- ✅ All 6 regression tests pass
- ✅ >= 95% success rate on 100-protein test set
- ✅ >= 75% average coverage
- ✅ >= 80% ECOD classification accuracy
- ✅ Performance benchmarks met
- ✅ No critical bugs or crashes
- ✅ Code coverage >= 80%
- ✅ All type checks pass (mypy)

**Optional (Nice to Have):**
- >= 90% exact domain count match
- >= 90% boundary accuracy (±5 residues)
- >= 85% classification accuracy
- Performance < 5s per protein

## Timeline

1. **Unit Tests**: 1 day
2. **Regression Tests**: 1 day (includes test case selection)
3. **Integration Tests**: 0.5 day
4. **Production Tests**: 1 day
5. **Report Generation**: 0.5 day

**Total**: 4 days for comprehensive validation

## Next Steps After Validation

1. If **PASS**: Proceed to production framework implementation
2. If **FAIL**: Debug issues, fix algorithm, re-validate

---

**Note**: This validation plan ensures the algorithm is production-ready before deploying at scale.
