# GitHub Actions Workflows

This directory contains CI/CD workflows for pyecod_mini.

## Workflows

### 1. `tests.yml` - Main Test Suite

**Triggers:** Push to main/develop, Pull Requests

**Jobs:**
- **test**: Runs unit tests on Python 3.9, 3.10, 3.11
  - Executes all unit tests with coverage
  - Generates coverage reports
  - Uploads to Codecov
  - Enforces 80% coverage minimum

- **integration-tests**: Runs integration tests (Python 3.11 only)
  - Runs after unit tests pass
  - Allows failure (may require external resources)

- **lint**: Code quality checks
  - Ruff linting
  - Black formatting
  - isort import sorting
  - mypy type checking

- **test-summary**: Aggregates results

### 2. `pr-checks.yml` - Pull Request Validation

**Triggers:** Pull Request events

**Jobs:**
- **quick-tests**: Fast validation for PRs
  - Runs unit tests (excluding slow tests)
  - Uses parallel execution (`pytest-xdist`)
  - Quick lint check
  - Validates test file naming
  - Warns about large PRs (>50 files)

- **new-test-check**: Ensures tests accompany code changes
  - Warns if source files changed without test updates

### 3. `coverage.yml` - Comprehensive Coverage

**Triggers:** Push to main, Weekly schedule

**Jobs:**
- **coverage**: Full coverage analysis
  - Runs complete test suite
  - Generates detailed coverage reports
  - Creates coverage badge
  - Uploads to Codecov
  - Publishes HTML reports as artifacts

## Running Tests Locally

### Quick Unit Tests
```bash
pytest tests/ -m "unit" -v
```

### Fast Tests Only (for rapid development)
```bash
pytest tests/ -m "unit and not slow" -n auto
```

### With Coverage
```bash
pytest tests/ -m "unit" --cov=src/pyecod_mini --cov-report=term-missing
```

### Integration Tests
```bash
pytest tests/ -m "integration" -v
```

### All Tests
```bash
pytest tests/ -v
```

## Test Markers

- `unit`: Fast, isolated unit tests
- `integration`: Tests requiring real data/external resources
- `regression`: Critical validation cases
- `slow`: Long-running tests
- `production`: Requires SLURM/database (not run in CI)
- `performance`: Performance benchmarks

## Coverage Requirements

- **Minimum:** 80% (enforced in CI)
- **Target:** 90%+
- **Measured:** `src/pyecod_mini/` only
- **Excluded:** Tests, `__main__.py`, type checking blocks

## Artifacts

### Coverage Reports
- **Location:** Workflow run → Artifacts
- **Retention:** 7 days (quick tests), 30 days (full coverage)
- **Contents:** HTML coverage report (open `index.html`)

### Codecov
- **Dashboard:** https://codecov.io/gh/YOUR_ORG/pyecod_mini
- **Badge:** Add to README.md:
  ```markdown
  [![codecov](https://codecov.io/gh/YOUR_ORG/pyecod_mini/branch/main/graph/badge.svg)](https://codecov.io/gh/YOUR_ORG/pyecod_mini)
  ```

## Setting Up Codecov (Optional)

1. Go to https://codecov.io/
2. Sign in with GitHub
3. Add `pyecod_mini` repository
4. Get your upload token
5. Add as GitHub secret: `CODECOV_TOKEN`
6. Update workflows to use token:
   ```yaml
   - uses: codecov/codecov-action@v4
     with:
       token: ${{ secrets.CODECOV_TOKEN }}
   ```

## Troubleshooting

### Tests Failing in CI but Passing Locally

**Common causes:**
- Missing dependencies in CI
- Different Python versions
- Relying on local files not in repo
- Environment-specific paths

**Solution:**
- Check workflow logs
- Ensure test data is committed or generated in CI
- Use relative paths
- Mock external dependencies

### Coverage Below Threshold

**Check:**
```bash
coverage report --show-missing
```

**Fix:**
- Add tests for uncovered lines
- Remove dead code
- Add `# pragma: no cover` for truly untestable code

### Slow CI Runs

**Optimizations:**
- Use `pytest-xdist` for parallelization: `-n auto`
- Skip slow tests in quick checks: `-m "not slow"`
- Cache pip dependencies (already configured)
- Use `--maxfail=3` to stop early

## Adding New Tests

1. Create test file: `tests/test_<module>.py`
2. Use appropriate markers:
   ```python
   @pytest.mark.unit
   def test_something():
       ...
   ```
3. Run locally first:
   ```bash
   pytest tests/test_<module>.py -v
   ```
4. Check coverage:
   ```bash
   pytest tests/test_<module>.py --cov=src/pyecod_mini/<module>.py --cov-report=term-missing
   ```
5. Commit and push

## Workflow Status Badges

Add to your README.md:

```markdown
[![Tests](https://github.com/YOUR_ORG/pyecod_mini/workflows/Tests/badge.svg)](https://github.com/YOUR_ORG/pyecod_mini/actions/workflows/tests.yml)
[![Coverage](https://github.com/YOUR_ORG/pyecod_mini/workflows/Coverage%20Report/badge.svg)](https://github.com/YOUR_ORG/pyecod_mini/actions/workflows/coverage.yml)
```
