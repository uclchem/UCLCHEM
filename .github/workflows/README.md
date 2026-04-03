# GitHub Actions Workflows

Workflows are organized by tier and use `[phase:category]` namespacing in step names for clarity.

## Naming Convention

- `ci-*.yml` — Continuous Integration (run on PR/push, must pass)
- `cd-*.yml` — Continuous Deployment (run on releases, generate artifacts)

## Step Naming

Steps use `[phase:category]` format:
- `[setup:vcs]` — Version control operations
- `[setup:env]` — Environment setup
- `[setup:deps]` — Dependency installation
- `[build:package]` — Package building
- `[test:*]` — Test execution
- `[lint:*]` — Linting operations
- `[deploy:*]` — Release operations

## Architecture

**TIER 0: Gatekeeping** (Independent)
- `ci-lint.yml` — Auto-fix style, always passes
- `ci-check-fortran-metadata.yml` — Validate Fortran source alignment

**TIER 1: Foundation** (Sequential)
- `ci-test-suite.yml` → `install` job — Package installation, gates all tests

**TIER 2: Parallel Tests** (All depend on install)
- `ci-test-suite.yml` → Multiple jobs:
  - `test-functional` — Core tests (3-5m)
  - `test-doctests` — Docstring examples (2-3m)
  - `test-heating-cooling` — Feature tests (1-2m)
  - `test-networks` — Custom network validation (1-2m)
  - `test-notebooks` — Notebook execution (10-20m, main PRs only)

**TIER 3: Release** (Independent)
- `cd-notebooks-release.yml` — Build release artifacts (tags/releases only)

## Workflows

### `ci-lint.yml`

Auto-fixes code style using pre-commit hooks. Always passes (auto-commits fixes).

**Triggers:** Every push/PR to main/develop

### `ci-check-fortran-metadata.yml`

Validates `fortran_metadata.yaml` matches Fortran source code.

**Triggers:** Push/PR when Fortran files change

**Fails if:** Metadata out of sync with source

### `ci-test-suite.yml`

Consolidated workflow containing all test jobs with explicit dependencies via `needs:`.

#### `install` job
Installs package and runs install tests. Required for all other tests.

**Duration:** 1-2 minutes

#### `test-functional` job
Runs core test suite with parallel execution: `pytest -n auto tests/ -k "not install and not heating_cooling"`

**Duration:** 3-5 minutes

#### `test-doctests` job
Tests docstring examples: `pytest --doctest-modules ...`

**Duration:** 2-3 minutes

#### `test-heating-cooling` job
Tests heating/cooling features: `pytest tests/test_heating_cooling_install.py`

**Duration:** 1-2 minutes

#### `test-networks` job
Tests custom networks in `tests/networks/`:
1. Run makerates to generate rates
2. Load and verify network
3. Run makerates again to verify consistency

**Duration:** 1-2 minutes

**Auto-discovery:** Any directory in `tests/networks/` with `user_settings.yaml` is tested.

#### `test-notebooks` job
Generates and executes notebooks: `pytest --nbmake notebooks`

**Duration:** 10-20 minutes

**Trigger:** Main PRs only (expensive operation)

### `cd-notebooks-release.yml`

Builds executable notebooks as release artifacts.

**Triggers:** Version tags (v*, rc*) or releases

**Output:** Attaches ZIP archive to release

## Dependency Graph

```
ci-lint.yml ──→ ci-test-suite.yml
                ├─ install (required)
                └─ TIER 2 (parallel, all need install):
                   ├─ test-functional
                   ├─ test-doctests
                   ├─ test-heating-cooling
                   ├─ test-networks
                   └─ test-notebooks (main only)

ci-check-fortran-metadata.yml (independent, only when Fortran changes)

cd-notebooks-release.yml (independent, tags/releases only)
```

## Execution Timeline

```
T+0s   ci-lint starts (parallel with ci-check-fortran-metadata if Fortran changed)
T+30s  lint completes, ci-test-suite triggered
T+1m   install job starts
T+3m   install completes, TIER 2 tests start (parallel):
       - test-functional [3-5m]
       - test-doctests [2-3m]
       - test-heating [1-2m]
       - test-networks [1-2m]
       - test-notebooks [10-20m, main only]
T+10-22m All tests complete
```

**Total time:** 10-22 minutes (50% faster than sequential)

## Reusable Action

`.github/actions/setup-uclchem/action.yml` standardizes environment setup:
- Installs Python, gfortran, system dependencies
- Installs test requirements
- Optionally installs the package

**Usage:**
```yaml
- uses: ./.github/actions/setup-uclchem
  with:
    python-version: "3.12"
    install-package: "true"
```

## Local Testing

Run tests locally to match CI behavior:

```bash
# All custom network tests
pytest tests/test_custom_networks.py -v

# Specific network
pytest tests/test_custom_networks.py -k "small_chemistry" -v

# Core tests
pytest tests/ -k "not install and not heating_cooling"

# Pre-commit
pre-commit run --all-files
```

## Common Issues

**Lint fails with style violations:**
```bash
pre-commit run --all-files
```
Fixes are auto-committed in CI.

**Install test fails:**
Check `pyproject.toml` dependencies and build system (meson.build).

**Network test fails:**
Verify `tests/networks/<name>/user_settings.yaml` points to valid files. See `CUSTOM_NETWORK_TESTS.md`.

**Notebook timeout:**
Notebooks run only on main (main PRs). Optimize expensive computations or split into smaller examples.
