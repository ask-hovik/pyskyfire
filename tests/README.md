# Test suite structure

The test suite is grouped by the kind of contract each test protects:

- `unit/` contains fast tests of formulas, geometry helpers, data models, and
  rendering semantics in isolation.
- `subsystem/` exercises several cooperating package objects or an external
  scientific backend without running a complete validation case.
- `integration/` runs complete public workflows such as the example scripts.
- `validation/` protects scientific reference data and numerical validation
  cases against regression.
- `packaging/` tests the installed distribution, optional dependencies, and
  packaged runtime resources.
- `browser/` protects behavior implemented by generated browser-side HTML,
  CSS, and JavaScript.

Run every test from the repository root with:

```bash
uv run --group test pytest -q tests
```

During development, a category can be run independently, for example:

```bash
uv run --group test pytest -q tests/unit
uv run --group test pytest -q tests/validation
```
