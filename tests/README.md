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

Run tests from the repository root using one of the three test tiers:

```bash
uv sync
source .venv/bin/activate
```

```bash
test-light   # Fast unit tests
test-medium  # Every test except the advanced example
test-heavy   # Every test, including the advanced example
```

Additional pytest arguments are passed through, for example:

```bash
test-light -x
test-medium -k serialization
```
