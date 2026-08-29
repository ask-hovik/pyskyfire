# Contributing to Pyskyfire

Thank you for considering a contribution to Pyskyfire. Contributions of all
sizes are welcome, including bug reports, documentation improvements, new test
cases, and code changes.

Pyskyfire is engineering software under active development. Changes to its
numerical models should therefore be clearly motivated, tested, and supported
by references or validation data where appropriate.

## Before you begin

- Search the [existing issues](https://github.com/ask-hovik/pyskyfire/issues)
  before opening a new one.
- For a substantial feature, API change, or numerical-model change, open an
  issue first so that the approach can be discussed before significant work is
  invested.
- Keep pull requests focused. Unrelated fixes are easier to review as separate
  contributions.

## Reporting bugs

Open an issue and include enough information to reproduce the problem:

- a concise description of the observed and expected behaviour;
- the smallest practical example that demonstrates the issue;
- the Pyskyfire and Python versions in use;
- the operating system and relevant dependency versions; and
- the complete traceback or error output, if applicable.

Please remove confidential or proprietary engineering data before sharing
inputs, logs, or results.

## Suggesting features

Feature requests are welcome. Describe the problem or engineering workflow the
feature would address, the proposed behaviour, and any alternatives you have
considered. For new physical models, include the source equations and relevant
literature or validation references where possible.

## Development setup

Pyskyfire requires Python 3.14 or newer and uses
[`uv`](https://docs.astral.sh/uv/) for dependency and environment management.

Fork the repository, clone your fork, and create a branch for the change:

```bash
git clone https://github.com/<your-username>/pyskyfire.git
cd pyskyfire
git switch -c <short-descriptive-branch-name>
```

Install the project and all development dependencies:

```bash
uv sync --group dev
```

If you only need the test dependencies, use:

```bash
uv sync --group test
```

## Making changes

- Follow the existing code structure and conventions.
- Use clear names and keep functions and classes focused.
- Add or update docstrings for public APIs.
- Add tests for bug fixes and new behaviour.
- Avoid changing public APIs unless the change has been discussed and the
  documentation and examples are updated accordingly.
- For numerical or physical-model changes, document assumptions, units,
  equations, and references. Add validation coverage where practical.
- Do not commit generated documentation, build artifacts, virtual environments,
  caches, or local simulation output.

## Running tests

Run the test suite from the repository root:

```bash
uv run --group test pytest -q tests
```

During development, you can run an individual test file or test by passing its
path to pytest:

```bash
uv run --group test pytest -q tests/test_package.py
uv run --group test pytest -q tests/test_package.py::<test-name>
```

Some example tests use optional visualisation dependencies and may require a
working display or an off-screen rendering setup. The GitHub Actions workflow
uses Xvfb for these tests on Linux.

## Documentation

Documentation is built with Sphinx. Install the documentation dependencies and
build the HTML documentation with:

```bash
rm -rf docs/_build/html
uv sync --group docs
uv run --group docs sphinx-build -b html docs docs/_build/html
```

The generated site is written to `docs/_build/html`. To force a completely
clean build, remove `docs/_build` before running Sphinx again.

Documentation changes should be clear, concise, and accompanied by runnable
examples where useful. Update relevant tutorials and examples when changing
user-facing behaviour.

## Commits and pull requests

Before opening a pull request:

1. Rebase or merge the latest `main` branch into your branch.
2. Run the relevant tests and documentation build locally.
3. Review the diff for accidental generated files or unrelated changes.
4. Write clear commit messages that describe the purpose of each change.

In the pull request description, explain:

- what changed and why;
- any user-visible or API changes;
- how the change was tested;
- relevant issues, references, or validation data; and
- any known limitations or follow-up work.

By submitting a contribution, you agree that it will be distributed under the
project's [MIT License](LICENSE).

## Other useful commands
```bash
.venv/bin/python validation/RL10_2/optimize_channel_heights.py --target heat_flux
``` 
