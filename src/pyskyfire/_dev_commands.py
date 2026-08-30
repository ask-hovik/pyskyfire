"""Console entry points for repository development commands."""

import os
from pathlib import Path
import sys


def _repository_root(command: str) -> Path:
    candidates = (Path.cwd(), Path(__file__).resolve().parents[2])

    for candidate in candidates:
        for directory in (candidate, *candidate.parents):
            if (directory / "pyproject.toml").is_file() and (
                directory / "tools" / command
            ).is_file():
                return directory

    raise SystemExit(
        f"{command} is a repository development command; run it from a "
        "pyskyfire checkout"
    )


def _run(command: str) -> None:
    script = _repository_root(command) / "tools" / command
    os.execv(script, (command, *sys.argv[1:]))


def build_docs() -> None:
    _run("build-docs")


def build_dist() -> None:
    _run("build-dist")


def lint() -> None:
    _run("lint")


def test_light() -> None:
    _run("test-light")


def test_medium() -> None:
    _run("test-medium")


def test_heavy() -> None:
    _run("test-heavy")
