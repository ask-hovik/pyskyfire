"""Report the composition of the Python source tree by physical lines."""

from __future__ import annotations

import argparse
import ast
import io
import tokenize
import warnings
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class LineCounts:
    """Mutually exclusive physical-line counts for one Python file."""

    path: Path
    total: int
    blank: int
    code: int
    docstrings: int
    comments: int

    @property
    def nonempty(self) -> int:
        return self.code + self.docstrings + self.comments


def _docstring_lines(tree: ast.AST) -> set[int]:
    """Return line numbers occupied by module, class, and function docstrings."""
    lines: set[int] = set()
    owners = (
        node
        for node in ast.walk(tree)
        if isinstance(
            node,
            (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef),
        )
    )

    for owner in owners:
        body = owner.body
        if not body or not isinstance(body[0], ast.Expr):
            continue
        value = body[0].value
        if isinstance(value, ast.Constant) and isinstance(value.value, str):
            lines.update(range(body[0].lineno, body[0].end_lineno + 1))

    return lines


def _comment_only_lines(source: str) -> set[int]:
    """Return lines containing a comment but no code before the comment."""
    lines: set[int] = set()
    source_lines = source.splitlines()

    try:
        tokens = tokenize.generate_tokens(io.StringIO(source).readline)
        for token in tokens:
            if token.type != tokenize.COMMENT:
                continue
            line_number, column = token.start
            if not source_lines[line_number - 1][:column].strip():
                lines.add(line_number)
    except tokenize.TokenError:
        # ast.parse reports the more useful syntax error when the file is analyzed.
        pass

    return lines


def analyze_file(path: Path) -> LineCounts:
    """Classify every physical line in *path*."""
    source = path.read_text(encoding="utf-8", errors="ignore")
    source_lines = source.splitlines()
    nonempty_lines = {
        number for number, line in enumerate(source_lines, start=1) if line.strip()
    }

    # Some embedded HTML/JavaScript strings contain escapes that Python warns
    # about while parsing. They do not affect the AST needed for this report.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SyntaxWarning)
        tree = ast.parse(source, filename=str(path))

    docstring_lines = nonempty_lines & _docstring_lines(tree)
    comment_lines = (
        nonempty_lines - docstring_lines
    ) & _comment_only_lines(source)
    code_lines = nonempty_lines - docstring_lines - comment_lines

    return LineCounts(
        path=path,
        total=len(source_lines),
        blank=len(source_lines) - len(nonempty_lines),
        code=len(code_lines),
        docstrings=len(docstring_lines),
        comments=len(comment_lines),
    )


def analyze_tree(source_root: Path) -> list[LineCounts]:
    """Analyze every Python file below *source_root*."""
    return [analyze_file(path) for path in sorted(source_root.rglob("*.py"))]


def _sum_counts(counts: Iterable[LineCounts]) -> dict[str, int]:
    items = list(counts)
    return {
        field: sum(getattr(item, field) for item in items)
        for field in ("total", "blank", "code", "docstrings", "comments", "nonempty")
    }


def _format_table(headers: list[str], rows: list[list[object]]) -> str:
    """Format a compact table, right-aligning numeric columns."""
    text_rows = [[str(cell) for cell in row] for row in rows]
    widths = [
        max([len(headers[index]), *(len(row[index]) for row in text_rows)])
        for index in range(len(headers))
    ]

    def format_row(row: list[str], header: bool = False) -> str:
        cells = []
        for index, cell in enumerate(row):
            if header or index == 0:
                cells.append(cell.ljust(widths[index]))
            else:
                cells.append(cell.rjust(widths[index]))
        return "  ".join(cells).rstrip()

    separator = ["-" * width for width in widths]
    return "\n".join(
        [format_row(headers, header=True), format_row(separator), *map(format_row, text_rows)]
    )


def _subpackage(path: Path, source_root: Path) -> str:
    relative = path.relative_to(source_root)
    return relative.parts[0] if len(relative.parts) > 1 else "(root)"


def render_report(source_root: Path, counts: list[LineCounts], top_n: int) -> str:
    """Render the complete line-count report."""
    totals = _sum_counts(counts)
    nonempty = totals["nonempty"]

    def percentage(value: int) -> str:
        return f"{value / nonempty:5.1%}" if nonempty else "  0.0%"

    sections = [
        f"Python line report for {source_root} ({len(counts)} files)",
        "",
        "Line types",
        _format_table(
            ["Type", "Lines", "% non-empty"],
            [
                ["Code", totals["code"], percentage(totals["code"])],
                ["Docstrings", totals["docstrings"], percentage(totals["docstrings"])],
                ["Comment-only", totals["comments"], percentage(totals["comments"])],
                ["Non-empty", totals["nonempty"], "100.0%" if nonempty else "0.0%"],
                ["Blank", totals["blank"], "-"],
                ["Physical total", totals["total"], "-"],
            ],
        ),
    ]

    by_subpackage: defaultdict[str, list[LineCounts]] = defaultdict(list)
    for item in counts:
        by_subpackage[_subpackage(item.path, source_root)].append(item)

    subpackage_rows = []
    for name, items in sorted(
        by_subpackage.items(),
        key=lambda entry: _sum_counts(entry[1])["nonempty"],
        reverse=True,
    ):
        subtotal = _sum_counts(items)
        subpackage_rows.append(
            [
                name,
                len(items),
                subtotal["code"],
                subtotal["docstrings"],
                subtotal["comments"],
                subtotal["nonempty"],
                percentage(subtotal["nonempty"]),
            ]
        )

    sections.extend(
        [
            "",
            "Subpackages",
            _format_table(
                ["Subpackage", "Files", "Code", "Docstrings", "Comments", "Non-empty", "Share"],
                subpackage_rows,
            ),
        ]
    )

    top_files = sorted(counts, key=lambda item: item.nonempty, reverse=True)[:top_n]
    file_rows = [
        [
            str(item.path.relative_to(source_root)),
            item.code,
            item.docstrings,
            item.comments,
            item.nonempty,
            item.total,
        ]
        for item in top_files
    ]
    sections.extend(
        [
            "",
            f"Top {len(top_files)} files by non-empty lines",
            _format_table(
                ["File", "Code", "Docstrings", "Comments", "Non-empty", "Physical"],
                file_rows,
            ),
        ]
    )

    return "\n".join(sections)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "source_root",
        nargs="?",
        type=Path,
        default=Path("src/pyskyfire"),
        help="Python package directory to analyze (default: src/pyskyfire)",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="number of largest files to display (default: 10)",
    )
    args = parser.parse_args()
    if args.top < 0:
        parser.error("--top must be zero or greater")
    if not args.source_root.is_dir():
        parser.error(f"source directory does not exist: {args.source_root}")
    return args


if __name__ == "__main__":
    arguments = parse_args()
    results = analyze_tree(arguments.source_root)
    print(render_report(arguments.source_root, results, arguments.top))
