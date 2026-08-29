"""Markdown prose support in generated reports."""

import ast
from pathlib import Path

import pytest

from pyskyfire.viz.report import Report, Tab, markdown_to_html


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def _render(tmp_path, markdown):
    report = Report("Markdown report")
    report.add_tab("Guide").add_markdown(markdown)
    output_path = tmp_path / "report.html"
    report.save_html(output_path)
    return output_path.read_text()


def test_add_markdown_renders_formatted_content_and_dedents(tmp_path) -> None:
    document = _render(
        tmp_path,
        """
            ## Reproduction

            Use **care** and read the [guide](https://example.com).

            1. Run `regen_sim.py`.
            2. Build the report.

            | File | Purpose |
            | --- | --- |
            | `result.pkl` | Saved result |
        """,
    )

    assert '<div class="psf-block psf-markdown">' in document
    assert "<h2>Reproduction</h2>" in document
    assert "Use <strong>care</strong>" in document
    assert '<a href="https://example.com">guide</a>' in document
    assert "<ol>" in document
    assert '<table class="psf-table psf-markdown-table">' in document
    assert "<code>result.pkl</code>" in document


def test_add_markdown_escapes_raw_html(tmp_path) -> None:
    document = _render(tmp_path, "Before <script>alert('x')</script> after")

    assert "<script>alert('x')</script>" not in document
    assert "&lt;script&gt;alert('x')&lt;/script&gt;" in document


def test_add_markdown_rejects_non_string_input() -> None:
    tab = Report("Markdown report").add_tab("Guide")

    with pytest.raises(TypeError, match="Markdown source must be a string"):
        tab.add_markdown(123)


def test_plain_text_api_has_been_replaced_by_markdown() -> None:
    assert not hasattr(Tab, "add_text")


def test_rl10_reproduction_guide_is_valid_markdown() -> None:
    source_path = REPOSITORY_ROOT / "validation" / "RL10" / "post_process.py"
    source = source_path.read_text(encoding="utf-8")
    module = ast.parse(source)
    guide_assignment = next(
        node
        for node in module.body
        if isinstance(node, ast.Assign)
        and any(
            isinstance(target, ast.Name)
            and target.id == "REPORT_REPRODUCTION_GUIDE"
            for target in node.targets
        )
    )
    guide = ast.literal_eval(guide_assignment.value)
    rendered = markdown_to_html(guide)

    assert "<h2>Files and reproduction</h2>" in rendered
    assert "<ol>" in rendered
    assert rendered.count('<table class="psf-table psf-markdown-table">') == 4
    assert "<code>regen_sim.py</code>" in rendered
    assert ".add_raw_html(" not in source
