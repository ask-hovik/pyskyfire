"""Responsive clustering checks for consecutive report tables."""

from pyskyfire.viz.report import Report


def _render(tmp_path, build):
    report = Report("Table layout")
    tab = report.add_tab("Tables")
    build(tab)
    output_path = tmp_path / "report.html"
    report.save_html(output_path)
    return output_path.read_text()


def test_consecutive_tables_share_a_responsive_cluster(tmp_path) -> None:
    document = _render(
        tmp_path,
        lambda tab: (
            tab.add_table({"a": 1}, caption="First")
            .add_table({"b": 2}, caption="Second")
        ),
    )

    assert document.count('class="psf-block psf-table-cluster"') == 1
    assert document.count('class="psf-table-item"') == 2
    assert "display: flex;" in document
    assert "flex-wrap: wrap;" in document
    assert "flex: 1 1 18rem;" in document
    assert "min-width: min(18rem, 100%);" in document
    assert "overflow-x: auto;" in document


def test_non_table_block_splits_table_clusters(tmp_path) -> None:
    document = _render(
        tmp_path,
        lambda tab: (
            tab.add_table({"a": 1})
            .add_table({"b": 2})
            .add_markdown("A boundary")
            .add_table({"c": 3})
            .add_table({"d": 4})
        ),
    )

    assert document.count('class="psf-block psf-table-cluster"') == 2
    first_cluster, boundary, second_cluster = (
        document.index("<td>1</td>"),
        document.index("A boundary"),
        document.index("<td>3</td>"),
    )
    assert first_cluster < boundary < second_cluster


def test_single_table_keeps_a_standalone_table_block(tmp_path) -> None:
    document = _render(tmp_path, lambda tab: tab.add_table({"a": 1}))

    assert 'class="psf-block"><table class=\'psf-table\'>' in document
    assert 'class="psf-block psf-table-cluster"' not in document
