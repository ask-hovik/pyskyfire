# pyskyfire/viz/report.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional, Any, Mapping
import html
import plotly.io as pio
import plotly.graph_objects as go
import numpy as _np
import base64, mimetypes, os

# -------- Helper functions ----------------------------------------------------

@dataclass
class _Block:
    kind: str  # "text" | "figure" | "raw" | "image" | "tikzjax"
    payload: str

def _file_to_data_uri(path: str) -> str:
    mime, _ = mimetypes.guess_type(path)
    if not mime:
        # default to png if unknown
        mime = "image/png"
    with open(path, "rb") as f:
        b64 = base64.b64encode(f.read()).decode("ascii")
    return f"data:{mime};base64,{b64}"

def _is_fluid(obj: Any) -> bool:
    # Duck-typing to avoid circular imports
    return hasattr(obj, "propellants") and hasattr(obj, "fractions")

def _is_material(obj: Any) -> bool:
    """
    Duck-typing to detect your Material objects without importing from solids.py.
    """
    return (
        hasattr(obj, "name") and isinstance(getattr(obj, "name"), str)
        # Optionally tighten the check to avoid false positives:
        # and any(hasattr(obj, m) for m in ("k", "E", "alpha", "nu", "rho"))
    )

def _format_material(mat: Any) -> str:
    """
    Render just the material's name, e.g., 'Inconel 718'.
    Keep it minimal per your request.
    """
    try:
        return getattr(mat, "name", str(mat))
    except Exception:
        return str(mat)

def _format_number(x: float, precision: int = 6) -> str:
    # Nice general-purpose number formatter (keeps integers neat, handles small/large)
    try:
        if x is None:
            return "—"
        # Show integers without decimal point
        if isinstance(x, (int,)) or (isinstance(x, float) and x.is_integer()):
            return f"{int(x)}"
        # Use g-format with a cap on significant digits
        return f"{x:.{precision}g}"
    except Exception:
        return str(x)

def _format_iterable(it, precision: int = 6) -> str:
    # Format list/tuple/np-array of scalars reasonably
    try:
        if isinstance(it, _np.ndarray):
            it = it.tolist()
        parts = []
        for v in it:
            if isinstance(v, (int, float)):
                parts.append(_format_number(v, precision))
            else:
                parts.append(str(v))
        return "[" + ", ".join(parts) + "]"
    except Exception:
        return str(it)

def _format_fluid(fluid: Any, precision_pct: int = 2) -> str:
    # Render: "C2H5OH (89%), H2O (10%), SiO2 (1%)"
    try:
        props = list(getattr(fluid, "propellants", []))
        fracs = list(getattr(fluid, "fractions", []))
        pairs = []
        for p, f in zip(props, fracs):
            pct = round(float(f) * 100.0, precision_pct)
            # drop trailing .0 if possible
            pct_str = f"{pct:.{precision_pct}f}".rstrip("0").rstrip(".")
            pairs.append(f"{p} ({pct_str}%)")
        return ", ".join(pairs) if pairs else "(empty fluid)"
    except Exception:
        return str(fluid)

def format_value(value: Any, precision: int = 6) -> str:
    """Human-friendly string for table cells."""
    # Fluid
    if _is_fluid(value):
        return _format_fluid(value)
    # Material
    if _is_material(value):
        return _format_material(value)
    # None
    if value is None:
        return "—"
    # Numpy scalars
    if isinstance(value, (_np.generic,)):
        value = value.item()
    # Numbers
    if isinstance(value, (int, float)):
        return _format_number(value, precision)
    # Strings
    if isinstance(value, str):
        return value
    # Iterables (lists/tuples/arrays)
    if isinstance(value, (list, tuple, _np.ndarray)):
        return _format_iterable(value, precision)
    # Mappings (dict) — one-line JSON-ish for compactness
    if isinstance(value, Mapping):
        try:
            import json as _json
            return _json.dumps(value, ensure_ascii=False)
        except Exception:
            return str(value)
    # Fallback
    return str(value)

def dict_to_table_html(
    data: Mapping[str, Any],
    *,
    col_key: str = "Key",
    col_val: str = "Value",
    caption: str | None = None,
    precision: int = 6,
) -> str:
    """Return HTML for a simple 2-column table of a dict."""
    # Escape keys and values (values are preformatted strings)
    rows = []
    for k, v in data.items():
        k_html = html.escape(str(k))
        v_html = html.escape(format_value(v, precision))
        rows.append(f"<tr><th>{k_html}</th><td>{v_html}</td></tr>")

    cap_html = f"<caption>{html.escape(caption)}</caption>" if caption else ""
    table = [
        f"<table class='psf-table'>",
        cap_html,
        f"<thead><tr><th>{html.escape(col_key)}</th><th>{html.escape(col_val)}</th></tr></thead>",
        "<tbody>",
        *rows,
        "</tbody></table>",
    ]
    return "\n".join(table)


# --- Content blocks ----------------------------------------------------------

@dataclass
class _Block:
    kind: str  # "text" | "figure" | "raw"
    payload: str  # already HTML for "raw" and "figure"; plain text for "text"

# --- Tab container -----------------------------------------------------------

@dataclass
class Tab:
    title: str
    blocks: List[_Block] = field(default_factory=list)

    def add_text(self, text: str) -> "Tab":
        # simple safe text (escape + preserve newlines)
        safe = html.escape(text).replace("\n", "<br>")
        self.blocks.append(_Block("text", safe))
        return self

    def add_figure(self, fig: go.Figure, caption: Optional[str] = None) -> "Tab":
        div = pio.to_html(
                        fig,
                        include_plotlyjs=False,  # you’re already loading plotly in <head>
                        include_mathjax='cdn',   # add MathJax script tag alongside the div
                        full_html=False
                    )
        if caption:
            cap = f"<div class='psf-caption'>{html.escape(caption)}</div>"
        else:
            cap = ""
        self.blocks.append(_Block("figure", f"<div class='psf-figure'>{div}{cap}</div>"))
        return self
    
    def add_table(self, data: Mapping[str, Any], *, caption: Optional[str] = None,
                  key_title: str = "Key", value_title: str = "Value",
                  precision: int = 6) -> "Tab":
        """Append a key/value HTML table made from a dict."""
        table_html = dict_to_table_html(
            data, col_key=key_title, col_val=value_title, caption=caption, precision=precision
        )
        return self.add_raw_html(table_html)


    def add_raw_html(self, raw_html: str) -> "Tab":
        self.blocks.append(_Block("raw", raw_html))
        return self
    
    def add_image(self, path: str, *, alt: str = "", caption: str | None = None,
                  style: str = "max-width:100%;height:auto;") -> "Tab":
        data_uri = _file_to_data_uri(path)
        cap = f"<div class='psf-caption'>{html.escape(caption)}</div>" if caption else ""
        img_html = f"<img src='{data_uri}' alt='{html.escape(alt)}' style='{style}'/>"
        self.blocks.append(_Block("image", f"<div class='psf-figure'>{img_html}{cap}</div>"))
        return self

    def add_svg(self, svg: str, *, caption: str | None = None) -> "Tab":
        cap = f"<div class='psf-caption'>{html.escape(caption)}</div>" if caption else ""
        self.blocks.append(_Block("raw", f"<div class='psf-figure'>{svg}{cap}</div>"))
        return self

# --- Report ------------------------------------------------------------------

class Report:
    def __init__(self, title: str = "pyskyfire Report"):
        self.title = title
        self.tabs: List[Tab] = []

    def add_tab(self, title: str) -> Tab:
        tab = Tab(title=title)
        self.tabs.append(tab)
        return tab

    def save_html(self, path: str):
        if not self.tabs:
            raise RuntimeError("Report has no tabs. Add at least one via add_tab().")

        # --- Head / styles / JS ---
        head = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>{html.escape(self.title)}</title>
<link rel="preconnect" href="https://cdn.plot.ly">
<style>
:root {{
  --sidebar-w: 260px;
  --bg: #ffffff;
  --fg: #111;
  --muted: #666;
  --border: #ddd;
  --accent: #0d6efd;
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  font-family: system-ui, -apple-system, Segoe UI, Roboto, Ubuntu, sans-serif;
  background: var(--bg);
  color: var(--fg);
}}
.header {{
  padding: 16px 20px;
  border-bottom: 1px solid var(--border);
  font-size: 20px;
  font-weight: 600;
}}
.container {{
  display: grid;
  grid-template-columns: var(--sidebar-w) 1fr;
  min-height: calc(100vh - 58px);
}}
.sidebar {{
  border-right: 1px solid var(--border);
  padding: 12px 0;
  overflow: auto;
}}
.tab-btn {{
  width: 100%;
  text-align: left;
  background: transparent;
  border: 0;
  outline: none;
  padding: 10px 16px;
  font-size: 14px;
  color: var(--fg);
  cursor: pointer;
  border-left: 3px solid transparent;
}}
.tab-btn:hover {{
  background: #f6f8fa;
}}
.tab-btn.active {{
  background: #eef4ff;
  border-left-color: var(--accent);
  color: #0a58ca;
  font-weight: 600;
}}
.content {{
  padding: 18px 22px 40px;
}}
.tab-panel {{ display: none; }}
.tab-panel.active {{ display: block; }}
/* Blocks */
.psf-block {{ margin: 16px 0; }}
.psf-text {{ line-height: 1.5; color: var(--fg); }}
.psf-figure {{ margin: 12px 0; }}
.psf-caption {{ font-size: 12px; color: var(--muted); margin-top: 4px; }}
/* Make Plotly plots use available width */
.psf-figure .plotly-graph-div {{ width: 100% !important; }}
/* --- Key/Value tables --- */
.psf-table {{
  width: 100%;
  border-collapse: collapse;
  margin: 12px 0;
  font-size: 14px;
}}
.psf-table th, .psf-table td {{
  border: 1px solid var(--border);
  padding: 8px 10px;
  vertical-align: top;
}}
.psf-table th {{
  background: #f7f9fb;
  text-align: left;
  width: 34%;
}}
.psf-table caption {{
  text-align: left;
  font-weight: 600;
  margin: 0 0 6px 0;
}}
</style>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<script>
function psfActivateTab(idx) {{
  const btns = document.querySelectorAll('.tab-btn');
  const panels = document.querySelectorAll('.tab-panel');
  btns.forEach((b, i) => {{
    if (i === idx) b.classList.add('active'); else b.classList.remove('active');
  }});
  panels.forEach((p, i) => {{
    if (i === idx) p.classList.add('active'); else p.classList.remove('active');
  }});
  // Persist selection
  try {{ localStorage.setItem('psf_active_tab', String(idx)); }} catch(e) {{}}
}}
window.addEventListener('DOMContentLoaded', () => {{
  // Restore last active tab
  let idx = 0;
  try {{
    const s = localStorage.getItem('psf_active_tab');
    if (s !== null) idx = Math.max(0, parseInt(s, 10) || 0);
  }} catch(e) {{}}
  psfActivateTab(idx);
}});
</script>
</head>
<body>
<div class="header">{html.escape(self.title)}</div>
<div class="container">
  <nav class="sidebar">
"""

        # --- Sidebar (tab buttons) ---
        sidebar_btns = []
        for idx, tab in enumerate(self.tabs):
            sidebar_btns.append(
                f'<button class="tab-btn" onclick="psfActivateTab({idx})">{html.escape(tab.title)}</button>'
            )
        sidebar_html = "\n".join(sidebar_btns) + "\n  </nav>\n"

        # --- Content panels ---
        content_start = '  <main class="content">\n'
        panels = []
        for idx, tab in enumerate(self.tabs):
            blocks_html = []
            for blk in tab.blocks:
                if blk.kind == "text":
                    blocks_html.append(f'<div class="psf-block psf-text">{blk.payload}</div>')
                elif blk.kind == "figure":
                    blocks_html.append(f'<div class="psf-block psf-figure">{blk.payload}</div>')
                elif blk.kind == "raw":
                    blocks_html.append(f'<div class="psf-block">{blk.payload}</div>')
                elif blk.kind == "image":
                    blocks_html.append(f'<div class="psf-block psf-figure">{blk.payload}</div>')
            panel_html = f'<section class="tab-panel" id="psf-tab-{idx}">\n' + "\n".join(blocks_html) + "\n</section>"
            panels.append(panel_html)
        content_html = content_start + "\n".join(panels) + "\n  </main>\n</div>\n</body>\n</html>"

        # --- Write file ---
        with open(path, "w", encoding="utf-8") as f:
            f.write(head + sidebar_html + content_html)
