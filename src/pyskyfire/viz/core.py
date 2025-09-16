# pyskyfire/viz/core.py
from typing import Any, Dict, Iterable, Optional, Callable, List
import plotly.graph_objects as go
import plotly.io as pio

class _Adder:
    def __init__(self, fig: go.Figure): self._fig = fig
    def annotation(self, **kw): self._fig.add_annotation(**kw); return self
    def shape(self, **kw):      self._fig.add_shape(**kw);      return self
    def hline(self, **kw):      self._fig.add_hline(**kw);      return self
    def vline(self, **kw):      self._fig.add_vline(**kw);      return self

class _Node:
    """Nested, chainable config: plot.layout(...), plot.xaxis(...), plot.traces(...)."""
    def __init__(self, fig: go.Figure, applier: Callable[[go.Figure, Dict[str, Any]], None], prefix: Dict[str, Any] = None):
        self._fig = fig
        self._applier = applier
        self._prefix = prefix or {}
    def __call__(self, **kw):
        d = {**self._prefix, **kw}
        self._applier(self._fig, d)
        return self
    # allow deeper nesting: plot.layout.legend(...), plot.traces.line(...)
    def __getattr__(self, name: str):
        nested_prefix = {**self._prefix, name: {}}
        def applier(fig: go.Figure, d: Dict[str, Any]):
            # merge nested dicts one level
            base = {**self._prefix}
            # graft d into base at the last key
            k = list(nested_prefix.keys())[-1]
            base[k] = {**nested_prefix[k], **d}
            self._applier(fig, base)
        return _Node(self._fig, applier, nested_prefix[name])

class PlotBase:
    """Common helpers + Plotly passthrough."""
    def __init__(self, fig: Optional[go.Figure] = None):
        self.fig = fig or go.Figure()
        # Namespaces
        self.layout  = _Node(self.fig, lambda f, d: f.update_layout(**d))
        self.traces  = _Node(self.fig, lambda f, d: f.update_traces(**d))
        self.xaxis   = _Node(self.fig, lambda f, d: f.update_layout(xaxis=d))
        self.yaxis   = _Node(self.fig, lambda f, d: f.update_layout(yaxis=d))
        self.add     = _Adder(self.fig)
        pio.renderers.default = "browser"

    # Fluent utilities
    def template(self, name: str): self.fig.update_layout(template=name); return self
    def config(self, **cfg):       self._config = cfg; return self
    def show(self):                return self.fig.show(config=getattr(self, "_config", None))
    def save_html(self, path: str): self.fig.write_html(path, include_plotlyjs="cdn"); return self
    def save_png(self, path: str, scale: float = 2): self.fig.write_image(path, scale=scale); return self

    # Plotly passthrough: call any go.Figure method directly if needed
    def __getattr__(self, name): return getattr(self.fig, name)
