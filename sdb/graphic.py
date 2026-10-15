"""Tiny scene model: a Graphic is SVG markup in local coordinates plus a bounding box."""
from dataclasses import dataclass, field
from html import escape

from . import style


@dataclass
class Graphic:
    body: str
    x0: float
    y0: float
    x1: float
    y1: float
    anchors: dict = field(default_factory=dict)   # named points, e.g. "base" (structure bottom)
    parts: dict = field(default_factory=dict)     # named sub-Graphics, exported as separate stickers

    @property
    def w(self):
        return self.x1 - self.x0

    @property
    def h(self):
        return self.y1 - self.y0

    @property
    def cx(self):
        return (self.x0 + self.x1) / 2

    def moved(self, dx, dy, attrs=""):
        """Return a copy translated by (dx, dy), optionally wrapping it in a <g attrs>."""
        body = f'<g transform="translate({dx:.2f},{dy:.2f})"{attrs}>{self.body}</g>'
        anchors = {k: (x + dx, y + dy) for k, (x, y) in self.anchors.items()}
        parts = {k: g.moved(dx, dy) for k, g in self.parts.items()}
        return Graphic(body, self.x0 + dx, self.y0 + dy, self.x1 + dx, self.y1 + dy, anchors, parts)


class BBox:
    def __init__(self):
        self.x0 = self.y0 = float("inf")
        self.x1 = self.y1 = float("-inf")

    def add(self, x0, y0, x1, y1):
        self.x0, self.y0 = min(self.x0, x0), min(self.y0, y0)
        self.x1, self.y1 = max(self.x1, x1), max(self.y1, y1)

    def add_g(self, g):
        self.add(g.x0, g.y0, g.x1, g.y1)


def line(x1, y1, x2, y2, stroke=None, w=None, cap="butt"):
    return (f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
            f'stroke="{style.color(stroke)}" stroke-width="{w or style.LINE_W:.2f}" stroke-linecap="{cap}"/>')


def path(d, stroke=None, w=None, fill="none"):
    return (f'<path d="{d}" fill="{fill}" stroke="{style.color(stroke)}" '
            f'stroke-width="{w or style.LINE_W:.2f}" stroke-linecap="round" stroke-linejoin="round"/>')


def document(g, title="", pad=None, data=None):
    """Wrap a Graphic into a standalone SVG document with a transparent background."""
    pad = style.S * 0.6 if pad is None else pad
    x0, y0 = g.x0 - pad, g.y0 - pad
    w, h = g.w + 2 * pad, g.h + 2 * pad
    attrs = "".join(f' data-{k}="{escape(str(v))}"' for k, v in (data or {}).items())
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{x0:.2f} {y0:.2f} {w:.2f} {h:.2f}" '
            f'width="{w:.0f}" height="{h:.0f}"{attrs}>'
            f'<title>{escape(title)}</title>{g.body}</svg>\n')
