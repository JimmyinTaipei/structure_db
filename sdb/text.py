"""Text → SVG outline paths, so every viewer (Typora, PDF, resvg PNG) renders identically.

Markup (used everywhere in YAML):
    CH_2OH        subscript single char        _{12}  subscript group
    NAD^+         superscript single char      ^{2-}  superscript group
    OPO_3^{2-}    a subscript directly followed by a superscript is stacked (² over ₃)
    O{P}          {P} = the phosphoryl symbol Ⓟ (orange disc, dark ring, serif P)
A '-' inside a superscript is drawn as a true minus sign.
"""
from functools import lru_cache

from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
from fontTools.ttLib import TTFont

from . import style

SCRIPT_SCALE = 0.68
PHOS = "\u24C5"            # Ⓟ — stands for a phosphoryl group, drawn by hand (not in the fonts)
PHOS_FILL = "#F3C388"


class Font:
    def __init__(self, path):
        self.tt = TTFont(str(path))
        self.upm = self.tt["head"].unitsPerEm
        self.cmap = self.tt.getBestCmap()
        self.glyphs = self.tt.getGlyphSet()
        self.hmtx = self.tt["hmtx"]
        os2 = self.tt["OS/2"]
        self.cap = getattr(os2, "sCapHeight", 0) or 0.7 * self.upm
        self.xh = getattr(os2, "sxHeight", 0) or 0.5 * self.upm
        self.desc = abs(self.tt["hhea"].descent)

    def gname(self, ch):
        name = self.cmap.get(ord(ch))
        if name is None:
            raise ValueError(f"glyph {ch!r} (U+{ord(ch):04X}) missing from font")
        return name

    def advance(self, ch):
        return self.hmtx[self.gname(ch)][0]

    def outline(self, ch, x, y, size):
        """Absolute SVG path data for one glyph with its baseline-left at (x, y)."""
        k = size / self.upm
        pen = SVGPathPen(self.glyphs, ntos=lambda v: f"{v:.2f}".rstrip("0").rstrip("."))
        self.glyphs[self.gname(ch)].draw(TransformPen(pen, (k, 0, 0, -k, x, y)))
        return pen.getCommands()


@lru_cache(maxsize=None)
def font(kind):
    return Font(style.SERIF_FONT if kind == "serif" else style.SANS_FONT)


def parse(markup):
    """Return a list of segments: (level, text) with level in {'n','sub','sup'}."""
    segs, i = [], 0
    markup = markup.replace("{P}", PHOS)
    while i < len(markup):
        c = markup[i]
        if c in "_^" and i + 1 < len(markup):
            level = "sub" if c == "_" else "sup"
            if markup[i + 1] == "{":
                j = markup.index("}", i)
                body, i = markup[i + 2:j], j + 1
            else:
                body, i = markup[i + 1], i + 2
            if level == "sup":
                body = body.replace("-", "−")
            segs.append((level, body))
        else:
            if segs and segs[-1][0] == "n":
                segs[-1] = ("n", segs[-1][1] + c)
            else:
                segs.append(("n", c))
            i += 1
    return segs


def plain(markup):
    return "".join(t for _, t in parse(markup))


def anchor_index(markup):
    """Index (in plain text) of the atom that chain bonds attach to: the first C, otherwise the
    first heavy-atom letter ('NH' → N, '^+NH_3' → N, 'SH' → S, 'OH' → O)."""
    p = plain(markup)
    if "C" in p:
        return p.index("C")
    for i, ch in enumerate(p):
        if ch.isascii() and ch.isupper() and ch != "H":
            return i
    return 0


class TextLine:
    """One line of rich text. Coordinates are relative to its baseline-left origin."""

    def __init__(self, markup, kind="serif", size=None):
        self.markup = markup
        self.kind = kind
        self.size = size or (style.S if kind == "serif" else style.NAME_SIZE)
        f = font(kind)
        k = self.size / f.upm
        ss = self.size * SCRIPT_SCALE
        sub_dy = self.size * 0.30
        sup_dy = -self.size * 0.42
        self.chars = []          # (char, x, dy, size)
        x = 0.0
        segs = parse(markup)
        i = 0
        while i < len(segs):
            level, body = segs[i]
            if level == "n":
                for ch in body:
                    if ch == PHOS:
                        x += self.size * 0.06
                        self.chars.append((ch, x, 0.0, self.size))
                        x += self._phos_d() + self.size * 0.06
                        continue
                    self.chars.append((ch, x, 0.0, self.size))
                    x += f.advance(ch) * k
                i += 1
                continue
            # script run; stack sub+sup if adjacent
            group = [segs[i]]
            if i + 1 < len(segs) and segs[i + 1][0] not in ("n", level):
                group.append(segs[i + 1])
            widths = []
            for lv, b in group:
                xx = x + self.size * 0.02
                for ch in b:
                    self.chars.append((ch, xx, sub_dy if lv == "sub" else sup_dy, ss))
                    xx += f.advance(ch) * ss / f.upm
                widths.append(xx - x)
            x += max(widths) + self.size * 0.03
            i += len(group)
        self.width = x
        cap = f.cap * k
        has_sup = any(dy < 0 for _, _, dy, _ in self.chars)
        has_sub = any(dy > 0 for _, _, dy, _ in self.chars)
        self.cap = cap
        self.top = -(cap + (self.size * 0.2 if has_sup else 0))
        self.bottom = self.size * (0.22 + (0.18 if has_sub else 0))

    def _phos_d(self):
        """Diameter of the Ⓟ disc: a bit taller than a capital letter."""
        f = font(self.kind)
        return f.cap * self.size / f.upm * 1.3

    def char_box(self, idx):
        """(x0, x1) of the idx-th rendered character."""
        f = font(self.kind)
        ch, x, dy, sz = self.chars[idx]
        if ch == PHOS:
            return x, x + self._phos_d()
        return x, x + f.advance(ch) * sz / f.upm

    def span(self, sub, start=0):
        """(start, end) character indices of the markup fragment `sub`, or None."""
        p, q = plain(self.markup), plain(sub)
        i = p.find(q, start)
        return None if i < 0 or not q else (i, i + len(q))

    def _phos(self, x, y, color):
        f = font("serif")
        d = self._phos_d()
        cap = f.cap * self.size / f.upm
        cx, cy = x + d / 2, y - cap / 2
        psize = d * 0.55 / (f.cap / f.upm)           # P's cap height = 55 % of the disc
        pw = f.advance("P") * psize / f.upm
        pcap = f.cap * psize / f.upm
        ring = (f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{d / 2 - d * 0.035:.2f}" fill="{PHOS_FILL}" '
                f'stroke="{style.color(color)}" stroke-width="{d * 0.07:.2f}"/>')
        letter = f'<path d="{f.outline("P", cx - pw / 2 + psize * 0.02, cy + pcap / 2, psize)}" fill="{style.color(color)}"/>'
        return ring + letter

    def anchor_x(self, idx):
        x0, x1 = self.char_box(idx)
        return (x0 + x1) / 2

    def svg(self, x, y, fill=None, extra="", spans=None):
        """spans: [(start, end, colour)] recolours those characters (later spans win)."""
        f = font(self.kind)
        cols = [fill] * len(self.chars)
        for a, b, c in spans or []:
            for i in range(a, b):
                cols[i] = c
        out, run, run_col = [], [], None

        def flush():
            if run:
                out.append(f'<path d="{"".join(run)}" fill="{style.color(run_col)}"{extra}/>')

        for (ch, cx, dy, sz), c in zip(self.chars, cols):
            if ch == PHOS:
                out.append(self._phos(x + cx, y + dy, c))
                continue
            if c != run_col and run:
                flush()
                run = []
            run_col = c
            run.append(f.outline(ch, x + cx, y + dy, sz))
        flush()
        return "".join(out)


class TextBlock:
    """Several lines, centred on x=0, first baseline at y=0."""

    def __init__(self, lines, kind="serif", size=None, leading=1.18, align="center"):
        if isinstance(lines, str):
            lines = lines.split("\n")
        self.lines = [TextLine(l, kind, size) for l in lines]
        self.size = self.lines[0].size if self.lines else 0
        self.lh = self.size * leading
        self.align = align
        self.width = max((l.width for l in self.lines), default=0)
        self.top = self.lines[0].top if self.lines else 0
        self.bottom = (len(self.lines) - 1) * self.lh + (self.lines[-1].bottom if self.lines else 0)

    def svg(self, cx, y, fill=None):
        out = []
        for i, l in enumerate(self.lines):
            if self.align == "center":
                x = cx - l.width / 2
            elif self.align == "left":
                x = cx
            else:
                x = cx - l.width
            out.append(l.svg(x, y + i * self.lh, fill))
        return "".join(out)
