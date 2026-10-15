"""Vertical-chain drawer: Lehninger-style condensed formulas AND Fischer projections.

A chain is a list of rows (top → bottom). Row keys:
    c       centre group, e.g. "C", "CH_2OH", "CH_2O{P}"
    l, r    left / right substituent (markup); lb, rb = bond order (1 default, 2 = double)
    u, d    substituent straight above / below the centre (ub, db = bond order); only on rows
            with no chain neighbour on that side, e.g. horizontal acetyl-CoA
            {c: C, l: CH_3, r: S-CoA, u: O, ub: 2}
    fork    only on the first row: {l, r, lb, rb} drawn up-left / up-right at 45°
            (O=C–H aldehyde, O=C–O⁻ carboxylate)
    dfork   only on the last row: the same, drawn down-left / down-right (Val CH(CH3)2, Asn amide)
    ring    only as the last row: "[*]c1ccccc1" — a ring hung below the chain (see ring.py)
    anchor  index of the character in `c` that bonds attach to (default: first 'C')
    down    bond order to the next row (default 1)

Decoration (comes from the reaction, so the same molecule can be plain elsewhere):
    colors  {row: colour}                            whole row
            {row: {c: colour, l: colour, ...}}       parts: c l r u d fork_l fork_r
            {row: {c: [colour, "{P}"]}}              only that fragment of the text
            {row: {c: [[cyan, "H"], [cyan, "H^+"]]}} several fragments
            {row: {ring: [{atoms: [7], color: magenta}]}}   ring row: only those SMILES atoms
                                                     (or {smarts: …}; see haworth.highlight_atoms)
            A bond is coloured when both of its ends have the same whole-part colour.
    marks   {row: {c: "*" | circle | tri | "1".."9", l: ...}}   isotope-tracing symbols;
            a list [tri, "2"] puts the first at the upper-left and the second at the upper-right
    boxes   [{rows: [2, 3], parts: ["4.c"], color: salmon, name: R}, …] background blocks,
            each exported as its own layer (box_<name>); a part box hugs just that text
    box     list of row indices, or "all"  → legacy single salmon box ("highlight" layer)
"""
import math

from .. import style
from ..graphic import BBox, Graphic, line
from ..text import TextLine, anchor_index, plain
from .ring import draw_ring

S = style.S
PITCH = 1.75 * S          # vertical distance between row centres
GAP_V = 0.20 * S          # gap between text and a vertical bond
GAP_H = 0.14 * S          # gap between text and a horizontal bond
BOND = 0.85 * S           # horizontal bond length
DBL = 0.085 * S           # half distance between the two lines of a double bond
FORK_DX, FORK_DY = 1.62 * S, 1.70 * S
FORK_TRIM = 0.58 * S
BOX_PAD = 0.32 * S


def _default_anchor(markup):
    return anchor_index(markup)


def _left_anchor(t):
    """Bond attaches to the last full-size letter of a left substituent ('HO' → O)."""
    idx = [i for i, (ch, _, dy, _) in enumerate(t.chars) if dy == 0 and ch.isascii() and ch.isalpha()]
    return idx[-1] if idx else len(t.chars) - 1


def _hbond(x_from, x_to, y, order, stroke):
    if order == 2:
        return line(x_from, y - DBL, x_to, y - DBL, stroke) + line(x_from, y + DBL, x_to, y + DBL, stroke)
    return line(x_from, y, x_to, y, stroke)


def _vbond(x, y_from, y_to, order, stroke):
    if order == 2:
        return line(x - DBL, y_from, x - DBL, y_to, stroke) + line(x + DBL, y_from, x + DBL, y_to, stroke)
    return line(x, y_from, x, y_to, stroke)


def _abond(x1, y1, x2, y2, order, stroke):
    if order != 2:
        return line(x1, y1, x2, y2, stroke)
    dx, dy = x2 - x1, y2 - y1
    n = math.hypot(dx, dy)
    ox, oy = -dy / n * DBL * 1.1, dx / n * DBL * 1.1
    return line(x1 + ox, y1 + oy, x2 + ox, y2 + oy, stroke) + line(x1 - ox, y1 - oy, x2 - ox, y2 - oy, stroke)


class Paint:
    """Resolves the `colors` decoration for (row, part)."""

    def __init__(self, colors):
        self.colors = {int(k): v for k, v in (colors or {}).items()}

    def _spec(self, i, part):
        v = self.colors.get(i)
        if isinstance(v, dict):
            return v.get(part)
        return v                                    # whole-row colour (or None)

    def whole(self, i, part):
        """Colour of the entire part, or None (fragment-only colouring counts as uncoloured)."""
        s = self._spec(i, part)
        return s if isinstance(s, str) else None

    def text(self, i, part, t):
        """(fill, spans) for drawing TextLine t."""
        s = self._spec(i, part)
        if isinstance(s, (list, tuple)):
            pairs = s if s and isinstance(s[0], (list, tuple)) else [s]
            spans, used = [], 0
            for col, frag in pairs:
                # search after the previous fragment so ["H", "H^+"] picks NADH's H, then H⁺
                sp = t.span(frag, start=used)
                if sp:
                    spans.append((sp[0], sp[1], col))
                    used = sp[1]
            return None, spans
        return s, None

    def bond(self, i, a, b):
        ca, cb = self.whole(i, a), self.whole(i, b)
        return ca if ca == cb else None


def _mark(sym, ax, ay, hw, right=False):
    """Small orange isotope symbol at the upper-left (or upper-right) of an atom centred on (ax, ay)."""
    col = style.color("orange")
    x, y = (ax + hw + 0.32 * S) if right else (ax - hw - 0.32 * S), ay - 0.62 * S
    w = style.LINE_W * 1.1
    if sym == "circle":
        return f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{0.17 * S:.2f}" fill="none" stroke="{col}" stroke-width="{w:.2f}"/>'
    if sym == "tri":
        r = 0.22 * S
        pts = " ".join(f"{x + r * math.cos(math.radians(a)):.2f},{y + r * math.sin(math.radians(a)):.2f}"
                       for a in (-90, 30, 150))
        return f'<polygon points="{pts}" fill="none" stroke="{col}" stroke-width="{w:.2f}" stroke-linejoin="round"/>'
    if sym.isdigit():
        t = TextLine(sym, "sans", 0.42 * S)
        circ = f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{0.25 * S:.2f}" fill="none" stroke="{col}" stroke-width="{w:.2f}"/>'
        return circ + t.svg(x - t.width / 2, y + t.cap / 2, "orange")
    t = TextLine(sym, "serif", 0.85 * S)            # "*"
    return t.svg(x - t.width / 2, y + t.cap * 0.9, "orange")


class _Placer:
    def __init__(self, paint, row_index, bb):
        self.paint, self.i, self.bb = paint, row_index, bb
        self.atoms = {}                              # part → (ax, ay, half width of anchor char)
        self.part_boxes = {}                         # part → tight (x0, y0, x1, y1) around its text

    def text(self, part, t, anchor_idx, ax, ay):
        """Place TextLine t so that its anchor character is centred on (ax, ay)."""
        x = ax - t.anchor_x(anchor_idx)
        base = ay + t.cap / 2
        self.bb.add(x, base + t.top, x + t.width, base + t.bottom)
        self.part_boxes[part] = (x, base - t.cap, x + t.char_box(len(t.chars) - 1)[1], base + t.size * 0.12)
        # isotope marks sit beside the atom that bears them: the bonding atom, except on left
        # substituents where it is the first letter (the C of CH_3, not its H)
        mi = 0 if part in ("l", "fork_l") else anchor_idx
        x0, x1 = t.char_box(mi)
        self.atoms[part] = (x + (x0 + x1) / 2, ay, (x1 - x0) / 2)
        fill, spans = self.paint.text(self.i, part, t)
        return t.svg(x, base, fill, spans=spans)


def draw_chain(rows, deco=None):
    deco = deco or {}
    paint = Paint(deco.get("colors"))
    marks = {int(k): v for k, v in (deco.get("marks") or {}).items()}
    boxes = list(deco.get("boxes") or [])
    if deco.get("box"):
        legacy = deco["box"]
        boxes.insert(0, {"rows": list(range(len(rows))) if legacy == "all" else legacy,
                         "color": "salmon", "name": "highlight"})

    parts, back = [], []
    total = BBox()
    row_boxes, row_svgs, centers, placers = [], [], [], []

    for i, row in enumerate(rows):
        y = i * PITCH
        rb = BBox()
        pl = _Placer(paint, i, rb)
        placers.append(pl)
        start = len(parts)
        if row.get("ring"):
            spec = paint._spec(i, "ring")
            hl = None
            if isinstance(spec, list):                # [{atoms: [7]} | {smarts: "c[OH]"}, …] → only those atoms
                from .haworth import highlight_atoms
                hl = highlight_atoms(row["ring"], spec)
            svg, rbb = draw_ring(row["ring"], 0, y, spec if isinstance(spec, str) else None, hl=hl)
            parts.append(svg)
            rb.add(rbb.x0, rbb.y0, rbb.x1, rbb.y1)
            row_boxes.append(rb)
            row_svgs.append("".join(parts[start:]))
            total.add(rb.x0, rb.y0, rb.x1, rb.y1)
            centers.append(None)
            continue
        ct = TextLine(row["c"])
        ai = row.get("anchor", _default_anchor(row["c"]))
        parts.append(pl.text("c", ct, ai, 0, y))
        # side bonds start at the edges of the whole centre group (CH_2—COO⁻), not just the anchor C
        cx_off = -ct.anchor_x(ai)
        c_left = ct.char_box(0)[0] + cx_off
        c_right = ct.char_box(len(ct.chars) - 1)[1] + cx_off
        centers.append(ct)

        if row.get("r"):
            t = TextLine(row["r"])
            xa = c_right + GAP_H
            xb = xa + BOND
            parts.append(_hbond(xa, xb, y, row.get("rb", 1), paint.bond(i, "c", "r")))
            ax = xb + GAP_H + (t.anchor_x(0) - t.char_box(0)[0])
            parts.append(pl.text("r", t, 0, ax, y))
        if row.get("l"):
            t = TextLine(row["l"])
            xa = c_left - GAP_H
            xb = xa - BOND
            parts.append(_hbond(xb, xa, y, row.get("lb", 1), paint.bond(i, "c", "l")))
            k = _left_anchor(t)
            # the text's right edge (including a trailing subscript, CH_3) sits GAP_H from the bond
            right_edge = t.char_box(len(t.chars) - 1)[1]
            ax = xb - GAP_H - (right_edge - t.anchor_x(k))
            parts.append(pl.text("l", t, k, ax, y))
        for side, sign in (("u", -1), ("d", 1)):
            if not row.get(side):
                continue
            t = TextLine(row[side])
            ay = y + sign * PITCH
            parts.append(pl.text(side, t, 0, 0, ay))
            ya = y + sign * (ct.cap / 2 + GAP_V)
            yb = ay - sign * (t.cap / 2 + GAP_V)
            parts.append(_vbond(0, ya, yb, row.get(side + "b", 1), paint.bond(i, "c", side)))
        for key, vdir in (("fork", -1), ("dfork", 1)):
            fork = row.get(key)
            if not fork:
                continue
            for side, sign in (("l", -1), ("r", 1)):
                if not fork.get(side):
                    continue
                t = TextLine(fork[side])
                k = _left_anchor(t) if side == "l" else 0
                fx, fy = sign * FORK_DX, y + vdir * FORK_DY
                parts.append(pl.text(key + "_" + side, t, k, fx, fy))
                n = math.hypot(fx, fy - y)
                ux, uy = fx / n, (fy - y) / n
                parts.append(_abond(ux * FORK_TRIM, y + uy * FORK_TRIM,
                                    fx - ux * FORK_TRIM, fy - uy * FORK_TRIM,
                                    fork.get(side + "b", 1), paint.bond(i, "c", key + "_" + side)))
        for part, sym in (marks.get(i) or {}).items():
            if part in pl.atoms:
                ax, ay, hw = pl.atoms[part]
                syms = sym if isinstance(sym, list) else [sym]
                # right-hand substituents get their mark on the right, clear of the centre's marks
                first_right = part in ("r", "fork_r")
                for k, one in enumerate(syms[:2]):
                    parts.append(_mark(str(one), ax, ay, hw, right=(k == 1) != first_right))
                rb.add(ax - hw - 0.6 * S, ay - 0.95 * S, ax + (hw + 0.6 * S if len(syms) > 1 else 0), ay)
        row_boxes.append(rb)
        row_svgs.append("".join(parts[start:]))
        total.add(rb.x0, rb.y0, rb.x1, rb.y1)

    # vertical bonds between consecutive rows
    for i in range(len(rows) - 1):
        ct, nt = centers[i], centers[i + 1]
        y_a = i * PITCH + ct.cap / 2 + GAP_V
        # a ring row starts at its (unlabelled) attachment carbon, so the bond runs right to it
        y_b = (i + 1) * PITCH - (nt.cap / 2 + GAP_V if nt else 0)
        c1, c2 = paint.whole(i, "c"), paint.whole(i + 1, "c")
        parts.append(_vbond(0, y_a, y_b, rows[i].get("down", 1), c1 if c1 == c2 else None))

    box_layers = {}
    for spec in boxes:
        hb = BBox()
        for i in spec.get("rows") or []:
            r = row_boxes[int(i)]
            hb.add(r.x0 - BOX_PAD, r.y0 - BOX_PAD, r.x1 + BOX_PAD, r.y1 + BOX_PAD)
        for ref in spec.get("parts") or []:
            ri, part = str(ref).split(".")
            pb = placers[int(ri)].part_boxes.get(part)
            if pb:
                pad = 0.08 * S
                hb.add(pb[0] - pad, pb[1] - pad, pb[2] + pad, pb[3] + pad)
        if hb.x0 == float("inf"):
            continue
        rect = (f'<rect x="{hb.x0:.2f}" y="{hb.y0:.2f}" width="{hb.x1 - hb.x0:.2f}" '
                f'height="{hb.y1 - hb.y0:.2f}" fill="{style.color(spec.get("color", "salmon"))}"/>')
        back.append(rect)
        total.add(hb.x0, hb.y0, hb.x1, hb.y1)
        name = spec.get("name", "R")
        box_layers["highlight" if name == "highlight" else f"box_{name}"] = Graphic(rect, hb.x0, hb.y0, hb.x1, hb.y1)

    last = row_boxes[-1]
    anchors = {"base": (0, last.y1), "axis": (0, 0)}
    sub = {f"row{i}": Graphic(svg, rb.x0, rb.y0, rb.x1, rb.y1) for i, (svg, rb) in enumerate(zip(row_svgs, row_boxes))}
    sub.update(box_layers)
    return Graphic("".join(back + parts), total.x0, total.y0, total.x1, total.y1, anchors, sub)


def draw_formula(markup, deco=None):
    """A molecule shown as plain formula text (e.g. CO_2, H_2O)."""
    deco = deco or {}
    t = TextLine(markup)
    bb = BBox()
    body = _Placer(Paint(deco.get("colors")), 0, bb).text("c", t, 0, t.anchor_x(0), 0)
    back = ""
    if deco.get("box"):
        x0, y0, x1, y1 = bb.x0 - BOX_PAD, bb.y0 - BOX_PAD, bb.x1 + BOX_PAD, bb.y1 + BOX_PAD
        back = (f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{x1 - x0:.2f}" height="{y1 - y0:.2f}" '
                f'fill="{style.BOX_COLOR}"/>')
        bb.add(x0, y0, x1, y1)
    sub = {"highlight": Graphic(back, x0, y0, x1, y1)} if back else {}
    return Graphic(back + body, bb.x0, bb.y0, bb.x1, bb.y1,
                   {"base": ((bb.x0 + bb.x1) / 2, bb.y1), "axis": (0, 0)}, sub)
