"""Illustration glyphs for `scene` figures (Alberts-style: pale flat fill + darker outline of the same
hue; `accent: true` adds a soft radial gradient for the protein the figure is about).

Every glyph takes the YAML item (positions and sizes in S units) and returns
    (svg, (x0, y0, x1, y1), (ax, ay))
in SVG units: markup, bounding box, and the anchor a label's leader line points at.
"""
import itertools
import math

from .. import style
from ..text import TextBlock

S = style.S
LW = style.LINE_W
_uid = itertools.count()

# Alberts-style flat palette: (fill, dark shade for loops / inner text, pale label-box background)
HUES = {
    "blue":   ("#5AA9DD", "#2E6E9E", "#C4E2F2"),
    "green":  ("#A9C96B", "#5E8A2E", "#DDE8C4"),
    "dkgreen": ("#4E9A51", "#2F6A33", "#C5E0C6"),
    "orange": ("#F0A35E", "#B5651F", "#F8D9BC"),
    "purple": ("#9C86C4", "#5F4A8C", "#DDD3EC"),
    "pink":   ("#E39A84", "#A9573F", "#F1C9BD"),
    "yellow": ("#F5D547", "#B08F10", "#FAEEB0"),
    "teal":   ("#63B9B3", "#2F7F7A", "#C8E6E3"),
    "gray":   ("#ABABAB", "#5F5F5F", "#E0E0E0"),
    "red":    ("#D23A3A", "#8E1F1F", "#F2C4C4"),
    "brown":  ("#B3966A", "#6F5634", "#E8DDCB"),
}
EDGE = "#FFFFFF"          # thin white seam between touching shapes (no dark outlines)
MEMBRANE_GRAY = "#CDD5CE"
BASES = {"A": "#5BB55B", "T": "#E2574C", "U": "#E2574C", "G": "#F2B632", "C": "#4C8BD6"}
PAIR = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A"}
STRANDS = {"parent": "#3F4F6B", "new": "#D9622B", "rna": "#C2338A", "primer": "#C2338A", "gray": "#8A8A8A"}

MEMBRANE_T = 2.6          # bilayer thickness in S (illustrations; schematics use 1.25)


def hue(name):
    """(fill, dark) for a palette name or a raw '#rrggbb'."""
    if (name or "").startswith("#"):
        return name, _mix(name, "#000000", 0.4)
    return HUES.get(name or "blue", HUES["blue"])[:2]


def label_bg(name):
    if (name or "").startswith("#"):
        return _mix(name, "#FFFFFF", 0.6)
    return HUES.get(name or "blue", HUES["blue"])[2]


def _mix(c1, c2, t):
    a = [int(c1[i:i + 2], 16) for i in (1, 3, 5)]
    b = [int(c2[i:i + 2], 16) for i in (1, 3, 5)]
    return "#" + "".join(f"{round(x + (y - x) * t):02X}" for x, y in zip(a, b))


def paint(it, default="blue"):
    """(fill, dark shade, defs): flat colour; `accent: true` deepens the fill a little."""
    fill, dark = hue(it.get("hue", default))
    if it.get("accent"):
        fill = _mix(fill, dark, 0.18)
    return fill, dark, ""


def P(v):
    return v * S


def _pt(it, key="at"):
    x, y = it[key]
    return P(x), P(y)


def _rot(svg, deg, cx, cy):
    return f'<g transform="rotate({deg:.1f},{cx:.2f},{cy:.2f})">{svg}</g>' if deg else svg


def _bb_rot(x0, y0, x1, y1, deg, cx, cy):
    if not deg:
        return x0, y0, x1, y1
    a = math.radians(deg)
    pts = []
    for x, y in ((x0, y0), (x1, y0), (x0, y1), (x1, y1)):
        dx, dy = x - cx, y - cy
        pts.append((cx + dx * math.cos(a) - dy * math.sin(a), cy + dx * math.sin(a) + dy * math.cos(a)))
    xs, ys = zip(*pts)
    return min(xs), min(ys), max(xs), max(ys)


def _inner_text(it, cx, cy, color=None, size=0.85):
    if not it.get("text"):
        return ""
    tb = TextBlock(it["text"], "sans", it.get("text_size", size) * S)
    return tb.svg(cx, cy - (tb.top + tb.bottom) / 2, color or style.INK)


def _shape_path(shape, cx, cy, w, h):
    """Outline path of a protein body centred on (cx, cy)."""
    rx, ry = w / 2, h / 2
    if shape == "round" or shape == "oval":
        return (f"M{cx - rx:.2f},{cy:.2f} a{rx:.2f},{ry:.2f} 0 1,0 {2 * rx:.2f},0 "
                f"a{rx:.2f},{ry:.2f} 0 1,0 {-2 * rx:.2f},0 Z")
    if shape == "box":
        r = min(rx, ry) * 0.45
        return (f"M{cx - rx + r:.2f},{cy - ry:.2f} H{cx + rx - r:.2f} Q{cx + rx:.2f},{cy - ry:.2f} {cx + rx:.2f},{cy - ry + r:.2f} "
                f"V{cy + ry - r:.2f} Q{cx + rx:.2f},{cy + ry:.2f} {cx + rx - r:.2f},{cy + ry:.2f} H{cx - rx + r:.2f} "
                f"Q{cx - rx:.2f},{cy + ry:.2f} {cx - rx:.2f},{cy + ry - r:.2f} V{cy - ry + r:.2f} "
                f"Q{cx - rx:.2f},{cy - ry:.2f} {cx - rx + r:.2f},{cy - ry:.2f} Z")
    if shape == "kidney":                       # bean with a dent on the top edge
        k = 0.55
        return (f"M{cx - rx:.2f},{cy:.2f} C{cx - rx:.2f},{cy - ry * 1.25:.2f} {cx - rx * 0.25:.2f},{cy - ry * 1.1:.2f} "
                f"{cx:.2f},{cy - ry * k:.2f} C{cx + rx * 0.25:.2f},{cy - ry * 1.1:.2f} {cx + rx:.2f},{cy - ry * 1.25:.2f} "
                f"{cx + rx:.2f},{cy:.2f} C{cx + rx:.2f},{cy + ry * 1.3:.2f} {cx - rx:.2f},{cy + ry * 1.3:.2f} {cx - rx:.2f},{cy:.2f} Z")
    if shape == "hand":                         # polymerase: palm with fingers / thumb around a cleft (opening right)
        return (f"M{cx + rx * 0.9:.2f},{cy - ry * 0.55:.2f} C{cx + rx * 0.2:.2f},{cy - ry * 1.25:.2f} {cx - rx * 1.05:.2f},{cy - ry * 1.15:.2f} "
                f"{cx - rx:.2f},{cy:.2f} C{cx - rx * 1.05:.2f},{cy + ry * 1.15:.2f} {cx + rx * 0.2:.2f},{cy + ry * 1.25:.2f} "
                f"{cx + rx * 0.9:.2f},{cy + ry * 0.55:.2f} C{cx + rx * 0.45:.2f},{cy + ry * 0.35:.2f} {cx + rx * 0.2:.2f},{cy + ry * 0.2:.2f} "
                f"{cx + rx * 0.2:.2f},{cy:.2f} C{cx + rx * 0.2:.2f},{cy - ry * 0.2:.2f} {cx + rx * 0.45:.2f},{cy - ry * 0.35:.2f} "
                f"{cx + rx * 0.9:.2f},{cy - ry * 0.55:.2f} Z")
    raise ValueError(f"未知的形狀 {shape}")


# ================================================================ proteins

def g_blob(it):
    """Generic protein: shape round | oval | box | kidney | hand | bilobe; w, h; rot; hue; accent; text."""
    cx, cy = _pt(it)
    w, h = P(it.get("w", 2.4)), P(it.get("h", 1.8))
    fill, line, defs = paint(it)
    shape = it.get("shape", "oval")
    sw = 1.15 * LW
    if shape == "bilobe":
        a = _shape_path("oval", cx - w * 0.22, cy, w * 0.62, h)
        b = _shape_path("oval", cx + w * 0.22, cy, w * 0.62, h * 0.9)
        body = (f'<path d="{a}" fill="{fill}" stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}"/>'
                f'<path d="{b}" fill="{fill}" stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}"/>')
    else:
        body = f'<path d="{_shape_path(shape, cx, cy, w, h)}" fill="{fill}" stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}"/>'
    deg = it.get("rot", 0)
    tx, ty = cx + P(it.get("text_dx", 0)), cy + P(it.get("text_dy", 0))
    svg = defs + _rot(body, deg, cx, cy) + _inner_text(it, tx, ty, line if it.get("text_color") == "hue" else None)
    return svg, _bb_rot(cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2, deg, cx, cy), (cx, cy)


def g_gpcr(it):
    """Seven transmembrane helices (rounded cylinders) with loops; `at` = centre of the membrane."""
    cx, cy = _pt(it)
    fill, line, defs = paint(it, "purple")
    n, cw, gap = 7, P(0.5), P(0.12)
    hgt = P(MEMBRANE_T + 0.5)
    total = n * cw + (n - 1) * gap
    x0 = cx - total / 2
    out = [defs]
    tops, bots = [], []
    for i in range(n):
        x = x0 + i * (cw + gap)
        out.append(f'<rect x="{x:.2f}" y="{cy - hgt / 2:.2f}" width="{cw:.2f}" height="{hgt:.2f}" rx="{cw / 2:.2f}" '
                   f'fill="{fill}" stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}"/>')
        tops.append((x + cw / 2, cy - hgt / 2 + P(0.1)))
        bots.append((x + cw / 2, cy + hgt / 2 - P(0.1)))
    loops = []
    for i in range(n - 1):              # extracellular loops between 1-2? no: N-term outside → loops alternate
        a, b = (bots[i], bots[i + 1]) if i % 2 == 0 else (tops[i], tops[i + 1])
        dy = P(0.55) * (1 if i % 2 == 0 else -1)
        loops.append(f"M{a[0]:.2f},{a[1]:.2f} C{a[0]:.2f},{a[1] + dy:.2f} {b[0]:.2f},{b[1] + dy:.2f} {b[0]:.2f},{b[1]:.2f}")
    nx, ny = tops[0]
    loops.append(f"M{nx:.2f},{ny:.2f} c{-P(0.2):.2f},{-P(0.6):.2f} {-P(0.8):.2f},{-P(0.5):.2f} {-P(1.0):.2f},{-P(0.9):.2f}")
    cx7, cy7 = bots[-1]
    loops.append(f"M{cx7:.2f},{cy7:.2f} c{P(0.2):.2f},{P(0.7):.2f} {P(0.9):.2f},{P(0.4):.2f} {P(1.1):.2f},{P(1.0):.2f}")
    out.insert(1, f'<path d="{" ".join(loops)}" fill="none" stroke="{fill}" stroke-width="{2.6 * LW:.2f}" stroke-linecap="round"/>')
    bb = (x0 - P(1.1), cy - hgt / 2 - P(1.0), x0 + total + P(1.2), cy + hgt / 2 + P(1.1))
    return "".join(out), bb, (cx, cy)


def g_channel(it):
    """Two subunits around a pore across the membrane; open: true/false; cargo drawn separately."""
    cx, cy = _pt(it)
    fill, line, defs = paint(it, "teal")
    w, h = P(it.get("w", 1.1)), P(it.get("h", MEMBRANE_T + 1.2))
    gap = P(0.55 if it.get("open", True) else 0.08)
    out = [defs]
    for side in (-1, 1):
        x = cx + side * (gap / 2 + w / 2)
        out.append(f'<path d="{_shape_path("box", x, cy, w, h)}" fill="{fill}" stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}"/>')
    out.append(_inner_text(it, cx, cy + h / 2 + P(0.7)))
    return "".join(out), (cx - gap / 2 - w, cy - h / 2, cx + gap / 2 + w, cy + h / 2), (cx, cy)


def g_transporter(it):
    """Carrier: two halves tilted into a V (open to `open: out|in`)."""
    cx, cy = _pt(it)
    fill, line, defs = paint(it, "blue")
    w, h = P(it.get("w", 1.25)), P(it.get("h", MEMBRANE_T + 1.0))
    tilt = 9 if it.get("open", "out") == "out" else -9
    out = [defs]
    for side in (-1, 1):
        x = cx + side * w * 0.52
        out.append(_rot(f'<path d="{_shape_path("box", x, cy, w, h)}" fill="{fill}" stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}"/>',
                        side * tilt, x, cy))
    return "".join(out), (cx - w * 1.2, cy - h / 2 - P(0.1), cx + w * 1.2, cy + h / 2 + P(0.1)), (cx, cy)


def g_pump(it):
    """P-type pump: transmembrane transporter body + cytosolic ATP-binding lobe (below)."""
    svg, bb, anchor = g_transporter({**it, "hue": it.get("hue", "orange"), "open": it.get("open", "in")})
    cx, cy = _pt(it)
    fill, line, defs = paint({**it, "hue": it.get("hue", "orange")})
    ly = bb[3] + P(0.55)
    lobe = f'<path d="{_shape_path("kidney", cx + P(0.3), ly, P(2.2), P(1.3))}" fill="{fill}" stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}"/>'
    return defs + lobe + svg, (min(bb[0], cx - P(0.8)), bb[1], max(bb[2], cx + P(1.4)), ly + P(0.9)), anchor


def g_gprotein(it):
    """Heterotrimeric G protein under the membrane: α (kidney, carries GDP/GTP) + β + γ.
    `split: true` draws α–GTP away from βγ."""
    cx, cy = _pt(it)
    fa, la, da = paint({**it, "hue": it.get("hue", "green")})
    fb, lb = hue("yellow")
    nuc = it.get("nucleotide", "GDP")
    split = P(1.4) if it.get("split") else 0
    out = [da]
    ax = cx - P(0.8) - split
    out.append(f'<path d="{_shape_path("kidney", ax, cy, P(2.0), P(1.5))}" fill="{fa}" stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}"/>')
    out.append(f'<path d="{_shape_path("oval", cx + P(0.85), cy, P(1.3), P(1.2))}" fill="{fb}" stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}"/>')
    out.append(f'<path d="{_shape_path("oval", cx + P(1.55), cy - P(0.55), P(0.55), P(0.9))}" fill="{fb}" stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}"/>')
    tb = TextBlock("α", "sans", 0.8 * S)
    out.append(tb.svg(ax, cy + P(0.15), la))
    tb = TextBlock("β", "sans", 0.75 * S)
    out.append(tb.svg(cx + P(0.85), cy + P(0.25), lb))
    tb = TextBlock("γ", "sans", 0.6 * S)
    out.append(tb.svg(cx + P(1.55), cy - P(0.42), lb))
    pill = TextBlock(nuc, "sans", 0.62 * S)
    px, py = ax, cy + P(0.95)
    pw = pill.width + P(0.4)
    col = "#D6177F" if nuc == "GTP" else "#777777"
    out.append(f'<rect x="{px - pw / 2:.2f}" y="{py - P(0.38):.2f}" width="{pw:.2f}" height="{P(0.68):.2f}" rx="{P(0.34):.2f}" '
               f'fill="#FFFFFF" stroke="{col}" stroke-width="{LW:.2f}"/>')
    out.append(pill.svg(px, py + P(0.2), col))
    return "".join(out), (ax - P(1.05), cy - P(1.0), cx + P(1.9), cy + P(1.4)), (ax, cy)


def g_ring(it):
    """Sliding clamp seen edge-on (kind: clamp) or a hexameric helicase ring (kind: hexamer)."""
    cx, cy = _pt(it)
    fill, line, defs = paint(it, "green")
    if it.get("kind", "clamp") == "hexamer":
        r, rr = P(it.get("r", 0.8)), P(0.5)
        out = [defs]
        for k in range(6):
            a = math.radians(60 * k - 90)
            out.append(f'<circle cx="{cx + r * math.cos(a):.2f}" cy="{cy + r * math.sin(a):.2f}" r="{rr:.2f}" '
                       f'fill="{fill}" stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}"/>')
        out.append(_inner_text(it, cx, cy, size=0.6))
        return "".join(out), (cx - r - rr, cy - r - rr, cx + r + rr, cy + r + rr), (cx, cy)
    rx, ry = P(it.get("w", 0.9)) / 2, P(it.get("h", 2.4)) / 2
    d = (f"M{cx - rx:.2f},{cy:.2f} a{rx:.2f},{ry:.2f} 0 1,0 {2 * rx:.2f},0 a{rx:.2f},{ry:.2f} 0 1,0 {-2 * rx:.2f},0 Z")
    return (defs + f'<path d="{d}" fill="{fill}" stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}" fill-opacity="0.85"/>',
            (cx - rx, cy - ry, cx + rx, cy + ry), (cx, cy))


# ================================================================ membranes

def g_membrane(it):
    """Straight lipid bilayer: y (centre), x0, x1; heads = circles, two tails each."""
    y, a0, a1 = P(it["y"]), P(it["x0"]), P(it["x1"])
    T = P(it.get("t", MEMBRANE_T))
    if not it.get("lipids"):
        return (f'<rect x="{a0:.2f}" y="{y - T / 2:.2f}" width="{a1 - a0:.2f}" height="{T:.2f}" fill="{it.get("color", MEMBRANE_GRAY)}"/>',
                (a0, y - T / 2, a1, y + T / 2), ((a0 + a1) / 2, y))
    r, step = P(0.24), P(0.56)
    head, head_line = it.get("head", "#F2C46D"), "#B8862B"
    out = [f'<rect x="{a0:.2f}" y="{y - T / 2 + r:.2f}" width="{a1 - a0:.2f}" height="{T - 2 * r:.2f}" fill="#FBF3DA"/>']
    n = max(2, int((a1 - a0) / step))
    tails = []
    heads = []
    for k in range(n + 1):
        x = a0 + (a1 - a0) * k / n
        for side in (-1, 1):
            hy = y + side * (T / 2 - r)
            ty = y + side * P(0.12)
            for dx in (-0.07, 0.07):
                tails.append(f"M{x + P(dx):.2f},{hy - side * r:.2f} L{x + P(dx) * 1.4:.2f},{ty:.2f}")
            heads.append(f'<circle cx="{x:.2f}" cy="{hy:.2f}" r="{r:.2f}" fill="{head}" stroke="{head_line}" stroke-width="{0.5 * LW:.2f}"/>')
    out.append(f'<path d="{" ".join(tails)}" stroke="#C9A85E" stroke-width="{0.55 * LW:.2f}" fill="none"/>')
    out += heads
    return "".join(out), (a0 - r, y - T / 2, a1 + r, y + T / 2), ((a0 + a1) / 2, y)


def g_region(it):
    """Compartment background: rounded box (x0, y0, x1, y1), pale fill, optional dashed outline."""
    x0, y0, x1, y1 = (P(it[k]) for k in ("x0", "y0", "x1", "y1"))
    fill, line = hue(it.get("hue", "gray"))
    op = it.get("opacity", 0.35)
    dash = ' stroke-dasharray="6 4"' if it.get("dashed") else ""
    svg = (f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{x1 - x0:.2f}" height="{y1 - y0:.2f}" rx="{P(0.8):.2f}" '
           f'fill="{fill}" fill-opacity="{op}" stroke="{line if it.get("outline", True) else "none"}" stroke-width="{LW:.2f}"{dash}/>')
    return svg, (x0, y0, x1, y1), (x0, y0)


# ================================================================ nucleic acids

def _smooth(pts):
    """Catmull–Rom through the points → cubic Bézier path."""
    if len(pts) == 2:
        (x0, y0), (x1, y1) = pts
        return f"M{x0:.2f},{y0:.2f} L{x1:.2f},{y1:.2f}"
    d = [f"M{pts[0][0]:.2f},{pts[0][1]:.2f}"]
    for i in range(len(pts) - 1):
        p0 = pts[i - 1] if i else pts[i]
        p1, p2 = pts[i], pts[i + 1]
        p3 = pts[i + 2] if i + 2 < len(pts) else p2
        c1 = (p1[0] + (p2[0] - p0[0]) / 6, p1[1] + (p2[1] - p0[1]) / 6)
        c2 = (p2[0] - (p3[0] - p1[0]) / 6, p2[1] - (p3[1] - p1[1]) / 6)
        d.append(f"C{c1[0]:.2f},{c1[1]:.2f} {c2[0]:.2f},{c2[1]:.2f} {p2[0]:.2f},{p2[1]:.2f}")
    return " ".join(d)


def _sample(pts, step):
    """Points every `step` along a polyline, with the unit tangent."""
    out = []
    segs = list(zip(pts, pts[1:]))
    carry = step / 2
    for (x0, y0), (x1, y1) in segs:
        L = math.hypot(x1 - x0, y1 - y0)
        if L == 0:
            continue
        ux, uy = (x1 - x0) / L, (y1 - y0) / L
        t = carry
        while t <= L:
            out.append((x0 + ux * t, y0 + uy * t, ux, uy))
            t += step
        carry = t - L
    return out


def _seq(it, n):
    s = (it.get("seq") or "ATGCGTACCGATTGCA").upper()
    return [s[i % len(s)] for i in range(n)]


def _end_label(text, x, y, ux, uy, color):
    tb = TextBlock(text, "sans", 0.62 * S)
    lx, ly = x + ux * P(0.55), y + uy * P(0.55)
    return tb.svg(lx, ly - (tb.top + tb.bottom) / 2, color), (lx - tb.width / 2, ly - P(0.4), lx + tb.width / 2, ly + P(0.4))


def g_strand(it):
    """Single nucleic-acid strand along `points`; kind parent | new | rna | primer | gray;
    bases: left | right (short coloured ticks); arrow: true (3′ end grows); ends: ["5'", "3'"]."""
    pts = [(P(x), P(y)) for x, y in it["points"]]
    col = STRANDS.get(it.get("kind", "parent"), it.get("kind"))
    w = P(it.get("width", 0.26))
    out = []
    xs, ys = zip(*pts)
    bb = [min(xs) - w, min(ys) - w, max(xs) + w, max(ys) + w]
    side = {"left": -1, "right": 1}.get(it.get("bases"), 0)
    if side:
        samples = _sample(pts, P(0.5))
        seq = _seq(it, len(samples))
        L = P(it.get("base_len", 0.6))
        for (x, y, ux, uy), b in zip(samples, seq):
            nx, ny = -uy * side, ux * side          # left of travel = (-uy, ux)
            out.append(f'<line x1="{x:.2f}" y1="{y:.2f}" x2="{x + nx * L:.2f}" y2="{y + ny * L:.2f}" '
                       f'stroke="{BASES[b]}" stroke-width="{P(0.13):.2f}" stroke-linecap="round"/>')
            bb = [min(bb[0], x + nx * L), min(bb[1], y + ny * L), max(bb[2], x + nx * L), max(bb[3], y + ny * L)]
    out.append(f'<path d="{_smooth(pts)}" fill="none" stroke="{col}" stroke-width="{w:.2f}" stroke-linecap="round" stroke-linejoin="round"/>')
    if it.get("arrow"):
        (x0, y0), (x1, y1) = pts[-2], pts[-1]
        L = math.hypot(x1 - x0, y1 - y0) or 1
        ux, uy = (x1 - x0) / L, (y1 - y0) / L
        hl, hw = P(0.55), P(0.32)
        tip = (x1 + ux * hl * 0.6, y1 + uy * hl * 0.6)
        out.append(f'<path d="M{tip[0]:.2f},{tip[1]:.2f} L{tip[0] - ux * hl - uy * hw:.2f},{tip[1] - uy * hl + ux * hw:.2f} '
                   f'L{tip[0] - ux * hl + uy * hw:.2f},{tip[1] - uy * hl - ux * hw:.2f} Z" fill="{col}"/>')
    ends = it.get("ends")
    if ends:
        for (txt, (x, y), (px, py)) in ((ends[0], pts[0], pts[1]), (ends[1], pts[-1], pts[-2])):
            if not txt:
                continue
            L = math.hypot(x - px, y - py) or 1
            svg, b = _end_label(txt, x, y, (x - px) / L, (y - py) / L, col)
            out.append(svg)
            bb = [min(bb[0], b[0]), min(bb[1], b[1]), max(bb[2], b[2]), max(bb[3], b[3])]
    mid = pts[len(pts) // 2]
    return "".join(out), tuple(bb), mid


def g_dsdna(it):
    """Double-stranded DNA from `from` to `to`: two backbones with coloured base-pair rungs.
    kinds: [top strand kind, bottom strand kind]; ends: true labels 5′/3′."""
    (x0, y0), (x1, y1) = _pt(it, "from"), _pt(it, "to")
    L = math.hypot(x1 - x0, y1 - y0)
    ux, uy = (x1 - x0) / L, (y1 - y0) / L
    nx, ny = -uy, ux
    half = P(it.get("width", 1.3)) / 2
    k1, k2 = it.get("kinds", ["parent", "parent"])
    out = []
    samples = _sample([(x0, y0), (x1, y1)], P(0.5))
    seq = _seq(it, len(samples))
    for (x, y, _, _), b in zip(samples, seq):
        ax, ay = x - nx * half, y - ny * half
        bx, by = x + nx * half, y + ny * half
        mx, my = (ax + bx) / 2, (ay + by) / 2
        out.append(f'<line x1="{ax:.2f}" y1="{ay:.2f}" x2="{mx:.2f}" y2="{my:.2f}" stroke="{BASES[b]}" stroke-width="{P(0.16):.2f}"/>')
        out.append(f'<line x1="{mx:.2f}" y1="{my:.2f}" x2="{bx:.2f}" y2="{by:.2f}" stroke="{BASES[PAIR[b]]}" stroke-width="{P(0.16):.2f}"/>')
    w = P(0.26)
    for sgn, k in ((-1, k1), (1, k2)):
        c = STRANDS.get(k, k)
        out.append(f'<line x1="{x0 + sgn * nx * half:.2f}" y1="{y0 + sgn * ny * half:.2f}" x2="{x1 + sgn * nx * half:.2f}" '
                   f'y2="{y1 + sgn * ny * half:.2f}" stroke="{c}" stroke-width="{w:.2f}" stroke-linecap="round"/>')
    bb = [min(x0, x1) - half - w, min(y0, y1) - half - w, max(x0, x1) + half + w, max(y0, y1) + half + w]
    if it.get("ends"):
        for sgn, k, a, b in ((-1, k1, "5'", "3'"), (1, k2, "3'", "5'")):
            c = STRANDS.get(k, k)
            for txt, (x, y), d in ((a, (x0, y0), -1), (b, (x1, y1), 1)):
                svg, lb = _end_label(txt, x + sgn * nx * half, y + sgn * ny * half, ux * d, uy * d, c)
                out.append(svg)
                bb = [min(bb[0], lb[0]), min(bb[1], lb[1]), max(bb[2], lb[2]), max(bb[3], lb[3])]
    return "".join(out), tuple(bb), ((x0 + x1) / 2, (y0 + y1) / 2)


# ================================================================ small things

def g_ion(it):
    """Coloured sphere with its symbol (Na⁺, K⁺, cAMP …); r in S."""
    cx, cy = _pt(it)
    fill, line, defs = paint(it, "yellow")
    r = P(it.get("r", 0.45))
    out = [defs, f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{r:.2f}" fill="{fill}" stroke="{line}" stroke-width="{0.8 * LW:.2f}"/>']
    if it.get("text"):
        tb = TextBlock(it["text"], "sans", it.get("text_size", 0.55) * S)
        out.append(tb.svg(cx, cy - (tb.top + tb.bottom) / 2, style.INK))
    return "".join(out), (cx - r, cy - r, cx + r, cy + r), (cx, cy)


def g_ligand(it):
    """Signal molecule: shape triangle | diamond | hexagon."""
    cx, cy = _pt(it)
    fill, line, defs = paint(it, "red")
    r = P(it.get("r", 0.5))
    n = {"triangle": 3, "diamond": 4, "hexagon": 6}[it.get("shape", "triangle")]
    pts = " ".join(f"{cx + r * math.cos(math.radians(360 * k / n - 90)):.2f},{cy + r * math.sin(math.radians(360 * k / n - 90)):.2f}"
                   for k in range(n))
    return (defs + f'<polygon points="{pts}" fill="{fill}" stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}" stroke-linejoin="round"/>',
            (cx - r, cy - r, cx + r, cy + r), (cx, cy))


def g_arrow(it):
    """from → to; curve (bend, + = left of travel); style thin | fat; color; inhibit: true (⊣)."""
    (x0, y0), (x1, y1) = _pt(it, "from"), _pt(it, "to")
    L = math.hypot(x1 - x0, y1 - y0)
    ux, uy = (x1 - x0) / L, (y1 - y0) / L
    bend = it.get("curve", 0) * L
    mx, my = (x0 + x1) / 2 - uy * bend, (y0 + y1) / 2 + ux * bend
    col = style.color(it.get("color", "black"))
    if it.get("style") == "fat":               # Alberts' broad process arrow (straight)
        w, hl, hw = P(0.5), P(1.1), P(0.95)
        bx, by = x1 - ux * hl, y1 - uy * hl
        pts = [(x0 - uy * w / 2, y0 + ux * w / 2), (bx - uy * w / 2, by + ux * w / 2), (bx - uy * hw, by + ux * hw), (x1, y1),
               (bx + uy * hw, by - ux * hw), (bx + uy * w / 2, by - ux * w / 2), (x0 + uy * w / 2, y0 - ux * w / 2)]
        fill = it.get("fill", "#D9D9D9")
        poly = " ".join(f"{x:.2f},{y:.2f}" for x, y in pts)
        return (f'<polygon points="{poly}" fill="{fill}" stroke="#8C8C8C" stroke-width="{0.8 * LW:.2f}" stroke-linejoin="round"/>',
                (min(x0, x1) - hw, min(y0, y1) - hw, max(x0, x1) + hw, max(y0, y1) + hw), (mx, my))
    # tangent at the end of a quadratic curve
    tx, ty = x1 - mx, y1 - my
    tl = math.hypot(tx, ty) or 1
    tx, ty = tx / tl, ty / tl
    w = LW * it.get("width", 1.3)
    if it.get("inhibit"):
        bar = P(0.4)
        head = f'<line x1="{x1 - ty * bar:.2f}" y1="{y1 + tx * bar:.2f}" x2="{x1 + ty * bar:.2f}" y2="{y1 - tx * bar:.2f}" stroke="{col}" stroke-width="{2 * w:.2f}"/>'
        ex, ey = x1, y1
    else:
        hl, hw = P(0.5), P(0.24)
        head = (f'<path d="M{x1:.2f},{y1:.2f} L{x1 - tx * hl - ty * hw:.2f},{y1 - ty * hl + tx * hw:.2f} '
                f'L{x1 - tx * hl + ty * hw:.2f},{y1 - ty * hl - tx * hw:.2f} Z" fill="{col}"/>')
        ex, ey = x1 - tx * hl * 0.8, y1 - ty * hl * 0.8
    dash = ' stroke-dasharray="5 4"' if it.get("dashed") else ""
    body = (f'<path d="M{x0:.2f},{y0:.2f} Q{2 * mx - (x0 + x1) / 2:.2f},{2 * my - (y0 + y1) / 2:.2f} {ex:.2f},{ey:.2f}" '
            f'fill="none" stroke="{col}" stroke-width="{w:.2f}" stroke-linecap="round"{dash}/>')
    return body + head, (min(x0, x1, mx) - P(0.3), min(y0, y1, my) - P(0.3), max(x0, x1, mx) + P(0.3), max(y0, y1, my) + P(0.3)), (mx, my)


def g_text(it):
    """Free text: font sans | serif, size (S), color, align."""
    x, y = _pt(it)
    font = it.get("font", "sans")
    tb = TextBlock(it["text"], font, it.get("size", 0.9) * S, align=it.get("align", "center"))
    dy = -(tb.top + tb.bottom) / 2
    w = tb.width
    x0 = x - w / 2 if tb.align == "center" else (x if tb.align == "left" else x - w)
    return tb.svg(x, y + dy, style.color(it.get("color"))), (x0, y + dy + tb.top, x0 + w, y + dy + tb.bottom), (x, y)


GLYPHS = {
    "blob": g_blob, "gpcr": g_gpcr, "channel": g_channel, "transporter": g_transporter, "pump": g_pump,
    "gprotein": g_gprotein, "ring": g_ring, "membrane": g_membrane, "region": g_region,
    "strand": g_strand, "dsdna": g_dsdna, "ion": g_ion, "ligand": g_ligand, "arrow": g_arrow, "text": g_text,
}


# ================================================================ DNA topology and chromatin

def g_helix(it):
    """B-DNA side view from `from` to `to`: two sinusoidal backbones offset so that the major and minor
    grooves differ; coloured rungs. pitch (S per turn), amp (half width), kinds [strand1, strand2]."""
    (x0, y0), (x1, y1) = _pt(it, "from"), _pt(it, "to")
    L = math.hypot(x1 - x0, y1 - y0)
    ux, uy = (x1 - x0) / L, (y1 - y0) / L
    nx, ny = -uy, ux
    pitch, amp = P(it.get("pitch", 4.2)), P(it.get("amp", 0.75))
    off = it.get("groove", 0.36)                       # fraction of a turn between the strands → minor groove
    k1, k2 = it.get("kinds", ["parent", "parent"])
    n = max(8, int(L / P(0.08)))

    def pt(t, phase):
        a = 2 * math.pi * (t / pitch - phase)
        return x0 + ux * t + nx * amp * math.sin(a), y0 + uy * t + ny * amp * math.sin(a), math.cos(a)

    A = [pt(L * i / n, 0) for i in range(n + 1)]
    B = [pt(L * i / n, off) for i in range(n + 1)]
    out = []
    step = pitch / 10.5                                # 10.5 bp per turn
    seq = _seq(it, int(L / step) + 1)
    for i, b in enumerate(seq):
        t = step * (i + 0.5)
        if t > L:
            break
        (ax, ay, _), (bx, by, _) = pt(t, 0), pt(t, off)
        mx, my = (ax + bx) / 2, (ay + by) / 2
        w = P(0.12)
        out.append(f'<line x1="{ax:.2f}" y1="{ay:.2f}" x2="{mx:.2f}" y2="{my:.2f}" stroke="{BASES[b]}" stroke-width="{w:.2f}"/>')
        out.append(f'<line x1="{mx:.2f}" y1="{my:.2f}" x2="{bx:.2f}" y2="{by:.2f}" stroke="{BASES[PAIR[b]]}" stroke-width="{w:.2f}"/>')
    W = P(0.3)

    def runs(pts, front):
        """split a backbone into segments where it is in front (cos > 0) or behind."""
        segs, cur = [], []
        for x, y, c in pts:
            if (c > 0) == front:
                cur.append((x, y))
            elif cur:
                segs.append(cur)
                cur = []
        if cur:
            segs.append(cur)
        return segs

    for front in (False, True):                        # back halves first, lighter
        for pts, k in ((A, k1), (B, k2)):
            col = STRANDS.get(k, k)
            col = col if front else _mix(col, "#FFFFFF", 0.45)
            for seg in runs(pts, front):
                if len(seg) > 1:
                    d = "M" + " L".join(f"{x:.2f},{y:.2f}" for x, y in seg)
                    out.append(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="{W:.2f}" stroke-linecap="round"/>')
    pad = amp + W
    bb = (min(x0, x1) - abs(nx) * pad - W, min(y0, y1) - abs(ny) * pad - W, max(x0, x1) + abs(nx) * pad + W, max(y0, y1) + abs(ny) * pad + W)
    return "".join(out), bb, ((x0 + x1) / 2, (y0 + y1) / 2)


def _ring_path(cx, cy, rx, ry):
    return f"M{cx - rx:.2f},{cy:.2f} a{rx:.2f},{ry:.2f} 0 1,0 {2 * rx:.2f},0 a{rx:.2f},{ry:.2f} 0 1,0 {-2 * rx:.2f},0 Z"


def g_plasmid(it):
    """Circular DNA: kind relaxed | supercoiled (twists: n crossings) | catenane | theta (replication bubble,
    `bubble: [deg0, deg1]` measured clockwise from the top)."""
    cx, cy = _pt(it)
    r = P(it.get("r", 1.6))
    kind = it.get("kind", "relaxed")
    par, new = STRANDS["parent"], STRANDS["new"]
    w = P(0.24)
    out = []
    if kind == "relaxed":
        out.append(f'<path d="{_ring_path(cx, cy, r, r)}" fill="none" stroke="{par}" stroke-width="{w:.2f}"/>')
        out.append(f'<path d="{_ring_path(cx, cy, r - P(0.32), r - P(0.32))}" fill="none" stroke="{par}" stroke-width="{w:.2f}"/>')
        return "".join(out), (cx - r - w, cy - r - w, cx + r + w, cy + r + w), (cx, cy - r)
    if kind == "supercoiled":                          # interwound (plectonemic) rod with `twists` crossings
        n, h, hw = it.get("twists", 2), 2 * r, r * 0.42
        m = 60
        A = [(cx + hw * math.sin(math.pi * (n + 1) * t / 1), cy - h / 2 + h * t) for t in (i / m for i in range(m + 1))]
        B = [(cx - hw * math.sin(math.pi * (n + 1) * t / 1), cy - h / 2 + h * t) for t in (i / m for i in range(m + 1))]
        d = _smooth(A[::3] + [A[-1]]) + " " + _smooth(B[::3] + [B[-1]])
        out.append(f'<path d="{d}" fill="none" stroke="{par}" stroke-width="{2.2 * w:.2f}" stroke-linecap="round"/>')
        out.append(f'<path d="{d}" fill="none" stroke="#FFFFFF" stroke-width="{0.5 * w:.2f}" stroke-linecap="round"/>')
        return "".join(out), (cx - hw - w, cy - h / 2 - w, cx + hw + w, cy + h / 2 + w), (cx, cy - h / 2)
    if kind == "catenane":
        d = 0.62 * r
        for dx in (-d, d):
            out.append(f'<path d="{_ring_path(cx + dx, cy, r, r)}" fill="none" stroke="{par}" stroke-width="{2 * w:.2f}"/>')
        # left ring passes over the right one at the upper crossing
        a = math.acos(d / r)
        x0, y0 = cx - d + r * math.cos(a - 0.35), cy - r * math.sin(a - 0.35)
        x1, y1 = cx - d + r * math.cos(a + 0.35), cy - r * math.sin(a + 0.35)
        arc = f"M{x0:.2f},{y0:.2f} A{r:.2f},{r:.2f} 0 0,0 {x1:.2f},{y1:.2f}"
        out.append(f'<path d="{arc}" fill="none" stroke="#FFFFFF" stroke-width="{3.4 * w:.2f}"/>')
        out.append(f'<path d="{arc}" fill="none" stroke="{par}" stroke-width="{2 * w:.2f}"/>')
        return "".join(out), (cx - d - r - w, cy - r - w, cx + d + r + w, cy + r + w), (cx, cy - r)
    if kind == "theta":
        a0, a1 = (math.radians(v) for v in it.get("bubble", [-50, 50]))
        out.append(f'<path d="{_ring_path(cx, cy, r, r)}" fill="none" stroke="{par}" stroke-width="{2 * w:.2f}"/>')
        ri = r - P(0.55)

        def pol(rr, a):
            return cx + rr * math.sin(a), cy - rr * math.cos(a)
        (sx, sy), (ex, ey) = pol(ri, a0), pol(ri, a1)
        large = 1 if a1 - a0 > math.pi else 0
        out.append(f'<path d="M{sx:.2f},{sy:.2f} A{ri:.2f},{ri:.2f} 0 {large},1 {ex:.2f},{ey:.2f}" fill="none" '
                   f'stroke="{par}" stroke-width="{2 * w:.2f}"/>')
        for rr, col in ((r - P(0.2), new), (ri + P(0.2), new)):
            (sx, sy), (ex, ey) = pol(rr, a0 + 0.05), pol(rr, a1 - 0.05)
            out.append(f'<path d="M{sx:.2f},{sy:.2f} A{rr:.2f},{rr:.2f} 0 {large},1 {ex:.2f},{ey:.2f}" fill="none" '
                       f'stroke="{col}" stroke-width="{1.2 * w:.2f}"/>')
        for a in (a0, a1):
            fx, fy = pol(r - P(0.28), a)
            out.append(f'<circle cx="{fx:.2f}" cy="{fy:.2f}" r="{P(0.22):.2f}" fill="{HUES["orange"][0]}"/>')
        return "".join(out), (cx - r - w, cy - r - w, cx + r + w, cy + r + w), pol(r, (a0 + a1) / 2)
    raise ValueError(f"未知的 plasmid kind {kind}")


def g_nucleosome(it):
    """Histone octamer (flat disk) with DNA wrapped ~1.7 turns; `at`, r (S)."""
    cx, cy = _pt(it)
    r = P(it.get("r", 0.9))
    fill, dark, _ = paint(it, "orange")
    dna = STRANDS["parent"]
    w = P(0.2)
    rx, ry = r, r * 0.62
    back = f'<path d="M{cx - rx:.2f},{cy:.2f} A{rx:.2f},{ry:.2f} 0 0,1 {cx + rx:.2f},{cy:.2f}" fill="none" stroke="{_mix(dna, "#FFFFFF", 0.45)}" stroke-width="{w:.2f}"/>'
    disk = (f'<path d="{_ring_path(cx, cy + P(0.12), rx * 0.86, ry * 1.25)}" fill="{dark}"/>'
            f'<path d="{_ring_path(cx, cy, rx * 0.86, ry * 1.25)}" fill="{fill}" stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}"/>')
    front = (f'<path d="M{cx - rx:.2f},{cy - P(0.15):.2f} C{cx - rx * 0.6:.2f},{cy + ry * 1.25:.2f} {cx + rx * 0.6:.2f},{cy + ry * 1.25:.2f} '
             f'{cx + rx:.2f},{cy - P(0.15):.2f}" fill="none" stroke="{dna}" stroke-width="{w:.2f}" stroke-linecap="round"/>'
             f'<path d="M{cx - rx:.2f},{cy + P(0.25):.2f} C{cx - rx * 0.6:.2f},{cy + ry * 1.7:.2f} {cx + rx * 0.6:.2f},{cy + ry * 1.7:.2f} '
             f'{cx + rx:.2f},{cy + P(0.25):.2f}" fill="none" stroke="{dna}" stroke-width="{w:.2f}" stroke-linecap="round"/>')
    return back + disk + front, (cx - rx - w, cy - ry * 1.25 - w, cx + rx + w, cy + ry * 1.8 + w), (cx, cy)


def g_chromosome(it):
    """Metaphase chromosome: two sister chromatids joined at the centromere; `at`, h (S)."""
    cx, cy = _pt(it)
    h = P(it.get("h", 5.0))
    fill, dark, _ = paint(it, "purple")
    w = h * 0.2
    c = h * 0.08                                       # centromere offset above centre
    out = []
    for side in (-1, 1):
        x = cx + side * w * 0.55
        for y0, y1 in ((cy - h / 2, cy - c), (cy - c, cy + h / 2)):
            tilt = side * (6 if y1 <= cy - c + 1e-6 else -6)
            ym = (y0 + y1) / 2
            seg = _shape_path("box", x, ym, w, (y1 - y0) * 1.02)
            out.append(_rot(f'<path d="{seg}" fill="{fill}" stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}"/>', tilt, x, y1 if y1 <= cy - c + 1e-6 else y0))
    out.append(f'<ellipse cx="{cx:.2f}" cy="{cy - c:.2f}" rx="{w * 0.55:.2f}" ry="{w * 0.3:.2f}" fill="{dark}"/>')
    return "".join(out), (cx - 1.4 * w, cy - h / 2, cx + 1.4 * w, cy + h / 2), (cx, cy - c)


GLYPHS.update({"helix": g_helix, "plasmid": g_plasmid, "nucleosome": g_nucleosome, "chromosome": g_chromosome})


def g_mark(it):
    """Small marks on DNA: kind lesion (red star) | cut (red wedge pointing at `at`, from `side`: up/down) |
    x (red ✕, blocked) | methyl (CH₃ tag on a stem; side up/down)."""
    cx, cy = _pt(it)
    kind = it.get("kind", "lesion")
    red = style.color(it.get("color", "red"))
    if kind == "lesion":
        r1, r2, n = P(it.get("r", 0.45)), P(it.get("r", 0.45)) * 0.45, 6
        pts = []
        for k in range(2 * n):
            r = r1 if k % 2 == 0 else r2
            a = math.pi * k / n - math.pi / 2
            pts.append(f"{cx + r * math.cos(a):.2f},{cy + r * math.sin(a):.2f}")
        return (f'<polygon points="{" ".join(pts)}" fill="{red}" stroke="{EDGE}" stroke-width="{0.6 * LW:.2f}"/>',
                (cx - r1, cy - r1, cx + r1, cy + r1), (cx, cy))
    if kind == "cut":
        s = -1 if it.get("side", "up") == "up" else 1
        h, w = P(0.7), P(0.32)
        tip = (cx, cy + s * P(0.15))
        poly = f"{tip[0]:.2f},{tip[1]:.2f} {cx - w:.2f},{cy + s * h:.2f} {cx + w:.2f},{cy + s * h:.2f}"
        return (f'<polygon points="{poly}" fill="{red}"/>', (cx - w, min(cy, cy + s * h), cx + w, max(cy, cy + s * h)), (cx, cy))
    if kind == "x":
        d = P(0.35)
        return (f'<path d="M{cx - d:.2f},{cy - d:.2f} L{cx + d:.2f},{cy + d:.2f} M{cx + d:.2f},{cy - d:.2f} L{cx - d:.2f},{cy + d:.2f}" '
                f'stroke="{red}" stroke-width="{2.2 * LW:.2f}" stroke-linecap="round"/>', (cx - d, cy - d, cx + d, cy + d), (cx, cy))
    if kind == "methyl":
        s = -1 if it.get("side", "up") == "up" else 1
        tb = TextBlock("CH_3", "serif", 0.62 * S)
        ty = cy + s * P(0.95)
        svg = (f'<line x1="{cx:.2f}" y1="{cy:.2f}" x2="{cx:.2f}" y2="{cy + s * P(0.45):.2f}" stroke="{style.INK}" stroke-width="{LW:.2f}"/>'
               + tb.svg(cx, ty - (tb.top + tb.bottom) / 2, style.color(it.get("color", "black"))))
        return svg, (cx - tb.width / 2, min(cy, ty - P(0.35)), cx + tb.width / 2, max(cy, ty + P(0.35))), (cx, cy)
    raise ValueError(f"未知的 mark kind {kind}")


GLYPHS["mark"] = g_mark


def _loop(ex, ey, dx, dy, g, r, n=26):
    """Points of a hairpin loop closing a stem whose end centre is (ex, ey), pointing (dx, dy), half-gap g."""
    h = math.sqrt(max(r * r - g * g, 0))
    cx, cy = ex + dx * h, ey + dy * h
    nx, ny = -dy, dx
    aA = math.atan2(ey + ny * g - cy, ex + nx * g - cx)
    aB = math.atan2(ey - ny * g - cy, ex - nx * g - cx)
    ad = math.atan2(dy, dx)
    delta = (aB - aA) % (2 * math.pi)
    span = delta if (ad - aA) % (2 * math.pi) < delta else -(2 * math.pi - delta)
    return [(cx + r * math.cos(aA + span * k / n), cy + r * math.sin(aA + span * k / n)) for k in range(n + 1)], (cx, cy)


def _trna_pts(it):
    """Backbone points of the cloverleaf (5′ → 3′, incl. the CCA tail) and the loop centres; SVG units."""
    cx, cy = _pt(it)
    u = P(it.get("s", 1.0))
    g, r = 0.22 * u, 0.6 * u
    top, yD = cy - 2.5 * u, cy - 0.15 * u
    pts, loops = [], {}
    pts += [(cx - g, top), (cx - g, cy - 0.55 * u), (cx - 0.55 * u, yD - g), (cx - 1.8 * u, yD - g)]
    lp, loops["D"] = _loop(cx - 1.8 * u, yD, -1, 0, g, r)
    pts += lp + [(cx - 0.55 * u, yD + g), (cx - g, cy + 0.55 * u), (cx - g, cy + 1.8 * u)]
    lp, loops["AC"] = _loop(cx, cy + 1.8 * u, 0, 1, g, r * 1.1)
    pts += lp + [(cx + g, cy + 0.55 * u), (cx + 0.55 * u, yD + g), (cx + 1.8 * u, yD + g)]
    lp, loops["T"] = _loop(cx + 1.8 * u, yD, 1, 0, g, r)
    pts += lp + [(cx + 0.55 * u, yD - g), (cx + g, cy - 0.55 * u), (cx + g, top)]
    if it.get("cca", True):
        pts.append((cx + g, top - 0.7 * u))
    return pts, loops


def g_trna(it):
    """tRNA cloverleaf (acceptor stem up): `at` = junction centre, `s` scale, `anticodon` text, `cca` (3′ CCA tail),
    `leader` / `trailer` (S, precursor ends in gray), `intron` (extra loop in the anticodon loop, red), `kind` strand colour.
    `shape: icon` = compact Alberts tRNA standing on the mRNA (`at` = anticodon foot): see _trna_icon."""
    if it.get("shape") == "icon":
        return _trna_icon(it)
    cx, cy = _pt(it)
    u = P(it.get("s", 1.0))
    col = STRANDS.get(it.get("kind", "rna"), it.get("kind"))
    g, w, r = 0.22 * u, 0.2 * u, 0.6 * u
    top, yD = cy - 2.5 * u, cy - 0.15 * u
    rungs = []

    def stem(x0, y0, x1, y1):                        # rungs between the two strands of a stem
        L = math.hypot(x1 - x0, y1 - y0)
        k = int(L / (0.38 * u))
        for i in range(1, k + 1):
            t = i / (k + 1)
            rungs.append((x0 + (x1 - x0) * t, y0 + (y1 - y0) * t))

    pts, loops = _trna_pts(it)
    stem(cx, top, cx, cy - 0.55 * u)
    stem(cx - 0.55 * u, yD, cx - 1.8 * u, yD)
    stem(cx, cy + 0.55 * u, cx, cy + 1.8 * u)
    stem(cx + 0.55 * u, yD, cx + 1.8 * u, yD)
    out = [f'<line x1="{x - (g if abs(x - cx) < 1e-6 else 0):.2f}" y1="{y - (0 if abs(x - cx) < 1e-6 else g):.2f}" '
           f'x2="{x + (g if abs(x - cx) < 1e-6 else 0):.2f}" y2="{y + (0 if abs(x - cx) < 1e-6 else g):.2f}" '
           f'stroke="#BBBBBB" stroke-width="{0.09 * u:.2f}"/>' for x, y in rungs]
    tail = pts[-1] if it.get("cca", True) else None
    out.append(f'<path d="M' + " L".join(f"{x:.2f},{y:.2f}" for x, y in pts) +
               f'" fill="none" stroke="{col}" stroke-width="{w:.2f}" stroke-linejoin="round" stroke-linecap="round"/>')
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    bb = [min(xs) - w, min(ys) - w, max(xs) + w, max(ys) + w]
    gray = STRANDS["gray"]
    if tail:
        tb = TextBlock("CCA", "sans", 0.55 * u)
        tx = tail[0] + 0.25 * u + tb.width / 2
        out.append(tb.svg(tx, tail[1] + 0.3 * u - (tb.top + tb.bottom) / 2, col))
        bb[2] = max(bb[2], tx + tb.width / 2)
    for key, (x0, y0), sx in (("leader", (cx - g, top), -1), ("trailer", tail or (cx + g, top), 1)):
        if it.get(key):
            L = it[key] * u
            x1, y1 = x0 + sx * L * 0.75, y0 - L * 0.66
            out.append(f'<path d="M{x0:.2f},{y0:.2f} Q{x0:.2f},{y1:.2f} {x1:.2f},{y1:.2f}" fill="none" stroke="{gray}" '
                       f'stroke-width="{w:.2f}" stroke-linecap="round"/>')
            bb = [min(bb[0], x1 - w), min(bb[1], y1 - w), max(bb[2], x1 + w), bb[3]]
    ax, ay = loops["AC"]
    if it.get("intron"):
        ib = ay + 1.15 * r
        out.append(f'<path d="M{ax - 0.3 * u:.2f},{ib - 0.1 * u:.2f} C{ax - 0.9 * u:.2f},{ib + 1.3 * u:.2f} {ax + 0.9 * u:.2f},{ib + 1.3 * u:.2f} '
                   f'{ax + 0.3 * u:.2f},{ib - 0.1 * u:.2f}" fill="none" stroke="{style.color("red")}" stroke-width="{w:.2f}" stroke-linecap="round"/>')
        bb[3] = max(bb[3], ib + 1.05 * u)
    elif it.get("anticodon"):
        tb = TextBlock(it["anticodon"], "sans", 0.62 * u)
        ty = ay + 1.1 * r + 0.5 * u
        out.append(tb.svg(ax, ty - (tb.top + tb.bottom) / 2, style.INK))
        bb[3] = max(bb[3], ty + 0.4 * u)
        bb[0], bb[2] = min(bb[0], ax - tb.width / 2), max(bb[2], ax + tb.width / 2)
    return "".join(out), tuple(bb), (cx, cy)


GLYPHS["trna"] = g_trna


AA_HUES = ["orange", "green", "blue", "pink", "yellow", "teal", "purple", "brown"]


def _trna_icon(it):
    """Compact tRNA: rounded body whose foot (`at`) sits on the mRNA, `anticodon` letters at the foot,
    `aa` (text) = one amino-acid ball on top, `peptide` (n) = a chain of n balls instead; `s` scale; `hue` body."""
    cx, cy = _pt(it)
    u = P(it.get("s", 1.0))
    fill, dark = hue(it.get("hue", "#E7A9C9"))
    w, h = 1.1 * u, 2.5 * u
    x0, y0 = cx - w / 2, cy - 0.15 * u - h
    out = [f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{w:.2f}" height="{h:.2f}" rx="{0.3 * u:.2f}" fill="{fill}" '
           f'stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}"/>']
    bb = [x0, y0, x0 + w, cy]
    if it.get("anticodon"):
        tb = TextBlock(it["anticodon"], "sans", 0.48 * u)
        ty = cy - 0.55 * u
        out.append(tb.svg(cx, ty - (tb.top + tb.bottom) / 2, style.INK))
    n = it.get("peptide", 1 if it.get("aa") else 0)
    r = 0.34 * u
    for k in range(n):                               # chain grows up (zig-zag) from the acceptor end
        bx = cx + 0.2 * u + (0.32 * u if k % 2 else 0)
        by = y0 - r * 0.9 - k * 1.55 * r
        hn = it.get("aa_hue") or AA_HUES[k % len(AA_HUES)]
        f, d = hue(hn)
        out.append(f'<circle cx="{bx:.2f}" cy="{by:.2f}" r="{r:.2f}" fill="{f}" stroke="{EDGE}" stroke-width="{0.6 * LW:.2f}"/>')
        bb = [min(bb[0], bx - r), min(bb[1], by - r), max(bb[2], bx + r), bb[3]]
    if it.get("aa") and n == 1:
        tb = TextBlock(it["aa"], "sans", 0.42 * u)
        out.append(tb.svg(cx + 0.2 * u, y0 - r * 0.9 - (tb.top + tb.bottom) / 2, style.INK))
    return "".join(out), tuple(bb), (cx, y0 + h / 2)


def g_ribosome(it):
    """Ribosome: large subunit above, small subunit below the mRNA line through `at` (= centre of the P site);
    `flip: true` puts the large subunit underneath (ER-bound ribosome, peptide going down). `w` width (S),
    `sites` draws the E / P / A slots, `large` / `small` subunit names, `hue` (default purple)."""
    cx, cy = _pt(it)
    w = P(it.get("w", 7.0))
    f = -1 if it.get("flip") else 1                  # +1: large subunit above the mRNA
    fill, dark = hue(it.get("hue", "purple"))
    light = _mix(fill, "#FFFFFF", 0.35)
    out = [f'<ellipse cx="{cx:.2f}" cy="{cy - f * 0.33 * w:.2f}" rx="{0.5 * w:.2f}" ry="{0.3 * w:.2f}" fill="{fill}" stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}"/>',
           f'<ellipse cx="{cx:.2f}" cy="{cy + f * 0.17 * w:.2f}" rx="{0.43 * w:.2f}" ry="{0.15 * w:.2f}" fill="{light}" stroke="{EDGE}" stroke-width="{0.8 * LW:.2f}"/>']
    slot = 0.27 * w
    if it.get("sites"):
        for k, name in zip((-1, 0, 1), "EPA"):
            sx = cx + k * slot
            y0 = cy - 0.47 * w if f == 1 else cy + 0.02 * w
            out.append(f'<rect x="{sx - 0.11 * w:.2f}" y="{y0:.2f}" width="{0.22 * w:.2f}" height="{0.45 * w:.2f}" '
                       f'rx="{0.04 * w:.2f}" fill="#FFFFFF" fill-opacity="0.45"/>')
            tb = TextBlock(name, "sans", 0.6 * S)          # site letters sit in the small subunit, under each slot
            out.append(tb.svg(sx, cy + f * 0.095 * w - (tb.top + tb.bottom) / 2, dark))
    for key, y in (("large", cy - f * 0.555 * w), ("small", cy + f * 0.245 * w)):
        if it.get(key):
            tb = TextBlock(it[key], "sans", 0.62 * S)
            out.append(tb.svg(cx, y - (tb.top + tb.bottom) / 2, style.INK))
    top, bot = (cy - 0.63 * w, cy + 0.32 * w) if f == 1 else (cy - 0.32 * w, cy + 0.63 * w)
    return "".join(out), (cx - 0.5 * w, top, cx + 0.5 * w, bot), (cx, cy - f * 0.33 * w)


GLYPHS["ribosome"] = g_ribosome
