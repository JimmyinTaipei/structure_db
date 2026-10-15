"""Free-coordinate drawings for the few molecules a vertical chain cannot express (proline).

    free:
      atoms: {a: ["COO^-", 0, 0], ca: [C, 0, 1], n: ["^+H_2N", -1.05, 1], ...}   # x, y in chain pitches
                                                       # optional 4th item: index of the anchor char
      bonds: [[a, ca], [ca, n], [n, cd, 1], ...]       # optional 3rd item: 1, 2, "bold" or "hb" (dashed H-bond)
      notes: [["Glycine", 2.4, -1.0, magenta], …]      # optional: sans labels (not atoms)
      schematic: true                                   # ring skeleton without H: skip the formula check
    decoration:
      colors: {atom_id: colour | [colour, "fragment"]}
      boxes:  [{atoms: [cb, cg], color: salmon, name: R}, …]   # one rectangle per entry; entries
                                                                # sharing a name form one layer
"""
import math
import re
from collections import Counter

from .. import style
from ..graphic import BBox, Graphic, line
from ..text import TextBlock, TextLine, anchor_index, plain
from .chain import BOX_PAD, DBL, GAP_H, GAP_V, PITCH

S = style.S


def draw_free(spec, deco=None):
    deco = deco or {}
    colors = deco.get("colors") or {}
    parts, back = [], []
    bb = BBox()
    pos, label_size, text_boxes = {}, {}, {}
    for aid, (text, gx, gy, *anchor) in spec["atoms"].items():
        x, y = gx * PITCH, gy * PITCH
        pos[aid] = (x, y)
        if not text:
            bb.add(x, y, x, y)                         # bare skeleton vertex still counts
            continue
        if text.startswith("@"):                       # SNFG monosaccharide symbol (glycan cartoons)
            svg, r = snfg(text[1:], x, y)
            parts.append(svg)
            label_size[aid] = (r - GAP_H, r - GAP_V)
            bb.add(x - r, y - r, x + r, y + r)
            continue
        t = TextLine(text)
        ai = anchor[0] if anchor else anchor_index(text)   # 4th item: anchor char ("CoA" left of S → 2)
        x0, x1 = t.char_box(ai)
        p = plain(text)
        if ai + 1 < len(p) and p[ai + 1].islower():     # two-letter element (Co, Fe, Mg): centre on both
            x1 = t.char_box(ai + 1)[1]
        tx, base = x - (x0 + x1) / 2, y + t.cap / 2
        c = colors.get(aid)
        if isinstance(c, list):
            sp = t.span(c[1])
            parts.append(t.svg(tx, base, None, spans=[(sp[0], sp[1], c[0])] if sp else None))
        else:
            parts.append(t.svg(tx, base, c))
        label_size[aid] = ((x1 - x0) / 2, t.cap / 2)
        text_boxes[aid] = (tx, base - t.cap, tx + t.char_box(len(t.chars) - 1)[1], base + t.size * 0.12)
        bb.add(tx, base + t.top, tx + t.width, base + t.bottom)

    def trim(aid, ux, uy):
        if aid not in label_size:
            return 0.0
        hw, hh = label_size[aid]
        return min((hw + GAP_H) / abs(ux) if ux else 1e9, (hh + GAP_V) / abs(uy) if uy else 1e9)

    for bond in spec.get("bonds", []):
        a, b = bond[0], bond[1]
        order = bond[2] if len(bond) > 2 else 1
        (x1, y1), (x2, y2) = pos[a], pos[b]
        L = math.hypot(x2 - x1, y2 - y1)
        ux, uy = (x2 - x1) / L, (y2 - y1) / L
        ta, tb = trim(a, ux, uy), trim(b, ux, uy)
        p1, p2 = (x1 + ux * ta, y1 + uy * ta), (x2 - ux * tb, y2 - uy * tb)
        if order == "hb":                                # hydrogen bond: grey dashes
            parts.append(line(*p1, *p2, stroke="gray").replace("/>", ' stroke-dasharray="3 3"/>', 1))
        elif order == "bold":                            # Haworth front edge
            parts.append(line(*p1, *p2, w=3.2 * style.LINE_W, cap="round"))
        elif order == 2:
            nx, ny = -uy * DBL, ux * DBL
            parts.append(line(p1[0] + nx, p1[1] + ny, p2[0] + nx, p2[1] + ny))
            parts.append(line(p1[0] - nx, p1[1] - ny, p2[0] - nx, p2[1] - ny))
        else:
            parts.append(line(*p1, *p2))

    layers = {}
    for spec_box in deco.get("boxes") or []:
        hb = BBox()
        for aid in spec_box.get("atoms", []):
            x0, y0, x1, y1 = text_boxes[aid]
            hb.add(x0 - BOX_PAD * 0.6, y0 - BOX_PAD * 0.6, x1 + BOX_PAD * 0.6, y1 + BOX_PAD * 0.6)
        if hb.x0 == float("inf"):
            continue
        rect = (f'<rect x="{hb.x0:.2f}" y="{hb.y0:.2f}" width="{hb.x1 - hb.x0:.2f}" height="{hb.y1 - hb.y0:.2f}" '
                f'fill="{style.color(spec_box.get("color", "salmon"))}"/>')
        back.append(rect)
        bb.add(hb.x0, hb.y0, hb.x1, hb.y1)
        key = f"box_{spec_box.get('name', 'R')}"
        prev = layers.get(key)
        if prev:
            layers[key] = Graphic(prev.body + rect, min(prev.x0, hb.x0), min(prev.y0, hb.y0),
                                  max(prev.x1, hb.x1), max(prev.y1, hb.y1))
        else:
            layers[key] = Graphic(rect, hb.x0, hb.y0, hb.x1, hb.y1)

    for text, gx, gy, *col in spec.get("notes") or []:
        tb = TextBlock(text, "sans", style.NAME_SIZE * 0.9)
        cx, cy = gx * PITCH, gy * PITCH
        base = cy - (tb.top + tb.bottom) / 2
        parts.append(tb.svg(cx, base, col[0] if col else None))
        bb.add(cx - tb.width / 2, base + tb.top, cx + tb.width / 2, base + tb.bottom)

    lowest = max(text_boxes.values(), key=lambda b: b[3])[3] if text_boxes else bb.y1
    if spec.get("notes"):                              # name goes under the notes too
        lowest = bb.y1
    return Graphic("".join(back + parts), bb.x0, bb.y0, bb.x1, bb.y1,
                   {"base": (0, lowest), "axis": (0, 0)}, layers)


# Symbol Nomenclature for Glycans: shape + colour per monosaccharide
SNFG = {"glc": ("circle", "#0072BC"), "gal": ("circle", "#FFD400"), "man": ("circle", "#00A651"),
        "glcnac": ("square", "#0072BC"), "galnac": ("square", "#FFD400"), "fru": ("pentagon", "#00A651"),
        "neu5ac": ("diamond", "#A54399"), "fuc": ("triangle", "#ED1C24")}


def snfg(name, x, y, r=0.42 * S):
    shape, fill = SNFG[name]
    st = f'fill="{fill}" stroke="{style.INK}" stroke-width="{style.LINE_W:.2f}"'
    if shape == "circle":
        return f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" {st}/>', r
    if shape == "square":
        return f'<rect x="{x - r:.2f}" y="{y - r:.2f}" width="{2 * r:.2f}" height="{2 * r:.2f}" {st}/>', r
    n, rot = {"pentagon": (5, -90), "diamond": (4, 0), "triangle": (3, -90)}[shape]
    pts = " ".join(f"{x + r * 1.1 * math.cos(math.radians(rot + 360 * i / n)):.2f},"
                   f"{y + r * 1.1 * math.sin(math.radians(rot + 360 * i / n)):.2f}" for i in range(n))
    return f'<polygon points="{pts}" {st}/>', r


def free_composition(spec):
    """Atom counts and net charge implied by the labels ('^+H_2N' → N + 2 H, charge +1)."""
    c, q = Counter(), 0
    for text, _x, _y, *_a in spec["atoms"].values():
        p = plain(text)
        q += p.count("+") - p.count("−") - p.count("-")
        for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", p):
            c[el] += int(n or 1)
    return c, q
