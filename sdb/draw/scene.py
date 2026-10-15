"""`scene` figures: illustrated glyphs (membranes, proteins, DNA …) placed by coordinates, with labels and
leader lines. Used as `layout: scene` or as a pathway variant `variants: {illus: {scene: …}}`.

    scene:
      size: [32, 18]                       # S units (informational; the bbox grows to fit)
      items:
        - {g: membrane, y: 6, x0: 0, x1: 32}
        - {g: gpcr, at: [8, 6], label: "β-adrenergic\\nreceptor", label_at: [8, 1.5]}
"""
from .. import style
from ..graphic import BBox, Graphic
from ..text import TextBlock
import math

from .glyphs import GLYPHS, _sample, _trna_pts, label_bg

S = style.S
LABEL_SIZE = style.NAME_SIZE * 0.9
PROTEINS = {"blob", "gpcr", "channel", "transporter", "pump", "gprotein", "ring", "ribosome"}
DEFAULT_HUE = {"ribosome": "purple", "gpcr": "purple", "channel": "teal", "pump": "orange", "gprotein": "green", "ring": "green", "ion": "yellow"}
PAD = 0.28 * S
# glyphs whose bbox a label must not cover (membranes / regions are backgrounds; arrows and text are thin)
SOLID = PROTEINS | {"ion", "ligand", "strand", "dsdna", "trna"}


def _label(it, bb, anchor):
    """Label (Alberts: tinted box for proteins, plain text otherwise) + thin black leader line.
    Returns (svg, label bbox)."""
    tb = TextBlock(it["label"], "sans", it.get("label_size", 0.9) * style.NAME_SIZE, align="center")
    boxed = it.get("label_box", it.get("g") in PROTEINS)
    pad = PAD if boxed else 0
    h = tb.bottom - tb.top + 2 * pad
    if it.get("label_at"):
        lx, ly = it["label_at"][0] * S, it["label_at"][1] * S
        leader = it.get("leader", True)
    else:                                          # default: centred under the glyph, no leader
        lx, ly = (bb[0] + bb[2]) / 2, bb[3] + 0.45 * S + h / 2
        leader = False
    dy = -(tb.top + tb.bottom) / 2
    if it.get("label_align") == "left":            # label_at = left edge (right-hand label columns)
        lx += tb.width / 2 + pad
    box = (lx - tb.width / 2 - pad, ly + dy + tb.top - pad, lx + tb.width / 2 + pad, ly + dy + tb.bottom + pad)
    col = style.color(it.get("label_color", "black"))
    svg = ""
    if boxed:
        svg += (f'<rect x="{box[0]:.2f}" y="{box[1]:.2f}" width="{box[2] - box[0]:.2f}" height="{box[3] - box[1]:.2f}" '
                f'fill="{label_bg(it.get("hue", DEFAULT_HUE.get(it.get("g"), "blue")))}"/>')
    svg += tb.svg(lx, ly + dy, col)
    if leader:
        ax, ay = (it["point"][0] * S, it["point"][1] * S) if it.get("point") else anchor
        if box[0] < ax < box[2]:                   # anchor above / below: leave from the top or bottom edge
            sx, sy = ax, (box[1] if ay < box[1] else box[3])
        else:
            sx, sy = (box[0] if ax < box[0] else box[2]), min(max(ay, box[1]), box[3])
        svg = (f'<line x1="{sx:.2f}" y1="{sy:.2f}" x2="{ax:.2f}" y2="{ay:.2f}" stroke="{style.INK}" '
               f'stroke-width="{0.6 * style.LINE_W:.2f}"/>') + svg
    return svg, box


def render_items(spec):
    """[(name, svg, glyph bbox, label bbox | None)] in drawing order."""
    out, seen = [], {}
    for i, it in enumerate(spec.get("items") or []):
        g = it.get("g")
        if g not in GLYPHS:
            raise ValueError(f"未知的插畫元件 g: {g}")
        svg, bb, anchor = GLYPHS[g](it)
        lb = None
        if it.get("label"):
            lsvg, lb = _label(it, bb, anchor)
            svg += lsvg
        name = it.get("id") or f"{g}{seen.get(g, 0) + 1}"
        seen[g] = seen.get(g, 0) + 1
        out.append((name, svg, bb, lb, g))
    return out


def compose_scene(spec):
    items = render_items(spec)
    bb = BBox()
    parts = {}
    body = []
    for name, svg, gb, lb, _ in items:
        x0, y0, x1, y1 = gb
        if lb:
            x0, y0, x1, y1 = min(x0, lb[0]), min(y0, lb[1]), max(x1, lb[2]), max(y1, lb[3])
        bb.add(x0, y0, x1, y1)
        body.append(svg)
        parts[name] = Graphic(svg, x0, y0, x1, y1)
    if spec.get("caption"):
        tb = TextBlock(spec["caption"], "sans", style.NAME_SIZE)
        cx, cy = (bb.x0 + bb.x1) / 2, bb.y1 + 0.9 * S - tb.top
        body.append(tb.svg(cx, cy))
        bb.add(cx - tb.width / 2, cy + tb.top, cx + tb.width / 2, cy + tb.bottom)
    return Graphic("".join(body), bb.x0, bb.y0, bb.x1, bb.y1, parts=parts)


def _overlap(a, b, pad=0.05 * S):
    return a[0] < b[2] - pad and b[0] < a[2] - pad and a[1] < b[3] - pad and b[1] < a[3] - pad


def _hits(it):
    """Sampled (x, y, r) discs along strands / duplexes (their bbox is mostly empty space)."""
    if it["g"] == "strand":
        pts = [(x * S, y * S) for x, y in it["points"]]
        r = (it.get("width", 0.26) / 2 + (it.get("base_len", 0.6) if it.get("bases") else 0)) * S
        return [(x, y, r) for x, y, _, _ in _sample(pts, 0.25 * S)]
    if it["g"] == "dsdna":
        (x0, y0), (x1, y1) = [(a * S, b * S) for a, b in (it["from"], it["to"])]
        r = it.get("width", 1.3) / 2 * S + 0.13 * S
        return [(x, y, r) for x, y, _, _ in _sample([(x0, y0), (x1, y1)], 0.25 * S)]
    if it["g"] == "trna" and it.get("shape") != "icon":   # backbone only (loops are open space)
        pts, _ = _trna_pts(it)
        return [(x, y, 0.15 * S) for x, y, _, _ in _sample(pts, 0.2 * S)]
    return None


def _box_hits(box, discs):
    return any(box[0] - r < x < box[2] + r and box[1] - r < y < box[3] + r and
               math.hypot(max(box[0] - x, 0, x - box[2]), max(box[1] - y, 0, y - box[3])) < r for x, y, r in discs)


def check_scene(spec):
    """(errors, warnings): unknown glyphs / bad parameters are errors; overlapping labels are warnings."""
    try:
        items = render_items(spec)
    except Exception as e:
        return [f"插畫設定錯誤 {e!r}"], []
    warns = []
    labels = [(n, lb) for n, _, _, lb, _ in items if lb]
    labels += [(n, gb) for n, _, gb, _, g in items if g == "text"]
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            if _overlap(labels[i][1], labels[j][1]):
                warns.append(f"插畫標籤重疊：{labels[i][0]} ↔ {labels[j][0]}")
    specs = spec.get("items") or []
    solids = [(n, gb, _hits(it)) for (n, _, gb, _, g), it in zip(items, specs) if g in SOLID]
    for ln, lb in labels:
        for gn, gb, discs in solids:
            if gn == ln:
                continue
            if (_box_hits(lb, discs) if discs is not None else _overlap(lb, gb, pad=0.1 * S)):
                warns.append(f"插畫文字蓋到圖形：{ln} → {gn}")
    return [], warns
