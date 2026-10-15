"""Plain text nodes for pathway / signalling diagrams."""
from .. import style
from ..graphic import Graphic
from ..text import TextBlock

S = style.S
KINDS = {
    #          font    text colour  box stroke   box fill
    "concept": ("sans", "black",    "#999999",   "none"),
    "enzyme":  ("sans", "cyan",     None,        "none"),
    "drug":    ("sans", "red",      "red",       "none"),
    "plain":   ("serif", "black",   None,        "none"),
    "heading": ("sans", "black",    None,        "none"),     # group titles in overview sheets
    "name":    ("sans", "black",    None,        "none"),     # intermediate without a molecule entry
    "protein": ("sans", "black",    "protein_line", "protein"),   # membrane complexes (I–V, transporters)
    "carrier": ("serif", "black",   "carrier_line", "carrier"),   # mobile carriers (Q, cyt c)
    "ion":     ("serif", "black",   None,        "none"),     # H⁺, e⁻ … (serif like cofactors)
}


def draw_textnode(text, kind="concept", min_w=0.0, min_h=0.0):
    """min_w / min_h (in S): a protein box tall enough to span the membrane, equal-width boxes."""
    font, fg, stroke, fill = KINDS.get(kind, KINDS["concept"])
    tb = TextBlock(text, font, style.NAME_SIZE if font == "sans" else S)
    pad_x, pad_y = 0.6 * S, 0.45 * S
    w = max(tb.width + 2 * pad_x, min_w * S)
    top, bottom = tb.top - pad_y, tb.bottom + pad_y
    if bottom - top < min_h * S:
        extra = (min_h * S - (bottom - top)) / 2
        top, bottom = top - extra, bottom + extra
    # vertically centre the block on y=0
    dy = -(top + bottom) / 2
    body = ""
    if stroke:
        body += (f'<rect x="{-w / 2:.2f}" y="{top + dy:.2f}" width="{w:.2f}" height="{bottom - top:.2f}" '
                 f'rx="{(0.9 if kind == "carrier" else 0.35) * S:.2f}" fill="{style.color(fill) if fill != "none" else "none"}" stroke="{style.color(stroke)}" stroke-width="{style.LINE_W:.2f}"/>')
    body += tb.svg(0, dy, fg)
    return Graphic(body, -w / 2, top + dy, w / 2, bottom + dy, {"base": (0, bottom + dy), "axis": (0, 0)})
