"""A molecule as a Graphic in any of its views, optionally with its name underneath."""
from .. import style
from ..draw.chain import draw_chain, draw_formula
from ..draw.free import draw_free
from ..draw.lipid import draw_lipid
from ..draw.haworth import draw_chair, draw_haworth, draw_ring_molecule
from ..draw.skeletal import draw_skeletal
from ..graphic import Graphic
from ..text import TextBlock

S = style.S
NAME_GAP = 1.45 * S       # structure bottom → first name baseline


def views_of(m):
    v = []
    if m.get("lipid"):
        v.append("lipid")
    if m.get("haworth"):
        v.append("haworth")
    if m.get("chair"):
        v.append("chair")
    if m.get("ring"):
        v.append("ring")
        if isinstance(m["ring"], dict) and m["ring"].get("mirror"):
            v.append("ring_mirror")
    if m.get("free"):
        v.append("free")
    if m.get("chain"):
        v.append("chain")
    if m.get("formula"):
        v.append("formula")
    v.append("skeletal")
    return v


def structure(m, view=None, deco=None, chain=None):
    """chain: name of an alternative layout in m["chains"] (None = the default `chain`)."""
    view = view or m["default_view"]
    if deco == "default":
        deco = m.get("decorate")
    if view == "lipid":
        g = draw_lipid(m["lipid"], deco)
    elif view == "haworth":
        g = draw_haworth(m["haworth"], deco)
    elif view == "chair":
        g = draw_chair(m["chair"], deco)
    elif view in ("ring", "ring_mirror"):
        g = draw_ring_molecule(m, deco, mirror=view == "ring_mirror")
    elif view == "free":
        g = draw_free(m["free"], deco)
    elif view == "chain":
        g = draw_chain(m["chains"][chain] if chain else m["chain"], deco)
    elif view == "formula":
        g = draw_formula(m["formula"], deco)
    elif view == "skeletal":
        g = draw_skeletal(m["smiles"], deco)
    elif view == "name_line":                         # one-line name (text versions of reactions)
        tb = TextBlock(m["name"], "sans", style.NAME_SIZE)
        g = Graphic(tb.svg(0, 0), -tb.width / 2, tb.top, tb.width / 2, tb.bottom, {"base": (0, tb.bottom), "axis": (0, -tb.lines[0].cap / 2)})
    elif view == "name":
        tb = TextBlock(m.get("label") or m["name"], "sans", style.NAME_SIZE)
        g = Graphic(tb.svg(0, 0), -tb.width / 2, tb.top, tb.width / 2, tb.bottom, {"base": (0, tb.bottom), "axis": (0, -tb.lines[0].cap / 2)})
    else:
        raise ValueError(f"unknown view {view}")
    return g


def name_block(m):
    return TextBlock(m.get("label") or m["name"], "sans", style.NAME_SIZE, leading=1.15)


def with_name(g, m, name_y=None):
    """Append the molecule name centred under the structure. name_y = first baseline."""
    tb = name_block(m)
    y = name_y if name_y is not None else g.anchors["base"][1] + NAME_GAP
    body = g.body + tb.svg(g.cx, y)
    return Graphic(body, min(g.x0, g.cx - tb.width / 2), g.y0, max(g.x1, g.cx + tb.width / 2),
                   y + tb.bottom, dict(g.anchors))


def show_name(m, view):
    return view not in ("formula", "name", "name_line") and m.get("show_name", True)
