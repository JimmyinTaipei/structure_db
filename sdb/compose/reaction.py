"""Single-step reaction figure (Lehninger layout).

    substrates  ──(curved cofactor arrows / prosthetic text)──▶  products
                        enzyme (cyan) / conditions (black)

Coordinates: every structure is bottom-aligned at y=0; the arrow sits at arrow_y < 0.
"""
from .. import style
from ..graphic import BBox, Graphic, line, path
from ..text import TextBlock, TextLine
from ..draw.chain import Paint, draw_formula
from .molecule import NAME_GAP, name_block, show_name, structure

S = style.S
PLUS_GAP = 0.75 * S
SIDE_GAP = 0.9 * S
HEAD = 0.55 * S
TIER = [2.15 * S, 4.05 * S, 5.9 * S]       # cofactor label baselines above the arrow
PLUS_R = 0.3 * S                           # half-size of the thin "+" between molecules


def _plus(cx, cy):
    return line(cx - PLUS_R, cy, cx + PLUS_R, cy) + line(cx, cy - PLUS_R, cx, cy + PLUS_R)


def _head(x, y, dx, dy):
    """Open chevron arrowhead with its tip at (x, y) pointing along (dx, dy)."""
    import math
    n = math.hypot(dx, dy)
    ux, uy = dx / n, dy / n
    a = math.radians(24)
    pts = []
    for s in (1, -1):
        c, si = math.cos(a), s * math.sin(a)
        bx, by = -(ux * c - uy * si), -(ux * si + uy * c)
        pts.append((x + bx * HEAD, y + by * HEAD))
    return path(f"M{pts[0][0]:.2f},{pts[0][1]:.2f} L{x:.2f},{y:.2f} L{pts[1][0]:.2f},{pts[1][1]:.2f}")


def _cof_label(db, name):
    if name in db.cofactors:
        return db.cofactors[name].get("label", name)
    return db.molecules[name].get("formula") or db.molecules[name]["name"]


def _cof_svg(t, x, y, deco):
    """A cofactor label, coloured per the reaction's decorate (e.g. the H of NADH in cyan)."""
    fill, spans = Paint((deco or {}).get("colors")).text(0, "c", t)
    return t.svg(x, y, fill, spans=spans)


def arrow_block(db, r, lay, deco=None):
    deco = deco or {}
    """Arrow + cofactors + enzyme text, arrow at y=0, starting at x=0."""
    parts, bb = [], BBox()
    x = 0.7 * S
    joins = []
    labels_right = 0.0                        # right edge of every label placed so far
    y0 = -0.13 * S if r["reversible"] else 0.0    # curves meet the upper harpoon of ⇌
    n_in = len(r["cofactors_in"])
    prev_right = None
    in_boxes = []
    cof, cof_bb, sub_extra = [], BBox(), {}          # finer sticker layers: cofactors / prosthetic
    for i, c in enumerate(r["cofactors_in"]):
        t = TextLine(_cof_label(db, c))
        tier = TIER[min(n_in - 1 - i, len(TIER) - 1)]
        cx = x + t.width / 2 if prev_right is None else max(prev_right + 0.45 * S + t.width / 2, x + t.width / 2)
        by = -tier
        parts.append(f'<g data-ref="cofactor:{c}">{_cof_svg(t, cx - t.width / 2, by, deco.get(c))}</g>')
        bb.add(cx - t.width / 2, by + t.top, cx + t.width / 2, by + t.bottom)
        cof.append(parts[-1])
        cof_bb.add(cx - t.width / 2, by + t.top, cx + t.width / 2, by + t.bottom)
        in_boxes.append((cx - t.width / 2, by + t.top, cx + t.width / 2, by + t.bottom))
        sy = by + 0.45 * S
        drop = y0 - sy
        jx = cx + 0.8 * drop
        parts.append(path(f"M{cx:.2f},{sy:.2f} Q{cx:.2f},{y0:.2f} {jx:.2f},{y0:.2f}"))
        cof.append(parts[-1])
        cof_bb.add(cx, sy, jx, y0)
        joins.append(jx)
        prev_right = cx + t.width / 2
        labels_right = max(labels_right, prev_right)
        x = cx - t.width / 2
    x = max(joins) + 0.5 * S if joins else x

    if r["prosthetic"]:
        kind = r.get("prosthetic_font", "serif")
        size = S if kind == "serif" else style.ENZYME_SIZE
        tb = TextBlock(r["prosthetic"], kind, size, leading=1.12)
        last = -0.5 * S - (0.0 if r["reversible"] else 0)
        top_y = last - (len(tb.lines) - 1) * tb.lh
        if r["reversible"]:
            top_y -= 0.15 * S
        cx = x + tb.width / 2 + (0.2 * S if joins else 0.6 * S)
        # a tall prosthetic block (FMN, Fe-S) must also clear the labels of the cofactors coming in
        for bx0, by0, bx1, by1 in in_boxes:
            if top_y + tb.top < by1 + 0.15 * S and by0 < top_y + tb.bottom and cx - tb.width / 2 < bx1 + 0.6 * S:
                cx = bx1 + 0.6 * S + tb.width / 2
        parts.append(tb.svg(cx, top_y))
        bb.add(cx - tb.width / 2, top_y + tb.top, cx + tb.width / 2, top_y + tb.bottom)
        sub_extra["prosthetic"] = Graphic(parts[-1], cx - tb.width / 2, top_y + tb.top, cx + tb.width / 2, top_y + tb.bottom)
        x = cx + tb.width / 2 + 0.5 * S
        labels_right = max(labels_right, cx + tb.width / 2)

    for j, c in enumerate(r["cofactors_out"]):
        t = TextLine(_cof_label(db, c))
        tier = TIER[min(j, len(TIER) - 1)]
        by = -tier
        ey = by + 0.45 * S
        drop = y0 - ey
        # the label must clear every label already placed (e.g. NAD⁺ next to NADH + H⁺)
        ex = max(x + 0.9 * drop, labels_right + 0.45 * S + t.width / 2)
        sx = ex - 0.9 * drop
        parts.append(path(f"M{sx:.2f},{y0:.2f} Q{ex:.2f},{y0:.2f} {ex:.2f},{ey:.2f}"))
        parts.append(_head(ex, ey, 0, -1))
        parts.append(f'<g data-ref="cofactor:{c}">{_cof_svg(t, ex - t.width / 2, by - 0.25 * S, deco.get(c))}</g>')
        bb.add(ex - t.width / 2, by - 0.25 * S + t.top, ex + t.width / 2, by)
        cof += parts[-3:]
        cof_bb.add(sx, ey - 0.3 * S, ex + t.width / 2, y0)
        cof_bb.add(ex - t.width / 2, by - 0.25 * S + t.top, ex + t.width / 2, by)
        x = max(ex, ex + t.width / 2 - 0.4 * S) + 0.3 * S
        labels_right = max(labels_right, ex + t.width / 2)

    below = []
    if r.get("enzyme"):
        below.append((TextBlock(r["enzyme"], "sans", style.ENZYME_SIZE, leading=1.15), "cyan"))
    if r.get("conditions"):
        below.append((TextBlock(r["conditions"], "sans", style.ENZYME_SIZE, leading=1.15), "black"))
    below_w = max((tb.width for tb, _ in below), default=0)

    length = max(x + 0.9 * S, below_w + 1.2 * S, lay.get("arrow_min", 4.2) * S)
    length = max(length, 3.2 * S)
    content_r = x
    shift = max(0.0, (length - 0.9 * S - content_r) / 2) if not (joins or r["cofactors_out"]) else 0
    if shift:
        parts = [f'<g transform="translate({shift:.2f},0)">{"".join(parts)}</g>']
        bb.add(bb.x0 + shift, bb.y0, bb.x1 + shift, bb.y1)

    n0 = len(parts)
    if r["reversible"]:
        g = 0.13 * S
        parts.append(line(0, -g, length, -g))
        parts.append(path(f"M{length - HEAD:.2f},{-g - 0.32 * S:.2f} L{length:.2f},{-g:.2f}"))
        parts.append(line(0, g, length, g))
        parts.append(path(f"M{HEAD:.2f},{g + 0.32 * S:.2f} L0,{g:.2f}"))
        bb.add(0, -g - 0.4 * S, length, g + 0.4 * S)
        line_g = Graphic("".join(parts[n0:]), 0, -g - 0.4 * S, length, g + 0.4 * S)
    else:
        parts.append(line(0, 0, length, 0))
        parts.append(_head(length, 0, 1, 0))
        bb.add(0, -0.3 * S, length, 0.3 * S)
        line_g = Graphic("".join(parts[n0:]), 0, -0.3 * S, length, 0.3 * S)

    arrow_only = Graphic("".join(parts), bb.x0, bb.y0, bb.x1, bb.y1)
    # "arrow" = the whole arrow block (as before); finer layers below
    sub = {"arrow": arrow_only}
    if cof or sub_extra:
        sub["arrow_line"] = line_g
        if cof:
            body = "".join(cof)
            if shift:
                body = f'<g transform="translate({shift:.2f},0)">{body}</g>'
            sub["cofactors"] = Graphic(body, cof_bb.x0 + shift, cof_bb.y0, cof_bb.x1 + shift, cof_bb.y1)
        if "prosthetic" in sub_extra:
            pg = sub_extra["prosthetic"]
            sub["prosthetic"] = pg.moved(shift, 0) if shift else pg
    y = 1.3 * S
    for tb, col in below:
        svg = tb.svg(length / 2, y, col)
        parts.append(svg)
        x0, y0, x1, y1 = length / 2 - tb.width / 2, y + tb.top, length / 2 + tb.width / 2, y + tb.bottom
        bb.add(x0, y0, x1, y1)
        sub["enzyme" if col == "cyan" else "conditions"] = Graphic(svg, x0, y0, x1, y1)
        y += tb.bottom + 1.05 * S
    return Graphic(f'<g data-ref="reaction:{r["id"]}" class="sdb-arrow">{"".join(parts)}</g>',
                   bb.x0, bb.y0, bb.x1, bb.y1, {}, sub)


def _side(db, r, items, layout, deco, lay, arrow_y):
    """Lay out one side. Returns (Graphic, list of (mol, structure Graphic placed))."""
    placed = []
    gs = []
    for p in items:
        if p.get("cofactor"):
            # a cofactor written into the equation (+ H2O, + NAD+): plain formula text, no name
            c = p["cofactor"]
            m = {"id": c, "name": c}
            view = "formula"
            g = draw_formula(_cof_label(db, c), deco.get(c))
        else:
            m = db.molecules[p["id"]]
            view = p.get("view") or m["default_view"]
            g = structure(m, view, deco.get(p["id"]), p.get("chain"))
        if p["n"] > 1:
            t = TextLine(str(p["n"]))
            if view in ("formula",) + NAME_VIEWS:    # text partners end up centred on the arrow line
                cy = (g.y0 + g.y1) / 2
            else:
                cy = arrow_y + g.anchors["base"][1]  # g is still in its own coordinates (base not yet at 0)
            body = t.svg(g.x0 - 0.35 * S - t.width, cy + t.cap / 2) + g.body
            g = Graphic(body, g.x0 - 0.35 * S - t.width, g.y0, g.x1, g.y1, g.anchors, g.parts)
        nd = lay.get("nudge", {}).get(p["id"])
        if nd:
            g = g.moved(nd[0] * S, nd[1] * S)
        gs.append((m, view, g))

    parts, bb = [], BBox()
    if layout == "stack":
        # first item on top, last item bottom-aligned at y=0
        y = 0.0
        stacked = []
        for k, (m, view, g) in enumerate(reversed(gs)):
            dy = y - g.anchors["base"][1] if k == 0 else y - g.y1
            stacked.append((m, view, g.moved(-g.cx, dy)))
            top = stacked[-1][2].y0
            if k < len(gs) - 1:
                py = top - 0.5 * S - PLUS_R
                parts.append(_plus(0, py))
                bb.add(-PLUS_R, py - PLUS_R, PLUS_R, py + PLUS_R)
                y = py - PLUS_R - 0.5 * S
        stacked.reverse()
        for m, view, g in stacked:
            parts.append(f'<g data-ref="molecule:{m["id"]}" class="sdb-mol">{g.body}</g>')
            bb.add_g(g)
            placed.append((m, view, g))
    else:
        x = 0.0
        for k, (m, view, g) in enumerate(gs):
            if k:
                parts.append(_plus(x + PLUS_GAP + PLUS_R, arrow_y))
                x += PLUS_GAP * 2 + 2 * PLUS_R
            from ..draw.chain import PITCH
            # the name under a narrow structure can be wider than the structure itself
            # (Dihydroxyacetone phosphate): keep neighbouring names apart as well
            left, right = g.x0, g.x1
            if show_name(m, view) and len(gs) > 1:
                half = name_block(m).width / 2 + 0.25 * S
                left, right = min(left, g.cx - half), max(right, g.cx + half)
            if view in NAME_VIEWS or (view == "formula" and not show_name(m, view) and len(gs) > 1) or \
                    (g.h < 2.2 * PITCH and max(o.h for _, _, o in gs) > 2.2 * PITCH):
                # small partners (+ H2O, + NAD+, one-row carbamoyl phosphate / urea) sit on the
                # arrow line, like the "+" signs, instead of dropping to the bottom
                g2 = g.moved(x - left, arrow_y - (g.y0 + g.y1) / 2)
            else:
                g2 = g.moved(x - left, -g.anchors["base"][1])
            parts.append(f'<g data-ref="molecule:{m["id"]}" class="sdb-mol">{g2.body}</g>')
            bb.add_g(g2)
            placed.append((m, view, g2))
            x = g2.x1 + (right - g.x1)
    return Graphic("".join(parts), bb.x0, bb.y0, bb.x1, bb.y1), placed


def merge_deco(base, extra):
    """Overlay a variant's decoration on the normal one, molecule by molecule and row by row."""
    out = {k: (dict(v) if isinstance(v, dict) else v) for k, v in (base or {}).items()}
    for mol, d in (extra or {}).items():
        if not isinstance(out.get(mol, {}), dict):
            out[mol] = {}
        tgt = out.setdefault(mol, {})
        for key, val in d.items():
            if isinstance(val, dict) and isinstance(tgt.get(key), dict):
                tgt[key] = {**tgt[key], **val}
            else:
                tgt[key] = val
    return out


TEXT_VARIANTS = ("text", "text_full")
NAME_VIEWS = ("name", "name_line")


def _text_view(db, p):
    """Small molecules keep their formula (CO₂, NH₄⁺), everything else its one-line name."""
    return "formula" if db.molecules[p["id"]].get("formula") else "name_line"


def text_reaction(db, r, variant):
    """Reaction as seen by a text variant: every molecule written as its name, no decoration;
    `text` also drops the enzyme, conditions and cofactors (just the substances and the arrow)."""
    q = {**r, "decorate": {}, "layout": {k: v for k, v in (r.get("layout") or {}).items() if k == "tiers"},
         "substrates": [{**p, "view": _text_view(db, p), "chain": None} for p in r["substrates"]],
         "products": [{**p, "view": _text_view(db, p), "chain": None} for p in r["products"]],
         "substrates_layout": "row", "products_layout": "row"}
    if variant == "text":
        q.update(cofactors_in=[], cofactors_out=[], prosthetic=[], enzyme="", conditions="",
                 cofactors_inline=False)
    return q


def compose_reaction(db, r, variant=None):
    if variant in TEXT_VARIANTS:
        r, variant = text_reaction(db, r, variant), None
    lay = r.get("layout", {})
    deco = r.get("decorate", {})
    if variant:
        deco = merge_deco(deco, r["variants"][variant].get("decorate"))
    # arrow height derives from the tallest substrate structure
    heights = []
    for p in r["substrates"]:
        m = db.molecules[p["id"]]
        g = structure(m, p.get("view") or m["default_view"], deco.get(p["id"]), p.get("chain"))
        heights.append(g.anchors["base"][1] - g.y0)
    arrow_y = -max(heights) * lay.get("arrow_at", 0.36)
    if all(p.get("view") in NAME_VIEWS + ("formula",) for p in r["substrates"]) and \
            any(p.get("view") in NAME_VIEWS for p in r["substrates"] + r["products"]):
        arrow_y = 0.0                                    # text version: names sit on the arrow line

    subs, prods, arrow_r = r["substrates"], r["products"], r
    if r.get("cofactors_inline"):
        subs = subs + [{"cofactor": c, "id": c, "n": 1} for c in r["cofactors_in"]]
        prods = prods + [{"cofactor": c, "id": c, "n": 1} for c in r["cofactors_out"]]
        arrow_r = {**r, "cofactors_in": [], "cofactors_out": []}
    left, lp = _side(db, r, subs, r.get("substrates_layout", "row"), deco, lay, arrow_y)
    right, rp = _side(db, r, prods, r.get("products_layout", "row"), deco, lay, arrow_y)
    arrow = arrow_block(db, arrow_r, lay, deco)

    off_l = -left.x0
    left = left.moved(off_l, 0)
    arrow = arrow.moved(left.x1 + SIDE_GAP, arrow_y)
    off_r = arrow.x1 + SIDE_GAP - right.x0
    right = right.moved(off_r, 0)

    bb = BBox()
    for g in (left, arrow, right):
        bb.add_g(g)
    parts = [left.body, arrow.body, right.body]
    sub = dict(arrow.parts)
    for placed, off in ((lp, off_l), (rp, off_r)):
        for m, view, g in placed:
            g = g.moved(off, 0)
            plain_g = Graphic(g.body, g.x0, g.y0, g.x1, g.y1)
            sub[f"mol_{m['id']}"] = plain_g
            if "highlight" in g.parts:
                sub[f"highlight_{m['id']}"] = g.parts["highlight"]

    # names on one common baseline below the lowest structure
    base_y = max(g.y1 for _, _, g in lp + rp)
    name_y = base_y + NAME_GAP + lay.get("name_dy", 0) * S
    for placed, off in ((lp, off_l), (rp, off_r)):
        for m, view, g in placed:
            if not show_name(m, view):
                continue
            tb = name_block(m)
            cx = g.cx + off
            svg = tb.svg(cx, name_y)
            parts.append(f'<g data-ref="molecule:{m["id"]}" class="sdb-name">{svg}</g>')
            nb = (cx - tb.width / 2, name_y + tb.top, cx + tb.width / 2, name_y + tb.bottom)
            bb.add(*nb)
            sub[f"name_{m['id']}"] = Graphic(svg, *nb)
    return Graphic("".join(parts), bb.x0, bb.y0, bb.x1, bb.y1, {}, sub)
