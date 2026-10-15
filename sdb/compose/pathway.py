"""Pathway / flow diagrams.

    linear  nodes in one row, left → right (follows the edges)
    free    every node has pos: [x, y] in "cells" (cell = size of the largest node)
    cycle   nodes without pos are spread clockwise on a circle (radius in cells);
            nodes with pos (e.g. acetyl-CoA entering) are placed explicitly

Edges: {from, to, step, label, side_text, reversible, bend, label_side, inhibit, color}
    color: orange (electron flow), magenta (proton flow) … recolours that arrow
Nodes: text nodes may set w / h (in S) as a minimum box size (protein complexes spanning a membrane)
Regions: membranes and compartment boxes drawn behind everything (see _regions)
    inhibit: true draws a red line ending in a bar (⊣) — drugs, feedback inhibition
    head: false draws a plain line without an arrowhead (family trees)
    resonance: true draws a straight double-headed arrow (↔) between resonance structures
    to: "edge:<from>-><to>"  ends on the midpoint of another edge (side entry)
    bend > 0 bows the arrow to the left of its travel direction (outwards on a clockwise ring)

Text versions (variants `text` / `text_full`, added automatically by the loader): every molecule
node is drawn as its name; `text` also drops enzyme labels and side texts and draws heavier
arrows (Lehninger overview style). `text_layout:` overrides layout keys (radius, cell, cell_size,
start_angle, bend …) and `pos: {node_id: [x, y]}` for the text versions only.

Graph layouts return parts: path (nodes + arrows), enzymes, cofactors — separable sticker layers.

Every node and arrow carries data-ref="molecule:<id>" / "reaction:<id>" so the web UI can
make them clickable.
"""
import math

from .. import style
from ..draw.textnode import draw_textnode
from ..graphic import BBox, Graphic, line, path
from ..text import TextBlock
from .molecule import name_block, show_name, structure
from .reaction import HEAD, _head

S = style.S
GAP = 0.8 * S
TEXT_VARIANTS = ("text", "text_full")
THICK = 2.0 * style.LINE_W


def _thick(svg):
    return svg.replace(f'stroke-width="{style.LINE_W:.2f}"', f'stroke-width="{THICK:.2f}"')


def text_view(p, variant):
    """Pathway dict as seen by a text variant: names for every molecule, text_layout applied."""
    tl = dict(p.get("text_layout") or {})
    pos = tl.pop("pos", {})
    nodes = []
    for n in p["nodes"]:
        n = dict(n)
        if n.get("ref"):
            n["view"] = "name"
            n.pop("chain", None)
        if n["id"] in pos:
            n["pos"] = pos[n["id"]]
        nodes.append(n)
    q = {**p, **tl, "nodes": nodes}
    if variant == "text":
        q["edges"] = [{k: v for k, v in e.items() if k not in ("label", "side_text")} for e in p.get("edges", [])]
    return q


def node_graphic(db, n, deco=None):
    if n.get("ref"):
        m = db.molecules[n["ref"]]
        view = n.get("view") or m["default_view"]
        g = structure(m, view, deco, n.get("chain"))
        # centre the structure vertically on y=0 and hang its name underneath
        g = g.moved(-g.cx, -(g.y0 + g.y1) / 2)
        if show_name(m, view):
            tb = name_block(m)
            y = g.y1 + 1.2 * S
            g = Graphic(g.body + tb.svg(0, y), min(g.x0, -tb.width / 2), g.y0,
                        max(g.x1, tb.width / 2), y + tb.bottom, g.anchors)
        return g, f"molecule:{m['id']}"
    return draw_textnode(n["text"], n.get("kind", "concept"), n.get("w", 0), n.get("h", 0)), None


def _order(p):
    nodes = {n["id"]: n for n in p["nodes"]}
    nxt = {e["from"]: e for e in p["edges"] if not e["to"].startswith("edge:")}
    targets = {e["to"] for e in nxt.values()}
    start = next(i for i in nodes if i not in targets)
    seq = [start]
    while seq[-1] in nxt:
        seq.append(nxt[seq[-1]]["to"])
    return [nodes[i] for i in seq], nxt


def compose_linear(db, p, variant=None):
    thick = variant == "text"
    if variant in TEXT_VARIANTS:
        p = text_view(p, variant)
    seq, nxt = _order(p)
    parts, bb = [], BBox()
    x = 0.0
    for i, n in enumerate(seq):
        g, ref = node_graphic(db, n)
        g = g.moved(x - g.x0, 0)
        attr = f' data-ref="{ref}" class="sdb-mol"' if ref else ""
        parts.append(f'<g{attr}>{g.body}</g>')
        bb.add_g(g)
        x = g.x1
        e = nxt.get(n["id"])
        if not e or i == len(seq) - 1:
            continue
        label = TextBlock(e.get("label", ""), "sans", style.ENZYME_SIZE) if e.get("label") else None
        length = max(4 * S, (label.width + 1.2 * S) if label else 0)
        x0 = x + GAP
        seg = line(x0, 0, x0 + length, 0) + _head(x0 + length, 0, 1, 0)
        if thick:
            seg = _thick(seg)
        if label:
            seg += label.svg(x0 + length / 2, -0.6 * S - (len(label.lines) - 1) * label.lh, "cyan")
            bb.add(x0, -0.6 * S + label.top - (len(label.lines) - 1) * label.lh, x0 + length, 0)
        ref = f' data-ref="reaction:{e["step"]}" class="sdb-step"' if e.get("step") else ""
        # a transparent hit area makes thin arrows easy to click in the web UI
        hit = (f'<rect x="{x0:.2f}" y="{-1.6 * S:.2f}" width="{length:.2f}" height="{2.2 * S:.2f}" '
               f'fill="transparent"/>') if ref else ""
        parts.append(f'<g{ref}>{hit}{seg}</g>')
        bb.add(x0, -0.4 * S, x0 + length, 0.4 * S)
        x = x0 + length + GAP
    return Graphic("".join(parts), bb.x0, bb.y0, bb.x1, bb.y1)


def _bez(p0, c, p1, t):
    return tuple((1 - t) ** 2 * a + 2 * (1 - t) * t * b + t * t * d for a, b, d in zip(p0, c, p1))


def _unit(dx, dy):
    n = math.hypot(dx, dy) or 1.0
    return dx / n, dy / n


def _exit(center, toward, box, gap):
    """Point where the ray center→toward leaves `box` (x0, y0, x1, y1) grown by `gap`."""
    ux, uy = _unit(toward[0] - center[0], toward[1] - center[1])
    hw, hh = (box[2] - box[0]) / 2 + gap, (box[3] - box[1]) / 2 + gap
    t = min(hw / abs(ux) if ux else 1e9, hh / abs(uy) if uy else 1e9)
    return center[0] + ux * t, center[1] + uy * t


def _half_head(tip, toward_dir):
    """Single barb on the traveller's left, like the harpoons of ⇌."""
    ux, uy = _unit(*toward_dir)
    lx, ly = uy, -ux
    bx, by = tip[0] - ux * HEAD + lx * 0.32 * S, tip[1] - uy * HEAD + ly * 0.32 * S
    return path(f"M{bx:.2f},{by:.2f} L{tip[0]:.2f},{tip[1]:.2f}")


def _hits(box, obstacles, points):
    """Number of obstacles / arrow sample points that fall inside box."""
    x0, y0, x1, y1 = box
    n = sum(1 for (a0, b0, a1, b1) in obstacles if x0 < a1 and a0 < x1 and y0 < b1 and b0 < y1)
    return n * 10 + sum(1 for px, py in points if x0 < px < x1 and y0 < py < y1)


def _place_block(tb, m, nv, tangent, obstacles, points, d=0.4 * S):
    """Centre for TextBlock tb beside point m, on the side of unit normal nv, avoiding obstacles.

    The block sits at distance d + (its support along nv) from m, so it never crosses the
    tangent line at m; if it collides with a node, another label or an arrow, it is pushed
    further out and/or slid along the arrow.
    """
    w, h = tb.width, tb.bottom - tb.top
    support = abs(nv[0]) * w / 2 + abs(nv[1]) * h / 2
    best = None
    for extra in [0, 0.3, 0.6, 1.0, 1.5, 2.1, 2.8, 3.6]:
        reach = w / S / 2 + 1.0              # sliding by half the block width clears a crossing arrow
        for slide in [0, 0.4, -0.4, 0.8, -0.8, 1.3, -1.3, 2.0, -2.0, 3.0, -3.0, reach, -reach]:
            t = d + support + extra * S
            cx = m[0] + nv[0] * t + tangent[0] * slide * S
            cy = m[1] + nv[1] * t + tangent[1] * slide * S
            box = (cx - w / 2 - 0.1 * S, cy - h / 2 - 0.1 * S, cx + w / 2 + 0.1 * S, cy + h / 2 + 0.1 * S)
            cost = _hits(box, obstacles, points)
            if cost == 0:
                return cx, cy, box
            if best is None or cost < best[0]:
                best = (cost, (cx, cy, box))
    return best[1]


def _draw_block(tb, cx, cy, fill, bb):
    base = cy - (tb.top + tb.bottom) / 2
    bb.add(cx - tb.width / 2, base + tb.top, cx + tb.width / 2, base + tb.bottom)
    return tb.svg(cx, base, fill)


def compose_graph(db, p, variant=None):
    from .reaction import merge_deco
    deco = {n["id"]: n["decorate"] for n in p["nodes"] if n.get("decorate")}
    if p.get("molecule_deco") == "default":           # every molecule node shows its own decoration
        deco = {n["id"]: deco.get(n["id"], "default") for n in p["nodes"] if n.get("ref")}
    text = variant in TEXT_VARIANTS
    if text:
        p, deco = text_view(p, variant), {}
    elif variant:
        deco = merge_deco(deco, p["variants"][variant].get("decorate"))
    nodes = {}
    for n in p["nodes"]:
        g, ref = node_graphic(db, n, deco.get(n["id"]))
        nodes[n["id"]] = (n, g.moved(-(g.x0 + g.x1) / 2, -(g.y0 + g.y1) / 2), ref)
    if p.get("cell_size"):                            # absolute cell size, in units of S
        cell = p["cell_size"] * S
    else:
        cell = max(max(g.w, g.h) for _, g, _ in nodes.values()) * p.get("cell", 1.0)

    pos = {}
    ring = [n for n in p["nodes"] if p["layout"] == "cycle" and "pos" not in n]
    if ring:
        radius = p.get("radius", max(1.6, len(ring) * 1.2 / (2 * math.pi)))
        a0 = p.get("start_angle", -90)
        for i, n in enumerate(ring):
            a = math.radians(a0 + 360 * i / len(ring))
            pos[n["id"]] = (radius * math.cos(a), radius * math.sin(a))
    for n in p["nodes"]:
        if "pos" in n:
            pos[n["id"]] = tuple(n["pos"])
    missing = [i for i in nodes if i not in pos]
    if missing:
        raise ValueError(f"節點沒有位置（free 版面需要 pos）：{missing}")

    parts, bb = [], BBox()
    layer = {k: ([], BBox()) for k in ("enzymes", "cofactors")}
    path_body = []
    boxes = {}
    region_svg, region_bb = _regions(p, cell)
    # finer sticker layers: what kind of node / arrow each piece is
    fine = {k: ([], BBox()) for k in ("regions", "molecules", "proteins", "labels", "drugs",
                                      "arrows", "electrons", "protons", "inhibitions")}
    if region_svg:
        fine["regions"][0].append(region_svg)
        fine["regions"][1].add(region_bb.x0, region_bb.y0, region_bb.x1, region_bb.y1)
    if region_svg:
        parts.append(f"<g>{region_svg}</g>")
        bb.add(region_bb.x0, region_bb.y0, region_bb.x1, region_bb.y1)
    for id_, (n, g, ref) in nodes.items():
        # node_anchor: top → pos is the top-centre of the node (so a row of amino acids lines
        # up on its COO⁻); otherwise pos is the node's centre. Headings are left-aligned.
        dx = g.w / 2 if n.get("kind") == "heading" else 0
        dy = g.h / 2 if p.get("node_anchor") == "top" else 0
        g = g.moved(pos[id_][0] * cell + dx, pos[id_][1] * cell + dy)
        boxes[id_] = (g.x0, g.y0, g.x1, g.y1)
        attr = f' data-ref="{ref}" class="sdb-mol"' if ref else ""
        parts.append(f"<g{attr}>{g.body}</g>")
        path_body.append(g.body)
        bb.add_g(g)
        kind = "molecules" if ref else NODE_LAYER.get(n.get("kind", "concept"), "labels")
        fine[kind][0].append(g.body)
        fine[kind][1].add_g(g)

    def center(i):
        b = boxes[i]
        return (b[0] + b[2]) / 2, (b[1] + b[3]) / 2

    geo = {}
    default_bend = p.get("bend", 0.12 if p["layout"] == "cycle" else 0.0)
    default_side = p.get("label_side", "right" if p["layout"] == "cycle" else "left")
    edges = sorted(p.get("edges", []), key=lambda e: e["to"].startswith("edge:") or e["from"].startswith("edge:"))
    drawn = []
    for e in edges:
        bend = e.get("bend", default_bend)
        if e["from"].startswith("edge:"):
            # branching off another arrow (e.g. fumarate leaving the ASL step)
            s0, c0, e0 = geo[e["from"][5:]]
            a = _bez(s0, c0, e0, 0.5)
            b = center(e["to"])
            ux, uy = _unit(b[0] - a[0], b[1] - a[1])
            start = (a[0] + ux * 0.25 * S, a[1] + uy * 0.25 * S)
            end = _exit(b, a, boxes[e["to"]], 0.4 * S)
        elif e["to"].startswith("edge:"):
            a = center(e["from"])
            s0, c0, e0 = geo[e["to"][5:]]
            target = _bez(s0, c0, e0, 0.5)
            ux, uy = _unit(target[0] - a[0], target[1] - a[1])
            end = (target[0] - ux * 0.25 * S, target[1] - uy * 0.25 * S)
            start = _exit(a, target, boxes[e["from"]], 0.4 * S)
        else:
            a = center(e["from"])
            b = center(e["to"])
            ux, uy = _unit(b[0] - a[0], b[1] - a[1])
            dist = math.hypot(b[0] - a[0], b[1] - a[1])
            c = ((a[0] + b[0]) / 2 + uy * bend * dist, (a[1] + b[1]) / 2 - ux * bend * dist)
            start = _exit(a, c, boxes[e["from"]], 0.4 * S)
            end = _exit(b, c, boxes[e["to"]], 0.4 * S)
        ux, uy = _unit(end[0] - start[0], end[1] - start[1])
        nlx, nly = uy, -ux                                   # traveller's left
        L = math.hypot(end[0] - start[0], end[1] - start[1])
        ctrl = ((start[0] + end[0]) / 2 + nlx * bend * L, (start[1] + end[1]) / 2 + nly * bend * L)
        geo[f"{e['from']}->{e['to']}"] = (start, ctrl, end)

        drawn.append((e, start, ctrl, end, (nlx, nly)))

    # every arrow is known now: place labels so they avoid nodes, arrows and each other
    points = []
    for _, s0, c0, e0, _ in drawn:                     # one sample every ~0.4 S along each arrow
        n = max(4, int(math.hypot(e0[0] - s0[0], e0[1] - s0[1]) / (0.4 * S)))
        points += [_bez(s0, c0, e0, t / n) for t in range(n + 1)]
    pad = 0.35 * S if text else 0.0                  # bare names need a margin around them
    obstacles = [(b[0] - pad, b[1] - pad, b[2] + pad, b[3] + pad) for b in boxes.values()]
    if p.get("center_text"):
        ctb = TextBlock(p["center_text"], "sans", style.NAME_SIZE * 1.1)
        h = ctb.bottom - ctb.top
        obstacles.append((-ctb.width / 2, -h / 2, ctb.width / 2, h / 2))
    for e, start, ctrl, end, (nlx, nly) in drawn:
        seg = []
        if e.get("reversible"):
            g = (0.2 if variant == "text" else 0.13) * S
            o = (nlx * g, nly * g)
            f = [(q[0] + o[0], q[1] + o[1]) for q in (start, ctrl, end)]
            r = [(q[0] - o[0], q[1] - o[1]) for q in (end, ctrl, start)]
            for (q0, qc, q1) in (f, r):
                seg.append(path(f"M{q0[0]:.2f},{q0[1]:.2f} Q{qc[0]:.2f},{qc[1]:.2f} {q1[0]:.2f},{q1[1]:.2f}"))
                seg.append(_half_head(q1, (q1[0] - qc[0], q1[1] - qc[1])))
        elif e.get("inhibit"):                       # drug / feedback inhibition: red ⊣
            seg.append(path(f"M{start[0]:.2f},{start[1]:.2f} Q{ctrl[0]:.2f},{ctrl[1]:.2f} {end[0]:.2f},{end[1]:.2f}",
                            stroke="red"))
            tx, ty = _unit(end[0] - ctrl[0], end[1] - ctrl[1])
            h = 0.4 * S
            seg.append(line(end[0] - ty * h, end[1] + tx * h, end[0] + ty * h, end[1] - tx * h, stroke="red",
                            w=2 * style.LINE_W))
        elif e.get("resonance"):                     # resonance structures: one line, a head at each end
            seg.append(line(*start, *end))
            seg.append(_head(end[0], end[1], end[0] - start[0], end[1] - start[1]))
            seg.append(_head(start[0], start[1], start[0] - end[0], start[1] - end[1]))
        else:
            seg.append(path(f"M{start[0]:.2f},{start[1]:.2f} Q{ctrl[0]:.2f},{ctrl[1]:.2f} {end[0]:.2f},{end[1]:.2f}"))
            if e.get("head", True):                   # head: false → plain line (family trees)
                seg.append(_head(end[0], end[1], end[0] - ctrl[0], end[1] - ctrl[1]))
        if e.get("color"):                           # e⁻ flow (orange), H⁺ flow (magenta) …
            seg = [x.replace(f'stroke="{style.INK}"', f'stroke="{style.color(e["color"])}"') for x in seg]
        for q in (start, ctrl, end):
            bb.add(q[0], q[1], q[0], q[1])

        m = _bez(start, ctrl, end, 0.5)
        tangent = _unit(end[0] - start[0], end[1] - start[1])
        side = e.get("label_side", default_side)
        nv = (nlx, nly) if side == "left" else (-nlx, -nly)
        if variant == "text":
            seg = [_thick(x) for x in seg]
        labels = []
        for key, kind, size, fill, normal, lay in (
                ("label", "sans", style.ENZYME_SIZE, "cyan", nv, "enzymes"),
                ("side_text", "serif", style.S * 0.85, None, (-nv[0], -nv[1]), "cofactors")):
            if not e.get(key):
                continue
            tb = TextBlock(e[key], kind, size)
            cx, cy, box = _place_block(tb, m, normal, tangent, obstacles, points)
            obstacles.append(box)
            svg = _draw_block(tb, cx, cy, fill, layer[lay][1])
            layer[lay][0].append(svg)
            labels.append(svg)
        arrow_only = "".join(seg)
        kind = "inhibitions" if e.get("inhibit") else EDGE_LAYER.get(e.get("color"), "arrows")
        fine[kind][0].append(arrow_only)
        for q in (start, ctrl, end):
            fine[kind][1].add(q[0] - 0.4 * S, q[1] - 0.4 * S, q[0] + 0.4 * S, q[1] + 0.4 * S)
        seg += labels

        if e.get("step"):
            d = f"M{start[0]:.2f},{start[1]:.2f} Q{ctrl[0]:.2f},{ctrl[1]:.2f} {end[0]:.2f},{end[1]:.2f}"
            hit = f'<path d="{d}" fill="none" stroke="transparent" stroke-width="{1.4 * S:.1f}"/>'
            parts.append(f'<g data-ref="reaction:{e["step"]}" class="sdb-step">{hit}{"".join(seg)}</g>')
        else:
            parts.append(f"<g>{''.join(seg)}</g>")
        path_body.append(arrow_only)

    if p.get("center_text"):
        tb = TextBlock(p["center_text"], "sans", style.NAME_SIZE * 1.1)
        parts.append(tb.svg(0, -(tb.top + tb.bottom) / 2))
        path_body.append(parts[-1])
    for _, lb in layer.values():
        if lb.x0 != float("inf"):
            bb.add(lb.x0, lb.y0, lb.x1, lb.y1)
    layers = {"path": Graphic(region_svg + "".join(path_body), bb.x0, bb.y0, bb.x1, bb.y1)}
    for k, (body, lb) in layer.items():
        if body:
            layers[k] = Graphic("".join(body), lb.x0, lb.y0, lb.x1, lb.y1)
    used = {k: v for k, v in fine.items() if v[0]}
    if len(used) > 1:                                 # one kind only → "path" already is that layer
        for k, (body, lb) in used.items():
            layers[k] = Graphic("".join(body), lb.x0, lb.y0, lb.x1, lb.y1)
    return Graphic("".join(parts), bb.x0, bb.y0, bb.x1, bb.y1, {}, layers)


NODE_LAYER = {"protein": "proteins", "carrier": "proteins", "drug": "drugs"}
EDGE_LAYER = {"orange": "electrons", "magenta": "protons"}


def _regions(p, cell):
    """Background drawn under the nodes (cells, like node positions):
        {kind: membrane, y: 2, x0: -5, x1: 5, label: "Inner membrane"}   horizontal bilayer
        {kind: membrane, x: 0, y0: -3, y1: 3}                              vertical bilayer
        {kind: box, x0, y0, x1, y1, label: "Liver"}                        compartment / organ
    """
    out, bb = [], BBox()
    T, R, STEP = 1.25 * S, 0.17 * S, 0.46 * S
    for rg in p.get("regions") or []:
        if rg["kind"] == "membrane":
            horiz = "y" in rg
            if horiz:
                a0, a1, c = rg["x0"] * cell, rg["x1"] * cell, rg["y"] * cell
                x0, y0, x1, y1 = a0, c - T / 2, a1, c + T / 2
            else:
                a0, a1, c = rg["y0"] * cell, rg["y1"] * cell, rg["x"] * cell
                x0, y0, x1, y1 = c - T / 2, a0, c + T / 2, a1
            out.append(f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{x1 - x0:.2f}" height="{y1 - y0:.2f}" '
                       f'fill="{style.color("membrane")}"/>')
            n = max(2, int((a1 - a0) / STEP))
            for k in range(n + 1):
                t = a0 + (a1 - a0) * k / n
                for side in (-1, 1):
                    hx, hy = (t, c + side * (T / 2 - R)) if horiz else (c + side * (T / 2 - R), t)
                    tx, ty = (t, c + side * 0.12 * S) if horiz else (c + side * 0.12 * S, t)
                    out.append(line(hx, hy, tx, ty, stroke="membrane_head", w=0.6 * style.LINE_W))
                    out.append(f'<circle cx="{hx:.2f}" cy="{hy:.2f}" r="{R:.2f}" fill="{style.color("membrane_head")}"/>')
            bb.add(x0, y0, x1, y1)
            if rg.get("label"):
                tb = TextBlock(rg["label"], "sans", style.ENZYME_SIZE)
                if horiz:
                    lx, ly = x1 + 0.4 * S + tb.width / 2, c - (tb.top + tb.bottom) / 2
                else:
                    lx, ly = c, y1 + 0.5 * S - tb.top
                out.append(tb.svg(lx, ly, "gray"))
                bb.add(lx - tb.width / 2, ly + tb.top, lx + tb.width / 2, ly + tb.bottom)
        elif rg["kind"] == "box":
            x0, y0, x1, y1 = rg["x0"] * cell, rg["y0"] * cell, rg["x1"] * cell, rg["y1"] * cell
            dash = ' stroke-dasharray="6 4"' if rg.get("dashed", True) else ""
            out.append(f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{x1 - x0:.2f}" height="{y1 - y0:.2f}" rx="{0.6 * S:.2f}" '
                       f'fill="none" stroke="{style.color("compartment")}" stroke-width="{style.LINE_W:.2f}"{dash}/>')
            bb.add(x0, y0, x1, y1)
            if rg.get("label"):
                tb = TextBlock(rg["label"], "sans", style.NAME_SIZE)
                lx, ly = x0 + 0.5 * S + tb.width / 2, y0 + 0.4 * S - tb.top
                out.append(tb.svg(lx, ly, "gray"))
    return "".join(out), bb


def compose_pathway(db, p, variant=None):
    scene = (p.get("variants") or {}).get(variant, {}).get("scene") if variant else None
    if scene or p["layout"] == "scene":            # illustrated version (glyphs placed by coordinates)
        from ..draw.scene import compose_scene
        return compose_scene(scene or p["scene"])
    if p["layout"] == "linear":
        return compose_linear(db, p, variant)
    if p["layout"] in ("free", "cycle"):
        return compose_graph(db, p, variant)
    if p["layout"] == "plot":                         # axes + curves (kinetics, O2 binding …)
        from ..draw.plot import compose_plot
        return compose_plot(db, p, variant)
    raise NotImplementedError(f"layout {p['layout']} 尚未實作")
