"""Sugar rings: Haworth projections (one ring or several linked rings) and chair conformations.

One ring (nucleotides, PRPP, monosaccharides):

    haworth:
      ring: furanose               # or pyranose
      start: 1                     # number of the first ring carbon (fructose: 2)
      subs:                        # carbon number → substituent above / below the ring plane
        1: {up: base, down: H}     # "base" = the ring given in `base:` (purine / pyrimidine)
        2: {up: H, down: OH}
        3: {up: H, down: OH}
        4: {up: "{P}OCH_2"}        # anchor atom (first C, else first heavy atom) sits on the bond;
                                   # ["^-O_3SO", 4] picks the anchor character explicitly
      base: "[*]n1cnc2c(N)ncnc21"  # attachment [*] = C1'
      base_right: 7                # atom index (in `base`) that must end up on the right
      base_len: 0.95               # C1′–base bond length in pitches (longer lifts the base clear)
      numbers: true                # small grey ring numbers (also exported as layer `numbers`)
      primes: true                 # 1′ 2′ … (default: only when a base is attached)
      mirror: false                # draw the ring left ↔ right (anomeric carbon on the left)
      flip: false                  # turned over (180° about the horizontal axis): ring O in front,
                                   # up ↔ down — the reducing-end glucose of lactose / cellobiose

Several rings (disaccharides): `units: [ring spec, …]`, each with `at: [dx, dy]` (pitches), and
`links: [{from: [unit, carbon, down], to: [unit, carbon, down], o: [x, y]}]` (o optional). A linked side is written
`link` in that unit's subs; the glycosidic O is drawn between the two carbons.

Chair (`chair:` instead of `haworth:`): `{ring: pyranose, subs: {1: {ax: H, eq: OH}, …}}`. The
ring is a real 3D chair (C1 the right tip, pointing down; C4 the left tip, pointing up; ring O
back right); axial bonds are vertical, equatorial bonds point outwards. The drawing is a
projection of the same 3D model the validator reads, so ax/eq mistakes are caught.

Front edges of a Haworth ring are drawn bold. Everything is built as a `free` drawing, so labels,
trimming and colours behave exactly like the other views.
"""
import math

from .. import style
from ..graphic import BBox, Graphic, line
from ..text import TextLine
from .chain import GAP_V, PITCH
from .free import draw_free
from .ring import draw_ring

S = style.S

# ring geometry in pitches: (x, y) with y growing towards the viewer (screen down)
RINGS = {
    "furanose": {"O": (0.0, -0.55), 1: (1.35, -0.02), 2: (0.68, 0.74), 3: (-0.68, 0.74), 4: (-1.35, -0.02)},
    # front corners sit further out than the back ones, so the H below C5 and the group above
    # C3 (both inside the ring) do not collide
    "pyranose": {"O": (0.55, -0.55), 1: (1.75, 0.05), 2: (1.0, 0.88), 3: (-1.0, 0.88), 4: (-1.75, 0.05),
                 5: (-0.55, -0.55)},
}
# substituent bond lengths (up, down) in pitches, per ring position
LENS = {
    "furanose": {1: (0.95, 0.66), 2: (0.55, 0.66), 3: (0.55, 0.66), 4: (0.95, 0.66)},
    "pyranose": {1: (0.75, 0.7), 2: (0.5, 0.66), 3: (0.5, 0.66), 4: (0.75, 0.7), 5: (0.95, 0.42)},
}


def ring_atoms(kind):
    """Ring order starting at O, then the ring carbons (used by the validator too)."""
    g = RINGS[kind]
    return ["O"] + [k for k in g if k != "O"]


def front(geo):
    """Edges nearest the viewer (drawn bold): those below the ring centre on screen."""
    keys = list(geo)
    cy = sum(p[1] for p in geo.values()) / len(geo)
    ring = ["O"] + [k for k in keys if k != "O"]
    return {(a, b) for a, b in zip(ring, ring[1:] + ring[:1]) if (geo[a][1] + geo[b][1]) / 2 > cy}


def units_of(spec):
    """Normalised list of ring units (a one-ring spec is a list of one)."""
    return spec["units"] if "units" in spec else [spec]


def unit_geo(u):
    """Ring positions (pitches) of one unit after mirror / flip and offset: {"O" | 1..n: (x, y)}."""
    geo = RINGS[u.get("ring", "furanose")]
    dx, dy = u.get("at", (0, 0))
    kx = -1 if u.get("mirror") else 1
    ky = -1 if u.get("flip") else 1
    return {key: (kx * x + dx, ky * y + dy) for key, (x, y) in geo.items()}


def sub_points(u):
    """[(position, side, text, (x, y))] — where each drawn substituent label sits."""
    kind = u.get("ring", "furanose")
    geo, start = unit_geo(u), u.get("start", 1)
    out = []
    for carbon, sub in (u.get("subs") or {}).items():
        pos = int(carbon) - start + 1
        x, y = geo[pos]
        up_len, down_len = LENS[kind][pos]
        if u.get("flip"):                              # front and back swap places
            up_len, down_len = down_len, up_len
        up_len, down_len = sub.get("up_len", up_len), sub.get("down_len", down_len)   # per-group override
        for side, dy in (("up", -up_len), ("down", down_len)):
            t = sub.get(side)
            if t is None or t == "":
                continue
            out.append((pos, side, t, (x + sub.get(f"{side}_dx", 0), y + dy)))   # *_dx: nudge a crowded label
        if sub.get("out"):                             # in the ring plane, pointing away (lactone C1=O)
            ring = list(geo.values())
            cx, cy = sum(p[0] for p in ring) / len(ring), sum(p[1] for p in ring) / len(ring)
            d = math.hypot(x - cx, y - cy) or 1.0
            out.append((pos, "out", sub["out"], (x + 0.85 * (x - cx) / d, y + 0.85 * (y - cy) / d)))
    return out


def _numbers(u, geo, n):
    nums, nb = [], BBox()
    primes = u.get("primes", bool(u.get("base")))
    start = u.get("start", 1)
    ring = [geo[i] for i in range(1, n + 1)] + [geo["O"]]
    cx0 = sum(p[0] for p in ring) / len(ring) * PITCH
    cy0 = sum(p[1] for p in ring) / len(ring) * PITCH
    for i in range(1, n + 1):
        x, y = geo[i][0] * PITCH, geo[i][1] * PITCH
        t = TextLine(f"{i + start - 1}{'′' if primes else ''}", "sans", 0.62 * S)
        # beside the vertex on the outer side, clear of its vertical substituent bonds
        out = 1 if x > cx0 else -1
        px = x + out * (0.5 * S + t.width / 2)
        py = y + (0.55 * S if y > cy0 else -0.45 * S)
        tx, base = px - t.width / 2, py + t.cap / 2
        nums.append(t.svg(tx, base, "gray"))
        nb.add(tx, base + t.top, tx + t.width, base + t.bottom)
    return nums, nb


def draw_haworth(spec, deco=None):
    deco = deco or {}
    units = units_of(spec)
    atoms, bonds = {}, []
    ends = {}                                     # (unit, pos, side) → ring-carbon atom id for links
    for ui, u in enumerate(units):
        kind = u.get("ring", "furanose")
        geo = unit_geo(u)
        n = len(geo) - 1
        p = f"u{ui}_" if len(units) > 1 else ""
        atoms[f"{p}O"] = ["O", *geo["O"]]
        for i in range(1, n + 1):
            atoms[f"{p}c{i}"] = ["", *geo[i]]
        keys = ["O"] + list(range(1, n + 1))
        fr = front(geo)
        for a, b in zip(keys, keys[1:] + keys[:1]):
            ida = f"{p}O" if a == "O" else f"{p}c{a}"
            idb = f"{p}O" if b == "O" else f"{p}c{b}"
            bonds.append([ida, idb, "bold" if (a, b) in fr or (b, a) in fr else 1])
        for pos, side, t, (x, y) in sub_points(u):
            if t == "base" and u.get("base"):
                continue
            if t == "link":
                ends[(ui, pos, side)] = f"{p}c{pos}"
                continue
            aid = f"{p}{side[0]}{pos}"
            # [text, i]: the bond meets character i (e.g. ["^-O_3SO", 4] → the last O)
            atoms[aid] = [t[0], x, y, t[1]] if isinstance(t, list) else [t, x, y]
            subs = u.get("subs") or {}
            sub = subs.get(pos + u.get("start", 1) - 1) or subs.get(str(pos + u.get("start", 1) - 1)) or {}
            bonds.append([f"{p}c{pos}", aid, sub.get(f"{side}_b", 1)])     # out_b: 2 → C=O

    for k, ln in enumerate(spec.get("links") or []):
        (ua, ca, sa), (ub, cb, sb) = ln["from"], ln["to"]
        pa = int(ca) - units[ua].get("start", 1) + 1
        pb = int(cb) - units[ub].get("start", 1) + 1
        ga, gb = unit_geo(units[ua])[pa], unit_geo(units[ub])[pb]
        drop = 0.75 if sa == "down" else -0.75
        oy = (max(ga[1], gb[1]) if sa == "down" else min(ga[1], gb[1])) + drop
        ox, oy = ln.get("o", ((ga[0] + gb[0]) / 2, oy))     # o: [x, y] places the glycosidic O by hand
        atoms[f"link{k}"] = [ln.get("label", "O"), ox, oy]
        bonds += [[ends[(ua, pa, sa)], f"link{k}"], [f"link{k}", ends[(ub, pb, sb)]]]

    colors = dict(deco.get("colors") or {})        # keys: "d2", "u4" (one ring) / "u1_d4" (units)
    g = draw_free({"atoms": atoms, "bonds": bonds}, {"colors": colors, "boxes": deco.get("boxes")})
    body, bb = g.body, BBox()
    bb.add(g.x0, g.y0, g.x1, g.y1)
    parts = dict(g.parts)

    for u in units:
        geo = unit_geo(u)
        base_at = next(((pos, side) for pos, side, t, _ in sub_points(u) if t == "base"), None)
        if base_at and u.get("base"):
            pos, side = base_at
            cx, cy = geo[pos][0] * PITCH, geo[pos][1] * PITCH
            up = side == "up"
            length = u.get("base_len", 0.95) * PITCH
            ay = cy - length if up else cy + length
            svg, rbb = draw_ring(u["base"], cx, ay, None, direction=90 if up else -90,
                                 flip=u.get("base_flip", False), right=u.get("base_right"),
                                 hl={**_base_colors(deco), **highlight_atoms(u["base"], deco.get("base_highlight"))})
            cap = TextLine("N").cap
            end = ay + (cap / 2 + GAP_V) * (1 if up else -1)
            body += svg + line(cx, cy, cx, end)
            bb.add(rbb.x0, rbb.y0, rbb.x1, rbb.y1)
            parts["base"] = Graphic(svg, rbb.x0, rbb.y0, rbb.x1, rbb.y1)

    nums, nb = [], BBox()
    for u in units:
        if u.get("numbers", spec.get("numbers")):
            geo = unit_geo(u)
            ns, b = _numbers(u, geo, len(geo) - 1)
            nums += ns
            nb.add(b.x0, b.y0, b.x1, b.y1)
    if nums:
        svg = "".join(nums)
        body += svg
        bb.add(nb.x0, nb.y0, nb.x1, nb.y1)
        parts["numbers"] = Graphic(svg, nb.x0, nb.y0, nb.x1, nb.y1)

    return Graphic(body, bb.x0, bb.y0, bb.x1, bb.y1, {"base": (0, bb.y1), "axis": (0, 0)}, parts)


def haworth_3d(u):
    """Validator coordinates of one unit: ring in the horizontal plane, up = +y, front = +z.
    Returns (ring {key: xyz}, {(pos, side): xyz of the substituent})."""
    geo = unit_geo(u)
    ring = {k: (x, 0.0, y) for k, (x, y) in geo.items()}
    cx = sum(p[0] for p in ring.values()) / len(ring)
    cz = sum(p[2] for p in ring.values()) / len(ring)
    subs = {}
    for pos in range(1, len(geo)):
        x, _y, z = ring[pos]
        d = math.hypot(x - cx, z - cz) or 1.0
        ox, oz = 0.35 * (x - cx) / d, 0.35 * (z - cz) / d   # tilt outwards: never collinear
        subs[(pos, "up")] = (x + ox, 1.0, z + oz)
        subs[(pos, "down")] = (x + ox, -1.0, z + oz)
    return ring, subs


# ---------------------------------------------------------------- chair

CHAIR_H = 0.30           # half the puckering height (ring radius 1; drawn a little exaggerated)
CHAIR_TILT = 0.38        # how much depth (y, back > 0) lifts a point: front edge low, back edge high


def chair_3d(spec):
    """Ring atoms of a 4C1 chair and their axial / equatorial directions, in a right-handed
    frame (x right, y back, z up). Order: C1 (right tip, down) C2 C3 C4 (left tip, up) C5 O.
    The drawing is an orthographic view of this model from slightly above."""
    keys = [1, 2, 3, 4, 5, "O"]
    ring, dirs = {}, {}
    for k, key in enumerate(keys):
        th = math.radians(-60 * k - 8)        # turned 8°: C2’s axial bond clears the ring O
        z = -CHAIR_H if k % 2 == 0 else CHAIR_H
        ring[key] = (math.cos(th), math.sin(th), z)        # th < 0 → y < 0 → front (C2, C3)
        sgn = 1 if z > 0 else -1
        ax = (0.0, 0.0, sgn)
        ex, ey = math.cos(th), math.sin(th)
        eq = (ex * 0.94, ey * 0.94, -sgn * 0.33)
        dirs[key] = {"ax": ax, "eq": eq}
    return ring, dirs


def _project(p):
    """Model units → pitches. Depth (y, back > 0) lifts a point; height z lifts it more."""
    x, y, z = p
    return x * 2.2, -2.2 * (z + y * CHAIR_TILT)


def draw_chair(spec, deco=None):
    deco = deco or {}
    ring, dirs = chair_3d(spec)
    start = spec.get("start", 1)
    atoms, bonds = {}, []
    for key, p in ring.items():
        X, Y = _project(p)
        atoms[str(key)] = ["O" if key == "O" else "", X, Y]
    order = ["O", "1", "2", "3", "4", "5"]
    for a, b in zip(order, order[1:] + order[:1]):
        bonds.append([a, b])
    for carbon, sub in (spec.get("subs") or {}).items():
        pos = int(carbon) - start + 1
        for kind in ("ax", "eq"):
            t = sub.get(kind)
            if not t:
                continue
            d = dirs[pos][kind]
            L = 0.45 if kind == "ax" else 0.75
            q = tuple(c + L * v for c, v in zip(ring[pos], d))
            aid = f"{kind}{pos}"
            atoms[aid] = [t, *_project(q)]
            bonds.append([str(pos), aid])
    g = draw_free({"atoms": atoms, "bonds": bonds}, {"colors": deco.get("colors") or {}})
    return Graphic(g.body, g.x0, g.y0, g.x1, g.y1, {"base": (0, g.y1), "axis": (0, 0)}, g.parts)


def chair_validation_3d(spec):
    """(ring {key: xyz}, {(pos, side): xyz}) with side "ax" / "eq", for the validator."""
    ring, dirs = chair_3d(spec)
    subs = {}
    for pos in range(1, 6):
        for kind in ("ax", "eq"):
            subs[(pos, kind)] = tuple(c + v for c, v in zip(ring[pos], dirs[pos][kind]))
    return ring, subs


# ---------------------------------------------------------------- standalone ring molecules

def _base_colors(deco):
    """`base_colors: {atom index in the base SMILES: colour}` (bonds between two atoms of the
    same colour follow); `base_highlight: [{smarts: …}]` does the same by SMARTS."""
    return {int(k): v for k, v in (deco.get("base_colors") or {}).items()}


def highlight_atoms(smiles, specs):
    """`highlight: [{smarts: "C(=O)N", color: magenta} | {atoms: [3, 4], color: cyan}]` →
    {atom index: colour}. Every SMARTS match is coloured; when the SMARTS carries atom maps
    ([N:1]), only the mapped atoms are."""
    from rdkit import Chem
    out = {}
    for h in specs or []:
        idx = list(h.get("atoms") or [])
        if h.get("smarts"):
            mol, patt = Chem.MolFromSmiles(smiles), Chem.MolFromSmarts(h["smarts"])
            keep = [a.GetIdx() for a in patt.GetAtoms() if a.GetAtomMapNum()] or range(patt.GetNumAtoms())
            idx += [match[k] for match in mol.GetSubstructMatches(patt) for k in keep]
        for i in idx:
            out[int(i)] = h.get("color", "magenta")
    return out


def draw_ring_molecule(m, deco=None, mirror=False):
    """Standalone ring molecule (adenine, uracil …):
    `ring: {smiles?, down | axis: [i, j], right, rotate, flip, explicit, nh_bond, h_down,
            numbers: {atom index: "9"}, mirror: true}`
    mirror: true adds the left/right twin as the view `ring_mirror`."""
    spec = m["ring"] if isinstance(m["ring"], dict) else {}
    smiles = spec.get("smiles", m["smiles"])
    layers = {}
    svg, bb = draw_ring(smiles, 0, 0, None, rotate=spec.get("rotate", 0), flip=spec.get("flip", False),
                        colors=_base_colors(deco or {}), down=spec.get("down"), right=spec.get("right"),
                        explicit=spec.get("explicit", False), numbers=spec.get("numbers"), layers=layers,
                        mirror=mirror, nh_bond=spec.get("nh_bond", False), axis=spec.get("axis"),
                        h_down=spec.get("h_down"), letters=spec.get("letters", False),
                        orient=spec.get("orient"), methyls=spec.get("methyls", True),
                        coordgen=spec.get("coordgen", True), seed=spec.get("seed"), hide_charges=tuple(spec.get("hide_charges", ())), hl=highlight_atoms(smiles, (deco or {}).get("highlight")))
    parts = {}
    if "numbers" in layers:
        nsvg, nb = layers["numbers"]
        parts["numbers"] = Graphic(nsvg, nb.x0, nb.y0, nb.x1, nb.y1)
    g = Graphic(svg, bb.x0, bb.y0, bb.x1, bb.y1, {}, parts)
    g = g.moved(-(g.x0 + g.x1) / 2, -(g.y0 + g.y1) / 2)
    g.anchors = {"base": (0, g.y1), "axis": (0, 0)}
    return g
