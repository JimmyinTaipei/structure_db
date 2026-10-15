"""Ring side chains (Phe, Tyr, Trp, His …) hung below a vertical chain.

The ring is given as SMILES with a dummy atom marking the attachment, e.g. "[*]c1ccc(O)cc1".
RDKit only supplies 2D coordinates; lines and labels are drawn here so that fonts and stroke
widths match the rest of the figure. The ring is rotated so that the attachment bond points
straight up, with the attachment atom at (0, 0).
"""
import math
import random

from rdkit import Chem
from rdkit.Chem import rdDepictor

from .. import style
from ..graphic import BBox, line
from ..text import TextLine

S = style.S
GAP_H = 0.14 * S
GAP_V = 0.20 * S


def _orient_steroid(mol, out):
    """Textbook steroid pose: rings A–B level at the bottom left, C and D raised to the right."""
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    five = [r for r in rings if len(r) == 5]
    six = [r for r in rings if len(r) == 6]
    if len(five) != 1 or len(six) != 3:
        return out
    D = five[0]
    C = next(r for r in six if r & D)
    B = next(r for r in six if r is not C and r & C)
    A = next(r for r in six if r is not C and r is not B)
    cen = lambda r: (sum(out[i][0] for i in r) / len(r), sum(out[i][1] for i in r) / len(r))  # noqa: E731
    a, b = cen(A), cen(B)
    rot = -math.atan2(b[1] - a[1], b[0] - a[0])           # A–B level, A on the left
    c, s_ = math.cos(rot), math.sin(rot)
    out = {i: ((x - a[0]) * c - (y - a[1]) * s_, (x - a[0]) * s_ + (y - a[1]) * c) for i, (x, y) in out.items()}
    if cen(C)[1] > 0:                                     # C must be raised (screen y < 0)
        out = {i: (x, -y) for i, (x, y) in out.items()}
    return out


def ring_geometry(smiles, bond_len, direction=-90, flip=False, rotate=0, down=None, right=None, mirror=False,
                  axis=None, orient=None, methyls=True, coordgen=True, seed=None):
    """direction: angle (degrees, screen coords, -90 = up) of the attachment bond attach → [*].
    Without a dummy atom the molecule is centred on its first atom; `down` (atom index) turns
    the molecule so that atom points straight down (N9 of a purine, N1 of a pyrimidine), else
    `rotate` (degrees) applies. `right` (atom index): mirror if needed so that atom ends up on
    the right-hand side (C2 of a pyrimidine), which overrides `flip`."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"ring SMILES 無法解析：{smiles}")
    Chem.Kekulize(mol, clearAromaticFlags=True)
    # seed: lay out a copy with shuffled atom order (RDKit's layout depends on it), then map the
    # coordinates back, so atom indices in the YAML (highlight, numbers, down …) stay valid
    order = list(range(mol.GetNumAtoms()))
    if seed:
        random.Random(seed).shuffle(order)
    work = Chem.RenumberAtoms(mol, order) if seed else mol
    rdDepictor.SetPreferCoordGen(coordgen)
    rdDepictor.Compute2DCoords(work)
    rdDepictor.SetPreferCoordGen(True)
    conf = work.GetConformer()
    pts = {order[k]: (conf.GetAtomPosition(k).x, -conf.GetAtomPosition(k).y) for k in range(work.GetNumAtoms())}
    # the attachment point is a plain [*]; numbered dummies ([16*] = rest of glutathione …) are
    # drawn as text labels (R_LABELS), like the R groups of skeletal drawings
    dummy = next((a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and not a.GetIsotope()), None)
    if dummy is None:
        attach = 0
        ax, ay = pts[attach]
        rot = math.radians(rotate)
        if down is not None:
            ring = next((r for r in mol.GetRingInfo().AtomRings() if down in r), list(pts))
            gx = sum(pts[i][0] for i in ring) / len(ring)
            gy = sum(pts[i][1] for i in ring) / len(ring)
            rot = math.radians(90) - math.atan2(pts[down][1] - gy, pts[down][0] - gx)
        if axis is not None:                       # bond axis[0] → axis[1] points straight up
            i, j = axis
            rot = math.radians(-90) - math.atan2(pts[j][1] - pts[i][1], pts[j][0] - pts[i][0])
    else:
        attach = mol.GetAtomWithIdx(dummy).GetNeighbors()[0].GetIdx()
        ax, ay = pts[attach]
        vx, vy = pts[dummy][0] - ax, pts[dummy][1] - ay
        # flipping mirrors x afterwards, so aim at the mirrored direction first
        target = math.radians(180 - direction if flip else direction)
        rot = target - math.atan2(vy, vx)
    lens = [math.dist(pts[b.GetBeginAtomIdx()], pts[b.GetEndAtomIdx()])
            for b in mol.GetBonds() if dummy not in (b.GetBeginAtomIdx(), b.GetEndAtomIdx())]
    k = bond_len / (sum(lens) / len(lens))
    c, s_ = math.cos(rot), math.sin(rot)
    out = {}
    for i, (x, y) in pts.items():
        dx, dy = x - ax, y - ay
        x, y = (dx * c - dy * s_) * k, (dx * s_ + dy * c) * k
        out[i] = (-x if flip else x, y)
    if right is not None:
        if dummy is None and axis is not None:
            ref = out[axis[0]][0]
        elif dummy is None and down is not None:
            ref = out[down][0]
        else:
            ref = out[attach][0]
        if out[right][0] < ref:
            out = {i: (-x, y) for i, (x, y) in out.items()}
    if orient == "steroid":
        out = _orient_steroid(mol, out)
    if mirror:                                     # the left/right twin of a base drawing
        out = {i: (-x, y) for i, (x, y) in out.items()}

    ri = mol.GetRingInfo()
    centres = [tuple(sum(out[i][d] for i in ring) / len(ring) for d in (0, 1)) for ring in ri.AtomRings()]
    # the centre each ring atom's number leans towards (a fused atom prefers its 6-membered ring)
    home = {}
    for ring, cen in sorted(zip(ri.AtomRings(), centres), key=lambda rc: -len(rc[0])):
        for i in ring:
            home.setdefault(i, cen)
    from .skeletal import R_LABELS
    atoms = [(a.GetIdx(), R_LABELS.get(a.GetIsotope(), "R") if a.GetAtomicNum() == 0 else a.GetSymbol(),
              a.GetTotalNumHs(), a.GetFormalCharge(),
              methyls and a.GetSymbol() == "C" and not a.IsInRing() and a.GetDegree() == 1,  # CH₃ label
              [n.GetIdx() for n in a.GetNeighbors()], home.get(a.GetIdx()))
             for a in mol.GetAtoms() if a.GetIdx() != dummy]
    bonds = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if dummy in (i, j):
            continue
        rings = [n for n, r in enumerate(ri.AtomRings()) if i in r and j in r]
        bonds.append((i, j, 2 if b.GetBondType() == Chem.BondType.DOUBLE else 1,
                      centres[rings[0]] if rings else None))
    ring_geometry.last_centres = sorted(centres)       # for ring letters (A–D of a steroid)
    return atoms, bonds, out


MIN_VISIBLE = 0.45 * S          # shortest visible stroke between two atom labels


def _stretch_label_bonds(atoms, bonds, pos):
    """A chain bond between two labelled atoms (P–O, O–P, S–S …) loses most of its length to the
    two labels; on a 45° bond nothing is left (phosphates). For every such acyclic bond whose
    visible part is too short, slide the smaller side of the molecule straight out along the
    bond, so every angle is kept and only that bond gets longer."""
    pos = dict(pos)
    nbrs = {idx: list(n) for idx, _s, _h, _c, _m, n, _home in atoms}
    half = {}
    for idx, sym, _nh, _ch, methyl, _n, _home in atoms:
        if sym != "C" or methyl:
            t = TextLine(sym)
            x0, x1 = t.char_box(0)
            half[idx] = ((x1 - x0) / 2 + GAP_H, t.cap / 2 + GAP_V)

    def side(start, block):
        seen, todo = {start}, [start]
        while todo:
            a = todo.pop()
            for b in nbrs[a]:
                if b not in seen and b != block:
                    seen.add(b)
                    todo.append(b)
        return seen

    for i, j, _order, centre in bonds:
        if centre is not None or i not in half or j not in half:
            continue
        (x1, y1), (x2, y2) = pos[i], pos[j]
        L = math.hypot(x2 - x1, y2 - y1)
        ux, uy = (x2 - x1) / L, (y2 - y1) / L
        cut = sum(min(hw / abs(ux) if ux else 1e9, hh / abs(uy) if uy else 1e9) for hw, hh in (half[i], half[j]))
        if L - cut >= MIN_VISIBLE:
            continue
        d = MIN_VISIBLE - (L - cut)
        far, near = side(j, i), side(i, j)
        if far & near:                               # in a ring after all: leave it
            continue
        move, sgn = (far, 1) if len(far) <= len(near) else (near, -1)
        for a in move:
            pos[a] = (pos[a][0] + sgn * ux * d, pos[a][1] + sgn * uy * d)
    return pos


def label_overlaps(smiles, spec):
    """Pairs of atom labels (not bonded to each other) whose boxes touch in the `ring:` drawing."""
    import itertools
    atoms, bonds, pos = ring_geometry(smiles, 0.82 * 1.75 * S, -90, spec.get("flip", False), spec.get("rotate", 0),
                                      spec.get("down"), spec.get("right"), False, spec.get("axis"), spec.get("orient"),
                                      spec.get("methyls", True), spec.get("coordgen", True), spec.get("seed"))
    pos = _stretch_label_bonds(atoms, bonds, pos)
    box = {}
    for idx, sym, nh, charge, methyl, nbrs, _home in atoms:
        if sym == "C" and not methyl:
            continue
        away = pos[idx][0] - sum(pos[n][0] for n in nbrs) / len(nbrs) if nbrs else 0.0
        if sym in spec.get("hide_charges", ()):
            charge = 0
        t = TextLine(_label(sym, nh, charge, away))
        ai = next(i for i, (ch, *_r) in enumerate(t.chars) if ch == sym[0])
        x0, x1 = t.char_box(ai)
        if len(sym) == 2:
            x1 = t.char_box(ai + 1)[1]
        tx = pos[idx][0] - (x0 + x1) / 2
        box[idx] = (tx, pos[idx][1] - t.cap / 2, tx + t.width, pos[idx][1] + t.cap / 2)
    bonded = {frozenset(b[:2]) for b in bonds}
    g = 0.08 * S
    return [(i, j) for i, j in itertools.combinations(box, 2) if frozenset((i, j)) not in bonded
            and box[i][0] < box[j][2] + g and box[j][0] < box[i][2] + g
            and box[i][1] < box[j][3] + g and box[j][1] < box[i][3] + g]


def _label(sym, nh, charge, x):
    h = "" if nh == 0 else ("H" if nh == 1 else f"H_{nh}")
    ch = "^+" if charge > 0 else ("^-" if charge < 0 else "")
    # heteroatoms left of the attachment axis read right-to-left: HN, HO
    text = f"{h}{sym}" if (h and x < -0.3 * S) else f"{sym}{h}"     # nearly vertical → NH₂
    return text + ch


def draw_ring(smiles, ox, oy, color=None, bond_len=None, direction=-90, flip=False, rotate=0, colors=None,
              down=None, right=None, explicit=False, numbers=None, layers=None, mirror=False,
              nh_bond=False, axis=None, h_down=None, letters=False, orient=None, methyls=True,
              coordgen=True, hl=None, seed=None, hide_charges=()):
    """Draw the ring with its attachment atom at (ox, oy). Returns (svg, BBox).
    colors: {atom index in the SMILES: colour} for single heteroatom labels.
    explicit: ring carbons are written as C and ring N–H / O–H get their H as a separate atom
    on a bond pointing outwards (textbook base drawings). numbers: {atom index: "9"} small grey
    ring numbers inside the ring; layers (dict) receives them as a separate "numbers" part.
    nh_bond: skeletal ring (carbons unlabelled) but ring N–H drawn as an N–H bond, as in
    textbook base drawings. hl: {atom index: colour} recolours those atoms' labels and every
    bond between two atoms of the same colour (the reactive part of a coenzyme)."""
    hl = {int(k): v for k, v in (hl or {}).items()}
    colors = {**{int(k): v for k, v in (colors or {}).items()}, **hl}
    bond_len = bond_len or 0.82 * 1.75 * S
    if explicit:                                   # every ring atom carries a letter: longer bonds
        bond_len *= 1.45
    atoms, bonds, pos = ring_geometry(smiles, bond_len, direction, flip, rotate, down, right, mirror, axis,
                                      orient, methyls, coordgen, seed)
    if not explicit:
        pos = _stretch_label_bonds(atoms, bonds, pos)
    parts, bb = [], BBox()
    labels = {}
    extra_h = []                                   # (atom idx, H position) for explicit ring N–H
    nums, nb = [], BBox()
    for idx, sym, nh, charge, methyl, nbrs, home in atoms:
        x, y = pos[idx][0] + ox, pos[idx][1] + oy
        bb.add(x, y, x, y)
        n = (numbers or {}).get(idx, (numbers or {}).get(str(idx)))
        if n is not None and home is not None:
            t = TextLine(str(n), "sans", 0.6 * S)
            px = x + (home[0] + ox - x) * 0.36
            py = y + (home[1] + oy - y) * 0.36
            nums.append(t.svg(px - t.width / 2, py + t.cap / 2, "gray"))
            nb.add(px - t.width / 2, py - t.cap / 2, px + t.width / 2, py + t.cap / 2)
        in_ring = home is not None
        if sym == "C" and not methyl and not (explicit and in_ring):
            continue
        if (explicit or nh_bond) and in_ring and nh and sym != "C":
            # N–H drawn as a bond: H sits outside the ring, away from the neighbours
            mx = sum(pos[n][0] for n in nbrs) / len(nbrs)
            my = sum(pos[n][1] for n in nbrs) / len(nbrs)
            d = math.hypot(pos[idx][0] - mx, pos[idx][1] - my) or 1.0
            ux, uy = (pos[idx][0] - mx) / d, (pos[idx][1] - my) / d
            if h_down is not None and idx == h_down:
                ux, uy = 0.0, 1.0                      # N9–H straight down (base cards)
            k = 0.72 if explicit else 1.2
            extra_h.append((idx, (pos[idx][0] + ux * k * bond_len, pos[idx][1] + uy * k * bond_len)))
            nh = 0
        if explicit and in_ring and sym == "C":
            nh = 0
        # H goes on the side facing away from the bonded neighbours (HN on the left of a ring)
        away = pos[idx][0] - sum(pos[n][0] for n in nbrs) / len(nbrs) if nbrs else 0.0
        if sym in hide_charges:                    # heme: plain "Fe" and "N", as in textbooks
            charge = 0
        t = TextLine(_label(sym, nh, charge, away))
        ai = next(i for i, (ch, *_r) in enumerate(t.chars) if ch == sym[0])
        x0, x1 = t.char_box(ai)
        if len(sym) == 2:                          # Fe, Co, Mg: centre on both letters
            x1 = t.char_box(ai + 1)[1]
        tx = x - (x0 + x1) / 2
        base = y + t.cap / 2
        parts.append(t.svg(tx, base, (colors or {}).get(idx, color)))
        bb.add(tx, base + t.top, tx + t.width, base + t.bottom)
        labels[idx] = ((x1 - x0) / 2, t.cap / 2)

    def trim(i, ux, uy):
        if i not in labels:
            return 0.0
        hw, hh = labels[i]
        return min((hw + GAP_H) / abs(ux) if ux else 1e9, (hh + GAP_V) / abs(uy) if uy else 1e9)

    for idx, (hx, hy) in extra_h:
        t = TextLine("H")
        x0, x1 = t.char_box(0)
        tx, base = hx + ox - (x0 + x1) / 2, hy + oy + t.cap / 2
        parts.append(t.svg(tx, base, (colors or {}).get(idx, color)))
        bb.add(tx, base + t.top, tx + t.width, base + t.bottom)
        key = f"h{idx}"
        pos[key] = (hx, hy)
        labels[key] = ((x1 - x0) / 2, t.cap / 2)
        bonds = list(bonds) + [(idx, key, 1, None)]
    if letters:                                     # A B C D … from left to right (steroid rings)
        for k, (cx, cy) in enumerate(ring_geometry.last_centres):
            t = TextLine("ABCDEFGH"[k], "sans", 0.75 * S)
            nums.append(t.svg(cx + ox - t.width / 2, cy + oy + t.cap / 2, "gray"))
            nb.add(cx + ox - t.width / 2, cy + oy - t.cap / 2, cx + ox + t.width / 2, cy + oy + t.cap / 2)
    if nums:
        svg = "".join(nums)
        parts.append(svg)
        if layers is not None:
            layers["numbers"] = (svg, nb)

    for i, j, order, centre in bonds:
        (x1, y1), (x2, y2) = pos[i], pos[j]
        L = math.hypot(x2 - x1, y2 - y1)
        ux, uy = (x2 - x1) / L, (y2 - y1) / L
        a, b = trim(i, ux, uy), trim(j, ux, uy)
        p1 = (x1 + ux * a + ox, y1 + uy * a + oy)
        p2 = (x2 - ux * b + ox, y2 - uy * b + oy)
        bc = hl.get(i) if isinstance(i, int) and hl.get(i) == hl.get(j if isinstance(j, int) else i) else None
        parts.append(line(*p1, *p2, bc))
        if order == 2:
            nx, ny = -uy, ux
            if centre is not None:
                # second line inside the ring, a little shorter (Kekulé style)
                mx, my = (x1 + x2) / 2, (y1 + y2) / 2
                if (centre[0] - mx) * nx + (centre[1] - my) * ny < 0:
                    nx, ny = -nx, -ny
                off, sh = 0.17 * bond_len, 0.16 * bond_len
                q1 = (p1[0] + nx * off + ux * sh * (a == 0), p1[1] + ny * off + uy * sh * (a == 0))
                q2 = (p2[0] + nx * off - ux * sh * (b == 0), p2[1] + ny * off - uy * sh * (b == 0))
                parts.append(line(*q1, *q2, bc))
            else:
                d = 0.085 * S
                parts.append(line(p1[0] + nx * d, p1[1] + ny * d, p2[0] + nx * d, p2[1] + ny * d, bc))
    return "".join(parts), bb
