"""Lipids built from a small spec: zigzag fatty acids, triacylglycerols, glycerophospholipids,
sphingolipids. The same spec produces the drawing and the SMILES, so the validator can compare
the molecule's SMILES (including every cis/trans double bond) with what is drawn.

    lipid:
      type: fatty_acid                     # fatty_acid | tag | phospholipid | sphingolipid
      chains: [{n: 18, cis: [9]}]          # carbons counted from the carboxyl C (Δ numbering)
      head: choline                        # phospholipid: choline | ethanolamine | serine | inositol |
                                           #   glycerol | phosphate;  sphingolipid: oh |
                                           #   phosphocholine | galactose | glucose | sulfogalactose
      numbers: true                        # grey C1, Δ carbons and the last carbon, plus ω

Zigzag rule: every bond turns ±60° from the previous one, alternating; a cis double bond makes
the two turns at its ends the same sign, so the chain kinks by 60°. With three or more cis bonds
the kinks are softened (total ≤ ~180°, each ≥ 36°) so the chain opens into a U as in textbook
arachidonate instead of closing on itself. trans keeps the plain zigzag.
"""
import math

from .. import style
from ..graphic import BBox, Graphic
from ..text import TextLine
from .free import draw_free

S = style.S
STEP = 0.62             # zigzag bond length, in pitches
ROW = 1.75              # distance between glycerol / sphingosine rows, in pitches

HEADS = {   # phospholipid head groups on the phosphate: (label, SMILES from the atom bonded to O)
    "choline": ("CH_2CH_2N^+(CH_3)_3", "CC[N+](C)(C)C"),
    "ethanolamine": ("CH_2CH_2NH_3^+", "CC[NH3+]"),
    "serine": ("CH_2CH(NH_3^+)COO^-", "CC([NH3+])C(=O)[O-]"),
    "inositol": ("Inositol", "C1C(O)C(O)C(O)C(O)C1O"),
    "glycerol": ("CH_2CH(OH)CH_2OH", "CC(O)CO"),
}
S_HEADS = {  # sphingolipid substituent on the C1 oxygen: (drawing kind, SMILES from that O)
    "oh": ("OH", "O"),
    "phosphocholine": ("phospho", "OP(=O)([O-])OCC[N+](C)(C)C"),
    "galactose": ("β-Galactose", "OC1OC(CO)C(O)C(O)C1O"),
    "glucose": ("β-Glucose", "OC1OC(CO)C(O)C(O)C1O"),
    "sulfogalactose": ("3-O-Sulfo-β-Galactose", "OC1OC(CO)C(O)C(OS(=O)(=O)[O-])C1O"),
}


# ---------------------------------------------------------------- chains

def _dbl(chain):
    d = {k: "cis" for k in chain.get("cis", [])}
    d.update({k: "trans" for k in chain.get("trans", [])})
    return d


def chain_smiles(start, end, dbl):
    """Atoms start+1 … end written after the caller has written atom `start` (Δk = k=k+1)."""
    out = []
    for j in range(start + 1, end + 1):
        if j in dbl:
            b = "/"
        elif j - 1 in dbl:
            b = "="
        elif j - 2 in dbl:
            b = "\\" if dbl[j - 2] == "cis" else "/"
        else:
            b = ""
        out.append(b + "C")
    return "".join(out)


def zigzag(x0, y0, start, end, dbl, heading=30.0, theta=None):
    """Positions (pitches) of atoms start+1 … end, the bond start→start+1 leaving at `heading`
    (degrees, screen coordinates: +30 = down-right). Returns {atom number: (x, y)} and the
    heading of every bond {k: degrees of bond k→k+1}.

    Each bond points along the chain direction θ ± 30° (the zigzag side s alternates). A cis
    double bond turns θ by the kink angle and keeps the next bond on the same side."""
    n_cis = sum(1 for v in dbl.values() if v == "cis")
    kink = 60.0 if n_cis < 3 else max(36.0, 180.0 / n_cis)
    if theta is None:                                   # horizontal chain
        side = 1 if heading > 0 else -1
        theta = heading - 30.0 * side
    else:
        side = 1 if heading > theta else -1
    pos, heads = {start: (x0, y0)}, {}
    x, y = x0, y0
    for j in range(start, end):
        if j > start:
            if dbl.get(j - 1) == "cis":
                theta += kink * side
            else:
                side = -side
        h = theta + 30.0 * side
        heads[j] = h
        x += STEP * math.cos(math.radians(h))
        y += STEP * math.sin(math.radians(h))
        pos[j + 1] = (x, y)
    return pos, heads


def _side(p, q, r):
    return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])


def check_geometry(pos, dbl):
    """Each drawn double bond must really look cis / trans."""
    bad = []
    for k, kind in dbl.items():
        if k - 1 in pos and k + 2 in pos:
            same = _side(pos[k], pos[k + 1], pos[k - 1]) * _side(pos[k], pos[k + 1], pos[k + 2]) > 0
            if same != (kind == "cis"):
                bad.append(f"Δ{k} 畫成 {'cis' if same else 'trans'}，但規格是 {kind}")
    return bad


# ---------------------------------------------------------------- builder

class Build:
    def __init__(self):
        self.atoms, self.bonds, self.nums, self.problems = {}, [], [], []

    def atom(self, aid, text, x, y, anchor=None):
        self.atoms[aid] = [text, x, y] + ([anchor] if anchor is not None else [])
        return aid

    def bond(self, a, b, order=1):
        self.bonds.append([a, b, order] if order != 1 else [a, b])

    def chain(self, prefix, at, start, end, dbl, heading=30.0, numbers=False, omega=False, down=False):
        """Zigzag carbons start+1…end hanging from atom id `at` (which sits at chain atom `start`).
        down: the chain runs straight down (θ = 90°); the zigzag starts on the side that makes
        the first cis kink bend to the right, into the free space."""
        x0, y0 = self.atoms[at][1], self.atoms[at][2]
        theta = None
        if down:
            theta = 90.0
            cis = sorted(k for k, v in dbl.items() if v == "cis")
            # side of bond k = side0 · (−1)^(k − start − 1) until the first kink; θ += kink·side,
            # and from straight down (θ = 90°) a larger θ swings the tail towards +x on screen
            side0 = -1 if not cis else (-1) ** (cis[0] - start - 1)
            heading = theta + 30.0 * side0
        pos, heads = zigzag(x0, y0, start, end, dbl, heading, theta)
        self.problems += check_geometry(pos, dbl)
        prev = at
        for j in range(start + 1, end + 1):
            aid = self.atom(f"{prefix}{j}", "", *pos[j])
            self.bond(prev, aid, 2 if (j - 1) in dbl else 1)
            prev = aid
        if numbers:
            show = {start, end} | {k for k in dbl} | {k + 1 for k in dbl}
            for j in sorted(show):
                if j == start or j not in pos:
                    continue
                hin = heads.get(j - 1)
                hout = heads.get(j, hin)
                ix, iy = math.cos(math.radians(hin)), math.sin(math.radians(hin))
                ox, oy = math.cos(math.radians(hout)), math.sin(math.radians(hout))
                bx, by = ix - ox, iy - oy
                d = math.hypot(bx, by)
                if d < 1e-6:                                 # chain end / straight: beside the bond
                    bx, by, d = -iy, ix, 1.0
                px, py = pos[j][0] + 0.42 * bx / d, pos[j][1] + 0.42 * by / d
                self.nums.append((str(j), px, py))
            if omega:
                ex, ey = pos[end]
                hx, hy = math.cos(math.radians(heads[end - 1])), math.sin(math.radians(heads[end - 1]))
                self.nums.append(("ω", ex + 0.5 * hx + 0.3 * hy, ey + 0.5 * hy - 0.3 * hx))
        return pos

    def ester(self, prefix, gly, x, chain, numbers, heading=60.0):
        """Glycerol carbon → O → C(=O) hanging below at column x; the acyl tail zigzags down."""
        o = self.atom(f"{prefix}o", "O", x, 1.0)
        self.bond(gly, o)
        c1 = self.atom(f"{prefix}c1", "", x, 1.95)
        self.bond(o, c1)
        self.bond(c1, self.atom(f"{prefix}oo", "O", x - 0.75, 1.95), 2)
        self.chain(f"{prefix}c", c1, 1, chain["n"], _dbl(chain), numbers=numbers, down=True)

    def phospho_head(self, carbon, x, head):
        """carbon → O → P(=O)(O⁻) → O → head, stacked upwards at column x."""
        o = self.atom("po", "O", x, -1.0)
        self.bond(carbon, o)
        p = self.atom("p", "P", x, -2.0)
        self.bond(o, p)
        self.bond(p, self.atom("pdo", "O", x - 0.75, -2.0), 2)
        self.bond(p, self.atom("pom", "O^-", x + 0.8, -2.0))
        if head == "phosphate":
            self.bond(p, self.atom("pe", "O^-", x, -3.0))
            return
        o2 = self.atom("pe", "O", x, -3.0)
        self.bond(p, o2)
        label = HEADS[head][0] if head in HEADS else head
        self.bond(o2, self.atom("head", label, x, -4.0))


GX = (0.0, 1.8, 3.6)     # x of glycerol C1, C2, C3 (sphingosine C3, C2, C1), in pitches


def build(spec):
    """→ (free-drawing spec, SMILES, problems, number marks)."""
    kind = spec.get("type", "fatty_acid")
    chains = spec.get("chains") or []
    nums = spec.get("numbers", False)
    b = Build()
    if kind == "fatty_acid":
        ch = chains[0]
        c1 = b.atom("c1", "", 0, 0)
        o1 = b.atom("o1", "O", 0, -0.72)
        end = spec.get("end", "carboxylate")           # carboxylate | coa (thioester) | carnitine (ester)
        b.bond(c1, o1, 2)
        if end == "carboxylate":
            b.bond(c1, b.atom("o2", "^-O", -0.66, 0.38))
            head = "[O-]C(=O)"
        else:
            x = b.atom("o2", "S" if end == "coa" else "O", -0.66, 0.38)
            b.bond(c1, x)
            label = "CoA" if end == "coa" else "Carnitine"
            b.bond(x, b.atom("rest", label, -1.62, 0.38, len(label) - 1))     # text runs to the left
            head = "[1*]SC(=O)" if end == "coa" else "C[N+](C)(C)C[C@@H](CC(=O)[O-])OC(=O)"
        b.chain("c", c1, 1, ch["n"], _dbl(ch), numbers=nums, omega=nums)
        if nums:
            b.nums.append(("1", 0.32, 0.4))
        smi = head + chain_smiles(1, ch["n"], _dbl(ch))
    elif kind in ("tag", "phospholipid"):
        # glycerol left → right (sn-1, sn-2, sn-3); acyl tails hang down, head group points up
        labels = ["CH_2", "CH", "CH_2"]
        g = [b.atom(f"g{i}", t, GX[i], 0) for i, t in enumerate(labels)]
        b.bond(g[0], g[1])
        b.bond(g[1], g[2])
        acyl = []
        n_acyl = 3 if kind == "tag" else 2
        for i in range(n_acyl):
            b.ester(f"a{i}", g[i], GX[i], chains[i], nums)
            acyl.append(chain_smiles(1, chains[i]["n"], _dbl(chains[i])))
        smi = f"C(OC(=O){acyl[0]})C(OC(=O){acyl[1]})CO"
        if kind == "tag":
            smi += f"C(=O){acyl[2]}"
        else:
            head = spec.get("head", "phosphate")
            b.phospho_head(g[2], GX[2], head)
            smi += "P(=O)([O-])[O-]" if head == "phosphate" else f"P(=O)([O-])O{HEADS[head][1]}"
    elif kind == "sphingolipid":
        # sphingosine C3 (HO–, trans-Δ4 tail down) – C2 (amide tail down) – C1 (head up)
        c3 = b.atom("s3", "CH", GX[0], 0)
        c2 = b.atom("s2", "CH", GX[1], 0)
        c1 = b.atom("s1", "CH_2", GX[2], 0)
        b.bond(c3, c2)
        b.bond(c2, c1)
        b.bond(c3, b.atom("s3o", "HO", GX[0] - 1.15, 0))
        c4 = b.atom("sph4", "", GX[0], 1.0)
        b.bond(c3, c4)
        b.chain("sph", c4, 4, 18, {4: "trans"}, down=True)
        nh = b.atom("n", "NH", GX[1], 1.0)
        b.bond(c2, nh)
        cc = b.atom("ac1", "", GX[1], 1.95)
        b.bond(nh, cc)
        b.bond(cc, b.atom("aco", "O", GX[1] + 0.75, 1.95), 2)
        ch = chains[0]
        b.chain("ac", cc, 1, ch["n"], _dbl(ch), numbers=nums, down=True)
        head = spec.get("head", "oh")
        draw, hsmi = S_HEADS[head]
        if draw == "OH":
            b.bond(c1, b.atom("h", "OH", GX[2], -1.0))
        elif draw == "phospho":
            b.phospho_head(c1, GX[2], "choline")
        else:
            o = b.atom("po", "O", GX[2], -1.0)
            b.bond(c1, o)
            b.bond(o, b.atom("head", draw, GX[2], -2.0))
        acyl = chain_smiles(1, ch["n"], _dbl(ch))
        smi = f"C({hsmi})C(NC(=O){acyl})C(O)" + chain_smiles(3, 18, {4: "trans"})
    else:
        raise ValueError(f"lipid type {kind}")
    return {"atoms": b.atoms, "bonds": b.bonds}, smi, b.problems, b.nums


def draw_lipid(spec, deco=None):
    free, _smi, _p, nums = build(spec)
    g = draw_free(free, {"colors": (deco or {}).get("colors") or {}})
    parts = dict(g.parts)
    body, bb = g.body, BBox()
    bb.add(g.x0, g.y0, g.x1, g.y1)
    if nums:
        from .chain import PITCH
        svgs, nb = [], BBox()
        for text, x, y in nums:
            t = TextLine(text, "sans", 0.6 * S)
            px, py = x * PITCH, y * PITCH
            svgs.append(t.svg(px - t.width / 2, py + t.cap / 2, "gray"))
            nb.add(px - t.width / 2, py - t.cap / 2, px + t.width / 2, py + t.cap / 2)
        svg = "".join(svgs)
        body += svg
        bb.add(nb.x0, nb.y0, nb.x1, nb.y1)
        parts["numbers"] = Graphic(svg, nb.x0, nb.y0, nb.x1, nb.y1)
    return Graphic(body, bb.x0, bb.y0, bb.x1, bb.y1, {"base": (0, bb.y1), "axis": (0, 0)}, parts)
