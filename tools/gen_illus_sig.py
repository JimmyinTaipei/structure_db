"""Round 25: illustrated (`illus`) versions of the membrane, signalling and energy figures, plus two new membrane scenes.

The boxed versions stay as they are; this script only adds `variants.illus` to the existing YAML files (keeping the review
status and every other field), so run it AFTER gen_sig.py / gen_notes.py / gen_energy.py whenever those are re-run.
Figures carry names and short labels only; explanations stay in notes_zh.

    python tools/gen_illus_sig.py data
"""
import math
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent))
from figure_text import tidy  # noqa: E402

H_COL, E_COL = "magenta", "#E07B1A"        # proton flow, electron flow


# ================================================================ helpers

def T(x, y, s, size=0.75, color=None, align="center"):
    d = {"g": "text", "at": [x, y], "text": s, "size": size, "align": align}
    if color:
        d["color"] = color
    return d


def B(x, y, w, h, hue, text="", shape="oval", ts=0.6, **kw):
    d = {"g": "blob", "at": [x, y], "w": w, "h": h, "shape": shape, "hue": hue, **kw}
    if text:
        d["text"], d["text_size"] = text, ts
    return d


def I(x, y, text="", hue="yellow", r=0.5, ts=0.45):
    d = {"g": "ion", "at": [x, y], "r": r, "hue": hue}
    if text:
        d["text"], d["text_size"] = text, ts
    return d


def A(x0, y0, x1, y1, color="black", curve=0, **kw):
    return {"g": "arrow", "from": [x0, y0], "to": [x1, y1], "color": color, "curve": curve, **kw}


def FAT(x0, y0, x1, y1):
    return {"g": "arrow", "from": [x0, y0], "to": [x1, y1], "style": "fat"}


def INH(x0, y0, x1, y1, curve=0):
    return {"g": "arrow", "from": [x0, y0], "to": [x1, y1], "inhibit": True, "color": "red", "width": 1.2, "curve": curve}


def MEM(y, x0, x1):
    return {"g": "membrane", "y": y, "x0": x0, "x1": x1, "lipids": True}


def TM(x, y, hue, h=3.6, w=0.6):
    """Single transmembrane helix (rounded rod) centred on the membrane."""
    return B(x, y, w, h, hue, shape="box")


def LINE(pts, w=0.12):
    return {"g": "strand", "kind": "gray", "points": pts, "width": w}


def NUC(x0, y0, x1, y1, name="nucleus"):
    return [{"g": "region", "x0": x0, "y0": y0, "x1": x1, "y1": y1, "hue": "blue", "opacity": 0.18, "dashed": True},
            T(x1 - 0.4, y0 + 0.6, name, 0.65, "gray", "right")]


def DNA(x0, x1, y, w=1.0):
    return {"g": "dsdna", "from": [x0, y], "to": [x1, y], "width": w}


def sides(m, x, top="Outside", bot="Cytosol"):
    """Compartment names at the right end of the membrane (x = right edge)."""
    return [T(x, m - 2.0, top, 0.7, "gray", "right"), T(x, m + 2.0, bot, 0.7, "gray", "right")]


def P_(x, y):
    """Phosphate on a protein (small yellow sphere)."""
    return I(x, y, "P", "yellow", 0.3, 0.4)


def ring7(cx, cy, r, hue, rr=0.45):
    """Heptameric ring (apoptosome) seen from above."""
    return [I(cx + r * math.cos(math.radians(360 * k / 7 - 90)), cy + r * math.sin(math.radians(360 * k / 7 - 90)), "", hue, rr)
            for k in range(7)]


# ================================================================ 1. receptor types (Lehninger Fig. 12–2)
M = 5.0
RECEPTORS = {"items": [
    *sides(M, 35.0), MEM(M, -1.0, 35.0),
    # GPCR
    {"g": "ligand", "at": [2.6, 1.0], "shape": "triangle", "hue": "red", "r": 0.45},
    {"g": "gpcr", "at": [2.6, M], "hue": "purple"}, T(2.6, 10.6, "GPCR\n(7TM)"),
    # RTK dimer
    {"g": "ligand", "at": [8.0, 0.6], "shape": "diamond", "hue": "red", "r": 0.45},
    B(7.4, 2.3, 1.0, 1.4, "blue"), B(8.6, 2.3, 1.0, 1.4, "blue"), TM(7.4, M, "blue"), TM(8.6, M, "blue"),
    B(7.2, 7.6, 1.2, 1.3, "blue", "Y", ts=0.5), B(8.8, 7.6, 1.2, 1.3, "blue", "Y", ts=0.5),
    P_(6.4, 8.3), P_(9.6, 8.3), T(8.0, 10.6, "RTK"),
    # receptor Ser/Thr kinase (TGF-β)
    {"g": "ligand", "at": [13.4, 0.6], "shape": "hexagon", "hue": "red", "r": 0.45},
    B(12.8, 2.3, 1.0, 1.4, "teal"), B(14.0, 2.3, 1.0, 1.4, "teal"), TM(12.8, M, "teal"), TM(14.0, M, "teal"),
    B(12.6, 7.6, 1.2, 1.3, "teal", "S/T", ts=0.42), B(14.2, 7.6, 1.2, 1.3, "teal", "S/T", ts=0.42),
    T(13.4, 10.6, "Ser/Thr\nkinase R"),
    # receptor guanylyl cyclase (ANF)
    {"g": "ligand", "at": [18.6, 0.8], "shape": "triangle", "hue": "red", "r": 0.4},
    B(18.6, 2.3, 1.4, 1.4, "green"), TM(18.6, M, "green"), B(18.6, 7.7, 1.8, 1.4, "green", "GC", ts=0.5),
    T(18.6, 10.6, "Guanylyl\ncyclase R"),
    # ligand-gated channel
    {"g": "ligand", "at": [23.6, 1.0], "shape": "diamond", "hue": "red", "r": 0.4},
    {"g": "channel", "at": [24.6, M], "hue": "orange", "open": True},
    I(24.6, 2.0, "Na^+", "purple", 0.4, 0.4), A(24.6, 2.6, 24.6, 7.6, "magenta"),
    T(24.6, 10.6, "Ligand-gated\nchannel"),
    # receptor without kinase → JAK
    {"g": "ligand", "at": [30.0, 0.6], "shape": "hexagon", "hue": "red", "r": 0.45},
    B(29.4, 2.3, 1.0, 1.4, "pink"), B(30.6, 2.3, 1.0, 1.4, "pink"), TM(29.4, M, "pink"), TM(30.6, M, "pink"),
    B(28.6, 7.5, 1.4, 1.1, "orange", "JAK", ts=0.42), B(31.4, 7.5, 1.4, 1.1, "orange", "JAK", ts=0.42),
    T(30.0, 10.6, "Cytokine R\n(→ JAK)"),
    # intracellular receptors
    T(-1.0, 13.0, "Intracellular receptors", 0.8, None, "left"),
    I(3.0, 14.4, "NO", "gray", 0.45, 0.4), A(3.6, 14.4, 5.0, 14.4),
    B(6.6, 14.4, 2.4, 1.4, "green", "sGC", ts=0.5), T(6.6, 16.2, "soluble GC"),
    {"g": "ligand", "at": [11.4, 13.6], "shape": "diamond", "hue": "orange", "r": 0.4},
    B(11.6, 14.9, 1.6, 1.2, "purple", "R", ts=0.5), B(13.0, 15.0, 1.2, 0.9, "gray", "Hsp90", ts=0.36),
    A(14.0, 14.4, 17.4, 14.4),
    *NUC(17.8, 12.6, 34.6, 18.0), DNA(19.0, 33.4, 15.6),
    B(24.0, 14.3, 1.6, 1.3, "purple", "R", ts=0.5), B(25.6, 14.3, 1.6, 1.3, "purple", "R", ts=0.5),
    T(26.2, 17.0, "nuclear receptor dimer on HRE", 0.65),
]}

# ================================================================ 2. GPCR → four G-protein classes
M = 4.0
GCLASSES = {"items": [
    *sides(M, 35.0), MEM(M, -1.0, 35.0),
    {"g": "ligand", "at": [2.4, 0.4], "shape": "triangle", "hue": "red", "r": 0.45},
    {"g": "gpcr", "at": [2.4, M], "hue": "purple"}, T(2.4, 9.0, "GPCR (GEF)"),
    FAT(4.6, 6.6, 7.0, 6.6),
    *[x for i, (gx, name, hue) in enumerate(((10.0, "G_s", "green"), (16.8, "G_i", "orange"), (23.6, "G_q", "blue"),
                                             (30.4, "G_{12/13}", "pink")))
      for x in ({"g": "gprotein", "at": [gx, M + 2.3], "nucleotide": "GTP", "hue": hue}, T(gx - 0.8, M + 4.6, name, 0.85))],
    # effectors
    A(9.2, 9.0, 9.2, 10.6, "green"), B(9.2, 11.6, 2.4, 1.3, "teal", "AC", ts=0.55), T(9.2, 13.4, "cAMP ↑"),
    INH(16.0, 9.0, 16.0, 10.7), B(16.0, 11.6, 2.4, 1.3, "teal", "AC", ts=0.55), T(16.0, 13.4, "cAMP ↓"),
    A(22.8, 9.0, 22.8, 10.6, "green"), B(22.8, 11.6, 2.4, 1.3, "blue", "PLC-β", ts=0.5), T(22.8, 13.4, "IP_3 + DAG"),
    A(29.6, 9.0, 29.6, 10.6, "green"), B(29.6, 11.6, 2.6, 1.3, "pink", "RhoGEF", ts=0.45), T(29.6, 13.4, "Rho"),
    T(9.2, 15.0, "cholera toxin: always on", 0.62, "red"), T(16.0, 15.0, "pertussis toxin: blocked", 0.62, "red"),
]}

# ================================================================ 3. G protein cycle
M = 4.0
GCYCLE = {"items": [
    *sides(M, 35.0), MEM(M, -1.0, 35.0),
    # 1 inactive
    {"g": "gpcr", "at": [2.4, M], "hue": "purple"},
    {"g": "gprotein", "at": [6.0, M + 2.3], "nucleotide": "GDP", "hue": "green"}, T(5.6, 9.6, "1. inactive"),
    FAT(9.0, 6.4, 11.6, 6.4),
    # 2 ligand bound, GDP → GTP
    {"g": "ligand", "at": [14.0, 0.4], "shape": "triangle", "hue": "red", "r": 0.45},
    {"g": "gpcr", "at": [14.0, M], "hue": "purple", "accent": True},
    {"g": "gprotein", "at": [17.6, M + 2.3], "nucleotide": "GTP", "hue": "green", "accent": True},
    T(17.2, 9.6, "2. GDP → GTP"),
    FAT(20.4, 6.4, 22.0, 6.4),
    # 3 α–GTP and βγ act on effectors
    {"g": "gprotein", "at": [29.0, M + 2.3], "nucleotide": "GTP", "hue": "green", "accent": True, "split": True},
    B(23.6, M, 1.8, 3.4, "teal", shape="box"), B(23.6, M + 2.3, 1.8, 1.0, "teal", "AC", ts=0.45),
    A(25.6, M + 2.4, 24.7, M + 2.4, "green"),
    T(27.0, 9.6, "3. α–GTP and βγ act"),
    # 4 GTPase resets
    A(26.0, 11.0, 6.0, 11.0, "black", -0.12),
    T(16.0, 15.0, "4. GTP → GDP + P_i (GTPase, RGS)"),
]}

# ================================================================ 4. GTPase switch (Ras)
GSWITCH = {"items": [
    {"g": "membrane", "y": 0.0, "x0": -1.0, "x1": 23.0, "lipids": True},
    LINE([[4.0, 1.3], [4.0, 2.6]]), LINE([[17.0, 1.3], [17.0, 2.6]]),
    B(4.0, 4.0, 3.6, 2.6, "gray", "Ras", ts=0.7), I(4.0, 5.6, "GDP", "gray", 0.55, 0.42), T(4.0, 7.6, "inactive"),
    B(17.0, 4.0, 3.6, 2.6, "green", "Ras", ts=0.7, accent=True), I(17.0, 5.6, "GTP", "pink", 0.55, 0.42), T(17.0, 7.6, "active"),
    A(6.4, 2.9, 14.6, 2.9, "black", -0.18), T(10.5, 3.4, "GEF (Sos)", 0.75, "green"),
    I(7.6, 4.6, "GDP", "gray", 0.45, 0.38), I(13.4, 4.6, "GTP", "pink", 0.45, 0.38),
    A(14.6, 6.4, 6.4, 6.4, "black", -0.18), T(10.5, 9.4, "GAP (RasGAP, NF1)", 0.75, "red"),
    I(10.5, 7.2, "P_i", "yellow", 0.4, 0.38),
    A(19.2, 4.0, 21.6, 4.0), T(22.0, 4.0, "effectors", 0.75, None, "left"),
    T(10.5, 11.2, "G12V: GTPase lost → always on", 0.7, "red"),
]}

# ================================================================ 5. PIP2: PLC and PI3K
M = 4.0
PIP = {"items": [
    *sides(M, 33.0), MEM(M, -1.0, 33.0),
    {"g": "ligand", "at": [2.4, 0.4], "shape": "triangle", "hue": "red", "r": 0.45},
    {"g": "gpcr", "at": [2.4, M], "hue": "purple"},
    {"g": "gprotein", "at": [6.0, M + 2.3], "nucleotide": "GTP", "hue": "blue", "split": True}, T(5.2, 9.0, "G_q"),
    B(10.6, 6.4, 2.6, 1.4, "blue", "PLC-β", ts=0.5, accent=True),
    {"g": "ligand", "at": [13.4, 5.6], "shape": "hexagon", "hue": "yellow", "r": 0.45}, T(14.3, 7.0, "PIP_2", 0.6, None, "left"),
    A(14.0, 5.4, 15.4, 4.6), T(16.4, 4.0, "DAG", 0.65, "brown"),
    A(13.6, 6.0, 13.6, 9.2), {"g": "ligand", "at": [13.6, 9.9], "shape": "hexagon", "hue": "yellow", "r": 0.45},
    T(14.8, 9.9, "IP_3", 0.65, None, "left"),
    B(18.0, 6.6, 2.0, 1.3, "green", "PKC", ts=0.55, accent=True), A(16.8, 4.9, 17.6, 5.8),
    # ER
    MEM(13.5, -1.0, 16.0), T(-1.0, 15.8, "ER lumen", 0.7, "gray", "left"),
    {"g": "channel", "at": [10.6, 13.5], "hue": "teal", "open": True}, A(13.0, 10.4, 11.6, 11.6),
    I(10.6, 15.6, "Ca^{2+}", "green", 0.45, 0.36), A(10.6, 15.0, 10.6, 11.2, "cyan"), I(10.6, 10.6, "Ca^{2+}", "green", 0.45, 0.36),
    A(11.2, 10.4, 17.2, 7.4, "cyan", -0.15),
    # PI3K branch
    B(22.0, 2.0, 1.0, 1.4, "blue"), B(23.2, 2.0, 1.0, 1.4, "blue"), TM(22.0, M, "blue"), TM(23.2, M, "blue"),
    T(22.6, -0.2, "RTK", 0.65),
    B(25.0, 6.8, 2.2, 1.3, "orange", "PI3K", ts=0.5),
    {"g": "ligand", "at": [27.4, 5.6], "shape": "hexagon", "hue": "orange", "r": 0.45}, T(27.4, 6.9, "PIP_3", 0.6),
    B(26.4, 8.6, 1.8, 1.1, "purple", "PDK1", ts=0.42), B(29.6, 7.2, 2.0, 1.2, "green", "Akt", ts=0.5, accent=True),
    A(27.3, 8.4, 28.6, 7.8), P_(30.6, 7.8),
    T(29.6, 9.6, "survival, GLUT4", 0.65),
]}

# ================================================================ 6. insulin receptor
M = 4.0
INSULIN = {"items": [
    *sides(M, 34.0), MEM(M, -1.0, 34.0),
    I(14.0, -0.4, "Ins", "red", 0.5, 0.4),
    B(13.4, 1.4, 1.2, 1.6, "blue"), B(14.6, 1.4, 1.2, 1.6, "blue"), T(16.4, 1.0, "α_2β_2", 0.65, None, "left"),
    TM(13.4, M, "blue"), TM(14.6, M, "blue"),
    B(13.2, 7.2, 1.3, 1.6, "blue", accent=True), B(14.8, 7.2, 1.3, 1.6, "blue", accent=True), P_(12.3, 7.8), P_(15.7, 7.8),
    B(14.0, 10.4, 2.0, 1.2, "gray", "IRS-1", ts=0.45), P_(15.2, 11.0), A(14.0, 8.2, 14.0, 9.6),
    # PI3K branch (left)
    A(13.0, 10.4, 10.6, 10.4),
    B(9.4, 10.4, 2.2, 1.3, "orange", "PI3K", ts=0.5), A(8.8, 9.8, 7.0, 6.5),
    {"g": "ligand", "at": [6.4, 5.8], "shape": "hexagon", "hue": "orange", "r": 0.42}, T(5.6, 6.4, "PIP_3", 0.6, None, "right"),
    B(4.6, 8.6, 1.8, 1.1, "purple", "PDK1", ts=0.42), A(4.6, 9.2, 4.6, 10.6),
    B(4.6, 11.6, 2.0, 1.3, "green", "Akt", ts=0.5, accent=True), P_(5.8, 12.2),
    T(4.6, 13.6, "GLUT4, glycogen ↑", 0.65),
    # MAPK branch (right)
    A(15.1, 10.4, 19.6, 10.4), B(20.6, 10.4, 1.8, 1.1, "teal", "Grb2", ts=0.42), B(23.0, 10.4, 1.8, 1.1, "pink", "Sos", ts=0.45),
    A(23.4, 9.8, 24.6, 7.4), B(25.4, 6.6, 1.8, 1.3, "green", "Ras", ts=0.45, accent=True), I(26.6, 7.4, "GTP", "pink", 0.35, 0.3),
    A(26.8, 6.6, 28.6, 6.6), B(29.8, 6.6, 1.9, 1.2, "orange", "Raf", ts=0.45),
    A(29.8, 7.3, 29.8, 8.3), B(29.8, 9.0, 1.9, 1.2, "orange", "MEK", ts=0.45),
    A(29.8, 9.7, 29.8, 10.7), B(29.8, 11.4, 1.9, 1.2, "orange", "ERK", ts=0.45, accent=True),
    A(29.8, 12.1, 29.8, 13.4),
    *NUC(21.0, 13.6, 34.0, 18.0), DNA(22.0, 33.0, 16.4), B(28.0, 15.2, 1.8, 1.1, "purple", "Elk-1", ts=0.4),
    B(10.2, 7.4, 2.0, 1.1, "gray", "PLC-γ", ts=0.42), A(12.5, 7.4, 11.3, 7.4),
]}

# ================================================================ 7. JAK–STAT
M = 4.0
JAKSTAT = {"items": [
    *sides(M, 30.0), MEM(M, -1.0, 30.0),
    {"g": "ligand", "at": [6.0, 0.0], "shape": "diamond", "hue": "red", "r": 0.5},
    B(5.4, 1.6, 1.0, 1.4, "pink"), B(6.6, 1.6, 1.0, 1.4, "pink"), TM(5.4, M, "pink"), TM(6.6, M, "pink"),
    LINE([[5.4, 5.8], [5.4, 9.0]], 0.22), LINE([[6.6, 5.8], [6.6, 9.0]], 0.22),
    B(4.4, 7.0, 1.6, 1.2, "orange", "JAK", ts=0.42, accent=True), B(7.6, 7.0, 1.6, 1.2, "orange", "JAK", ts=0.42, accent=True),
    P_(5.4, 8.6), P_(6.6, 8.6),
    B(10.4, 8.6, 2.2, 1.3, "blue", "STAT", ts=0.45), A(9.2, 8.6, 7.2, 8.6), P_(11.6, 9.2),
    FAT(12.4, 9.4, 15.0, 10.6),
    B(16.6, 11.2, 2.0, 1.2, "blue", "STAT", ts=0.42, accent=True), B(16.6, 12.4, 2.0, 1.2, "blue", "STAT", ts=0.42, accent=True), P_(15.4, 11.2), P_(15.4, 12.4),
    A(17.8, 11.8, 19.6, 11.8),
    *NUC(19.8, 8.6, 30.0, 16.2), DNA(20.6, 29.4, 14.4),
    B(24.0, 13.1, 2.0, 1.2, "blue", "STAT", ts=0.42), B(26.0, 13.1, 2.0, 1.2, "blue", "STAT", ts=0.42), P_(23.2, 12.3), P_(26.8, 12.3),
    T(25.0, 10.0, "gene expression", 0.65),
    T(6.0, 11.4, "EPO, leptin, GH, IFN", 0.65, "gray"),
]}

# ================================================================ 8. NF-κB
M = 4.0
NFKB = {"items": [
    *sides(M, 32.0), MEM(M, -1.0, 32.0),
    {"g": "ligand", "at": [4.0, 0.0], "shape": "triangle", "hue": "red", "r": 0.5}, T(5.4, 0.0, "TNF-α", 0.65, None, "left"),
    B(3.4, 1.6, 0.9, 1.4, "purple"), B(4.6, 1.6, 0.9, 1.4, "purple"), TM(3.4, M, "purple"), TM(4.6, M, "purple"),
    B(4.0, 7.4, 2.0, 1.3, "purple", "IKK", ts=0.5, accent=True), A(4.0, 6.0, 4.0, 6.6),
    A(5.2, 7.6, 8.6, 8.6),
    B(10.4, 9.2, 2.6, 1.6, "blue", "NF-κB", ts=0.45), B(11.6, 10.6, 1.6, 0.9, "red", "IκB", ts=0.4), P_(12.6, 11.2),
    A(11.8, 11.4, 12.6, 13.0),
    *[I(12.4 + 0.5 * k, 13.6 - 0.1 * k, "", "yellow", 0.26) for k in range(4)], T(14.8, 13.2, "Ub", 0.6, None, "left"),
    B(13.6, 15.6, 1.4, 2.4, "gray", shape="box"), B(13.6, 15.6, 1.8, 0.5, "gray", shape="box"),
    T(15.0, 15.8, "proteasome", 0.62, "gray", "left"),
    A(12.0, 8.8, 18.6, 8.8),
    *NUC(19.0, 7.0, 31.0, 15.0), DNA(20.0, 30.2, 12.8), B(24.4, 11.2, 2.6, 1.4, "blue", "NF-κB", ts=0.45, accent=True),
    T(25.0, 8.8, "inflammation, anti-apoptosis", 0.62),
]}

# ================================================================ 9. ACh → NO → cGMP
NO = {"items": [
    {"g": "region", "x0": -1.0, "y0": 2.0, "x1": 13.6, "y1": 16.0, "hue": "yellow", "opacity": 0.15},
    {"g": "region", "x0": 15.6, "y0": 2.0, "x1": 33.0, "y1": 16.0, "hue": "pink", "opacity": 0.15},
    T(-0.6, 15.2, "Endothelial cell", 0.75, "gray", "left"), T(32.6, 15.2, "Smooth muscle cell", 0.75, "gray", "right"),
    MEM(2.0, -1.0, 13.6),
    {"g": "ligand", "at": [3.0, -1.2], "shape": "triangle", "hue": "red", "r": 0.45}, T(4.2, -1.2, "ACh", 0.65, None, "left"),
    {"g": "gpcr", "at": [3.0, 2.0], "hue": "purple"},
    {"g": "gprotein", "at": [6.6, 4.3], "nucleotide": "GTP", "hue": "blue"}, T(5.8, 6.4, "M_3 · G_q", 0.6),
    A(7.6, 5.6, 7.6, 7.6), I(7.6, 8.4, "Ca^{2+}", "green", 0.45, 0.36), B(9.6, 8.4, 1.8, 1.0, "yellow", "CaM", ts=0.42),
    A(9.0, 9.0, 7.6, 10.4),
    B(6.6, 11.2, 2.4, 1.4, "orange", "eNOS", ts=0.55, accent=True), T(6.6, 13.2, "Arg → citrulline", 0.6, "gray"),
    I(10.6, 11.2, "NO", "gray", 0.45, 0.4), A(11.2, 11.2, 17.0, 11.2, "gray", dashed=True),
    I(17.6, 11.2, "NO", "gray", 0.45, 0.4),
    B(20.4, 11.2, 2.8, 1.6, "green", "sGC", ts=0.55, accent=True), I(21.6, 12.0, "Fe", "red", 0.3, 0.32),
    T(20.4, 9.4, "GTP → cGMP", 0.6),
    I(24.4, 11.2, "cGMP", "yellow", 0.6, 0.38), A(25.1, 11.2, 26.6, 11.2),
    B(28.4, 11.2, 2.4, 1.4, "teal", "PKG", ts=0.55), T(28.4, 13.2, "Ca^{2+} ↓ → relaxation", 0.62),
    A(24.4, 10.4, 24.4, 7.8, "gray"), I(24.4, 7.2, "GMP", "gray", 0.5, 0.36), T(23.8, 9.0, "PDE5", 0.6, "cyan", "right"),
    INH(28.6, 8.4, 24.7, 9.0), T(28.8, 8.2, "sildenafil", 0.62, "red", "left"),
    T(17.6, 6.0, "nitroglycerin", 0.62, "red"), A(17.6, 6.8, 17.6, 10.4, "red"),
]}

# ================================================================ 10. two-component system
M = 4.0
TWOCOMP = {"items": [
    *sides(M, 24.0, "Periplasm", "Cytoplasm"), MEM(M, -1.0, 24.0),
    I(6.0, -0.4, "", "red", 0.45), T(7.0, -0.4, "signal", 0.65, None, "left"),
    B(5.4, 1.4, 1.0, 1.4, "teal"), B(6.6, 1.4, 1.0, 1.4, "teal"), TM(5.4, M, "teal"), TM(6.6, M, "teal"),
    B(5.4, 7.4, 1.3, 2.0, "teal", accent=True), B(6.6, 7.4, 1.3, 2.0, "teal", accent=True),
    I(7.8, 7.0, "His~P", "yellow", 0.62, 0.32), T(3.4, 7.4, "sensor\nkinase", 0.62, None, "right"),
    A(8.5, 7.6, 11.4, 8.6), B(12.8, 8.8, 2.2, 1.3, "orange", "RR", ts=0.5), I(14.4, 9.6, "Asp~P", "yellow", 0.62, 0.32),
    A(13.2, 9.8, 15.0, 12.0), DNA(13.0, 23.0, 13.2), B(17.4, 12.0, 1.8, 1.1, "orange", "RR", ts=0.45, accent=True),
    T(20.6, 11.6, "transcription", 0.65, None, "left"),
    A(14.8, 8.6, 18.6, 7.2), T(19.0, 7.0, "flagellar motor", 0.65, None, "left"),
]}

# ================================================================ 11. ADP-ribosylating toxins
M = 3.0
TOXINS = {"items": [
    MEM(M, -1.0, 21.0), T(20.8, M - 2.0, "Outside", 0.7, "gray", "right"), T(20.8, M + 2.0, "Cytosol", 0.7, "gray", "right"),
    # cholera
    T(4.0, -0.2, "Cholera toxin", 0.8, "red"),
    {"g": "gprotein", "at": [4.0, M + 2.3], "nucleotide": "GTP", "hue": "green", "accent": True},
    {"g": "ligand", "at": [3.2, 7.4], "shape": "hexagon", "hue": "red", "r": 0.36}, T(3.2, 8.4, "ADP-ribose (Arg)", 0.55, "red"),
    T(4.0, 9.8, "G_sα: GTPase lost", 0.65), A(4.0, 10.6, 4.0, 11.6),
    T(4.0, 12.4, "cAMP ↑ → Cl^- out\n(watery diarrhea)", 0.65),
    # pertussis
    T(14.6, -0.2, "Pertussis toxin", 0.8, "red"),
    {"g": "gprotein", "at": [14.6, M + 2.3], "nucleotide": "GDP", "hue": "orange"},
    {"g": "ligand", "at": [13.8, 7.4], "shape": "hexagon", "hue": "red", "r": 0.36}, T(13.8, 8.4, "ADP-ribose (Cys)", 0.55, "red"),
    T(14.6, 9.8, "G_iα stuck: no inhibition", 0.65), A(14.6, 10.6, 14.6, 11.6), T(14.6, 12.4, "cAMP ↑", 0.65),
    # diphtheria
    T(28.0, -0.2, "Diphtheria toxin", 0.8, "red"),
    {"g": "strand", "kind": "rna", "points": [[23.0, 6.0], [33.0, 6.0]], "ends": ["5'", "3'"]},
    {"g": "ribosome", "at": [28.0, 6.0], "w": 4.4},
    B(31.2, 3.0, 1.8, 1.2, "teal", "eEF2", ts=0.45),
    {"g": "ligand", "at": [32.4, 2.0], "shape": "hexagon", "hue": "red", "r": 0.36},
    T(28.0, 9.8, "eEF2 (diphthamide) blocked", 0.65), A(28.0, 10.6, 28.0, 11.6), T(28.0, 12.4, "translation stops", 0.65),
    T(16.0, 14.6, "NAD^+ → ADP-ribose–protein + nicotinamide", 0.7, "gray"),
]}

# ================================================================ 12. apoptosis
M = 3.0
APOP = {"items": [
    MEM(M, -1.0, 12.0), T(-1.0, 0.0, "Extrinsic", 0.8, None, "left"),
    {"g": "ligand", "at": [4.0, -0.8], "shape": "triangle", "hue": "red", "r": 0.45}, T(5.2, -0.8, "FasL, TNF", 0.62, None, "left"),
    B(3.4, 0.6, 0.8, 1.2, "purple"), B(4.6, 0.6, 0.8, 1.2, "purple"), TM(3.4, M, "purple"), TM(4.6, M, "purple"),
    B(4.0, 6.0, 2.6, 1.3, "orange", "caspase 8", ts=0.45, accent=True), A(4.0, 4.9, 4.0, 5.3),
    # mitochondrion
    T(15.0, 0.0, "Intrinsic", 0.8, None, "left"),
    {"g": "region", "x0": 15.0, "y0": 1.4, "x1": 29.0, "y1": 7.0, "hue": "orange", "opacity": 0.25},
    T(28.6, 2.0, "mitochondrion", 0.62, "gray", "right"),
    T(12.6, 10.2, "DNA damage, ROS", 0.62), A(14.6, 9.8, 18.4, 8.4),
    {"g": "channel", "at": [20.0, 7.0], "hue": "purple", "open": True, "h": 2.2}, T(20.0, 9.4, "Bax/Bak", 0.62),
    B(24.6, 8.2, 1.8, 1.1, "green", "Bcl-2", ts=0.42), INH(23.6, 8.0, 21.4, 7.6),
    I(20.0, 4.8, "c", "red", 0.35, 0.4), A(20.0, 5.4, 20.0, 10.4, "red"), I(20.0, 11.0, "c", "red", 0.35, 0.4),
    A(6.0, 6.2, 15.8, 6.2, "gray", 0.12), T(10.6, 7.6, "tBid", 0.6, "gray"),
    # apoptosome
    *ring7(20.0, 14.2, 1.3, "purple"), I(20.0, 14.2, "c", "red", 0.35, 0.36), T(22.0, 14.2, "apoptosome", 0.62, None, "left"),
    T(22.0, 15.2, "(Apaf-1 · cyt c · dATP)_7", 0.55, "gray", "left"),
    A(20.0, 11.6, 20.0, 12.4),
    B(14.6, 15.4, 2.6, 1.3, "orange", "caspase 9", ts=0.45, accent=True), A(18.4, 14.8, 16.0, 15.2),
    # executioner
    B(6.6, 13.4, 3.0, 1.4, "red", "caspase 3, 7", ts=0.45, accent=True),
    A(4.0, 6.8, 6.0, 12.6), A(13.2, 15.2, 8.2, 13.8),
    T(6.6, 16.0, "shrinks · DNA ladder\nPS out · no inflammation", 0.6),
]}

# ================================================================ 13. LDL receptor-mediated endocytosis
M = 3.0
LDL = {"items": [
    *sides(M, 31.0, "Outside", "Cytosol"), MEM(M, -1.0, 31.0),
    I(3.0, -1.6, "LDL", "yellow", 0.9, 0.42), B(3.0, 0.2, 1.4, 1.0, "blue"), TM(3.0, M, "blue", 3.0, 0.5),
    T(3.0, 5.6, "LDL receptor", 0.6),
    # coated vesicle
    I(10.6, 9.4, "", "blue", 1.5),
    *[I(10.6 + 1.85 * math.cos(math.radians(a)), 9.4 + 1.85 * math.sin(math.radians(a)), "", "gray", 0.26) for a in range(0, 360, 30)],
    I(10.6, 9.4, "LDL", "yellow", 0.75, 0.38), T(10.6, 12.2, "clathrin-coated vesicle", 0.6), T(13.0, 6.6, "dynamin", 0.6, "gray", "left"),
    A(4.4, 5.6, 8.8, 7.8),
    # endosome
    A(12.6, 9.6, 15.4, 9.6),
    I(18.0, 9.6, "", "teal", 2.2), I(17.4, 9.8, "LDL", "yellow", 0.7, 0.36), T(20.6, 10.4, "endosome (pH ↓)", 0.6, None, "left"),
    A(19.6, 7.8, 23.6, 4.6, "blue", -0.1), T(23.8, 6.2, "receptor recycles", 0.6, "blue", "left"),
    # lysosome
    A(17.4, 11.8, 17.4, 13.6), I(17.4, 15.0, "", "pink", 1.3), T(19.2, 15.0, "lysosome", 0.62, None, "left"),
    T(17.4, 17.2, "cholesterol inhibits HMG-CoA reductase", 0.62),
]}

# ================================================================ 14. membrane asymmetry
M = 6.0
OUT_Y, IN_Y = M - 1.06, M + 1.06
ASYM = {"items": [
    T(-1.0, M - 3.4, "Outer leaflet", 0.75, "gray", "left"), T(-1.0, M + 3.4, "Inner leaflet (cytosol)", 0.75, "gray", "left"),
    MEM(M, -1.0, 31.0),
    *[I(x, OUT_Y, t, h, 0.34, 0.3) for x, t, h in ((1.0, "PC", "blue"), (2.4, "SM", "teal"), (4.6, "PC", "blue"), (12.6, "SM", "teal"),
                                                   (19.4, "PC", "blue"), (21.4, "GL", "green"))],
    *[I(x, IN_Y, t, h, 0.34, 0.3) for x, t, h in ((1.6, "PE", "orange"), (3.4, "PS", "red"), (5.4, "PI", "yellow"), (12.0, "PE", "orange"),
                                                  (19.0, "PS", "red"), (20.6, "PE", "orange"))],
    # flippase: PS out → in (ATP)
    {"g": "pump", "at": [8.6, M], "hue": "orange"}, T(8.6, 1.2, "flippase", 0.65),
    I(7.2, 2.6, "PS", "red", 0.34, 0.3), A(7.2, 3.2, 7.2, 8.6, "magenta"),
    # floppase: in → out (ABC)
    {"g": "transporter", "at": [15.6, M], "hue": "pink", "open": "in"}, T(15.6, 1.2, "floppase (ABC)", 0.65),
    A(17.2, 8.8, 17.2, 3.4, "cyan"),
    # scramblase: both ways, Ca2+
    {"g": "channel", "at": [24.6, M], "hue": "teal", "open": True}, T(24.6, 1.2, "scramblase", 0.65),
    A(23.6, 3.6, 23.6, 8.6, "gray"), A(25.6, 8.6, 25.6, 3.6, "gray"), I(26.6, 9.6, "Ca^{2+}", "green", 0.4, 0.32),
    # anchored proteins
    B(28.6, 2.6, 2.0, 1.3, "purple", "GPI", ts=0.45), LINE([[28.6, 3.3], [28.6, OUT_Y]]),
    B(29.0, 10.0, 1.8, 1.2, "green", "Ras", ts=0.42), LINE([[29.0, 9.4], [29.0, IN_Y]]),
    T(8.6, 11.8, "ATP", 0.6, "gray"), T(15.6, 11.8, "ATP", 0.6, "gray"), T(24.6, 11.8, "no ATP", 0.6, "gray"),
]}

# ================================================================ 15. electron transport chain
M = 7.0
ETC = {"items": [
    T(35.8, 3.4, "IMS (H^+ high)", 0.7, "gray", "right"), T(35.8, 14.6, "Matrix", 0.7, "gray", "right"),
    MEM(M, -1.0, 36.0),
    # complex I (L-shaped)
    B(3.0, M, 2.8, 3.2, "blue", shape="box"), B(2.4, M + 3.0, 1.8, 3.2, "blue", "I", ts=0.8, shape="box"),
    T(2.4, 12.6, "NADH → NAD^+", 0.62),
    # complex II
    B(8.0, M + 2.2, 2.0, 2.2, "teal", "II", ts=0.7), T(8.0, 12.6, "succinate → fumarate", 0.62),
    # Q
    I(11.0, M, "Q", "yellow", 0.55, 0.5),
    # complex III
    B(15.0, M, 3.2, 3.4, "purple", "III", ts=0.8, shape="box"),
    # cyt c
    I(18.2, M - 2.4, "c", "red", 0.5, 0.5),
    # complex IV
    B(21.4, M, 3.0, 3.4, "green", "IV", ts=0.8, shape="box"), T(21.4, 12.6, "½O_2 → H_2O", 0.62),
    # ATP synthase
    B(29.0, M, 3.0, 3.0, "orange", "F_o", ts=0.6, shape="box"), B(29.0, M + 2.6, 0.7, 2.0, "orange", shape="box"),
    B(29.0, M + 4.4, 3.4, 2.2, "orange", "F_1", ts=0.6, accent=True),
    T(31.0, M + 4.4, "ADP + P_i → ATP", 0.62, None, "left"),
    # electrons
    A(3.4, M + 1.6, 10.4, M + 0.4, E_COL), A(8.6, M + 1.2, 10.4, M + 0.5, E_COL), A(11.6, M, 13.4, M, E_COL),
    A(16.0, M - 1.9, 17.6, M - 2.4, E_COL), A(18.8, M - 2.4, 20.6, M - 1.9, E_COL),
    # protons
    *[x for xh, n in ((3.6, "4H^+"), (15.0, "4H^+"), (21.4, "2H^+"))
      for x in (A(xh + 2.0, M + 1.0, xh + 2.0, M - 3.0, H_COL), T(xh + 2.0, M - 3.6, n, 0.62, H_COL))],
    A(31.2, M - 3.0, 31.2, M + 1.2, H_COL), T(31.2, M - 3.6, "H^+", 0.62, H_COL),
    # inhibitors
    T(3.0, M - 5.0, "rotenone", 0.58, "red"), T(9.4, M + 3.6, "malonate", 0.58, "red", "left"), T(15.0, M - 5.0, "antimycin A", 0.58, "red"),
    T(21.4, M - 5.0, "CN^-, CO", 0.58, "red"), T(27.0, M - 5.0, "oligomycin", 0.58, "red"),
]}

# ================================================================ 16. ATP synthase
M = 4.0
ATPS = {"items": [
    T(-1.0, 0.4, "IMS (H^+ high)", 0.75, "gray", "left"), T(-1.0, 16.4, "Matrix", 0.75, "gray", "left"),
    MEM(M, -1.0, 31.0),
    # F_o: a + c ring
    B(5.0, M, 1.6, 3.0, "orange", "a", ts=0.6, shape="box"),
    *[B(6.6 + 0.62 * k, M, 0.56, 2.8, "yellow", shape="box") for k in range(8)], T(8.8, M - 2.0, "c ring (8 c)", 0.6),
    # b2 peripheral stalk and δ
    B(4.4, 9.8, 0.6, 8.8, "orange", shape="box"), B(6.6, 14.2, 4.8, 0.9, "orange", "δ", ts=0.5),
    # central stalk γε
    B(8.8, M + 3.6, 0.9, 4.6, "purple", "γ", ts=0.5, shape="box"), B(9.8, M + 1.9, 0.8, 0.7, "purple", "ε", ts=0.42),
    # F_1 head α3β3 (side view)
    B(7.0, 11.8, 1.8, 3.6, "pink", "α", ts=0.55), B(8.8, 11.8, 1.8, 3.6, "green", "β", ts=0.55, accent=True),
    B(10.6, 11.8, 1.8, 3.6, "pink", "α", ts=0.55),
    # proton path
    I(6.0, 0.2, "H^+", "pink", 0.4, 0.36), A(6.0, 0.8, 6.0, 2.4, H_COL),
    A(11.6, 5.2, 12.8, 7.6, H_COL), I(13.2, 8.2, "H^+", "pink", 0.4, 0.36),
    A(7.0, M + 1.9, 11.2, M + 1.9, "gray", -0.3), T(8.0, 7.4, "rotor", 0.6, "gray", "right"),
    T(12.6, 11.2, "ADP + P_i → ATP", 0.65, None, "left"),
    T(-1.0, 7.0, "oligomycin", 0.6, "red", "left"),
    # binding change: O → L → T
    T(26.0, 8.4, "binding change (β)", 0.7),
    I(26.0, 10.4, "O", "gray", 0.8, 0.6), I(29.0, 14.6, "L", "yellow", 0.8, 0.6), I(23.0, 14.6, "T", "green", 0.8, 0.6),
    A(27.0, 10.9, 28.6, 13.4, "cyan", -0.2), A(28.0, 15.2, 24.0, 15.2, "cyan", -0.2), A(23.4, 13.6, 25.0, 10.9, "cyan", -0.2),
    T(30.0, 12.0, "ADP + P_i bind", 0.55, None, "left"), T(22.0, 12.0, "ATP forms", 0.55, None, "right"),
    T(26.0, 17.0, "ATP released (needs energy)", 0.6),
    # ANT / PiC
    {"g": "transporter", "at": [17.0, M], "hue": "blue"}, T(17.0, M - 3.0, "ANT", 0.6),
    A(16.2, 6.2, 16.2, 1.8, "black"), A(17.8, 1.8, 17.8, 6.2, "gray"), T(17.0, 7.0, "ATP out · ADP in", 0.55),
    {"g": "transporter", "at": [23.0, M], "hue": "teal"}, T(23.0, M - 3.0, "PiC", 0.6), A(23.0, 1.6, 23.0, 6.2, "gray"),
    T(23.0, 7.0, "P_i + H^+ in", 0.55),
]}

# ================================================================ 17. malate–aspartate shuttle
M = 7.0
MAL = {"items": [
    T(-1.0, 0.8, "Cytosol", 0.75, "gray", "left"), T(-1.0, 13.4, "Matrix", 0.75, "gray", "left"),
    MEM(M, -1.0, 27.0), T(26.8, M + 2.0, "inner membrane", 0.55, "gray", "right"),
    I(4.0, 2.6, "Asp", "purple", 0.7, 0.42), I(12.0, 2.6, "OAA", "orange", 0.7, 0.42), I(20.0, 2.6, "Mal", "green", 0.7, 0.42),
    I(4.0, 11.4, "Asp", "purple", 0.7, 0.42), I(12.0, 11.4, "OAA", "orange", 0.7, 0.42), I(20.0, 11.4, "Mal", "green", 0.7, 0.42),
    A(12.8, 2.6, 19.2, 2.6), T(16.0, 1.4, "NADH → NAD^+", 0.6, "cyan"), T(16.0, 3.6, "MDH", 0.55, "gray"),
    A(4.8, 2.6, 11.2, 2.6), T(8.0, 1.4, "AST", 0.55, "gray"),
    {"g": "transporter", "at": [20.0, M], "hue": "green"}, A(20.0, 3.4, 20.0, 10.6, "black"),
    I(22.0, 11.4, "αKG", "yellow", 0.55, 0.36), A(22.0, 10.8, 22.0, 3.4, "gray"), I(22.0, 2.6, "αKG", "yellow", 0.55, 0.36),
    A(19.2, 11.4, 12.8, 11.4), T(16.0, 12.6, "NAD^+ → NADH", 0.6, "cyan"), T(16.0, 10.4, "MDH", 0.55, "gray"),
    A(11.2, 11.4, 4.8, 11.4), T(8.0, 12.6, "AST", 0.55, "gray"),
    {"g": "transporter", "at": [4.0, M], "hue": "purple"}, A(4.0, 10.6, 4.0, 3.4, "black"),
    I(2.0, 2.6, "Glu", "blue", 0.55, 0.36), A(2.0, 3.2, 2.0, 10.8, "gray"), I(2.0, 11.4, "Glu", "blue", 0.55, 0.36),
    T(16.0, 14.6, "NADH → complex I (2.5 ATP)", 0.62),
]}

# ================================================================ 18. glycerol 3-phosphate shuttle
M = 7.0
G3P = {"items": [
    T(-1.0, 0.6, "Cytosol / IMS", 0.75, "gray", "left"), T(-1.0, 12.4, "Matrix", 0.75, "gray", "left"),
    MEM(M, -1.0, 27.0),
    I(4.0, 2.6, "DHAP", "orange", 0.8, 0.38), I(14.0, 2.6, "G3P", "yellow", 0.8, 0.38),
    A(4.9, 2.2, 13.1, 2.2, "black", -0.1), T(9.0, 0.6, "NADH → NAD^+", 0.6, "cyan"),
    B(15.4, 5.0, 2.8, 1.3, "teal", "mGPDH", ts=0.45, accent=True), I(17.2, 5.0, "FAD", "yellow", 0.45, 0.32),
    A(13.6, 3.4, 13.0, 4.4, "black"), A(13.0, 4.4, 4.6, 3.4, "black", -0.1),
    I(19.6, M, "Q", "yellow", 0.55, 0.5), A(17.8, 5.6, 19.0, 6.6, E_COL),
    B(24.0, M, 3.0, 3.4, "purple", "III", ts=0.8, shape="box"), A(20.2, M, 22.4, M, E_COL),
    T(13.0, 11.0, "bypasses complex I → 1.5 ATP", 0.62),
]}

# ================================================================ 19. membrane protein classes (new)
M = 5.0
MPCLASS = {"items": [
    *sides(M, 33.0), MEM(M, -1.0, 33.0),
    # single-pass
    B(2.0, 2.6, 1.8, 1.4, "blue"), TM(2.0, M, "blue"), B(2.0, 7.4, 1.6, 1.2, "blue"), T(2.0, 10.4, "single-pass"),
    # multipass
    {"g": "gpcr", "at": [7.6, M], "hue": "purple"}, T(7.6, 10.4, "multipass"),
    # β barrel / channel
    {"g": "channel", "at": [13.0, M], "hue": "teal", "open": True}, T(13.0, 10.4, "channel / barrel"),
    T(7.6, 12.0, "integral: detergent to remove", 0.6, "gray"),
    # peripheral
    B(18.6, M, 1.0, 3.2, "orange", shape="box"), B(18.6, 7.1, 2.2, 1.2, "yellow", "P", ts=0.5),
    T(20.6, 10.4, "peripheral"), T(20.6, 12.0, "salt / pH releases", 0.6, "gray"),
    # amphitropic on headgroups
    B(22.6, 6.9, 2.0, 1.2, "yellow", "P", ts=0.5),
    # lipid-anchored
    B(26.0, 1.4, 2.0, 1.3, "green", "GPI", ts=0.42), LINE([[26.0, 2.1], [26.0, M - 1.06]]),
    B(28.8, 8.6, 1.8, 1.2, "pink", "Ras", ts=0.42), LINE([[28.8, 8.0], [28.8, M + 1.06]]),
    T(27.4, 10.4, "lipid-anchored"), T(27.4, 12.0, "GPI: PI-PLC · prenyl, acyl", 0.6, "gray"),
]}

# ================================================================ 20. membrane permeability (new)
M = 6.0
PERM = {"items": [
    *sides(M, 31.0), MEM(M, -1.0, 31.0),
    # gases / hydrophobic: free
    *[x for k, t in enumerate(("O_2", "CO_2", "N_2")) for x in (I(1.4 + 2.2 * k, 1.6, t, "gray", 0.55, 0.36),)],
    A(3.6, 2.4, 3.6, 10.4, "black", width=1.8), T(3.6, 12.2, "gases, hydrophobic"), T(3.6, 13.4, "free", 0.62, "green"),
    # small uncharged polar: slow
    *[I(9.6 + 2.2 * k, 1.6, t, "blue", 0.6, 0.32) for k, t in enumerate(("H_2O", "urea", "glyc"))],
    A(11.8, 2.4, 11.8, 10.4, "black", dashed=True), T(11.8, 12.2, "small polar"), T(11.8, 13.4, "slight", 0.62, "orange"),
    # large polar: blocked
    I(19.0, 1.6, "Glc", "orange", 0.65, 0.4), A(19.0, 2.4, 19.0, M - 1.6), {"g": "mark", "kind": "x", "at": [19.0, M - 2.4]},
    T(19.0, 12.2, "large polar"), T(19.0, 13.4, "no", 0.62, "red"),
    # ions and charged
    *[I(24.0 + 1.8 * k, 1.6, t, h, 0.5, 0.34) for k, (t, h) in enumerate((("Na^+", "purple"), ("K^+", "green"), ("Cl^-", "teal"),
                                                                           ("ATP", "pink")))],
    A(26.6, 2.4, 26.6, M - 1.6), {"g": "mark", "kind": "x", "at": [26.6, M - 2.4]},
    T(26.6, 12.2, "ions, charged"), T(26.6, 13.4, "no", 0.62, "red"),
]}

ILLUS = {
    "receptor_types": RECEPTORS, "gpcr_g_protein_classes": GCLASSES, "g_protein_cycle": GCYCLE, "gtpase_switch": GSWITCH,
    "phosphoinositide_signaling": PIP, "insulin_rtk_signaling": INSULIN, "jak_stat_pathway": JAKSTAT, "nfkb_pathway": NFKB,
    "no_cgmp_relaxation": NO, "two_component_system": TWOCOMP, "bacterial_toxins": TOXINS, "apoptosis_pathways": APOP,
    "ldl_receptor_endocytosis": LDL, "membrane_asymmetry": ASYM, "electron_transport_chain": ETC, "atp_synthase_fof1": ATPS,
    "malate_aspartate_shuttle": MAL, "glycerol_phosphate_shuttle": G3P,
}

NEW = [
    {"id": "membrane_protein_classes", "title": "Classes of membrane proteins", "layout": "scene", "scene": MPCLASS,
     "tags": ["membrane", "membrane_lipid"],
     "notes_zh": "Integral（內在）：穿過疏水核心，單次跨膜（α 螺旋約 20 個疏水胺基酸）、多次跨膜（GPCR 7 次）、β-barrel（porin、外膜）；"
                 "需界面活性劑（SDS、Triton X-100）才能取出。Peripheral（周邊）：以離子鍵、氫鍵結合在膜表面或其他膜蛋白上（如 cytochrome c、"
                 "spectrin），高鹽、改變 pH 或 EDTA 即可洗下。Amphitropic：可逆地結合在膜上。Lipid-anchored：GPI（膜外側，PI-PLC 可切下）；"
                 "內側的 myristoyl（N 端 Gly）、palmitoyl（Cys）、prenyl（farnesyl／geranylgeranyl，C 端 CaaX，如 Ras）。"
                 "Hydropathy plot 找出約 20 個疏水胺基酸的區段可預測跨膜螺旋。",
     "parts_png": False, "status": "pending"},
    {"id": "membrane_permeability", "title": "Permeability of the lipid bilayer", "layout": "scene", "scene": PERM,
     "tags": ["membrane", "membrane_transport"],
     "notes_zh": "脂雙層的通透性取決於大小與極性：氣體（O₂、CO₂、N₂、NO）與疏水小分子（steroid、乙醇）自由擴散；小的不帶電極性分子"
                 "（H₂O、urea、glycerol）緩慢通過（水另經 aquaporin）；較大的極性分子（glucose、fructose）幾乎不能通過，需 transporter；"
                 "所有離子（Na⁺、K⁺、Ca²⁺、Cl⁻、H⁺）與帶電分子（胺基酸、ATP、核苷酸）不能通過，需 channel、transporter 或 pump。",
     "parts_png": False, "status": "pending"},
]


# ================================================================ writer

def dump(path, d):
    if path.exists() and "status" in d:
        old = (yaml.safe_load(path.read_text(encoding="utf-8")) or {}).get("status")
        if old in ("approved", "changes_requested"):
            d["status"] = old
    path.write_text(yaml.safe_dump(d, allow_unicode=True, sort_keys=False, width=200), encoding="utf-8")


def with_parent(tags):
    parent = {"membrane_transport": "membrane", "membrane_lipid": "membrane"}
    out = []
    for x in tags:
        if x in parent:
            out.append(parent[x])
        out.append(x)
    return list(dict.fromkeys(out))


def main(data):
    pdir = Path(data) / "pathways"
    for pid, scene in ILLUS.items():
        path = pdir / f"{pid}.yaml"
        d = yaml.safe_load(path.read_text(encoding="utf-8"))
        d.setdefault("variants", {})["illus"] = {"note_zh": "圖片版（插畫）", "scene": scene}
        path.write_text(yaml.safe_dump(d, allow_unicode=True, sort_keys=False, width=200), encoding="utf-8")
    for p in NEW:
        d = tidy({**p, "tags": with_parent(p["tags"])})
        dump(pdir / f"{d['id']}.yaml", d)
    print(f"{len(ILLUS)} 個圖片版、{len(NEW)} 個新插畫")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data")
