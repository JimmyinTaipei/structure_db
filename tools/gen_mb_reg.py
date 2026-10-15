"""Round 22: molecular biology V — gene regulation (operons, attenuation, bacterial RNA regulators, eukaryotic activation,
DNA-binding motifs, chromatin / epigenetics, RNAi, nuclear receptors; DNA methylation and lac-inducer structures).
Figures carry names and short labels only; explanations live in notes_zh.

    python tools/gen_mb_reg.py data
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from gen_mb_repair import GAP, down, dump, duplex, table, txt, with_parent  # noqa: E402
from gen_mb_tx import ball, box, hairpin, prot, rna, rungs  # noqa: E402
from gen_mb_tl import ribo, tr  # noqa: E402
from figure_text import tidy  # noqa: E402
from reaction_colors import COLORS, hl  # noqa: E402

M, C = "magenta", "cyan"
ON, OFF = "#2E9E4F", "#D62828"


def region(x0, x1, y, hue, text="", op=0.85, ts=0.62):
    """Coloured stretch of a duplex (y = top strand) with its name between the strands."""
    out = [{"g": "region", "x0": x0, "y0": y - 0.2, "x1": x1, "y1": y + GAP + 0.2, "hue": hue, "opacity": op, "outline": False}]
    if text:
        out.append(txt((x0 + x1) / 2, y + GAP / 2, text, ts))
    return out


def on(x0, x1, y, weak=False):
    return {"g": "arrow", "from": [x0, y], "to": [x1, y], "color": ON, "width": 1.0 if weak else 1.8, **({"dashed": True} if weak else {})}


def off(x, y):
    return {"g": "mark", "kind": "x", "at": [x, y]}


def inhibit(x0, y0, x1, y1, curve=0):
    return {"g": "arrow", "from": [x0, y0], "to": [x1, y1], "inhibit": True, "color": OFF, "width": 1.2, "curve": curve}


def head(x, y, text, size=0.9):
    return txt(x, y, text, size, None, "left")


# ================================================================ control modes (notes p.74 right; lecture 4-216/217)
CX, CY = [0.0, 15.0], [2.5, 10.5]


def cell(x, y, kind, bound, sig, result):
    out = [*duplex(y, x, x + 12, top="parent", bot="parent", ends=False),
           *region(x + 1.0, x + 2.6, y, "gray", "P"), *region(x + 2.6, x + 4.2, y, "pink" if kind == "R" else "green", "O" if kind == "R" else "A"),
           *region(x + 4.8, x + 11.4, y, "blue", "genes", 0.5)]
    hue_ = "purple" if kind == "R" else "dkgreen"
    px, py = (x + 3.4, y - 0.75) if bound else (x + 3.4, y - 2.9)
    out.append(prot(px, py, 2.0, 1.3, hue_, kind, "kidney", 0.7))
    if sig:
        out.append(ball(px + 1.25, py - 0.55, "", "orange" if sig == "ind" else "red", 0.35))
    out.append(on(x + 5.0, x + 11.0, y - 1.2) if result else off(x + 8.0, y - 1.3))
    return out


CONTROL = {"items": [
    txt(CX[0] + 6.0, -2.2, "Negative control", 0.95), txt(CX[1] + 6.0, -2.2, "Positive control", 0.95),
    txt(-1.0, CY[0] + 0.5, "Induction", 0.9, None, "right"), txt(-1.0, CY[1] + 0.5, "Repression", 0.9, None, "right"),
    *cell(CX[0], CY[0], "R", False, "ind", True), txt(CX[0] + 6.0, CY[0] + 2.4, "lac: LacI + allolactose", 0.7, "gray"),
    *cell(CX[1], CY[0], "A", True, "ind", True), txt(CX[1] + 6.0, CY[0] + 2.4, "ara: AraC + arabinose", 0.7, "gray"),
    *cell(CX[0], CY[1], "R", True, "co", False), txt(CX[0] + 6.0, CY[1] + 2.4, "trp: TrpR + Trp", 0.7, "gray"),
    *cell(CX[1], CY[1], "A", False, "co", False), txt(CX[1] + 6.0, CY[1] + 2.4, "(rare)", 0.7, "gray"),
    prot(1.0, 15.2, 1.4, 0.9, "purple", "R", "kidney", 0.55), txt(2.0, 15.2, "repressor", 0.7, None, "left"),
    prot(6.6, 15.2, 1.4, 0.9, "dkgreen", "A", "kidney", 0.55), txt(7.6, 15.2, "activator", 0.7, None, "left"),
    ball(12.0, 15.2, "", "orange", 0.33), txt(12.6, 15.2, "inducer", 0.7, None, "left"),
    ball(16.6, 15.2, "", "red", 0.33), txt(17.2, 15.2, "corepressor", 0.7, None, "left"),
    on(21.4, 23.0, 15.2), txt(23.3, 15.2, "on", 0.7, None, "left"), off(25.4, 15.2), txt(26.0, 15.2, "off", 0.7, None, "left"),
]}

# ================================================================ lac operon (notes p.75; lecture 4-218, Lehninger 28-18)
LY = [6.2, 10.4, 14.6, 18.8]


def lac_row(y, glc, lac):
    cap = not glc
    rep = not lac
    out = [*duplex(y, 4.0, 22.0, top="parent", bot="parent", ends=False),
           *region(5.0, 6.4, y, "green", "CAP"), *region(6.4, 8.0, y, "gray", "P"), *region(8.0, 9.2, y, "pink", "O"),
           *region(9.8, 21.6, y, "blue", "lacZ  lacY  lacA", 0.45)]
    if lac:
        out += [prot(7.6, y - 0.95, 2.6, 1.5, "blue", "pol", "oval", 0.55)]
    if cap:
        out += [prot(5.7, y - 0.8, 1.6, 1.2, "dkgreen", "CAP", "kidney", 0.5), ball(5.0, y - 1.6, "", "yellow", 0.26), ball(6.4, y - 1.6, "", "yellow", 0.26)]
    if rep:
        out.append(prot(8.6, y - 0.75, 1.9, 1.3, "purple", "LacI", "kidney", 0.55))
    else:
        out += [prot(17.6, y - 1.75, 1.9, 1.3, "purple", "LacI", "kidney", 0.55), ball(18.7, y - 2.25, "", "orange", 0.3)]
    if lac and cap:
        out.append(on(10.0, 21.0, y - 0.55))
    elif lac:
        out.append(on(10.0, 14.0, y - 0.55, weak=True))
    else:
        out.append(off(15.0, y - 1.0))
    out.append(txt(2.8, y + 0.5, ("glucose +" if glc else "glucose −") + "\n" + ("lactose +" if lac else "lactose −"), 0.68, None, "right"))
    out.append(txt(23.0, y + 0.5, "high" if (lac and cap) else "low" if lac else "off", 0.8, ON if lac else OFF, "left"))
    return out


LAC = {"items": [
    *duplex(1.0, 0.0, 26.0, top="parent", bot="parent", ends=True),
    *region(0.6, 3.6, 1.0, "purple", "lacI", 0.6), *region(5.0, 6.4, 1.0, "green", "CAP"), *region(6.4, 8.0, 1.0, "gray", "P"),
    *region(8.0, 9.2, 1.0, "pink", "O_1"), *region(9.8, 16.6, 1.0, "blue", "lacZ"), *region(16.6, 20.6, 1.0, "teal", "lacY"),
    *region(20.6, 24.6, 1.0, "orange", "lacA"),
    txt(13.2, 3.0, "β-galactosidase", 0.62, "gray"), txt(18.6, 3.0, "permease", 0.62, "gray"), txt(22.6, 3.0, "transacetylase", 0.62, "gray"),
    *lac_row(LY[0], True, False), *lac_row(LY[1], False, False), *lac_row(LY[2], True, True), *lac_row(LY[3], False, True),
    ball(5.0, 21.6, "", "yellow", 0.26), txt(5.5, 21.6, "cAMP", 0.68, None, "left"),
    ball(9.0, 21.6, "", "orange", 0.3), txt(9.5, 21.6, "allolactose", 0.68, None, "left"),
]}

# ================================================================ trp attenuation (notes p.76; lecture Lehninger 28-20)
AY = 7.5


def attenuation(x, low):
    out = [*duplex(AY + 2.2, x, x + 13.0, top="parent", bot="parent", ends=False)]
    if low:
        hp = hairpin(x + 6.6, AY, 2.0, 0.45, 0.42)
        out += [rna([(x + 0.4, AY), (x + 5.2, AY), (x + 6.39, AY)] + hp + [(x + 6.81, AY), (x + 9.4, AY), (x + 10.6, AY + 1.9)], ("5'", "")),
                *rungs(x + 6.6, AY, 2.0, 0.42),
                ribo(x + 3.0, AY, 3.6, False), txt(x + 0.9, AY + 0.65, "1", 0.7, "#B5651F"),
                txt(x + 5.4, AY - 1.2, "2", 0.7, "#B5651F"), txt(x + 7.8, AY - 1.2, "3", 0.7, "#B5651F"), txt(x + 9.0, AY + 0.6, "4", 0.7, "#B5651F"),
                prot(x + 11.0, AY + 2.7, 2.6, 1.8, "blue", "pol", "oval", 0.55)]
    else:
        hp = hairpin(x + 7.8, AY, 2.0, 0.45, 0.42)
        out += [rna([(x + 0.4, AY), (x + 6.4, AY), (x + 7.59, AY)] + hp + [(x + 8.01, AY), (x + 9.4, AY)], ("5'", "3'")),
                *rungs(x + 7.8, AY, 2.0, 0.42),
                ribo(x + 4.6, AY, 3.6, False), txt(x + 6.6, AY - 1.2, "3", 0.7, "#B5651F"), txt(x + 9.0, AY - 1.2, "4", 0.7, "#B5651F"),
                txt(x + 9.0, AY + 0.65, "UUUU", 0.6, "#C2338A"),
                prot(x + 11.6, AY + 0.3, 2.6, 1.8, "blue", "pol", "oval", 0.55)]
    return out


TRP = {"items": [
    *duplex(1.0, 0.0, 27.0, top="parent", bot="parent", ends=True),
    *region(0.6, 2.2, 1.0, "gray", "P O"), *region(2.2, 5.0, 1.0, "yellow", "trpL"),
    *[x for k, g in enumerate("EDCBA") for x in region(5.4 + 4.2 * k, 9.4 + 4.2 * k, 1.0, "teal", "trp" + g, 0.6)],
    txt(3.6, 3.0, "leader + attenuator", 0.62, "gray"),
    txt(0.0, 4.6, "Trp low", 0.9, None, "left"),
    txt(16.0, 4.6, "Trp high", 0.9, None, "left"),
    *attenuation(0.0, True), *attenuation(16.0, False),
    txt(6.6, AY + 4.4, "2:3 pair → read-through", 0.75, ON),
    txt(22.6, AY + 4.4, "3:4 terminator → stop", 0.75, OFF),
]}

# ================================================================ stringent response + r-protein feedback (notes p.77)
SY = 5.0
STRINGENT = {"items": [
    rna([(0.0, SY), (10.0, SY)], ("5'", "3'")), ribo(5.0, SY, 7.0, True),
    tr(5.0, SY, "", pep=3), tr(5.0 + 1.89, SY, "", hue="#D9D9D9"),
    prot(7.6, SY - 4.4, 2.2, 1.3, "orange", "RelA", "kidney", 0.55),
    txt(5.0 + 1.89, SY + 2.9, "uncharged tRNA", 0.65, "gray"),
    {"g": "arrow", "from": [9.4, SY - 4.6], "to": [13.2, SY - 4.8], "color": "black", "width": 1.1},
    txt(11.6, SY - 5.7, "GTP + ATP", 0.62, "gray"),
    ball(14.4, SY - 4.8, "ppGpp", "yellow", 0.75, 0.45),
    inhibit(15.2, SY - 4.0, 18.6, SY - 1.6, 0.2),
    *duplex(SY, 15.0, 27.0, top="parent", bot="parent", ends=False),
    *region(16.0, 17.4, SY, "gray", "P"), *region(17.8, 26.4, SY, "green", "rRNA, tRNA genes", 0.5),
    prot(19.4, SY - 0.9, 2.6, 1.6, "blue", "pol", "oval", 0.55), off(23.0, SY - 1.4),
    head(0.0, 10.4, "r-protein as translational repressor"),
    rna([(0.0, 13.0), (12.0, 13.0)], ("5'", "3'")), box(1.2, 2.6, 13.0, "yellow", "SD", 0.75),
    # L11 operon = rplK (L11) + rplA (L1) only; L10 / L7/L12 form a separate unit regulated by L10
    box(2.8, 6.2, 13.0, "teal", "L11"), box(6.2, 9.6, 13.0, "teal", "L1"),
    prot(1.9, 12.0, 1.6, 1.2, "brown", "L1", "round", 0.55), off(4.4, 11.7), off(7.9, 11.7),
]}

# ================================================================ sRNA / Hfq and riboswitch (notes p.77, p.79)
RY = 4.0
RIBOSW = {"items": [
    head(0.0, -1.0, "sRNA + Hfq (trans)"),
    rna([(0.0, RY), (11.0, RY)], ("5'", "3'")), box(3.0, 4.6, RY, "yellow", "SD", 0.75), box(5.4, 7.0, RY, "gray", "AUG", 0.75, 0.5),
    rna([(2.6, RY - 0.75), (5.2, RY - 0.75)], ("3'", "5'"), kind="#E0A100"),
    {"g": "ring", "kind": "hexamer", "at": [3.9, RY - 2.0], "r": 0.55, "hue": "teal", "text": "Hfq"},
    prot(9.0, RY - 2.6, 2.6, 1.5, "#C9BEDF", "30S", "oval", 0.55), off(6.8, RY - 2.6),
    txt(3.9, RY + 1.3, "SD blocked", 0.7, OFF),
    head(14.0, -1.0, "Riboswitch (cis): TPP"),
    rna([(14.0, RY), (16.4, RY), (16.4, RY - 2.2), (17.4, RY - 2.8), (18.4, RY - 2.2), (18.4, RY),
         (20.6, RY), (20.6, RY - 2.0), (21.1, RY - 2.5), (21.6, RY - 2.0), (21.6, RY), (25.0, RY)], ("5'", "3'")),
    {"g": "ligand", "at": [17.4, RY - 1.3], "shape": "diamond", "hue": "orange", "r": 0.42},
    {"g": "blob", "at": [21.6, RY - 1.0], "w": 1.4, "h": 0.55, "shape": "box", "hue": "yellow", "rot": 90, "text": ""}, txt(22.4, RY - 1.0, "SD", 0.6, "#B08F10", "left"),
    txt(17.4, RY + 1.0, "aptamer", 0.62, "gray"), 
    txt(17.4, RY - 3.4, "TPP", 0.62, "#B5651F"),
    txt(19.5, RY + 2.2, "metabolite bound → translation off", 0.7, OFF),
]}

# ================================================================ Salmonella phase variation (notes p.79; lecture 4-215)
PV = [2.5, 9.5]


def phase(y, forward):
    out = [*duplex(y, 0.0, 26.0, top="parent", bot="parent", ends=False),
           {"g": "ligand", "at": [1.2, y + 0.5], "shape": "triangle", "hue": "purple", "r": 0.35},
           {"g": "ligand", "at": [8.2, y + 0.5], "shape": "triangle", "hue": "purple", "r": 0.35},
           *region(1.8, 6.0, y, "orange", "hin"), *region(6.2, 7.6, y, "gray", "P"),
           *region(8.8, 13.0, y, "teal", "fljB"), *region(13.0, 16.4, y, "pink", "fljA"), *region(20.0, 25.0, y, "green", "fliC")]
    if forward:
        out += [on(7.0, 16.2, y - 0.9), txt(11.0, y + 2.2, "FljB flagellin (H2)", 0.68, ON), txt(22.5, y + 2.2, "repressed by FljA", 0.68, OFF),
                inhibit(16.2, y - 1.8, 21.0, y - 0.8, -0.2)]
    else:
        out += [on(6.8, 1.6, y - 0.9), off(11.0, y - 1.0), on(20.2, 24.8, y - 0.9), txt(22.5, y + 2.2, "FliC flagellin (H1)", 0.68, ON)]
    return out


PHASE = {"items": [
    *phase(PV[0], True),
    {"g": "arrow", "from": [4.6, PV[0] + 2.3], "to": [4.6, PV[1] - 1.8], "color": "gray", "width": 1.0}, txt(5.2, (PV[0] + PV[1]) / 2 + 0.6, "Hin inverts (hix)", 0.7, None, "left"),
    *phase(PV[1], False),
]}

# ================================================================ eukaryotic activation (notes p.73, p.80; lecture 4-211)
DNA = [(0.0, 12.0), (6.0, 12.0), (8.6, 11.6), (9.4, 8.2), (10.2, 5.4), (12.6, 4.6), (15.0, 5.4), (15.8, 8.2), (16.6, 11.6), (19.0, 12.0), (30.0, 12.0)]
EUK = {"items": [
    {"g": "strand", "kind": "parent", "points": [list(p) for p in DNA], "width": 0.42},
    {"g": "region", "x0": 11.4, "y0": 4.2, "x1": 13.8, "y1": 5.0, "hue": "yellow", "opacity": 0.9, "outline": False},
    prot(12.6, 3.0, 2.4, 1.6, "dkgreen", "activator", "kidney", 0.5),
    prot(17.8, 6.8, 6.4, 4.6, "teal", "Mediator", "oval", 0.7),
    prot(20.8, 10.8, 4.6, 3.2, "blue", "Pol II", "oval", 0.7, accent=True),
    prot(17.6, 11.8, 1.8, 1.2, "green", "TFIID", "oval", 0.45),
    {"g": "nucleosome", "at": [2.2, 11.6], "r": 0.9}, {"g": "nucleosome", "at": [26.8, 11.6], "r": 0.9},
    ball(1.4, 10.4, "Ac", "yellow", 0.3, 0.4), ball(2.3, 10.2, "Ac", "yellow", 0.3, 0.4),
    prot(3.9, 10.1, 1.9, 1.4, "orange", "HAT", "round", 0.55),
    prot(26.8, 10.0, 2.8, 1.3, "purple", "SWI/SNF", "oval", 0.5),
    txt(12.6, 5.7, "enhancer", 0.62, "gray"), txt(17.6, 13.2, "TATA", 0.62, "gray"),
    txt(3.9, 8.8, "(SAGA)", 0.6, "gray"), txt(26.8, 8.8, "remodeler (ATP)", 0.6, "gray"),
]}

# ================================================================ DNA-binding motifs (lecture 4-211; notes p.72)
MX, MYY = [0.0, 14.0], [4.0, 12.5]


def motif_dna(x, y):
    return {"g": "helix", "from": [x, y], "to": [x + 11, y], "pitch": 3.6, "amp": 0.7}


def helix_bar(x, y, rot, hue="purple", h=2.2):
    return {"g": "blob", "at": [x, y], "w": 0.75, "h": h, "shape": "box", "hue": hue, "rot": rot}


MOTIFS = {"items": [
    motif_dna(MX[0], MYY[0]), helix_bar(MX[0] + 4.4, MYY[0] - 2.2, 50), helix_bar(MX[0] + 5.6, MYY[0] - 0.6, 90, "pink"),
    txt(MX[0] + 5.5, MYY[0] - 4.3, "Helix-turn-helix", 0.85), txt(MX[0] + 5.5, MYY[0] + 1.8, "CAP, LacI, Cro, homeodomain", 0.65, "gray"),
    motif_dna(MX[1], MYY[0]), ball(MX[1] + 5.0, MYY[0] - 2.0, "Zn", "gray", 0.42, 0.45),
    {"g": "arrow", "from": [MX[1] + 3.4, MYY[0] - 0.9], "to": [MX[1] + 4.0, MYY[0] - 2.6], "color": "#5F5F5F", "width": 1.6},
    rna([(MX[1] + 4.0, MYY[0] - 1.8), (MX[1] + 4.6, MYY[0] - 2.0)], kind="gray", w=0.08),
    rna([(MX[1] + 5.3, MYY[0] - 1.7), (MX[1] + 5.6, MYY[0] - 0.9)], kind="gray", w=0.08),
    helix_bar(MX[1] + 5.8, MYY[0] - 0.6, 90, "pink"),
    txt(MX[1] + 5.5, MYY[0] - 4.3, "Zinc finger", 0.85), txt(MX[1] + 5.5, MYY[0] + 1.8, "steroid / thyroid receptors, TFIIIA, Sp1", 0.65, "gray"),
    motif_dna(MX[0], MYY[1]), helix_bar(MX[0] + 4.6, MYY[1] - 2.0, 22, "pink", 4.2), helix_bar(MX[0] + 6.4, MYY[1] - 2.0, -22, "teal", 4.2),
    txt(MX[0] + 5.5, MYY[1] - 4.8, "Leucine zipper (bZIP)", 0.85), txt(MX[0] + 5.5, MYY[1] + 1.8, "c-Fos/c-Jun (AP-1), CREB, C/EBP", 0.65, "gray"),
    motif_dna(MX[1], MYY[1]), helix_bar(MX[1] + 4.4, MYY[1] - 1.4, 22, "pink", 3.0), helix_bar(MX[1] + 6.6, MYY[1] - 1.4, -22, "teal", 3.0),
    helix_bar(MX[1] + 4.6, MYY[1] - 3.4, -50, "pink", 1.8), helix_bar(MX[1] + 6.4, MYY[1] - 3.4, 50, "teal", 1.8),
    txt(MX[1] + 5.5, MYY[1] - 4.8, "Helix-loop-helix (bHLH)", 0.85), txt(MX[1] + 5.5, MYY[1] + 1.8, "MyoD, Myc/Max", 0.65, "gray"),
    txt(12.5, 16.0, "pink = recognition helix in the major groove", 0.65, "gray"),
]}

# ================================================================ chromatin and epigenetics (notes p.71; lecture 4-214)
EY = 4.0
EPI = {"items": [
    head(0.0, -0.8, "Euchromatin: on"), head(16.0, -0.8, "Heterochromatin: off"),
    {"g": "strand", "kind": "parent", "points": [[0.0, EY], [13.0, EY]], "width": 0.3},
    *[x for k in range(3) for x in ({"g": "nucleosome", "at": [2.0 + 4.4 * k, EY - 0.3], "r": 0.9},
                                     ball(1.6 + 4.4 * k, EY - 1.7, "Ac", "yellow", 0.3, 0.4), ball(2.5 + 4.4 * k, EY - 1.8, "Ac", "yellow", 0.3, 0.4))],
    on(4.0, 11.0, EY + 1.6),
    {"g": "strand", "kind": "parent", "points": [[16.0, EY], [26.0, EY]], "width": 0.3},
    *[x for k in range(5) for x in ({"g": "nucleosome", "at": [17.0 + 2.0 * k, EY - 0.3], "r": 0.9},
                                     ball(17.0 + 2.0 * k, EY - 1.8, "Me", "purple", 0.3, 0.38))],
    *[{"g": "mark", "kind": "methyl", "at": [17.9 + 2.0 * k, EY + 0.15], "side": "down", "color": M} for k in range(4)],
    off(27.2, EY + 0.2),
    {"g": "arrow", "from": [13.6, EY - 1.4], "to": [15.6, EY - 1.4], "color": "black", "width": 1.1},
    {"g": "arrow", "from": [15.6, EY + 0.6], "to": [13.6, EY + 0.6], "color": "black", "width": 1.1},
    txt(14.6, EY - 2.3, "HDAC, HMT, DNMT", 0.62), txt(14.6, EY + 1.5, "HAT", 0.62),
    txt(6.5, EY + 2.6, "H3/H4 Lys acetylated · CpG unmethylated", 0.65, "gray"),
    txt(21.0, EY + 2.6, "H3K9 / K27 methylated · CpG 5mC", 0.65, "gray"),
    head(0.0, 9.6, "Genomic imprinting"),
    # Igf2: the paternal ICR (insulator) is methylated → CTCF cannot bind → the enhancer reaches IGF2 (Alberts 7e, p. 438)
    *duplex(11.4, 0.0, 12.0, top="parent", bot="parent", ends=False), *region(0.6, 4.4, 11.4, "blue", "IGF2"),
    *region(5.6, 7.4, 11.4, "gray", "ICR"), *region(9.4, 11.4, 11.4, "yellow", "enh"), on(0.8, 4.2, 10.4),
    *[{"g": "mark", "kind": "methyl", "at": [x, 11.4], "side": "up", "color": M} for x in (6.5,)],
    txt(-0.4, 11.9, "paternal", 0.7, None, "right"),
    *duplex(14.6, 0.0, 12.0, top="parent", bot="parent", ends=False), *region(0.6, 4.4, 14.6, "blue", "IGF2"),
    *region(5.6, 7.4, 14.6, "gray", "ICR"), *region(9.4, 11.4, 14.6, "yellow", "enh"), off(2.5, 13.6),
    prot(6.5, 13.7, 1.6, 0.8, "teal", "CTCF", ts=0.5),
    txt(-0.4, 15.1, "maternal", 0.7, None, "right"),
    txt(15.0, 12.0, "Prader–Willi: paternal 15q11 lost", 0.7, None, "left"),
    txt(15.0, 13.4, "Angelman: maternal 15q11 (UBE3A) lost", 0.7, None, "left"),
]}

# ================================================================ miRNA / siRNA (notes p.80; lecture 4-218)


def hp_h(x0, x1, y, ends=("5'", "3'"), kind="rna", gap=0.5):
    """Horizontal hairpin: top strand left → right, loop on the right, back along y + gap."""
    r = gap / 2
    loop = [(x1 + r * math.cos(math.radians(t)), y + r + r * math.sin(math.radians(t))) for t in range(-90, 91, 20)]
    return [rna([(x0, y), (x1, y)] + loop + [(x1, y + gap), (x0, y + gap)], ends, kind),
            *[rna([(x, y + 0.12), (x, y + gap - 0.12)], kind="gray", w=0.06) for x in [x0 + 0.4 + 0.45 * k for k in range(int((x1 - x0 - 0.5) / 0.45))]]]


def duplex_rna(x0, y, w=2.2, guide="#C2338A", passenger="#E7A9C9"):
    return [rna([(x0, y), (x0 + w, y)], kind=guide), rna([(x0 + 0.15, y + 0.4), (x0 + w + 0.15, y + 0.4)], kind=passenger)]


QY = [3.0, 12.0]
RNAI = {"items": [
    {"g": "region", "x0": -0.6, "y0": -1.6, "x1": 9.6, "y1": 6.6, "hue": "blue", "opacity": 0.12, "outline": True, "dashed": True},
    txt(0.0, -1.0, "nucleus", 0.65, "gray", "left"),
    head(10.4, -1.0, "miRNA"), head(0.0, QY[1] - 2.4, "siRNA"),
    *hp_h(0.2, 3.2, QY[0], kind="rna"), txt(1.7, QY[0] - 0.9, "pri-miRNA", 0.62, "gray"),
    prot(4.2, QY[0] + 0.25, 1.6, 1.4, "orange", "Drosha", "kidney", 0.42),
    {"g": "arrow", "from": [5.2, QY[0] + 0.25], "to": [6.4, QY[0] + 0.25], "color": "black", "width": 1.0},
    *hp_h(6.6, 8.4, QY[0], ("", "")), txt(7.5, QY[0] - 0.9, "pre-miRNA", 0.62, "gray"),
    prot(9.6, QY[0] + 0.25, 1.6, 1.2, "teal", "Exp5", "oval", 0.42),
    {"g": "arrow", "from": [10.6, QY[0] + 0.25], "to": [11.8, QY[0] + 0.25], "color": "black", "width": 1.0},
    prot(13.0, QY[0] + 0.25, 1.9, 1.4, "purple", "Dicer", "kidney", 0.45),
    {"g": "arrow", "from": [14.2, QY[0] + 0.25], "to": [15.2, QY[0] + 0.25], "color": "black", "width": 1.0},
    *duplex_rna(15.4, QY[0]), txt(16.5, QY[0] - 0.8, "21–23 nt", 0.6, "gray"),
    {"g": "arrow", "from": [18.0, QY[0] + 0.25], "to": [19.0, QY[0] + 0.25], "color": "black", "width": 1.0},
    prot(20.6, QY[0] + 0.1, 2.6, 1.8, "green", "RISC (Ago)", "oval", 0.45),
    rna([(19.0, QY[0] + 2.6), (30.0, QY[0] + 2.6)], ("5'", "3'")), box(26.0, 29.0, QY[0] + 2.6, "gray", "3′UTR", 0.7, 0.5),
    rna([(26.4, QY[0] + 2.0), (28.4, QY[0] + 2.0)], kind="#C2338A"),
    *[rna([(x, QY[0] + 2.1), (x, QY[0] + 2.5)], kind="gray", w=0.06) for x in (26.7, 27.1, 27.5)],
    txt(27.5, QY[0] + 3.6, "partial pairing (seed)", 0.62, "gray"),
    txt(24.0, QY[0] + 5.0, "translation ↓, deadenylation", 0.7, OFF),
    {"g": "arrow", "from": [21.6, QY[0] + 1.0], "to": [25.8, QY[0] + 1.9], "color": "gray", "width": 1.0, "curve": -0.2},
    *duplex_rna(0.2, QY[1], 6.0, "#3F4F6B", "#8A8A8A"), txt(3.2, QY[1] - 0.8, "long dsRNA (virus, experiment)", 0.62, "gray"),
    {"g": "arrow", "from": [7.0, QY[1] + 0.25], "to": [11.8, QY[1] + 0.25], "color": "black", "width": 1.0},
    prot(13.0, QY[1] + 0.25, 1.9, 1.4, "purple", "Dicer", "kidney", 0.45),
    {"g": "arrow", "from": [14.2, QY[1] + 0.25], "to": [15.2, QY[1] + 0.25], "color": "black", "width": 1.0},
    *duplex_rna(15.4, QY[1]),
    {"g": "arrow", "from": [18.0, QY[1] + 0.25], "to": [19.0, QY[1] + 0.25], "color": "black", "width": 1.0},
    prot(20.6, QY[1] + 0.1, 2.6, 1.8, "green", "RISC (Ago2)", "oval", 0.45),
    rna([(19.0, QY[1] + 2.6), (30.0, QY[1] + 2.6)], ("5'", "3'")),
    rna([(23.4, QY[1] + 2.0), (25.8, QY[1] + 2.0)], kind="#C2338A"),
    *[rna([(x, QY[1] + 2.1), (x, QY[1] + 2.5)], kind="gray", w=0.06) for x in (23.7, 24.1, 24.5, 24.9, 25.3)],
    {"g": "mark", "kind": "cut", "at": [24.6, QY[1] + 2.6], "side": "down"},
    txt(24.6, QY[1] + 4.0, "perfect pairing → mRNA cleaved", 0.7, OFF),
    {"g": "arrow", "from": [21.6, QY[1] + 1.0], "to": [23.2, QY[1] + 1.9], "color": "gray", "width": 1.0, "curve": -0.2},
]}

# ================================================================ nuclear (steroid) receptors (lecture 4-211)
NR = {"items": [
    {"g": "membrane", "y": 1.0, "x0": 0.0, "x1": 28.0, "t": 1.6},
    {"g": "ligand", "at": [3.0, -1.4], "shape": "diamond", "hue": "orange", "r": 0.45},
    {"g": "arrow", "from": [3.0, -0.6], "to": [3.0, 3.6], "color": "gray", "width": 1.0},
    txt(4.0, -1.4, "steroid", 0.7, None, "left"),
    prot(6.0, 6.0, 2.6, 1.8, "purple", "GR", "kidney", 0.6), prot(7.9, 7.2, 2.0, 1.4, "gray", "HSP90", "oval", 0.45),
    {"g": "arrow", "from": [9.4, 6.0], "to": [12.2, 6.0], "style": "fat"},
    prot(14.4, 5.6, 2.6, 1.8, "purple", "GR", "kidney", 0.6), {"g": "ligand", "at": [14.8, 4.4], "shape": "diamond", "hue": "orange", "r": 0.38},
    prot(15.6, 8.6, 2.0, 1.4, "gray", "HSP90", "oval", 0.45),
    {"g": "region", "x0": 18.6, "y0": 3.0, "x1": 29.0, "y1": 14.6, "hue": "blue", "opacity": 0.12, "dashed": True},
    txt(28.6, 3.6, "nucleus", 0.65, "gray", "right"),
    {"g": "arrow", "from": [16.6, 5.6], "to": [19.8, 7.6], "color": "gray", "width": 1.0},
    *duplex(11.4, 19.2, 28.4, top="parent", bot="parent", ends=False), *region(20.4, 23.6, 11.4, "yellow", "HRE"),
    prot(21.0, 10.4, 1.9, 1.5, "purple", "", "kidney"), prot(23.0, 10.4, 1.9, 1.5, "purple", "", "kidney"),
    {"g": "ligand", "at": [21.0, 9.4], "shape": "diamond", "hue": "orange", "r": 0.32}, {"g": "ligand", "at": [23.0, 9.4], "shape": "diamond", "hue": "orange", "r": 0.32},
    on(24.2, 28.0, 10.0),
]}

# ================================================================ small RNA table (lecture 4-232; notes p.80)
SRNA_TABLE = table([
    ("miRNA", "own genes (Pol II) · hairpin → Drosha → Dicer · 21–23 nt · Ago · partial pairing (seed 2–8) in 3′UTR → translation ↓, decay"),
    ("siRNA", "long dsRNA (virus, transposon, experiment) · Dicer · 21–23 nt · Ago2 · perfect pairing → mRNA cleaved (RNAi)"),
    ("piRNA", "single-stranded precursors, no Dicer · 24–30 nt · PIWI · silences transposons in the germ line"),
    ("sRNA (bacteria)", "50–300 nt · Hfq chaperone · pairs with SD / mRNA → translation ↓ or ↑, decay (OxyS, DsrA)"),
    ("Riboswitch", "part of the mRNA 5′UTR · aptamer binds a metabolite (TPP, FMN, SAM, lysine) → terminator or SD hidden"),
], w_name=6.2, x_desc=6.8, hue="green")

# ================================================================ molecules and reactions
MOLECULES = [
    {"id": "ppgpp", "name": "Guanosine tetraphosphate (ppGpp)", "label": "ppGpp", "aliases": ["(p)ppGpp", "alarmone", "magic spot"],
     "smiles": "Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](OP(=O)([O-])OP(=O)([O-])[O-])[C@H]1O",
     "ring": {"coordgen": True}, "tags": ["gene_regulation"],
     "notes": "Alarmone：guanosine 5′-二磷酸 3′-二磷酸（pppGpp 在 5′ 端多一個磷酸）。胺基酸缺乏時，核糖體 A 位出現未帶胺基酸的 tRNA → RelA 以 ATP 的焦磷酸接到 GTP/GDP 的 3′-OH。"
              "結合 RNA polymerase → rRNA、tRNA 轉錄下降（stringent response）。"},
    {"id": "iptg", "name": "Isopropyl β-D-1-thiogalactopyranoside (IPTG)", "label": "IPTG", "aliases": ["isopropylthiogalactoside"],
     "smiles": "CC(C)S[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O",
     "haworth": {"ring": "pyranose", "numbers": True,
                 "subs": {1: {"up": ["S–CH(CH_3)_2", 0], "down": "H"}, 2: {"up": "H", "down": "OH"}, 3: {"up": "OH", "down": "H"},
                          4: {"up": "HO", "down": "H"}, 5: {"up": "CH_2OH", "down": "H"}}},
     "tags": ["gene_regulation"],
     "notes": "Gratuitous inducer：與 LacI 結合使其離開 operator，但 C1 是硫醚（S-糖苷鍵），β-galactosidase 不能水解 → 濃度不會下降；實驗與藍白篩選常用。"},
    {"id": "allolactose", "name": "Allolactose", "label": "Allolactose", "aliases": ["Gal β(1→6) Glc"],
     "smiles": "OC[C@H]1O[C@@H](OC[C@H]2O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O", "ring": {"coordgen": True},
     "tags": ["gene_regulation"],
     "notes": "Lac operon 的天然 inducer：β-galactosidase 的轉糖苷副反應把 lactose 的 β(1→4) 改成 β(1→6)。結合 LacI 四聚體 → 構形改變、離開 operator。"},
    {"id": "cytosine_dr", "name": "Cytosine (in DNA)", "label": "Cytosine (in DNA)", "aliases": [],
     "smiles": "[17*]n1ccc(N)nc1=O", "ring": {"coordgen": True}, "tags": ["gene_regulation"],
     "notes": "DNA 中的 cytosine（dR = 去氧核糖與骨架）。"},
    {"id": "methylcytosine_dr", "name": "5-Methylcytosine (in DNA)", "label": "5-Methylcytosine (in DNA)", "aliases": ["5mC"],
     "smiles": "[17*]n1cc(C)c(N)nc1=O", "ring": {"coordgen": True}, "tags": ["gene_regulation"],
     "notes": "CpG 中 C 的 C5 甲基化（位於大溝，不影響配對）。啟動子 CpG island 甲基化 → 基因沉默（招募 MBD 蛋白與 HDAC）。自發去胺 → T（CpG 突變熱點）。"},
]

REACTIONS = [
    {"id": "dna_methyltransferase", "title": "CpG methylation (DNA methyltransferase)", "substrates": ["cytosine_dr"], "products": ["methylcytosine_dr"],
     "cofactors_in": ["SAM"], "cofactors_out": ["SAH"], "extra_out": ["H+"], "enzyme": "DNA methyltransferase\n(DNMT1, DNMT3A/B)", "tags": ["gene_regulation"],
     "card": {"type": "SAM 的甲基（洋紅）轉到 cytosine 的 C5；酵素 Cys 先暫時加成到 C6 活化 C5，H⁺ 從 C5 離開",
              "exam_points": ["DNMT3A/B：de novo 甲基化（胚胎、imprinting）；DNMT1：維持型，辨識半甲基化 CpG，複製後照樣甲基化",
                              "去甲基化：TET 氧化 5mC → 5hmC → … → BER", "5-azacytidine / decitabine 抑制 DNMT（MDS 治療）"]}},
]

COLORS.update({
    "dna_methyltransferase": {"SAM": {"colors": {0: M}}, "methylcytosine_dr": hl("[CH3:1]c")},
})

PATHWAYS = [
    {"id": "control_modes", "title": "Modes of transcriptional control", "layout": "scene", "scene": CONTROL, "tags": ["gene_regulation"],
     "notes": "依調節蛋白結合 DNA 後的效果：結合即關 = negative control（repressor），結合即開 = positive control（activator）。"
              "依小分子：存在使基因開 = induction（inducer），存在使基因關 = repression（corepressor）。"
              "Negative inducible：lac（allolactose 使 LacI 離開）；negative repressible：trp（Trp 使 TrpR 結合）；positive inducible：CAP–cAMP、AraC–arabinose。"
              "沒有調節 = constitutive（如 lac O^c 突變、管家基因）。"},
    {"id": "lac_operon", "title": "The lac operon: glucose and lactose", "layout": "scene", "scene": LAC, "tags": ["gene_regulation"],
     "notes": "lacI（repressor，四聚體）、CAP site、P、O、lacZ（β-galactosidase）、lacY（permease）、lacA（transacetylase）。"
              "Lactose → allolactose（β1→6）結合 LacI → 離開 operator（O1；O2、O3 形成 DNA 環加強抑制）。"
              "葡萄糖低 → cAMP 高 → CAP–cAMP 結合 CAP site、彎曲 DNA、招募 RNA pol（catabolite repression / positive control）。"
              "只有「葡萄糖低 + 乳糖有」時大量表現。突變：lacI⁻、O^c 為組成性表現；I^s（不受 inducer 影響）恆關。IPTG 為 gratuitous inducer。"},
    {"id": "trp_attenuation", "title": "trp operon: attenuation", "layout": "scene", "scene": TRP, "tags": ["gene_regulation"],
     "notes": "trp operon 有兩層控制：TrpR + Trp（corepressor）結合 operator（negative repressible，約 70 倍）；attenuation（約 10 倍）。"
              "Leader（trpL）14 aa 的胜肽含兩個連續 Trp 密碼子，mRNA 有序列 1–4。Trp 少：核糖體停在序列 1 → 2:3 配對（antiterminator）→ 繼續轉錄。"
              "Trp 多：核糖體蓋住序列 2 → 3:4 形成終止莖環 + UUUU → 提前終止。需要轉錄與轉譯偶聯（只有原核）。his、phe、leu operon 也用 attenuation。"},
    {"id": "stringent_response", "title": "Stringent response and r-protein feedback", "layout": "scene", "scene": STRINGENT, "tags": ["gene_regulation"],
     "notes": "胺基酸缺乏 → 未帶胺基酸的 tRNA 進入 A 位 → 核糖體上的 RelA（stringent factor）以 GTP + ATP 合成 (p)ppGpp（alarmone）→ 結合 RNA pol，rRNA、tRNA 轉錄下降。"
              "核糖體蛋白的轉譯回饋：r-protein（如 L1）優先結合 rRNA；rRNA 不足時多出的 L1 結合自己 operon 的 mRNA（SD 附近）→ 抑制整個 operon（L11 operon = L11 + L1）的轉譯（translational repressor）；L10 operon（L10、L7/L12）則由 L10 自我調控。"},
    {"id": "bacterial_riboswitch_srna", "title": "Bacterial sRNA and riboswitch", "layout": "scene", "scene": RIBOSW, "tags": ["gene_regulation"],
     "notes": "sRNA（50–300 nt，反式）：Hfq 協助與 mRNA 配對，蓋住 SD → 抑制轉譯或促使分解（OxyS 抑制 rpoS；DsrA、RprA 打開 rpoS 的髮夾反而促進）。"
              "Riboswitch（順式，mRNA 5′UTR 本身）：aptamer 結合代謝物（TPP、FMN、SAM、lysine）→ expression platform 改變構形：形成終止子（轉錄）或把 SD 藏進髮夾（轉譯）；"
              "真菌、植物的 TPP riboswitch 也可改變剪接。"},
    {"id": "phase_variation", "title": "Salmonella phase variation (Hin inversion)", "layout": "scene", "scene": PHASE, "tags": ["gene_regulation", "recombination"],
     "notes": "Hin recombinase 在兩個 hix 反向重複間倒轉一段 DNA（含 hin 與 fljB 的啟動子）。順向：fljB（H2 flagellin）與 fljA（fliC 的 repressor）表現 → fliC 關。"
              "倒轉：fljB/fljA 不表現 → fliC（H1 flagellin）開。約每 1000 代切換一次，逃避宿主免疫。其他：phage Mu Gin、酵母 mating-type（HO）、錐蟲 VSG（基因轉換）。"},
    {"id": "eukaryotic_activation", "title": "Eukaryotic transcription activation", "layout": "scene", "scene": EUK, "tags": ["gene_regulation", "transcription"],
     "notes": "Activator 的 DNA-binding domain 結合 enhancer（可遠在數千 bp 外、方向不拘），activation domain（acidic：Gal4；Pro-rich：CTF1；Gln-rich：Sp1）"
              "經 DNA 彎成環，透過 coactivators 招募 Pol II：Mediator（> 20 次單元）、HAT（SAGA，組蛋白乙醯化）、chromatin remodeling（SWI/SNF、ISWI，ATP）、TFIID（TAFs）。"
              "Architectural protein（HMG）協助彎曲 DNA。酵母 GAL 系統：Gal4（activator）、Gal80（抑制）、Gal3 + galactose 解除抑制。"},
    {"id": "dna_binding_motifs", "title": "DNA-binding motifs", "layout": "scene", "scene": MOTIFS, "tags": ["gene_regulation"],
     "notes": "多數以 α-helix 嵌入 major groove（氫鍵多），以 Arg、Gln、Asn 辨識鹼基。HTH：CAP、LacI、Cro、homeodomain。Zinc finger：Cys₂His₂（TFIIIA、Sp1）、"
              "Cys₄（核受體：steroid、thyroid、retinoic acid、vitamin D receptor）。Leucine zipper（bZIP）：每 7 個胺基酸一個 Leu 形成 coiled-coil 二聚，鹼性區抓 DNA（c-Fos/c-Jun = AP-1、CREB）。"
              "bHLH：MyoD、Myc/Max。例外：TBP 以 β-sheet 結合 TATA 的 minor groove。"},
    {"id": "chromatin_epigenetics", "title": "Chromatin states, DNA methylation and imprinting", "layout": "scene", "scene": EPI, "tags": ["gene_regulation"],
     "notes": "Euchromatin（開）：H3/H4 Lys 乙醯化（HAT，中和正電、鬆開），CpG island 低甲基化。Heterochromatin（關）：去乙醯化（HDAC）、H3K9me3 / H3K27me3（HMT；HP1、Polycomb）、"
              "CpG 5mC（DNMT；MeCP2 結合）。H3K4me3 則與活化相關。Imprinting：親源特異的甲基化使單一 allele 表現（IGF2 父源、H19 母源）。IGF2/H19：父源 ICR（insulator）被甲基化 → CTCF 無法結合 → 下游 enhancer 可活化 IGF2；母源 ICR 未甲基化 → CTCF 結合隔絕 → IGF2 關、H19 開（甲基化在此反而促成表現；Alberts 7e p.438）。"
              "15q11–13 缺失：父源 → Prader–Willi；母源（UBE3A）→ Angelman。X 不活化（Xist）也屬表觀遺傳。"},
    {"id": "rnai_mirna_sirna", "title": "RNA interference: miRNA and siRNA", "layout": "scene", "scene": RNAI, "tags": ["gene_regulation"],
     "notes": "miRNA：Pol II 轉錄 pri-miRNA（有 cap、poly(A)）→ 核內 Drosha–DGCR8 切成 ~70 nt pre-miRNA → Exportin-5（Ran-GTP）出核 → Dicer 切成 21–23 nt 雙股 → "
              "guide strand 進 RISC（Argonaute）→ seed（2–8）與 3′UTR 部分配對 → 抑制轉譯、去 poly(A)、分解（P bodies）。"
              "siRNA：外來長 dsRNA → Dicer → RISC（Ago2 slicer）→ 完全配對 → 切斷 mRNA（RNAi，基因 knock-down；Fire & Mello）。piRNA：PIWI，生殖細胞抑制轉位子。"},
    {"id": "nuclear_receptor_signaling", "title": "Nuclear (steroid) receptors", "layout": "scene", "scene": NR, "tags": ["gene_regulation", "signal_transduction"],
     "notes": "Type I（glucocorticoid、mineralocorticoid、sex steroids）：細胞質中與 HSP90 結合；激素穿膜結合後 HSP90 離開 → 同源二聚 → 入核 → 以 Zn-finger DBD 結合 HRE（反向重複）→ 招募 coactivator（HAT）。"
              "Type II（thyroid、retinoic acid、vitamin D、PPAR）：已在核內與 RXR 異源二聚、無激素時結合 corepressor（HDAC），激素結合後換成 coactivator。結構：AF-1、DBD（兩個 Cys₄ 鋅指）、LBD/AF-2。"},
    {"id": "small_regulatory_rnas", "title": "Small regulatory RNAs compared", "layout": "scene", "scene": SRNA_TABLE, "tags": ["gene_regulation"],
     "notes": "miRNA（內生、部分配對、抑制轉譯）、siRNA（外來 dsRNA、完全配對、切斷）、piRNA（PIWI、生殖細胞轉位子）、細菌 sRNA（Hfq）、riboswitch（順式、代謝物）。"},
    {"id": "lac_inducers", "title": "Lactose, allolactose and IPTG", "layout": "free", "cell_size": 13.0,
     "nodes": [{"id": "lac", "ref": "lactose", "pos": [0, 0]}, {"id": "allo", "ref": "allolactose", "pos": [1.2, 0]},
               {"id": "iptg", "ref": "iptg", "pos": [2.2, 0]}],
     "edges": [], "text": False, "tags": ["gene_regulation"],
     "notes": "Lactose（Gal β1→4 Glc）經 β-galactosidase 轉糖苷成 allolactose（Gal β1→6 Glc），後者才是真正的 inducer。IPTG 的 C1 為 S-異丙基（硫糖苷），可誘導但不被水解（gratuitous inducer）。"},
    {"id": "epigenetic_marks", "title": "Epigenetic marks", "layout": "free", "cell_size": 12.0,
     "nodes": [{"id": "mc", "ref": "methylcytosine_dr", "pos": [0, 0], "decorate": {"highlight": [{"smarts": "[CH3]c", "color": "magenta"}]}},
               {"id": "ac", "ref": "acetyllysine", "pos": [1.05, 0]}, {"id": "me3", "ref": "trimethyllysine", "pos": [2.0, 0]}],
     "edges": [], "text": False, "tags": ["gene_regulation"],
     "notes": "DNA：CpG 的 5-methylcytosine（DNMT，SAM）。組蛋白：Lys 乙醯化（HAT，acetyl-CoA；去除正電 → 開）、Lys 三甲基化（HMT，SAM；保留正電，依位置開或關：H3K4me3 開、H3K9me3 / H3K27me3 關）。"},
]


# ================================================================ writer

def main(data):
    data = Path(data)
    for m in MOLECULES:
        d = {"id": m["id"], "name": m["name"], "label": m["label"], "aliases": m["aliases"], "smiles": m["smiles"]}
        for k in ("ring", "haworth"):
            if k in m:
                d[k] = m[k]
        d.update({"tags": with_parent(m["tags"]), "notes_zh": m["notes"], "parts_png": False, "status": "pending"})
        dump(data / "molecules" / f"{m['id']}.yaml", d)
    for r in REACTIONS:
        d = {"id": r["id"], "title": r["title"], "substrates": r["substrates"], "products": r["products"],
             "cofactors_in": r["cofactors_in"], "cofactors_out": r["cofactors_out"]}
        for k in ("extra_in", "extra_out"):
            if r.get(k):
                d[k] = r[k]
        d.update({"enzyme": r["enzyme"], "reversible": False})
        if COLORS.get(r["id"]):
            d["decorate"] = COLORS[r["id"]]
        d.update({"card_zh": r["card"], "tags": with_parent(r["tags"]), "parts_png": False, "status": "pending"})
        dump(data / "reactions" / f"{r['id']}.yaml", d)
    for p in PATHWAYS:
        d = {k: v for k, v in p.items() if k not in ("tags", "notes")}
        d["tags"] = with_parent(p["tags"])
        d["notes_zh"] = p["notes"]
        d = tidy(d)
        d["parts_png"] = False
        d["status"] = "pending"
        dump(data / "pathways" / f"{p['id']}.yaml", d)
    print(f"{len(MOLECULES)} molecules, {len(REACTIONS)} reactions, {len(PATHWAYS)} pathways")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data")
