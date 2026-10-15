"""Round 23: molecular biology VI — techniques (PCR / qPCR, Sanger, blotting, cloning and selection, cDNA, CRISPR,
DNA–protein assays, yeast two-hybrid, RFLP; vector / technique tables; X-gal chemistry).
Figures carry names and short labels only; explanations live in notes_zh.

    python tools/gen_mb_tech.py data
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from figure_text import tidy  # noqa: E402
from gen_mb_repair import GAP, down, dump, duplex, seg, table, txt, with_parent  # noqa: E402
from gen_mb_reg import inhibit, off, on, region  # noqa: E402
from gen_mb_tx import ball, box, prot, rna  # noqa: E402
from reaction_colors import COLORS, hl  # noqa: E402

M, C = "magenta", "cyan"
BASE = {"A": "#5BB55B", "C": "#4C8BD6", "G": "#E0A100", "T": "#E2574C"}


def head(x, y, text, size=0.85, color=None):
    return txt(x, y, text, size, color, "left")


def strand(x0, x1, y, kind="parent", ends=None, arrow=False, w=0.3):
    d = {"g": "strand", "kind": kind, "points": [[x0, y], [x1, y]], "width": w}
    if ends:
        d["ends"] = list(ends)
    if arrow:
        d["arrow"] = True
    return d


def band(x, y, w=1.6, hue="#3F4F6B", h=0.28):
    return {"g": "region", "x0": x - w / 2, "y0": y - h / 2, "x1": x + w / 2, "y1": y + h / 2, "hue": hue, "opacity": 1.0, "outline": False}


def gel(x0, y0, x1, y1):
    return {"g": "region", "x0": x0, "y0": y0, "x1": x1, "y1": y1, "hue": "gray", "opacity": 0.18}


def fat(x0, y0, x1, y1):
    return {"g": "arrow", "from": [x0, y0], "to": [x1, y1], "style": "fat"}


def arr(x0, y0, x1, y1, color="black"):
    return {"g": "arrow", "from": [x0, y0], "to": [x1, y1], "color": color, "width": 1.0}


# ================================================================ PCR (notes p.83; lecture 4-224)
PX = [0.0, 11.5, 23.0, 34.5]
PCR = {"items": [
    head(PX[0], -0.8, "template"), *duplex(1.0, PX[0], PX[0] + 8.0, top="parent", bot="parent"),
    head(PX[1], -0.8, "1. denature 95 °C"),
    strand(PX[1], PX[1] + 8.0, 0.6, ends=("5'", "3'")), strand(PX[1], PX[1] + 8.0, 3.2, ends=("3'", "5'")),
    head(PX[2], -0.8, "2. anneal 55 °C"),
    strand(PX[2], PX[2] + 8.0, 0.6, ends=("5'", "3'")), strand(PX[2] + 6.2, PX[2] + 8.0, 1.4, "new"),
    strand(PX[2], PX[2] + 8.0, 3.2, ends=("3'", "5'")), strand(PX[2], PX[2] + 1.8, 2.4, "new"),
    txt(PX[2] + 4.0, 4.4, "primers", 0.68, "#D9622B"),
    head(PX[3], -0.8, "3. extend 72 °C"),
    strand(PX[3], PX[3] + 8.0, 0.6, ends=("5'", "3'")),
    {"g": "strand", "kind": "new", "points": [[PX[3] + 8.0, 1.4], [PX[3] + 0.8, 1.4]], "width": 0.3, "arrow": True},
    strand(PX[3], PX[3] + 8.0, 3.2, ends=("3'", "5'")), strand(PX[3], PX[3] + 7.2, 2.4, "new", arrow=True),
    txt(PX[3] + 4.0, 4.4, "Taq polymerase", 0.68, "#D9622B"),
    *[fat(PX[k] + 9.1, 1.9, PX[k + 1] - 0.9, 1.9) for k in range(3)],
    {"g": "arrow", "from": [PX[3] + 4.0, 5.2], "to": [PX[1] + 4.0, 5.2], "curve": -0.12, "color": "gray", "width": 1.0},
    txt(PX[2] + 2.0, 9.0, "25–35 cycles → 2^n copies", 0.8, "gray"),
]}

# ================================================================ qPCR (notes p.84; lecture 4-237)
QPCR = {"size": [15, 9.5], "params": {"A": 22, "B": 30, "n": 9, "th": 0.15, "k": 0.8248},
        "x": {"min": 0, "max": 40, "label": "cycle", "ticks": [0, 10, 20, 30, 40]},
        "y": {"min": 0, "max": 1.15, "label": "fluorescence"},
        "hlines": [{"y": "th", "label": "threshold"}],
        "curves": [{"id": "a", "y": "hill(x, A/k, n)", "color": "cyan"}, {"id": "b", "y": "hill(x, B/k, n)", "color": "red"}],
        "points": [{"on": "a", "y": "th", "drop": "x", "dot": True, "label": "C_t = 22"},
                   {"on": "b", "y": "th", "drop": "x", "dot": True, "label": "C_t = 30"}],
        "legend": {"at": [2, 1.05], "items": [{"curve": "a", "text": "sample A"}, {"curve": "b", "text": "sample B"}]}}

# ================================================================ Sanger (notes p.81–82; lecture 4-236)
SEQ = "ATGCCAGT"                     # new strand 5'→3'
GY0, GSTEP = 13.6, 1.05
SANGER = {"items": [
    head(0.0, -1.0, "ddNTP chain termination"),
    strand(0.0, 11.0, 0.4, "parent", ("3'", "5'")),
    *[txt(2.4 + 0.95 * k, 1.05, "TACGGTCA"[k], 0.6, "gray") for k in range(8)],
    *[x for k in range(8) for x in (strand(0.0, 2.0 + 0.95 * k, 2.2 + 0.95 * k, "new"),
                                     ball(2.35 + 0.95 * k, 2.2 + 0.95 * k, "", BASE[SEQ[k]] and {"A": "green", "C": "blue", "G": "yellow", "T": "red"}[SEQ[k]], 0.3))],
    txt(1.0, 10.6, "primer", 0.65, "#D9622B"),
    head(14.0, -1.0, "gel"),
    gel(14.0, 0.2, 22.6, GY0 + 0.8),
    *[txt(15.2 + 2.0 * j, 0.9, "dd" + b, 0.68) for j, b in enumerate("ACGT")],
    *[band(15.2 + 2.0 * "ACGT".index(SEQ[k]), GY0 - GSTEP * k, 1.5, BASE[SEQ[k]]) for k in range(8)],
    {"g": "arrow", "from": [23.6, GY0], "to": [23.6, GY0 - GSTEP * 7.4], "color": "black", "width": 1.0},
    *[txt(24.4, GY0 - GSTEP * k, SEQ[k], 0.75, BASE[SEQ[k]], "left") for k in range(8)],
    txt(25.4, GY0 + 0.1, "5'", 0.62, "gray", "left"), txt(25.4, GY0 - GSTEP * 7, "3'", 0.62, "gray", "left"),
    txt(18.3, GY0 + 1.6, "short fragments run farthest", 0.65, "gray"),
]}

# ================================================================ blotting (lecture 4-237; notes p.81)
BY = 3.0


def antibody(x, y, hue="#7B4FC0"):
    return [rna([(x, y + 0.9), (x, y)], kind=hue, w=0.22), rna([(x, y), (x - 0.5, y - 0.6)], kind=hue, w=0.22),
            rna([(x, y), (x + 0.5, y - 0.6)], kind=hue, w=0.22)]


BLOT = {"items": [
    head(0.0, -0.6, "1. gel"),
    gel(0.0, 0.2, 5.4, 6.0), *[band(1.3 + 1.4 * j, 1.4 + 1.1 * k, 1.0, "#9A9A9A", 0.22) for j in range(3) for k in range(4)],
    fat(5.9, BY, 7.7, BY),
    head(8.2, -0.6, "2. transfer"),
    {"g": "region", "x0": 8.2, "y0": 0.2, "x1": 13.6, "y1": 6.0, "hue": "yellow", "opacity": 0.25},
    *[band(9.5 + 1.4 * j, 1.4 + 1.1 * k, 1.0, "#BDBDBD", 0.22) for j in range(3) for k in range(4)],
    fat(14.1, BY, 15.9, BY),
    head(16.4, -0.6, "3. probe"),
    strand(16.6, 19.0, 2.0, "new"), {"g": "mark", "kind": "lesion", "at": [19.4, 2.0], "r": 0.3, "color": "orange"},
    *antibody(17.6, 4.0), {"g": "mark", "kind": "lesion", "at": [17.6, 5.2], "r": 0.3, "color": "orange"},
    fat(20.3, BY, 22.1, BY),
    head(22.6, -0.6, "4. detect"),
    {"g": "region", "x0": 22.6, "y0": 0.2, "x1": 28.0, "y1": 6.0, "hue": "gray", "opacity": 0.12},
    band(23.9, 2.5, 1.0, "#222222", 0.3), band(25.3, 2.5, 1.0, "#222222", 0.3), band(26.7, 4.7, 1.0, "#222222", 0.3),
    box(0.0, 6.0, 8.4, "blue", "Southern", 1.0, 0.7), txt(6.6, 8.4, "DNA ← DNA probe", 0.75, None, "left"),
    box(0.0, 6.0, 9.9, "pink", "Northern", 1.0, 0.7), txt(6.6, 9.9, "RNA ← DNA probe", 0.75, None, "left"),
    box(0.0, 6.0, 11.4, "purple", "Western", 1.0, 0.7), txt(6.6, 11.4, "protein ← antibody (SDS-PAGE)", 0.75, None, "left"),
    box(16.0, 22.0, 8.4, "teal", "Southwestern", 1.0, 0.7), txt(22.6, 8.4, "protein ← DNA probe", 0.75, None, "left"),
]}

# ================================================================ restriction enzymes and cloning (notes p.85–86; lecture 4-243)
RY = 1.2
STAIR = [(1.95, RY - 0.6), (1.95, RY + 0.5), (5.0, RY + 0.5), (5.0, RY + 1.6)]
CLONE = {"items": [
    head(0.0, -1.0, "EcoRI: sticky ends"),
    *[txt(0.6 + 1.0 * k, RY, b, 0.85, "#3F4F6B") for k, b in enumerate("GAATTC")],
    *[txt(0.6 + 1.0 * k, RY + 1.0, b, 0.85, "#3F4F6B") for k, b in enumerate("CTTAAG")],
    *[{"g": "strand", "kind": "red", "points": [list(a), list(b)], "width": 0.08} for a, b in (((1.1, RY - 0.6), (1.1, RY + 0.5)), ((1.1, RY + 0.5), (5.1, RY + 0.5)), ((5.1, RY + 0.5), (5.1, RY + 1.6)))],
    txt(-0.5, RY, "5'", 0.6, "gray"), txt(6.6, RY, "3'", 0.6, "gray"), txt(-0.5, RY + 1.0, "3'", 0.6, "gray"), txt(6.6, RY + 1.0, "5'", 0.6, "gray"),
    fat(7.6, RY + 0.5, 9.4, RY + 0.5),
    strand(10.0, 11.4, RY), strand(10.0, 14.0, RY + GAP), strand(15.4, 19.4, RY), strand(18.0, 19.4, RY + GAP),
    txt(12.7, RY + 1.9, "TTAA", 0.6, "#D62828"), txt(16.7, RY - 0.8, "AATT", 0.6, "#D62828"),
    head(0.0, 5.0, "Cloning"),
    {"g": "plasmid", "at": [3.0, 9.0], "kind": "relaxed", "r": 1.8},
    {"g": "mark", "kind": "cut", "at": [3.0, 7.2], "side": "up"}, txt(3.0, 11.6, "vector (MCS cut)", 0.68, "gray"),
    *duplex(8.4, 7.0, 12.0, top="new", bot="new", ends=False), txt(9.5, 10.6, "insert (EcoRI ends)", 0.68, "gray"),
    fat(13.0, 9.0, 15.0, 9.0), txt(14.0, 7.7, "DNA ligase", 0.7),
    {"g": "plasmid", "at": [18.4, 9.0], "kind": "relaxed", "r": 1.8},
    rna([(17.2, 7.65), (18.4, 7.2), (19.6, 7.65)], kind="new", w=0.42), txt(18.4, 11.6, "recombinant plasmid", 0.68, "gray"),
    fat(21.0, 9.0, 23.0, 9.0), txt(22.0, 7.7, "transform", 0.7),
    {"g": "region", "x0": 23.8, "y0": 7.4, "x1": 27.4, "y1": 10.6, "hue": "blue", "opacity": 0.15},
    {"g": "plasmid", "at": [25.6, 9.0], "kind": "relaxed", "r": 0.6}, txt(25.6, 11.4, "E. coli", 0.68, "gray"),
]}

# ================================================================ blue-white and replica plating (notes p.86; lecture 4-238)
WY = 4.5


def plate(x, y, colonies, r=2.4):
    out = [{"g": "region", "x0": x - r, "y0": y - r, "x1": x + r, "y1": y + r, "hue": "yellow", "opacity": 0.25}]
    out += [ball(x + dx, y + dy, "", hue_, 0.32) for dx, dy, hue_ in colonies]
    return out


SPOTS = [(-1.2, -1.0), (0.4, -1.4), (1.3, -0.2), (-0.6, 0.6), (0.8, 1.2), (-1.5, 0.9), (0.0, -0.1)]
BW = {"items": [
    head(0.0, -0.6, "Blue-white screening (X-gal)"),
    {"g": "plasmid", "at": [2.6, WY], "kind": "relaxed", "r": 1.8},
    rna([(1.4, WY - 1.35), (2.6, WY - 1.8), (3.8, WY - 1.35)], kind="#4C8BD6", w=0.42), txt(2.6, WY - 2.8, "lacZα (MCS)", 0.65, "#2E6E9E"),
    rna([(1.2, WY + 1.2), (2.6, WY + 1.8), (4.0, WY + 1.2)], kind="#E0A100", w=0.42), txt(2.6, WY + 2.8, "amp^R", 0.65, "#B08F10"),
    *plate(10.0, WY, [(dx, dy, "#F2F2F2" if k in (1, 4) else "blue") for k, (dx, dy) in enumerate(SPOTS)]),
    txt(13.2, WY - 1.2, "blue: no insert", 0.72, "#2E6E9E", "left"), txt(13.2, WY + 0.2, "white: insert", 0.72, None, "left"),
    txt(13.2, WY + 1.6, "(lacZα disrupted)", 0.65, "gray", "left"),
    head(0.0, 9.6, "Replica plating (pBR322)"),
    *plate(3.0, 13.6, [(dx, dy, "brown") for dx, dy in SPOTS]), txt(3.0, 16.6, "ampicillin", 0.7),
    fat(6.0, 13.6, 8.0, 13.6),
    *plate(11.0, 13.6, [(dx, dy, "brown") for k, (dx, dy) in enumerate(SPOTS) if k not in (1, 4)]), txt(11.0, 16.6, "tetracycline", 0.7),
    *[{"g": "region", "x0": 11.0 + dx - 0.32, "y0": 13.6 + dy - 0.32, "x1": 11.0 + dx + 0.32, "y1": 13.6 + dy + 0.32,
       "hue": "red", "opacity": 0.0, "dashed": True} for k, (dx, dy) in enumerate(SPOTS) if k in (1, 4)],
    txt(14.4, 13.0, "missing on Tet = insert in tet^R", 0.72, "#D62828", "left"),
]}

# ================================================================ cDNA synthesis (lecture 4-240)
CY = [1.0, 4.2, 7.4, 10.6]
CDNA = {"items": [
    rna([(0.0, CY[0]), (12.0, CY[0])], ("5'", "")), txt(12.2, CY[0], "AAAA", 0.7, "#C2338A", "left"), ball(-0.6, CY[0], "", "green", 0.3),
    strand(12.2, 14.0, CY[0] + 0.8, "new", ("3'", "5'")), txt(15.3, CY[0] + 0.8, "oligo(dT)", 0.65, "#D9622B", "left"),
    down(6.0, CY[0] + 0.8, CY[1] - 0.8), txt(6.6, (CY[0] + CY[1]) / 2, "reverse transcriptase", 0.72, None, "left"),
    rna([(0.0, CY[1]), (12.0, CY[1])], ("5'", "")), txt(12.2, CY[1], "AAAA", 0.7, "#C2338A", "left"),
    {"g": "strand", "kind": "new", "points": [[14.0, CY[1] + 0.8], [0.6, CY[1] + 0.8]], "width": 0.3, "arrow": True, "ends": ["5'", ""]},
    down(6.0, CY[1] + 1.6, CY[2] - 0.8), txt(6.6, (CY[1] + CY[2]) / 2 + 0.3, "RNase H nicks RNA", 0.72, None, "left"),
    *[rna([(a, CY[2]), (b, CY[2])]) for a, b in ((1.0, 2.6), (5.0, 6.0), (9.0, 10.2))],
    strand(0.6, 14.0, CY[2] + 0.8, "new", ("", "5'")),
    down(6.0, CY[2] + 1.6, CY[3] - 0.8), txt(6.6, (CY[2] + CY[3]) / 2 + 0.3, "DNA Pol I + ligase", 0.72, None, "left"),
    strand(0.6, 14.0, CY[3], "parent", ("5'", "3'")), strand(0.6, 14.0, CY[3] + 0.8, "new", ("3'", "5'")),
    txt(7.3, CY[3] + 1.9, "ds cDNA (no introns)", 0.72),
]}

# ================================================================ CRISPR–Cas9 (notes p.82; lecture 4-236)
KY = 3.2
CRISPR = {"items": [
    prot(8.0, KY - 0.3, 9.0, 4.4, "#C9C3DD", "", "bilobe", 0.7),
    *duplex(KY, 0.0, 16.0, top="parent", bot="parent", ends=True),
    box(11.2, 12.8, KY + 0.5, "yellow", "NGG", 1.5, 0.6),
    rna([(5.0, KY + 0.5), (10.6, KY + 0.5)], kind="#C2338A", w=0.26),
    rna([(10.6, KY + 0.5), (11.6, KY - 1.0), (10.8, KY - 2.0), (9.6, KY - 1.6), (8.4, KY - 2.1)], kind="#C2338A", w=0.26),
    txt(-0.2, KY - 4.3, "Cas9 + sgRNA", 0.8, "#5F4A8C", "left"), txt(7.8, KY + 2.9, "20-nt guide", 0.68, "#C2338A"),
    txt(12.0, KY + 2.9, "PAM", 0.68, "#B08F10"),
    {"g": "mark", "kind": "cut", "at": [10.2, KY], "side": "up"}, {"g": "mark", "kind": "cut", "at": [10.2, KY + GAP], "side": "down"},
    txt(0.0, 8.2, "NHEJ", 0.85, None, "left"),
    *duplex(9.4, 0.0, 7.0, top="parent", bot="parent", ends=False), {"g": "mark", "kind": "lesion", "at": [3.5, 9.9], "r": 0.35},
    txt(3.5, 11.4, "indel → knock-out", 0.72, "#D62828"),
    txt(9.0, 8.2, "HDR + donor", 0.85, None, "left"),
    *duplex(9.4, 9.0, 16.0, top="parent", bot="parent", ends=False), *region(11.6, 13.4, 9.4, "green", "new"),
    txt(12.5, 11.4, "precise edit / knock-in", 0.72, "#2E9E4F"),
]}

# ================================================================ EMSA, footprinting, ChIP (lecture 4-240; notes p.85)
EM = {"items": [
    head(0.0, -0.8, "EMSA"),
    gel(0.0, 0.0, 7.6, 8.0), *[txt(1.2 + 1.7 * j, 0.7, t, 0.62) for j, t in enumerate(["probe", "+P", "+P+Ab", "+cold"])],
    *[band(1.2 + 1.7 * j, 7.0, 1.2, "#3F4F6B") for j in (0, 3)], band(4.6, 7.0, 1.2, "#9A9A9A"),
    band(2.9, 4.4, 1.2, "#3F4F6B"), band(4.6, 2.4, 1.2, "#3F4F6B"), band(6.3, 4.4, 1.2, "#BDBDBD"),
    txt(8.0, 4.4, "shift", 0.65, None, "left"), txt(8.0, 2.4, "supershift", 0.65, None, "left"), txt(8.0, 7.0, "free", 0.65, None, "left"),
    head(12.0, -0.8, "DNase I footprint"),
    gel(12.0, 0.0, 17.0, 8.0), txt(13.2, 0.7, "−P", 0.62), txt(15.8, 0.7, "+P", 0.62),
    *[band(13.2, 1.6 + 0.45 * k, 1.2, "#3F4F6B", 0.2) for k in range(14)],
    *[band(15.8, 1.6 + 0.45 * k, 1.2, "#3F4F6B", 0.2) for k in range(14) if k not in (5, 6, 7, 8)],
    {"g": "region", "x0": 17.2, "y0": 3.75, "x1": 17.5, "y1": 5.65, "hue": "red", "opacity": 0.9, "outline": False},
    txt(17.8, 4.7, "binding site", 0.65, "#D62828", "left"),
    head(0.0, 10.2, "ChIP"),
    *duplex(12.4, 0.0, 4.4, top="parent", bot="parent", ends=False), prot(2.2, 11.7, 1.4, 1.1, "teal", "TF"),
    txt(2.2, 14.4, "cross-link", 0.65, "gray"), fat(4.8, 12.6, 6.0, 12.6),
    *duplex(12.4, 6.4, 8.6, top="parent", bot="parent", ends=False), prot(7.5, 11.7, 1.4, 1.1, "teal", "TF"),
    txt(7.5, 14.4, "shear", 0.65, "gray"), fat(9.0, 12.6, 10.2, 12.6),
    *duplex(12.4, 10.6, 12.8, top="parent", bot="parent", ends=False), prot(11.7, 11.7, 1.4, 1.1, "teal", "TF"),
    {"g": "strand", "kind": "#7B4FC0", "points": [[11.7, 11.1], [11.7, 9.9]], "width": 0.22},
    txt(11.7, 14.4, "antibody IP", 0.65, "gray"), fat(13.4, 12.6, 15.2, 12.6),
    *duplex(12.4, 15.8, 18.0, top="parent", bot="parent", ends=False), txt(17.4, 14.4, "PCR / sequencing", 0.65, "gray"),
]}

# ================================================================ yeast two-hybrid (notes p.84–85; lecture 4-238)
YY = [4.0, 11.0]


def y2h(y, bind):
    out = [*duplex(y, 0.0, 16.0, top="parent", bot="parent", ends=False), *region(1.0, 3.6, y, "gray", "UAS"),
           *region(7.0, 15.0, y, "blue", "reporter (lacZ, HIS3)", 0.45),
           prot(2.3, y - 0.8, 2.2, 1.4, "orange", "BD", "box", 0.6), prot(2.3, y - 2.3, 1.8, 1.4, "pink", "X", "round", 0.6)]
    if bind:
        out += [prot(3.8, y - 3.4, 1.8, 1.4, "green", "Y", "round", 0.6), prot(5.4, y - 3.4, 1.8, 1.3, "teal", "AD", "box", 0.6),
                prot(6.4, y - 1.1, 2.4, 1.6, "blue", "Pol II", "oval", 0.5), on(7.6, 14.6, y - 0.6)]
    else:
        out += [prot(8.6, y - 3.4, 1.8, 1.4, "purple", "Z", "round", 0.6), prot(10.2, y - 3.4, 1.8, 1.3, "teal", "AD", "box", 0.6), off(10.0, y - 0.9)]
    return out


Y2H = {"items": [
    *y2h(YY[0], True), txt(17.0, YY[0] + 0.5, "X binds Y → ON", 0.8, "#2E9E4F", "left"),
    *y2h(YY[1], False), txt(17.0, YY[1] + 0.5, "no binding → OFF", 0.8, "#D62828", "left"),
    txt(2.3, YY[1] + 2.3, "bait: BD–X", 0.68, "gray"), txt(9.4, YY[1] + 2.3, "prey: AD–Y / AD–Z", 0.68, "gray"),
]}

# ================================================================ RFLP: sickle cell (notes p.81)
RF = [2.0, 6.5]


def rflp_row(y, sites, label):
    out = [*duplex(y, 0.0, 16.0, top="parent", bot="parent", ends=False), *region(5.4, 10.4, y, "teal", "β-globin", 0.4),
           {"g": "region", "x0": 2.0, "y0": y + GAP + 0.5, "x1": 6.0, "y1": y + GAP + 0.8, "hue": "orange", "opacity": 0.9, "outline": False},
           txt(-0.6, y + 0.5, label, 0.75, None, "right")]
    out += [{"g": "mark", "kind": "cut", "at": [x, y], "side": "up"} for x in sites]
    return out


RFLP = {"items": [
    head(0.0, -1.0, "MstII sites (CCTNAGG)"),
    *rflp_row(RF[0], (1.0, 6.6, 8.6), "β^A"), txt(3.8, RF[0] - 1.2, "1.15 kb", 0.65, "gray"), txt(7.6, RF[0] - 1.2, "0.2 kb", 0.65, "gray"),
    *rflp_row(RF[1], (1.0, 8.6), "β^S"), txt(4.8, RF[1] - 1.2, "1.35 kb (site lost: GAG → GTG)", 0.65, "gray"),
    txt(4.0, RF[1] + 2.6, "probe", 0.65, "#B5651F"),
    head(19.0, -1.0, "Southern blot"),
    gel(19.0, 0.0, 27.0, 9.0), *[txt(20.4 + 2.6 * j, 0.7, t, 0.65) for j, t in enumerate(["AA", "AS", "SS"])],
    txt(18.6, 3.0, "1.35", 0.62, "gray", "right"), txt(18.6, 4.4, "1.15", 0.62, "gray", "right"),
    band(20.4, 4.4, 1.6), band(23.0, 3.0, 1.6), band(23.0, 4.4, 1.6), band(25.6, 3.0, 1.6),
]}

# ================================================================ tables
VECTOR_TABLE = table([
    ("Plasmid", "E. coli · up to ~15 kb · subcloning, cDNA library; ori + selectable marker + MCS"),
    ("λ phage", "E. coli · up to ~25 kb · genomic or cDNA library"),
    ("Cosmid", "E. coli · 30–45 kb · λ cos sites + plasmid ori · genomic library"),
    ("BAC", "E. coli (F plasmid) · 100–300 kb · genome projects"),
    ("YAC", "yeast · 250–2000 kb · ARS + centromere + telomeres"),
    ("Shuttle / expression", "two origins (E. coli + yeast) / promoter + RBS for protein production; baculovirus–insect cells"),
], w_name=7.8, x_desc=8.4, hue="teal")

TECH_TABLE = table([
    ("DNA", "Southern, PCR, Sanger / NGS (pyrosequencing), RFLP, FISH (location), microarray"),
    ("RNA", "Northern, RT-PCR, qRT-PCR, nuclease S1 / RNase protection, microarray, RNA-seq"),
    ("Protein", "Western, ELISA, immunoprecipitation; GST / His-tag purification"),
    ("DNA–protein", "filter binding, EMSA, DNase I footprinting, ChIP, yeast one-hybrid"),
    ("Protein–protein", "co-IP, yeast two-hybrid, FRET (GFP variants), pull-down"),
    ("RNA–protein", "yeast three-hybrid, CLIP"),
    ("Gene function", "knock-in (HR), knock-out (Cre–loxP, CRISPR), knock-down (siRNA / miRNA)"),
    ("DNA delivery", "transformation (CaCl_2), electroporation, microinjection, lipofection, viral vectors, Ti plasmid (plants)"),
], w_name=6.8, x_desc=7.4, hue="purple")

# ================================================================ molecules and reactions
MOLECULES = [
    {"id": "xgal", "name": "X-gal (5-bromo-4-chloro-3-indolyl β-D-galactoside)", "label": "X-gal", "aliases": ["X-gal"],
     "smiles": "OC[C@H]1O[C@@H](Oc2c[nH]c3ccc(Br)c(Cl)c23)[C@H](O)[C@@H](O)[C@H]1O", "ring": {"coordgen": True}, "tags": ["techniques"],
     "notes": "β-galactosidase 的人工受質：水解出 galactose 與 5-bromo-4-chloro-indoxyl；後者在空氣中兩分子氧化二聚成藍色靛藍（不溶）→ 藍色菌落。"},
    {"id": "bromochloroindoxyl", "name": "5-Bromo-4-chloro-3-hydroxyindole", "label": "5-Bromo-4-chloro-indoxyl", "aliases": ["indoxyl"],
     "smiles": "Oc1c[nH]c2ccc(Br)c(Cl)c12", "ring": {"coordgen": True}, "tags": ["techniques"],
     "notes": "X-gal 水解產物（無色），兩分子經 O₂ 氧化二聚為 5,5′-dibromo-4,4′-dichloro-indigo（藍）。"},
    {"id": "dibromodichloroindigo", "name": "5,5′-Dibromo-4,4′-dichloro-indigo", "label": "Blue dye (dibromo-dichloro-indigo)", "aliases": [],
     "smiles": "O=C1/C(=C2\\Nc3ccc(Br)c(Cl)c3C2=O)Nc2ccc(Br)c(Cl)c21", "ring": {"coordgen": True}, "tags": ["techniques"],
     "notes": "藍白篩選的藍色來源：不溶於水，沉積在菌落中。"},
]

REACTIONS = [
    {"id": "xgal_hydrolysis", "title": "X-gal hydrolysis (β-galactosidase)", "substrates": ["xgal"], "products": ["b_galactopyranose", "bromochloroindoxyl"],
     "cofactors_in": ["H2O"], "cofactors_out": [], "enzyme": "β-galactosidase (lacZ)", "tags": ["techniques"],
     "card": {"type": "水解 β-糖苷鍵：indoxyl 的 O（青）留在 indoxyl 上；H₂O（洋紅）的 OH 接到 galactose C1",
              "exam_points": ["藍白篩選：lacZα 完整 → 有 β-gal 活性 → 藍；插入片段破壞 lacZα → 白（含重組質體）",
                              "需 IPTG 解除 LacI 抑制；α-complementation：宿主提供 lacZΔM15"]}},
]

COLORS.update({
    "xgal_hydrolysis": {"H2O": {"colors": {0: M}}, "b_galactopyranose": {"colors": {"u1": M}},
                        "xgal": hl("[O:1]c1c[nH]c2ccc(Br)c(Cl)c12", col=C), "bromochloroindoxyl": hl("[OH:1]c", col=C)},
})

PATHWAYS = [
    {"id": "pcr_cycle", "title": "Polymerase chain reaction (PCR)", "layout": "scene", "scene": PCR, "tags": ["techniques"],
     "notes": "一對 primer（15–25 nt，分別結合兩股的 3′ 端）、Taq DNA polymerase（Thermus aquaticus，耐熱、無 3′→5′ 校對）、dNTP、Mg²⁺。"
              "每循環：95 °C 變性 → 50–65 °C 黏合（約 Tm − 5；Tm ≈ 4(G+C) + 2(A+T)）→ 72 °C 延長。n 個循環產物約 2ⁿ（嚴格長度的 amplicon = 2ⁿ − 2n）。"
              "應用：診斷、STR／DNA 指紋、RAPD、RT-PCR（先反轉錄）。高保真 PCR 用 Pfu（有校對）。"},
    {"id": "qpcr_ct", "title": "Real-time PCR and Ct", "layout": "plot", "plot": QPCR, "tags": ["techniques"],
     "notes": "Real-time（quantitative）PCR：以 SYBR Green 或 TaqMan 探針即時偵測產物。Ct（threshold cycle）＝螢光越過閾值的循環數；起始模板越多 Ct 越小。"
              "每差 1 個 Ct ≈ 2 倍：A 的 Ct 22、B 的 Ct 30 → ΔCt = 8 → A 的起始量約為 B 的 2⁸ = 256 倍。mRNA 定量用 qRT-PCR（先以 reverse transcriptase 做 cDNA）。"},
    {"id": "sanger_sequencing", "title": "Sanger (dideoxy) sequencing", "layout": "scene", "scene": SANGER, "tags": ["techniques"],
     "notes": "Primer + DNA polymerase + 4 種 dNTP + 少量 ddNTP（缺 3′-OH）→ 隨機在各位置終止 → 一組長度差 1 nt 的片段。"
              "四管（ddA、ddC、ddG、ddT）各跑一道膠，由下（短）往上讀 = 新股 5′→3′；與模板互補。現代用四色螢光 ddNTP + 毛細管電泳。"
              "Maxam–Gilbert 為化學切割法；NGS（如 pyrosequencing：PPi → ATP sulfurylase → luciferase 發光）。"},
    {"id": "blotting_methods", "title": "Southern, Northern and Western blotting", "layout": "scene", "scene": BLOT, "tags": ["techniques"],
     "notes": "電泳分離 → 轉漬到膜（nitrocellulose / nylon）→ 以標記探針或抗體雜交 → 偵測（放射自顯影、化學發光）。"
              "Southern（DNA；E. Southern）：限制酶切 → 電泳 → 鹼變性 → 轉漬 → DNA 探針；考題順序：切割 → 電泳 → 轉漬 → 雜交 → 自顯影。"
              "Northern：RNA（不需限制酶）。Western：SDS-PAGE 分離蛋白質 → 一級、二級抗體。Southwestern：蛋白質以 DNA 探針偵測（DNA 結合蛋白）。"},
    {"id": "restriction_cloning", "title": "Restriction enzymes and gene cloning", "layout": "scene", "scene": CLONE, "tags": ["techniques"],
     "notes": "Type II 限制酶：辨識 4–8 bp 回文、只需 Mg²⁺、不需 ATP、切割與甲基化分開 → 實驗用。命名：EcoRI = E. coli、R 株、第一個。"
              "細菌的 restriction–modification：自己的位點被 methylase 甲基化保護，外來（噬菌體）DNA 被切。EcoRI（G↓AATTC）產生 5′ 突出的黏端；SmaI 產生鈍端。"
              "4ⁿ：EcoRI 平均每 4⁶ ≈ 4 kb 切一次。Isoschizomer（同位點同切點）、neoschizomer（同位點不同切點）。同一酶切載體與插入片段 → DNA ligase 接合 → 轉形。"},
    {"id": "blue_white_selection", "title": "Clone selection: blue-white and replica plating", "layout": "scene", "scene": BW, "tags": ["techniques"],
     "notes": "質體帶 amp^R（只讓轉形菌生長）與 lacZα（MCS 位於其中）。宿主 lacZΔM15 + 質體 lacZα → α-complementation 形成 β-gal → X-gal 變藍；"
              "插入片段破壞 lacZα → 白色菌落 = 含重組質體。需 IPTG 誘導。pBR322：amp^R 與 tet^R，插入 tet^R → 可長在 ampicillin、不能長在 tetracycline（replica plating 挑出）。"},
    {"id": "cdna_synthesis", "title": "cDNA synthesis", "layout": "scene", "scene": CDNA, "tags": ["techniques"],
     "notes": "mRNA 的 poly(A) 以 oligo(dT) 為引子 → reverse transcriptase 合成第一股（RNA–DNA hybrid）→ RNase H 切 RNA 留下小段當引子 → DNA Pol I 合成第二股、ligase 接合。"
              "（另法：鹼／S1 nuclease 去除髮夾。）cDNA 沒有 intron，可在細菌表現真核蛋白；cDNA library 代表某組織當下表現的基因，genomic library 代表全基因體。"},
    {"id": "crispr_cas9", "title": "CRISPR–Cas9 genome editing", "layout": "scene", "scene": CRISPR, "tags": ["techniques"],
     "notes": "細菌的適應性免疫（spacer 記住噬菌體序列）。Cas9 + sgRNA（20 nt 引導 + scaffold）找到互補 DNA，緊鄰 PAM（NGG）時，在 PAM 上游 3 bp 切斷雙股（HNH、RuvC）。"
              "NHEJ → 插入／缺失 → frameshift 剔除（knock-out）；有 donor 模板時 HDR → 精確修改或嵌入（knock-in）。dCas9（無切割）+ 抑制子 = CRISPRi（knock-down 轉錄）、+ 活化子 = CRISPRa。"
              "比較：Cre–loxP（條件式剔除）、siRNA / miRNA（mRNA 層次 knock-down）。"},
    {"id": "dna_protein_assays", "title": "EMSA, DNase I footprinting and ChIP", "layout": "scene", "scene": EM, "tags": ["techniques"],
     "notes": "EMSA（gel shift）：標記 DNA 與蛋白質結合後跑得較慢 → 上移；加抗體 → supershift；加過量未標記 DNA（cold competitor）→ 競爭掉。"
              "DNase I footprinting：結合處受保護，ladder 出現空缺（footprint）→ 定出結合序列；methylation interference / protection 亦可。"
              "ChIP：甲醛交聯 → 切碎染色質 → 抗體沉澱蛋白質 → 解交聯 → PCR / 定序（ChIP-seq）找體內結合位點。其他：filter binding、yeast one-hybrid、reporter（luciferase）assay。"},
    {"id": "yeast_two_hybrid", "title": "Yeast two-hybrid", "layout": "scene", "scene": Y2H, "tags": ["techniques"],
     "notes": "Gal4 分成 DNA-binding domain（BD，結合 UAS）與 activation domain（AD）。BD–X（bait）與 AD–Y（prey，可為 cDNA library）共表現；"
              "X 與 Y 結合才把 AD 帶到啟動子 → reporter（lacZ、HIS3）表現。用於蛋白–蛋白交互作用；one-hybrid 測 DNA–蛋白、three-hybrid 測 RNA–蛋白。偽陽性：自活化的 bait。"},
    {"id": "rflp_sickle_cell", "title": "RFLP: sickle cell (MstII)", "layout": "scene", "scene": RFLP, "tags": ["techniques"],
     "notes": "HbS：β-globin 第 6 個密碼子 GAG（Glu）→ GTG（Val），破壞 MstII 位點（CCTNAGG）。以 β-globin 探針做 Southern：正常 1.15 kb（+ 0.2 kb），鐮刀 1.35 kb；"
              "AA 一條 1.15、SS 一條 1.35、AS 兩條。RFLP 也用於連鎖分析與親子鑑定；現多以 PCR 加 STR 取代。"},
    {"id": "cloning_vectors", "title": "Cloning vectors", "layout": "scene", "scene": VECTOR_TABLE, "tags": ["techniques"],
     "notes": "載體至少需複製起點、選擇標記、多選殖位點。容量：plasmid < λ < cosmid < BAC < YAC。Shuttle vector 兩種起點；expression vector 帶啟動子與 RBS；"
              "真核表現另有 CMV promoter（哺乳類），E. coli 載體（如 pET）則不需要。"},
    {"id": "techniques_by_purpose", "title": "Molecular techniques by purpose", "layout": "scene", "scene": TECH_TABLE, "tags": ["techniques"],
     "notes": "記法：Southern-DNA、Northern-RNA、Western-蛋白質（SNoW DRoP）。交互作用：DNA–蛋白（EMSA、footprint、ChIP、one-hybrid）、蛋白–蛋白（co-IP、Y2H、FRET）。"
              "基因功能：knock-in（同源重組）、knock-out（Cre–loxP、CRISPR）、knock-down（RNAi）。"},
    {"id": "xgal_blue_dye", "title": "X-gal → blue dye", "layout": "free", "cell_size": 13.0,
     "nodes": [{"id": "ind", "ref": "bromochloroindoxyl", "pos": [0, 0]}, {"id": "dye", "ref": "dibromodichloroindigo", "pos": [1.2, 0]}],
     "edges": [], "text": False, "tags": ["techniques"],
     "notes": "β-galactosidase 水解 X-gal 得到無色 indoxyl；兩分子 indoxyl 經空氣氧化二聚成不溶的藍色靛藍染料（5,5′-dibromo-4,4′-dichloro-indigo）。"},
]


# ================================================================ writer

def main(data):
    data = Path(data)
    for m in MOLECULES:
        d = {"id": m["id"], "name": m["name"], "label": m["label"], "aliases": m["aliases"], "smiles": m["smiles"], "ring": m["ring"],
             "tags": with_parent(m["tags"]), "notes_zh": m["notes"], "parts_png": False, "status": "pending"}
        dump(data / "molecules" / f"{m['id']}.yaml", d)
    for r in REACTIONS:
        d = {"id": r["id"], "title": r["title"], "substrates": r["substrates"], "products": r["products"],
             "cofactors_in": r["cofactors_in"], "cofactors_out": r["cofactors_out"], "enzyme": r["enzyme"], "reversible": False}
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
