"""Round 19: molecular biology II — DNA damage, repair and recombination (illustrated scenes, tables,
damage structures and the MGMT reaction).

    python tools/gen_mb_repair.py data
"""
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent))
from figure_text import tidy  # noqa: E402
from reaction_colors import COLORS, frag, hl  # noqa: E402

M, C = "magenta", "cyan"
PATCH = "#2E9E4F"        # newly synthesised repair patch
W = 0.3                  # strand width in repair diagrams (no bases drawn: like the lecture notes)
GAP = 1.0                # distance between the two strands of a duplex


def seg(x0, x1, y, kind="parent", ends=None):
    d = {"g": "strand", "kind": kind, "points": [[x0, y], [x1, y]], "width": W}
    if ends:
        d["ends"] = list(ends)
    return d


def duplex(y, x0, x1, top_gaps=(), bot_gaps=(), top="parent", bot="new", ends=True, patches=()):
    """Two strands (top 5'→3' left to right); *_gaps = [(a, b)] holes; patches = [(strand 'top'|'bot', a, b)]."""
    out = []
    for yy, gaps, kind, e in ((y, top_gaps, top, ("5'", "3'")), (y + GAP, bot_gaps, bot, ("3'", "5'"))):
        cuts = sorted(gaps)
        xs = [x0] + [v for g in cuts for v in g] + [x1]
        pieces = [(xs[i], xs[i + 1]) for i in range(0, len(xs), 2)]
        for k, (a, b) in enumerate(pieces):
            if b - a < 0.05:
                continue
            ends_k = None
            if ends:                       # True | "left" | "right": which outer ends get 5'/3'
                left = ends in (True, "left") and k == 0
                right = ends in (True, "right") and k == len(pieces) - 1
                ends_k = (e[0] if left else "", e[1] if right else "")
            out.append(seg(a, b, yy, kind, ends_k))
    for which, a, b in patches:
        out.append(seg(a, b, y if which == "top" else y + GAP, PATCH))
    return out


def txt(x, y, text, size=0.85, color=None, align="center"):
    d = {"g": "text", "at": [x, y], "text": text, "size": size, "align": align}
    if color:
        d["color"] = color
    return d


def step(y, n):
    """Step number at the left margin + down arrow to the next row."""
    return [txt(-1.6, y + 0.5, f"{n}.", 0.95, "gray")]


def down(x, y0, y1):
    return {"g": "arrow", "from": [x, y0], "to": [x, y1], "color": "gray", "width": 1.0}


def note(y, text, x=0.0, size=0.78):
    """Explanatory paragraph: kept out of the figure (user: short labels only) and moved into notes_zh by with_notes()."""
    return {**txt(x, y, text, size, None, "left"), "_note": True}


def with_notes(d):
    """Strip note() paragraphs from a pathway's scene and append their text to notes_zh (【補充】)."""
    sc = d.get("scene")
    if not sc:
        return d
    extra = [" ".join(it["text"].split()) for it in sc["items"] if it.get("_note")]
    d["scene"] = {**sc, "items": [it for it in sc["items"] if not it.get("_note")]}
    if extra:
        d["notes_zh"] = d.get("notes_zh", "") + "\n【補充】" + " ".join(extra)
    return d


# ================================================================ mismatch repair (lecture notes p.30)
Y = [1.5, 6.0, 10.5, 14.5]
MMR = {"items": [
    *step(Y[0], 1), *duplex(Y[0], 0, 18),
    {"g": "mark", "kind": "methyl", "at": [3.0, Y[0]], "side": "up"},
    txt(3.0, Y[0] - 1.75, "GATC (hemimethylated)", 0.72, "gray"),
    txt(12.0, Y[0] - 0.55, "G", 0.8), txt(12.0, Y[0] + GAP + 0.8, "T", 0.8, "red"),
    {"id": "muts", "g": "blob", "at": [12.0, Y[0] + 0.5], "w": 1.6, "h": 1.0, "shape": "kidney", "hue": "blue",
     "label": "MutS: finds the mismatch", "label_align": "left", "label_at": [19.4, Y[0] - 1.2]},
    {"id": "mutl", "g": "blob", "at": [7.5, Y[0] + 0.5], "w": 2.2, "h": 1.0, "shape": "box", "hue": "teal",
     "label": "MutL: links MutS–MutH", "label_align": "left", "label_at": [19.4, Y[0] + 0.5], "point": [8.6, Y[0] + 0.5]},
    {"id": "muth", "g": "blob", "at": [3.0, Y[0] + 0.5], "w": 1.4, "h": 1.2, "shape": "round", "hue": "orange",
     "label": "MutH nicks the\nunmethylated (new) strand", "label_align": "left", "label_at": [19.4, Y[0] + 2.5], "point": [3.6, Y[0] + 1.0]},
    {"g": "mark", "kind": "cut", "at": [3.7, Y[0] + GAP], "side": "down"},
    down(16.5, Y[0] + 2.6, Y[1] - 0.9),
    *step(Y[1], 2), *duplex(Y[1], 0, 18, bot_gaps=[(3.7, 13.2)]),
    {"g": "mark", "kind": "methyl", "at": [3.0, Y[1]], "side": "up"},
    {"g": "strand", "kind": "gray", "points": [[4.4, Y[1] + 2.2], [12.6, Y[1] + 2.2]], "width": 0.22},
    txt(8.5, Y[1] + 3.0, "excised: UvrD (helicase II) + exonuclease, SSB", 0.75, "gray"),
    down(16.5, Y[1] + 3.5, Y[2] - 0.9),
    *step(Y[2], 3), *duplex(Y[2], 0, 18, patches=[("bot", 3.7, 13.2)]),
    {"g": "mark", "kind": "methyl", "at": [3.0, Y[2]], "side": "up"},
    txt(12.0, Y[2] + GAP + 0.6, "C", 0.8, PATCH),
    {"id": "pol3", "g": "blob", "at": [13.4, Y[2] + 0.5], "w": 1.8, "h": 1.4, "shape": "hand", "hue": "blue",
     "label": "DNA Pol III fills the gap;\nDNA ligase (NAD^+) seals", "label_align": "left", "label_at": [19.4, Y[2] + 0.5]},
    down(16.5, Y[2] + 2.0, Y[3] - 0.9),
    *step(Y[3], 4), *duplex(Y[3], 0, 18),
    {"g": "mark", "kind": "methyl", "at": [3.0, Y[3]], "side": "up"},
    {"g": "mark", "kind": "methyl", "at": [3.0, Y[3] + GAP], "side": "down"},
    txt(12.0, Y[3] - 0.55, "G", 0.8), txt(12.0, Y[3] + GAP + 0.6, "C", 0.8, PATCH),
    txt(19.4, Y[3] + 0.5, "Dam methylase (SAM) methylates\nthe new GATC → fully methylated", 0.8, None, "left"),
    note(Y[3] + 3.4, "Eukaryotes: MSH2–MSH6 (MutSα) / MSH2–MSH3, MLH1–PMS2 (MutLα, endonuclease); strand signal = nicks / Okazaki ends (no MutH).\n"
         "Defect → Lynch syndrome (HNPCC): microsatellite instability, colorectal / endometrial cancer."),
]}

# ================================================================ base-excision repair
YB = [1.5, 5.6, 9.7, 13.8]
BER = {"items": [
    *step(YB[0], 1), *duplex(YB[0], 0, 16),
    txt(8.0, YB[0] + GAP + 0.6, "U", 0.85, "red"),
    {"id": "udg", "g": "blob", "at": [8.0, YB[0] + 0.45], "w": 1.8, "h": 1.4, "shape": "kidney", "hue": "orange",
     "label": "DNA glycosylase (e.g. uracil-DNA glycosylase)\ncleaves the N-glycosidic bond", "label_align": "left", "label_at": [17.4, YB[0] + 0.5]},
    down(8.0, YB[0] + 2.4, YB[1] - 0.9),
    *step(YB[1], 2), *duplex(YB[1], 0, 16),
    txt(8.0, YB[1] + GAP + 0.6, "AP site", 0.75, "red"),
    {"g": "mark", "kind": "cut", "at": [7.7, YB[1] + GAP], "side": "up"},
    {"id": "ape", "g": "blob", "at": [7.0, YB[1] + 0.45], "w": 1.4, "h": 1.0, "shape": "round", "hue": "purple",
     "label": "AP endonuclease (APE1)\nnicks 5' of the AP site", "label_align": "left", "label_at": [17.4, YB[1] + 0.5]},
    down(8.0, YB[1] + 2.4, YB[2] - 0.9),
    *step(YB[2], 3), *duplex(YB[2], 0, 16, patches=[("bot", 7.7, 8.3)]),
    {"id": "polb", "g": "blob", "at": [8.9, YB[2] + 0.45], "w": 1.8, "h": 1.4, "shape": "hand", "hue": "blue", "rot": 180,
     "label": "Pol β (eukaryotes) / Pol I (E. coli):\ninserts 1 nt, removes dRP", "label_align": "left", "label_at": [17.4, YB[2] + 0.5]},
    down(8.0, YB[2] + 2.4, YB[3] - 0.9),
    *step(YB[3], 4), *duplex(YB[3], 0, 16, patches=[("bot", 7.7, 8.3)]),
    txt(8.0, YB[3] + GAP + 0.6, "C", 0.85, PATCH),
    txt(17.4, YB[3] + 0.5, "DNA ligase III–XRCC1 seals", 0.8, None, "left"),
    note(YB[3] + 3.2, "Small, non-distorting lesions: deaminated bases (U from C, hypoxanthine from A), 8-oxoG (OGG1), alkylated bases, AP sites (depurination).\n"
         "5-methyl-C → T gives a T·G mismatch: removed by TDG (thymine-DNA glycosylase) — CpG sites are mutation hot spots."),
]}

# ================================================================ nucleotide-excision repair
YN = [1.5, 5.8, 10.1, 14.4]
NER = {"items": [
    *step(YN[0], 1), *duplex(YN[0], 0, 18),
    {"g": "mark", "kind": "lesion", "at": [9.0, YN[0] - 0.05]},
    {"id": "uvra", "g": "blob", "at": [8.0, YN[0] + 0.5], "w": 2.0, "h": 1.5, "shape": "bilobe", "hue": "green",
     "label": "UvrA_2B: scans, recognises\nthe distortion (ATP)", "label_align": "left", "label_at": [19.4, YN[0] + 0.5], "point": [7.5, YN[0] + 1.1]},
    txt(9.0, YN[0] - 1.1, "T=T dimer (UV), bulky adduct", 0.72, "red"),
    down(16.5, YN[0] + 2.4, YN[1] - 0.9),
    *step(YN[1], 2), *duplex(YN[1], 0, 18),
    {"g": "mark", "kind": "lesion", "at": [9.0, YN[1] - 0.05]},
    {"g": "mark", "kind": "cut", "at": [6.0, YN[1]], "side": "up"}, {"g": "mark", "kind": "cut", "at": [11.0, YN[1]], "side": "up"},
    {"id": "uvrc", "g": "blob", "at": [9.0, YN[1] + 0.5], "w": 1.6, "h": 1.2, "shape": "kidney", "hue": "pink",
     "label": "UvrC: dual incision\n(8 nt 5', 4–5 nt 3' of the lesion)", "label_align": "left", "label_at": [19.4, YN[1] + 0.5], "point": [9.6, YN[1] + 0.8]},
    down(16.5, YN[1] + 2.4, YN[2] - 0.9),
    *step(YN[2], 3), *duplex(YN[2], 0, 18, top_gaps=[(6.0, 11.0)]),
    {"g": "strand", "kind": "parent", "points": [[6.3, YN[2] - 1.0], [10.7, YN[2] - 1.0]], "width": W},
    {"g": "mark", "kind": "lesion", "at": [9.0, YN[2] - 1.05], "r": 0.35},
    txt(19.4, YN[2] + 0.3, "UvrD (helicase II) releases the\n12–13-mer oligonucleotide", 0.8, None, "left"),
    down(16.5, YN[2] + 2.4, YN[3] - 0.9),
    *step(YN[3], 4), *duplex(YN[3], 0, 18, patches=[("top", 6.0, 11.0)]),
    txt(19.4, YN[3] + 0.5, "DNA Pol I fills, ligase seals", 0.8, None, "left"),
    note(YN[3] + 3.2, "Eukaryotes: XPC (global genome) or stalled RNA Pol II + CSA/CSB (transcription-coupled) recognise the lesion; TFIIH helicases XPB/XPD open it;\n"
         "XPG (3') and XPF–ERCC1 (5') cut out 24–32 nt; Pol δ/ε fill. Defects: xeroderma pigmentosum (UV-sensitive skin cancer), Cockayne syndrome."),
]}

# ================================================================ homologous recombination (RecBCD, notes p.32)
def recA(x0, x1, y):
    n = int((x1 - x0) / 0.62)
    return [{"g": "blob", "at": [x0 + 0.31 + 0.62 * k, y], "w": 0.6, "h": 0.6, "shape": "round", "hue": "orange"} for k in range(n)]


H1, H2, H3, H4, H5 = 1.5, 8.4, 14.4, 20.4, 25.9
HR = {"items": [
    *step(H1, 1), *duplex(H1, 0, 7.5, top="parent", bot="parent", ends="left"), *duplex(H1, 8.5, 16, top="parent", bot="parent", ends="right"),
    *duplex(H1 + 2.6, 0, 16, top="new", bot="new"),
    txt(8.0, H1 - 1.0, "double-strand break", 0.75, "red"),
    {"id": "recbcd", "g": "blob", "at": [7.3, H1 + 0.5], "w": 1.7, "h": 1.7, "shape": "bilobe", "hue": "purple",
     "label": "RecBCD (helicase + nuclease)\nresects until chi (GCTGGTGG)", "label_at": [19.4, H1 + 0.5], "label_align": "left"},
    txt(8.0, H1 + 4.6, "homologous sister duplex", 0.75, "gray"),
    down(14.0, H1 + 5.0, H2 - 1.8),
    *step(H2, 2), seg(0, 7.0, H2), seg(0, 4.0, H2 + GAP), *duplex(H2 + 2.6, 0, 16, top="new", bot="new"),
    *recA(4.2, 7.0, H2 - 0.05),
    {"id": "reca", "g": "text", "at": [5.6, H2 - 1.1], "text": "RecA filament on the 3' overhang", "size": 0.75},
    txt(19.4, H2 + 0.5, "3' single-strand tail coated by RecA\n(eukaryotes: Rad51, loaded with BRCA2)", 0.8, None, "left"),
    down(14.0, H2 + 4.4, H3 - 1.0),
    *step(H3, 3), seg(0, 5.0, H3), seg(0, 4.0, H3 + GAP),
    {"g": "strand", "kind": "new", "points": [[0, H3 + 2.6], [5.4, H3 + 2.6], [6.4, H3 + 1.7], [11.6, H3 + 1.7], [12.6, H3 + 2.6], [16, H3 + 2.6]], "width": W},
    {"g": "strand", "kind": "new", "points": [[0, H3 + 3.6], [16, H3 + 3.6]], "width": W},
    {"g": "strand", "kind": "parent", "points": [[5.0, H3], [6.0, H3 + 2.4], [10.0, H3 + 2.6]], "width": W, "arrow": True},
    txt(19.4, H3 + 2.0, "strand invasion → D-loop;\nthe invading 3' end is extended", 0.8, None, "left"),
    down(14.0, H3 + 4.3, H4 - 1.2),
    *step(H4, 4),
    {"g": "strand", "kind": "parent", "points": [[0, H4], [7.0, H4], [8.0, H4 + 1.3], [16, H4 + 1.3]], "width": W},
    {"g": "strand", "kind": "new", "points": [[0, H4 + 1.3], [7.0, H4 + 1.3], [8.0, H4], [16, H4]], "width": W},
    seg(0, 16, H4 - 1.0), {"g": "strand", "kind": "new", "points": [[0, H4 + 2.3], [16, H4 + 2.3]], "width": W},
    {"id": "ruvab", "g": "blob", "at": [7.5, H4 + 0.65], "w": 1.1, "h": 1.1, "shape": "round", "hue": "teal",
     "label": "Holliday junction: RuvA (binds)\n+ RuvB (branch migration, ATP)", "label_at": [19.4, H4 + 0.65], "label_align": "left"},
    down(14.0, H4 + 3.2, H5 - 1.0),
    *step(H5, 5), seg(0, 8.0, H5), seg(8.0, 16, H5, "new"), seg(0, 8.0, H5 + GAP), seg(8.0, 16, H5 + GAP, "new"),
    seg(0, 8.0, H5 + 2.6, "new"), seg(8.0, 16, H5 + 2.6), seg(0, 8.0, H5 + 3.6, "new"), seg(8.0, 16, H5 + 3.6),
    txt(19.4, H5 + 1.8, "RuvC (resolvase) cuts the junction\n→ crossover (or non-crossover) products", 0.8, None, "left"),
    note(H5 + 5.4, "Error-free repair of DSBs and stalled forks (S/G2: needs the sister chromatid). Eukaryotes: MRN (Mre11–Rad50–Nbs1) + CtIP resect, BRCA1, BRCA2–Rad51.\n"
         "Defects: BRCA1/2 (breast / ovarian cancer; PARP inhibitors), ataxia telangiectasia (ATM), Bloom (BLM helicase), Werner (WRN)."),
]}

# ================================================================ NHEJ (notes p.35)
NJ = 2.0
NHEJ = {"items": [
    *duplex(NJ, 0, 5.4, top="parent", bot="parent", ends="left"), *duplex(NJ, 6.6, 12, top="parent", bot="parent", ends="right"),
    txt(6.0, NJ - 1.2, "DSB (radiation, ROS)", 0.75, "red"),
    {"g": "arrow", "from": [13.0, NJ + 0.5], "to": [15.0, NJ + 0.5], "style": "fat"},
    *duplex(NJ, 16, 21.0, top="parent", bot="parent", ends="left"), *duplex(NJ, 23.0, 28, top="parent", bot="parent", ends="right"),
    {"g": "ring", "at": [20.6, NJ + 0.5], "w": 0.9, "h": 2.4, "hue": "dkgreen"},
    {"id": "ku", "g": "ring", "at": [23.4, NJ + 0.5], "w": 0.9, "h": 2.4, "hue": "dkgreen",
     "label": "Ku70/80 rings\nbind the ends", "label_at": [22.0, NJ + 4.4], "point": [22.4, NJ + 1.5]},
    {"id": "pkcs", "g": "blob", "at": [22.0, NJ - 1.4], "w": 3.4, "h": 1.3, "shape": "box", "hue": "yellow",
     "label": "DNA-PKcs bridges the ends", "label_at": [22.0, NJ - 3.4]},
    {"g": "arrow", "from": [22.0, NJ + 6.0], "to": [22.0, NJ + 7.6], "style": "fat"},
    *duplex(NJ + 9.0, 16, 21.4, top="parent", bot="parent", ends="left"), *duplex(NJ + 9.0, 22.6, 28, top="parent", bot="parent", ends="right"),
    {"id": "artemis", "g": "blob", "at": [20.2, NJ + 8.0], "w": 1.3, "h": 1.0, "shape": "kidney", "hue": "red",
     "label": "Artemis trims", "label_at": [17.0, NJ + 7.0]},
    {"id": "polmu", "g": "blob", "at": [22.0, NJ + 10.9], "w": 1.6, "h": 1.0, "shape": "hand", "hue": "blue",
     "label": "Pol μ / λ fill", "label_at": [18.0, NJ + 12.6]},
    {"id": "lig4", "g": "blob", "at": [24.0, NJ + 10.9], "w": 1.8, "h": 1.0, "shape": "box", "hue": "teal",
     "label": "Ligase IV–XRCC4", "label_at": [26.8, NJ + 12.6]},
    {"g": "arrow", "from": [15.0, NJ + 9.5], "to": [13.0, NJ + 9.5], "style": "fat"},
    *duplex(NJ + 9.0, 0, 12, top="parent", bot="parent"),
    {"g": "mark", "kind": "lesion", "at": [6.0, NJ + 8.6], "r": 0.3},
    txt(6.0, NJ + 7.2, "joined — a few nt may be lost\n(error-prone)", 0.75),
    note(NJ + 14.6, "NHEJ: main DSB pathway in G1 and in mammals; no homologous template. Also joins RAG-cut ends in V(D)J recombination and class switching.\n"
         "Defects in Artemis, DNA-PKcs or ligase IV → SCID (radiosensitive)."),
]}

# ================================================================ SOS response and translesion synthesis (notes p.33)
SOS = {"items": [
    *duplex(2.0, 0, 14, top="parent", bot="parent", ends=False),
    {"g": "region", "x0": 6.0, "y0": 1.4, "x1": 9.0, "y1": 3.6, "hue": "orange", "opacity": 0.45, "outline": False},
    {"id": "lexa", "g": "blob", "at": [7.5, 1.0], "w": 2.4, "h": 1.3, "shape": "bilobe", "hue": "purple",
     "label": "LexA dimer (repressor)\non SOS boxes", "label_at": [3.0, -1.2], "point": [7.0, 0.7]},
    {"g": "blob", "at": [11.6, 1.0], "w": 2.6, "h": 1.8, "shape": "oval", "hue": "gray"},
    {"g": "mark", "kind": "x", "at": [10.0, 0.4]},
    txt(12.2, -1.3, "RNA polymerase blocked", 0.75, "gray"),
    txt(7.0, 5.0, "SOS genes OFF:\nrecA, lexA, uvrA/B, umuDC (Pol V), dinB (Pol IV), sulA", 0.75, "gray"),
    {"g": "arrow", "from": [15.5, 2.5], "to": [17.5, 2.5], "style": "fat"},
    txt(16.5, 1.0, "DNA damage", 0.8, "red"),
    {"g": "strand", "kind": "parent", "points": [[19.0, 0.6], [27.0, 0.6]], "width": W},
    *recA(19.2, 27.0, 0.55),
    txt(23.0, -0.8, "ssDNA + RecA → RecA* (activated)", 0.8),
    {"g": "blob", "at": [21.4, 3.2], "w": 1.2, "h": 1.2, "shape": "round", "hue": "purple"},
    {"g": "blob", "at": [23.4, 3.2], "w": 1.2, "h": 1.2, "shape": "round", "hue": "purple"},
    {"g": "mark", "kind": "cut", "at": [22.4, 3.2], "side": "down"},
    txt(23.0, 5.0, "RecA* acts as a co-protease:\nLexA self-cleaves → SOS genes ON", 0.8),
    # translesion synthesis
    # new strand on top (5'→3', growing rightwards towards the lesion); the damaged template is the bottom strand
    *duplex(9.0, 0, 13, top="new", bot="parent", top_gaps=[(7.0, 13.0)]),
    {"g": "mark", "kind": "lesion", "at": [7.6, 10.0]},
    {"g": "blob", "at": [6.4, 9.0], "w": 1.8, "h": 1.4, "shape": "hand", "hue": "blue",
     "label": "Pol III stalls", "label_at": [3.0, 12.4]},
    {"g": "arrow", "from": [14.0, 9.5], "to": [16.0, 9.5], "style": "fat"},
    txt(15.0, 8.0, "polymerase\nswitch", 0.72, "gray"),
    *duplex(9.0, 17, 30, top="new", bot="parent", top_gaps=[(25.0, 30.0)]),
    {"g": "mark", "kind": "lesion", "at": [24.6, 10.0]},
    {"g": "blob", "at": [24.6, 9.0], "w": 2.0, "h": 1.5, "shape": "hand", "hue": "orange",
     "label": "Pol V (UmuD'_2C) / Pol IV inserts across\nthe lesion: no proofreading → mutations", "label_at": [23.5, 12.8]},
    note(16.2, "SOS = last-resort, error-prone survival response. Eukaryotic TLS polymerases: Pol η (accurate across T=T; defect → XP variant), Pol ζ, Rev1, Pol κ, Pol ι.\n"
         "Stalled forks can also restart by recombination (Pol II, RecA; \"chicken-foot\" fork regression)."),
]}

# ================================================================ site-specific recombination (notes p.37)
SITE = {"items": [
    txt(0.0, -0.6, "Integration / excision (direct repeats) — phage λ", 0.95, None, "left"),
    {"g": "plasmid", "at": [4.0, 3.5], "kind": "relaxed", "r": 1.8},
    {"g": "region", "x0": 3.3, "y0": 5.0, "x1": 4.7, "y1": 5.8, "hue": "orange", "opacity": 0.9, "outline": False},
    txt(4.0, 6.6, "attP (λ)", 0.75),
    *duplex(8.5, 0, 9.0, top="parent", bot="parent", ends=False),
    {"g": "region", "x0": 3.9, "y0": 8.3, "x1": 5.1, "y1": 9.7, "hue": "orange", "opacity": 0.9, "outline": False},
    txt(4.5, 10.6, "attB (E. coli)", 0.75),
    {"g": "arrow", "from": [10.0, 5.0], "to": [12.6, 5.0], "style": "fat"},
    txt(11.3, 3.4, "Int + IHF", 0.75, "#2E6E9E"),
    {"g": "arrow", "from": [12.6, 7.0], "to": [10.0, 7.0], "style": "fat"},
    txt(11.3, 8.6, "Int + Xis\n(+ IHF)", 0.75, "#2E6E9E"),
    *duplex(5.5, 13.5, 30.0, top="parent", bot="parent", ends=False),
    {"g": "region", "x0": 15.5, "y0": 5.3, "x1": 16.7, "y1": 6.7, "hue": "orange", "opacity": 0.9, "outline": False},
    {"g": "region", "x0": 16.7, "y0": 5.3, "x1": 26.8, "y1": 6.7, "hue": "teal", "opacity": 0.35, "outline": False},
    {"g": "region", "x0": 26.8, "y0": 5.3, "x1": 28.0, "y1": 6.7, "hue": "orange", "opacity": 0.9, "outline": False},
    txt(16.1, 4.4, "attL", 0.75), txt(27.4, 4.4, "attR", 0.75), txt(21.8, 7.6, "prophage (lysogeny)", 0.8),
    txt(0.0, 13.0, "Inversion (inverted repeats) — e.g. Salmonella flagellar phase variation (Hin), Cre–loxP in opposite orientation", 0.95, None, "left"),
    *duplex(15.5, 0, 12, top="parent", bot="parent", ends=False),
    {"g": "ligand", "at": [3.0, 16.0], "shape": "triangle", "hue": "red", "r": 0.55},
    {"g": "region", "x0": 4.0, "y0": 15.3, "x1": 8.0, "y1": 16.7, "hue": "green", "opacity": 0.7, "outline": False},
    txt(6.0, 16.0, "A → B", 0.75),
    {"g": "ligand", "at": [9.0, 16.0], "shape": "triangle", "hue": "red", "r": 0.55},
    {"g": "arrow", "from": [13.0, 16.0], "to": [15.6, 16.0], "style": "fat"},
    txt(14.3, 14.5, "recombinase", 0.75, "#2E6E9E"),
    *duplex(15.5, 16.6, 28.6, top="parent", bot="parent", ends=False),
    {"g": "ligand", "at": [19.6, 16.0], "shape": "triangle", "hue": "red", "r": 0.55},
    {"g": "region", "x0": 20.6, "y0": 15.3, "x1": 24.6, "y1": 16.7, "hue": "green", "opacity": 0.7, "outline": False},
    txt(22.6, 16.0, "B ← A", 0.75),
    {"g": "ligand", "at": [25.6, 16.0], "shape": "triangle", "hue": "red", "r": 0.55},
    note(19.4, "Recombinase (Tyr or Ser family) cuts and rejoins at short specific sites, no DNA synthesis, no RecA.\n"
         "Direct repeats → excision / integration; inverted repeats → inversion. Cre–loxP (P1 phage) is used for conditional knock-outs."),
]}

# ================================================================ transposition (notes p.39–40)
def tn(x0, x1, y, fill="pink", label=None, ir=True):
    out = [{"g": "region", "x0": x0, "y0": y - 0.2, "x1": x1, "y1": y + GAP + 0.2, "hue": fill, "opacity": 0.75, "outline": False}]
    if ir:
        out += [{"g": "ligand", "at": [x0 + 0.35, y + 0.5], "shape": "triangle", "hue": "purple", "r": 0.35},
                {"g": "ligand", "at": [x1 - 0.35, y + 0.5], "shape": "triangle", "hue": "purple", "r": 0.35}]
    if label:
        out.append(txt((x0 + x1) / 2, y + 0.5, label, 0.7))
    return out


def tsd(x, y):
    return {"g": "region", "x0": x, "y0": y - 0.2, "x1": x + 0.8, "y1": y + GAP + 0.2, "hue": "yellow", "opacity": 1.0, "outline": False}


TRANSPOSON = {"items": [
    txt(0.0, -0.8, "Bacterial transposons (inverted repeats ▲ at both ends)", 0.95, None, "left"),
    *duplex(0.5, 0, 9, top="parent", bot="parent", ends=False), *tn(2.0, 7.0, 0.5, "pink", "transposase"),
    txt(4.5, 2.6, "IS element (insertion sequence)", 0.75),
    *duplex(0.5, 10.5, 29, top="parent", bot="parent", ends=False), *tn(11.5, 16.0, 0.5, "pink", "IS"),
    *tn(16.0, 23.5, 0.5, "orange", "drug resistance", ir=False), *tn(23.5, 28.0, 0.5, "pink", "IS"),
    txt(19.8, 2.6, "composite transposon (e.g. Tn10, Tn5)", 0.75),
    txt(0.0, 4.8, "Two modes of transposition", 0.95, None, "left"),
    {"g": "plasmid", "at": [3.0, 8.6], "kind": "relaxed", "r": 1.6},
    {"g": "region", "x0": 2.2, "y0": 9.9, "x1": 3.8, "y1": 10.6, "hue": "pink", "opacity": 0.9, "outline": False},
    txt(3.0, 11.4, "donor", 0.75),
    *duplex(8.0, 6.0, 12.5, top="parent", bot="parent", ends=False),
    txt(9.2, 10.1, "target", 0.75),
    {"g": "arrow", "from": [13.4, 7.0], "to": [15.6, 6.0], "style": "fat"},
    {"g": "arrow", "from": [13.4, 10.0], "to": [15.6, 11.0], "style": "fat"},
    *duplex(5.5, 16.4, 27, top="parent", bot="parent", ends=False), tsd(18.4, 5.5), *tn(19.2, 23.6, 5.5, "pink"), tsd(23.6, 5.5),
    {"g": "plasmid", "at": [29.6, 6.0], "kind": "relaxed", "r": 1.2}, {"g": "mark", "kind": "x", "at": [29.6, 4.8]},
    txt(22.0, 3.9, "cut-and-paste (non-replicative): donor left broken", 0.75),
    *duplex(11.0, 16.4, 27, top="parent", bot="parent", ends=False), tsd(18.4, 11.0), *tn(19.2, 23.6, 11.0, "pink"), tsd(23.6, 11.0),
    {"g": "plasmid", "at": [29.6, 11.5], "kind": "relaxed", "r": 1.2},
    {"g": "region", "x0": 28.9, "y0": 12.5, "x1": 30.3, "y1": 13.1, "hue": "pink", "opacity": 0.9, "outline": False},
    txt(22.0, 13.6, "copy-and-paste (replicative): donor keeps a copy", 0.75),
    txt(21.0, 15.0, "yellow = target-site duplication (short direct repeats)", 0.72, "gray"),
    note(16.8, "Retrotransposons (eukaryotes) copy via RNA → reverse transcriptase → cDNA: LTR type (retrovirus-like) and non-LTR LINEs (L1, autonomous, encode RT)\n"
         "and SINEs (Alu ~300 bp, non-autonomous, use L1's RT). ~45 % of the human genome is transposon-derived; only ~1–2 % codes for protein."),
]}

# ================================================================ V(D)J recombination (notes p.38)
def gene(x0, x1, y, hue, text):
    return [{"g": "region", "x0": x0, "y0": y - 0.2, "x1": x1, "y1": y + GAP + 0.2, "hue": hue, "opacity": 0.85, "outline": False},
            txt((x0 + x1) / 2, y + 0.5, text, 0.68)]


VY = [1.0, 6.0, 11.0]
VDJ = {"items": [
    txt(-0.2, VY[0] - 1.3, "Germline (heavy-chain locus)", 0.85, "gray", "left"),
    *duplex(VY[0], 0, 31, top="parent", bot="parent", ends=False),
    *gene(1.0, 2.6, VY[0], "blue", "V1"), *gene(3.4, 5.0, VY[0], "blue", "V2"), *gene(5.8, 7.4, VY[0], "blue", "Vn"),
    *gene(10.0, 11.0, VY[0], "green", "D1"), *gene(11.6, 12.6, VY[0], "green", "D2"), *gene(13.2, 14.2, VY[0], "green", "Dn"),
    *gene(17.0, 18.0, VY[0], "orange", "J1"), *gene(18.6, 19.6, VY[0], "orange", "J2"), *gene(20.2, 21.2, VY[0], "orange", "J6"),
    *gene(24.0, 29.0, VY[0], "purple", "C (μ, δ, γ …)"),
    {"g": "arrow", "from": [15.0, VY[0] + 2.0], "to": [15.0, VY[1] - 0.8], "color": "gray"},
    txt(15.4, VY[0] + 2.9, "D–J joining (RAG1/2 cut at RSS, NHEJ joins)", 0.75, None, "left"),
    *duplex(VY[1], 0, 31, top="parent", bot="parent", ends=False),
    *gene(1.0, 2.6, VY[1], "blue", "V1"), *gene(3.4, 5.0, VY[1], "blue", "V2"), *gene(5.8, 7.4, VY[1], "blue", "Vn"),
    # downstream J segments (J3–J6) stay in the DNA and are removed later by RNA splicing
    *gene(11.0, 12.0, VY[1], "green", "D2"), *gene(12.0, 13.0, VY[1], "orange", "J2"), *gene(14.6, 15.6, VY[1], "orange", "J6"),
    *gene(20.0, 25.0, VY[1], "purple", "C"),
    {"g": "arrow", "from": [15.0, VY[1] + 2.0], "to": [15.0, VY[2] - 0.8], "color": "gray"},
    txt(15.4, VY[1] + 2.9, "V–DJ joining; intervening DNA lost as a circle", 0.75, None, "left"),
    *duplex(VY[2], 0, 31, top="parent", bot="parent", ends=False),
    *gene(1.0, 2.6, VY[2], "blue", "V1"), *gene(3.4, 4.4, VY[2], "blue", "V2"), *gene(4.4, 5.4, VY[2], "green", "D2"),
    *gene(5.4, 6.4, VY[2], "orange", "J2"), *gene(8.0, 9.0, VY[2], "orange", "J6"), *gene(13.0, 18.0, VY[2], "purple", "C"),
    {"g": "plasmid", "at": [26.0, VY[2] + 0.5], "kind": "relaxed", "r": 1.3},
    txt(26.0, VY[2] + 2.8, "excised signal circle", 0.72, "gray"),
    txt(9.6, VY[2] + 2.6, "transcribed: V–D–J exon + C → one heavy chain", 0.75),
    txt(-0.2, 16.3, "RSS = heptamer – spacer (12 or 23 bp) – nonamer; 12/23 rule: a 12-RSS pairs only with a 23-RSS. RAG1/2 make hairpin coding ends,\n"
        "Artemis opens them, TdT adds N nucleotides → junctional diversity. RAG or Artemis defect → SCID (no B and T cells).", 0.75, None, "left"),
]}

# ================================================================ tables
def table(rows, w_name=6.2, x_desc=6.8, hue="blue", head_cols=None):
    items, y = [], 0.0
    if head_cols:
        for x, t in head_cols:
            items.append({"g": "text", "at": [x, y], "text": t, "size": 0.9, "align": "left", "color": "gray"})
        y += 1.4
    for name, desc in rows:
        items.append({"g": "blob", "at": [w_name / 2, y], "w": w_name, "h": 1.3, "shape": "box", "hue": hue,
                      "text": name, "text_size": 0.85})
        items.append({"g": "text", "at": [x_desc, y], "text": desc, "size": 0.82, "align": "left"})
        y += 1.6 + 0.62 * desc.count("\n")
    return {"items": items}


REPAIR_TABLE = table([
    ("Mismatch repair", "replication errors, small loops · MutS/L/H (E. coli), MSH2/MLH1 · Lynch syndrome (HNPCC)"),
    ("Base excision", "deaminated / oxidised / alkylated single bases, AP sites · glycosylase → APE1 → Pol β → ligase III"),
    ("Nucleotide excision", "bulky, helix-distorting: T=T dimers, benzo[a]pyrene adducts · UvrABC; XPA–XPG, TFIIH · xeroderma pigmentosum"),
    ("Direct reversal", "photolyase splits T=T with light (not in humans); MGMT removes O^6-methyl-G (suicide enzyme); AlkB demethylates"),
    ("Homologous recomb.", "double-strand breaks, stalled forks (S/G2) · RecBCD, RecA; MRN, BRCA1/2, Rad51 · error-free"),
    ("NHEJ", "double-strand breaks (G1) · Ku70/80, DNA-PKcs, Artemis, ligase IV · error-prone; also V(D)J"),
    ("Translesion (SOS)", "lesion blocks the replicative polymerase · Pol IV, V (E. coli); Pol η, ζ · error-prone"),
], w_name=8.4, x_desc=9.0)

MUTATION_TABLE = table([
    ("Transition", "purine ↔ purine or pyrimidine ↔ pyrimidine (A·T ↔ G·C); e.g. deamination, tautomers, 5-BU"),
    ("Transversion", "purine ↔ pyrimidine (e.g. 8-oxoG pairs with A → G·C to T·A)"),
    ("Silent", "codon changes, same amino acid (often the 3rd / wobble position)"),
    ("Missense", "one amino acid replaced (e.g. HbS: β6 Glu → Val)"),
    ("Nonsense", "codon → stop (UAA, UAG, UGA): truncated protein, nonsense-mediated decay"),
    ("Frameshift", "insertion / deletion not a multiple of 3 (e.g. intercalators: acridine, ethidium)"),
    ("Repeat expansion", "trinucleotide repeats grow by slippage: Huntington (CAG), fragile X (CGG), myotonic dystrophy (CTG)"),
], w_name=7.4, x_desc=8.0, hue="purple")


# ================================================================ damaged bases and MGMT
MOLECULES = [
    {"id": "thymine_dimer", "name": "Cyclobutane thymine dimer", "label": "Thymine dimer (CPD)", "aliases": ["T=T", "CPD", "胸腺嘧啶雙體"],
     "smiles": "[17*]N1C(=O)NC(=O)C2(C)C1C1N([17*])C(=O)NC(=O)C12C", "ring": {"coordgen": True, "rotate": 90, "methyls": False},
     "tags": ["repair"], "notes": "UV（UVB）使相鄰兩個嘧啶的 C5–C6 雙鍵形成環丁烷（cis-syn CPD，約 80 %）；另一種是 6-4 photoproduct（約 20 %）。"
              "扭曲雙螺旋、阻擋複製與轉錄。修復：NER（人類唯一途徑；XP 缺陷）；細菌與植物另有 photolyase（FADH⁻ 吸光，直接切開）。"},
    {"id": "oxoguanine", "name": "8-Oxoguanine", "label": "8-oxoG", "aliases": ["8-oxo-dG", "8-oxoG"],
     "smiles": "[17*]n1c2nc(N)[nH]c(=O)c2[nH]c1=O", "ring": {"coordgen": True},
     "tags": ["repair"], "notes": "ROS 氧化 G 的 C8。syn 構形可與 A 配對 → G·C 變 T·A（transversion）。由 OGG1（glycosylase）以 BER 移除；MUTYH 移除對面誤配的 A。"},
    {"id": "o6_methylguanine", "name": "O6-Methylguanine", "label": "O^6-Methylguanine", "aliases": ["O6-meG"],
     "smiles": "[17*]n1cnc2c(OC)nc(N)nc21", "ring": {"coordgen": True},
     "tags": ["repair"], "notes": "烷化劑（nitrosamine、temozolomide）在 G 的 O6 加甲基 → 與 T 配對 → G·C 變 A·T（transition）。由 MGMT 直接修復。"},
    {"id": "guanine_dr", "name": "Guanine (in DNA)", "label": "Guanine (in DNA)", "aliases": [],
     "smiles": "[17*]n1cnc2c(=O)[nH]c(N)nc21", "ring": {"coordgen": True},
     "tags": ["repair"], "notes": "DNA 中的 guanine（dR = 去氧核糖與骨架）。"},
]

REACTIONS = [
    {"id": "mgmt", "title": "O6-Methylguanine → Guanine (MGMT, direct repair)",
     "substrates": ["o6_methylguanine"], "products": ["guanine_dr"], "cofactors_in": ["MGMT-SH"], "cofactors_out": ["MGMT-SCH3"],
     "enzyme": "O^6-methylguanine-DNA methyltransferase (MGMT)", "reversible": False,
     "card": {"type": "直接修復：O6 上的甲基（洋紅）轉移到酵素自身的 Cys–S；N1 由 Cys 的 H 補上",
              "mechanism": "Suicide enzyme：每個 MGMT 只能用一次（甲基化後被泛素化分解），不需要切除或合成 DNA",
              "exam_points": ["MGMT 啟動子甲基化（基因沉默）的膠質母細胞瘤對 temozolomide 較敏感", "E. coli 對應：Ada 蛋白（亦為轉錄活化子）",
                              "其他直接修復：photolyase（T=T，人類沒有）、AlkB 家族（ALKBH2/3，氧化去甲基）"]},
     "tags": ["repair"]},
]

COLORS.update({
    "mgmt": {"o6_methylguanine": hl("[CH3:1]Oc"), "MGMT-SCH3": frag((M, "CH_3")),
             "MGMT-SH": frag((C, "H")), "guanine_dr": hl("[nH:1]1c(N)nc2c(c1=O)nc[n]2", col=C)},
})

NEW_COFACTORS = {
    "MGMT-SH": {"label": "MGMT–Cys–SH", "smiles": "S[18*]"},
    "MGMT-SCH3": {"label": "MGMT–Cys–S–CH_3", "smiles": "CS[18*]"},
}

# ================================================================ pathways

PATHWAYS = [
    {"id": "mismatch_repair", "title": "Mismatch repair (E. coli)", "layout": "scene", "scene": MMR, "tags": ["repair"],
     "notes": "複製錯誤（錯配、小插入／缺失）。E. coli 以 Dam 甲基化區分新舊股：新合成股的 GATC 暫時未甲基化（半甲基化）。"
              "MutS 認錯配 → MutL 連接 → MutH 切開最近的半甲基化 GATC 的未甲基化（新）股 → UvrD + 外切酶從切口移除到錯配之後 → Pol III 補、ligase 封口 → Dam 甲基化新股。"
              "真核：MSH2–MSH6、MLH1–PMS2，以切口辨識新股；缺陷 → Lynch syndrome（HNPCC，微衛星不穩定）。"},
    {"id": "base_excision_repair", "title": "Base-excision repair", "layout": "scene", "scene": BER, "tags": ["repair"],
     "notes": "單一受損鹼基（不扭曲雙螺旋）：DNA glycosylase 切斷 N-醣苷鍵 → AP（無嘌呤／無嘧啶）位 → APE1 在 5′ 端切開 → Pol β 補 1 nt 並移除 dRP → ligase III–XRCC1。"
              "Glycosylase 例：UNG（uracil，來自 C 去胺）、OGG1（8-oxoG）、TDG（T·G，來自 5mC 去胺）、AAG（烷化）。自發去嘌呤每細胞每天約上萬次。"
              "DNA 用 T 而不用 U：C 去胺成 U 才能被辨認為錯誤。"},
    {"id": "nucleotide_excision_repair", "title": "Nucleotide-excision repair", "layout": "scene", "scene": NER, "tags": ["repair"],
     "notes": "扭曲雙螺旋的大型損傷（UV 嘧啶雙體、6-4 光產物、苯并芘加成物、cisplatin 交聯）。E. coli：UvrA₂B 偵測 → UvrC 雙切（5′ 端 8 nt、3′ 端 4–5 nt）→ UvrD 移走 12–13-mer → Pol I、ligase。"
              "真核：GG-NER（XPC）與 TC-NER（停住的 Pol II + CSA/CSB）；TFIIH（XPB、XPD 解旋）打開，XPG 與 XPF–ERCC1 切除 24–32 nt。"
              "XP（xeroderma pigmentosum）：日光敏感、皮膚癌；Cockayne syndrome：TC-NER 缺陷。TFIIH 同時參與轉錄（CTD kinase 只用於轉錄）。"},
    {"id": "homologous_recombination", "title": "Homologous recombination repair (RecBCD → Holliday junction)", "layout": "scene", "scene": HR,
     "tags": ["repair", "recombination"],
     "notes": "雙股斷裂（或停住的複製叉）以同源姊妹染色分體為模板，無誤修復。E. coli：RecBCD（helicase + nuclease）處理斷端直到 chi 序列 → 留下 3′ 單股尾端 → RecA 包覆並搜尋同源 → 股入侵形成 D-loop → 延伸 → Holliday junction"
              " → RuvA（認 HJ）＋RuvB（ATP，branch migration）→ RuvC（resolvase）切開。真核：MRN、CtIP、BRCA1、BRCA2–Rad51。"
              "RecA 也是 SOS 的共蛋白酶。減數分裂交換也走 HR。"},
    {"id": "nhej", "title": "Non-homologous end joining (NHEJ)", "layout": "scene", "scene": NHEJ, "tags": ["repair", "recombination"],
     "notes": "雙股斷裂直接接回，不需同源模板（G1 期、哺乳類主要途徑），可能丟失幾個核苷酸 → 易錯。"
              "Ku70/80 環套在斷端 → 招募 DNA-PKcs（kinase，與 Ku 合稱 DNA-PK）→ Artemis（核酸酶）修整 → Pol μ/λ 補 → ligase IV–XRCC4（–XLF）連接。"
              "也用於 V(D)J 重組與類別轉換；缺陷 → 放射線敏感型 SCID。"},
    {"id": "sos_response", "title": "SOS response and translesion synthesis", "layout": "scene", "scene": SOS, "tags": ["repair"],
     "notes": "平時 LexA 二聚體壓住 SOS 基因。DNA 損傷產生單股 DNA → RecA 結合成 RecA*（活化）→ 作為共蛋白酶促使 LexA 自我切割 → SOS 基因表現："
              "recA、lexA、uvrA/B（NER）、sulA（暫停分裂）、Pol II、Pol IV（dinB）、Pol V（UmuD′₂C；活性形式為 UmuD′₂C–RecA–ATP，UmuD 需 RecA* 協助自切成 UmuD′）。"
              "TLS：Pol III 在損傷前停住 → 換成 Pol IV/V 越過損傷（無校對，易錯）→ 再換回 Pol III。真核 TLS：Pol η（可正確越過 T=T；缺陷 → XP variant）、Pol ζ、Rev1。"},
    {"id": "site_specific_recombination", "title": "Site-specific recombination", "layout": "scene", "scene": SITE, "tags": ["recombination"],
     "notes": "Recombinase 在短的特定序列切開再接回，不需 DNA 合成與 RecA。直接重複 → 整合或切除（λ phage：attP × attB 由 integrase + IHF 整合成溶源，attL/attR；切除需 Xis）；"
              "反向重複 → 倒轉（Salmonella 鞭毛相變異 Hin）。Cre–loxP（P1 phage）用於條件式基因剔除。"},
    {"id": "transposition", "title": "Transposable elements", "layout": "scene", "scene": TRANSPOSON, "tags": ["recombination"],
     "notes": "Transposon（跳躍基因、mobile element）。細菌：IS（只有 transposase，兩端反向重複）、composite（IS–抗藥基因–IS，Tn10、Tn5）、Tn3（transposase + resolvase + 抗藥基因）、phage Mu。"
              "非複製型（cut-and-paste）：轉位子離開原處；複製型（copy-and-paste）：原處保留一份。插入處兩側產生短的直接重複（target-site duplication）。"
              "真核以 retrotransposon 為主：RNA → 反轉錄酶 → cDNA 再插入；LTR 型（似反轉錄病毒）、LINE（L1，自主，有 RT）、SINE（Alu，約 300 bp，需借 L1 的 RT）。"},
    {"id": "vdj_recombination", "title": "V(D)J recombination", "layout": "scene", "scene": VDJ, "tags": ["recombination"],
     "notes": "B、T 細胞發育時，抗體重鏈先 D–J、再 V–DJ 接合（輕鏈只有 V–J）。RAG1/2 在 RSS（heptamer–12 或 23 bp spacer–nonamer）切開，遵守 12/23 規則；"
              "編碼端形成髮夾，由 Artemis 打開，TdT 加 N 核苷酸，NHEJ 接合 → 接合多樣性。中間 DNA 以環狀切除（signal joint）。下游未用到的 J（如 J3–J6）留在 DNA 上，轉錄後由 RNA splicing 去除。RAG、Artemis 缺陷 → SCID。"},
    {"id": "dna_repair_pathways", "title": "DNA repair pathways at a glance", "layout": "scene", "scene": REPAIR_TABLE, "tags": ["repair"],
     "notes": "依損傷選擇修復途徑：錯配 → MMR；單一鹼基 → BER；大型扭曲 → NER；甲基 → MGMT（直接）；雙股斷裂 → HR（S/G2，無誤）或 NHEJ（G1，易錯）；無法修復時以 TLS 越過。"
              "相關疾病：XP（NER）、Lynch（MMR）、BRCA1/2（HR）、ataxia telangiectasia（ATM）、Bloom（BLM）、Werner（WRN）、Fanconi anemia（交聯修復）。"},
    {"id": "mutation_types", "title": "Types of point mutation", "layout": "scene", "scene": MUTATION_TABLE, "tags": ["repair"],
     "notes": "鹼基替換：transition（同類嘌呤或嘧啶互換）、transversion（嘌呤↔嘧啶）。對蛋白質的影響：silent、missense、nonsense；插入／缺失非 3 的倍數 → frameshift。"
              "三核苷酸重複擴增（滑動）：Huntington（CAG，多麩醯胺）、fragile X（CGG，FMR1 甲基化沉默）、myotonic dystrophy（CTG），具 anticipation。"},
    {"id": "dna_damage_structures", "title": "Damaged bases", "layout": "free", "cell_size": 12.0,
     "nodes": [{"id": "td", "ref": "thymine_dimer", "pos": [0, 0], "decorate": {"highlight": [{"smarts": "C1(C)C(N)C(N)C1(C)", "color": "magenta"}]}},
               {"id": "og", "ref": "oxoguanine", "pos": [1.05, 0], "decorate": {"highlight": [{"smarts": "[nH]c(=O)n", "color": "magenta"}]}},
               {"id": "mg", "ref": "o6_methylguanine", "pos": [2.0, 0], "decorate": {"highlight": [{"smarts": "[CH3]O", "color": "magenta"}]}}],
     "edges": [], "text": False, "tags": ["repair"],
     "notes": "常見 DNA 損傷：UV → 環丁烷嘧啶雙體（NER；photolyase）；ROS → 8-oxoG（BER，OGG1）；烷化 → O⁶-methyl-G（MGMT）；"
              "自發 → 去嘌呤（AP 位）、去胺（C→U、A→hypoxanthine、G→xanthine、5mC→T）。洋紅為受損部位；dR = DNA 中的去氧核糖。"},
]


# ================================================================ writer

def dump(path, d):
    """Write the YAML, keeping a review status already given in the web UI (re-running a generator must not reset it)."""
    if path.exists() and "status" in d:
        old = (yaml.safe_load(path.read_text(encoding="utf-8")) or {}).get("status")
        if old in ("approved", "changes_requested"):
            d["status"] = old
    path.write_text(yaml.safe_dump(d, allow_unicode=True, sort_keys=False, width=200), encoding="utf-8")


PARENT = {t: "molecular_biology" for t in ("dna_structure", "replication", "repair", "recombination", "transcription",
                                           "rna_processing", "translation", "gene_regulation", "techniques")}


def with_parent(tags):
    out = []
    for x in tags:
        if x in PARENT:
            out.append(PARENT[x])
        out.append(x)
    return list(dict.fromkeys(out))


def add_cofactors(data):
    p = data / "cofactors.yaml"
    text = p.read_text(encoding="utf-8")
    new = [k for k in NEW_COFACTORS if f"\n{k}:" not in text]
    if new:
        lines = ["", "# --- 第十九輪（DNA 修復）：MGMT 的 Cys（[18*] = 蛋白質其餘部分；由 tools/gen_mb_repair.py 產生）"]
        for k in new:
            v = NEW_COFACTORS[k]
            lines.append(f'{k}: {{label: "{v["label"]}", smiles: "{v["smiles"]}"}}')
        p.write_text(text.rstrip("\n") + "\n" + "\n".join(lines) + "\n", encoding="utf-8")


def main(data):
    data = Path(data)
    add_cofactors(data)
    for m in MOLECULES:
        d = {"id": m["id"], "name": m["name"], "label": m["label"], "aliases": m["aliases"], "smiles": m["smiles"], "ring": m["ring"],
             "tags": with_parent(m["tags"]), "notes_zh": m["notes"], "parts_png": False, "status": "pending"}
        dump(data / "molecules" / f"{m['id']}.yaml", d)
    for r in REACTIONS:
        d = {"id": r["id"], "title": r["title"], "substrates": r["substrates"], "products": r["products"],
             "cofactors_in": r["cofactors_in"], "cofactors_out": r["cofactors_out"], "enzyme": r["enzyme"], "reversible": r["reversible"]}
        if COLORS.get(r["id"]):
            d["decorate"] = COLORS[r["id"]]
        d.update({"card_zh": r["card"], "tags": with_parent(r["tags"]), "parts_png": False, "status": "pending"})
        dump(data / "reactions" / f"{r['id']}.yaml", d)
    for p in PATHWAYS:
        d = {k: v for k, v in p.items() if k not in ("tags", "notes")}
        d["tags"] = with_parent(p["tags"])
        d["notes_zh"] = p["notes"]
        d = tidy(with_notes(d))
        d["parts_png"] = False
        d["status"] = "pending"
        dump(data / "pathways" / f"{p['id']}.yaml", d)
    print(f"{len(MOLECULES)} molecules, {len(REACTIONS)} reactions, {len(PATHWAYS)} pathways")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data")
