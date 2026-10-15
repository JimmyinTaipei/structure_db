"""Round 18: molecular biology I — DNA structure and replication (illustrated `scene` figures).

    python tools/gen_mb_rep.py data
"""
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent))
from figure_text import tidy  # noqa: E402

L = 11.8          # lagging new-strand row (y)
T_LAG = 12.8      # lagging template row
T_LEAD = 5.2      # leading template row
N_LEAD = 6.2      # leading new-strand row


def frag(x0, x1, primer=None, arrow=False, ends=("", "")):
    """Okazaki fragment drawn 5′ (right) → 3′ (left); primer = (x0, x1) RNA at its 5′ end."""
    out = []
    if primer:
        out.append({"g": "strand", "kind": "primer", "points": [[primer[1], L], [primer[0], L]], "bases": "left",
                    "base_len": 0.45, "seq": "GCAU"})
    out.append({"g": "strand", "kind": "new", "points": [[x1, L], [x0, L]], "bases": "left", "base_len": 0.45,
                "arrow": arrow, "ends": list(ends), "seq": "TACGGATC"})
    return out


REPLICATION_FORK = {
    "size": [34, 15],
    "items": [
        # parental duplex ahead of the fork + gyrase
        {"g": "dsdna", "from": [23, 9], "to": [33.5, 9], "width": 1.3, "seq": "GATCCGTAAGCT"},
        # the top strand continues the leading-strand template (3' at the far left), so its right end is 5'
        {"g": "text", "at": [34.4, 8.35], "text": "5'", "size": 0.62, "color": "#3F4F6B"},
        {"g": "text", "at": [34.4, 9.65], "text": "3'", "size": 0.62, "color": "#3F4F6B"},
        {"id": "gyrase", "g": "blob", "at": [29, 9], "w": 1.9, "h": 2.6, "shape": "oval", "hue": "purple",
         "label": "Topoisomerase II\n(gyrase): relieves\nsupercoils", "label_at": [29, 3.4]},
        # template arms
        {"g": "strand", "kind": "parent", "points": [[23, 8.35], [18, T_LEAD], [1, T_LEAD]], "bases": "left",
         "base_len": 0.45, "ends": ["", "3'"], "seq": "CTAGGCATTCGA"},
        {"g": "strand", "kind": "parent", "points": [[23, 9.65], [18, T_LAG], [1, T_LAG]], "bases": "right",
         "base_len": 0.45, "ends": ["", "5'"], "seq": "ATGCCTAGGATC"},
        # leading strand
        {"g": "strand", "kind": "primer", "points": [[1.2, N_LEAD], [2.2, N_LEAD]], "bases": "left", "base_len": 0.45, "seq": "GAUC"},
        {"g": "strand", "kind": "new", "points": [[2.2, N_LEAD], [16.6, N_LEAD]], "bases": "left", "base_len": 0.45,
         "arrow": True, "ends": ["", ""], "seq": "GATCCGTAAGCT"},
        {"g": "text", "at": [0.25, N_LEAD], "text": "5'", "size": 0.62, "color": "#C2338A"},
        {"id": "clamp_lead", "g": "ring", "at": [16.2, 5.7], "w": 0.8, "h": 2.4, "hue": "green"},
        {"id": "pol_lead", "g": "blob", "at": [18.0, 6.3], "w": 2.6, "h": 2.1, "shape": "hand", "hue": "blue", "accent": True,
         "label": "DNA polymerase III\n+ β sliding clamp", "label_at": [17.6, 1.6], "point": [18.0, 5.4]},
        # helicase, SSB, primase on the lagging template
        {"id": "helicase", "g": "ring", "kind": "hexamer", "at": [21.7, 10.3], "r": 0.75, "hue": "orange", "accent": True,
         "label": "Helicase (DnaB)\nATP", "label_at": [25.5, 13.6]},
        {"id": "ssb1", "g": "blob", "at": [19.7, 11.7], "w": 0.9, "h": 0.75, "shape": "round", "hue": "yellow"},
        {"id": "ssb2", "g": "blob", "at": [18.9, 12.35], "w": 0.9, "h": 0.75, "shape": "round", "hue": "yellow",
         "label": "SSB", "label_at": [19.3, 14.6], "point": [19.2, 12.6]},
        {"id": "primase", "g": "blob", "at": [17.6, 11.3], "w": 1.6, "h": 1.1, "shape": "kidney", "hue": "pink",
         "label": "Primase\n(DnaG)", "label_at": [18.4, 8.6], "point": [17.6, 10.9], "label_size": 0.8},
        # lagging strand: oldest (left) → newest (right)
        *frag(1.0, 4.4, ends=("", "3'")),
        *frag(4.55, 8.6, primer=(8.6, 9.3)),
        *frag(9.75, 12.6, primer=(12.6, 13.4)),
        *frag(15.0, 16.0, primer=(16.0, 16.8), arrow=True),
        {"id": "ligase", "g": "blob", "at": [4.5, 10.85], "w": 1.3, "h": 0.95, "shape": "round", "hue": "teal",
         "label": "DNA ligase", "label_at": [3.6, 8.9], "point": [4.4, 10.5], "label_size": 0.8},
        {"id": "pol1", "g": "blob", "at": [9.5, 10.85], "w": 1.5, "h": 1.05, "shape": "kidney", "hue": "gray",
         "label": "Pol I removes\nRNA primer", "label_at": [8.9, 8.9], "point": [9.4, 10.4], "label_size": 0.8},
        {"id": "pol_lag", "g": "blob", "at": [13.9, 11.45], "w": 2.0, "h": 1.7, "shape": "hand", "rot": 180, "hue": "blue",
         "accent": True, "label": "Pol III\n+ clamp", "label_at": [13.8, 8.9], "point": [13.9, 10.7], "label_size": 0.8},
        {"id": "clamp_lag", "g": "ring", "at": [15.25, 11.7], "w": 0.6, "h": 1.8, "hue": "green"},
        # strand names and fork direction
        {"g": "text", "at": [7.5, 3.8], "text": "Leading strand (continuous)", "color": "#D9622B", "size": 0.95},
        {"g": "text", "at": [8, 14.6], "text": "Lagging strand: Okazaki fragments (1000–2000 nt)", "color": "#D9622B", "size": 0.95},
        {"g": "arrow", "from": [26.5, 15.8], "to": [32.5, 15.8], "style": "fat"},
        {"g": "text", "at": [29.5, 17.0], "text": "Fork movement", "size": 0.85, "color": "gray"},
    ],
}


def lbl(x, y, text, size=0.85, color=None, align="center"):
    d = {"g": "text", "at": [x, y], "text": text, "size": size, "align": align}
    if color:
        d["color"] = color
    return d


def fat(x0, y0, x1, y1):
    return {"g": "arrow", "from": [x0, y0], "to": [x1, y1], "style": "fat"}


# ---------------------------------------------------------------- B-DNA
P_TURN, AMP = 8.4, 1.7                       # S per turn (34 Å), half width (2 nm ↔ 3.4 S)
DOUBLE_HELIX = {"items": [
    {"id": "helix", "g": "helix", "from": [0, 5], "to": [25.2, 5], "pitch": P_TURN, "amp": AMP, "groove": 0.36},
    {"g": "text", "at": [13.0, 0.6], "text": "Minor groove", "size": 0.9},
    {"g": "arrow", "from": [13.0, 1.2], "to": [P_TURN * 1.43, 5 - AMP + 0.3], "color": "gray", "width": 0.8},
    {"g": "text", "at": [20.6, 0.6], "text": "Major groove", "size": 0.9},
    {"g": "arrow", "from": [20.6, 1.2], "to": [P_TURN * 1.93, 5 - AMP + 0.3], "color": "gray", "width": 0.8},
    {"g": "arrow", "from": [P_TURN * 0.25, 8.4], "to": [P_TURN * 1.25, 8.4], "color": "black", "width": 0.8},
    {"g": "arrow", "from": [P_TURN * 1.25, 8.4], "to": [P_TURN * 0.25, 8.4], "color": "black", "width": 0.8},
    lbl(P_TURN * 0.75, 9.3, "one turn = 3.4 nm (10.5 bp)", 0.85),
    {"g": "arrow", "from": [26.6, 5 - AMP], "to": [26.6, 5 + AMP], "color": "black", "width": 0.8},
    {"g": "arrow", "from": [26.6, 5 + AMP], "to": [26.6, 5 - AMP], "color": "black", "width": 0.8},
    lbl(27.2, 5, "2 nm", 0.85, align="left"),
    lbl(0.0, 11.3, "B-DNA: right-handed, antiparallel; bases stacked 0.34 nm apart, perpendicular to the axis;\n"
        "A–T 2 H-bonds, G≡C 3 H-bonds; proteins read the sequence mainly in the major groove.\n"
        "A-DNA: right-handed, wider, 11 bp/turn (dehydrated; RNA–DNA hybrids, dsRNA).   "
        "Z-DNA: left-handed zigzag, 12 bp/turn (alternating GC, methylated).", 0.8, align="left"),
]}

# ---------------------------------------------------------------- chromatin levels
def _beads(x0, y, n, gap=2.2):
    pts = [[x0 - 0.8, y + 0.3]]
    items = []
    for i in range(n):
        x = x0 + i * gap
        pts += [[x - 0.5, y + 0.15], [x + 0.5, y - 0.15]] if i % 2 == 0 else [[x - 0.5, y - 0.15], [x + 0.5, y + 0.15]]
        items.append({"g": "nucleosome", "at": [x, y - 0.1], "r": 0.8})
    pts.append([x0 + (n - 1) * gap + 1.0, y])
    return [{"g": "strand", "kind": "parent", "points": pts, "width": 0.2}] + items


CHROMATIN = {"items": [
    {"g": "helix", "from": [0, 4], "to": [5.5, 4], "pitch": 2.6, "amp": 0.45},
    lbl(2.75, 6.6, "DNA double helix\n2 nm"),
    fat(6.0, 4, 7.6, 4),
    *_beads(8.8, 4, 3),
    {"g": "text", "at": [11.0, 1.4], "text": "H1 binds linker DNA", "size": 0.75, "color": "gray"},
    lbl(11.0, 6.6, "\"Beads on a string\"\nnucleosome 11 nm"),
    fat(14.6, 4, 16.2, 4),
    *[{"g": "nucleosome", "at": [17.6 + 1.2 * (i % 2), 1.8 + 0.9 * i], "r": 0.7} for i in range(4)],
    *[{"g": "nucleosome", "at": [20.0 + 1.2 * (i % 2), 1.8 + 0.9 * i], "r": 0.7} for i in range(4)],
    lbl(19.4, 6.6, "30 nm fiber"),
    fat(22.6, 4, 24.2, 4),
    {"g": "blob", "at": [27.0, 4.0], "w": 0.7, "h": 5.0, "shape": "box", "hue": "brown"},
    *[{"g": "strand", "kind": "parent", "width": 0.28, "points": [[27.0, y], [25.0 + 4.0 * (k % 2), y - 0.6], [25.3 + 3.4 * (k % 2), y + 0.4], [27.0, y + 0.5]]}
      for k, y in enumerate((1.8, 2.9, 4.0, 5.1))],
    lbl(27.0, 7.4, "Loops on a protein\nscaffold (~300 nm)"),
    fat(30.0, 4, 31.6, 4),
    {"g": "chromosome", "at": [33.6, 4.0], "h": 5.6},
    lbl(33.6, 8.6, "Metaphase\nchromosome\n(~1400 nm)"),
    lbl(0.0, 10.8, "Nucleosome core: histone octamer (H2A, H2B, H3, H4) × 2 + ~147 bp DNA wrapped 1.7 turns; histones rich in Lys/Arg (+).\n"
        "Euchromatin: open, acetylated, active.  Heterochromatin: condensed, H3K9me3, inactive.", 0.8, align="left"),
]}

# ---------------------------------------------------------------- oriC and initiation
def _box(x0, x1, y, hue, text):
    return [{"g": "region", "x0": x0, "y0": y - 0.95, "x1": x1, "y1": y + 0.95, "hue": hue, "opacity": 0.55, "outline": False},
            ] + ([{"g": "text", "at": [(x0 + x1) / 2, y - 1.55], "text": text, "size": 0.7}] if text else [])


ORIC = {"items": [
    *_box(2.0, 3.6, 3, "orange", ""), *_box(4.0, 5.6, 3, "orange", ""), *_box(6.0, 7.6, 3, "orange", ""),
    *[b for k in range(5) for b in _box(11.0 + 2.2 * k, 12.4 + 2.2 * k, 3, "green", "")],
    lbl(4.8, 1.2, "3 × 13-mer", 0.75), lbl(15.4, 1.2, "DnaA boxes: 9-mer × 5", 0.75),
    {"g": "dsdna", "from": [0.5, 3], "to": [23.5, 3], "width": 1.1, "seq": "GATCAATTATAT"},
    lbl(4.8, 5.2, "DUE (AT-rich): melts first", 0.85), lbl(16.4, 5.2, "DnaA–ATP binding sites", 0.85),
    lbl(12.0, 17.0, "oriC (245 bp): 11 GATC sites — fully methylated (Dam) = ready; hemimethylated → SeqA blocks re-initiation",
        0.75, "gray"),
    # step 1: DnaA–ATP oligomer wraps DNA
    {"g": "helix", "from": [0.5, 11], "to": [8.0, 11], "pitch": 2.4, "amp": 0.45},
    *[{"g": "blob", "at": [3.0 + 1.0 * k, 10.0 + (0.0 if k % 2 else 0.25)], "w": 1.0, "h": 0.9, "shape": "round", "hue": "green",
       **({"label": "DnaA–ATP\noligomer", "label_at": [4.5, 7.9]} if k == 2 else {})} for k in range(4)],
    fat(9.0, 11, 10.8, 11),
    # step 2: DUE melts; DnaC loads DnaB on each single strand
    {"g": "strand", "kind": "parent", "points": [[11.5, 10.7], [14, 10.7], [16, 9.2], [19, 9.2], [21, 10.7], [24, 10.7]]},
    {"g": "strand", "kind": "parent", "points": [[11.5, 11.3], [14, 11.3], [16, 12.8], [19, 12.8], [21, 11.3], [24, 11.3]]},
    {"id": "dnab1", "g": "ring", "kind": "hexamer", "at": [16.6, 9.2], "r": 0.55, "hue": "orange",
     "label": "DnaB helicase\n(loaded by DnaC)", "label_at": [20.5, 7.4]},
    {"g": "ring", "kind": "hexamer", "at": [18.4, 12.8], "r": 0.55, "hue": "orange"},
    {"id": "ssb", "g": "blob", "at": [20.2, 12.2], "w": 0.8, "h": 0.65, "shape": "round", "hue": "yellow",
     "label": "SSB", "label_at": [22.8, 13.9]},
    {"id": "dnag", "g": "blob", "at": [15.0, 13.2], "w": 1.3, "h": 0.9, "shape": "kidney", "hue": "pink",
     "label": "DnaG primase", "label_at": [13.4, 15.0]},
    lbl(4.2, 13.8, "IHF, HU, FIS bend DNA;\ngyrase supplies\nnegative supercoils", 0.75, "gray"),
]}

# ---------------------------------------------------------------- bacterial chromosome cycle
CHROMOSOME_CYCLE = {"items": [
    {"g": "plasmid", "at": [2.0, 3], "kind": "relaxed", "r": 1.7},
    lbl(2.0, 6.0, "Relaxed circle"),
    fat(4.2, 3, 6.0, 3), lbl(5.1, 1.0, "Gyrase\n(ATP)", 0.75, "#5F4A8C"),
    {"g": "plasmid", "at": [7.6, 3], "kind": "supercoiled", "r": 1.8, "twists": 3},
    lbl(7.6, 7.4, "Negatively\nsupercoiled", 0.8),
    fat(9.3, 3, 11.1, 3), lbl(10.2, 1.0, "DnaA\nat oriC", 0.75, "#2F6A33"),
    {"g": "plasmid", "at": [13.4, 3], "kind": "theta", "r": 1.8, "bubble": [-70, 70]},
    lbl(13.4, 6.0, "θ replication:\n2 forks", 0.8),
    fat(15.7, 3, 17.5, 3), lbl(16.6, -0.2, "Ter / Tus\n(forks meet)", 0.75),
    {"g": "plasmid", "at": [20.4, 3], "kind": "catenane", "r": 1.5},
    lbl(20.4, 6.0, "Catenanes\n(interlinked)"),
    fat(23.0, 3, 24.8, 3), lbl(23.9, 1.0, "Topo IV", 0.75, "#5F4A8C"),
    {"g": "plasmid", "at": [26.8, 3], "kind": "relaxed", "r": 1.3},
    {"g": "plasmid", "at": [30.0, 3], "kind": "relaxed", "r": 1.3},
    lbl(28.4, 6.0, "2 daughter\nchromosomes"),
]}

# ---------------------------------------------------------------- topoisomerases
TOPO = {"items": [
    lbl(5.5, 0.4, "Topoisomerase I", 1.0),
    {"g": "strand", "kind": "parent", "points": [[0.5, 3.0], [4.8, 3.0]], "width": 0.28},
    {"g": "strand", "kind": "parent", "points": [[6.2, 3.0], [10.5, 3.0]], "width": 0.28},
    {"g": "strand", "kind": "parent", "points": [[0.5, 4.3], [10.5, 4.3]], "width": 0.28},
    {"id": "topo1", "g": "blob", "at": [5.5, 2.4], "w": 2.2, "h": 1.5, "shape": "kidney", "hue": "purple",
     "label": "Tyr–O–P (covalent)\ncuts ONE strand", "label_at": [5.5, 7.6], "point": [5.5, 3.0]},
    {"g": "arrow", "from": [7.4, 5.4], "to": [3.6, 5.4], "curve": 0.35, "color": "magenta"},
    lbl(5.5, 9.6, "strand rotates → ΔLk ±1, no ATP", 0.8),
    lbl(20.5, -0.6, "Topoisomerase II (gyrase in bacteria)", 1.0),
    {"g": "dsdna", "from": [14.5, 3.6], "to": [19.4, 3.6], "width": 0.9, "seq": "GC"},
    {"g": "dsdna", "from": [21.6, 3.6], "to": [26.5, 3.6], "width": 0.9, "seq": "GC"},
    {"g": "dsdna", "from": [20.5, 0.9], "to": [20.5, 7.8], "width": 0.9, "seq": "AT", "kinds": ["new", "new"]},
    {"id": "topo2a", "g": "blob", "at": [18.7, 3.6], "w": 1.8, "h": 1.6, "shape": "kidney", "hue": "purple", "rot": 90},
    {"id": "topo2b", "g": "blob", "at": [22.3, 3.6], "w": 1.8, "h": 1.6, "shape": "kidney", "hue": "purple", "rot": -90,
     "label": "cuts BOTH strands (G segment);\nT segment passes through", "label_at": [28.6, 6.6], "point": [22.9, 4.2]},
    lbl(20.5, 9.6, "ΔLk ±2, ATP; gyrase adds negative supercoils", 0.8),
    lbl(0.0, 11.6, "Inhibitors: camptothecin, irinotecan, topotecan (topo I); etoposide, doxorubicin (topo II);\n"
        "quinolones / fluoroquinolones (bacterial gyrase, topo IV); novobiocin (gyrase ATPase).", 0.8, align="left"),
]}

# ---------------------------------------------------------------- telomere
TEL_Y = 3.0
TELOMERE = {"items": [
    {"g": "strand", "kind": "parent", "points": [[0.5, TEL_Y], [13.0, TEL_Y]], "ends": ["", "3'"], "bases": "right", "base_len": 0.45,
     "seq": "TTAGGG"},
    {"g": "strand", "kind": "parent", "points": [[8.8, TEL_Y + 1.0], [0.5, TEL_Y + 1.0]], "ends": ["5'", ""], "bases": "right",
     "base_len": 0.45, "seq": "AATCCC"},
    lbl(6.0, 0.9, "G-rich strand (TTAGGG)_n, 3' overhang", 0.8),
    lbl(4.5, 6.0, "C-rich strand: lagging strand cannot be finished\n(last RNA primer removed) → ends shorten each division", 0.75, "gray"),
    fat(14.2, TEL_Y + 0.5, 16.0, TEL_Y + 0.5),
    {"g": "strand", "kind": "parent", "points": [[17.0, TEL_Y], [27.0, TEL_Y]], "bases": "right", "base_len": 0.45, "seq": "TTAGGG"},
    {"g": "strand", "kind": "new", "points": [[27.0, TEL_Y], [31.2, TEL_Y]], "bases": "right", "base_len": 0.45, "seq": "TTAGGG",
     "arrow": True, "ends": ["", "3'"]},
    {"g": "strand", "kind": "parent", "points": [[25.3, TEL_Y + 1.0], [17.0, TEL_Y + 1.0]], "bases": "right", "base_len": 0.45,
     "seq": "AATCCC"},
    {"id": "telomerase", "g": "blob", "at": [28.6, TEL_Y + 0.85], "w": 3.6, "h": 2.4, "shape": "oval", "hue": "teal", "accent": True,
     "label": "Telomerase (TERT reverse\ntranscriptase + RNA template)", "label_at": [26.0, 8.0]},
    {"g": "strand", "kind": "rna", "points": [[27.4, TEL_Y + 0.75], [30.0, TEL_Y + 0.75]], "width": 0.2, "bases": "left",
     "base_len": 0.35, "seq": "CCCAAU", "ends": ["3'", "5'"]},
    lbl(24.0, 0.9, "extends the 3' end using its own RNA template (3'-AAUCCC-5')", 0.8),
    lbl(0.0, 10.4, "Then Pol α–primase fills in the C-rich strand; the end folds into a T-loop (D-loop) or G-quadruplex, protected by shelterin (TRF1/2, POT1).\n"
        "Telomerase active in germ cells, stem cells and ~90 % of cancers; somatic cells lack it → senescence. Telomere 10–20 kb, overhang > 100 nt.", 0.75,
        align="left"),
]}

# ---------------------------------------------------------------- eukaryotic origin licensing
LICENSING = {"items": [
    lbl(6.5, 0.3, "G1: licensing (helicase loaded, inactive)", 0.95),
    {"g": "dsdna", "from": [0.5, 4], "to": [13, 4], "width": 1.0, "seq": "ATTTAT"},
    {"id": "orc", "g": "blob", "at": [3.2, 2.7], "w": 2.4, "h": 1.5, "shape": "kidney", "hue": "blue", "label": "ORC", "label_at": [3.2, 7.4]},
    {"id": "cdc6", "g": "blob", "at": [5.0, 2.6], "w": 1.2, "h": 1.0, "shape": "round", "hue": "teal", "label": "Cdc6, Cdt1",
     "label_at": [5.8, 9.0], "point": [5.0, 3.0]},
    {"g": "ring", "at": [8.6, 4], "w": 1.0, "h": 2.6, "hue": "orange"},
    {"g": "ring", "at": [9.7, 4], "w": 1.0, "h": 2.6, "hue": "orange"},
    {"id": "mcm", "g": "text", "at": [9.15, 7.4], "text": "MCM2–7\ndouble hexamer", "size": 0.85},
    fat(14.0, 4, 16.0, 4),
    lbl(23.5, 0.3, "S phase: firing (CDK + DDK)", 0.95),
    {"g": "dsdna", "from": [17.0, 4], "to": [19.6, 4], "width": 1.0, "seq": "ATTTAT"},
    {"g": "dsdna", "from": [27.4, 4], "to": [30.0, 4], "width": 1.0, "seq": "ATTTAT"},
    {"g": "strand", "kind": "parent", "points": [[19.6, 3.5], [21.0, 2.4], [26.0, 2.4], [27.4, 3.5]]},
    {"g": "strand", "kind": "parent", "points": [[19.6, 4.5], [21.0, 5.6], [26.0, 5.6], [27.4, 4.5]]},
    {"id": "cmg1", "g": "ring", "kind": "hexamer", "at": [20.6, 2.6], "r": 0.6, "hue": "orange"},
    {"id": "cmg2", "g": "ring", "kind": "hexamer", "at": [26.4, 5.4], "r": 0.6, "hue": "orange",
     "label": "CMG helicase\n(Cdc45–MCM–GINS)", "label_at": [26.6, 8.8]},
    lbl(23.5, 4.0, "bubble", 0.75, "gray"),
    lbl(0.0, 11.2, "Licensing only in G1 (low CDK); in S/G2 high CDK and geminin (blocks Cdt1) prevent reloading → each origin fires once per cycle.\n"
        "Eukaryotes: many origins (ARS in yeast ≈ 11 bp core), replicons; pols α (primase), ε (leading), δ (lagging); PCNA clamp, RFC clamp loader, primers removed by Pol δ strand displacement + FEN1 (Dna2 for long flaps); RNase H1 cannot cut the last rN.",
        0.75, align="left"),
]}

# ---------------------------------------------------------------- polymerase active site
POL_MECH = {"items": [
    {"id": "pol", "g": "blob", "at": [8.0, 5.0], "w": 11.0, "h": 8.0, "shape": "hand", "hue": "blue", "accent": False},
    {"g": "strand", "kind": "parent", "points": [[0.0, 6.1], [10.6, 6.1]], "bases": "left", "base_len": 0.5, "ends": ["3'", "5'"], "seq": "TACGGC"},
    {"g": "strand", "kind": "new", "points": [[0.0, 4.6], [8.2, 4.6]], "bases": "right", "base_len": 0.5, "ends": ["5'", ""], "seq": "ATGCCG"},
    {"id": "dntp", "g": "text", "at": [14.3, 5.0], "text": "dNTP: base–sugar–{P}α–{P}β–{P}γ", "size": 0.8, "align": "left"},
    {"g": "arrow", "from": [9.3, 4.2], "to": [15.6, 4.4], "curve": -0.18, "color": "magenta"},
    lbl(14.3, 2.8, "3'-OH attacks Pα", 0.75, "magenta", "left"),
    {"g": "ion", "at": [10.4, 5.0], "text": "Mg^{2+}", "hue": "green", "r": 0.5, "text_size": 0.4},
    {"g": "ion", "at": [11.6, 5.0], "text": "Mg^{2+}", "hue": "green", "r": 0.5, "text_size": 0.4},
    {"g": "arrow", "from": [19.2, 5.7], "to": [20.4, 7.4], "color": "gray"},
    lbl(20.6, 8.4, "PP_i → 2 P_i\n(drives reaction)", 0.75, "gray"),
    lbl(0.0, 10.6, "Palm (catalytic Asp + 2 Mg^{2+}) · fingers (bind incoming dNTP) · thumb (holds DNA)", 0.8, align="left"),
    lbl(4.0, 0.4, "Polymerase site: 5'→3'", 0.85, "#2E6E9E"),
    lbl(0.0, 12.2, "New strand grows 5'→3' only; template read 3'→5'. A primer with a free 3'-OH is required (no de novo start; primase/RNA pol can start).\n"
        "Proofreading: mismatch → 3' end shifts to the 3'→5' exonuclease site (Pol I, III, δ, ε; absent in Pol α, β, RT, Taq). Error ~10^{-9} with MMR.",
        0.75, align="left"),
]}

# ---------------------------------------------------------------- polymerase table
_ROWS = [("E. coli", None), ("Pol I", "removes RNA primers (5'→3' exo), fills gaps, repair; Klenow fragment lacks 5'→3' exo"),
         ("Pol III", "main replicase (α polymerase, ε 3'→5' exo, β clamp): fast, highly processive"),
         ("Pol II, IV, V", "repair / SOS translesion synthesis (IV, V error-prone, no proofreading)"),
         ("Eukaryotes", None), ("Pol α", "with primase: makes the RNA–DNA primer (no proofreading)"),
         ("Pol ε", "leading strand; proofreading"), ("Pol δ", "lagging strand; proofreading; PCNA clamp"),
         ("Pol β", "base-excision repair gap filling (no 3'→5' exo)"), ("Pol γ", "mitochondrial DNA replication")]


def _table():
    items, y = [], 0.0
    for name, desc in _ROWS:
        if desc is None:
            y += 0.6
            items.append({"g": "text", "at": [0, y], "text": name, "size": 1.0, "align": "left", "color": "gray"})
            y += 1.5
            continue
        items.append({"g": "blob", "at": [2.2, y], "w": 4.4, "h": 1.3, "shape": "box", "hue": "blue", "text": name, "text_size": 0.9})
        items.append({"g": "text", "at": [5.0, y], "text": desc, "size": 0.85, "align": "left"})
        y += 1.6
    return {"items": items}


POL_TABLE = _table()

# ---------------------------------------------------------------- chain terminators (structures)
ATP = "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O"
DATP = ATP.replace("[C@@H](O)[C@H]1O", "[C@@H](O)C1")
DDATP = ATP.replace("[C@@H](O)[C@H]1O", "CC1")
DTMP = "Cc1cn([C@@H]2O[C@H](COP(=O)([O-])[O-])[C@@H](O)C2)c(=O)[nH]c1=O"
AZT = DTMP.replace("COP(=O)([O-])[O-]", "CO").replace("[C@@H](O)C2", "[C@@H](N=[N+]=[N-])C2")
TP = "{P}–{P}–{P}–OCH_2"
ADE = "[*]n1cnc2c(N)ncnc21"


def hw(c3, c4=TP, base=ADE, right=5):
    return {"ring": "furanose", "numbers": True, "base": base, "base_right": right,
            "subs": {1: {"up": "base", "down": "H"}, 2: {"up": "H", "down": "H"}, 3: c3, 4: {"up": c4, "down": "H"}}}


MOLECULES = [
    {"id": "datp", "name": "dATP", "aliases": ["deoxyadenosine triphosphate"], "smiles": DATP, "haworth": hw({"up": "H", "down": "OH"}),
     "tags": ["replication", "nucleotide"], "notes": "DNA polymerase 的受質之一；3′-OH 讓下一個核苷酸可以接上。"},
    {"id": "ddatp", "name": "ddATP", "aliases": ["dideoxy ATP", "雙去氧核苷酸"], "smiles": DDATP, "haworth": hw({"up": "H", "down": "H"}),
     "tags": ["replication", "techniques"], "notes": "2′,3′-dideoxy：沒有 3′-OH → 接上後鏈就停止。Sanger 定序用螢光標記的 ddNTP。"},
    {"id": "azt", "name": "AZT (zidovudine)", "label": "AZT (zidovudine)", "aliases": ["zidovudine", "azidothymidine"], "smiles": AZT,
     "haworth": hw({"up": "H", "down": "N_3"}, "HOCH_2", "[*]n1cc(C)c(=O)[nH]c1=O", 8),
     "tags": ["replication"], "notes": "3′-azido-thymidine：細胞內磷酸化成 AZT-TP，被 HIV 反轉錄酶接上後因沒有 3′-OH 而終止鏈；反轉錄酶比人類 DNA pol 更容易誤用它。"},
    {"id": "acyclovir", "name": "Acyclovir", "aliases": ["aciclovir"], "smiles": "Nc1nc2c(ncn2COCCO)c(=O)[nH]1",
     "ring": {"down": 6},
     "tags": ["replication"], "notes": "無環的 guanosine 類似物（沒有完整的糖、沒有 3′-OH）。只有皰疹病毒的 thymidine kinase 能做第一次磷酸化 → 選擇性；acyclo-GTP 抑制病毒 DNA pol 並終止鏈。"},
]

PATHWAYS = [
    {"id": "dna_ligase_mechanism", "title": "DNA ligase: sealing a nick in three steps", "layout": "linear",
     "nodes": [
         {"id": "s0", "text": "Ligase–Lys–NH_2\n+ ATP (eukaryotes)\nor NAD^+ (E. coli)", "kind": "protein"},
         {"id": "s1", "text": "Ligase–Lys–AMP\n(adenylylated enzyme)", "kind": "protein"},
         {"id": "s2", "text": "AMP–{P}–5' of the nick\n(activated 5'-phosphate)", "kind": "name"},
         {"id": "s3", "text": "3'-OH attacks →\nphosphodiester bond", "kind": "name"},
     ],
     "edges": [{"from": "s0", "to": "s1", "label": "1. PP_i (or NMN) out"},
               {"from": "s1", "to": "s2", "label": "2. AMP moves\nto 5'-P"},
               {"from": "s2", "to": "s3", "label": "3. AMP out"}],
     "tags": ["replication"],
     "notes": "Ligase 封閉 nick：需要相鄰的 3′-OH 與 5′-磷酸（雙股 DNA 中）。① 酵素的 Lys 與 ATP（真核、T4）或 NAD⁺（E. coli）反應形成 Lys–AMP，放出 PPi（或 NMN）；"
              "② AMP 轉到 nick 的 5′-磷酸上（5′-5′ 焦磷酸鍵，活化）；③ 3′-OH 攻擊活化的磷酸形成磷酸二酯鍵，放出 AMP。"
              "用途：延遲股 Okazaki 片段、修復最後一步、重組 DNA（T4 ligase 可接平端）。"},
    {"id": "chain_terminators", "title": "Chain terminators lack a 3'-OH", "layout": "free", "cell_size": 14.5,
     "nodes": [{"id": "d", "ref": "datp", "pos": [0, 0], "decorate": {"colors": {"d3": "cyan"}}},
               {"id": "dd", "ref": "ddatp", "pos": [1, 0], "decorate": {"colors": {"d3": "magenta"}}},
               {"id": "z", "ref": "azt", "pos": [2, 0], "decorate": {"colors": {"d3": "magenta"}}},
               {"id": "a", "ref": "acyclovir", "pos": [2.85, 0], "decorate": {"highlight": [{"smarts": "n[CH2]OCC[OH]", "color": "magenta"}]}}],
     "edges": [], "text": False,
     "tags": ["replication", "techniques"],
     "notes": "DNA 聚合需要引子末端的 3′-OH 攻擊下一個 dNTP。沒有 3′-OH 的類似物一旦接上就終止延長："
              "ddNTP（Sanger 定序）、AZT（HIV 反轉錄酶）、acyclovir（皰疹，需病毒 TK 活化）、ganciclovir（CMV）。"
              "對照：cytarabine（ara-C，阿拉伯糖）、gemcitabine 屬於抗代謝物；5-FU、methotrexate 抑制 dTMP 合成。"},

    {"id": "dna_double_helix", "title": "B-DNA double helix", "layout": "scene", "scene": DOUBLE_HELIX,
     "tags": ["dna_structure"],
     "notes": "Watson–Crick B-DNA：右手、反平行雙股；糖–磷酸骨架在外、鹼基在內堆疊（stacking 是主要穩定力，H-bond 給專一性）。"
              "每圈 3.4 nm、10.5 bp，直徑 2 nm；大溝與小溝交替，序列專一的 DNA 結合蛋白（HTH、zinc finger、leucine zipper）多讀大溝。"
              "Chargaff：A = T、G = C（雙股）；GC 越多 T_m 越高。A-DNA 較寬短（11 bp/turn，RNA–DNA hybrid、dsRNA）；Z-DNA 左手鋸齒（12 bp/turn，嘌呤嘧啶交替、甲基化 CpG）。"},
    {"id": "chromatin_levels", "title": "Levels of chromatin packing", "layout": "scene", "scene": CHROMATIN,
     "tags": ["dna_structure"],
     "notes": "DNA（2 nm）→ 核小體（11 nm，「串珠」）→ 30 nm 纖維（H1 協助）→ 掛在蛋白質骨架上的環（~300 nm）→ 中期染色體（~1400 nm）。"
              "核小體：H2A、H2B、H3、H4 各兩個組成八聚體，纏繞 ~147 bp DNA 1.7 圈；H1 結合 linker DNA。組蛋白富含 Lys、Arg（正電）。"
              "原核：類核（nucleoid）以 HU、IHF 等蛋白壓縮，無組蛋白（古菌有）。"},
    {"id": "oric_initiation", "title": "Initiation at E. coli oriC", "layout": "scene", "scene": ORIC,
     "tags": ["replication"],
     "notes": "oriC 245 bp：3 個 13-mer（DUE，AT-rich，先解開）＋ 多個 9-mer DnaA box。訊號：完全甲基化的 GATC（11 個，Dam methylase）。"
              "DnaA–ATP 多聚化並纏繞 DNA（IHF、HU、FIS 協助彎曲，gyrase 提供負超螺旋）→ DUE 打開 → DnaC（helicase loader）把 DnaB 六聚體裝到兩條單股上 → SSB → DnaG primase 合成引子 → Pol III 全酶 → replisome。"
              "防止重複起始：SeqA 結合半甲基化 GATC；DnaA–ATP 被水解成 DnaA–ADP。primosome = DnaB + DnaG。"},
    {"id": "bacterial_chromosome_replication", "title": "Replicating a circular bacterial chromosome", "layout": "scene",
     "scene": CHROMOSOME_CYCLE, "tags": ["replication", "dna_structure"],
     "notes": "Gyrase（topo II，耗 ATP）引入負超螺旋（ΔLk = −2），讓 oriC 容易解開。雙向複製形成 θ 結構，兩個複製叉在終點區（Ter 序列＋Tus 蛋白阻擋）會合。"
              "複製完成的兩個環互相套住（catenane），由 topo IV 切開分離。E. coli 複製叉 ~1000 nt/s，真核 ~50 nt/s。"},
    {"id": "topoisomerase_types", "title": "Topoisomerase I and II", "layout": "scene", "scene": TOPO,
     "tags": ["replication", "dna_structure"],
     "notes": "Topo I：切斷一股（活性 Tyr 與磷酸形成共價中間物，能量保留所以不耗 ATP），另一股旋轉後接回 → 連結數 Lk 改變 1。真核 topo I 可鬆弛正、負超螺旋。"
              "Topo II：二聚體切斷雙股（G 段），讓另一段雙股（T 段）穿過再接回 → ΔLk = 2，耗 ATP。細菌 gyrase 是唯一能主動引入負超螺旋的 topo II；topo IV 分離連環體。"
              "藥物：camptothecin、irinotecan、topotecan（topo I）；etoposide、doxorubicin（topo II，抗癌）；quinolone（ciprofloxacin）抑制 gyrase 與 topo IV。"},
    {"id": "telomere_replication", "title": "Telomeres and telomerase", "layout": "scene", "scene": TELOMERE,
     "tags": ["replication"],
     "notes": "末端複製問題：延遲股最後的 RNA 引子移除後無法補齊 → 每次分裂端粒變短（並留下 3′ overhang）。"
              "Telomerase：TERT（反轉錄酶）＋ 自帶 RNA 模板（人 3′-AAUCCC-5′），把 3′ 端延長 TTAGGG 重複；之後 Pol α–primase 補上 C-rich 股。"
              "端粒 10–20 kb，3′ overhang > 100 nt，折成 T-loop（overhang 侵入雙股形成 D-loop）或 G-quadruplex（K⁺ 穩定），由 shelterin（TRF1、TRF2、POT1）保護。"
              "生殖細胞、幹細胞與多數癌細胞有 telomerase 活性；體細胞沒有 → 老化。"},
    {"id": "eukaryotic_origin_licensing", "title": "Eukaryotic origin licensing and firing", "layout": "scene", "scene": LICENSING,
     "tags": ["replication"],
     "notes": "G1（CDK 低）：ORC 結合起點 → Cdc6、Cdt1 把 MCM2–7 以雙六聚體套在雙股 DNA 上（licensed，但不活化）。"
              "S 期：CDK 與 DDK（Cdc7–Dbf4）磷酸化 → Cdc45 + GINS 加入成 CMG helicase，兩組朝相反方向移動 → 雙向複製。"
              "S/G2 時 CDK 高、geminin 抑制 Cdt1 → 不能再裝 MCM → 每個起點每週期只發一次（防止重複複製）。"
              "真核聚合酶：α（含 primase，起始 RNA–DNA 引子）、ε（領先股）、δ（延遲股）、γ（粒線體）、β（BER 修補）；PCNA = clamp，RFC = clamp loader；引子移除：Pol δ 股置換把引子推成 flap → FEN1 切除（長 flap 加 Dna2）；RNase H1 可降解 RNA 但切不掉最後一個與 DNA 相連的 rN。"},
    {"id": "dna_polymerase_mechanism", "title": "DNA polymerase: adding a nucleotide", "layout": "scene", "scene": POL_MECH,
     "tags": ["replication"],
     "notes": "DNA polymerase 外形如右手（palm、fingers、thumb）。Palm 上的 Asp 抓住兩個 Mg²⁺：一個活化引子的 3′-OH，另一個穩定 dNTP 的三磷酸。"
              "3′-OH 攻擊進來 dNTP 的 α-磷酸 → 形成 3′,5′-磷酸二酯鍵，釋出 PPi；pyrophosphatase 水解 PPi 使反應不可逆。"
              "需要引子（3′-OH），新股只能 5′→3′。校對：3′→5′ exonuclease；DNA pol α、β、反轉錄酶、Taq 沒有。"
              "E. coli Pol I 另有 5′→3′ exonuclease（去除引子、nick translation）；Klenow fragment 去掉了這段。"},
    {"id": "dna_polymerases", "title": "DNA polymerases: E. coli vs eukaryotes", "layout": "scene", "scene": POL_TABLE,
"tags": ["replication"],
     "notes": "所有 DNA polymerase：5′→3′ 合成、需要模板與引子、兩個 Mg²⁺。3′→5′ exo（校對）例外：反轉錄酶、Taq、Pol α、Pol β、Pol IV/V、RNA 依賴的 RNA pol。"
              "5′→3′ exo：E. coli Pol I（移除引子、修復），Taq 也有。真核引子主要由 Pol δ 股置換 + FEN1（長 flap 加 Dna2）移除，RNase H1 為輔。"},
    {"id": "replication_fork", "title": "The bacterial replication fork", "layout": "scene", "scene": REPLICATION_FORK,
     "tags": ["replication"],
     "notes": "E. coli 複製叉：helicase（DnaB，六聚體，ATP）在 5′→3′ 方向沿延遲股模板解開雙股；SSB 保護單股；gyrase（topo II）在前方解除正超螺旋。"
              "Primase（DnaG）合成 RNA 引子 → Pol III（β clamp 增加持續性）延伸。領先股連續合成；延遲股以 Okazaki 片段（原核 1000–2000 nt，真核 100–200 nt）不連續合成。"
              "Pol I 以 5′→3′ exonuclease 移除引子並補上 DNA（nick translation），DNA ligase（E. coli 用 NAD⁺，真核用 ATP）封閉切口。"
              "新股一律 5′→3′ 延伸（3′-OH 攻擊 dNTP 的 α-P）；模板 3′→5′ 讀取。"},
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


def main(data):
    data = Path(data)
    for m in MOLECULES:
        d = {"id": m["id"], "name": m["name"]}
        if m.get("label"):
            d["label"] = m["label"]
        d["aliases"] = m.get("aliases", [])
        d["smiles"] = m["smiles"]
        for k in ("haworth", "ring"):
            if m.get(k):
                d[k] = m[k]
        d["tags"] = with_parent(m["tags"])
        d["notes_zh"] = m["notes"]
        d["parts_png"] = False
        d["status"] = "pending"
        dump(data / "molecules" / f"{m['id']}.yaml", d)
    for p in PATHWAYS:
        d = {k: v for k, v in p.items() if k not in ("tags", "notes")}
        d["tags"] = with_parent(p["tags"])
        d["notes_zh"] = p["notes"]
        d = tidy(d)
        d["parts_png"] = False
        d["status"] = "pending"
        dump(data / "pathways" / f"{p['id']}.yaml", d)
    print(f"{len(PATHWAYS)} pathways")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data")
