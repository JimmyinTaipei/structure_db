"""Round 20: molecular biology III — transcription and RNA processing (illustrated scenes, an RNA-polymerase
table, the 5′ cap, modified nucleosides and the two RNA-editing deaminations).

    python tools/gen_mb_tx.py data
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from gen_mb_repair import with_notes, GAP, W, down, dump, duplex, note, seg, table, txt, with_parent  # noqa: E402
from figure_text import tidy  # noqa: E402
from reaction_colors import COLORS, hl  # noqa: E402

M, C = "magenta", "cyan"
YEL, TEAL = "#E8B400", "#3FA7A0"          # promoter / signal elements drawn as coloured DNA or RNA stretches


def rna(pts, ends=None, kind="rna", w=W):
    d = {"g": "strand", "kind": kind, "points": [list(p) for p in pts], "width": w}
    if ends:
        d["ends"] = list(ends)
    return d


def box(x0, x1, y, hue, text="", h=0.85, ts=0.62):
    """Exon / RNA segment drawn as a flat box on the RNA line (text inside)."""
    return {"g": "blob", "at": [(x0 + x1) / 2, y], "w": x1 - x0, "h": h, "shape": "box", "hue": hue,
            "text": text, "text_size": ts}


def ball(x, y, text, hue="yellow", r=0.38, ts=0.5):
    return {"g": "ion", "at": [x, y], "r": r, "hue": hue, "text": text, "text_size": ts}


def prot(x, y, w, h, hue, text, shape="round", ts=0.62, **kw):
    return {"g": "blob", "at": [x, y], "w": w, "h": h, "shape": shape, "hue": hue, "text": text, "text_size": ts, **kw}


def both(x0, x1, y, kind):
    """Colour a stretch of both DNA strands (promoter element)."""
    return [seg(x0, x1, y, kind), seg(x0, x1, y + GAP, kind)]


def bubble(y, x0, x1, a, b, h=0.5, ends=True):
    """Duplex with a melted bubble between a and b (top strand bows up, template bows down)."""
    t = [(x0, y), (a, y), (a + 0.6, y - h), (b - 0.6, y - h), (b, y), (x1, y)]
    u = [(x0, y + GAP), (a, y + GAP), (a + 0.6, y + GAP + h), (b - 0.6, y + GAP + h), (b, y + GAP), (x1, y + GAP)]
    return [rna(t, ("5'", "3'") if ends else None, "parent"), rna(u, ("3'", "5'") if ends else None, "parent")]


def hairpin(x, y, up=1.7, r=0.32, gap=0.3):
    """Points of an RNA stem-loop rising from (x, y); the stem is `gap` wide."""
    a, b = x - gap / 2, x + gap / 2
    top = y - up
    h = math.sqrt(max(r * r - (gap / 2) ** 2, 0))
    cy = top - h
    loop = [(x + r * math.cos(math.radians(t)), cy + r * math.sin(math.radians(t)))
            for t in range(int(math.degrees(math.atan2(h, -gap / 2))), int(math.degrees(math.atan2(h, gap / 2))) - 360, -20)]
    return [(a, y), (a, top)] + loop + [(b, top), (b, y)]


def rungs(x, y, up, gap=0.3, n=4):
    return [rna([(x - gap / 2, y - up * (k + 1) / (n + 1)), (x + gap / 2, y - up * (k + 1) / (n + 1))], kind="gray", w=0.08)
            for k in range(n)]


def lariat(cx, cy, r=1.1, tail=None, a0=25):
    """Lariat loop (centre, radius) closed at the branch A (angle 0); optional tail from A to `tail`."""
    pts = [(cx + r * math.cos(math.radians(t)), cy + r * math.sin(math.radians(t))) for t in range(a0, 361, 15)]
    if tail:
        pts += tail
    return [rna(pts), ball(cx + r, cy, "A", "red", 0.3, 0.42)]


def arc(x0, x1, y, h=1.1, above=True):
    """Thin splice line joining two exon ends."""
    s = -1 if above else 1
    return rna([(x0, y + s * 0.45), ((x0 + x1) / 2, y + s * h), (x1, y + s * 0.45)], kind="gray", w=0.08)


# ================================================================ prokaryotic vs eukaryotic mRNA (lecture 4-114)
MRNA = {"items": [
    txt(0.0, -0.9, "Prokaryote (E. coli): polycistronic mRNA, translated while it is still being made", 0.9, None, "left"),
    *duplex(0.6, 0, 26, top="parent", bot="parent", ends=True),
    *both(0.5, 2.3, 0.6, "#8A8A8A"),
    {"g": "region", "x0": 3.0, "y0": 0.4, "x1": 9.0, "y1": 1.8, "hue": "blue", "opacity": 0.85, "outline": False}, txt(6.0, 1.1, "gene A", 0.7),
    {"g": "region", "x0": 9.6, "y0": 0.4, "x1": 15.6, "y1": 1.8, "hue": "green", "opacity": 0.85, "outline": False}, txt(12.6, 1.1, "gene B", 0.7),
    {"g": "region", "x0": 16.2, "y0": 0.4, "x1": 22.2, "y1": 1.8, "hue": "orange", "opacity": 0.85, "outline": False}, txt(19.2, 1.1, "gene C", 0.7),
    txt(1.4, 2.5, "promoter", 0.65, "gray"), txt(23.6, 2.5, "terminator", 0.65, "gray"),
    *both(22.8, 24.4, 0.6, "#8A8A8A"),
    down(27.4, 1.0, 4.0), txt(27.9, 2.6, "transcription", 0.7, "gray", "left"),
    rna([(2.2, 5.0), (24.0, 5.0)], ("5'", "3'")),
    box(3.6, 8.4, 5.0, "blue", "ORF A"), box(10.2, 15.0, 5.0, "green", "ORF B"), box(16.8, 21.6, 5.0, "orange", "ORF C"),
    prot(4.6, 4.1, 1.5, 0.9, "purple", "", "oval"), prot(11.3, 4.1, 1.5, 0.9, "purple", "", "oval"),
    txt(13.0, 6.5, "each ORF has its own Shine–Dalgarno + AUG; ribosomes (purple) load before transcription ends", 0.72, "gray"),
    txt(0.0, 9.0, "Eukaryote: monocistronic; processed in the nucleus, translated in the cytoplasm", 0.9, None, "left"),
    *duplex(10.4, 0, 26, top="parent", bot="parent", ends=True),
    *both(0.5, 2.3, 10.4, "#8A8A8A"),
    {"g": "region", "x0": 3.0, "y0": 10.2, "x1": 6.0, "y1": 11.6, "hue": "blue", "opacity": 0.85, "outline": False}, txt(4.5, 10.9, "exon 1", 0.7),
    {"g": "region", "x0": 6.0, "y0": 10.2, "x1": 9.0, "y1": 11.6, "hue": "gray", "opacity": 0.3, "outline": False}, txt(7.5, 10.9, "intron", 0.65, "gray"),
    {"g": "region", "x0": 9.0, "y0": 10.2, "x1": 12.5, "y1": 11.6, "hue": "green", "opacity": 0.85, "outline": False}, txt(10.75, 10.9, "exon 2", 0.7),
    {"g": "region", "x0": 12.5, "y0": 10.2, "x1": 16.0, "y1": 11.6, "hue": "gray", "opacity": 0.3, "outline": False}, txt(14.25, 10.9, "intron", 0.65, "gray"),
    {"g": "region", "x0": 16.0, "y0": 10.2, "x1": 20.0, "y1": 11.6, "hue": "orange", "opacity": 0.85, "outline": False}, txt(18.0, 10.9, "exon 3", 0.7),
    down(27.4, 10.8, 13.8), txt(27.9, 11.6, "transcription +\ncapping, splicing,\npoly(A)", 0.7, "gray", "left"),
    ball(2.6, 15.0, "m^7G", "green", 0.55, 0.42),
    rna([(3.15, 15.0), (19.8, 15.0)], ("", "3'")),
    box(3.4, 5.0, 15.0, "gray", "5'UTR", ts=0.5), box(5.0, 8.0, 15.0, "blue", "1"), box(8.0, 11.5, 15.0, "green", "2"),
    box(11.5, 15.5, 15.0, "orange", "3"), box(15.5, 17.4, 15.0, "gray", "3'UTR", ts=0.5),
    txt(5.0, 13.9, "AUG", 0.62), txt(15.5, 13.9, "stop", 0.62),
    txt(17.8, 15.85, "AAAA…A_n (~250)", 0.7, "#C2338A", "left"),
    note(17.6, "Prokaryotic mRNA: no cap, no introns, no poly(A) tail; 5′-ppp end; half-life of minutes. Both have 5′ and 3′ UTRs.\n"
         "Eukaryotic mRNA: 5′ m⁷G cap, introns removed, 3′ poly(A); exported through nuclear pores before translation.", 0.78),
]}

# ================================================================ E. coli initiation: the σ cycle (notes p.49–50)
SY = [2.5, 8.0, 13.5, 19.0]


def holo(y, x=9.5, sigma=True, bub=None, ends=True):
    out = [prot(x, y + 0.5, 7.0, 3.6, "blue", "RNA pol core", "oval", 0.6, text_dy=1.3)]
    if sigma:
        out.append(prot(x - 3.4, y - 1.05, 3.2, 1.2, "orange", "σ^{70}", "oval", 0.65))
    out += bubble(y, 0, 20, *bub, ends=ends) if bub else duplex(y, 0, 20, top="parent", bot="parent", ends=ends)
    return out


SIGMA = {"items": [
    *holo(SY[0]), *both(3.6, 4.8, SY[0], YEL), *both(7.4, 8.6, SY[0], YEL),
    txt(4.2, SY[0] + 1.75, "−35", 0.7),
    {**seg(7.9, 8.1, SY[0] + GAP, YEL), "label": "−10", "label_at": [7.4, SY[0] + 3.1], "point": [8.0, SY[0] + GAP]},
    {**seg(10.1, 10.3, SY[0] + GAP, "parent"), "label": "+1", "label_at": [11.0, SY[0] + 3.1], "point": [10.2, SY[0] + GAP]},
    txt(21.2, SY[0] + 0.5, "1. Closed complex: holoenzyme (α_2ββ'ωσ)\n    σ reads −35 (TTGACA) and −10 (TATAAT)", 0.8, None, "left"),
    *holo(SY[1], bub=(7.4, 10.4)), *both(3.6, 4.8, SY[1], YEL),
    txt(21.2, SY[1] + 0.5, "2. Open complex: DNA melted from −10 to +2\n    (~13 bp bubble; no helicase needed)", 0.8, None, "left"),
    *holo(SY[2], bub=(7.4, 10.4)), *both(3.6, 4.8, SY[2], YEL),
    rna([(9.6, SY[2] + 1.1), (10.6, SY[2] + 1.1)]),
    rna([(13.4, SY[2] - 1.7), (14.4, SY[2] - 1.5)]), rna([(13.9, SY[2] - 0.9), (14.9, SY[2] - 0.8)]),
    txt(21.2, SY[2] + 0.5, "3. Abortive initiation: short RNAs (2–9 nt)\n    are made and released again and again", 0.8, None, "left"),
    prot(13.5, SY[3] + 0.5, 7.0, 3.6, "blue", "RNA pol core", "oval", 0.6, text_dy=1.3),
    prot(2.6, SY[3] - 2.0, 3.4, 1.2, "orange", "σ^{70}", "oval", 0.65),
    *bubble(SY[3], 0, 20, 11.8, 15.0), *both(3.6, 4.8, SY[3], YEL),
    rna([(10.4, SY[3] - 1.7), (11.4, SY[3] - 0.9), (12.4, SY[3] + 0.95), (14.6, SY[3] + 1.05)], ("5'", "")),
    txt(21.2, SY[3] + 0.5, "4. Promoter escape (RNA ≥ ~10 nt): σ leaves;\n    core (+ NusA) elongates; σ is reused", 0.8, None, "left"),
    {"g": "arrow", "from": [-1.9, SY[3] - 1.4], "to": [-1.9, SY[0] - 0.9], "color": "#E8731A", "width": 1.2},
    txt(-2.4, (SY[0] + SY[3]) / 2, "σ\ncycle", 0.75, "#E8731A", "right"),
    note(SY[3] + 3.6, "σ factors choose the promoter: σ^{70} housekeeping (spacer 16–19 bp), σ^{32} heat shock, σ^{54} nitrogen genes (−24/−12; needs an ATP-using enhancer-binding activator).\n"
         "σ region 4.2 binds −35, region 2.4 binds −10, region 2.3 melts DNA. −10 and −35 lie on the same face of B-DNA (~10.5 bp/turn).\n"
         "Rifampicin binds β and blocks the exit channel (abortive stage); actinomycin D intercalates in DNA. Footprinting maps the promoter (RNA pol covers 60–80 bp, including −35 and −10).", 0.76),
]}

# ================================================================ elongation: bubble, hybrid, twin supercoils (notes p.49, 4-122)
EY = 4.0
ELONG = {"items": [
    prot(12.2, EY + 0.5, 9.6, 5.2, "blue", "RNA polymerase", "oval", 0.7, text_dy=2.0, accent=True),
    *duplex(EY, 0, 8.2, top="parent", bot="parent", ends="left"),
    *duplex(EY, 16.2, 26, top="parent", bot="parent", ends="right"),
    rna([(8.2, EY), (9.0, EY), (9.9, EY - 1.0), (14.5, EY - 1.0), (15.4, EY), (16.2, EY)], kind="parent"),
    rna([(8.2, EY + GAP), (9.0, EY + GAP), (9.9, EY + GAP + 1.0), (14.5, EY + GAP + 1.0), (15.4, EY + GAP), (16.2, EY + GAP)], kind="parent"),
    rna([(10.0, 0.8), (10.5, 2.0), (11.1, EY + 0.9), (11.9, EY + 1.5), (14.6, EY + 1.5)], ("5'", "3'")),
    prot(3.2, EY + 0.5, 1.6, 1.5, "green", "topo I", "round", 0.5),
    prot(21.6, EY + 0.5, 2.4, 1.6, "purple", "gyrase", "bilobe", 0.5),
    txt(3.6, 0.6, "behind: underwound\n(−) supercoils", 0.75), txt(21.6, 0.6, "ahead: overwound\n(+) supercoils", 0.75),
    txt(-0.4, EY - 0.75, "coding strand", 0.62, "gray", "left"), txt(-0.4, EY + GAP + 0.75, "template strand", 0.62, "gray", "left"),
    {"id": "hybrid", "g": "strand", "kind": "rna", "points": [[12.0, EY + 1.5], [14.6, EY + 1.5]], "width": W,
     "label": "RNA–DNA hybrid ~8 bp", "label_at": [18.0, EY + 3.3], "label_align": "left", "point": [13.6, EY + 1.6], "label_size": 0.8},
    txt(12.2, EY + 4.0, "transcription bubble ~17 bp", 0.72),
    {"g": "arrow", "from": [0.5, EY + 4.6], "to": [6.5, EY + 4.6], "style": "fat"},
    txt(3.5, EY + 6.2, "direction of transcription", 0.72, "gray"),
    note(EY + 8.9, "RNA is made 5′→3′, antiparallel to the template; its sequence equals the coding (sense) strand with U for T. No primer is needed.\n"
         "E. coli ~50–90 nt/s. No 3′→5′ proofreading exonuclease (error ~10^{-4}–10^{-5}); a misincorporation makes the enzyme backtrack and\n"
         "GreA/B (E. coli) or TFIIS (eukaryotes) cleave the 3′ end. Twin-domain model: topo I removes (−) behind, gyrase removes (+) ahead.\n"
         "Elongation factors: NusA (E. coli); TFIIS, P-TEFb (CTD Ser2 kinase), elongin (eukaryotes).", 0.76),
]}

# ================================================================ termination (notes p.53, lecture 4-135)
TY = [4.5, 13.0]
HP1 = hairpin(5.3, 2.5, 2.0, 0.45, 0.42)
HP2 = hairpin(4.6, TY[1] - 1.9, 2.0, 0.45, 0.42)
TERM = {"items": [
    txt(0.0, -1.8, "ρ-independent (intrinsic) terminator", 0.95, None, "left"),
    prot(9.6, TY[0] + 0.5, 5.6, 3.4, "blue", "RNA pol", "oval", 0.6, text_dy=-1.15),
    *duplex(TY[0], 0, 14, top="parent", bot="parent", ends=True),
    rna([(0.4, 2.5), (3.4, 2.5), (4.9, 2.5)] + HP1 + [(5.7, 2.5), (6.9, 2.9), (8.0, TY[0] + 0.3), (8.9, TY[0] + 0.55), (10.4, TY[0] + 0.55)], ("5'", "")),
    *rungs(5.3, 2.5, 2.0, 0.42),
    txt(4.6, 0.6, "GC-rich stem-loop", 0.72, None, "right"),
    {"id": "urun", "g": "strand", "kind": "rna", "points": [[8.9, TY[0] + 0.55], [10.4, TY[0] + 0.55]], "width": W,
     "label": "U-run: weak rU·dA hybrid", "label_at": [11.0, 1.2], "point": [9.8, TY[0] + 0.4], "label_size": 0.75},
    txt(10.6, TY[0] + 3.4, "pause at the hairpin + weak\nhybrid → RNA released", 0.72, "gray"),
    down(-1.0, TY[0] + 2.4, TY[1] - 1.8),
    *duplex(TY[1], 0, 14, top="parent", bot="parent", ends=True),
    rna([(0.4, TY[1] - 1.9), (4.2, TY[1] - 1.9)] + HP2 + [(5.0, TY[1] - 1.9), (7.8, TY[1] - 1.9)], ("5'", "3'")),
    txt(6.8, TY[1] - 1.15, "…UUUU-OH", 0.65, "#C2338A"),
    *rungs(4.6, TY[1] - 1.9, 2.0, 0.42),
    prot(12.2, TY[1] - 1.9, 3.2, 1.9, "blue", "RNA pol", "oval", 0.55),
    txt(0.0, TY[1] + 2.4, "inverted repeat (palindrome) in DNA", 0.7, "gray", "left"),

    txt(17.0, -1.8, "ρ-dependent terminator", 0.95, None, "left"),
    prot(27.4, TY[0] + 0.5, 5.6, 3.4, "blue", "RNA pol", "oval", 0.6, text_dy=-1.15),
    *duplex(TY[0], 17, 31, top="parent", bot="parent", ends=True),
    rna([(17.6, 1.6), (24.2, 1.6), (25.4, 2.4), (26.4, TY[0] + 0.55), (28.2, TY[0] + 0.55)], ("5'", "")),
    rna([(18.4, 1.6), (20.8, 1.6)], kind=TEAL, w=0.55),
    {"g": "ring", "kind": "hexamer", "at": [22.4, 1.6], "r": 0.62, "hue": "purple", "text": "ρ"},
    txt(19.6, -0.25, "rut site (C-rich, G-poor)", 0.72),
    {"g": "arrow", "from": [23.6, 0.2], "to": [26.0, 0.2], "color": "#E8731A", "width": 1.1},
    txt(24.0, TY[0] + 3.4, "ρ (hexameric RNA helicase) moves\n5′→3′ along the RNA (ATP)", 0.72, "gray"),
    down(16.0, TY[0] + 2.4, TY[1] - 1.8),
    *duplex(TY[1], 17, 31, top="parent", bot="parent", ends=True),
    rna([(17.6, TY[1] - 1.9), (25.0, TY[1] - 1.9)], ("5'", "3'")),
    rna([(18.4, TY[1] - 1.9), (20.8, TY[1] - 1.9)], kind=TEAL, w=0.55),
    {"g": "ring", "kind": "hexamer", "at": [24.2, TY[1] - 1.9], "r": 0.62, "hue": "purple", "text": "ρ"},
    prot(28.6, TY[1] - 1.9, 3.2, 1.9, "blue", "RNA pol", "oval", 0.55),
    txt(17.0, TY[1] + 2.4, "ρ unwinds the RNA–DNA hybrid", 0.7, "gray", "left"),
    note(TY[1] + 4.6, "Intrinsic: GC-rich stem (stable) followed by 6–8 U; the polymerase stalls and the A=U hybrid melts. ρ-dependent: no hairpin needed; inhibited by bicyclomycin;\n"
         "NusA helps intrinsic termination, NusG helps ρ. Eukaryotes: Pol II — cleavage at the poly(A) site, then the 5′→3′ exonuclease Xrn2 ('torpedo') catches the polymerase;\n"
         "Pol III — a run of T in the template (no hairpin); Pol I — terminator protein (TTF-I).", 0.76),
]}

# ================================================================ Pol II preinitiation complex (notes p.49, 51; lecture 4-123, 4-134)
PY = 11.0
PIC = {"items": [
    txt(0.0, -1.2, "Core promoter (class II gene)", 0.9, None, "left"),
    *duplex(1.0, 0, 26, top="parent", bot="parent", ends=True),
    *both(3.2, 4.2, 1.0, "#9C86C4"), *both(4.6, 6.0, 1.0, YEL), *both(10.6, 11.8, 1.0, "#4E9A51"), *both(17.4, 18.6, 1.0, TEAL),
    txt(3.7, 0.1, "BRE", 0.65), txt(5.3, 3.0, "TATA box (−30)\nTATAAA", 0.68), txt(11.2, 3.0, "Inr (+1)", 0.68), txt(18.0, 3.0, "DPE (+30)", 0.68),
    {"g": "arrow", "from": [11.2, 0.6], "to": [13.2, -0.3], "curve": -0.3, "color": "black", "width": 1.1},
    txt(0.0, 4.9, "Assembly: TFIID (TBP binds TATA and bends DNA) → TFIIA → TFIIB → Pol II·TFIIF → TFIIE → TFIIH = preinitiation complex (PIC)", 0.78, None, "left"),
    prot(12.2, PY + 0.6, 6.0, 5.4, "blue", "Pol II", "oval", 0.85, text_dy=-1.6, accent=True),
    *duplex(PY, 0, 26, top="parent", bot="parent", ends=True),
    prot(2.3, PY - 0.75, 1.6, 1.4, "pink", "IIA"),
    prot(4.9, PY - 0.55, 2.4, 1.1, "dkgreen", "TBP", "kidney", 0.58, rot=180),
    prot(4.9, PY - 2.5, 3.4, 2.6, "green", "TFIID\n(TAFs)", "oval", 0.6),
    prot(7.7, PY - 0.85, 1.6, 1.6, "orange", "IIB"),
    prot(15.2, PY - 1.3, 1.7, 1.6, "teal", "IIF"),
    prot(16.9, PY - 2.4, 1.7, 1.6, "purple", "IIE"),
    prot(19.3, PY - 2.7, 3.0, 2.2, "yellow", "TFIIH", "bilobe", 0.6),
    rna([(14.6, PY + 2.6), (16.2, PY + 3.6), (18.4, PY + 3.9), (20.6, PY + 4.4), (23.6, PY + 4.4)], kind="gray", w=0.24),
    ball(18.4, PY + 3.9, "P", "yellow", 0.32, 0.45), ball(20.6, PY + 4.4, "P", "yellow", 0.32, 0.45), ball(22.6, PY + 4.4, "P", "yellow", 0.32, 0.45),
    {"g": "arrow", "from": [20.6, PY - 1.6], "to": [22.4, PY + 3.7], "curve": -0.3, "color": "#E8731A", "width": 1.1},
    txt(21.4, PY - 2.4, "CDK7 (kinase)", 0.68, "#E8731A", "left"),
    txt(24.2, PY + 4.4, "CTD: 52 × YSPTSPS\n(Rpb1 tail)", 0.68, None, "left"),
    note(PY + 8.0, "TFIID = TBP + TAFs; reads TATA, Inr, DPE (TBP is also used by Pol I and Pol III). TFIIA stabilises TBP–DNA. TFIIB bridges TBP and Pol II and picks +1 (σ-like).\n"
         "TFIIF escorts Pol II. TFIIE recruits TFIIH. TFIIH: helicase (XPB) opens the promoter (ATP); CDK7 phosphorylates CTD Ser5 → promoter clearance; TFIIH also works in NER.\n"
         "Ser5-P CTD recruits the capping enzymes; P-TEFb (CDK9) phosphorylates Ser2 → productive elongation and loading of splicing and poly(A) factors.\n"
         "Unphosphorylated Pol II (IIA) joins the PIC; phosphorylated Pol II (IIO) elongates. α-Amanitin blocks Pol II translocation.", 0.76),
]}

# ================================================================ RNA polymerases (lecture 4-121, 4-131)
POL_TABLE = table([
    ("E. coli RNA pol", "one enzyme for all RNAs · core α_2ββ'ω + σ = holoenzyme · rifampicin binds β"),
    ("Pol I (nucleolus)", "45S pre-rRNA → 28S, 18S, 5.8S rRNA · α-amanitin resistant · ~50 % of transcription"),
    ("Pol II (nucleoplasm)", "pre-mRNA, most snRNAs, miRNA, lncRNA · α-amanitin very sensitive · CTD heptad repeats"),
    ("Pol III (nucleoplasm)", "tRNA, 5S rRNA, U6 snRNA, 7SL RNA · α-amanitin moderately sensitive · internal promoter (A, B boxes)"),
    ("Mitochondrial", "single subunit, T7-phage-like · transcribes mtDNA (both strands, polycistronic)"),
    ("Shared features", "no primer, no proofreading exonuclease; 5′→3′; Mg^{2+}; actinomycin D (intercalator) blocks all"),
], w_name=8.6, x_desc=9.2, hue="teal")

# ================================================================ 5′ capping (notes p.57; lecture 4-145)
CY = [1.0, 4.6, 8.2, 11.8, 15.4]


def capped(y, n_p=3, g=None, m7=False, m2=False):
    out = [ball(6.0, y, "N", "blue", 0.4, 0.5), rna([(6.45, y), (13.5, y)], ("", "3'"))]
    xs = [5.1, 4.15, 3.2][:n_p]
    out += [ball(x, y, "P", "yellow") for x in xs]
    if g:
        out += [ball(3.2, y, "P", "teal"), ball(2.1, y, "m^7G" if m7 else "G", "green", 0.55, 0.42 if m7 else 0.5)]
    if m7:
        out.append({"g": "mark", "kind": "methyl", "at": [2.1, y - 0.55], "side": "up", "color": M})
    if m2:
        out.append({"g": "mark", "kind": "methyl", "at": [6.0, y - 0.4], "side": "up", "color": M})
    return out


def enz(i, text):
    y = (CY[i] + CY[i + 1]) / 2
    return [down(8.4, CY[i] + 0.7, CY[i + 1] - 0.7), txt(9.0, y, text, 0.78, None, "left")]


CAP = {"items": [
    *capped(CY[0]), txt(3.2, CY[0] - 0.85, "γ", 0.6, "gray"), txt(4.15, CY[0] - 0.85, "β", 0.6, "gray"), txt(5.1, CY[0] - 0.85, "α", 0.6, "gray"),
    *enz(0, "RNA triphosphatase removes the γ-phosphate (→ P_i)"),
    *capped(CY[1], 2),
    *enz(1, "guanylyltransferase: GTP → GMP joined 5′→5′ (→ PP_i)"),
    *capped(CY[2], 2, g=True), txt(4.15, CY[2] + 0.95, "5′–5′ triphosphate", 0.62, "gray"),
    *enz(2, "guanine-N7-methyltransferase (SAM → SAH)"),
    *capped(CY[3], 2, g=True, m7=True), txt(14.6, CY[3], "cap 0", 0.8, None, "left"),
    *enz(3, "2′-O-methyltransferase (SAM) on nucleotide 1 (→ cap 1; nt 2 too → cap 2)"),
    *capped(CY[4], 2, g=True, m7=True, m2=True), txt(14.6, CY[4], "cap 1", 0.8, None, "left"),
    txt(15.2, CY[0] + 0.05, "pppN… (transcript 5′ end)", 0.72, "gray", "left"),
    note(CY[4] + 2.7, "Added co-transcriptionally (RNA ~25 nt) by enzymes bound to the Ser5-phosphorylated CTD of Pol II; only Pol II transcripts are capped.\n"
         "Functions: blocks 5′→3′ exonucleases, nuclear export (cap-binding complex CBC), translation initiation (eIF4E binds the cap), splicing of the first intron.\n"
         "Decapping (Dcp1/Dcp2) → Xrn1 degrades the mRNA 5′→3′. teal = from GTP.", 0.76),
]}

# ================================================================ cleavage and polyadenylation (notes p.59; lecture 4-145, 4-160)
AY = [3.0, 8.6, 13.8]
POLYA = {"items": [
    rna([(0.0, AY[0]), (24.0, AY[0])], ("5'", "3'")),
    rna([(7.0, AY[0]), (9.2, AY[0])], kind=YEL, w=0.55), rna([(15.2, AY[0]), (18.2, AY[0])], kind=TEAL, w=0.55),
    prot(8.1, AY[0] - 1.45, 3.0, 1.3, "orange", "CPSF", "oval"),
    prot(16.7, AY[0] - 1.45, 2.8, 1.3, "teal", "CstF", "oval"),
    prot(11.4, AY[0] + 1.3, 1.6, 1.2, "purple", "PAP"), prot(13.4, AY[0] + 1.3, 2.0, 1.1, "gray", "CF I/II", "oval", 0.55),
    {"g": "mark", "kind": "cut", "at": [12.2, AY[0]], "side": "up"},
    txt(12.2, AY[0] - 1.6, "cleavage\n(CA)", 0.68), txt(8.1, AY[0] + 0.95, "AAUAAA", 0.72, "#B08F10"),
    txt(16.7, AY[0] + 0.95, "GU-rich", 0.72, "#2F7F7A"),
    txt(2.0, AY[0] + 1.2, "pre-mRNA 3′ end", 0.7, "gray"),
    down(4.0, AY[0] + 1.8, AY[1] - 0.8),
    rna([(0.0, AY[1]), (12.2, AY[1])], ("5'", "3'")), rna([(7.0, AY[1]), (9.2, AY[1])], kind=YEL, w=0.55),
    rna([(13.4, AY[1]), (24.0, AY[1])], kind="gray"), rna([(15.2, AY[1]), (18.2, AY[1])], kind=TEAL, w=0.55),
    prot(14.0, AY[1], 1.4, 1.1, "pink", "Xrn2", "round", 0.5),
    txt(19.0, AY[1] + 1.2, "downstream RNA degraded by Xrn2 → Pol II terminates", 0.68, "gray"),
    down(4.0, AY[1] + 1.1, AY[2] - 0.8),
    rna([(0.0, AY[2]), (12.2, AY[2])], ("5'", "")), rna([(7.0, AY[2]), (9.2, AY[2])], kind=YEL, w=0.55),
    rna([(12.2, AY[2]), (20.6, AY[2])], ("", "3'"), kind=YEL),
    prot(13.7, AY[2], 1.3, 0.95, "teal", "", "oval"), prot(15.5, AY[2], 1.3, 0.95, "teal", "", "oval"),
    prot(17.3, AY[2], 1.3, 0.95, "teal", "", "oval"),
    prot(21.6, AY[2] - 1.3, 1.6, 1.2, "purple", "PAP"),
    txt(15.2, AY[2] - 1.3, "poly(A) tail ~200–250 A + PABP", 0.72, "#B08F10"),
    txt(10.0, AY[2] + 1.4, "poly(A) polymerase: ATP → AMP + PP_i, no template", 0.72, None, "left"),
    note(AY[2] + 3.4, "CPSF binds AAUAAA; CstF binds the GU/U-rich element ~30 nt downstream; CPSF-73 cuts 10–30 nt after AAUAAA (usually after CA).\n"
         "Poly(A): stability, export, translation (PABP–eIF4G closes the loop); shortening (deadenylation) starts decay. Histone mRNAs have a stem-loop instead.\n"
         "Alternative poly(A) sites change the 3′ UTR (miRNA sites) or the C terminus (IgM membrane vs secreted). Bacterial poly(A) instead marks RNA for decay.", 0.76),
]}

# ================================================================ spliceosome (notes p.58; lecture 4-151)
PX = 21.6
SPY = [2.2, 8.6, 14.6, 19.8]
ARCH = [(4.0, SPY[1]), (4.9, SPY[1] - 0.6), (6.8, SPY[1] - 2.3), (10.0, SPY[1] - 2.8), (13.2, SPY[1] - 2.3), (15.1, SPY[1] - 0.6), (16.0, SPY[1])]
SPLICE = {"items": [
    rna([(0, SPY[0]), (20, SPY[0])], ("5'", "3'")),
    box(0, 4.0, SPY[0], "blue", "exon 1"), box(16.0, 20.0, SPY[0], "green", "exon 2"),
    ball(11.6, SPY[0], "A", "red", 0.32, 0.42),
    txt(4.6, SPY[0] + 0.85, "GU", 0.72), txt(13.6, SPY[0] + 0.85, "(Py)_n", 0.68), txt(15.3, SPY[0] + 0.85, "AG", 0.72),
    prot(4.9, SPY[0] - 1.2, 1.5, 1.2, "teal", "U1"), prot(11.6, SPY[0] - 1.2, 1.5, 1.2, "purple", "U2"),
    prot(14.6, SPY[0] - 1.15, 2.1, 1.0, "gray", "U2AF", "oval", 0.55),
    txt(PX, SPY[0], "1. U1 snRNP pairs with the 5′ splice site (GU);\n    U2 pairs at the branch point (A bulges out);\n    U2AF binds the Py tract and the 3′ AG", 0.78, None, "left"),
    {"g": "region", "x0": 3.3, "y0": SPY[1] - 3.6, "x1": 16.7, "y1": SPY[1] + 0.9, "hue": "gray", "opacity": 0.18, "dashed": True},
    rna([(0, SPY[1]), (4.0, SPY[1])]), rna(ARCH), rna([(16.0, SPY[1]), (20, SPY[1])]),
    box(0, 4.0, SPY[1], "blue", "exon 1"), box(16.0, 20.0, SPY[1], "green", "exon 2"),
    ball(12.4, SPY[1] - 2.55, "A", "red", 0.32, 0.42),
    {"g": "arrow", "from": [11.9, SPY[1] - 2.2], "to": [4.6, SPY[1] - 0.5], "curve": 0.25, "color": "#E8731A", "width": 1.2},
    ball(6.2, SPY[1] - 2.9, "1", "orange", 0.3, 0.45),
    txt(10.0, SPY[1] - 4.1, "spliceosome: U2, U5, U6 snRNPs + ~150 proteins", 0.68, "gray"),
    txt(PX, SPY[1] - 0.6, "2. U4/U6·U5 join; U1 and U4 leave; U6 pairs with\n    the 5′ splice site and with U2 → RNA catalytic core.\n"
        "    (1) the branch-A 2′-OH attacks the 5′ splice site", 0.78, None, "left"),
    rna([(0, SPY[2]), (4.0, SPY[2])]), box(0, 4.0, SPY[2], "blue", "exon 1"), txt(4.9, SPY[2] - 0.75, "3'-OH", 0.68, "#C2338A"),
    *lariat(9.4, SPY[2] - 1.5, 1.2, tail=[(11.4, SPY[2] - 0.8), (12.8, SPY[2]), (16.0, SPY[2]), (20, SPY[2])]),
    box(16.0, 20.0, SPY[2], "green", "exon 2"),
    {"g": "arrow", "from": [4.6, SPY[2] + 0.6], "to": [15.8, SPY[2] + 0.6], "curve": 0.12, "color": "#E8731A", "width": 1.2},
    ball(10.2, SPY[2] + 2.3, "2", "orange", 0.3, 0.45),
    txt(PX, SPY[2] - 0.4, "3. Lariat intermediate (2′–5′ link at the branch A);\n    (2) the exon-1 3′-OH attacks the 3′ splice site", 0.78, None, "left"),
    rna([(0, SPY[3]), (8.0, SPY[3])], ("5'", "3'")), box(0, 4.0, SPY[3], "blue", "exon 1"), box(4.0, 8.0, SPY[3], "green", "exon 2"),
    *lariat(13.0, SPY[3] - 0.4, 1.1, tail=[(15.4, SPY[3] - 0.4)]),
    txt(PX, SPY[3] - 0.4, "4. Exons joined (mRNA); lariat released,\n    debranched and degraded", 0.78, None, "left"),
    note(SPY[3] + 2.9, "Two transesterifications: the number of phosphodiester bonds does not change (no ATP for the chemistry; ATP drives the RNA helicases that remodel snRNPs).\n"
         "Major (U2-type) GU–AG introns > 99 %; minor (U12-type) AU–AC introns use U11, U12, U4atac, U6atac + U5. snRNAs are U-rich; U3 is a snoRNA (rRNA processing), not a splicing snRNA.\n"
         "Splice-site mutations: β-thalassaemia; SMN2 exon-7 skipping (spinal muscular atrophy). Anti-Sm (snRNP) antibodies in SLE.", 0.76),
]}

# ================================================================ group I / group II / spliceosome (notes p.60)
MX = [0.0, 11.0, 22.0]
MY = [2.8, 8.6, 14.2]


def pre(x, y, a=False):
    arch = [(x + 2.4, y), (x + 2.9, y - 0.6), (x + 3.6, y - 2.0), (x + 4.5, y - 2.4), (x + 5.4, y - 2.0), (x + 6.1, y - 0.6), (x + 6.6, y)]
    out = [rna([(x, y), (x + 2.4, y)]), rna(arch), rna([(x + 6.6, y), (x + 9.0, y)]),
           box(x, x + 2.4, y, "blue", "1", 0.75), box(x + 6.6, x + 9.0, y, "green", "2", 0.75)]
    if a:
        out.append(ball(x + 5.45, y - 1.95, "A", "red", 0.3, 0.42))
    return out


def step1(x, y, group1=False):
    out = [rna([(x, y), (x + 2.4, y)]), box(x, x + 2.4, y, "blue", "1", 0.75), txt(x + 2.3, y - 0.8, "3'-OH", 0.6, "#C2338A")]
    if group1:
        out += [rna([(x + 4.0, y - 0.6), (x + 4.1, y - 1.9), (x + 4.9, y - 2.4), (x + 5.7, y - 2.0), (x + 6.2, y - 0.6), (x + 6.6, y), (x + 9.0, y)]),
                ball(x + 3.85, y - 0.25, "G", "green", 0.36, 0.48)]
    else:
        out += lariat(x + 4.5, y - 1.5, 0.9, tail=[(x + 5.9, y - 0.8), (x + 6.6, y), (x + 9.0, y)])
    out += [box(x + 6.6, x + 9.0, y, "green", "2", 0.75),
            {"g": "arrow", "from": [x + 2.5, y + 0.55], "to": [x + 6.7, y + 0.55], "curve": 0.3, "color": "#E8731A", "width": 1.1}]
    return out


def products(x, y, group1=False):
    out = [rna([(x, y), (x + 4.8, y)]), box(x, x + 2.4, y, "blue", "1", 0.75), box(x + 2.4, x + 4.8, y, "green", "2", 0.75)]
    if group1:
        out += [ball(x + 6.0, y - 1.3, "G", "green", 0.36, 0.48), rna([(x + 6.4, y - 1.3), (x + 9.0, y - 1.3)], ("", "3'"))]
    else:
        out += lariat(x + 7.2, y - 1.0, 0.8, tail=[(x + 9.0, y - 1.0)])
    return out


MECH = {"items": [
    txt(4.5, -1.6, "Group I intron", 0.95), txt(15.5, -1.6, "Group II intron", 0.95), txt(26.5, -1.6, "Spliceosomal intron", 0.95),
    *pre(MX[0], MY[0]), ball(1.3, MY[0] - 2.2, "G", "green", 0.36, 0.48), txt(1.3, MY[0] - 3.1, "G-3'OH", 0.6),
    {"g": "arrow", "from": [1.7, MY[0] - 1.9], "to": [2.6, MY[0] - 0.45], "curve": 0.2, "color": "#E8731A", "width": 1.1},
    *pre(MX[1], MY[0], True),
    {"g": "arrow", "from": [MX[1] + 5.1, MY[0] - 1.7], "to": [MX[1] + 2.6, MY[0] - 0.45], "curve": -0.2, "color": "#E8731A", "width": 1.1},
    *pre(MX[2], MY[0], True),
    {"g": "arrow", "from": [MX[2] + 5.1, MY[0] - 1.7], "to": [MX[2] + 2.6, MY[0] - 0.45], "curve": -0.2, "color": "#E8731A", "width": 1.1},
    txt(MX[2] + 4.5, MY[0] + 1.25, "+ U1, U2, U4, U5, U6 snRNPs", 0.62, "gray"),
    txt(MX[0] + 4.5, MY[0] + 1.25, "free guanosine (GTP) attacks", 0.62, "gray"), txt(MX[1] + 4.5, MY[0] + 1.25, "branch-point A 2′-OH attacks", 0.62, "gray"),
    *[down(x + 4.5, MY[0] + 1.8, MY[1] - 3.0) for x in MX],
    *step1(MX[0], MY[1], True), *step1(MX[1], MY[1]), *step1(MX[2], MY[1]),
    *[down(x + 4.5, MY[1] + 2.2, MY[2] - 2.0) for x in MX],
    *products(MX[0], MY[2], True), *products(MX[1], MY[2]), *products(MX[2], MY[2]),
    txt(MX[0], MY[2] + 1.6, "linear intron (5′ G added)\nself-splicing ribozyme:\nTetrahymena pre-rRNA, organelles", 0.68, None, "left"),
    txt(MX[1], MY[2] + 1.6, "lariat intron\nself-splicing ribozyme:\nfungal / plant organelles", 0.68, None, "left"),
    txt(MX[2], MY[2] + 1.6, "lariat intron\nsnRNPs; catalytic core = U2/U6 RNA:\nnuclear pre-mRNA", 0.68, None, "left"),
    note(MY[2] + 4.6, "All three use two transesterifications. Group II introns and the spliceosome share the lariat chemistry (group II is the likely ancestor; also LINE-like retrohoming).\n"
         "Nuclear pre-tRNA introns are removed by proteins instead (endonuclease + ligase, uses ATP). Cech (self-splicing) and Altman (RNase P) — ribozymes, Nobel 1989.", 0.76),
]}

# ================================================================ alternative splicing (notes p.59–60; lecture 4-145)


def pre_mrna(x, y, exons, intron=1.2):
    """Pre-mRNA: exons [(text, width, hue)] separated by introns; returns items and exon (x0, x1) list."""
    out, spans, cx = [], [], x
    for i, (t, w, h) in enumerate(exons):
        spans.append((cx, cx + w))
        cx += w + (intron if i < len(exons) - 1 else 0)
    out.append(rna([(x - 0.3, y), (cx + 0.3, y)]))
    out += [box(a, b, y, h, t, 0.75) for (a, b), (t, _, h) in zip(spans, exons)]
    return out, spans


def mature(x, y, exons, tail=True):
    out, cx = [ball(x, y, "", "green", 0.3)], x + 0.35
    end = cx + sum(w for _, w, _ in exons)
    out.append(rna([(cx, y), (end + 0.2, y)]))
    for t, w, h in exons:
        out.append(box(cx, cx + w, y, h, t, 0.75))
        cx += w
    if tail:
        out.append(txt(end + 0.4, y, "A_n", 0.7, "#C2338A", "left"))
    return out


E1, E2, E3 = ("1", 1.8, "blue"), ("2", 1.8, "orange"), ("3", 1.8, "green")
ALT_A, SA = pre_mrna(0.4, 1.8, [E1, E2, E3])
ALT_B, SB = pre_mrna(16.4, 1.8, [E1, E2, E3])
CT = [("1", 1.3, "blue"), ("2", 1.3, "blue"), ("3", 1.3, "blue"), ("4", 1.6, "orange"), ("5", 1.6, "purple"), ("6", 1.6, "purple")]
ALT_C, SC = pre_mrna(0.4, 11.2, CT, 0.9)
ALT = {"items": [
    txt(0.0, -0.6, "Constitutive splicing", 0.9, None, "left"),
    *ALT_A, arc(SA[0][1], SA[1][0], 1.8), arc(SA[1][1], SA[2][0], 1.8),
    down(4.6, 2.6, 4.0), *mature(0.6, 4.9, [E1, E2, E3]),
    txt(16.0, -0.6, "Alternative splicing (exon skipping)", 0.9, None, "left"),
    *ALT_B, arc(SB[0][1], SB[1][0], 1.8), arc(SB[1][1], SB[2][0], 1.8), arc(SB[0][1], SB[2][0], 1.8, 1.25, above=False),
    *mature(16.6, 4.9, [E1, E2, E3]), txt(24.4, 4.9, "isoform A", 0.7, None, "left"),
    *mature(16.6, 6.4, [E1, E3]), txt(24.4, 6.4, "isoform B (exon 2 skipped)", 0.7, None, "left"),
    txt(0.0, 9.0, "Calcitonin / CGRP gene:\nalternative poly(A) site + splicing", 0.9, None, "left"),
    *ALT_C, txt(SC[3][1] + 0.45, 12.1, "pA1", 0.6, "#C2338A"), txt(SC[5][1] + 0.45, 12.1, "pA2", 0.6, "#C2338A"),
    *mature(0.6, 14.2, CT[:4]), txt(7.6, 14.2, "thyroid C cells → calcitonin", 0.72, None, "left"),
    *mature(0.6, 15.8, CT[:3] + CT[4:]), txt(9.2, 15.8, "neurons → CGRP (vasodilator)", 0.72, None, "left"),
    txt(19.0, 9.0, "trans-splicing (trypanosomes, nematodes)", 0.9, None, "left"),
    rna([(19.2, 10.8), (22.6, 10.8)]), box(19.2, 20.8, 10.8, "pink", "SL", 0.75), txt(23.0, 10.8, "spliced-leader RNA", 0.65, "gray", "left"),
    rna([(20.8, 12.4), (25.0, 12.4)]), box(25.0, 27.4, 12.4, "green", "exon", 0.75), txt(28.0, 12.4, "pre-mRNA", 0.65, "gray", "left"),
    down(23.2, 13.2, 14.7), *mature(19.6, 15.6, [("SL", 1.6, "pink"), ("exon", 2.4, "green")]),
    txt(24.6, 16.7, "+ Y-shaped intron (not a lariat)", 0.62, "gray"),
    note(18.4, "~95 % of human multi-exon genes are alternatively spliced: exon skipping, mutually exclusive exons, alternative 5′ / 3′ sites, intron retention.\n"
         "Chosen by SR proteins (bind exonic splicing enhancers) vs hnRNPs (silencers). Drosophila Dscam: > 38 000 possible isoforms. Sex lethal / transformer: sex determination.", 0.76),
]}

# ================================================================ RNA editing (notes p.60; lecture 4-157)
ED = {"items": [
    txt(0.0, -0.8, "C → U: apolipoprotein B mRNA (APOBEC-1 cytidine deaminase)", 0.9, None, "left"),
    rna([(3.0, 1.2), (17.0, 1.2)], ("5'", "3'")), box(9.0, 11.0, 1.2, "yellow", "CAA", 0.85, 0.65), txt(10.0, 2.25, "codon 2153 (Gln)", 0.65, "gray"),
    {"g": "arrow", "from": [8.6, 2.6], "to": [5.0, 4.2], "color": "gray", "width": 1.0},
    {"g": "arrow", "from": [11.4, 2.6], "to": [15.0, 4.2], "color": "gray", "width": 1.0},
    txt(3.6, 2.7, "liver: unedited", 0.7, None, "center"), txt(16.6, 2.7, "intestine: edited", 0.7, None, "center"),
    rna([(0.0, 5.0), (9.0, 5.0)]), box(3.5, 5.5, 5.0, "yellow", "CAA", 0.85, 0.65),
    rna([(11.0, 5.0), (20.0, 5.0)]), box(14.5, 16.5, 5.0, "red", "UAA", 0.85, 0.65),
    {"g": "region", "x0": 0.0, "y0": 6.4, "x1": 9.0, "y1": 7.0, "hue": "purple", "opacity": 0.9, "outline": False},
    {"g": "region", "x0": 11.0, "y0": 6.4, "x1": 15.4, "y1": 7.0, "hue": "purple", "opacity": 0.9, "outline": False},
    txt(4.5, 7.8, "apoB-100 (4536 aa): VLDL, LDL", 0.72), txt(15.5, 7.8, "apoB-48 (2152 aa, ~48 %): chylomicron", 0.72),
    txt(0.0, 10.0, "A → I: ADAR (adenosine deaminase acting on RNA) in double-stranded RNA; I is read as G", 0.9, None, "left"),
    rna([(0.6, 11.6), (7.4, 11.6)] + [(7.4 + 0.5 * math.cos(math.radians(t)), 12.1 + 0.5 * math.sin(math.radians(t))) for t in range(-90, 91, 20)]
        + [(7.4, 12.6), (0.6, 12.6)], ("5'", "3'")),
    *[rna([(x, 11.75), (x, 12.45)], kind="gray", w=0.08) for x in (1.2, 1.9, 2.6, 3.3, 4.6, 5.3, 6.0, 6.7)],
    ball(4.0, 11.6, "A", "red", 0.33, 0.45),
    {"g": "arrow", "from": [9.0, 12.1], "to": [11.4, 12.1], "style": "fat"}, txt(10.2, 10.95, "ADAR", 0.7),
    rna([(12.6, 11.6), (19.4, 11.6)] + [(19.4 + 0.5 * math.cos(math.radians(t)), 12.1 + 0.5 * math.sin(math.radians(t))) for t in range(-90, 91, 20)]
        + [(19.4, 12.6), (12.6, 12.6)], ("5'", "3'")),
    *[rna([(x, 11.75), (x, 12.45)], kind="gray", w=0.08) for x in (13.2, 13.9, 14.6, 15.3, 16.6, 17.3, 18.0, 18.7)],
    ball(16.0, 11.6, "I", "teal", 0.33, 0.45),
    txt(0.0, 14.0, "GluR-B (AMPA receptor) Q/R site: CAG (Gln) → CIG = CGG (Arg) → Ca^{2+}-impermeable channel; also 5-HT_{2C} receptor; tRNA wobble I (ADAT)", 0.72, None, "left"),
    txt(0.0, 16.0, "U insertion / deletion: trypanosome mitochondria (kinetoplastid); a guide RNA (gRNA) is the template", 0.9, None, "left"),
    rna([(0.0, 17.6), (19.0, 17.6)], ("5'", "3'")),
    *[txt(x, 16.95, "U", 0.7, "red") for x in (6.0, 7.2, 9.4, 11.8)],
    rna([(4.0, 18.6), (15.0, 18.6)], ("3'", "5'"), kind="gray"),
    txt(16.6, 18.6, "gRNA", 0.7, "gray", "left"),
    note(20.2, "Editing changes the RNA sequence after transcription (the gene still reads CAA). APOBEC family also includes APOBEC3G (antiviral; HIV Vif degrades it)\n"
         "and AID (somatic hypermutation, class switching). The first editing found was U insertion in trypanosomes. Compare CRISPR–Cas9: a guide RNA edits DNA.", 0.76),
]}

# ================================================================ tRNA processing (notes p.54; lecture 4-138)
TR = {"items": [
    {"g": "trna", "at": [5.0, 6.2], "leader": 2.4, "trailer": 2.4, "cca": False, "intron": True},
    {"g": "mark", "kind": "cut", "at": [4.78, 3.7], "side": "up"}, {"g": "mark", "kind": "cut", "at": [5.22, 3.7], "side": "up"},
    txt(0.6, 0.2, "RNase P (ribozyme)\ncuts off the 5′ leader", 0.72),
    txt(7.0, 0.2, "3′ trailer: RNase Z (endo, eukaryotes)\nor RNase D / PH / T (exo, E. coli)", 0.72, None, "left"),
    txt(5.0, 11.4, "intron in the anticodon loop:\nsplicing endonuclease + ligase", 0.72),
    txt(2.6, 2.3, "5′ leader", 0.62, "gray", "right"), txt(7.4, 2.5, "3′ trailer (…UUU, La)", 0.62, "gray", "left"),
    {"g": "arrow", "from": [10.4, 6.6], "to": [13.4, 6.6], "style": "fat"},
    {"g": "trna", "at": [18.0, 6.2], "anticodon": "anticodon"},
    txt(15.6, 4.65, "D loop", 0.65, "gray"), txt(20.5, 4.65, "TΨC loop", 0.65, "gray"), txt(17.4, 3.2, "5′", 0.62, "#C2338A"),
    txt(20.0, 2.6, "3′ CCA added by tRNA nucleotidyl-\ntransferase (no template; encoded in E. coli)", 0.72, None, "left"),
    txt(21.5, 8.2, "modified bases: D (dihydro-U),\nΨ, T (ribothymidine), I, m^1G …", 0.72, None, "left"),
    note(14.4, "Pol III transcribes tRNA genes (internal promoter: A and B boxes). Order: RNase P (5′) and RNase Z / exonucleases (3′) → CCA → splicing → modifications → export.\n"
         "E. coli: tRNAs also come out of rRNA operons; CCA is encoded (nucleotidyltransferase only repairs it). RNase P = protein + catalytic M1 RNA.\n"
         "The amino acid is attached to the 3′-terminal A (2′- or 3′-OH) by aminoacyl-tRNA synthetase (round 21).", 0.76),
]}

# ================================================================ rRNA processing (notes p.57; lecture 4-142)


def bar(y, parts):
    """Pre-rRNA: parts [(x0, x1, hue | None, text)]; None = spacer (gray line only)."""
    x0, x1 = parts[0][0], parts[-1][1]
    out = [rna([(x0, y), (x1, y)], kind="gray")]
    out += [box(a, b, y, h, t, 0.8, 0.62) for a, b, h, t in parts if h]
    return out


RY = [1.5, 4.6, 10.6, 13.8]
RRNA = {"items": [
    txt(0.0, -0.6, "E. coli: rrn operon (7 copies) → one 30S pre-rRNA, one RNA polymerase", 0.9, None, "left"),
    *bar(RY[0], [(0.0, 0.8, None, ""), (0.8, 6.8, "blue", "16S"), (6.8, 7.3, None, ""), (7.3, 9.1, "yellow", "tRNA"),
                 (9.1, 9.6, None, ""), (9.6, 18.6, "green", "23S"), (18.6, 19.4, None, ""), (19.4, 20.8, "orange", "5S"), (20.8, 21.6, None, "")]),
    *[{"g": "mark", "kind": "cut", "at": [x, RY[0] - 0.4], "side": "up"} for x in (0.8, 6.8, 9.6, 18.6)],
    txt(22.4, RY[0], "RNase III, E, P\n+ trimming", 0.68, "gray", "left"),
    down(10.5, RY[0] + 0.7, RY[1] - 0.7),
    box(0.8, 6.8, RY[1], "blue", "16S"), txt(3.8, RY[1] + 1.0, "→ 30S subunit (+21 proteins)", 0.68),
    box(9.6, 18.6, RY[1], "green", "23S"), box(19.0, 20.4, RY[1], "orange", "5S"), txt(15.0, RY[1] + 1.0, "→ 50S subunit (+ ~33 proteins)", 0.68),
    txt(22.4, RY[1], "30S + 50S = 70S", 0.75, None, "left"),
    txt(0.0, 7.9, "Eukaryotes: Pol I makes 45S pre-rRNA in the nucleolus; the 5S rRNA gene is transcribed separately by Pol III", 0.9, None, "left"),
    *bar(RY[2], [(0.0, 1.6, None, ""), (1.6, 7.0, "blue", "18S"), (7.0, 8.2, None, ""), (8.2, 9.4, "yellow", "5.8S"),
                 (9.4, 10.6, None, ""), (10.6, 20.0, "green", "28S"), (20.0, 21.4, None, "")]),
    *[{"g": "mark", "kind": "cut", "at": [x, RY[2] - 0.4], "side": "up"} for x in (1.6, 7.0, 8.2, 9.4, 10.6, 20.0)],
    txt(0.8, RY[2] + 0.85, "5′ETS", 0.58, "gray"), txt(7.6, RY[2] + 0.85, "ITS1", 0.58, "gray"),
    txt(10.0, RY[2] + 0.85, "ITS2", 0.58, "gray"), txt(20.7, RY[2] + 0.85, "3′ETS", 0.58, "gray"),
    down(10.5, RY[2] + 1.3, RY[3] - 0.7),
    box(1.6, 7.0, RY[3], "blue", "18S"), txt(4.3, RY[3] + 1.0, "→ 40S (+ ~33 proteins)", 0.68),
    box(8.2, 9.4, RY[3], "yellow", "5.8S"), box(10.6, 20.0, RY[3], "green", "28S"), box(20.4, 21.8, RY[3], "orange", "5S"),
    txt(15.0, RY[3] + 1.0, "→ 60S (+ ~47 proteins; 5S from Pol III)", 0.68),
    txt(22.4, RY[3], "40S + 60S = 80S", 0.75, None, "left"),
    note(RY[3] + 3.0, "Modification guided by snoRNPs (snoRNA + enzyme): C/D-box snoRNA → 2′-O-methylation (fibrillarin); H/ACA-box snoRNA → pseudouridylation (dyskerin;\n"
         "mutated in dyskeratosis congenita). E. coli modifies rRNA with enzymes alone (no snoRNA). U3 snoRNA guides 5′ETS cleavage; RNase MRP cuts ITS1.\n"
         "Ribosomes are assembled in the nucleolus and exported as subunits. The rRNA (23S / 28S) catalyses peptide bonds: a ribozyme.", 0.76),
]}


# ================================================================ molecules and reactions
MOLECULES = [
    {"id": "mrna_cap", "name": "5′ Cap (m7GpppNm, cap 1)", "label": "5' Cap (m^7GpppAm)", "aliases": ["m7G cap", "5′ cap", "帽子結構"],
     "smiles": "C[n+]1cn([C@@H]2O[C@H](COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OC)[C@@H]3O[19*])[C@@H](O)[C@H]2O)c2nc(N)[nH]c(=O)c21",
     "ring": {"coordgen": True}, "tags": ["rna_processing"],
     "notes": "7-甲基鳥苷（m⁷G，N7 帶正電）以 5′–5′ 三磷酸鍵接在 mRNA 第一個核苷酸（此處為 A）上。Cap 0 只有 m⁷G；cap 1 第一個核苷酸 2′-O-甲基化（如圖）；cap 2 第二個也甲基化。"
              "甲基都來自 SAM。功能：抗 5′ 外切酶、出核（CBC）、轉譯起始（eIF4E）、第一個 intron 的剪接。Rib = 其餘 RNA。"},
    {"id": "uridine", "name": "Uridine", "label": "Uridine", "aliases": ["尿苷"], "tags": ["rna_processing"],
     "smiles": "O=c1ccn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)[nH]1",
     "haworth": {"ring": "furanose", "subs": {1: {"up": "base", "down": "H"}, 2: {"up": "H", "down": "OH"}, 3: {"up": "H", "down": "OH"},
                                               4: {"up": "HOCH_2", "down": "H"}},
                 "base": "[*]n1ccc(=O)[nH]c1=O", "base_right": 7, "numbers": True},
     "notes": "Uracil 以 N1 接核糖 C1′（N-醣苷鍵）。"},
    {"id": "pseudouridine", "name": "Pseudouridine", "label": "Pseudouridine (Ψ)", "aliases": ["Ψ", "假尿苷"], "tags": ["rna_processing"],
     "smiles": "O=c1[nH]cc([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)[nH]1",
     "haworth": {"ring": "furanose", "subs": {1: {"up": "base", "down": "H"}, 2: {"up": "H", "down": "OH"}, 3: {"up": "H", "down": "OH"},
                                               4: {"up": "HOCH_2", "down": "H"}},
                 "base": "[*]c1c[nH]c(=O)[nH]c1=O", "base_right": 6, "numbers": True},
     "notes": "Uridine 的異構物：鹼基轉 180°，改以 C5 接 C1′（C-醣苷鍵），多出一個 N1–H 氫鍵供體。H/ACA snoRNA + dyskerin（真核）或 pseudouridine synthase 催化，"
              "不需要能量、分子量不變。最常見的 RNA 修飾：tRNA 的 TΨC 環、rRNA、snRNA。"},
    {"id": "cytosine_rib", "name": "Cytosine (in RNA)", "label": "Cytosine (in RNA)", "aliases": [], "tags": ["rna_processing"],
     "smiles": "[19*]n1ccc(N)nc1=O", "ring": {"coordgen": True}, "notes": "RNA 中的 cytosine（Rib = 核糖與骨架）。"},
    {"id": "uracil_rib", "name": "Uracil (in RNA)", "label": "Uracil (in RNA)", "aliases": [], "tags": ["rna_processing"],
     "smiles": "[19*]n1ccc(=O)[nH]c1=O", "ring": {"coordgen": True}, "notes": "RNA 中的 uracil（Rib = 核糖與骨架）。"},
    {"id": "adenine_rib", "name": "Adenine (in RNA)", "label": "Adenine (in RNA)", "aliases": [], "tags": ["rna_processing"],
     "smiles": "[19*]n1cnc2c(N)ncnc21", "ring": {"coordgen": True}, "notes": "RNA 中的 adenine（Rib = 核糖與骨架）。"},
    {"id": "hypoxanthine_rib", "name": "Hypoxanthine (inosine in RNA)", "label": "Hypoxanthine (I, in RNA)", "aliases": ["inosine in RNA"],
     "tags": ["rna_processing"], "smiles": "[19*]n1cnc2c(=O)[nH]cnc21", "ring": {"coordgen": True},
     "notes": "RNA 中的 inosine 鹼基。I 與 C 配對（像 G），所以 A→I 編輯在轉譯時被讀成 G；tRNA 反密碼子第一位的 I 可與 U、C、A 擺動配對。"},
]

REACTIONS = [
    {"id": "apobec1_editing", "title": "C → U RNA editing (APOBEC-1)", "substrates": ["cytosine_rib"], "products": ["uracil_rib"],
     "cofactors_in": ["H2O"], "cofactors_out": ["ammonium"], "enzyme": "APOBEC-1\n(cytidine deaminase)", "tags": ["rna_processing"],
     "card": {"type": "水解脫胺：C4 的 NH₂（洋紅）離開成 NH₄⁺，H₂O 的 O（青）成為 C4=O",
              "exam_points": ["apoB mRNA codon 2153 CAA → UAA（終止）：小腸 apoB-48，肝臟不編輯 → apoB-100",
                              "需輔助蛋白 ACF；同家族 APOBEC3G 為抗病毒蛋白、AID 負責體細胞高度突變"]}},
    {"id": "adar_editing", "title": "A → I RNA editing (ADAR)", "substrates": ["adenine_rib"], "products": ["hypoxanthine_rib"],
     "cofactors_in": ["H2O"], "cofactors_out": ["ammonium"], "enzyme": "ADAR\n(adenosine deaminase acting on RNA)", "tags": ["rna_processing"],
     "card": {"type": "水解脫胺：C6 的 NH₂（洋紅）離開，H₂O 的 O（青）成為 C6=O → inosine 鹼基 hypoxanthine",
              "exam_points": ["只作用於雙股 RNA；I 被讀成 G（CAG Gln → CIG 讀作 Arg：GluR-B Q/R 位）",
                              "與 adenosine deaminase（ADA，游離核苷，SCID）同一化學反應、不同受質"]}},
]

COLORS.update({
    "apobec1_editing": {"cytosine_rib": hl("[NH2:1]c"), "uracil_rib": hl("[O:1]=c1ccn(*)c(=O)[nH]1", col=C),
                        "H2O": {"colors": {0: C}}, "ammonium": {"colors": {0: M}}},
    "adar_editing": {"adenine_rib": hl("[NH2:1]c"), "hypoxanthine_rib": hl("[O:1]=c", col=C),
                     "H2O": {"colors": {0: C}}, "ammonium": {"colors": {0: M}}},
})

# ================================================================ pathways

PATHWAYS = [
    {"id": "mrna_prok_vs_euk", "title": "Prokaryotic vs eukaryotic mRNA", "layout": "scene", "scene": MRNA, "tags": ["transcription", "rna_processing"],
     "notes": "原核：polycistronic（一條 mRNA 多個 ORF，各有 SD + AUG）、無核膜 → 轉錄與轉譯偶聯、無 intron、無 5′ cap、無 poly(A)。"
              "真核：monocistronic、核內加工（5′ cap、剪接、3′ poly(A) ~250 A）後出核轉譯。兩者都有 5′UTR、3′UTR。"},
    {"id": "bacterial_transcription_initiation", "title": "E. coli transcription initiation (σ cycle)", "layout": "scene", "scene": SIGMA, "tags": ["transcription"],
     "notes": "Holoenzyme（α₂ββ′ωσ）的 σ 認 −35（TTGACA）與 −10（TATAAT，Pribnow box）→ closed complex → 不需 helicase 自行解開 −10 至 +2（open complex）"
              " → abortive initiation（反覆合成 2–9 nt 短 RNA）→ RNA 超過約 10 nt 時 promoter escape，σ 離開並重複使用（σ cycle），core + NusA 延長。"
              "σ⁷⁰ housekeeping、σ³² heat shock、σ⁵⁴（−24/−12，需 enhancer 活化子）。Rifampicin 結合 β 次單元。"},
    {"id": "transcription_elongation", "title": "Transcription elongation: bubble, hybrid and twin supercoils", "layout": "scene", "scene": ELONG, "tags": ["transcription"],
     "notes": "Transcription bubble ~17 bp，RNA–DNA hybrid ~8 bp。RNA 5′→3′ 合成，與 template 反向平行、序列同 coding strand（U 取代 T）；不需 primer。"
              "無 3′→5′ 外切校對，但可 backtrack 由 GreA/B（原核）或 TFIIS（真核）切除錯誤 3′ 端。"
              "Twin-supercoiling domain：前方 (+) 超螺旋由 gyrase 解、後方 (−) 由 topo I 解。延長因子：NusA；TFIIS、P-TEFb、elongin。"},
    {"id": "transcription_termination", "title": "Transcription termination in E. coli", "layout": "scene", "scene": TERM, "tags": ["transcription"],
     "notes": "ρ-independent（intrinsic）：DNA 的反向重複（palindrome）使 RNA 形成 GC-rich 莖環，接著 6–8 個 U；聚合酶停頓、rU·dA 雜合弱 → 釋放。"
              "ρ-dependent：ρ 六聚體（RNA–DNA helicase，ATP）結合 rut site（C 多 G 少、無二級結構），沿 RNA 5′→3′ 追上聚合酶、拆開雜合。"
              "Bicyclomycin 抑制 ρ。真核 Pol II：poly(A) 位切開後 Xrn2 torpedo；Pol III：模板上一串 T。"},
    {"id": "pol2_preinitiation_complex", "title": "RNA polymerase II promoter and preinitiation complex", "layout": "scene", "scene": PIC, "tags": ["transcription"],
     "notes": "Core promoter：BRE、TATA box（−30，TATAAA，Hogness box）、Inr（+1）、DPE（+30）。"
              "組裝：TFIID（TBP 認 TATA、使 DNA 彎曲）→ IIA → IIB → Pol II·IIF → IIE → IIH。TFIIB 似 σ（決定起點），TFIIF 唯一先與 Pol II 結合，"
              "TFIIH 有 helicase（XPB，也參與 NER）與 CTD kinase（CDK7，Ser5）→ promoter clearance；P-TEFb（CDK9）磷酸化 Ser2 → 延長。"
              "CTD = 52 個 YSPTSPS 重複，作為 cap、剪接、poly(A) 因子的平台。"},
    {"id": "rna_polymerases", "title": "RNA polymerases compared", "layout": "scene", "scene": POL_TABLE, "tags": ["transcription"],
     "notes": "E. coli 一種 RNA pol 合成所有 RNA。真核：Pol I（核仁，45S pre-rRNA → 28S、18S、5.8S）、Pol II（mRNA、多數 snRNA、miRNA）、"
              "Pol III（tRNA、5S rRNA、U6 snRNA、7SL）。α-amanitin（毒蕈）：Pol II 最敏感、Pol III 中度、Pol I 不敏感。粒線體 RNA pol 為單一次單元（似 T7）。"},
    {"id": "mrna_capping", "title": "5′ capping of pre-mRNA", "layout": "scene", "scene": CAP, "tags": ["rna_processing"],
     "notes": "RNA triphosphatase 去 γ-P → guanylyltransferase 由 GTP 加上 GMP（5′–5′ 三磷酸鍵，放出 PPi）→ guanine-N7-methyltransferase（SAM）→ cap 0；"
              "2′-O-methyltransferase 甲基化第 1（cap 1）、第 2（cap 2）個核苷酸的核糖 2′-OH。酵素結合在 Ser5 磷酸化的 CTD，RNA 約 25 nt 時加上。"},
    {"id": "polyadenylation", "title": "Cleavage and polyadenylation", "layout": "scene", "scene": POLYA, "tags": ["rna_processing"],
     "notes": "CPSF 認 AAUAAA，CstF 認下游 GU-rich；在 AAUAAA 後 10–30 nt（常在 CA 之後）切開，poly(A) polymerase（不需模板，ATP）加上 ~200–250 A，PABP 包覆。"
              "下游 RNA 被 Xrn2 分解並使 Pol II 終止（torpedo）。Poly(A)：穩定、出核、轉譯（PABP–eIF4G）；histone mRNA 無 poly(A)（莖環）。"},
    {"id": "spliceosome_splicing", "title": "Pre-mRNA splicing by the spliceosome", "layout": "scene", "scene": SPLICE, "tags": ["rna_processing"],
     "notes": "U1 認 5′ss（GU）、U2 認 branch point A（A 突出）、U2AF 認 Py tract 與 3′ AG → U4/U6·U5 加入 → U1、U4 離開 → U6 與 5′ss、U2 配對成 RNA 催化中心（ribozyme）。"
              "第一次轉酯：branch A 的 2′-OH 攻擊 5′ss → exon 1 3′-OH + lariat（2′–5′ 鍵）；第二次：exon 1 3′-OH 攻擊 3′ss → exon 接合、lariat 釋出。"
              "Major（U2 型，GU–AG，>99 %）；minor（U12 型，AU–AC：U11、U12、U4atac、U6atac）。"},
    {"id": "splicing_mechanisms", "title": "Group I, group II and spliceosomal introns", "layout": "scene", "scene": MECH, "tags": ["rna_processing"],
     "notes": "三種都是兩次轉酯。Group I：游離 G（GTP/guanosine）的 3′-OH 攻擊 5′ss，intron 以線形（5′ 端帶 G）釋出（Tetrahymena rRNA，Cech）。"
              "Group II：intron 內 A 的 2′-OH 攻擊，產生 lariat（胞器）。Spliceosome：同 group II 化學，但需 snRNP。核內 pre-tRNA intron 由酵素（endonuclease + ligase）切除。"},
    {"id": "alternative_splicing", "title": "Alternative splicing and trans-splicing", "layout": "scene", "scene": ALT, "tags": ["rna_processing"],
     "notes": "Constitutive：所有 exon 都保留。Alternative：exon skipping、互斥 exon、替代 5′/3′ 剪接位、intron retention；calcitonin 基因在甲狀腺 C 細胞產生 calcitonin，"
              "在神經產生 CGRP（替代 poly(A) 位 + 剪接）。Trans-splicing：兩條不同 RNA 接合（錐蟲 SL RNA），intron 為 Y 形。SR protein 促進、hnRNP 抑制。"},
    {"id": "rna_editing", "title": "RNA editing", "layout": "scene", "scene": ED, "tags": ["rna_processing"],
     "notes": "C→U：APOBEC-1 使 apoB mRNA codon 2153 CAA → UAA，小腸產生 apoB-48（chylomicron），肝臟不編輯產生 apoB-100（VLDL/LDL）。"
              "A→I：ADAR 作用於雙股 RNA，I 讀成 G（GluR-B Q/R 位，CAG→CIG：Gln→Arg）。U 插入／刪除：錐蟲粒線體以 guide RNA 為模板（最早發現的編輯）。"},
    {"id": "trna_processing", "title": "tRNA processing", "layout": "scene", "scene": TR, "tags": ["rna_processing"],
     "notes": "Pol III 轉錄（內部啟動子 A、B box）。RNase P（核酶）切 5′ leader；3′ trailer：E. coli 由 RNase D 等外切酶、真核由 RNase Z（內切）；"
              "CCA 由 tRNA nucleotidyltransferase 加上（E. coli 基因已編碼）；真核部分 pre-tRNA 在反密碼子環有 intron（endonuclease + ligase）；"
              "鹼基修飾：D、Ψ、T、I、m¹G。"},
    {"id": "rrna_processing", "title": "rRNA processing", "layout": "scene", "scene": RRNA, "tags": ["rna_processing"],
     "notes": "E. coli：rrn operon 產生 30S pre-rRNA（16S、tRNA、23S、5S），由 RNase III 等切開 → 30S、50S → 70S。"
              "真核：Pol I 在核仁產生 45S pre-rRNA（18S、5.8S、28S，夾 ETS 與 ITS），5S 由 Pol III 另外轉錄 → 40S、60S → 80S。"
              "修飾：C/D box snoRNA → 2′-O-甲基化；H/ACA box snoRNA → pseudouridylation（dyskerin，dyskeratosis congenita）。"},
    {"id": "mrna_cap_structure", "title": "Structure of the 5′ cap", "layout": "free", "cell_size": 18.0,
     "nodes": [{"id": "cap", "ref": "mrna_cap", "pos": [0, 0],
                "decorate": {"highlight": [{"smarts": "[CH3][n+]", "color": "magenta"},
                                           {"smarts": "COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])OC", "color": "cyan"},
                                           {"smarts": "[CH3]O[CH1]", "color": "orange"}]}}],
     "edges": [], "text": False, "tags": ["rna_processing"],
     "notes": "洋紅：m⁷G 的 N7-甲基（帶正電）；青：5′–5′ 三磷酸鍵；橘：第一個核苷酸的 2′-O-甲基（cap 1）。"},
    {"id": "rna_modified_nucleosides", "title": "Uridine vs pseudouridine", "layout": "free", "cell_size": 12.0,
     "nodes": [{"id": "u", "ref": "uridine", "pos": [0, 0], "decorate": {"base_colors": {1: "magenta"}}},
               {"id": "psi", "ref": "pseudouridine", "pos": [1.05, 0], "decorate": {"base_colors": {3: "cyan"}}}],
     "edges": [], "text": False, "tags": ["rna_processing"],
     "notes": "Uridine 以 N1（洋紅）接 C1′（N-醣苷鍵）；pseudouridine 鹼基轉 180°，改以 C5 接 C1′（C–C 醣苷鍵），原本的 N1 變成多出的 N–H（青），可多形成一個氫鍵。"},
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
             "cofactors_in": r["cofactors_in"], "cofactors_out": r["cofactors_out"], "cofactors_inline": True, "extra_in": ["H+"],
             "enzyme": r["enzyme"], "reversible": False}
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
