"""Round 21: molecular biology IV — translation and protein fate (ribosome scenes, genetic code, aminoacylation and
peptidyl-transfer chemistry, targeting, ubiquitin–proteasome, factor / inhibitor / modification tables).

    python tools/gen_mb_tl.py data
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from gen_mb_repair import with_notes, down, dump, note, table, txt, with_parent  # noqa: E402
from gen_mb_tx import ball, box, prot, rna  # noqa: E402
from figure_text import tidy  # noqa: E402
from reaction_colors import COLORS, frag, hl  # noqa: E402

M, C = "magenta", "cyan"
OR = "#E8731A"


def ribo(x, y, w=7.0, sites=True, large=None, small=None, **kw):
    d = {"g": "ribosome", "at": [x, y], "w": w, "sites": sites, **kw}
    if large:
        d["large"] = large
    if small:
        d["small"] = small
    return d


def tr(x, y, ac="", aa=None, pep=0, s=0.8, **kw):
    """Compact tRNA icon standing on the mRNA at (x, y)."""
    d = {"g": "trna", "shape": "icon", "at": [x, y], "s": s, "anticodon": ac, **kw}
    if aa:
        d["aa"] = aa
    if pep:
        d["peptide"] = pep
    return d


def mrna(x0, x1, y, ends=("5'", "3'")):
    return rna([(x0, y), (x1, y)], ends)


def pointer(x, y, label, lx, ly, align=None, size=0.75):
    """Invisible anchor carrying a leader-line label (several labels on one big glyph)."""
    d = {"g": "region", "x0": x, "y0": y, "x1": x, "y1": y, "opacity": 0, "outline": False,
         "label": label, "label_at": [lx, ly], "point": [x, y], "label_size": size}
    if align:
        d["label_align"] = align
    return d


def fat(x0, y0, x1, y1):
    return {"g": "arrow", "from": [x0, y0], "to": [x1, y1], "style": "fat"}


def arr(x0, y0, x1, y1, curve=0, color="gray"):
    return {"g": "arrow", "from": [x0, y0], "to": [x1, y1], "curve": curve, "color": color, "width": 1.1}


SLOT = 0.27 * 7.0                    # E / P / A spacing for the default ribosome width

# ================================================================ genetic code (notes p.62–63; lecture 4-159, 4-161)
CODE = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu", "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met", "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser", "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr", "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "stop", "UAG": "stop", "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys", "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "UGU": "Cys", "UGC": "Cys", "UGA": "stop", "UGG": "Trp", "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg", "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}
CLASS = {**{a: "yellow" for a in ("Gly", "Ala", "Val", "Leu", "Ile", "Pro", "Phe", "Trp", "Met")},
         **{a: "green" for a in ("Ser", "Thr", "Cys", "Tyr", "Asn", "Gln")},
         **{a: "pink" for a in ("Asp", "Glu")}, **{a: "blue" for a in ("Lys", "Arg", "His")}, "stop": "gray"}
CW, CH, LH = 6.4, 3.0, 0.68


def code_table():
    items = [txt(2 * CW + 1.2, -2.0, "second base", 0.85, "gray")]
    for j, b in enumerate("UCAG"):
        items.append(txt(1.2 + j * CW + CW / 2, -1.0, b, 1.0))
    items.append(txt(-1.6, 2 * CH, "first\nbase", 0.8, "gray"))
    items.append(txt(4 * CW + 4.0, 2 * CH, "third\nbase", 0.8, "gray"))
    for i, b1 in enumerate("UCAG"):
        y0 = i * CH
        items.append(txt(0.2, y0 + CH / 2 - 0.15, b1, 1.0))
        for j, b2 in enumerate("UCAG"):
            x0 = 1.2 + j * CW
            for k, b3 in enumerate("UCAG"):
                c = b1 + b2 + b3
                aa = CODE[c]
                y = y0 + 0.45 + k * LH
                hue_ = "dkgreen" if c == "AUG" else CLASS[aa]
                items.append({"g": "region", "x0": x0 + 0.1, "y0": y - 0.31, "x1": x0 + CW - 0.1, "y1": y + 0.31,
                              "hue": hue_, "opacity": 0.45 if c != "AUG" else 0.6, "outline": False})
                label = {"stop": "Stop", "Met": "Met (start)"}.get(aa, aa)
                items.append(txt(x0 + 0.4, y, f"{c}   {label}", 0.66, None, "left"))
                if j == 3:
                    items.append(txt(4 * CW + 2.0, y, b3, 0.7, "gray"))
    items.append(txt(0.0, 4 * CH + 0.9, "yellow nonpolar · green polar · pink acidic · blue basic · gray stop · dark green start (AUG)",
                     0.72, "gray", "left"))
    items.append(note(4 * CH + 2.9, "Triplet, non-overlapping, comma-free, read 5′→3′; degenerate (61 sense codons for 20 aa); nearly universal (mitochondria: UGA = Trp, AUA = Met,\n"
                      "AGA/AGG = stop in vertebrates). Single-codon aa: Met (AUG), Trp (UGG); six-codon aa: Leu, Ser, Arg. 21st aa selenocysteine (UGA + SECIS), 22nd pyrrolysine (UAG).\n"
                      "Stop codons are read by release factors, not tRNAs. Codon usage is biased; in E. coli UAA > UGA > UAG.", 0.76))
    return {"items": items}


GENETIC_CODE = code_table()

# ================================================================ wobble (lecture 4-161; notes p.63)
WB = {"items": [
    {"g": "trna", "at": [4.5, 1.6], "s": 0.9},
    *[box(x - 0.42, x + 0.42, 5.4, "orange" if b == "I" else "teal", b, 0.8, 0.62) for x, b in zip((3.6, 4.5, 5.4), "CGI")],
    *[rna([(x, 5.85), (x, 6.75)], kind="gray", w=0.08) for x in (3.6, 4.5, 5.4)],
    mrna(0.6, 8.4, 7.2),
    *[box(x - 0.42, x + 0.42, 7.2, "orange" if b == "U" else "teal", b, 0.8, 0.62) for x, b in zip((3.6, 4.5, 5.4), "GCU")],
    txt(2.9, 5.4, "anticodon 3′", 0.62, "gray", "right"), txt(6.1, 5.4, "5′", 0.62, "gray", "left"),
    txt(6.3, 6.3, "wobble", 0.68, OR, "left"),
    txt(4.5, 8.6, "anticodon IGC reads GCU, GCC and GCA (Ala)", 0.68),
    txt(11.5, -0.6, "Wobble rules (Crick)", 0.9, None, "left"),
    *[x for i, (a, b) in enumerate([("anticodon 5′ base", "pairs with codon 3′ base"), ("C", "G only"), ("A", "U only"),
                                     ("G", "C or U"), ("U", "A or G"), ("I (inosine)", "U, C or A")])
      for x in ([{"g": "blob", "at": [14.0, 0.9 + i * 1.25], "w": 4.6, "h": 1.0, "shape": "box",
                  "hue": "gray" if i == 0 else "teal", "text": a, "text_size": 0.62},
                 txt(17.0, 0.9 + i * 1.25, b, 0.72, None, "left")])],
    note(10.6, "The first two codon positions pair strictly (Watson–Crick); the third is loose, so one tRNA reads several codons. ≥ 32 tRNAs (31 + initiator) read 61 codons.\n"
         "Isoaccepting tRNAs: different tRNAs for the same amino acid (e.g. tRNA^{Ala}_1 anticodon IGC reads GCU/GCC/GCA; tRNA^{Ala}_2 CGC reads GCG).\n"
         "Inosine at the wobble position comes from A → I editing by ADAT. Two codons differing only at the third base (pyrimidine-ending) always code the same amino acid.", 0.76),
]}

# ================================================================ aminoacyl-tRNA synthetase (notes p.64; lecture 4-164)
AY = 3.0
AARS = {"items": [
    prot(3.2, AY, 3.6, 2.6, "teal", "aaRS", "kidney", 0.7),
    ball(2.4, AY - 2.2, "aa", "orange", 0.42, 0.5), txt(4.4, AY - 2.2, "+ ATP", 0.78),
    fat(6.0, AY, 8.4, AY), txt(7.2, AY - 1.2, "PP_i → 2 P_i", 0.68, "gray"),
    prot(11.6, AY, 3.6, 2.6, "teal", "", "kidney", 0.7),
    ball(10.9, AY - 0.3, "aa", "orange", 0.42, 0.5), ball(11.9, AY - 0.3, "AMP", "yellow", 0.42, 0.42),
    txt(11.4, AY + 2.1, "aminoacyl-AMP (enzyme-bound)", 0.68),
    fat(14.4, AY, 16.8, AY), txt(15.6, AY - 1.2, "+ tRNA", 0.75),
    prot(21.0, AY + 2.4, 4.4, 2.4, "teal", "aaRS", "kidney", 0.7, text_dx=1.3, text_dy=0.3),
    {"g": "trna", "at": [20.4, AY + 0.3], "s": 0.7},
    ball(20.55, AY - 2.4, "aa", "orange", 0.4, 0.5), txt(23.2, AY - 1.0, "+ AMP", 0.78, None, "left"),
    txt(23.2, AY - 2.4, "aminoacyl-tRNA (charged)", 0.75, None, "left"),
    note(AY + 6.8, "Step 1 (activation): aa + ATP → aa~AMP + PP_i (PP_i hydrolysed → 2 high-energy bonds per amino acid). Step 2: aa~AMP + tRNA → aa-tRNA + AMP.\n"
         "Class I aaRS (Rossmann fold, monomeric) attach the aa to the 2′-OH of A76, then it moves to 3′-OH; class II (dimeric, antiparallel β) attach directly to 3′-OH.\n"
         "Specificity ('second genetic code'): anticodon + acceptor stem / discriminator base (tRNA^{Ala}: one G3·U70 pair is enough).\n"
         "Proofreading (double sieve): IleRS activates Val by mistake (~1/200) and hydrolyses Val-AMP / Val-tRNA^{Ile} in an editing site → error ~1/50 000.\n"
         "Indirect routes: Glu-tRNA^{Gln} → Gln-tRNA^{Gln} (amidotransferase, many bacteria); Ser-tRNA^{Sec} → Sec-tRNA^{Sec} (selenocysteine).", 0.76),
]}

# ================================================================ ribosome structure (notes p.66; lecture 4-181)
RB = {"items": [
    mrna(0.0, 11.0, 7.0), ribo(5.5, 7.0, 7.0, True, "50S", "30S"),
    tr(5.5 - SLOT, 7.0, "", s=0.8, hue="#D9D9D9"), tr(5.5, 7.0, "UAC", pep=3), tr(5.5 + SLOT, 7.0, "CGA", aa="Ala"),
    txt(13.0, 1.6, "E site: deacylated tRNA leaves", 0.78, None, "left"),
    txt(13.0, 2.8, "P site: peptidyl-tRNA (initiator tRNA starts here)", 0.78, None, "left"),
    txt(13.0, 4.0, "A site: incoming aminoacyl-tRNA", 0.78, None, "left"),
    txt(13.0, 5.6, "large subunit: peptidyl transferase centre (23S / 28S rRNA = ribozyme), exit tunnel", 0.78, None, "left"),
    txt(13.0, 6.8, "small subunit: binds mRNA first, decoding centre (16S / 18S rRNA)", 0.78, None, "left"),
    txt(13.0, 8.0, "mRNA read 5′→3′; protein made N → C terminus", 0.78, None, "left"),
    {"g": "blob", "at": [17.0, 11.4], "w": 8.0, "h": 1.1, "shape": "box", "hue": "purple", "text": "70S (bacteria)", "text_size": 0.7},
    txt(21.5, 11.4, "30S (16S rRNA + 21 proteins) + 50S (23S + 5S rRNA + 33 proteins)", 0.75, None, "left"),
    {"g": "blob", "at": [17.0, 12.9], "w": 8.0, "h": 1.1, "shape": "box", "hue": "pink", "text": "80S (eukaryotes)", "text_size": 0.7},
    txt(21.5, 12.9, "40S (18S rRNA + ~33 proteins) + 60S (28S, 5.8S, 5S rRNA + ~47 proteins)", 0.75, None, "left"),
    txt(17.0, 14.3, "mitochondria: 55S (12S + 16S rRNA)", 0.72, "gray"),
    note(16.2, "How the sites were found: puromycin (looks like Tyr-tRNA) reacts only with peptidyl-tRNA in the P site, not in the A site; three deacylated tRNAs fit on one ribosome\n"
         "but only two aminoacyl-tRNAs → a third (E) site. Ribosome ~2.5 MDa (bacteria), two-thirds rRNA. Polysome: many ribosomes on one mRNA.", 0.76),
]}

# ================================================================ initiation (notes p.65–66, 4-170, 4-179)
IY = [4.0, 15.0]
INIT = {"items": [
    txt(0.0, -1.6, "Bacteria (E. coli)", 0.95, None, "left"),
    mrna(0.0, 11.0, IY[0]), box(1.6, 3.6, IY[0], "yellow", "SD", 0.75), txt(2.6, IY[0] - 0.75, "AGGAGG", 0.62, "#B08F10"),
    prot(5.6, IY[0] + 1.0, 6.2, 1.5, "#C9BEDF", "30S", "oval", 0.55),
    txt(5.6, IY[0] + 2.4, "16S rRNA 3′ end pairs with SD", 0.62, "gray"),
    tr(5.4, IY[0], "UAC", aa="fMet"),
    prot(7.2, IY[0] - 1.6, 1.8, 1.0, "green", "IF2·GTP", "oval", 0.45),
    prot(9.0, IY[0] - 0.55, 1.1, 0.9, "orange", "IF1"), prot(3.7, IY[0] - 1.8, 1.1, 0.9, "yellow", "IF3"),
    fat(12.0, IY[0] - 0.5, 14.4, IY[0] - 0.5), txt(13.4, IY[0] - 1.7, "+ 50S", 0.75),
    txt(13.4, IY[0] + 1.2, "IF2: GTP → GDP\nIF1, 2, 3 leave", 0.65, "gray"),
    mrna(15.0, 26.5, IY[0]), ribo(20.6, IY[0], 7.0, True, "50S", "30S"), tr(20.6, IY[0], "UAC", aa="fMet"),
    txt(20.6, IY[0] + 3.2, "70S initiation complex: fMet-tRNA in P, A site empty", 0.72),
    txt(0.0, IY[1] - 4.6, "Eukaryotes (cap-dependent scanning)", 0.95, None, "left"),
    ball(0.5, IY[1], "m^7G", "green", 0.5, 0.4), mrna(1.0, 27.0, IY[1], ("", "")),
    txt(27.3, IY[1], "AAAA", 0.7, "#C2338A", "left"), prot(26.6, IY[1] - 0.95, 1.6, 0.9, "teal", "PABP", "oval", 0.45),
    prot(1.6, IY[1] - 0.95, 1.2, 0.95, "pink", "4E"), prot(3.2, IY[1] - 1.2, 2.2, 0.9, "gray", "4G", "box", 0.55),
    prot(3.4, IY[1] + 0.85, 1.3, 1.0, "orange", "4A"),
    txt(2.6, IY[1] - 2.3, "eIF4F = 4E + 4G + 4A", 0.65, "gray"),
    prot(8.8, IY[1] + 1.0, 4.6, 1.4, "#E6C3D6", "40S", "oval", 0.55), tr(8.6, IY[1], "UAC", aa="Met"),
    prot(10.6, IY[1] - 1.3, 2.2, 1.0, "green", "eIF2·GTP", "oval", 0.45), prot(6.9, IY[1] - 0.9, 1.4, 0.9, "yellow", "eIF3", "oval", 0.45),
    arr(11.4, IY[1] + 2.1, 14.4, IY[1] + 2.1, color=OR), txt(12.9, IY[1] + 2.9, "scanning (ATP)", 0.65, OR),
    txt(14.8, IY[1] - 0.8, "gccRccAUGG", 0.65, "#2F6A33"), txt(14.8, IY[1] + 0.9, "Kozak", 0.6, "gray"),
    fat(17.4, IY[1] - 0.5, 18.8, IY[1] - 0.5),
    ribo(22.6, IY[1], 6.6, True, "60S", "40S", hue="pink"), tr(22.6, IY[1], "UAC", aa="Met", s=0.75),
    txt(22.6, IY[1] + 3.0, "80S: eIF5 → GTP hydrolysis, eIF5B joins 60S", 0.68),
    note(IY[1] + 6.0, "Bacteria: IF3 keeps 30S apart from 50S and checks the start codon; IF1 blocks the A site; IF2·GTP brings fMet-tRNA_{fMet} to the P site. Start codons AUG (GUG, UUG).\n"
         "Methionine is formylated on the tRNA (N10-formyl-THF transformylase); deformylase and Met aminopeptidase later trim the N terminus.\n"
         "Eukaryotes: 43S = 40S + eIF1, 1A, 3 + eIF2·GTP·Met-tRNA_i; eIF4F + PABP activate the mRNA (closed loop) → 48S scans to the first AUG in a Kozak context.\n"
         "eIF2-GDP is recycled by eIF2B (GEF); eIF2α phosphorylation (PKR, PERK, GCN2, HRI) blocks this — global shut-off. IRESs (viruses) bypass the cap.", 0.76),
]}

# ================================================================ elongation cycle (notes p.67; lecture 4-181, 4-168)
EX, EYY = [5.5, 21.5], [6.0, 19.5]


def state(x, y, e=None, p=None, a=None):
    out = [mrna(x - 5.0, x + 5.0, y), ribo(x, y, 7.0, True)]
    for dx, spec in ((-SLOT, e), (0, p), (SLOT, a)):
        if spec:
            out.append(tr(x + dx, y, spec[0], pep=spec[1], **({"hue": "#D9D9D9"} if spec[1] == 0 else {})))
    return out


ELONG = {"items": [
    *state(EX[0], EYY[0], p=("UAC", 3)),
    *state(EX[1], EYY[0], p=("UAC", 3), a=("CGA", 1)),
    *state(EX[1], EYY[1], p=("UAC", 0), a=("CGA", 4)),
    *state(EX[0], EYY[1], e=("UAC", 0), p=("CGA", 4)),
    fat(11.0, EYY[0] - 1.2, 15.6, EYY[0] - 1.2), txt(13.3, EYY[0] - 5.6, "1. EF-Tu·GTP·aa-tRNA\nenters the A site", 0.72),
    tr(12.6, EYY[0] - 2.6, "CGA", aa="Ala", s=0.7), prot(14.2, EYY[0] - 3.1, 1.9, 1.0, "green", "EF-Tu·GTP", "oval", 0.42),
    txt(13.3, EYY[0] + 2.7, "codon–anticodon OK →\nGTP hydrolysis;\nEF-Tu·GDP leaves", 0.62, "gray"),
    fat(EX[1], EYY[0] + 3.0, EX[1], EYY[1] - 5.2), txt(EX[1] + 1.4, (EYY[0] + EYY[1]) / 2 - 1.0, "2. peptidyl transfer:\nα-NH_2 of aa-tRNA (A)\nattacks the peptidyl\nester (P); 23S rRNA", 0.72, None, "left"),
    fat(15.6, EYY[1] - 1.2, 11.0, EYY[1] - 1.2), txt(13.3, EYY[1] + 1.4, "3. EF-G·GTP:\ntranslocation by\none codon", 0.72),
    prot(13.3, EYY[1] - 3.0, 1.9, 1.0, "orange", "EF-G·GTP", "oval", 0.42),
    fat(EX[0], EYY[1] - 5.2, EX[0], EYY[0] + 3.0), txt(EX[0] - 1.4, (EYY[0] + EYY[1]) / 2 - 1.0, "4. tRNA leaves\nfrom E; A site\nfree for the\nnext aa-tRNA", 0.72, None, "right"),
    note(EYY[1] + 4.6, "Each cycle: 2 GTP (EF-Tu, EF-G) per peptide bond, ~15–20 aa/s in E. coli. EF-Ts (eEF1B) exchanges GDP for GTP on EF-Tu (eEF1A). EF-G (eEF2) mimics the EF-Tu·tRNA shape.\n"
         "Peptidyl transfer needs no protein or energy: the ester bond of peptidyl-tRNA pays for it (23S rRNA A2451 positions the substrates). The growing chain is handed from the P-site tRNA onto the A-site aminoacyl-tRNA (head growth).\n"
         "Fidelity: kinetic proofreading after GTP hydrolysis (accommodation). Error ~1/10^4. Diphtheria toxin ADP-ribosylates eEF2 (diphthamide).", 0.76),
]}

# ================================================================ termination, recycling, tmRNA (notes p.68, 70; lecture 4-190, 4-193)
TY = [5.5, 18.5]
TERM = {"items": [
    mrna(0.0, 10.0, TY[0]), ribo(5.0, TY[0], 7.0, True), tr(5.0, TY[0], "CGA", pep=5),
    {"g": "blob", "at": [5.0 + SLOT, TY[0] - 1.35], "w": 0.95, "h": 2.4, "shape": "box", "hue": "orange", "text": "RF1", "text_size": 0.42},
    txt(5.0 + SLOT, TY[0] + 3.0, "UAA in A site", 0.68, "red"),
    fat(10.6, TY[0] - 1.0, 12.6, TY[0] - 1.0),
    mrna(13.0, 23.0, TY[0]), ribo(18.0, TY[0], 7.0, True), tr(18.0, TY[0], "CGA", hue="#D9D9D9"),
    {"g": "blob", "at": [18.0 + SLOT, TY[0] - 1.35], "w": 0.95, "h": 2.4, "shape": "box", "hue": "orange", "text": "RF1", "text_size": 0.42},
    *[ball(20.4 + 0.55 * k, TY[0] - 5.4 + (0.3 if k % 2 else 0), "", h, 0.27) for k, h in enumerate(["orange", "green", "blue", "pink", "yellow"])],
    txt(23.4, TY[0] - 5.3, "polypeptide (COO⁻) released", 0.68, None, "left"),
    txt(18.0, TY[0] + 3.0, "H_2O hydrolyses peptidyl-tRNA", 0.68, "gray"),
    fat(23.6, TY[0] - 1.0, 25.6, TY[0] - 1.0),
    mrna(26.0, 35.0, TY[0] + 0.6),
    prot(30.0, TY[0] - 2.4, 6.4, 3.2, "purple", "50S", "oval", 0.6), prot(30.5, TY[0] + 1.6, 5.4, 1.6, "#C9BEDF", "30S", "oval", 0.55),
    txt(30.0, TY[0] + 3.3, "RRF + EF-G·GTP split the ribosome", 0.68, "gray"),
    txt(0.0, TY[0] + 4.4, "RF1 reads UAA, UAG; RF2 reads UAA, UGA (eRF1 reads all three). Its GGQ motif puts water into the peptidyl transferase centre. RF3·GTP (eRF3) removes RF1/2.", 0.72, None, "left"),
    txt(0.0, TY[1] - 6.0, "Rescue of stalled ribosomes: trans-translation by tmRNA (bacteria)", 0.95, None, "left"),
    mrna(0.0, 6.2, TY[1], ("5'", "3'")), ribo(5.0, TY[1], 7.0, True), tr(5.0, TY[1], "", pep=4),
    tr(5.0 + SLOT, TY[1], "", aa="Ala", hue="#F2C46D"), txt(8.8, TY[1] - 3.4, "tmRNA·SmpB\n(Ala-charged)", 0.68, "#B5651F", "left"),
    txt(4.0, TY[1] + 2.9, "mRNA without stop codon", 0.62, "gray"),
    fat(10.6, TY[1] - 1.0, 12.6, TY[1] - 1.0),
    rna([(13.0, TY[1]), (19.0, TY[1]), (23.0, TY[1])], ("", "")), rna([(19.0, TY[1]), (23.0, TY[1])], kind="#E0A100"),
    ribo(18.5, TY[1], 7.0, True), tr(18.5, TY[1], "", pep=6, aa_hue="red"),
    txt(23.6, TY[1], "ribosome switches to the tmRNA reading frame:\n11-aa ssrA tag (AANDENYALAA) + stop", 0.68, None, "left"),
    txt(23.6, TY[1] - 2.6, "tagged protein → ClpXP protease\nfaulty mRNA → RNase R", 0.68, "gray", "left"),
    note(TY[1] + 4.6, "Energy for a protein of N amino acids: 2N (aminoacylation, ATP → AMP + PP_i) + 1 (IF2) + (N−1) (EF-Tu) + (N−1) (EF-G) = 4N − 1 high-energy bonds;\n"
         "counting RF3 as well gives 4N. Eukaryotes: nonsense-mediated decay (premature stop before an exon junction complex), no-go and non-stop decay.", 0.76),
]}

# ================================================================ antibiotics / toxins acting on the ribosome (lecture 4-207)
RX, RYY = 14.0, 10.0
INHIB = {"items": [
    mrna(7.0, 21.0, RYY), ribo(RX, RYY, 9.0, True, "50S", "30S"),
    tr(RX - 0.27 * 9, RYY, "", hue="#D9D9D9"), tr(RX, RYY, "UAC", pep=3), tr(RX + 0.27 * 9, RYY, "CGA", aa="Ala"),
    pointer(RX, RYY - 5.5, "Macrolides (erythromycin): 50S exit\ntunnel, block elongation / translocation", 23.0, RYY - 6.8, "left"),
    pointer(RX + 1.0, RYY - 3.3, "Chloramphenicol, linezolid, clindamycin\n(50S peptidyl transferase centre)", 23.0, RYY - 3.9, "left"),
    pointer(RX + 0.27 * 9 + 0.3, RYY - 2.3, "Puromycin: mimics Tyr-tRNA end, enters A,\naccepts the peptide → premature release", 23.0, RYY - 1.0, "left"),
    pointer(RX + 0.27 * 9 + 0.6, RYY - 0.6, "Tetracyclines (30S A site):\nblock aa-tRNA binding", 23.0, RYY + 1.9, "left"),
    pointer(RX - 4.3, RYY - 2.2, "Fusidic acid: EF-G·GDP stuck;\ncycloheximide (80S): translocation", -2.0, RYY - 6.0, "left"),
    pointer(RX - 2.4, RYY + 2.0, "Streptomycin, aminoglycosides (30S):\nmisreading; high dose blocks initiation", -2.0, RYY + 4.0, "left"),
    note(RYY + 6.6, "Eukaryotic-only: cycloheximide (60S); diphtheria toxin and Pseudomonas exotoxin A ADP-ribosylate diphthamide (modified His) of eEF2 (NAD⁺ → nicotinamide);\n"
         "ricin and Shiga toxin are N-glycosidases that depurinate one A in 28S rRNA (sarcin–ricin loop). Chloramphenicol also hits mitochondrial ribosomes (gray baby syndrome).\n"
         "Aminoglycosides bind 16S rRNA near A1492/A1493 (decoding); tetracycline binds 16S; macrolides, chloramphenicol, linezolid bind 23S.", 0.76),
]}

# ================================================================ co-translational targeting to the ER (notes p.69; lecture 4-200, 4-201)
MY = 9.0
SRP = {"items": [
    {"g": "region", "x0": -0.6, "y0": MY + 1.6, "x1": 30.0, "y1": MY + 6.0, "hue": "yellow", "opacity": 0.18, "outline": False},
    txt(0.0, MY + 5.2, "ER lumen", 0.75, "gray", "left"), txt(29.6, -1.2, "cytosol", 0.75, "gray", "right"),
    {"g": "membrane", "y": MY, "x0": -0.6, "x1": 30.0, "t": 2.4},
    mrna(0.2, 5.8, 2.0, ("5'", "")), ribo(3.0, 2.0, 5.0, False, flip=True),
    rna([(3.0, 5.2), (3.4, 5.8), (2.8, 6.4), (3.2, 7.0)], kind="#D23A3A", w=0.32),
    txt(3.6, 6.9, "signal sequence", 0.68, "red", "left"),
    mrna(7.7, 13.3, 2.0, ("", "")), ribo(10.5, 2.0, 5.0, False, flip=True),
    rna([(10.5, 5.2), (10.9, 5.8), (10.3, 6.4)], kind="#D23A3A", w=0.32),
    prot(11.5, 6.6, 2.6, 1.2, "pink", "SRP", "oval", 0.55),
    txt(10.5, -0.6, "SRP binds the signal;\ntranslation pauses", 0.68, "gray"),
    mrna(14.7, 20.3, 4.6, ("", "")), ribo(17.5, 4.6, 5.0, False, flip=True),
    prot(14.2, MY - 1.75, 1.8, 1.0, "pink", "SRP", "oval", 0.45),
    {"g": "blob", "at": [14.2, MY], "w": 1.6, "h": 2.8, "shape": "box", "hue": "purple", "text": "SR", "text_size": 0.5},
    {"g": "channel", "at": [17.5, MY], "open": True, "hue": "teal"},
    txt(17.0, MY + 3.4, "SRP·SR (both GTP) dock on Sec61;\nGTP hydrolysis releases SRP", 0.66, "gray"),
    mrna(22.2, 27.8, 4.6, ("", "3'")), ribo(25.0, 4.6, 5.0, False, flip=True),
    {"g": "channel", "at": [25.0, MY], "open": True, "hue": "teal"},
    rna([(25.0, MY - 1.4), (25.0, MY + 2.0), (25.6, MY + 3.0), (26.4, MY + 2.6), (27.2, MY + 3.4)], kind="#5F4A8C", w=0.3),
    rna([(22.4, MY + 2.4), (22.9, MY + 3.2), (23.4, MY + 2.7)], kind="#D23A3A", w=0.32),
    prot(23.6, MY + 1.75, 1.6, 1.0, "orange", "SPase", "oval", 0.42),
    {"g": "ligand", "at": [27.6, MY + 4.4], "shape": "hexagon", "hue": "blue", "r": 0.35},
    {"g": "ligand", "at": [28.3, MY + 4.9], "shape": "hexagon", "hue": "green", "r": 0.35},
    {"g": "ligand", "at": [27.6, MY + 5.3], "shape": "hexagon", "hue": "green", "r": 0.35},
    prot(29.0, MY + 2.9, 1.6, 1.2, "brown", "BiP", "round", 0.45),
    txt(24.6, MY + 7.0, "translation resumes; signal peptidase cuts; OST adds the N-glycan (Asn-X-Ser/Thr) from dolichol-PP; BiP folds", 0.68),
    fat(5.9, 2.0, 7.5, 2.0), fat(13.4, 2.4, 14.9, 3.9), fat(20.4, 4.6, 22.0, 4.6),
    note(MY + 9.4, "SRP = 7SL RNA (Alu-like) + 6 proteins; SRP54 binds the hydrophobic signal (+ charge, 10–15 hydrophobic aa) and is a GTPase, as is the SRP receptor; GTP hydrolysis releases both.\n"
         "Proteins made this way: secreted, plasma-membrane, ER / Golgi, lysosomal. Tunicamycin blocks the first step of dolichol-PP-oligosaccharide synthesis (N-glycosylation).\n"
         "Bacteria: SRP (4.5S RNA) + FtsY → SecYEG co-translationally; post-translational SecB (chaperone) + SecA (ATPase) → SecYEG; folded proteins via Tat (twin-Arg signal).", 0.76),
]}

# ================================================================ ubiquitin–proteasome (lecture 4-197, 4-210)
UY = [3.0, 11.0]
UB = {"items": [
    ball(0.6, UY[0], "Ub", "yellow", 0.45, 0.45),
    prot(3.6, UY[0], 2.4, 1.6, "blue", "E1", "kidney", 0.6), txt(3.6, UY[0] - 1.7, "ATP → AMP + PP_i", 0.62, "gray"),
    ball(4.4, UY[0] - 0.75, "Ub", "yellow", 0.38, 0.4),
    fat(5.3, UY[0], 6.7, UY[0]),
    prot(8.4, UY[0], 2.0, 1.4, "teal", "E2", "round", 0.6), ball(9.1, UY[0] - 0.7, "Ub", "yellow", 0.38, 0.4),
    fat(9.8, UY[0], 11.2, UY[0]),
    prot(13.6, UY[0] + 0.2, 3.0, 2.0, "purple", "E3", "bilobe", 0.6), prot(16.4, UY[0] + 0.3, 2.4, 1.6, "gray", "substrate", "oval", 0.5),
    txt(13.6, UY[0] - 1.8, "E3 picks the substrate", 0.62, "gray"),
    fat(18.1, UY[0], 19.5, UY[0]),
    prot(22.0, UY[0] + 0.5, 2.6, 1.7, "gray", "substrate", "oval", 0.5),
    *[ball(23.0 + 0.55 * k, UY[0] - 0.7 - 0.6 * k, "Ub", "yellow", 0.34, 0.36) for k in range(4)],
    txt(25.8, UY[0] - 1.1, "K48-linked chain (≥ 4 Ub)", 0.68, None, "left"),
    down(22.0, UY[0] + 1.6, UY[1] - 2.8),
    prot(22.0, UY[1] - 2.0, 3.6, 1.2, "orange", "19S cap", "bilobe", 0.52),
    *[prot(22.0, UY[1] - 0.9 + 0.95 * k, 3.0, 0.9, "green" if k in (1, 2) else "dkgreen", ("α", "β", "β", "α")[k], "box", 0.5) for k in range(4)],
    prot(22.0, UY[1] + 3.0, 3.6, 1.2, "orange", "19S cap", "bilobe", 0.52, rot=180),
    txt(19.6, UY[1] - 1.9, "19S: binds the Ub chain, removes Ub (DUB),\nunfolds the protein (ATP)", 0.68, None, "right"),
    txt(19.6, UY[1] + 1.0, "20S core (α7β7β7α7): β1 caspase-like,\nβ2 trypsin-like, β5 chymotrypsin-like (N-terminal Thr)", 0.68, None, "right"),
    *[rna([(25.0 + 1.3 * k, UY[1] + 1.4 + 0.5 * (k % 2)), (25.9 + 1.3 * k, UY[1] + 1.6 + 0.5 * (k % 2))], kind="gray") for k in range(4)],
    txt(29.0, UY[1] + 3.2, "peptides (7–9 aa) → amino acids", 0.68),
    *[ball(25.4 + 0.9 * k, UY[1] - 0.6, "Ub", "yellow", 0.34, 0.36) for k in range(3)], txt(26.3, UY[1] - 1.6, "Ub recycled", 0.68, "gray"),
    note(UY[1] + 5.8, "E1 activates Ub (thioester on Cys, via Ub-AMP); E2 carries it; E3 (hundreds: RING / HECT; e.g. MDM2 → p53, SCF / APC-C → cyclins, parkin, VHL → HIF-1α) gives specificity.\n"
         "N-end rule: half-life depends on the N-terminal residue (Met, Ser, Ala long; Arg, Lys, Phe, Leu short). PEST sequences shorten half-life. K63 chains → signalling (NF-κB), not degradation.\n"
         "Bortezomib inhibits β5 (multiple myeloma). Lysosomes (autophagy) degrade long-lived and membrane proteins; the proteasome degrades short-lived, misfolded (ERAD) proteins.", 0.76),
]}

# ================================================================ tables
FACTOR_TABLE = table([
    ("Initiation", "bacteria IF1, IF2·GTP, IF3  ·  eukaryotes eIF1, 1A, 2·GTP, 2B (GEF), 3, 4F (4E cap-binding, 4G scaffold, 4A helicase), 5, 5B"),
    ("Elongation", "EF-Tu·GTP brings aa-tRNA (eEF1A) · EF-Ts exchanges GDP (eEF1B) · EF-G·GTP translocation (eEF2)"),
    ("Termination", "RF1 (UAA, UAG), RF2 (UAA, UGA), RF3·GTP  ·  eRF1 (all three stops), eRF3·GTP"),
    ("Recycling", "RRF + EF-G  ·  ABCE1"),
    ("Start signal", "Shine–Dalgarno (AGGAGG) pairs with 16S rRNA  ·  5′ cap + scanning to a Kozak AUG (gccRccAUGG)"),
    ("Initiator", "fMet-tRNA_{fMet}; start AUG (GUG, UUG)  ·  Met-tRNA_i; AUG only"),
    ("mRNA", "polycistronic, coupled to transcription  ·  monocistronic, exported first; circular polysome (cap–eIF4G–PABP)"),
    ("Ribosome", "70S = 30S + 50S  ·  80S = 40S + 60S"),
], w_name=6.4, x_desc=7.0, hue="purple")

ENERGY_TABLE = table([
    ("Aminoacylation", "aa + ATP + tRNA → aa-tRNA + AMP + PP_i (2 high-energy bonds)  ×  N  =  2N"),
    ("Initiation", "IF2·GTP (bacteria) / eIF2·GTP  =  1"),
    ("Elongation: EF-Tu", "one GTP per aa-tRNA delivered  =  N − 1"),
    ("Elongation: EF-G", "one GTP per translocation  =  N − 1"),
    ("Total", "4N − 1; 4N if RF3·GTP is counted (lecture). Lehninger: ≥ 4 per peptide bond"),
], w_name=7.4, x_desc=8.0, hue="orange")

SORT_TABLE = table([
    ("ER / secretion", "N-terminal signal sequence: 1–2 basic aa + 10–15 hydrophobic core; cleaved by signal peptidase (preproinsulin)"),
    ("ER retention", "C-terminal KDEL: escaped proteins are returned from the Golgi by the KDEL receptor"),
    ("Lysosome", "mannose-6-phosphate on N-glycans (GlcNAc-phosphotransferase in cis-Golgi; missing in I-cell disease)"),
    ("Mitochondrion", "N-terminal amphipathic α-helix (Arg / Lys on one face); TOM / TIM import, cleaved in the matrix"),
    ("Nucleus", "internal NLS, basic cluster (PKKKRKV, SV40 T); importin + Ran-GTP; not cleaved"),
    ("Peroxisome", "C-terminal PTS1 (SKL, most) or N-terminal PTS2; folded proteins imported (Zellweger syndrome)"),
    ("Membrane", "stop-transfer / signal-anchor hydrophobic segments; GPI anchor added in the ER"),
], w_name=6.8, x_desc=7.4, hue="teal")

PTM_TABLE = table([
    ("Proteolysis", "signal peptide removal; proinsulin → insulin (C-peptide); zymogens: trypsinogen → trypsin, fibrinogen → fibrin"),
    ("N / C terminus", "deformylation, Met removal; N-acetylation; N-myristoylation (Gly); C-terminal prenylation of CaaX (Ras)"),
    ("Phosphorylation", "Ser, Thr, Tyr (kinases / phosphatases); His, Asp in bacterial two-component systems"),
    ("Glycosylation", "N-linked on Asn-X-Ser/Thr (ER, dolichol-PP) · O-linked on Ser / Thr (Golgi)"),
    ("Hydroxylation", "Pro, Lys in collagen (vitamin C); Pro of HIF-1α (O_2 sensor → VHL → proteasome)"),
    ("Carboxylation", "γ-carboxy-Glu in factors II, VII, IX, X, protein C (vitamin K; warfarin blocks)"),
    ("Methyl / acetyl", "Lys, Arg of histones; trimethyl-Lys in calmodulin, cytochrome c"),
    ("Disulfides, folding", "Cys–Cys in the ER (PDI); Pro cis–trans (PPIase); chaperones Hsp70, chaperonin GroEL / TRiC"),
    ("Ub, SUMO, ADP-ribose", "ubiquitin K48 → proteasome; SUMO → localisation; toxins ADP-ribosylate eEF2, Gs α, Gi α"),
    ("Other", "iodination of Tyr (thyroglobulin → T_3, T_4); biotin / lipoate / heme attachment; protein splicing (inteins)"),
], w_name=8.6, x_desc=9.2, hue="pink")

# ================================================================ molecules and reactions
MOLECULES = [
    {"id": "alanyl_adenylate", "name": "Alanyl-adenylate (Ala-AMP)", "label": "Alanyl-AMP", "aliases": ["aminoacyl-AMP", "aminoacyl adenylate"],
     "smiles": "C[C@H]([NH3+])C(=O)OP(=O)([O-])O[4*]", "ring": {"coordgen": True}, "tags": ["translation"],
     "notes": "胺醯-AMP：胺基酸的羧基與 AMP 的磷酸形成混合酸酐（高能），由 aaRS 結合在活性中心。Ado = adenosine。"},
    {"id": "alanyl_trna", "name": "Alanyl-tRNA", "label": "Ala-tRNA^{Ala}", "aliases": ["aminoacyl-tRNA", "charged tRNA"],
     "smiles": "C[C@H]([NH3+])C(=O)O[20*]", "ring": {"coordgen": True}, "tags": ["translation"],
     "notes": "胺醯-tRNA：胺基酸以酯鍵接在 tRNA 3′ 端 A76 的 3′-OH（class I 先接 2′-OH 再轉位）。tRNA = 其餘 tRNA。"},
    {"id": "fmet_trna", "name": "N-Formylmethionyl-tRNA", "label": "fMet-tRNA^{fMet}", "aliases": ["fMet-tRNA", "initiator tRNA"],
     "smiles": "O=CN[C@@H](CCSC)C(=O)O[20*]", "ring": {"coordgen": True}, "tags": ["translation"],
     "notes": "細菌起始 tRNA：Met-tRNA^{fMet} 由 transformylase 以 N10-formyl-THF 甲醯化。甲醯基使它只能進 P 位、不被 EF-Tu 接受。"},
    {"id": "fmet_ala_trna", "name": "fMet-Ala-tRNA (dipeptidyl-tRNA)", "label": "fMet-Ala-tRNA^{Ala}", "aliases": ["dipeptidyl-tRNA", "peptidyl-tRNA"],
     "smiles": "O=CN[C@@H](CCSC)C(=O)N[C@@H](C)C(=O)O[20*]", "ring": {"coordgen": True}, "tags": ["translation"],
     "notes": "第一個胜肽鍵形成後的 dipeptidyl-tRNA，位於 A 位，之後由 EF-G 轉位到 P 位。"},
    {"id": "puromycin", "name": "Puromycin", "label": "Puromycin", "aliases": ["嘌呤黴素"],
     "smiles": "COc1ccc(C[C@H](N)C(=O)N[C@@H]2[C@H](O)[C@@H](n3cnc4c(N(C)C)ncnc43)O[C@H]2CO)cc1", "ring": {"coordgen": True, "rotate": 180},
     "tags": ["translation"],
     "notes": "模仿 tyrosyl-tRNA 的 3′ 端（O-甲基酪胺酸 + N6,N6-二甲基腺苷），但以醯胺鍵（不可水解）取代酯鍵。進入 A 位接下胜肽後脫落 → 提早終止；原核、真核都抑制。"},
    {"id": "tyrosyl_trna_end", "name": "Tyrosyl-tRNA 3′ end", "label": "Tyrosyl-tRNA (3′ end)", "aliases": ["Tyr-tRNA"],
     "smiles": "N[C@@H](Cc1ccc(O)cc1)C(=O)O[C@@H]1[C@H](O)[C@@H](n2cnc3c(N)ncnc32)O[C@H]1COP(=O)([O-])O[20*]", "ring": {"coordgen": True, "rotate": 180},
     "tags": ["translation"], "notes": "tRNA 3′ 端的 A76 與酯鍵連接的 Tyr（tRNA = 其餘 tRNA）。與 puromycin 對照。"},
]

REACTIONS = [
    {"id": "aars_activation", "title": "Amino acid activation (aminoacyl-tRNA synthetase, step 1)", "substrates": ["alanine"], "products": ["alanyl_adenylate"],
     "cofactors_in": ["ATP"], "cofactors_out": ["PPi"], "extra_in": ["H+"], "enzyme": "alanyl-tRNA synthetase", "tags": ["translation"],
     "card": {"type": "羧酸根的 O（洋紅）攻擊 ATP 的 α-P → 混合酸酐 aa~AMP；β、γ 磷酸以 PPi 離開（青）",
              "exam_points": ["PPi 再被 pyrophosphatase 水解 → 每個胺基酸耗 2 個高能鍵", "aa~AMP 不離開酵素，直接進行第二步"]}},
    {"id": "aars_transfer", "title": "Aminoacyl transfer to tRNA (aminoacyl-tRNA synthetase, step 2)", "substrates": ["alanyl_adenylate"], "products": ["alanyl_trna"],
     "cofactors_in": ["tRNA"], "cofactors_out": ["AMP"], "extra_out": ["H+"], "enzyme": "alanyl-tRNA synthetase (class II)", "tags": ["translation"],
     "card": {"type": "tRNA 3′-OH 的 O（青）攻擊醯基碳，形成胺醯-tRNA 酯鍵；AMP 離開",
              "exam_points": ["class I（如 Ile、Val、Leu、Met）接 2′-OH，class II（如 Ala、Gly、Ser、Phe）接 3′-OH",
                              "tRNA^{Ala} 的 G3·U70 是 AlaRS 的辨識要素（second genetic code）"]}},
    {"id": "peptidyl_transfer", "title": "Peptidyl transfer (first peptide bond)", "substrates": ["fmet_trna", "alanyl_trna"], "products": ["fmet_ala_trna"],
     "cofactors_in": [], "cofactors_out": ["tRNA"], "extra_out": ["H+"], "enzyme": "peptidyl transferase\n(23S rRNA, ribozyme)", "tags": ["translation"],
     "card": {"type": "A 位 aa-tRNA 的 α-NH₂（洋紅 N）攻擊 P 位 fMet-tRNA 的酯羰基（青 C=O）→ 胜肽轉移到 A 位 tRNA，P 位留下去醯化 tRNA",
              "mechanism": "不需蛋白質與額外能量：酯鍵的能量推動；23S rRNA（A2451）定位受質，P 位 tRNA 的 2′-OH 參與質子轉移",
              "exam_points": ["胜肽鏈由 P 轉到 A（大分子轉給小分子），與脂肪酸合成相同模式", "Chloramphenicol 抑制此步（50S）"]}},
]

AMP_PART = "C(=O)[O]-[P:1](=[O:2])([O-:3])[O:4][#0]"      # from ATP (the bridging O came from the carboxylate)
AMP_LEAVING = "C(=O)[O:5]-[P:1](=[O:2])([O-:3])[O:4][#0]"  # leaves as AMP in step 2
COLORS.update({
    "aars_activation": {"alanine": {"colors": {0: {"c": [["black", "CO"], [M, "O^-"]]}}},
                        "alanyl_adenylate": {"highlight": [{"smarts": "C(=O)[O:1]P", "color": M}, {"smarts": AMP_PART, "color": C}]},
                        "ATP": frag((C, "ATP")), "PPi": frag((C, "PP_i"))},
    "aars_transfer": {"alanyl_adenylate": hl(AMP_LEAVING), "AMP": frag((M, "AMP")),
                      "tRNA": frag((C, "tRNA")), "alanyl_trna": hl("C(=O)[O:1]*", col=C)},
    "peptidyl_transfer": {"alanyl_trna": hl("[NH3+:1]"), "fmet_trna": hl("[C:1](=[O:2])O*", col=C),
                          "fmet_ala_trna": {"highlight": [{"smarts": "C(=O)[NH:1]C(C)C(=O)O*", "color": M},
                                                          {"smarts": "[C:1](=[O:2])[NH]C(C)C(=O)O*", "color": C}]},
                          "tRNA": frag((C, "tRNA"))},
})

NEW_COFACTORS = {"tRNA": {"label": "tRNA", "smiles": "O[20*]"}}

PATHWAYS = [
    {"id": "genetic_code", "title": "The genetic code", "layout": "scene", "scene": GENETIC_CODE, "tags": ["translation"],
     "notes": "三聯體、不重疊、無標點、5′→3′、簡併（61 個有義密碼子對 20 種胺基酸）、近乎通用（粒線體 UGA = Trp、AUA = Met、AGA/AGG = stop）。"
              "AUG = 起始兼 Met；UAA、UAG、UGA = 終止。單一密碼子：Met、Trp；六個：Leu、Ser、Arg。Selenocysteine 用 UGA（+ SECIS）、pyrrolysine 用 UAG。"},
    {"id": "wobble_pairing", "title": "Codon–anticodon pairing and wobble", "layout": "scene", "scene": WB, "tags": ["translation"],
     "notes": "密碼子第 3 位與反密碼子第 1 位（5′ 端）可擺動：G–C/U、U–A/G、I–U/C/A；C、A 只認一種。前兩位嚴格配對。"
              "≥ 32 種 tRNA 即可讀 61 個密碼子；同一胺基酸的不同 tRNA = isoaccepting tRNA。"},
    {"id": "aminoacyl_trna_synthetase", "title": "Aminoacyl-tRNA synthetase: charging tRNA", "layout": "scene", "scene": AARS, "tags": ["translation"],
     "notes": "兩步：aa + ATP → aa~AMP + PPi；aa~AMP + tRNA → aa-tRNA + AMP（每個胺基酸 2 個高能鍵）。Class I 接 2′-OH、class II 接 3′-OH。"
              "辨識 anticodon 與 acceptor stem / discriminator base（second genetic code，tRNA^{Ala} 的 G3·U70）。Double-sieve 校對（IleRS 剔除 Val）。"},
    {"id": "ribosome_structure", "title": "The ribosome: subunits and E / P / A sites", "layout": "scene", "scene": RB, "tags": ["translation"],
     "notes": "70S = 30S（16S + 21 蛋白）+ 50S（23S + 5S + 33 種蛋白，L7/L12 有 4 份共 36 個）；80S = 40S（18S）+ 60S（28S、5.8S、5S）。A = aminoacyl-tRNA、P = peptidyl-tRNA、E = exit。"
              "小次單元先結合 mRNA 並負責解碼；大次單元的 23S/28S rRNA 形成 peptidyl transferase center（ribozyme）。"
              "證明：puromycin 只與 P 位的 peptidyl-tRNA 反應；一個核糖體可結合 3 個去醯化 tRNA 但只能 2 個胺醯 tRNA → E 位。"},
    {"id": "translation_initiation", "title": "Translation initiation: bacteria vs eukaryotes", "layout": "scene", "scene": INIT, "tags": ["translation"],
     "notes": "原核：30S + IF3 + IF1 + IF2·GTP·fMet-tRNA + mRNA（SD 與 16S rRNA 3′ 端配對）→ 50S 加入、IF2 水解 GTP、IF 離開 → 70S 起始複合體（fMet-tRNA 在 P）。"
              "真核：43S（40S + eIF1/1A/3 + eIF2·GTP·Met-tRNAi）+ eIF4F（4E 認 cap、4G 支架、4A helicase）與 PABP 活化的 mRNA → 48S 以 ATP 掃描至 Kozak 序列中的 AUG"
              " → eIF5 使 eIF2 水解 GTP、eIF5B 帶 60S → 80S。eIF2α 磷酸化 → 抑制 eIF2B → 整體轉譯關閉。"},
    {"id": "translation_elongation", "title": "The elongation cycle", "layout": "scene", "scene": ELONG, "tags": ["translation"],
     "notes": "① EF-Tu·GTP·aa-tRNA 進 A 位，密碼子正確 → GTP 水解、EF-Tu·GDP 離開（EF-Ts 換 GTP）；② peptidyl transfer：A 位胺基攻擊 P 位酯鍵（23S rRNA）；"
              "③ EF-G·GTP 轉位：P→E、A→P、mRNA 前進 3 nt；④ E 位 tRNA 離開。每個胜肽鍵 2 GTP。真核 eEF1A、eEF1B、eEF2（白喉毒素 ADP-ribosylation）。"},
    {"id": "translation_termination", "title": "Termination, ribosome recycling and tmRNA rescue", "layout": "scene", "scene": TERM, "tags": ["translation"],
     "notes": "終止密碼子進 A 位 → RF1（UAA、UAG）或 RF2（UAA、UGA）結合（eRF1 讀全部），GGQ 使水攻擊 → 胜肽釋出；RF3·GTP（eRF3）移除 RF1/2；RRF + EF-G 拆開核糖體（真核 ABCE1）。"
              "RF 與 tRNA、EF-G 與 EF-Tu·tRNA 為分子擬態。tmRNA·SmpB 拯救停在無終止密碼子 mRNA 上的核糖體：加 11-aa ssrA tag → ClpXP 分解；mRNA 被分解。"
              "合成 N 個胺基酸至少 4N−1 個高能鍵。"},
    {"id": "translation_inhibitors", "title": "Antibiotics and toxins that act on the ribosome", "layout": "scene", "scene": INHIB, "tags": ["translation"],
     "notes": "原核：streptomycin（30S，誤讀、高濃度阻斷起始）、tetracycline（30S A 位）、chloramphenicol（50S peptidyl transferase，也抑制粒線體）、erythromycin（50S 出口通道）、"
              "linezolid、clindamycin、fusidic acid（EF-G）。兩者：puromycin（仿 Tyr-tRNA，提早終止）。真核：cycloheximide（60S）、白喉毒素（eEF2 ADP-ribosylation）、ricin（28S rRNA 去嘌呤）。"},
    {"id": "srp_er_targeting", "title": "Co-translational targeting to the ER (SRP cycle)", "layout": "scene", "scene": SRP, "tags": ["translation"],
     "notes": "N 端 signal sequence 冒出 → SRP（7SL RNA + 蛋白，GTPase）結合、轉譯暫停 → SRP receptor（GTP）把核糖體送到 Sec61 translocon → SRP 水解 GTP 離開、轉譯恢復"
              " → 胜肽穿膜進入內腔 → signal peptidase 切除訊號 → OST 由 dolichol-PP 轉移 N-glycan 到 Asn-X-Ser/Thr（tunicamycin 抑制）→ BiP 協助摺疊 → Golgi。"
              "細菌：SRP（4.5S RNA）/FtsY、SecB + SecA → SecYEG、Tat（已摺疊）。"},
    {"id": "ubiquitin_proteasome", "title": "Ubiquitin–proteasome pathway", "layout": "scene", "scene": UB, "tags": ["translation"],
     "notes": "E1（ATP → AMP + PPi，Ub 與 Cys 形成硫酯）→ E2 → E3（受質專一性）→ K48 多聚泛素（≥ 4）→ 26S proteasome（19S 認 Ub、去泛素、ATP 解摺疊；20S 核心 β1、β2、β5 蛋白酶）"
              " → 7–9 aa 胜肽，Ub 回收。N-end rule、PEST 決定半衰期。K63 鏈為訊號用途。Bortezomib 抑制 β5。"},
    {"id": "translation_factors", "title": "Translation factors: bacteria vs eukaryotes", "layout": "scene", "scene": FACTOR_TABLE, "tags": ["translation"],
     "notes": "起始 IF1-3 vs eIF1、1A、2、2B、3、4F、5、5B；延長 EF-Tu/EF-Ts/EF-G vs eEF1A/eEF1B/eEF2；終止 RF1/RF2/RF3 vs eRF1/eRF3；回收 RRF vs ABCE1。"
              "SD vs cap + Kozak；fMet vs Met；70S vs 80S。"},
    {"id": "translation_energy_cost", "title": "Energy cost of protein synthesis", "layout": "scene", "scene": ENERGY_TABLE, "tags": ["translation"],
     "notes": "N 個胺基酸的蛋白質：活化 2N + 起始 1 + EF-Tu (N−1) + EF-G (N−1) = 4N − 1 個高能磷酸鍵；講義另計 RF3 為 4N。Lehninger 只寫「每個胜肽鍵至少 4 個高能磷酸當量」（6e p.1139）。"},
    {"id": "protein_sorting_signals", "title": "Protein sorting signals", "layout": "scene", "scene": SORT_TABLE, "tags": ["translation"],
     "notes": "ER/分泌：N 端疏水 signal sequence（切除）；KDEL 留 ER；M6P → 溶體（I-cell disease）；粒線體：N 端兩性 α 螺旋；細胞核：內部 NLS（不切除）；"
              "過氧化體：C 端 SKL（PTS1）。共轉譯（SRP）vs 轉譯後（粒線體、核、過氧化體）運送。"},
    {"id": "post_translational_modifications", "title": "Post-translational modifications", "layout": "scene", "scene": PTM_TABLE, "tags": ["translation"],
     "notes": "切割（signal peptide、proinsulin、zymogen）、N/C 端修飾、磷酸化、醣基化（N：Asn-X-Ser/Thr 在 ER；O：Ser/Thr 在 Golgi）、羥化（collagen，維生素 C）、"
              "γ-羧化（凝血因子，維生素 K）、甲基化／乙醯化、雙硫鍵（PDI）、prenylation、泛素化、ADP-ribosylation、碘化。"},
    {"id": "puromycin_vs_trna", "title": "Puromycin mimics the 3′ end of tyrosyl-tRNA", "layout": "free", "cell_size": 15.0,
     "nodes": [{"id": "pur", "ref": "puromycin", "pos": [0, 0], "decorate": {"highlight": [{"smarts": "[C:1](=[O:2])[N:3][C;R]", "color": "magenta"}]}},
               {"id": "tyr", "ref": "tyrosyl_trna_end", "pos": [1.25, 0], "decorate": {"highlight": [{"smarts": "[C:1](=[O:2])[O:3][C;R]", "color": "magenta"}]}}],
     "edges": [], "text": False, "tags": ["translation"],
     "notes": "洋紅：puromycin 的醯胺鍵 vs tyrosyl-tRNA 的酯鍵。Puromycin 進 A 位接受胜肽後，醯胺鍵不能再轉移 → 胜肽-puromycin 脫落（提早終止）；只與 P 位的 peptidyl-tRNA 反應，因此用來定義 A、P 位。"},
]


def add_cofactors(data):
    p = data / "cofactors.yaml"
    text = p.read_text(encoding="utf-8")
    new = [k for k in NEW_COFACTORS if f"\n{k}:" not in text]
    if new:
        lines = ["", "# --- 第二十一輪（轉譯）：tRNA 3′ 端 A76 的 OH（[20*] = 其餘 tRNA；由 tools/gen_mb_tl.py 產生）"]
        for k in new:
            v = NEW_COFACTORS[k]
            lines.append(f'{k}: {{label: "{v["label"]}", smiles: "{v["smiles"]}"}}')
        p.write_text(text.rstrip("\n") + "\n" + "\n".join(lines) + "\n", encoding="utf-8")


def main(data):
    data = Path(data)
    add_cofactors(data)
    for m in MOLECULES:
        d = {"id": m["id"], "name": m["name"], "label": m["label"], "aliases": m["aliases"], "smiles": m["smiles"]}
        if "ring" in m:
            d["ring"] = m["ring"]
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
        d = tidy(with_notes(d))
        d["parts_png"] = False
        d["status"] = "pending"
        dump(data / "pathways" / f"{p['id']}.yaml", d)
    print(f"{len(MOLECULES)} molecules, {len(REACTIONS)} reactions, {len(PATHWAYS)} pathways")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data")
