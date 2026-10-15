"""Round 26: gaps found in the biochemistry notes (typora_file/生化 B1/B2 code blocks): serine-protease catalysis,
chaperonin cycle, phospholipase cleavage sites, protein purification / sequencing, PDH complex, RBC glutathione
cycle, THF one-carbon forms and the starvation fuel curve (Lehninger 6e Fig. 23–31).
Figures carry names and short labels only; explanations live in notes_zh.

    python tools/gen_gap.py data
"""
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent))
from figure_text import tidy  # noqa: E402
from gen_illus_sig import A, B, FAT, I, INH, LINE, T  # noqa: E402
from gen_mb_repair import table  # noqa: E402

M, C = "magenta", "cyan"

# ================================================================ cofactors, molecules

NEW_COFACTORS = {
    "EnzSer": {"label": "Ser^{195}–OH (enzyme)", "smiles": "OC[22*]"},
}

PHE = "N[C@@H](Cc1ccccc1)"
MOLECULES = [
    {"id": "chymo_substrate", "name": "Peptide substrate (Phe–X bond)", "label": "Peptide (…Phe–X…)",
     "aliases": ["chymotrypsin substrate"], "smiles": f"[23*]C(=O){PHE}C(=O)N[24*]", "ring": {"coordgen": False},
     "tags": ["catalysis"],
     "notes": "Chymotrypsin 切在芳香族胺基酸（Phe、Trp、Tyr）羧基端的胜肽鍵。Pep_N、Pep_C = 受質胜肽 N 端側、C 端側的其餘部分。"},
    {"id": "acyl_enzyme", "name": "Acyl-enzyme intermediate (chymotrypsin)", "label": "Acyl-enzyme",
     "aliases": ["acyl-enzyme", "醯基-酵素中間物"], "smiles": f"[23*]C(=O){PHE}C(=O)OC[22*]", "ring": {"coordgen": False},
     "tags": ["catalysis"],
     "notes": "受質的醯基以酯鍵接在 Ser195 的 O 上（共價催化）；胺基端產物已離開。Enz = 酵素其餘部分。"},
    {"id": "peptide_c_fragment", "name": "C-terminal fragment (new amino end)", "label": "H_2N–Pep_C",
     "aliases": ["first product"], "smiles": "N[24*]", "ring": {"coordgen": True}, "tags": ["catalysis"],
     "notes": "第一個產物：帶新胺基端的 C 端片段（醯化步驟釋出）。"},
    {"id": "peptide_n_fragment", "name": "N-terminal fragment (new carboxyl end)", "label": "Pep_N–Phe–COOH",
     "aliases": ["second product"], "smiles": f"[23*]C(=O){PHE}C(=O)O", "ring": {"coordgen": True}, "tags": ["catalysis"],
     "notes": "第二個產物：帶新羧基端的 N 端片段（去醯化步驟由水解釋出）。"},
    {"id": "formyl_thf", "name": "N10-Formyl-THF", "label": "N^{10}-Formyl-THF", "aliases": ["10-formyltetrahydrofolate", "CHO-THF"],
     "smiles": "Nc1nc2NC[C@H](CN(C=O)c8ccc(cc8)C(=O)N[C@@H](CCC(=O)[O-])C(=O)[O-])Nc2c(=O)[nH]1",
     "ring": {"coordgen": False}, "highlight": "[N:1][CH:2]=[O:3]", "tags": ["one_carbon", "coenzyme"],
     "notes": "最氧化的一碳單位（甲酸層級）。提供嘌呤環的 C2 與 C8；細菌起始 tRNA 的 fMet 也由它甲醯化。"},
    {"id": "methenyl_thf", "name": "N5,N10-Methenyl-THF", "label": "N^5,N^{10}-Methenyl-THF", "aliases": ["5,10-methenyltetrahydrofolate"],
     "smiles": "Nc1nc2NC[C@@H]3C[N+](c8ccc(cc8)C(=O)N[C@@H](CCC(=O)[O-])C(=O)[O-])=CN3c2c(=O)[nH]1",
     "ring": {"coordgen": False}, "highlight": "[N+:1]=[CH:2][N:3]", "tags": ["one_carbon", "coenzyme"],
     "notes": "與 N¹⁰-formyl-THF 同一氧化層級（cyclohydrolase 互換）；由 methylene-THF 經 NADP⁺ 依賴的去氫酶氧化而來。"},
]

# ================================================================ reactions

REACTIONS = [
    {"id": "chymotrypsin_acylation", "title": "Serine protease step 1: acylation (chymotrypsin)",
     "substrates": ["chymo_substrate"], "products": ["acyl_enzyme", "peptide_c_fragment"],
     "cofactors_in": ["EnzSer"], "cofactors_out": [], "enzyme": "chymotrypsin (Ser195, His57, Asp102)", "tags": ["catalysis"],
     "decorate": {"EnzSer": {"colors": {0: {"c": [[C, "Ser^{195}–OH (enzyme)"]]}}},
                  "chymo_substrate": {"highlight": [{"smarts": "C(=O)[NH:1][#0]", "color": M}]},
                  "acyl_enzyme": {"highlight": [{"smarts": "C(=O)[O:1]C[#0]", "color": C}]},
                  "peptide_c_fragment": {"highlight": [{"smarts": "[NH2:1]", "color": M}]}},
     "card": {"type": "Ser195 的 O（青）親核攻擊胜肽鍵羰基 → 四面體中間物（氧陰離子洞穩定）→ 醯基-酵素；醯胺 N（洋紅）帶著 C 端片段離開",
              "exam_points": ["催化三聯體：Asp102 定位並穩定 His57，His57 當廣義鹼拔 Ser195 的 H（之後當廣義酸給離去的 N）",
                              "氧陰離子洞：Gly193 與 Ser195 主鏈 NH 以氫鍵穩定四面體中間物的 O⁻",
                              "Burst kinetics：第一個產物（胺）快速釋出，去醯化較慢（速率限制）"]}},
    {"id": "chymotrypsin_deacylation", "title": "Serine protease step 2: deacylation (chymotrypsin)",
     "substrates": ["acyl_enzyme"], "products": ["peptide_n_fragment"],
     "cofactors_in": ["H2O"], "cofactors_out": ["EnzSer"], "enzyme": "chymotrypsin (His57 activates water)", "tags": ["catalysis"],
     "decorate": {"H2O": {"colors": {0: {"c": [[M, "H_2O"]]}}},
                  "acyl_enzyme": {"highlight": [{"smarts": "C(=O)[O:1]C[#0]", "color": C}]},
                  "peptide_n_fragment": {"highlight": [{"smarts": "C(=O)[OH:1]", "color": M}]},
                  "EnzSer": {"colors": {0: {"c": [[C, "Ser^{195}–OH (enzyme)"]]}}}},
     "card": {"type": "His57 活化水，水的 O（洋紅）攻擊醯基-酵素酯鍵 → 第二個四面體中間物 → 羧酸產物；Ser195–OH（青）復原",
              "exam_points": ["DFP（diisopropyl fluorophosphate）共價修飾 Ser195 → 不可逆抑制（找出活性 Ser）",
                              "TPCK 烷化 His57（chymotrypsin），TLCK 對 trypsin"]}},
]

# ================================================================ scenes

# --- serine protease: catalytic triad + oxyanion hole, and S1 pockets
SERPROT = {"items": [
    T(0.0, -0.4, "Catalytic triad", 0.85, None, "left"),
    B(2.0, 6.0, 2.8, 1.5, "red", "Asp^{102}", ts=0.55), B(7.0, 6.0, 2.8, 1.5, "blue", "His^{57}", ts=0.55),
    B(12.0, 6.0, 2.8, 1.5, "green", "Ser^{195}", ts=0.55, accent=True),
    A(3.5, 6.0, 5.5, 6.0, "gray", dashed=True), T(4.5, 5.0, "H-bond", 0.55, "gray"),
    A(10.6, 6.6, 8.5, 6.6, "cyan", 0.0), T(9.6, 7.6, "H^+", 0.6, "cyan"),
    B(12.0, 2.0, 3.4, 1.4, "gray", "C=O (peptide)", ts=0.5), A(12.0, 5.2, 12.0, 2.8, M), T(12.6, 3.9, "O^- attacks", 0.55, M, "left"),
    B(16.6, 0.8, 2.2, 1.0, "yellow", "Gly^{193} NH", ts=0.4), B(16.6, 3.2, 2.2, 1.0, "yellow", "Ser^{195} NH", ts=0.4),
    A(15.4, 1.0, 13.8, 1.8, "gray", dashed=True), A(15.4, 3.0, 13.8, 2.3, "gray", dashed=True),
    T(16.6, 4.6, "oxyanion hole", 0.6),
    # S1 specificity pockets
    T(21.0, -0.4, "S1 pocket", 0.85, None, "left"),
    B(23.0, 3.0, 3.4, 3.0, "teal"), I(23.0, 3.6, "", "yellow", 0.75), T(23.0, 6.0, "Chymotrypsin", 0.65),
    T(23.0, 7.0, "Phe, Trp, Tyr", 0.6, "gray"),
    B(28.0, 3.0, 3.4, 3.0, "blue"), I(28.0, 3.6, "", "yellow", 0.55), I(28.0, 4.5, "–", "red", 0.3, 0.45), T(28.0, 6.0, "Trypsin", 0.65),
    T(28.0, 7.0, "Lys, Arg (Asp^{189})", 0.6, "gray"),
    B(33.0, 3.0, 3.4, 3.0, "purple"), I(33.0, 3.4, "", "yellow", 0.3), T(33.0, 6.0, "Elastase", 0.65),
    T(33.0, 7.0, "Ala, Gly (Val, Thr)", 0.6, "gray"),
]}

# --- GroEL / GroES chaperonin cycle
def groel(x, y, cap=False, closed=False):
    """Side view of the GroEL double ring (two walls per ring); GroES lid on top when `cap`."""
    w = 2.6 if cap else 2.0          # the capped (cis) ring opens wider
    out = []
    for dy in (0.0, 2.8):
        ww = w if dy == 0 else 2.0
        out += [B(x - ww, y + dy, 1.1, 2.6, "green", shape="box", accent=(cap and dy == 0)),
                B(x + ww, y + dy, 1.1, 2.6, "green", shape="box", accent=(cap and dy == 0))]
    if cap:
        out.append(B(x, y - 1.8, 5.0, 1.1, "yellow", "GroES", ts=0.45))
    return out


def unfolded(x, y):
    return {"g": "strand", "kind": "rna", "width": 0.16,
            "points": [[x - 1.0, y], [x - 0.6, y - 0.5], [x - 0.2, y + 0.3], [x + 0.2, y - 0.5], [x + 0.6, y + 0.3], [x + 1.0, y - 0.2]]}


GROEL = {"items": [
    *groel(3.0, 4.0), unfolded(3.0, 1.0), A(3.0, 1.6, 3.0, 3.0, "black"),
    T(3.0, 9.4, "1. hydrophobic rim\nbinds unfolded protein", 0.6),
    FAT(6.6, 5.4, 8.6, 5.4),
    *groel(12.4, 4.0, cap=True), unfolded(12.4, 4.0), *[I(9.2 + 0.0 * k, 2.4 + 0.6 * k, "", "pink", 0.22) for k in range(0)],
    I(15.8, 2.2, "7 ATP", "pink", 0.6, 0.32),
    T(12.4, 9.4, "2. 7 ATP + GroES lid:\nhydrophilic cavity", 0.6),
    FAT(16.8, 5.4, 18.8, 5.4),
    *groel(22.6, 4.0, cap=True), B(22.6, 4.0, 1.6, 1.4, "red", ts=0.4), I(26.0, 2.2, "ADP", "gray", 0.55, 0.34),
    T(22.6, 9.4, "3. folds inside (~10 s);\nATP hydrolysed", 0.6),
    FAT(27.0, 5.4, 29.0, 5.4),
    *groel(32.4, 4.0), B(32.4, 0.6, 1.6, 1.4, "red"), A(32.4, 3.0, 32.4, 1.6, "black"),
    T(32.4, 9.4, "4. ATP on trans ring:\nGroES and protein leave", 0.6),
]}

# --- PDH complex
PDH = {"items": [
    B(3.0, 6.0, 4.0, 3.0, "blue", "E1\n(TPP)", ts=0.6, accent=True),
    B(12.0, 9.0, 5.0, 3.0, "green", "E2\n(lipoamide)", ts=0.6),
    B(21.0, 6.0, 4.0, 3.0, "orange", "E3\n(FAD)", ts=0.6),
    LINE([[12.0, 7.5], [12.0, 4.6]], 0.18), I(12.0, 4.0, "Lip", "yellow", 0.6, 0.4),
    A(11.2, 3.8, 5.4, 5.0, "gray", 0.25), A(12.8, 3.8, 18.6, 5.0, "gray", -0.25), T(12.0, 2.2, "swinging arm (Lys–lipoate)", 0.6, "gray"),
    I(3.0, 1.2, "pyruvate", "orange", 0.9, 0.36), A(3.0, 2.1, 3.0, 4.4), I(-0.6, 3.4, "CO_2", "gray", 0.55, 0.36), A(1.2, 4.6, -0.2, 3.8, "gray"),
    T(3.0, 8.4, "1. decarboxylate", 0.55),
    T(12.8, 5.8, "2. acetyl → lipoamide", 0.55, None, "left"),
    I(9.0, 12.6, "CoA", "teal", 0.55, 0.36), A(9.4, 12.0, 10.4, 10.6), A(14.0, 10.6, 15.2, 12.0), I(15.8, 12.6, "Ac-CoA", "teal", 0.7, 0.32),
    T(12.0, 13.6, "3. acetyl-CoA", 0.55),
    T(21.0, 8.4, "4. reoxidise lipoamide", 0.55),
    I(19.6, 1.2, "NAD^+", "gray", 0.6, 0.34), A(20.2, 1.6, 21.0, 4.4), A(21.8, 4.4, 22.6, 1.6), I(23.2, 1.2, "NADH", "pink", 0.6, 0.34),
    T(21.0, 0.0, "5. FADH_2 → NADH", 0.55),
    # regulation
    T(27.0, 3.0, "Regulation", 0.75, None, "left"),
    B(28.4, 5.6, 2.6, 1.3, "blue", "PDH", ts=0.5, accent=True), B(33.0, 5.6, 2.6, 1.3, "gray", "PDH–P", ts=0.45),
    A(29.8, 5.0, 31.6, 5.0, "red", -0.3), T(30.7, 3.6, "PDK", 0.6, "red"),
    A(31.6, 6.2, 29.8, 6.2, "green", -0.3), T(30.7, 7.6, "PDP", 0.6, "green"),
    T(30.7, 9.0, "PDK ↑: ATP, NADH, acetyl-CoA", 0.55, "gray"), T(30.7, 9.8, "PDP ↑: Ca^{2+}, insulin", 0.55, "gray"),
]}

# --- RBC glutathione / NADPH cycle
GSH = {"items": [
    I(2.0, 2.0, "G6P", "orange", 0.75, 0.4), A(2.8, 2.0, 6.6, 2.0), I(7.6, 2.0, "6-PGL", "orange", 0.85, 0.34),
    T(4.7, 0.8, "G6PD", 0.6, "red"),
    I(3.4, 5.6, "NADP^+", "gray", 0.75, 0.32), I(6.2, 5.6, "NADPH", "pink", 0.75, 0.32),
    A(3.8, 4.8, 4.6, 2.8, "gray", 0.2), A(5.4, 2.8, 5.8, 4.8, "black", 0.2),
    # glutathione reductase
    A(6.2, 6.4, 7.4, 10.4), A(6.0, 10.4, 3.6, 6.4, "gray"),
    B(7.0, 11.2, 4.6, 1.2, "teal", "GSH reductase (FAD)", ts=0.42),
    I(9.6, 9.4, "GSSG", "yellow", 0.75, 0.36), I(9.6, 13.0, "2 GSH", "green", 0.8, 0.36),
    A(9.6, 10.2, 9.6, 12.2, "black"),
    # glutathione peroxidase
    B(14.4, 11.2, 4.6, 1.2, "purple", "GSH peroxidase (Se)", ts=0.42),
    A(10.4, 12.8, 12.6, 11.9, "black"), A(12.6, 10.5, 10.4, 9.6, "black"),
    I(18.6, 9.4, "H_2O_2", "red", 0.7, 0.36), A(18.0, 9.8, 16.6, 10.6), I(18.6, 13.0, "2 H_2O", "blue", 0.75, 0.34), A(16.6, 11.8, 18.0, 12.6),
    T(18.6, 7.6, "oxidant drugs, fava beans", 0.58, "gray"),
    T(9.0, 15.0, "G6PD deficiency → GSH ↓ → hemolysis", 0.6, "red"),
]}

# --- THF one-carbon forms (molecules)
THF_P = {
    "layout": "free", "cell_size": 15.0,
    "nodes": [
        {"id": "f", "ref": "formyl_thf", "pos": [0, 0], "decorate": {"highlight": [{"smarts": "[N:1][CH:2]=[O:3]", "color": M}]}},
        {"id": "m", "ref": "methenyl_thf", "pos": [1.75, 0], "decorate": {"highlight": [{"smarts": "[N+:1]=[CH:2][N:3]", "color": M}]}},
        {"id": "ch2", "ref": "methylene_thf", "pos": [0, 1.15], "decorate": {"highlight": [{"smarts": "[#7:1][CH2;r5:2][#7:3]", "color": M}]}},
        {"id": "ch3", "ref": "methyl_thf", "pos": [0, 2.3], "decorate": {"highlight": [{"smarts": "[CH3:1][N:2]", "color": M}]}},
        {"id": "pur", "text": "purine C2, C8\n(fMet-tRNA)", "kind": "concept", "pos": [0, -0.6]},
        {"id": "shmt", "text": "THF + Ser → Gly\n(SHMT)", "kind": "concept", "pos": [1.6, 0.92]},
        {"id": "dtmp", "text": "dUMP → dTMP\n(thymidylate synthase)", "kind": "concept", "pos": [1.6, 1.42]},
        {"id": "met", "text": "Hcy → Met, THF\n(methionine synthase, B_{12})", "kind": "concept", "pos": [1.6, 2.3]},
        {"id": "lv1", "text": "formate level", "kind": "heading", "pos": [2.95, 0]},
        {"id": "lv2", "text": "formaldehyde level", "kind": "heading", "pos": [2.95, 1.15]},
        {"id": "lv3", "text": "methanol level", "kind": "heading", "pos": [2.95, 2.3]},
    ],
    "edges": [
        {"from": "f", "to": "m", "label": "cyclohydrolase", "reversible": True},
        {"from": "m", "to": "ch2", "label": "MTHFD", "side_text": "NADPH", "reversible": True},
        {"from": "ch2", "to": "ch3", "label": "MTHFR (irreversible)", "side_text": "NADPH"},
        {"from": "shmt", "to": "ch2"}, {"from": "ch2", "to": "dtmp"}, {"from": "ch3", "to": "met"}, {"from": "f", "to": "pur"},
    ],
    "text": False,
}

# --- starvation plot (Lehninger 6e Fig. 23–31)
STARVE = {
    "size": [16, 10],
    "x": {"min": 0, "max": 40, "label": "days of starvation", "ticks": [0, 10, 20, 30, 40]},
    "y": {"min": 0, "max": 10, "label": "plasma concentration (mM)", "ticks": [0, 2, 4, 6, 8, 10]},
    "curves": [
        {"id": "bhb", "y": "6.1*(1-exp(-x/8))", "color": "red"},
        {"id": "glc", "y": "3.0 + 1.0*exp(-x/25) + 0.7*exp(-x/0.8)", "color": "cyan"},
        {"id": "ffa", "y": "0.5 + 0.8*(1-exp(-x/1.5)) + 0.35*x/40", "color": "orange"},
        {"id": "aca", "y": "1.35*(1-exp(-x/6))", "color": "magenta"},
    ],
    "legend": {"at": [2, 9.7], "items": [{"curve": "bhb", "text": "β-hydroxybutyrate"}, {"curve": "glc", "text": "glucose"},
                                           {"curve": "ffa", "text": "free fatty acids"}, {"curve": "aca", "text": "acetoacetate"}]},
}

PURIFY = table([
    ("Charge", "ion exchange (CM = cation, DEAE = anion) · isoelectric focusing · isoelectric precipitation"),
    ("Size", "gel filtration (large elute first) · dialysis · SDS-PAGE (mass) · ultracentrifugation (S)"),
    ("Solubility", "salting out ((NH_4)_2SO_4) · organic solvents · pH at pI"),
    ("Affinity", "ligand / antibody columns · His-tag (Ni^{2+}) · GST (glutathione)"),
    ("Hydrophobicity", "hydrophobic interaction · reverse-phase HPLC"),
    ("Check purity", "specific activity ↑ (units/mg) · one band on SDS-PAGE · 2-D gel (pI × mass)"),
], w_name=6.4, x_desc=7.0)

SEQUENCE = table([
    ("1. Purify", "> 98 % pure; count chains (end-group analysis)"),
    ("2. Composition", "6 M HCl, 110 °C, 24 h (Trp lost; Asn, Gln → Asp, Glu) · ninhydrin / HPLC"),
    ("3. Ends", "N-term: Sanger (FDNB), dansyl, Edman (PITC) · C-term: carboxypeptidase"),
    ("4. S–S", "performic acid, or DTT / β-ME + iodoacetate (blocks Cys)"),
    ("5. Cleave", "trypsin: after Lys, Arg · chymotrypsin: Phe, Trp, Tyr · CNBr: after Met"),
    ("6. Sequence", "Edman degradation of each fragment (≤ 50 residues)"),
    ("7. Overlap", "a second cleavage gives overlapping fragments → order"),
    ("8. Today", "MS (MALDI, ESI, MS/MS) · or deduce from the DNA / cDNA sequence"),
], w_name=6.0, x_desc=6.6, hue="teal")

# --- phospholipase sites on PC (ring view, R groups)
PC_R = "[25*]C(=O)OC[C@@H](OC(=O)[26*])COP(=O)([O-])OCC[N+](C)(C)C"
MOLECULES.append({"id": "phosphatidylcholine_r", "name": "Phosphatidylcholine (R groups)", "label": "Phosphatidylcholine",
                  "aliases": ["PC (generic)", "glycerophospholipid"], "smiles": PC_R, "ring": {"coordgen": True},
                  "tags": ["phospholipid"],
                  "notes": "以 R₁（sn-1）、R₂（sn-2）代表脂肪酸鏈的 PC，用來標示各種 phospholipase 的切點。"})
PLASE_P = {
    "layout": "free", "cell_size": 16.0,
    "nodes": [
        {"id": "pc", "ref": "phosphatidylcholine_r", "pos": [0, 0], "decorate": {"highlight": [
            {"smarts": "[25#0][C:1](=O)[O:2]", "color": "cyan"},
            {"smarts": "[26#0][C:1](=O)[O:2]", "color": "magenta"},
            {"smarts": "[CH2][O:1][P:2]", "color": "green"},
            {"smarts": "P[O:1][CH2:2]C[N+]", "color": "orange"}]}},
        {"id": "a1", "text": "PLA_1 (cyan): sn-1 ester → lyso-PC", "kind": "concept", "pos": [1.25, -0.45]},
        {"id": "a2", "text": "PLA_2 (magenta): sn-2 ester → arachidonate", "kind": "concept", "pos": [1.25, -0.15]},
        {"id": "c", "text": "PLC (green): → DAG + phosphocholine", "kind": "concept", "pos": [1.25, 0.15]},
        {"id": "d", "text": "PLD (orange): → phosphatidate + choline", "kind": "concept", "pos": [1.25, 0.45]},
    ],
    "edges": [], "text": False,
}

PATHWAYS = [
    {"id": "serine_protease_mechanism", "title": "Serine proteases: catalytic triad and S1 pocket", "layout": "scene", "scene": SERPROT,
     "tags": ["enzyme", "catalysis"],
     "notes": "催化三聯體（chymotrypsin 編號）：Asp102 以氫鍵（低障壁氫鍵）定位並穩定 His57 → His57 當廣義鹼拔取 Ser195 –OH 的質子 → "
              "Ser195 O⁻ 親核攻擊胜肽鍵的羰基碳 → 四面體中間物，O⁻ 嵌入氧陰離子洞（Gly193 與 Ser195 主鏈 NH 兩個氫鍵穩定過渡態）→ "
              "His57（此時為廣義酸）把質子給離去的胺 → 醯基-酵素；第二步由 His57 活化水進行去醯化。兩步反應見 chymotrypsin_acylation、"
              "chymotrypsin_deacylation。S1 口袋決定專一性：chymotrypsin 深而疏水（Phe、Trp、Tyr）；trypsin 底部 Asp189（Lys、Arg）；"
              "elastase 被 Val216、Thr226 填淺（Ala、Gly 等小側鏈）。共同點：都由 zymogen 活化（trypsinogen 由 enteropeptidase 切）；"
              "DFP 共價修飾 Ser195 而不可逆抑制。Subtilisin 為趨同演化的例子（不同摺疊、同樣三聯體）。"},
    {"id": "chaperonin_groel_cycle", "title": "GroEL / GroES chaperonin cycle", "layout": "scene", "scene": GROEL,
     "tags": ["protein", "protein_structure"],
     "notes": "GroEL（Hsp60）為兩個七聚體環背對背疊成的圓柱；GroES（Hsp10）為七聚體蓋子。① 未摺疊蛋白暴露的疏水區結合在開放環（cis）"
              "頂端的疏水內壁；② 7 ATP 與 GroES 結合 → 構形大幅改變，腔室變大、內壁變親水 → 蛋白質被推入孤立空間自行摺疊（約 10 秒）；"
              "③ ATP 水解；④ 另一端（trans）環結合 ATP 與新受質 → GroES 與摺好的蛋白質離開。尚未摺好可再進一次循環。"
              "真核胞質對應 TRiC（CCT），粒線體為 Hsp60/Hsp10。Hsp70（DnaK）則在轉譯中結合新生鏈的疏水片段。"},
    {"id": "pdh_complex_mechanism", "title": "Pyruvate dehydrogenase complex (E1–E2–E3)", "layout": "scene", "scene": PDH,
     "tags": ["TCA", "pyruvate"],
     "notes": "PDH 複合體三個酵素、五個輔因子（TPP、lipoamide、CoA、FAD、NAD⁺；對應 B1、B5 以外的硫辛酸、B5、B2、B3）。"
              "① E1（pyruvate dehydrogenase，TPP）脫羧放出 CO₂，形成 hydroxyethyl-TPP；② 羥乙基被氧化並轉到 E2 lipoyl-Lys 的硫上"
              "（acetyl-dihydrolipoamide）；③ E2（dihydrolipoyl transacetylase）把乙醯基轉給 CoA → acetyl-CoA；④ E3（dihydrolipoyl "
              "dehydrogenase，FAD）把還原態 lipoamide 再氧化；⑤ FADH₂ 把電子交給 NAD⁺ → NADH。Lipoyl-Lys 長臂在三個活性位之間擺動"
              "（substrate channeling）。α-KG dehydrogenase、branched-chain α-keto acid dehydrogenase 同樣結構（E3 共用）。"
              "調控：PDK 磷酸化 E1 → 失活（ATP、NADH、acetyl-CoA 活化 PDK；pyruvate、ADP 抑制 PDK）；PDP 去磷酸 → 活化"
              "（Ca²⁺、insulin）。砷（arsenite）結合 lipoamide；thiamine 缺乏（Wernicke、酒精）→ 乳酸堆積。"},
    {"id": "rbc_glutathione_cycle", "title": "NADPH and glutathione in the red cell", "layout": "scene", "scene": GSH,
     "tags": ["carb_metabolism", "PPP"],
     "notes": "紅血球沒有粒線體，NADPH 只能來自 PPP 的 G6PD。NADPH → glutathione reductase（FAD）把 GSSG 還原成 2 GSH → glutathione "
              "peroxidase（含 selenocysteine）用 GSH 把 H₂O₂ 還原成水。G6PD 缺乏（X 性聯、最常見酵素缺陷；抗瘧疾）：遇氧化壓力（primaquine、"
              "sulfonamide、dapsone、蠶豆的 divicine、感染）→ GSH 不足 → Hb 氧化沉澱（Heinz bodies）、膜受損 → 溶血（bite cells）。"
              "NADPH 其他用途：脂肪酸與固醇合成、P450、NOS、吞噬細胞 NADPH oxidase（呼吸爆發；缺陷 → 慢性肉芽腫病）。"},
    {"id": "thf_one_carbon_forms", "title": "One-carbon units on tetrahydrofolate", **THF_P,
     "tags": ["amino_acid_metabolism", "one_carbon", "vitamin", "coenzyme"],
     "notes": "洋紅為一碳單位，由上到下氧化程度遞減：formyl／methenyl（甲酸層級）→ methylene（甲醛層級）→ methyl（甲醇層級）。"
              "N¹⁰-formyl-THF：嘌呤 C2、C8，fMet-tRNA；N⁵,N¹⁰-methylene-THF：Ser → Gly（SHMT）產生，供 thymidylate synthase "
              "（dUMP → dTMP，同時給出還原力，THF → DHF，需 DHFR 再生；methotrexate、5-FU 作用點）；N⁵-methyl-THF：MTHFR 不可逆生成，"
              "唯一去路為 methionine synthase（B12）把甲基給 homocysteine。B12 缺乏 → methyl-folate trap。"},
    {"id": "phospholipase_sites", "title": "Where phospholipases cut", **PLASE_P, "tags": ["lipid", "phospholipid", "eicosanoid"],
     "notes": "PLA₁ 水解 sn-1 酯鍵；PLA₂ 水解 sn-2 酯鍵（釋出 arachidonate → eicosanoid；corticosteroid 經 lipocortin/annexin 抑制；"
              "蛇毒、胰臟 PLA₂）；PLC 切甘油 O–P 之間 → DAG + 磷酸頭基（PIP₂ → DAG + IP₃，Gq 訊號；細菌 α-toxin 為 PLC）；"
              "PLD 切 P–O（頭基一側）→ phosphatidate + choline。Lysophospholipid 有清潔劑性質（溶血）。"},
    {"id": "protein_purification", "title": "Protein purification: what each method separates by", "layout": "scene", "scene": PURIFY,
     "tags": ["protein", "protein_structure"],
     "notes": "依性質選方法：電荷（離子交換、等電聚焦、等電點沉澱）、大小（膠體過濾大分子先出、透析、SDS-PAGE 依質量、超高速離心）、"
              "溶解度（硫酸銨鹽析）、親和力（配體、抗體、His-tag–Ni²⁺、GST–glutathione）、疏水性（HIC、RP-HPLC）。"
              "純化表：總蛋白質下降、比活性（units/mg）上升，回收率下降。"},
    {"id": "protein_sequencing_strategy", "title": "Determining a protein's primary structure", "layout": "scene", "scene": SEQUENCE,
     "tags": ["protein", "peptide"],
     "notes": "Edman 降解：PITC（phenyl isothiocyanate）在弱鹼下接上 N 端胺基，酸處理把 N 端殘基以 PTH-胺基酸切下，其餘胜肽保留 → "
              "每次少一個，可重複約 50 次。Sanger 試劑 FDNB 標 N 端但會水解整條鏈（只能測一次）。不同切法產生重疊片段以排序；"
              "雙硫鍵位置：不先還原、切開後以 diagonal electrophoresis 比對。"},
    {"id": "starvation_fuel_plot", "title": "Plasma fuels during starvation", "layout": "plot", "plot": STARVE,
     "tags": ["energy_metabolism", "integration"],
     "notes": "依 Lehninger 6e 圖 23–31 的趨勢重畫（非精確數值）。血糖 2 天內即下降並維持約 3.5 mM；游離脂肪酸上升；酮體（以 "
              "β-hydroxybutyrate 為主）在 2–4 天後大幅上升到約 6–7 mM，腦改用酮體供應約 60–70 % 能量 → 糖質新生需求下降、肌肉蛋白分解減緩"
              "（尿素排泄由第 1–3 天約 30 g/天降到數週後約 10 g/天）。脂肪耗盡後才再分解蛋白質（含心肌、呼吸肌）→ 死亡。"},
]


# ================================================================ writer

def dump(path, d):
    if path.exists() and "status" in d:
        old = (yaml.safe_load(path.read_text(encoding="utf-8")) or {}).get("status")
        if old in ("approved", "changes_requested"):
            d["status"] = old
    path.write_text(yaml.safe_dump(d, allow_unicode=True, sort_keys=False, width=200), encoding="utf-8")


TREE = yaml.safe_load((Path(__file__).parent.parent / "data" / "tags.yaml").read_text(encoding="utf-8"))["tree"]
PARENT = {kid: top for top, kids in TREE.items() for kid in kids}


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
        lines = ["", "# --- 第二十六輪：酵素的 Ser195–OH（[22*] = 酵素其餘部分；由 tools/gen_gap.py 產生）"]
        lines += [f'{k}: {{label: "{NEW_COFACTORS[k]["label"]}", smiles: "{NEW_COFACTORS[k]["smiles"]}"}}' for k in new]
        p.write_text(text.rstrip("\n") + "\n" + "\n".join(lines) + "\n", encoding="utf-8")


def main(data):
    data = Path(data)
    add_cofactors(data)
    for m in MOLECULES:
        d = {"id": m["id"], "name": m["name"], "label": m["label"], "aliases": m["aliases"], "smiles": m["smiles"]}
        if m.get("ring"):
            d["ring"] = m["ring"]
        if m.get("highlight"):
            d["decorate"] = {"highlight": [{"smarts": m["highlight"], "color": M}]}
        d.update({"tags": with_parent(m["tags"]), "notes_zh": m["notes"], "parts_png": False, "status": "pending"})
        dump(data / "molecules" / f"{m['id']}.yaml", d)
    for r in REACTIONS:
        d = {k: r[k] for k in ("id", "title", "substrates", "products", "cofactors_in", "cofactors_out")}
        d.update({"enzyme": r["enzyme"], "reversible": False, "decorate": r["decorate"], "card_zh": r["card"],
                  "tags": with_parent(r["tags"]), "parts_png": False, "status": "pending"})
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
