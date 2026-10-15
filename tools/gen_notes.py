"""Round 15: gaps found by comparing with the Book 1 lecture notes (醣類／酵素／脂質／核酸／蛋白質).

    python tools/gen_notes.py data

Plots (mixed / substrate inhibition, burst kinetics, reaction coordinate, lysozyme pH profile, Hb T/R
states, His & Glu titration, DNA melting), small sugars (trehalose, gluconate, glucuronate, enediol),
2,3-BPG shunt, PTMs (Hyp, Hyl, Gla, P-Ser, acetyl-/trimethyl-Lys), cystine reduction, base pairs and
tautomeric mispairs (dashed H-bonds), deamination, and text overviews (PUFA, eicosanoid series,
isoprenoids, phospholipases, PIP2 signalling, enzyme classes, catalytic triad, Bohr effect, G protein).
"""
import math
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent))
from gen_misc import disaccharide  # noqa: E402
from reaction_colors import COLORS  # noqa: E402

M, C = "magenta", "cyan"
RED, GREEN, ORANGE, PURPLE, GRAY = "red", "green", "orange", "purple", "gray"


def n(id_, ref, x, y, **kw):
    return {"id": id_, "ref": ref, "pos": [x, y], "view": "name", **kw}


def t(id_, text, x, y, kind="concept", **kw):
    return {"id": id_, "text": text, "pos": [x, y], "kind": kind, **kw}


# ================================================================ plots

MM_LB = {"x": {"min": -1.0, "max": 1.4, "label": "1/[S]"}, "y": {"min": -0.6, "max": 6, "label": "1/v_0"}}


def mixed_panel(title, km, vm, note, xmin=-1.0, ymin=-0.6, lab=(0.25, [0.4, 0.2], "left")):
    return {
        "title": title, "size": [12, 8.5],
        "params": {"Vmax": 1, "Km": 2, "KmI": km, "VmI": vm},
        "x": {"min": xmin, "max": 1.4, "label": "1/[S]"}, "y": {"min": ymin, "max": 6, "label": "1/v_0"},
        "curves": [
            {"id": "no", "y": "Km/Vmax*x + 1/Vmax", "from": xmin, "to": 1.4, "label": "No inhibitor",
             "label_at": 1.15, "label_off": [0.3, 0.8]},
            {"id": "inh", "y": "KmI/VmI*x + 1/VmI", "from": xmin, "to": 1.4, "color": RED,
             "label": "+ Inhibitor", "label_at": lab[0], "label_off": lab[1], "label_align": lab[2]},
        ],
        "points": [{"on": "no", "x": 0, "y": "1/Vmax"}, {"on": "inh", "x": 0, "y": "1/VmI", "color": RED},
                   {"on": "no", "x": "-1/Km", "y": 0}, {"on": "inh", "x": "-1/KmI", "y": 0, "color": RED}],
        "texts": [{"at": [xmin + 0.05, 5.0], "text": note, "color": RED, "size": 0.9, "align": "left"}],
    }


INHIBITION_MIXED = [
    mixed_panel("Mixed, α > α' (common)", "1.8*Km", "Vmax/1.8*0.9", "K_m ↑  V_{max} ↓\nlines cross left of y axis,\nabove x axis"),
    mixed_panel("Mixed, α < α' (rare)", "Km/1.6", "Vmax/2.6", "K_m ↓  V_{max} ↓\nlines cross below the x axis",
                xmin=-1.7, ymin=-2.1, lab=(0.55, [0.35, 0.45], "left")),
]

SUBSTRATE_INHIBITION = [
    {"title": "Substrate inhibition", "size": [12, 8],
     "params": {"Vmax": 1, "Km": 1, "Ki": 12},
     "x": {"min": 0, "max": 30, "label": "[S]"}, "y": {"min": 0, "max": 1.05, "label": "v_0"},
     "curves": [{"id": "si", "y": "Vmax*x/(Km + x + x*x/Ki)", "width": 1.1},
                {"id": "mm", "y": "mm(x, Vmax, Km)", "color": GRAY, "dash": True}],
     "vlines": [{"x": "sqrt(Km*Ki)"}],
     "points": [{"on": "si", "x": "sqrt(Km*Ki)", "label": "optimum [S]", "label_off": [0.4, -0.55]}],
     "texts": [{"at": [20, 0.17], "text": "excess S binds a second\n(inhibitory) site", "color": RED, "size": 0.95},
               {"at": [21, 0.83], "text": "M–M (no inhibition)", "color": GRAY, "size": 0.9}]},
    {"title": "PFK-1 vs [ATP]", "size": [12, 8],
     "params": {"K": 0.4, "Ki": 3.2},
     "x": {"min": 0, "max": 6, "label": "[ATP]"}, "y": {"min": 0, "max": 1.05, "label": "PFK-1 activity"},
     "curves": [{"id": "pfk", "y": "1.25*x/(K + x) / (1 + (x/Ki)**4)", "width": 1.1, "color": M}],
     "texts": [{"at": [1.15, 0.3], "text": "ATP as substrate\n(active site)", "size": 0.85, "align": "left", "color": C},
               {"at": [4.6, 0.88], "text": "ATP as allosteric\ninhibitor (high energy)", "size": 0.9, "color": RED}]},
]

BURST = {
    "size": [14, 9],
    "params": {"A": 1, "k": 2.2, "s": 0.18},
    "x": {"min": 0, "max": 10, "label": "Time"}, "y": {"min": 0, "max": 3.2, "label": "[P_1] (p-nitrophenolate)"},
    "curves": [{"id": "p", "y": "A*(1 - exp(-k*x)) + s*x", "width": 1.1},
               {"id": "ss", "y": "A + s*x", "color": GRAY, "dash": True, "from": 0, "to": 10}],
    "points": [{"on": "ss", "x": 0, "y": "A", "color": GRAY}],
    "texts": [{"at": [0.3, 2.7], "text": "Burst phase:\nfast acylation\n(1 P_1 per enzyme)", "align": "left", "size": 0.9, "color": M},
              {"at": [7.0, 1.15], "text": "Steady state:\nslow deacylation\n(rate-limiting)", "size": 0.9, "color": C},
              {"at": [2.0, 0.55], "text": "y intercept = [E]_T", "size": 0.85, "color": GRAY, "align": "left"}],
}

REACTION_COORD = {
    "size": [15, 10],
    "params": {"G0": 0, "G1": -1.4, "Ea": 4.2, "Ec": 2.0},
    "x": {"min": 0, "max": 10, "label": "Reaction coordinate", "at": -2.6},
    "y": {"min": -2.6, "max": 5.2, "label": "Free energy, G", "label_pos": "center", "at": 0},
    "curves": [
        {"id": "un", "y": "G0 + (G1-G0)/(1+exp(-(x-5)*1.6)) + (Ea - (G1-G0)/2)*exp(-((x-5)/1.25)**2)", "width": 1.1},
        {"id": "cat", "y": "G0 + (G1-G0)/(1+exp(-(x-5)*1.6)) + (Ec - (G1-G0)/2)*exp(-((x-5)/1.25)**2)", "color": C},
    ],
    "hlines": [{"y": "G0", "to": 4.2}, {"y": "G1", "from": 5.8, "to": 10}],
    "arrows": [{"from": [3.3, "G0"], "to": [3.3, "Ea + (G1-G0)/2 + G0 - 0.05"], "color": RED,
                "label": "ΔG^‡ (uncatalyzed)", "label_off": [0.2, -2.9], "label_align": "right"},
               {"from": [4.2, "G0"], "to": [4.2, "Ec + (G1-G0)/2 + G0 - 0.05"], "color": C},
               {"from": [8.6, "G0"], "to": [8.6, "G1"], "color": M, "label": "ΔG (unchanged)",
                "label_off": [0.35, 0], "label_align": "left"}],
    "texts": [{"at": [0.6, -0.5], "text": "S", "size": 1.2}, {"at": [9.7, -1.95], "text": "P", "size": 1.2},
              {"at": [5, 4.55], "text": "Transition state (‡)", "size": 0.9},
              {"at": [0.3, 3.0], "text": "ΔG^‡ (enzyme)", "color": C, "align": "left"}],
}

LYSOZYME_PH = {
    "size": [14, 9],
    "params": {"pKa": 4.5, "pKb": 5.9},
    "x": {"min": 2, "max": 9, "label": "pH", "ticks": [2, 3, 4, 5, 6, 7, 8, 9]},
    "y": {"min": 0, "max": 1.15, "label": "Fraction / activity"},
    "curves": [
        {"id": "asp", "y": "1/(1+10**(pKa-x))", "color": C, "dash": True},
        {"id": "glu", "y": "1/(1+10**(x-pKb))", "color": M, "dash": True},
        {"id": "act", "y": "(1/(1+10**(pKa-x)))*(1/(1+10**(x-pKb)))/0.7", "width": 1.2},
    ],
    "points": [{"on": "asp", "x": "pKa", "y": 0.5, "color": C, "drop": "x"},
               {"on": "glu", "x": "pKb", "y": 0.5, "color": M, "drop": "x"}],
    "vlines": [{"x": "(pKa+pKb)/2", "label": "pH_{opt} ≈ (4.5 + 5.9)/2 = 5.2", "label_off": [3.2, 0]}],
    "legend": {"at": [6.3, 0.72], "items": [
        {"curve": "asp", "text": "Asp52 as COO^- (pK_a 4.5)"}, {"curve": "glu", "text": "Glu35 as COOH (pK_a 5.9)"},
        {"curve": "act", "text": "activity ~ product of both"}]},
}

HB_TR = {
    "size": [16, 10],
    "params": {"KR": 1.5, "KT": 250, "L": 8.5e4},
    "x": {"min": 0, "max": 110, "label": "pO_2 (mmHg)", "ticks": [0, 20, 40, 60, 80, 100]},
    "y": {"min": 0, "max": 1.12, "label": "Fractional saturation (Y)", "label_pos": "center", "ticks": [0, 0.5, 1]},
    "bands": [{"x0": 20, "x1": 40, "label": "Tissues", "color": C}, {"x0": 90, "x1": 105, "label": "Lungs", "color": M}],
    "curves": [
        {"id": "r", "y": "hill(x, KR, 1)", "color": GRAY, "dash": True},
        {"id": "t", "y": "hill(x, KT, 1)", "color": GRAY, "dash": True},
        {"id": "hb", "y": "mwc(x, KR, KT, L)", "color": RED, "width": 1.2},
    ],
    "arrows": [{"from": [24, 0.12], "to": [34, 0.62], "color": RED, "label": "T → R switch", "label_off": [1.6, 0.3]}],
    "legend": {"at": [54, 0.5], "items": [
        {"curve": "r", "text": "R state (strong binding, ≈ Mb)"}, {"curve": "t", "text": "T state (weak binding)"},
        {"curve": "hb", "text": "Hb: starts on T, ends on R → sigmoid"}]},
}


def titration(name, pks, pi_pair, labels, species, ymax=13.5):
    env = {f"pK{i + 1}": v for i, v in enumerate(pks)}
    env["pI"] = f"(pK{pi_pair[0]}+pK{pi_pair[1]})/2"
    xs = " + ".join(f"1/(1+10**(pK{i + 1}-t))" for i in range(len(pks)))
    pts = [{"on": "c", "x": i + 0.5, "y": f"pK{i + 1}", "drop": "y", "label": labels[i], "label_off": [0.4, 0.55]}
           for i in range(len(pks))]
    pts.append({"on": "c", "y": "pI", "color": M, "label": f"pI = {sum(pks[j - 1] for j in pi_pair) / 2:.2f}",
                "label_off": [-0.4, -0.5], "label_align": "right"})
    return {"title": name, "size": [12, 9.5], "params": env,
            "x": {"min": 0, "max": len(pks) + 0.15, "label": "OH^- (equivalents)", "ticks": list(range(len(pks) + 1))},
            "y": {"min": 0, "max": ymax, "label": "pH", "ticks": [0, 3, 6, 9, 12]},
            "curves": [{"id": "c", "t": [0.3, 13.2], "x": xs, "y": "t", "width": 1.1}],
            "points": pts,
            "texts": [{"at": [0.1, 12.6], "text": species, "align": "left", "size": 0.85, "color": GRAY}]}


TITRATION_HIS_GLU = [
    titration("Histidine", [1.82, 6.0, 9.17], (2, 3), ["pK_1 1.8 (α-COOH)", "pK_R 6.0 (imidazole)", "pK_2 9.2 (α-NH_3^+)"],
              "+2 → +1 → 0 → −1"),
    titration("Glutamate", [2.19, 4.25, 9.67], (1, 2), ["pK_1 2.2 (α-COOH)", "pK_R 4.3 (γ-COOH)", "pK_2 9.7 (α-NH_3^+)"],
              "+1 → 0 → −1 → −2"),
]

MELTING = {
    "size": [15, 9.5],
    "params": {"Tm1": 70, "Tm2": 82},
    "x": {"min": 40, "max": 100, "label": "Temperature (°C)", "ticks": [40, 60, 80, 100], "at": 0.7},
    "y": {"min": 0, "max": 1.5, "label": "A_{260}", "ticks": [], "at": 40},
    "curves": [
        {"id": "at", "y": "1 + 0.37/(1+exp(-(x-Tm1)/1.4))", "width": 1.1, "color": C},
        {"id": "gc", "y": "1 + 0.37/(1+exp(-(x-Tm2)/1.4))", "width": 1.1, "color": M},
        {"id": "ss", "y": "1.0 + 0.0042*(x-40)", "color": ORANGE, "dash": True},
    ],
    "points": [{"on": "at", "x": "Tm1", "y": 1.185, "color": C, "drop": "x"},
               {"on": "gc", "x": "Tm2", "y": 1.185, "color": M, "drop": "x"}],
    "x_marks": True,
    "legend": {"at": [62, 0.93], "items": [
        {"curve": "at", "text": "AT-rich dsDNA (T_m lower)"}, {"curve": "gc", "text": "GC-rich dsDNA (T_m higher; G≡C)"},
        {"curve": "ss", "text": "ssRNA: gradual, non-cooperative"}]},
    "texts": [{"at": [92, 1.44], "text": "single strands\n(hyperchromic)", "size": 0.85, "color": GRAY},
              {"at": [41, 1.3], "text": "double helix\n(stacked bases\nabsorb less)", "size": 0.85, "color": GRAY,
               "align": "left"}],
}
MELTING.pop("x_marks")
MELTING["x"]["marks"] = [{"at": "Tm1", "label": "T_m", "color": C}, {"at": "Tm2", "label": "T_m", "color": M}]
MELTING["y"]["min"] = 0.7


# ================================================================ sugars

GLC_A = "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O"
ALPHA = {1: {"up": "H", "down": "link"}, 2: {"up": "H", "down": "OH"}, 3: {"up": "OH", "down": "H"},
         4: {"up": "H", "down": "HO"}, 5: {"up": "CH_2OH", "down": "H"}}
# right-hand glucose drawn mirrored (anomeric C on the left); up/down as checked by the validator.
# CH2OH sits inside the ring, so it is written HOCH_2 with the bond on its C (text runs left).
ALPHA_M = {1: {"up": "link", "down": "H"}, 2: {"up": "OH", "down": "H"}, 3: {"up": "H", "down": "OH"},
           4: {"up": "OH", "down": "H"}, 5: {"up": "H", "down": ["HOCH_2", 2], "down_len": 0.7}}


def fischer(top, rows, bottom):
    return [top] + [{"c": "C", "l": l, "r": r} for l, r in rows] + [bottom]


GLC_ROWS = [("H", "OH"), ("HO", "H"), ("H", "OH"), ("H", "OH")]

MOLECULES = [
    {"id": "trehalose", "name": "Trehalose", "aliases": ["海藻糖", "α,α-trehalose"],
     "smiles": disaccharide(GLC_A, GLC_A, 1),
     "haworth": {"units": [{"ring": "pyranose", "subs": ALPHA},
                           {"ring": "pyranose", "subs": ALPHA_M, "at": [4.7, 0.9], "mirror": True}],
                 "links": [{"from": [0, 1, "down"], "to": [1, 1, "up"], "o": [2.25, 0.75]}]},
     "tags": ["carbohydrate", "disaccharide"],
     "notes": "Glc(α1↔α1)Glc：兩個 anomeric C 互相連結 → 沒有游離 anomeric OH → 非還原糖（同 sucrose）。昆蟲血糖、真菌、抗乾燥。右側葡萄糖左右翻（mirror）畫，anomeric C 在左，以 α 鍵與左側相接。"},
    {"id": "d_gluconate", "name": "D-Gluconate", "aliases": ["gluconic acid", "葡萄糖酸"],
     "smiles": "[O-]C(=O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO",
     "chain": fischer({"c": "COO^-"}, GLC_ROWS, {"c": "CH_2OH"}),
     "tags": ["carbohydrate", "sugar_derivative"],
     "notes": "C1（醛基）氧化成羧酸：aldonic acid。glucose oxidase 的產物（血糖試紙）；PPP 的 6-phosphogluconate 即其 6-磷酸。"},
    {"id": "d_glucuronate_open", "name": "D-Glucuronate (open chain)", "label": "D-Glucuronate",
     "aliases": ["glucuronic acid (Fischer)", "葡萄糖醛酸（開鏈）"],
     "smiles": "O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)C(=O)[O-]",
     "chain": fischer({"c": "C", "fork": {"l": "O", "lb": 2, "r": "H"}}, GLC_ROWS, {"c": "COO^-"}),
     "tags": ["carbohydrate", "sugar_derivative"],
     "notes": "C6（一級醇）氧化成羧酸：uronic acid，C1 醛基保留。UDP-glucuronate → 膽紅素、藥物 glucuronidation；GAG 的組成；植物與多數動物可由它合成 L-ascorbate（人類缺 gulonolactone oxidase）。"},
    {"id": "enediol", "name": "1,2-Enediol", "label": "cis-1,2-Enediol",
     "aliases": ["enediol intermediate", "烯二醇"],
     "smiles": "O/C=C(/O)[C@@H](O)[C@H](O)[C@H](O)CO",
     "chain": [{"c": "C", "l": "H", "r": "OH", "down": 2}, {"c": "C", "r": "OH"}] +
              [{"c": "C", "l": l, "r": r} for l, r in GLC_ROWS[1:]] + [{"c": "CH_2OH"}],
     "tags": ["carbohydrate", "monosaccharide"],
     "notes": "鹼性條件下（或 phosphohexose isomerase）aldose ⇌ ketose 的共同中間物：C1=C2 雙鍵、兩個 OH 在同側（cis）。"
              "由 enediol 重新質子化可得 glucose、mannose（C2 epimer）或 fructose（Lobry de Bruyn–van Ekenstein 轉換）。"},
    {"id": "bpg_23", "name": "2,3-Bisphosphoglycerate", "label": "2,3-Bisphospho-\nglycerate",
     "aliases": ["2,3-BPG", "2,3-DPG", "2,3-二磷酸甘油酸"],
     "smiles": "[O-]C(=O)[C@H](OP(=O)([O-])[O-])COP(=O)([O-])[O-]",
     "chain": [{"c": "C", "fork": {"l": "O", "lb": 2, "r": "O^-"}}, {"c": "C", "l": "H", "r": "O{P}"},
               {"c": "CH_2O{P}"}],
     "tags": ["carb_metabolism", "glycolysis", "hemoglobin"],
     "notes": "RBC 內約 5 mM，由 1,3-BPG 經 bisphosphoglycerate mutase 形成（Rapoport–Luebering shunt，少產 1 ATP）。"
              "結合 deoxy-Hb 中央空腔（β 鏈正電殘基）→ 穩定 T 態、O₂ 曲線右移。"},
    {"id": "cystine", "name": "Cystine", "aliases": ["Cys-S-S-Cys", "胱胺酸", "disulfide"],
     "smiles": "[NH3+][C@@H](CSSC[C@H]([NH3+])C(=O)[O-])C(=O)[O-]",
     "free": {"atoms": {"n1": ["^+H_3N", -3.95, 0.55], "ca1": ["CH", -2.75, 0], "a1": ["COO^-", -2.75, -1.1],
                        "cb1": ["CH_2", -1.6, 0.55], "s1": ["S", -0.55, 0], "s2": ["S", 0.55, 0.55],
                        "cb2": ["CH_2", 1.6, 0], "ca2": ["CH", 2.75, 0.55], "a2": ["COO^-", 2.75, 1.65],
                        "n2": ["NH_3^+", 3.95, 0]},
              "bonds": [["n1", "ca1"], ["ca1", "a1"], ["ca1", "cb1"], ["cb1", "s1"], ["s1", "s2"], ["s2", "cb2"],
                        ["cb2", "ca2"], ["ca2", "a2"], ["ca2", "n2"]],
              "notes": [["disulfide bond", 0, -0.75, "orange"]]},
     "decorate": {"colors": {"s1": "orange", "s2": "orange"}},
     "tags": ["amino_acid", "protein_structure"],
     "notes": "兩個 Cys 的 –SH 氧化成 –S–S–（disulfide bond，共價，穩定三級／四級結構，多見於分泌蛋白）。"
              "β-mercaptoethanol、DTT 還原打斷（可逆）；performic acid 氧化成 cysteic acid（不可逆）；iodoacetate 封住 –SH 防止重接。"},
    {"id": "selenocysteine", "name": "Selenocysteine", "aliases": ["Sec", "U", "硒半胱胺酸", "21st amino acid"],
     "smiles": "[SeH]C[C@H]([NH3+])C(=O)[O-]",
     "chain": [{"c": "COO^-"}, {"c": "C", "l": "^+H_3N", "r": "H"}, {"c": "CH_2"}, {"c": "SeH"}],
     "decorate": {"boxes": [{"rows": [2, 3], "name": "R"}]},
     "tags": ["amino_acid", "polar"],
     "notes": "第 21 個胺基酸（Sec, U）：Cys 的 S 換成 Se，pK_a ≈ 5.2（pH 7 帶負電，更強的親核基）。由 UGA（原為終止密碼）配合 SECIS 元件重新編碼。"
              "存在於 glutathione peroxidase、thioredoxin reductase、deiodinase（甲狀腺素）。第 22 個為 pyrrolysine（UAG，僅古菌）。"},
    {"id": "hydroxyproline", "name": "4-Hydroxyproline", "aliases": ["Hyp", "trans-4-hydroxy-L-proline", "羥脯胺酸"],
     "smiles": "[O-]C(=O)[C@@H]1C[C@@H](O)C[NH2+]1",
     "free": {"atoms": {"a": ["COO^-", 0, -0.15], "ca": ["C", 0, 1], "h": ["H", 0.72, 0.5], "n": ["^+H_2N", -1.15, 1],
                        "cb": ["CH_2", 1.15, 1], "cg": ["CH", 1.15, 2], "oh": ["OH", 2.3, 2], "cd": ["H_2C", -1.15, 2]},
              "bonds": [["a", "ca"], ["ca", "h"], ["ca", "n"], ["ca", "cb"], ["cb", "cg"], ["cg", "oh"], ["cg", "cd"],
                        ["cd", "n"]]},
     "decorate": {"colors": {"oh": "magenta"}},
     "tags": ["amino_acid", "protein_structure"],
     "notes": "Collagen 的 Gly-X-Y 中 Y 位常為 Hyp：prolyl 4-hydroxylase 在 ER 內修飾（轉譯後），需 Fe²⁺、α-KG、O₂ 與 vitamin C（維持 Fe²⁺）。"
              "Hyp 的 OH 形成氫鍵穩定三股螺旋；缺 vitamin C → scurvy（牙齦出血、傷口難癒）。HIF-1α 的 Pro 羥化則標記它被降解。"},
    {"id": "hydroxylysine", "name": "5-Hydroxylysine", "aliases": ["Hyl", "羥離胺酸"],
     "smiles": "[NH3+]C[C@H](O)CC[C@H]([NH3+])C(=O)[O-]",
     "chain": [{"c": "COO^-"}, {"c": "C", "l": "^+H_3N", "r": "H"}, {"c": "CH_2"}, {"c": "CH_2"},
               {"c": "C", "l": "HO", "r": "H"}, {"c": "CH_2"}, {"c": "^+NH_3"}],
     "decorate": {"colors": {4: {"l": "magenta"}}},
     "tags": ["amino_acid", "protein_structure"],
     "notes": "Collagen 中 lysyl hydroxylase（同樣需 vitamin C、Fe²⁺、α-KG）產物；Hyl 的 OH 可接 Gal／Glc（O-linked glycosylation），"
              "且 lysyl oxidase（含 Cu²⁺）把 Lys／Hyl 氧化成醛 → 交聯。Ehlers–Danlos（部分型）為 lysyl hydroxylase 缺陷。"},
    {"id": "gla", "name": "γ-Carboxyglutamate", "aliases": ["Gla", "γ-羧基麩胺酸"],
     "smiles": "[NH3+][C@@H](CC(C(=O)[O-])C(=O)[O-])C(=O)[O-]",
     "chain": [{"c": "COO^-"}, {"c": "C", "l": "^+H_3N", "r": "H"}, {"c": "CH_2"},
               {"c": "C", "l": "H", "r": "COO^-"}, {"c": "COO^-"}],
     "decorate": {"colors": {3: {"r": "magenta"}}},
     "tags": ["amino_acid", "protein_structure"],
     "notes": "凝血因子 II、VII、IX、X（與 protein C、S）N 端的 Glu 由 γ-glutamyl carboxylase 加上第二個 COO⁻（需 vitamin K hydroquinone）。"
              "兩個相鄰 COO⁻ 夾住 Ca²⁺，使因子能結合磷脂膜。Warfarin 抑制 VKOR（vitamin K epoxide reductase）→ 無法再生 KH₂。"},
    {"id": "phosphoserine", "name": "Phosphoserine", "aliases": ["P-Ser", "pSer", "磷酸絲胺酸"],
     "smiles": "[O-]P(=O)([O-])OC[C@H]([NH3+])C(=O)[O-]",
     "chain": [{"c": "COO^-"}, {"c": "C", "l": "^+H_3N", "r": "H"}, {"c": "CH_2O{P}"}],
     "decorate": {"colors": {2: [M, "{P}"]}},
     "tags": ["amino_acid", "protein_structure"],
     "notes": "Ser／Thr kinase（如 PKA、PKC）把 ATP 的 γ-磷酸接到 –OH；Tyr kinase（如 insulin receptor）作用在 Tyr。phosphatase 水解移除。最常見的可逆調控型 PTM。"},
    {"id": "acetyllysine", "name": "N^ε-Acetyllysine", "label": "N^ε-Acetyllysine",
     "aliases": ["acetyl-Lys", "Kac", "乙醯離胺酸"],
     "smiles": "CC(=O)NCCCC[C@H]([NH3+])C(=O)[O-]",
     "chain": [{"c": "COO^-"}, {"c": "C", "l": "^+H_3N", "r": "H"}, {"c": "CH_2"}, {"c": "CH_2"}, {"c": "CH_2"},
               {"c": "CH_2"}, {"c": "NH"}, {"c": "C", "r": "O", "rb": 2}, {"c": "CH_3"}],
     "decorate": {"colors": {7: M, 8: M}},
     "tags": ["amino_acid", "protein_structure"],
     "notes": "組蛋白 Lys 的 ε-NH₃⁺ 被 HAT（histone acetyltransferase）乙醯化 → 失去正電、與 DNA 磷酸的吸引減弱 → 染色質變鬆、基因表現（gene on）。HDAC 移除。"},
    {"id": "trimethyllysine", "name": "N^ε-Trimethyllysine", "label": "N^ε-Trimethyllysine",
     "aliases": ["Kme3", "三甲基離胺酸"],
     "smiles": "C[N+](C)(C)CCCC[C@H]([NH3+])C(=O)[O-]",
     "chain": [{"c": "COO^-"}, {"c": "C", "l": "^+H_3N", "r": "H"}, {"c": "CH_2"}, {"c": "CH_2"}, {"c": "CH_2"},
               {"c": "CH_2"}, {"c": "N^+", "l": "H_3C", "r": "CH_3"}, {"c": "CH_3"}],
     "decorate": {"colors": {6: {"l": M, "r": M}, 7: M}},
     "tags": ["amino_acid", "protein_structure"],
     "notes": "HMT（histone methyltransferase）以 SAM 為甲基供體，可做 mono-/di-/tri-methyl。甲基化不改變電荷；依位置可活化（H3K4me3）或抑制（H3K9me3、H3K27me3）轉錄。"},
    {"id": "methylcytosine", "name": "5-Methylcytosine", "aliases": ["5mC", "5-甲基胞嘧啶"],
     "smiles": "Cc1c[nH]c(=O)nc1N",
     "ring": {"down": 4, "right": 5, "nh_bond": True, "mirror": True},
     "tags": ["nucleotide", "pyrimidine"],
     "notes": "CpG 中 C 的 5 位甲基化（DNMT，SAM）→ 基因沉默（epigenetic）。自發去胺基直接變成 thymine，修復系統不易辨認 → CpG 為突變熱點。"},
]


# ---------------------------------------------------------------- base pairs (free drawings, dashed H-bonds)

B = 0.95          # bond length (pitches)
R6 = B            # hexagon radius


def pt(c, ang, r):
    a = math.radians(ang)
    return (round(c[0] + r * math.cos(a), 3), round(c[1] - r * math.sin(a), 3))


def purine(cx, kind):
    """atoms + bonds of a purine on the left; N1 faces right. kind: A, Aimino, G."""
    c = (cx, 0)
    ang = {"n1": 30, "c6": 90, "c5": 150, "c4": 210, "n3": 270, "c2": 330}
    P = {k: pt(c, a, R6) for k, a in ang.items()}
    # imidazole fused on C4–C5 (left edge): regular pentagon outward
    mid = ((P["c4"][0] + P["c5"][0]) / 2, (P["c4"][1] + P["c5"][1]) / 2)
    apo = B / (2 * math.tan(math.radians(36)))
    pc = (mid[0] - apo, mid[1])
    rp = B / (2 * math.sin(math.radians(36)))
    P["n7"] = pt(pc, 108, rp)
    P["c8"] = pt(pc, 180, rp)
    P["n9"] = pt(pc, 252, rp)
    atoms = {"n1": ["N", *P["n1"]], "c6": ["", *P["c6"]], "c5": ["", *P["c5"]], "c4": ["", *P["c4"]],
             "n3": ["N", *P["n3"]], "c2": ["", *P["c2"]], "n7": ["N", *P["n7"]], "c8": ["", *P["c8"]],
             "n9": ["N", *P["n9"]]}
    bonds = [["c6", "n1", 2 if kind == "A" else 1], ["n1", "c2"], ["c2", "n3", 2], ["n3", "c4"], ["c4", "c5", 2],
             ["c5", "c6"], ["c5", "n7"], ["n7", "c8", 2], ["c8", "n9"], ["n9", "c4"]]
    hb_from = {}
    up6 = pt(P["c6"], 90, B)
    if kind in ("A", "Aimino"):
        atoms["n6"] = ["N", *up6]
        if kind == "A":
            atoms["h6a"] = ["H", *pt(up6, 30, 0.8)]
            atoms["h6b"] = ["H", *pt(up6, 150, 0.8)]
            bonds += [["c6", "n6"], ["n6", "h6a"], ["n6", "h6b"]]
            hb_from["top"] = "h6a"
            hb_from["mid"] = "n1"
        else:                                           # imino tautomer: C6=N6–H (H away), N1–H
            atoms["h6b"] = ["H", *pt(up6, 150, 0.8)]
            atoms["h1"] = ["H", *pt(P["n1"], 30, 0.8)]
            bonds += [["c6", "n6", 2], ["n6", "h6b"], ["n1", "h1"]]
            hb_from["top"] = "n6"
            hb_from["mid"] = "h1"
    else:                                               # G: C6=O6, N1–H, C2–N2H2
        atoms["o6"] = ["O", *up6]
        atoms["h1"] = ["H", *pt(P["n1"], 30, 0.8)]
        n2 = pt(P["c2"], 330, B)
        atoms["n2"] = ["N", *n2]
        atoms["h2a"] = ["H", *pt(n2, 10, 0.75)]
        atoms["h2b"] = ["H", *pt(n2, 270, 0.8)]
        bonds += [["c6", "o6", 2], ["n1", "h1"], ["c2", "n2"], ["n2", "h2a"], ["n2", "h2b"]]
        hb_from["top"] = "o6"
        hb_from["mid"] = "h1"
        hb_from["low"] = "h2a"
    sugar = pt(P["n9"], 230, B * 1.1)
    atoms["s1"] = ["Sugar", *sugar, 0]
    bonds.append(["n9", "s1"])
    return atoms, bonds, hb_from, P


def pyrimidine(cx, kind):
    """atoms + bonds of a pyrimidine on the right; N3 faces left. kind: T, Tenol, C."""
    c = (cx, 0)
    ang = {"n3": 150, "c4": 90, "c5": 30, "c6": 330, "n1": 270, "c2": 210}
    P = {k: pt(c, a, R6) for k, a in ang.items()}
    atoms = {k: [("N" if k[0] == "n" else ""), *v] for k, v in P.items()}
    up4 = pt(P["c4"], 90, B)
    o2 = pt(P["c2"], 210, B)
    atoms["o2"] = ["O", *o2]
    bonds = [["n1", "c2"], ["c2", "o2", 2], ["c5", "c6", 2], ["c6", "n1"]]
    hb_to = {}
    if kind == "T":
        atoms["o4"] = ["O", *up4]
        atoms["h3"] = ["H", *pt(P["n3"], 150, 0.8)]
        atoms["me"] = ["CH_3", *pt(P["c5"], 30, B), 0]
        bonds += [["c2", "n3"], ["n3", "c4"], ["c4", "c5"], ["c4", "o4", 2], ["n3", "h3"], ["c5", "me"]]
        hb_to.update(top="o4", mid="h3")
    elif kind == "Tenol":                               # enol: C4–O4–H, N3=C4
        atoms["o4"] = ["O", *up4]
        atoms["h4"] = ["H", *pt(up4, 150, 0.8)]
        atoms["me"] = ["CH_3", *pt(P["c5"], 30, B), 0]
        bonds += [["c2", "n3"], ["n3", "c4", 2], ["c4", "c5"], ["c4", "o4"], ["o4", "h4"], ["c5", "me"]]
        hb_to.update(top="h4", mid="n3", low="o2")
    else:                                               # C: C4–N4H2, N3=C4 … wait: C has C4=N3
        atoms["n4"] = ["N", *up4]
        atoms["h4a"] = ["H", *pt(up4, 150, 0.8)]
        atoms["h4b"] = ["H", *pt(up4, 30, 0.8)]
        bonds += [["c2", "n3"], ["n3", "c4", 2], ["c4", "c5"], ["c4", "n4"], ["n4", "h4a"], ["n4", "h4b"]]
        hb_to.update(top="h4a", mid="n3", low="o2")
    sugar = pt(P["n1"], 300, B * 1.1)
    atoms["s2"] = ["Sugar", *sugar, 0]
    bonds.append(["n1", "s2"])
    return atoms, bonds, hb_to


def rotated(atoms, c, deg):
    """Rotate every atom position about c by deg (counter-clockwise on screen)."""
    a = math.radians(deg)
    out = {}
    for k, (txt, x, y, *rest) in atoms.items():
        dx, dy = x - c[0], -(y - c[1])
        rx, ry = dx * math.cos(a) - dy * math.sin(a), dx * math.sin(a) + dy * math.cos(a)
        out[k] = [txt, round(c[0] + rx, 3), round(c[1] - ry, 3), *rest]
    return out


def base_pair(left, right, pairs, notes, cx=2.0, rot=15):
    la, lb, lhb, _ = purine(-cx, left)
    ra, rb, rhb = pyrimidine(cx, right)
    la, ra = rotated(la, (-cx, 0), -rot), rotated(ra, (cx, 0), rot)
    atoms = {f"l_{k}": v for k, v in la.items()}
    atoms.update({f"r_{k}": v for k, v in ra.items()})
    bonds = [[f"l_{b[0]}", f"l_{b[1]}", *b[2:]] for b in lb] + [[f"r_{b[0]}", f"r_{b[1]}", *b[2:]] for b in rb]
    for key in pairs:
        bonds.append([f"l_{lhb[key]}", f"r_{rhb[key]}", "hb"])
    return {"schematic": True, "atoms": atoms, "bonds": bonds, "notes": notes}


GROOVES = [["Major groove", 0, -2.85, "gray"], ["Minor groove", 0, 2.75, "gray"]]
BASE_PAIRS = [
    {"id": "base_pair_at", "name": "A=T base pair", "label": "A=T (2 H-bonds)", "aliases": ["A-T pair", "Watson-Crick AT"],
     "smiles": "Nc1ncnc2[nH]cnc12.Cc1c[nH]c(=O)[nH]c1=O", "free": base_pair("A", "T", ["top", "mid"], GROOVES),
     "notes": "Watson–Crick A=T：2 個氫鍵（A 的 N6–H⋯O4，A 的 N1⋯H–N3）。糖接在同一側 → 上方大溝、下方小溝。A=U 在 RNA 中相同。"},
    {"id": "base_pair_gc", "name": "G≡C base pair", "label": "G≡C (3 H-bonds)", "aliases": ["G-C pair", "Watson-Crick GC"],
     "smiles": "Nc1nc2[nH]cnc2c(=O)[nH]1.Nc1cc[nH]c(=O)n1", "free": base_pair("G", "C", ["top", "mid", "low"], GROOVES),
     "notes": "Watson–Crick G≡C：3 個氫鍵（G O6⋯H–N4 C、G N1–H⋯N3 C、G N2–H⋯O2 C）→ GC 含量越高 T_m 越高。"},
    {"id": "mispair_a_imino_c", "name": "A* (imino) : C mispair", "label": "A* (imino) pairs with C",
     "aliases": ["tautomeric mispair", "imino adenine"],
     "smiles": "N=c1nc[nH]c2[nH]cnc12.Nc1cc[nH]c(=O)n1", "free": base_pair("Aimino", "C", ["top", "mid"],
                                                                         [["rare imino tautomer (~1/10^4)", 0, -2.85, "red"]]),
     "notes": "互變異構：A 的 amino 型罕見地變成 imino 型（A*，像 keto）→ 氫鍵供體/受體位置對調，改與 C 配對 → 複製後 A:T → G:C 轉換突變（transition）。"},
    {"id": "mispair_t_enol_g", "name": "T* (enol) : G mispair", "label": "T* (enol) pairs with G",
     "aliases": ["enol thymine", "tautomeric mispair"],
     "smiles": "Nc1nc2[nH]cnc2c(=O)[nH]1.Cc1c[nH]c(=O)nc1O", "free": base_pair("G", "Tenol", ["top", "mid", "low"],
                                                                             [["rare enol tautomer (~1/10^4)", 0, -2.85, "red"]]),
     "notes": "T 的 keto 型罕見地變成 enol 型（T*，像 amino）→ 與 G 形成 3 個氫鍵 → 複製後 T:A → C:G 轉換突變。"},
]


# ---------------------------------------------------------------- schematic overview molecules

PHOSPHOLIPASE = {
    "id": "phospholipase_sites", "name": "Phospholipase cleavage sites", "label": "Phospholipase sites",
    "aliases": ["PLA1", "PLA2", "PLC", "PLD", "磷脂酶"],
    "smiles": "CCC(=O)OCC(COP(=O)([O-])OC)OC(=O)CC",
    "free": {"schematic": True,
             "atoms": {"r1": ["R_1", -3.4, -1.0], "c1o": ["C", -2.3, -1.0], "o1d": ["O", -2.3, -2.05], "o1": ["O", -1.2, -1.0],
                       "g1": ["CH_2", 0, -1.0], "g2": ["CH", 0, 0.05], "o2": ["O", -1.2, 0.05], "c2o": ["C", -2.3, 0.05],
                       "o2d": ["O", -2.3, 1.1], "r2": ["R_2", -3.4, 0.05], "g3": ["CH_2", 0, 1.1], "o3": ["O", 1.15, 1.1],
                       "p": ["P", 2.3, 1.1], "pd": ["O", 2.3, 0.05], "pm": ["O^-", 2.3, 2.15], "o4": ["O", 3.45, 1.1],
                       "x": ["X", 4.5, 1.1]},
             "bonds": [["r1", "c1o"], ["c1o", "o1d", 2], ["c1o", "o1"], ["o1", "g1"], ["g1", "g2"], ["g2", "o2"],
                       ["o2", "c2o"], ["c2o", "o2d", 2], ["c2o", "r2"], ["g2", "g3"], ["g3", "o3"], ["o3", "p"],
                       ["p", "pd", 2], ["p", "pm"], ["p", "o4"], ["o4", "x"]],
             "notes": [["PLA_1", -1.55, -1.5, "red"], ["PLA_2", -1.55, 0.6, "red"], ["PLC", 1.65, 1.6, "red"],
                       ["PLD", 2.95, 0.6, "red"], ["sn-1", 1.3, -1.0, "gray"], ["sn-2", 1.3, 0.05, "gray"],
                       ["sn-3", -0.95, 1.55, "gray"], ["X = head group", 4.6, 1.75, "gray"]]},
    "tags": ["lipid", "phospholipid"],
    "notes": "磷脂酶切點：PLA₁ 切 sn-1 酯鍵；PLA₂ 切 sn-2 酯鍵（釋出 arachidonate → 類二十烷酸；glucocorticoid 經 lipocortin 抑制）；"
             "PLC 切甘油側的磷酸酯（PIP₂ → DAG + IP₃）；PLD 切頭基側（→ phosphatidate + 頭基）。PLB 兼具 A₁、A₂ 活性。",
}

TRIAD = {
    "id": "catalytic_triad", "name": "Catalytic triad (serine protease)", "label": "Catalytic triad",
    "aliases": ["Ser-His-Asp", "chymotrypsin active site", "催化三元組"],
    "smiles": "CC(=O)[O-].Cc1c[nH]cn1.CCO",
    "free": {"schematic": True,
             "atoms": {"acb": ["CH_2", -5.2, 0.55], "acg": ["C", -4.1, 0], "ao1": ["O", -4.1, -1.05], "ao2": ["O^-", -3.0, 0.55],
                       "hc": ["CH_2", -0.55, 2.1], "hcg": ["", -0.55, 1.15], "hnd": ["N", -1.45, 0.6], "hd": ["H", -2.15, 0.6],
                       "hce": ["", -1.0, -0.35], "hne": ["N", 0.05, -0.35], "hcd": ["", 0.35, 0.65],
                       "sh": ["H", 1.05, -0.35], "so": ["O", 1.95, -0.35], "scb": ["CH_2", 2.9, 0.2]},
             "bonds": [["acb", "acg"], ["acg", "ao1", 2], ["acg", "ao2"], ["ao2", "hd", "hb"], ["hd", "hnd"],
                       ["hc", "hcg"], ["hcg", "hnd"], ["hnd", "hce", 2], ["hce", "hne"], ["hne", "hcd"], ["hcd", "hcg", 2],
                       ["hne", "sh", "hb"], ["sh", "so"], ["so", "scb"]],
             "notes": [["Asp102", -4.4, 1.5, "red"], ["His57", -0.45, 2.95, "cyan"], ["Ser195", 2.9, 1.15, "magenta"],
                       ["general base", -0.1, -1.3, "cyan"], ["nucleophile", 2.95, -1.3, "magenta"],
                       ["orients His", -4.1, -1.9, "red"]]},
    "decorate": {"colors": {"so": "magenta", "sh": "magenta", "hne": "cyan"}},
    "tags": ["enzyme", "catalysis"],
    "notes": "Chymotrypsin／trypsin／elastase 的催化三元組：Asp102 的 COO⁻ 固定 His57 方向並穩定其正電；His57 當 general base 拉走 Ser195 的 H⁺；"
             "Ser195 的 O⁻ 親核攻擊受質羰基 C → 四面體中間物（oxyanion hole：Gly193、Ser195 主鏈 NH 穩定 O⁻）→ acyl-enzyme（burst）→ 水解（慢）。"
             "DIFP 共價修飾 Ser195、TPCK 修飾 His57。專一性口袋：trypsin（Asp189，切 Lys/Arg 後）、chymotrypsin（疏水大，切 Phe/Tyr/Trp 後）、elastase（Val/Thr 填塞，切小 AA 後）。",
}


# ================================================================ reactions

REACTIONS = [
    {"id": "bpg_mutase", "title": "Bisphosphoglycerate mutase", "substrates": ["bpg_13"], "products": ["bpg_23"],
     "cofactors_out": ["H+"],
     "enzyme": "Bisphosphoglycerate mutase", "reversible": False,
     "card": {"location": "RBC 細胞質", "type": "分子內磷酸轉移：C1 的 acyl-phosphate 移到 C2（Rapoport–Luebering shunt）",
              "exam_points": ["繞過 PGK → 每經此路徑少產 1 ATP", "產物 2,3-BPG 使 Hb 曲線右移；高海拔、貧血、慢性缺氧時增加",
                              "儲存血 2,3-BPG 下降 → 輸血後 O₂ 釋放較差"]},
     "tags": ["glycolysis", "hemoglobin"]},
    {"id": "bpg_phosphatase", "title": "2,3-BPG phosphatase", "substrates": ["bpg_23"], "products": ["pg_3"],
     "cofactors_in": ["H2O"], "cofactors_out": ["Pi"], "enzyme": "2,3-BPG phosphatase", "reversible": False,
     "card": {"location": "RBC", "type": "水解 C2 磷酸 → 回到糖解的 3-PG（與 mutase 為同一雙功能酵素）"},
     "tags": ["glycolysis", "hemoglobin"]},
    {"id": "disulfide_reduction", "title": "Disulfide reduction (β-mercaptoethanol)",
     "substrates": ["cystine"], "products": [{"id": "cysteine", "n": 2}],
     "cofactors_in": ["2BME"], "cofactors_out": ["BME_ox"], "enzyme": "(chemical) β-mercaptoethanol / DTT",
     "reversible": True,
     "card": {"type": "硫醇–雙硫交換：還原劑的 –SH 把 –S–S– 還原成兩個 –SH（可逆）",
              "exam_points": ["SDS-PAGE 前加 β-ME/DTT 才會把以 –S–S– 相連的次單元分開",
                              "Anfinsen（RNase A）：urea + β-ME 變性，移除後自發復性 → 一級結構決定三級結構",
                              "performic acid 則把 –S–S– 氧化成 2 個 cysteic acid（–SO₃⁻，不可逆）"]},
     "tags": ["protein_structure"]},
    {"id": "prolyl_hydroxylase", "title": "Prolyl 4-hydroxylase",
     "substrates": ["proline", "alpha_kg"], "products": ["hydroxyproline", "succinate"],
     "cofactors_in": ["O2"], "cofactors_out": ["CO2"], "prosthetic": ["Fe^{2+}", "vitamin C"],
     "enzyme": "Prolyl 4-hydroxylase", "reversible": False,
     "card": {"location": "ER lumen（轉譯後修飾，作用於 collagen 鏈上的 Pro）",
              "type": "α-KG–dependent dioxygenase：O₂ 一個 O 進入 Hyp（洋紅），另一個進入 succinate；α-KG 的 C1 以 CO₂ 離開（青）",
              "cofactors": "Fe²⁺；ascorbate 把氧化成 Fe³⁺ 的鐵還原回 Fe²⁺",
              "deficiency": "Vitamin C 缺乏 → Hyp 不足 → collagen 三股螺旋不穩 → scurvy",
              "exam_points": ["Collagen：(Gly-X-Y)ₙ，X 常為 Pro、Y 常為 Hyp", "同類酵素：lysyl hydroxylase、HIF-1α prolyl hydroxylase（氧感測）"]},
     "tags": ["protein_structure"]},
    {"id": "lysyl_hydroxylase", "title": "Lysyl hydroxylase",
     "substrates": ["lysine", "alpha_kg"], "products": ["hydroxylysine", "succinate"],
     "cofactors_in": ["O2"], "cofactors_out": ["CO2"], "prosthetic": ["Fe^{2+}", "vitamin C"],
     "enzyme": "Lysyl hydroxylase", "reversible": False,
     "card": {"location": "ER lumen", "type": "同 prolyl hydroxylase（α-KG–dependent dioxygenase）",
              "exam_points": ["Hyl 可再接 Gal-Glc（O-linked）", "Lysyl oxidase（Cu²⁺）氧化 Lys/Hyl → 醛 → 共價交聯；Menkes（缺 Cu）交聯不良"]},
     "tags": ["protein_structure"]},
    {"id": "gamma_glutamyl_carboxylase", "title": "γ-Glutamyl carboxylase",
     "substrates": ["glutamate", "phylloquinol"], "products": ["gla", "vitamin_k_epoxide"],
     "cofactors_in": ["CO2", "O2"], "cofactors_out": ["H2O", "H+"],
     "enzyme": "γ-Glutamyl carboxylase", "reversible": False,
     "card": {"location": "肝細胞 ER（凝血因子 II、VII、IX、X，protein C、S）",
              "type": "Vitamin K hydroquinone 被氧化成 epoxide，提供能量拔除 Glu γ-H；CO₂（洋紅）接上成第二個 COO⁻",
              "regulation": "Warfarin 抑制 VKOR → KO 無法還原回 KH₂",
              "exam_points": ["Gla 的兩個 COO⁻ 螯合 Ca²⁺ → 因子結合磷脂膜", "新生兒缺 vitamin K → 出血（出生時注射）"]},
     "tags": ["protein_structure"]},
    {"id": "protein_kinase_ser", "title": "Protein kinase (Ser phosphorylation)",
     "substrates": ["serine"], "products": ["phosphoserine"], "cofactors_in": ["ATP"], "cofactors_out": ["ADP", "H+"],
     "enzyme": "Ser/Thr protein kinase", "reversible": False,
     "card": {"type": "ATP 的 γ-磷酸（洋紅）轉移到 Ser 側鏈 –OH（Ser 的 O 為親核基）",
              "regulation": "phosphatase 水解移除 → 可逆開關",
              "exam_points": ["Ser/Thr kinase：PKA、PKC、CaMK…；Tyr kinase：insulin receptor、EGFR、Src", "Mg²⁺–ATP 才是真正受質"]},
     "tags": ["protein_structure"]},
    {"id": "histone_acetyltransferase", "title": "Histone acetyltransferase (HAT)",
     "substrates": ["lysine", "acetyl_coa"], "products": ["acetyllysine"], "cofactors_out": ["CoA-SH", "H+"],
     "enzyme": "Histone acetyltransferase (HAT)", "reversible": False,
     "card": {"type": "acetyl（洋紅）由 acetyl-CoA 轉到組蛋白 Lys 的 ε-NH₂", "regulation": "HDAC 移除（HDAC inhibitor 為抗癌藥）",
              "exam_points": ["Lys 正電消失 → 與 DNA 結合變鬆 → euchromatin、gene on"]},
     "tags": ["protein_structure"]},
    {"id": "histone_methyltransferase", "title": "Histone methyltransferase (HMT)",
     "substrates": ["lysine"], "products": ["trimethyllysine"], "cofactors_in": ["3SAM"], "cofactors_out": ["3SAH", "3H+"],
     "enzyme": "Histone methyltransferase (HMT)", "reversible": False,
     "card": {"type": "SAM 的甲基（洋紅）轉到 Lys ε-N；可 mono/di/tri（圖為 tri）；離開的 H（青）以 H⁺ 釋出",
              "exam_points": ["甲基化不改變電荷", "H3K4me3 活化；H3K9me3、H3K27me3 抑制（heterochromatin）"]},
     "tags": ["protein_structure"]},
    {"id": "cytosine_deamination", "title": "Cytosine deamination", "substrates": ["cytosine"], "products": ["uracil", "ammonium"],
     "cofactors_in": ["H2O", "H+"], "enzyme": "spontaneous (hydrolytic) deamination", "reversible": False,
     "card": {"type": "C4 的 –NH₂（洋紅）以 NH₄⁺ 離開，H₂O 的 O（青）成為 C4=O",
              "exam_points": ["DNA 中出現 U → uracil-DNA glycosylase 辨認移除（BER）；這也是 DNA 用 T 而不用 U 的原因",
                              "未修復 → C:G → T:A 轉換突變"]},
     "tags": ["nucleotide", "degradation"]},
    {"id": "methylcytosine_deamination", "title": "5-Methylcytosine deamination",
     "substrates": ["methylcytosine"], "products": ["thymine", "ammonium"], "cofactors_in": ["H2O", "H+"],
     "enzyme": "spontaneous deamination", "reversible": False,
     "card": {"type": "5mC 去胺基直接得到 thymine（洋紅 NH₂ → NH₄⁺；青 O 來自 H₂O）",
              "exam_points": ["T 是正常鹼基，只能靠 T:G mismatch 修復辨認 → CpG 為突變熱點"]},
     "tags": ["nucleotide", "degradation"]},
]

COLORS.update({
    "bpg_mutase": {"bpg_13": {"colors": {0: {"fork_r": M}}}, "bpg_23": {"colors": {1: {"r": M}}}},
    "bpg_phosphatase": {"bpg_23": {"colors": {1: {"r": [M, "{P}"]}}}, "H2O": {"colors": {0: C}},
                        "Pi": {"colors": {0: M}}},
    "disulfide_reduction": {"cystine": {"colors": {"s1": C, "s2": C}},
                            "cysteine": {"colors": {3: {"c": [[C, "S"], [M, "H"]]}}},
                            "2BME": {"colors": {0: {"c": [[M, "H"]]}}}},
    "prolyl_hydroxylase": {"O2": {"colors": {0: M}}, "hydroxyproline": {"colors": {"oh": M}},
                           "alpha_kg": {"colors": {4: C}}, "CO2": {"colors": {0: C}}},
    "lysyl_hydroxylase": {"O2": {"colors": {0: M}}, "hydroxylysine": {"colors": {4: {"l": M}}},
                          "alpha_kg": {"colors": {4: C}}, "CO2": {"colors": {0: C}}},
    "gamma_glutamyl_carboxylase": {"CO2": {"colors": {0: M}}, "gla": {"colors": {3: {"r": M}}}},
    "protein_kinase_ser": {"ATP": {"colors": {0: M}}, "phosphoserine": {"colors": {2: [M, "{P}"]}}},
    "histone_acetyltransferase": {"acetyl_coa": {"colors": {0: {"c": M, "l": M, "u": M}}},
                                  "acetyllysine": {"colors": {7: M, 8: M}}},
    "histone_methyltransferase": {"3SAM": {"colors": {0: M}}, "trimethyllysine": {"colors": {6: {"l": M, "r": M}, 7: M}},
                                  "lysine": {"colors": {6: [C, "H_3"]}}, "3H+": {"colors": {0: C}}},
    "cytosine_deamination": {"cytosine": {"highlight": [{"smarts": "[NH2:1]c1nc(=O)[nH]cc1", "color": M}]},
                             "ammonium": {"colors": {0: M}}, "H2O": {"colors": {0: C}},
                             "uracil": {"highlight": [{"smarts": "[O:1]=c1[cH][cH][nH]c(=O)[nH]1", "color": C}]}},
    "methylcytosine_deamination": {"methylcytosine": {"highlight": [{"smarts": "[NH2:1]c1nc(=O)[nH]cc1C", "color": M}]},
                                   "ammonium": {"colors": {0: M}}, "H2O": {"colors": {0: C}},
                                   "thymine": {"highlight": [{"smarts": "[O:1]=c1[nH]c(=O)[nH]cc1C", "color": C}]}},
})

NEW_COFACTORS = {
    "2BME": {"label": "2 HSCH_2CH_2OH", "smiles": "OCCS.OCCS"},
    "BME_ox": {"label": "HOCH_2CH_2S–SCH_2CH_2OH", "smiles": "OCCSSCCO"},
    "3SAM": {"label": "3 SAM", "smiles": "C[S+]([14*])[14*].C[S+]([14*])[14*].C[S+]([14*])[14*]", "molecule": "sam"},
    "3SAH": {"label": "3 SAH", "smiles": "[14*]S[14*].[14*]S[14*].[14*]S[14*]", "molecule": "sah"},
}


# ================================================================ pathways

PATHWAYS = [
    {"id": "inhibition_mixed", "title": "Mixed inhibition", "layout": "plot", "panels": INHIBITION_MIXED,
     "caption": {"text": "Mixed inhibition: I binds E (α) and ES (α') with different affinity → V_{max} ↓, K_m ↑ or ↓", "color": RED},
     "tags": ["enzyme_kinetics", "inhibition"],
     "notes": "Mixed inhibition：抑制劑可結合 E 與 ES，但親和力不同（α ≠ α'）。V_max^app = V_max/α' 一定下降；"
              "α > α'（較偏好 E）→ K_m 上升，雙倒數直線交在 y 軸左側、x 軸上方（常見）；α < α'（較偏好 ES）→ K_m 下降，交點在 x 軸下方（少見）。"
              "α = α' 即為 pure noncompetitive。"},
    {"id": "substrate_inhibition", "title": "Substrate inhibition", "layout": "plot", "panels": SUBSTRATE_INHIBITION,
     "tags": ["enzyme_kinetics", "inhibition"],
     "notes": "受質抑制：[S] 過高時第二個 S 結合在抑制位（或形成無生產力的 ESS）→ 速率下降，v₀ 有最佳 [S]（≈ √(K_m K_i)）。"
              "PFK-1：ATP 在低濃度是受質（活性位），高濃度時結合別位當抑制劑（能量充足時減慢糖解）。"},
    {"id": "burst_kinetics", "title": "Burst kinetics (chymotrypsin)", "layout": "plot", "plot": BURST,
     "tags": ["enzyme_kinetics", "catalysis"],
     "notes": "Chymotrypsin 水解 p-nitrophenyl acetate：開始時快速釋出一當量 p-nitrophenolate（burst，acylation 快），"
              "之後速率由慢的 deacylation（水解 acyl-enzyme）決定 → 雙相曲線。外推直線的 y 截距 = 活性酵素濃度 [E]_T，證明有 covalent acyl-enzyme 中間物（ping-pong）。"},
    {"id": "reaction_coordinate", "title": "Reaction coordinate: enzyme lowers ΔG‡", "layout": "plot",
     "plot": REACTION_COORD, "tags": ["enzyme", "catalysis"],
     "notes": "酵素降低活化能 ΔG‡（穩定過渡態、binding energy）→ 正逆反應都加快；但不改變 ΔG、K_eq 或平衡方向。"
              "ΔG < 0 自發 S → P；ΔG°′ = −RT ln K′_eq（ΔG°′ < 0 ⇔ K_eq > 1）。"},
    {"id": "lysozyme_ph_profile", "title": "Lysozyme activity vs pH", "layout": "plot", "plot": LYSOZYME_PH,
     "tags": ["enzyme_kinetics", "catalysis"],
     "notes": "Lysozyme：Glu35 需保持質子化（COOH，general acid 給 H⁺ 給醣苷鍵 O），Asp52 需去質子化（COO⁻，穩定 oxocarbenium 離子）。"
              "兩者同時成立的比例在兩個 pK_a 中點最大 → 最適 pH ≈ (4.5 + 5.9)/2 ≈ 5.2。同理：活性位兩個可解離基團 → 鐘形 pH 曲線。"},
    {"id": "hb_t_r_states", "title": "Hb switches from T to R state", "layout": "plot", "plot": HB_TR,
     "tags": ["hemoglobin", "allostery"],
     "notes": "若 Hb 只有 R 態（≈ Mb）：肺能吸 O₂ 但組織放不出；只有 T 態：組織放得出但肺吸不到。"
              "實際 Hb 在低 pO₂ 時沿 T 曲線，結合 O₂ 後逐步轉成 R 態 → S 形曲線，同時兼顧肺部攝取與組織卸氧（MWC 模型計算）。"},
    {"id": "titration_his_glu", "title": "Titration curves: histidine and glutamate", "layout": "plot",
     "panels": TITRATION_HIS_GLU, "tags": ["amino_acid"],
     "notes": "三個可解離基團 → 三個緩衝區。pI = 淨電荷 0 那一型兩側的 pK 平均：His（鹼性）pI = (6.0 + 9.2)/2 ≈ 7.6；"
              "Glu（酸性）pI = (2.2 + 4.3)/2 ≈ 3.2。His 的 pK_R ≈ 6 接近生理 pH → 蛋白中的緩衝（如 Hb His146 的 Bohr effect）、酵素 general acid/base。"},
    {"id": "dna_melting_curve", "title": "DNA melting curve", "layout": "plot", "plot": MELTING,
     "tags": ["nucleotide"],
     "notes": "加熱 dsDNA → 協同性解旋（S 形），A₂₆₀ 上升約 40 %（hyperchromic effect：疊在一起的鹼基吸光較少）。T_m = 一半解開的溫度。"
              "T_m 上升：G+C 含量高（3 個氫鍵、堆疊較強）、鹽濃度高（遮蔽磷酸負電）；T_m 下降：低鹽、極端 pH、formamide/urea。"
              "ssRNA 只有局部配對 → 曲線平緩、無協同性。冷卻復性（renaturation）→ hypochromic。"},
    {"id": "aldose_ketose_isomerization", "title": "Aldose ⇌ enediol ⇌ ketose", "layout": "free", "cell_size": 7.5,
     "nodes": [{"id": "glc", "ref": "glucose", "pos": [0, 0]}, {"id": "ene", "ref": "enediol", "pos": [1.35, 0]},
               {"id": "fru", "ref": "d_fructose", "pos": [2.7, 0]}, {"id": "man", "ref": "d_mannose", "pos": [1.35, 1.95]}],
     "edges": [{"from": "glc", "to": "ene", "reversible": True, "label": "base / isomerase"},
               {"from": "ene", "to": "fru", "reversible": True},
               {"from": "ene", "to": "man", "reversible": True, "label": "C2 epimer"}],
     "tags": ["carbohydrate", "monosaccharide"],
     "notes": "互變異構（tautomerism）：aldose 的 C1 醛基與 C2–H 互換 → cis-enediol → ketose。鹼性下 glucose、mannose、fructose 互相轉換"
              "（Lobry de Bruyn–van Ekenstein）；酵素版即 phosphoglucose isomerase（G6P ⇌ F6P）、triose phosphate isomerase。"
              "Fructose 在鹼性 Fehling 試劑中也呈還原性，就是因為能轉成 enediol／aldose。"},
    {"id": "glucose_oxidation_products", "title": "Glucose: oxidation and reduction products", "layout": "free",
     "cell_size": 7.5,
     "nodes": [{"id": "glc", "ref": "glucose", "pos": [1.3, 0]}, {"id": "gna", "ref": "d_gluconate", "pos": [-0.1, 1.95]},
               {"id": "gua", "ref": "d_glucuronate_open", "pos": [1.3, 1.95]}, {"id": "sor", "ref": "sorbitol", "pos": [2.7, 1.95]}],
     "edges": [{"from": "glc", "to": "gna", "label": "C1 oxidation\n(glucose oxidase)", "label_side": "right"},
               {"from": "glc", "to": "gua", "label": "C6 oxidation\n(via UDP-Glc DH)", "label_side": "right"},
               {"from": "glc", "to": "sor", "label": "C1 reduction\n(aldose reductase)"}],
     "molecule_deco": "default",
     "tags": ["carbohydrate", "sugar_derivative"],
     "notes": "Aldonic acid（C1 → COOH，如 gluconate）；uronic acid（C6 → COOH，如 glucuronate）；alditol（C1 醛 → CH₂OH，如 sorbitol，"
              "注意是「還原」）。Glucuronate → 膽紅素結合、L-ascorbate（人類不能）。"},
    {"id": "rapoport_luebering_shunt", "title": "2,3-BPG shunt (RBC)", "layout": "free", "cell_size": 7.0,
     "nodes": [{"id": "b13", "ref": "bpg_13", "pos": [0, 0]}, {"id": "pg3", "ref": "pg_3", "pos": [2.2, 0]},
               {"id": "b23", "ref": "bpg_23", "pos": [1.1, 1.35]}],
     "edges": [{"from": "b13", "to": "pg3", "step": "phosphoglycerate_kinase", "label": "PGK", "side_text": "ADP → ATP"},
               {"from": "b13", "to": "b23", "step": "bpg_mutase", "label": "BPG mutase", "label_side": "right"},
               {"from": "b23", "to": "pg3", "step": "bpg_phosphatase", "label": "2,3-BPG\nphosphatase", "side_text": "H_2O → P_i", "label_side": "right"}],
     "tags": ["glycolysis", "hemoglobin"],
     "notes": "RBC 糖解的支路：1,3-BPG → 2,3-BPG → 3-PG，繞過 PGK（犧牲 1 ATP）換取調節 Hb 的 2,3-BPG。"
              "Pyruvate kinase 缺乏 → 上游 2,3-BPG 堆積 → 曲線右移；hexokinase 缺乏 → 2,3-BPG 下降 → 左移。"},
    {"id": "bohr_effect_co2_transport", "title": "Bohr effect and CO2 transport", "layout": "free", "cell_size": 5.0,
     "regions": [{"kind": "box", "x0": -0.75, "y0": -0.6, "x1": 2.45, "y1": 2.6, "label": "Tissue capillary"},
                 {"kind": "box", "x0": 2.6, "y0": -0.6, "x1": 5.7, "y1": 2.6, "label": "Lung"}],
        "nodes": [
            t("co2t", "CO_2 + H_2O", 0.3, 0.15, "ion"), t("h2co3", "H_2CO_3", 0.3, 1.05, "ion"),
            t("hco3", "HCO_3^- + H^+", 0.3, 1.95, "ion"),
            t("hbo2", "HbO_2 (R)", 1.5, 0.15, "protein"), t("hbh", "H^+–Hb (T)\n+ O_2 → tissue", 1.45, 1.95, "protein"),
            t("hbh2", "H^+–Hb (T)", 3.35, 1.95, "protein"), t("hbo2l", "HbO_2 (R) + H^+", 3.45, 0.15, "protein"),
            t("co2l", "CO_2 + H_2O\n(exhaled)", 5.0, 0.15, "ion"), t("hco3l", "HCO_3^- + H^+", 5.0, 1.95, "ion"),
        ],
        "edges": [
            {"from": "co2t", "to": "h2co3", "label": "carbonic\nanhydrase"}, {"from": "h2co3", "to": "hco3"},
            {"from": "hco3", "to": "hbh", "color": "magenta", "label": "H^+ binds\nHis146 (β)"},
            {"from": "hbo2", "to": "hbh", "label": "O_2 released"},
            {"from": "hbh", "to": "hbh2", "label": "blood"},
            {"from": "hbh2", "to": "hbo2l", "label": "O_2 binds\n→ H^+ released"},
            {"from": "hbo2l", "to": "hco3l", "color": "magenta"}, {"from": "hco3l", "to": "co2l", "label": "carbonic\nanhydrase"},
        ],
     "tags": ["hemoglobin"],
     "notes": "Bohr effect：組織 CO₂ 多 → carbonic anhydrase 產生 H⁺ → H⁺ 結合 Hb（β His146 與 Asp94 形成鹽橋）→ 穩定 T 態、放出 O₂。"
              "肺部 O₂ 高 → 轉成 R 態釋出 H⁺ → 與 HCO₃⁻ 變回 CO₂ 呼出。CO₂ 也直接接在 N 端胺基形成 carbamate（–NH–COO⁻），同樣穩定 T 態。"
              "運輸形式：HCO₃⁻（~70 %，chloride shift）、carbamino-Hb（~20 %）、溶解（~10 %）。"},
    {"id": "pufa_synthesis", "title": "PUFA synthesis: desaturases and elongases", "layout": "free", "cell_size": 5.4,
     "nodes": [
         n("ole", "oleate", 0, 0), n("la", "linoleate", 0, 1.2), n("ala", "alpha_linolenate", 0, 2.4),
         t("gla", "γ-Linolenate\n18:3 (ω6)", 1.55, 1.2, "name"), t("dgla", "Dihomo-γ-\nlinolenate\n20:3 (ω6)", 3.1, 1.2, "name"),
         n("aa", "arachidonate", 4.75, 1.2), t("sda", "18:4 (ω3)", 1.55, 2.4, "name"), t("eta", "20:4 (ω3)", 3.1, 2.4, "name"),
         n("epa", "epa", 4.75, 2.4), n("dha", "dha", 4.75, 3.6),
         t("plants", "plants only", -1.1, 0.6, "name"),
     ],
     "edges": [
         {"from": "ole", "to": "la", "label": "Δ12", "color": "green"}, {"from": "la", "to": "ala", "label": "Δ15", "color": "green"},
         {"from": "la", "to": "gla", "label": "Δ6\ndesaturase"}, {"from": "gla", "to": "dgla", "label": "elongase"},
         {"from": "dgla", "to": "aa", "label": "Δ5\ndesaturase"},
         {"from": "ala", "to": "sda", "label": "Δ6"}, {"from": "sda", "to": "eta", "label": "elongase"},
         {"from": "eta", "to": "epa", "label": "Δ5"}, {"from": "epa", "to": "dha", "label": "elongase, Δ4\n(+ β-oxidation)", "label_side": "right"},
     ],
     "tags": ["lipid", "fatty_acid"],
     "notes": "人體去飽和酶只能在 C9 與羧基之間加雙鍵（Δ9、Δ6、Δ5、Δ4），沒有 Δ12、Δ15（綠色，只有植物）→ linoleate（ω6）與 α-linolenate（ω3）為必需脂肪酸。"
              "延長與去飽和保留 ω 位置：ω6 → arachidonate；ω3 → EPA → DHA（魚油）。ω 編號：從甲基端數第一個雙鍵；n−x 同義。"},
    {"id": "eicosanoid_series", "title": "Eicosanoid series 1, 2, 3", "layout": "free", "cell_size": 5.0,
     "nodes": [
         t("dgla", "DGLA\n20:3 (ω6)", 0, 1.2, "name"), n("aa", "arachidonate", 2.3, 1.2), n("epa", "epa", 4.6, 1.2),
         t("pg1", "PGE_1, TXA_1\n(series 1)", 0, 0, "concept"), t("pg2", "PGE_2, PGI_2, TXA_2\n(series 2)", 2.3, 0, "concept"),
         t("pg3", "PGI_3, TXA_3\n(series 3)", 4.6, 0, "concept"),
         t("lt3", "LTA_3\n(series 3)", 0, 2.4, "concept"), t("lt4", "LTA_4, LTB_4, LTC_4\n(series 4)", 2.3, 2.4, "concept"),
         t("lt5", "LTA_5\n(series 5)", 4.6, 2.4, "concept"),
     ],
     "edges": [
         {"from": "dgla", "to": "pg1", "label": "COX"}, {"from": "aa", "to": "pg2", "label": "COX\n(aspirin\ninhibits)"},
         {"from": "epa", "to": "pg3", "label": "COX"},
         {"from": "dgla", "to": "lt3", "label": "LOX"}, {"from": "aa", "to": "lt4", "label": "5-LOX"},
         {"from": "epa", "to": "lt5", "label": "LOX"},
     ],
     "tags": ["lipid", "eicosanoid"],
     "notes": "前列腺素／血栓素的下標 = 前驅物雙鍵數 − 2（COX 環化用掉 2 個）；白三烯的下標 = 雙鍵數（LOX 不形成環）。"
              "DGLA（20:3）→ series 1 PG、LT₃；AA（20:4）→ series 2 PG/TX、LT₄；EPA（20:5）→ series 3、LT₅（魚油：TXA₃ 促凝效果弱 → 抗血栓）。"},
    {"id": "isoprenoid_family", "title": "Isoprenoids derived from IPP", "layout": "free", "cell_size": 4.6,
     "nodes": [
         n("ipp", "ipp", 3.2, 0), n("chol", "cholesterol", 1.0, 1.2), n("bile", "cholate", -0.3, 2.4),
         n("vd", "calcitriol", 1.0, 3.2), n("ster", "cortisol", 2.3, 2.4),
         t("vit", "Vitamins A, E, K", 4.6, 1.2, "name"), n("ret", "retinol", 3.6, 2.4), n("toc", "tocopherol", 4.6, 3.2),
         n("vk", "phylloquinone", 5.7, 2.4),
         t("other", "CoQ, dolichol,\nfarnesyl / geranylgeranyl\n(protein prenylation),\nheme A, carotenoids", 6.9, 1.2, "name"),
     ],
     "edges": [
         {"from": "ipp", "to": "chol", "head": False}, {"from": "ipp", "to": "vit", "head": False},
         {"from": "ipp", "to": "other", "head": False},
         {"from": "chol", "to": "bile", "head": False}, {"from": "chol", "to": "vd", "head": False},
         {"from": "chol", "to": "ster", "head": False},
         {"from": "vit", "to": "ret", "head": False}, {"from": "vit", "to": "toc", "head": False},
         {"from": "vit", "to": "vk", "head": False},
     ],
     "tags": ["lipid", "steroid"],
     "notes": "Isoprene（C₅）單位 = IPP／DMAPP（由 mevalonate 途徑；植物與細菌另有 MEP 途徑）。"
              "衍生：cholesterol（→ 膽酸、vitamin D、固醇類荷爾蒙）、脂溶性維生素 A、E、K（D 由 cholesterol 而來）、"
              "CoQ（ubiquinone 的側鏈）、dolichol（N-linked glycosylation 載體）、farnesyl／geranylgeranyl（Ras 等蛋白質 prenylation）、heme A、carotenoid／lycopene。"},
    {"id": "phosphoinositide_signaling", "title": "PIP2 signalling: PLC and PI3K", "layout": "free", "cell_size": 5.0,
     "nodes": [
         t("pip2", "PIP_2\n(PI-4,5-P_2)", 0, 0, "name"), t("dag", "DAG", 1.8, -0.6, "name"),
         t("ip3", "IP_3", 1.8, 0.6, "name"), t("pkc", "PKC", 3.4, -0.6, "protein"),
         t("er", "ER IP_3 receptor\n(Ca^{2+} channel)", 3.4, 0.6, "protein"), t("ca", "Ca^{2+} ↑", 4.9, 0.6, "ion"),
         t("pip3", "PIP_3\n(PI-3,4,5-P_3)", 0, 1.8, "name"), t("pdk", "PDK1", 1.8, 1.8, "protein"),
         t("akt", "PKB / Akt", 3.4, 1.8, "protein"), t("sig", "survival, GLUT4,\nglycogen synthesis", 5.3, 1.8, "name"),
     ],
     "edges": [
         {"from": "pip2", "to": "dag", "label": "PLC"}, {"from": "pip2", "to": "ip3", "label": "PLC"},
         {"from": "dag", "to": "pkc"}, {"from": "ip3", "to": "er"}, {"from": "er", "to": "ca"},
         {"from": "ca", "to": "pkc", "label": "+"},
         {"from": "pip2", "to": "pip3", "label": "PI3K\n(insulin)"}, {"from": "pip3", "to": "pdk"},
         {"from": "pdk", "to": "akt"}, {"from": "akt", "to": "sig"},
     ],
     "tags": ["lipid", "phospholipid"],
     "notes": "Gq → PLC-β 切 PIP₂：DAG 留在膜上活化 PKC；IP₃ 開啟 ER 的 IP₃ receptor 釋出 Ca²⁺（也協同活化 PKC）。"
              "Insulin／RTK → PI3K 把 PIP₂ 的 3 位磷酸化成 PIP₃ → 招募 PDK1、Akt（PKB）→ Akt 被磷酸化活化。PTEN 把 PIP₃ 去磷酸回 PIP₂。"},
    {"id": "enzyme_classes", "title": "Seven classes of enzymes (EC)", "layout": "free", "cell_size": 4.0,
     "nodes": [
         t("c1", "1 Oxidoreductase", 0, 0, "concept"), t("c2", "2 Transferase", 0, 0.8, "concept"),
         t("c3", "3 Hydrolase", 0, 1.6, "concept"), t("c4", "4 Lyase", 0, 2.4, "concept"),
         t("c5", "5 Isomerase", 0, 3.2, "concept"), t("c6", "6 Ligase", 0, 4.0, "concept"),
         t("c7", "7 Translocase", 0, 4.8, "concept"),
         t("e1", "Lactate DH (NADH)", 3.2, 0, "name"), t("e2", "Hexokinase (Ⓟ from ATP)", 3.2, 0.8, "name"),
         t("e3", "Arginase (H_2O)", 3.2, 1.6, "name"), t("e4", "Fumarase (adds H_2O to C=C)", 3.2, 2.4, "name"),
         t("e5", "Phosphoglucomutase", 3.2, 3.2, "name"), t("e6", "Pyruvate carboxylase (ATP)", 3.2, 4.0, "name"),
         t("e7", "ATP synthase / Na^+,K^+-ATPase", 3.2, 4.8, "name"),
     ],
     "edges": [
         {"from": "c1", "to": "e1", "step": "ldh"}, {"from": "c2", "to": "e2", "step": "hexokinase"},
         {"from": "c3", "to": "e3", "step": "arginase"}, {"from": "c4", "to": "e4", "step": "fumarase"},
         {"from": "c5", "to": "e5", "step": "phosphoglucomutase"}, {"from": "c6", "to": "e6", "step": "pyruvate_carboxylase"},
         {"from": "c7", "to": "e7", "step": "atp_synthase"},
     ],
     "tags": ["enzyme"],
     "notes": "EC 分類（口訣 OTHLIL＋T）：1 氧化還原（dehydrogenase、oxidase、oxygenase）；2 轉移官能基（kinase、transaminase）；3 水解（protease、lipase、nuclease、phosphatase）；"
              "4 Lyase：非水解、非氧化地加成或移除基團形成雙鍵（decarboxylase、aldolase、hydratase；synthase 不耗 ATP，如 citrate synthase 屬 lyase）；"
              "5 異構（mutase、epimerase、racemase）；6 Ligase：耗 ATP 形成鍵（synthetase、carboxylase）；7 Translocase：跨膜運送（2018 新增）。點箭頭可跳到例子反應。"},
    {"id": "g_protein_cycle", "title": "G protein cycle", "layout": "cycle", "radius": 2.1, "cell_size": 5.2,
     "start_angle": -90, "center_text": "Gα\nGTPase\nswitch",
     "nodes": [
         {"id": "gdp", "text": "Gαβγ–GDP\n(inactive)", "kind": "protein"},
         {"id": "ex", "text": "Receptor (GPCR) + ligand\nGDP → GTP exchange", "kind": "name"},
         {"id": "gtp", "text": "Gα–GTP + βγ\n(active)", "kind": "protein"},
         {"id": "eff", "text": "Effector: adenylyl cyclase\n(cAMP), PLC-β …", "kind": "name"},
         {"id": "hyd", "text": "GTP → GDP + P_i\n(intrinsic GTPase, RGS)", "kind": "name"},
     ],
     "edges": [{"from": "gdp", "to": "ex"}, {"from": "ex", "to": "gtp"}, {"from": "gtp", "to": "eff"},
               {"from": "eff", "to": "hyd"}, {"from": "hyd", "to": "gdp"}],
     "tags": ["enzyme", "allostery"],
     "notes": "GPCR 活化 → 促使 Gα 釋出 GDP、結合 GTP（GEF 作用）→ Gα–GTP 與 βγ 分開並各自作用（Gs 活化、Gi 抑制 adenylyl cyclase；Gq 活化 PLC-β）。"
              "Gα 本身的 GTPase 把 GTP 水解成 GDP → 回到三聚體（關閉）。霍亂毒素：Gs ADP-ribosylation → GTPase 失效 → 持續活化；百日咳毒素：Gi 無法活化。"
              "PKA：4 cAMP 結合 2 個 R 次單元 → 釋出 2 個活性 C 次單元。"},
]


# ================================================================ writer

def dump(path, d):
    """Write the YAML, keeping a review status given in the web UI and an `illus` variant added by gen_illus_sig.py
    (re-running this generator must not reset either)."""
    if path.exists():
        old = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        if "status" in d and old.get("status") in ("approved", "changes_requested"):
            d["status"] = old["status"]
        illus = (old.get("variants") or {}).get("illus")
        if illus and "illus" not in (d.get("variants") or {}):
            d["variants"] = {**(d.get("variants") or {}), "illus": illus}
    path.write_text(yaml.safe_dump(d, allow_unicode=True, sort_keys=False, width=200), encoding="utf-8")


PARENT = {"enzyme_kinetics": "enzyme", "inhibition": "enzyme", "allostery": "enzyme", "catalysis": "enzyme",
          "hemoglobin": "protein", "peptide": "protein", "protein_structure": "protein",
          "glycolysis": "carb_metabolism", "fatty_acid": "lipid", "eicosanoid": "lipid", "steroid": "lipid",
          "phospholipid": "lipid", "disaccharide": "carbohydrate", "sugar_derivative": "carbohydrate",
          "monosaccharide": "carbohydrate", "polar": "amino_acid", "pyrimidine": "nucleotide",
          "degradation": "nucleotide", "base_pairing": "nucleotide"}


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
        lines = ["", "# --- 第十五輪（筆記補齊）：β-mercaptoethanol、3 SAM／SAH（由 tools/gen_notes.py 產生）"]
        for k in new:
            v = NEW_COFACTORS[k]
            extra = f', molecule: {v["molecule"]}' if v.get("molecule") else ""
            lines.append(f'{k}: {{label: "{v["label"]}", smiles: "{v["smiles"]}"{extra}}}')
        p.write_text(text.rstrip("\n") + "\n" + "\n".join(lines) + "\n", encoding="utf-8")


def main(data):
    data = Path(data)
    add_cofactors(data)
    mols = MOLECULES + [dict(b, tags=["nucleotide", "base_pairing"]) for b in BASE_PAIRS] + [PHOSPHOLIPASE, TRIAD]
    for m in mols:
        d = {"id": m["id"], "name": m["name"]}
        if m.get("label"):
            d["label"] = m["label"]
        d["aliases"] = m.get("aliases", [])
        d["smiles"] = m["smiles"]
        for k in ("chain", "free", "haworth", "ring", "decorate"):
            if m.get(k):
                d[k] = m[k]
        d["tags"] = with_parent(m["tags"])
        d["notes_zh"] = m["notes"]
        d["status"] = "pending"
        dump(data / "molecules" / f"{m['id']}.yaml", d)
    for r in REACTIONS:
        d = {"id": r["id"], "title": r["title"], "substrates": r["substrates"], "products": r["products"]}
        for k in ("cofactors_in", "cofactors_out", "prosthetic"):
            if r.get(k):
                d[k] = r[k]
        d["enzyme"] = r["enzyme"]
        d["reversible"] = r["reversible"]
        if COLORS.get(r["id"]):
            d["decorate"] = COLORS[r["id"]]
        d["card_zh"] = r["card"]
        d["tags"] = with_parent(r["tags"])
        d["status"] = "pending"
        dump(data / "reactions" / f"{r['id']}.yaml", d)
    for p in PATHWAYS:
        d = {k: v for k, v in p.items() if k not in ("tags", "notes")}
        d["tags"] = with_parent(p["tags"])
        d["notes_zh"] = p["notes"]
        d["status"] = "pending"
        dump(data / "pathways" / f"{p['id']}.yaml", d)
    print(f"{len(mols)} molecules, {len(REACTIONS)} reactions, {len(PATHWAYS)} pathways")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data")
