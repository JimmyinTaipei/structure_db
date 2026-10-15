"""Round 17: signal transduction and plasma membrane (Book 3 handout 生化_B3_訊號傳遞.pdf).

    python tools/gen_sig.py data

Molecules (cAMP, cGMP, GTP), reactions (adenylyl / guanylyl cyclase, cAMP / cGMP phosphodiesterase),
overviews (receptor types, GPCR classes, cAMP–PKA cascade, NO–cGMP relaxation, insulin / RTK → MAPK and
PI3K, JAK-STAT, NF-κB, apoptosis, GTPase switch, p53–Rb, two-component system, bacterial toxins,
second messengers), membrane figures (transporters, lipid asymmetry, LDL endocytosis) and a fluidity plot.
"""
import sys
from pathlib import Path

import yaml
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).parent))
from reaction_colors import COLORS, hl, hw, whole  # noqa: E402

M, C = "magenta", "cyan"
RED, GREEN, ORANGE, GRAY = "red", "green", "orange", "gray"


def n(id_, ref, x, y, **kw):
    return {"id": id_, "ref": ref, "pos": [x, y], "view": "name", **kw}


def t(id_, text, x, y, kind="concept", **kw):
    return {"id": id_, "text": text, "pos": [x, y], "kind": kind, **kw}


# ================================================================ molecules

AMP = "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]1O"
GMP = "O=c1[nH]c(N)nc2c1ncn2[C@@H]1O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]1O"
GTP = GMP.replace("COP(=O)([O-])[O-]", "COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-]")


def cyclic(smi):
    """3',5'-cyclic nucleotide: bond the 3'-O to the 5'-phosphate P and drop one P–O⁻ (sugar stereo untouched)."""
    mol = Chem.RWMol(Chem.MolFromSmiles(smi))
    patt = Chem.MolFromSmarts("[O-]P(=O)([O-])OC[C;R]1[O;R][C;R][C;R][C;R]1[OX2H1]")
    o_minus, p, _, _, _, _, _, _, _, _, _, o3 = mol.GetSubstructMatch(patt)
    mol.AddBond(o3, p, Chem.BondType.SINGLE)
    mol.GetAtomWithIdx(o3).SetNumExplicitHs(0)
    mol.RemoveAtom(o_minus)
    m = mol.GetMol()
    Chem.SanitizeMol(m)
    return Chem.MolToSmiles(m)


NUC_HW = {1: {"up": "base", "down": "H"}, 2: {"up": "H", "down": "OH"}, 3: {"up": "H", "down": "OH"}}

MOLECULES = [
    {"id": "camp", "name": "cAMP", "label": "cAMP (3',5'-cyclic AMP)", "aliases": ["cyclic AMP", "環腺苷酸"],
     "smiles": cyclic(AMP), "ring": {"coordgen": True},
     "tags": ["second_messenger", "nucleotide"],
     "notes": "第二傳訊者。Adenylyl cyclase 由 ATP 生成（3′-OH 攻擊 α-P，放出 PPi），phosphodiesterase 水解 3′ 酯鍵 → 5′-AMP。"
              "4 cAMP 結合 PKA 的兩個 R 次單元 → 釋出 2 個活性 C 次單元。Caffeine、theophylline 抑制 PDE → cAMP 上升。"},
    {"id": "cgmp", "name": "cGMP", "label": "cGMP (3',5'-cyclic GMP)", "aliases": ["cyclic GMP", "環鳥苷酸"],
     "smiles": cyclic(GMP), "ring": {"coordgen": True, "rotate": 180},
     "tags": ["second_messenger", "nucleotide"],
     "notes": "Guanylyl cyclase 產生：膜型（ANF 受體）與可溶型（NO 受體，含 heme）。活化 PKG → 平滑肌放鬆；視網膜中 cGMP 打開 Na⁺ 通道，"
              "光 → transducin → PDE6 分解 cGMP → 通道關閉（過極化）。PDE5 抑制劑 sildenafil 延長 cGMP 作用。"},
    {"id": "gtp", "name": "GTP", "aliases": ["guanosine triphosphate", "鳥苷三磷酸"], "smiles": GTP,
     "haworth": {"ring": "furanose", "subs": {**NUC_HW, 4: {"up": "{P}–{P}–{P}–OCH_2", "down": "H"}},
                 "base": "[*]n1cnc2c(=O)[nH]c(N)nc21", "base_right": 5, "numbers": True},
     "tags": ["nucleotide", "purine"],
     "notes": "G 蛋白與小 GTPase 的開關分子（結合 GTP = on）；轉譯（IF2、EF-Tu、EF-G、RF3）、糖質新生 PEPCK、TCA succinyl-CoA synthetase 的產物；guanylyl cyclase 的受質。"},
]

# ================================================================ reactions

REACTIONS = [
    {"id": "adenylyl_cyclase", "title": "ATP → cAMP + PPi (adenylyl cyclase)",
     "substrates": ["atp_n"], "products": ["camp"], "cofactors_out": ["PPi"], "prosthetic": ["Mg^{2+}"],
     "enzyme": "adenylyl cyclase", "reversible": False,
     "card": {"location": "細胞膜內側（12 次跨膜的膜蛋白）",
              "type": "分子內磷酸轉移：3′-OH（青）攻擊 α-P，β、γ 以 PPi（洋紅）離開，形成 3′,5′-環狀磷酸二酯",
              "regulation": "Gsα–GTP 活化、Giα 抑制；forskolin 直接活化；cholera toxin 使 Gsα 持續活化 → cAMP 持續上升",
              "exam_points": ["PPi 被 pyrophosphatase 水解 → 推動反應", "cAMP → PKA（R₂C₂ + 4 cAMP → 2 C）、Epac、CNG 通道",
                              "β-adrenergic、glucagon、ACTH、TSH、PTH、LH/FSH、ADH(V₂) 都經 Gs–cAMP"]},
     "tags": ["second_messenger", "gpcr"]},
    {"id": "camp_phosphodiesterase", "title": "cAMP → AMP (cAMP phosphodiesterase)",
     "substrates": ["camp"], "products": ["amp"], "cofactors_in": ["H2O"], "extra_out": ["H+"],
     "enzyme": "cAMP phosphodiesterase", "reversible": False,
     "card": {"type": "水解 3′ 磷酸酯鍵（H₂O 的 O 進入磷酸，青）→ 5′-AMP；終止 cAMP 訊號",
              "regulation": "insulin 經 PKB 活化 PDE3B → 降低 cAMP（抗 glucagon 作用）",
              "exam_points": ["Caffeine、theophylline 抑制 PDE → cAMP ↑（心跳、支氣管擴張）"]},
     "tags": ["second_messenger", "gpcr"]},
    {"id": "guanylyl_cyclase", "title": "GTP → cGMP + PPi (guanylyl cyclase)",
     "substrates": ["gtp"], "products": ["cgmp"], "cofactors_out": ["PPi"], "prosthetic": ["Mg^{2+}"],
     "enzyme": "guanylyl cyclase", "reversible": False,
     "card": {"location": "膜型（ANF / natriuretic peptide 受體本身）；可溶型（細胞質，含 heme，NO 受體）",
              "type": "同 adenylyl cyclase：3′-OH 攻擊 α-P，放出 PPi",
              "regulation": "NO 結合可溶型 GC 的 heme Fe²⁺ → 活化（nitroglycerin、nitroprusside 放出 NO）",
              "exam_points": ["NO–cGMP–PKG → 平滑肌放鬆、血管擴張", "ANF → 膜型 GC → 排鈉利尿、血管擴張"]},
     "tags": ["second_messenger"]},
    {"id": "cgmp_phosphodiesterase", "title": "cGMP → GMP (cGMP phosphodiesterase, PDE5)",
     "substrates": ["cgmp"], "products": ["gmp"], "cofactors_in": ["H2O"], "extra_out": ["H+"],
     "enzyme": "cGMP phosphodiesterase (PDE5)", "reversible": False,
     "card": {"type": "水解 cGMP 的 3′ 酯鍵 → 5′-GMP；終止 NO 訊號",
              "exam_points": ["Sildenafil（Viagra）抑制 PDE5 → cGMP 持續 → 海綿體平滑肌放鬆", "視網膜 PDE6 由 transducin 活化（光訊號）"]},
     "tags": ["second_messenger"]},
]

COLORS.update({
    "adenylyl_cyclase": {"atp_n": hw(d3=[C, "O"], u4=[M, "{P}–{P}"]), "PPi": whole(),
                         "camp": hl("P1(=O)([O-])OCC2OC(n)C(O)C2[O:1]1", col=C)},
    "camp_phosphodiesterase": {"H2O": whole(C), "amp": hw(u4=[C, "{P}"])},
    "guanylyl_cyclase": {"gtp": hw(d3=[C, "O"], u4=[M, "{P}–{P}"]), "PPi": whole(),
                         "cgmp": hl("P1(=O)([O-])OCC2OC(n)C(O)C2[O:1]1", col=C)},
    "cgmp_phosphodiesterase": {"H2O": whole(C), "gmp": hw(u4=[C, "{P}"])},
})


# ================================================================ plot

FLUIDITY = {
    "size": [15, 9.5],
    "x": {"min": 0, "max": 70, "label": "Temperature (°C)", "ticks": [0, 20, 40, 60]},
    "y": {"min": 0, "max": 1.15, "label": "Fluidity (disorder)", "label_pos": "center", "ticks": []},
    "curves": [
        {"id": "sat", "y": "1/(1+exp(-(x-41)*0.7))", "width": 1.2, "color": M},
        {"id": "uns", "y": "1/(1+exp(-(x-12)*0.7))", "width": 1.2, "color": C},
        {"id": "chol", "y": "1/(1+exp(-(x-41)*0.12))", "width": 1.2, "dash": True},
    ],
    "points": [{"on": "sat", "x": 41, "y": 0.5, "drop": "x", "color": M}],
    "texts": [{"at": [15, 0.62], "text": "Unsaturated (cis kinks):\nT_m low → fluid", "align": "left", "color": C, "size": 0.85},
              {"at": [46, 0.18], "text": "Saturated (16:0):\nT_m ~ 41 °C", "align": "left", "color": M, "size": 0.85},
              {"at": [50, 0.58], "text": "+ Cholesterol:\nsame T_m, broad\n(fluidity buffer)", "align": "left", "size": 0.85},
              {"at": [26, 0.42], "text": "Gel (ordered)", "size": 0.85, "color": GRAY},
              {"at": [60, 1.08], "text": "Liquid-disordered", "size": 0.85, "color": GRAY}],
}


# ================================================================ pathways

PATHWAYS = [
    {"id": "receptor_types", "title": "Receptor types and their ligands", "layout": "free", "cell_size": 6.0,
     "nodes": [
         t("mem", "Plasma membrane receptors", 1.5, 0, "heading"),
         t("gpcr", "GPCR (7TM)\nepinephrine, glucagon,\nlight, odorants", -0.6, 0.9, "protein"),
         t("rtk", "RTK\ninsulin, EGF, PDGF", 1.5, 0.9, "protein"),
         t("ser", "Receptor Ser/Thr kinase\nTGF-β", 3.6, 0.9, "protein"),
         t("rgc", "Receptor guanylyl cyclase\nANF", -0.6, 2.1, "protein"),
         t("ion", "Ligand-gated channel\nnicotinic AChR", 1.5, 2.1, "protein"),
         t("jak", "No enzyme activity (→ JAK)\ncytokines, leptin, EPO", 3.6, 2.1, "protein"),
         t("in", "Intracellular receptors", 1.5, 3.2, "heading"),
         t("sgc", "Cytosolic: soluble GC\n(NO)", -0.6, 4.0, "protein"),
         t("nr", "Nuclear (zinc finger):\nsteroids (Hsp90 in cytosol),\nvit. D, thyroid, retinoic acid", 2.3, 4.0, "protein"),
     ],
     "edges": [],
     "tags": ["signal_transduction"],
     "notes": "第一傳訊者分兩類：親水性（peptide、catecholamine）→ 膜受體 → 第二傳訊者，半衰期短；親脂性（steroid、thyroxine、vit. D、retinoic acid）→ 細胞內受體 → 調控轉錄，作用慢而久。"
              "膜受體：GPCR、RTK、受體 Ser/Thr kinase（TGF-β）、受體 guanylyl cyclase（ANF）、配體門控離子通道（nAChR）、本身無酵素活性的受體（經 JAK）；另有 adhesion receptor（integrin、cadherin、selectin、Ig-CAM）。"
              "Steroid 受體平時與 Hsp90 結合在細胞質，結合激素後入核成 homodimer；vit. D、thyroid、RA 受體已在核內，與 RXR 形成 heterodimer。"
              "訊號系統特徵：specificity、sensitivity、amplification、modularity、desensitization、integration、divergence、localized response。"},
    {"id": "gpcr_g_protein_classes", "title": "GPCR → four classes of G proteins", "layout": "free", "cell_size": 5.2,
     "nodes": [
         t("lig", "Ligand + GPCR (GEF)", 1.8, 0, "protein"),
         t("gs", "G_s (G_{olf})", -0.6, 1.2, "protein"), t("gi", "G_i (G_t, G_{gust})", 1.0, 1.2, "protein"),
         t("gq", "G_{q/11}", 2.6, 1.2, "protein"), t("g12", "G_{12/13}", 4.2, 1.2, "protein"),
         t("ac", "Adenylyl cyclase", -0.3, 2.4, "protein"), t("camp", "cAMP → PKA", -0.3, 3.5, "name"),
         t("plc", "PLC-β", 2.6, 2.4, "protein"), t("ip3", "IP_3 → Ca^{2+}\nDAG → PKC", 2.6, 3.5, "name"),
         t("rho", "RhoGEF → Rho", 4.2, 2.4, "protein"),
         t("ctx", "Cholera toxin\n(Gsα always on)", -1.9, 1.2, "drug"), t("ptx", "Pertussis toxin\n(Gi blocked)", 1.55, 2.4, "drug"),
     ],
     "edges": [
         {"from": "lig", "to": "gs"}, {"from": "lig", "to": "gi"}, {"from": "lig", "to": "gq"}, {"from": "lig", "to": "g12"},
         {"from": "gs", "to": "ac", "label": "+", "color": GREEN}, {"from": "gi", "to": "ac", "inhibit": True},
         {"from": "ac", "to": "camp", "step": "adenylyl_cyclase", "label": "ATP → cAMP"},
         {"from": "gq", "to": "plc"}, {"from": "plc", "to": "ip3", "label": "PIP_2"}, {"from": "g12", "to": "rho"},
         {"from": "ctx", "to": "gs", "color": RED}, {"from": "ptx", "to": "gi", "inhibit": True},
     ],
     "tags": ["signal_transduction", "gpcr"],
     "notes": "GPCR = 7 次跨膜（舊名 serpentine receptor），本身是 G 蛋白的 GEF：促使 Gα 釋出 GDP、結合 GTP。"
              "四類：Gs（β-adrenergic、glucagon；嗅覺 Golf）活化 AC → cAMP ↑；Gi（α₂、M₂、M₄；視覺 transducin Gt 與苦味 gustducin 經 PDE 使 cyclic nucleotide 下降）抑制 AC；"
              "Gq/11（α₁、M₁/₃/₅、angiotensin II、vasopressin V₁）活化 PLC-β → IP₃ + DAG；G12/13 → RhoGEF → Rho（細胞骨架）。"
              "毒素：cholera toxin ADP-ribosylate Gsα 的 Arg → GTPase 失效 → cAMP 持續；pertussis toxin ADP-ribosylate Giα 的 Cys → 無法抑制 AC → cAMP 上升。"},
    {"id": "camp_pka_cascade", "title": "Epinephrine → cAMP → PKA cascade (amplification)", "layout": "free", "cell_size": 4.4,
     "nodes": [
         t("epi", "Epinephrine / glucagon", 0, 0, "name"), t("rec", "β-Adrenergic receptor (GPCR)", 0, 0.9, "protein"),
         t("gs", "G_sα–GTP", 0, 1.8, "protein"), t("ac", "Adenylyl cyclase", 0, 2.7, "protein"),
         n("camp", "camp", 0, 3.6), t("pka", "PKA: R_2C_2 + 4 cAMP\n→ 2 C (active)", 0, 4.6, "protein"),
         t("phk", "Phosphorylase kinase", -1.8, 5.7, "protein"), t("gp", "Glycogen phosphorylase\nb → a", -1.8, 6.7, "protein"),
         t("glc", "Glycogen → G1P\n(10^4–10^6 ×)", -1.8, 7.7, "name"),
         t("creb", "CREB-P + CBP/p300\n→ gene expression\n(e.g. PEPCK)", 1.8, 5.8, "name"),
         n("amp", "amp", 2.6, 3.6),
         t("gs2", "Glycogen synthase → inactive", 1.8, 7.0, "name"),
     ],
     "edges": [
         {"from": "epi", "to": "rec"}, {"from": "rec", "to": "gs", "label": "GDP → GTP"}, {"from": "gs", "to": "ac"},
         {"from": "ac", "to": "camp", "step": "adenylyl_cyclase", "label": "ATP → PP_i"},
         {"from": "camp", "to": "pka"}, {"from": "pka", "to": "phk"}, {"from": "phk", "to": "gp"}, {"from": "gp", "to": "glc"},
         {"from": "pka", "to": "creb", "label": "nucleus"}, {"from": "pka", "to": "gs2"},
         {"from": "camp", "to": "amp", "step": "camp_phosphodiesterase", "label": "PDE\n(+ H_2O)"},
     ],
     "tags": ["signal_transduction", "gpcr", "second_messenger"],
     "notes": "每一層都是酵素 → 訊號放大（enzyme cascade）：一個 epinephrine 可使數百萬個 glucose 釋出。"
              "PKA = R₂C₂；每個 R 結合 2 個 cAMP（共 4 個）→ 釋出 2 個 C；AKAP 把 PKA 錨定在特定位置（localized response）。"
              "PKA 兩條路：① 直接磷酸化代謝酵素（phosphorylase kinase、glycogen synthase、PFK-2/FBPase-2、HSL、ACC）② 入核磷酸化 CREB，結合 CBP/p300 → 轉錄。"
              "終止：Gα 水解 GTP、PDE 分解 cAMP、phosphatase、受體被 GRK 磷酸化後由 β-arrestin 去敏感化。"},
    {"id": "no_cgmp_relaxation", "title": "Acetylcholine → NO → cGMP: vascular relaxation", "layout": "free", "cell_size": 4.4,
     "nodes": [
         t("endo", "Endothelial cell", -1.2, -0.7, "heading"), t("smc", "Smooth muscle cell", 2.6, -0.7, "heading"),
         t("ach", "Acetylcholine → M_3 (GPCR, G_q)", 0, 0, "protein"), t("plc", "PLC-β: PIP_2 → IP_3 + DAG", 0, 1.0, "protein"),
         t("ca", "IP_3 receptor (ER)\n→ Ca^{2+} ↑ → Ca^{2+}–calmodulin", 0, 2.0, "protein"),
         t("nos", "eNOS", -0.3, 3.0, "protein"), n("no", "nitric_oxide", 1.3, 3.0),
         t("sgc", "Soluble guanylyl cyclase\n(heme)", 3.5, 3.0, "protein"),
         n("cgmp", "cgmp", 3.5, 4.1), t("pkg", "PKG → Ca^{2+} ↓\n→ relaxation", 3.5, 5.2, "name"),
         n("gmp", "gmp", 5.4, 5.2),
         t("vi", "Sildenafil", 5.6, 4.1, "drug"), t("ntg", "Nitroglycerin\n(NO donor)", 5.4, 2.0, "drug"),
         t("arg", "Arg → citrulline", -0.3, 3.75, "name"),
     ],
     "edges": [
         {"from": "ach", "to": "plc"}, {"from": "plc", "to": "ca"}, {"from": "ca", "to": "nos", "label": "+"},
         {"from": "nos", "to": "no", "step": "nitric_oxide_synthase", "label": "O_2, NADPH"},
         {"from": "no", "to": "sgc", "label": "diffuses"},
         {"from": "sgc", "to": "cgmp", "step": "guanylyl_cyclase", "label": "GTP → PP_i"},
         {"from": "cgmp", "to": "pkg"}, {"from": "cgmp", "to": "gmp", "step": "cgmp_phosphodiesterase", "label": "PDE5"},
         {"from": "vi", "to": "edge:cgmp->gmp", "inhibit": True}, {"from": "ntg", "to": "sgc"},
     ],
     "regions": [{"kind": "membrane", "x": 2.15, "y0": -1.0, "y1": 5.8, "label": ""}],
     "tags": ["signal_transduction", "second_messenger"],
     "notes": "內皮細胞：ACh → M₃（Gq）→ PLC-β → IP₃ → ER 釋出 Ca²⁺ → Ca²⁺–calmodulin 活化 eNOS（Arg + O₂ + NADPH → citrulline + NO）。"
              "NO 為氣體，擴散進平滑肌，結合可溶型 guanylyl cyclase 的 heme → cGMP ↑ → PKG → 細胞內 Ca²⁺ 下降、myosin light chain 去磷酸 → 放鬆、血管擴張。"
              "PDE5 把 cGMP → GMP 終止訊號；sildenafil 抑制 PDE5。Nitroglycerin、nitroprusside 為 NO donor。"
              "氣體傳訊分子（gasotransmitters）：NO（NOS）、CO（heme oxygenase）、H₂S（cystathionine γ-lyase）。"},
    {"id": "insulin_rtk_signaling", "title": "Insulin receptor (RTK): MAPK, PI3K–Akt and PLCγ branches", "layout": "free", "cell_size": 4.4,
     "nodes": [
         t("ins", "Insulin / growth factor", 1.6, 0, "name"),
         t("ir", "Receptor Tyr kinase\n(trans-autophosphorylation, Tyr-P)", 1.6, 0.95, "protein"),
         t("irs", "IRS-1 (Tyr-P)", 1.6, 1.95, "protein"),
         t("grb", "Grb2 (SH2 → P-Tyr;\nSH3 → Pro of Sos)", -0.4, 2.95, "protein"),
         t("sos", "Sos (GEF)", -0.4, 3.9, "protein"), t("ras", "Ras–GDP → Ras–GTP", -0.4, 4.8, "protein"),
         t("raf", "Raf-1 (MAPKKK, Ser/Thr)", -0.4, 5.7, "protein"), t("mek", "MEK (MAPKK, dual: Thr + Tyr)", -0.4, 6.6, "protein"),
         t("erk", "ERK (MAPK)", -0.4, 7.5, "protein"), t("elk", "Nucleus: Elk-1 → growth genes", -0.4, 8.4, "name"),
         t("pi3k", "PI3K: PIP_2 → PIP_3", 2.2, 2.95, "protein"), t("pdk", "PDK1", 2.2, 3.9, "protein"),
         t("akt", "PKB / Akt", 2.2, 4.8, "protein"),
         t("eff", "GLUT4 translocation,\nGSK3 inhibited → glycogen synthesis,\nPDE3B → cAMP ↓, mTOR", 2.9, 5.95, "name"),
         t("plcg", "PLC-γ: PIP_2 → IP_3 + DAG", 4.3, 1.95, "protein"),
     ],
     "edges": [
         {"from": "ins", "to": "ir"}, {"from": "ir", "to": "irs", "label": "Tyr-P"}, {"from": "irs", "to": "grb"},
         {"from": "grb", "to": "sos"}, {"from": "sos", "to": "ras", "label": "GDP → GTP"}, {"from": "ras", "to": "raf"},
         {"from": "raf", "to": "mek", "label": "2 Ser-P"}, {"from": "mek", "to": "erk", "label": "Thr-P, Tyr-P"},
         {"from": "erk", "to": "elk"},
         {"from": "irs", "to": "pi3k"}, {"from": "pi3k", "to": "pdk"}, {"from": "pdk", "to": "akt", "label": "P"},
         {"from": "akt", "to": "eff"}, {"from": "ir", "to": "plcg"},
     ],
     "tags": ["signal_transduction", "rtk"],
     "notes": "Insulin receptor 為 α₂β₂ 四聚體（已是二聚體）；生長因子受體（EGFR、PDGFR）要配體誘導二聚化才能 trans-autophosphorylation。"
              "MAPK 路徑：IRS-1-P —SH2-Grb2-SH3— Pro-Sos → Ras（GEF 換 GTP）→ Raf-1 → MEK（磷酸化 ERK 的 Thr 與 Tyr，雙特異性 kinase）→ ERK 入核 → Elk-1。IRS-1、Grb2、Sos 本身沒有酵素活性（adaptor）。"
              "PI3K–Akt：PIP₃ 招募 PDK1 與 Akt 到膜 → PDK1 磷酸化 Akt → GLUT4、肝醣合成（抑制 GSK3）。PTEN 去 PIP₃ 的 3-磷酸。"
              "PLC-γ（由 RTK 活化；GPCR 活化的是 PLC-β）。Ras 突變（GTPase 失效）常見於癌症。"},
    {"id": "jak_stat_pathway", "title": "JAK–STAT pathway", "layout": "free", "cell_size": 4.4,
     "nodes": [
         t("lig", "Cytokine, leptin, EPO, GH", 0, 0, "name"),
         t("rec", "Receptor dimerizes\n(no intrinsic kinase)", 0, 1.0, "protein"),
         t("jak", "JAK (Janus kinase)\ntrans-autophosphorylation", 0, 2.0, "protein"),
         t("ty", "Receptor Tyr-P", 0, 3.0, "protein"), t("stat", "STAT binds via SH2", 0, 4.0, "protein"),
         t("statp", "STAT-P (by JAK)\n→ dimer", 0, 5.0, "protein"), t("nuc", "Nucleus: transcription factor", 0, 6.0, "name"),
     ],
     "edges": [{"from": "lig", "to": "rec"}, {"from": "rec", "to": "jak"}, {"from": "jak", "to": "ty", "label": "Tyr-P"},
               {"from": "ty", "to": "stat"}, {"from": "stat", "to": "statp", "label": "Tyr-P"}, {"from": "statp", "to": "nuc"}],
     "tags": ["signal_transduction", "rtk"],
     "notes": "受體本身無酵素活性，配體使受體二聚化 → 結合的 JAK 互相磷酸化 → 磷酸化受體的 Tyr → STAT 以 SH2 結合 P-Tyr → JAK 磷酸化 STAT → STAT-P 脫離、二聚化、入核當轉錄因子。"
              "例：EPO（紅血球生成）、leptin（飽足）、GH、interferon、IL。JAK2 V617F → polycythemia vera。SOCS 回饋抑制。"},
    {"id": "nfkb_pathway", "title": "NF-κB pathway", "layout": "free", "cell_size": 4.4,
     "nodes": [
         t("sig", "TNF-α → TNF receptor\nLPS → Toll-like receptor\nUV, virus", 0, 0, "name"),
         t("ikk", "IκB kinase (IKK)", 0, 1.15, "protein"),
         t("cx", "NF-κB · IκB (cytosol, inactive)", 0, 2.15, "protein"),
         t("ib", "IκB-P → ubiquitin\n→ proteasome", -1.6, 3.2, "name"),
         t("nf", "Free NF-κB → nucleus", 1.6, 3.2, "protein"),
         t("gene", "Inflammatory cytokines,\nanti-apoptotic genes", 1.6, 4.2, "name"),
     ],
     "edges": [{"from": "sig", "to": "ikk"}, {"from": "ikk", "to": "cx", "label": "IκB-P"},
               {"from": "cx", "to": "ib"}, {"from": "cx", "to": "nf"}, {"from": "nf", "to": "gene"}],
     "tags": ["signal_transduction"],
     "notes": "NF-κB 平時被 IκB 留在細胞質。TNF-α、TLR（細菌 LPS、病毒）、UV → IKK 磷酸化 IκB → 泛素化 → proteasome 分解 → NF-κB 露出 NLS 入核 → 發炎細胞激素、抗凋亡基因。"
              "Glucocorticoid 誘導 IκB 合成 → 抗發炎；aspirin（高劑量）抑制 IKK。"},
    {"id": "apoptosis_pathways", "title": "Apoptosis: extrinsic and intrinsic pathways", "layout": "free", "cell_size": 4.4,
     "nodes": [
         t("ext", "Extrinsic", -1.4, -0.65, "heading"), t("int", "Intrinsic (mitochondrial)", 1.8, -0.65, "heading"),
         t("fas", "FasL, TNF → death receptor", -1.4, 0.2, "protein"), t("c8", "Caspase 8", -1.4, 1.4, "protein"),
         t("ros", "DNA damage, ROS,\ngrowth-factor loss", 1.8, 0.2, "name"),
         t("bax", "Bax / Bak pores\n(outer mito. membrane)", 1.8, 1.4, "protein"), t("bcl", "Bcl-2", 3.5, 1.4, "drug"),
         t("cyt", "Cytochrome c release\n(cardiolipin peroxidation)", 1.8, 2.5, "carrier"),
         t("apo", "Apoptosome\n(Apaf-1 · cyt c · dATP)_7", 1.8, 3.6, "protein"), t("c9", "Caspase 9", 1.8, 4.6, "protein"),
         t("c3", "Caspase 3, 7 (executioner)", 0.2, 5.6, "protein"),
         t("death", "Proteolysis, DNA laddering, PS exposed,\ncell shrinks → phagocytosed (no inflammation)", 0.2, 6.6, "name"),
     ],
     "edges": [
         {"from": "fas", "to": "c8"}, {"from": "c8", "to": "c3"}, {"from": "c8", "to": "bax", "label": "Bid"},
         {"from": "ros", "to": "bax"}, {"from": "bcl", "to": "bax", "inhibit": True}, {"from": "bax", "to": "cyt"},
         {"from": "cyt", "to": "apo"}, {"from": "apo", "to": "c9"}, {"from": "c9", "to": "c3"}, {"from": "c3", "to": "death"},
     ],
     "tags": ["signal_transduction", "apoptosis"],
     "notes": "Caspase：活性位 Cys，切在 Asp 之後。Initiator：2、8、9、10；executioner：3、6、7；發炎型：1、4、5、12。"
              "外源：FasL／TNF → death receptor → caspase 8（可經 Bid 放大到粒線體）。內源：DNA damage（p53 → Bax）、ROS → Bax/Bak 在外膜形成孔，cardiolipin 過氧化使 cyt c 釋出 → 與 Apaf-1、dATP 組成 apoptosome → caspase 9 → caspase 3/7。Bcl-2 抗凋亡。"
              "C. elegans 對應：Ced-9 = Bcl-2（抗），Ced-4 = Apaf-1，Ced-3 = caspase 9。"
              "凋亡 vs 壞死：凋亡細胞縮小、DNA ladder、PS 外翻（scramblase，Ca²⁺）、無發炎；壞死細胞腫脹、破裂、發炎。"},
    {"id": "gtpase_switch", "title": "GTPase switch: GEF and GAP", "layout": "cycle", "radius": 1.7, "cell_size": 5.0,
     "start_angle": 180, "center_text": "GTPase\nswitch",
     "nodes": [
         {"id": "off", "text": "GTPase–GDP\n(inactive)", "kind": "protein"},
         {"id": "on", "text": "GTPase–GTP\n(active) → effectors", "kind": "protein"},
     ],
     "edges": [{"from": "off", "to": "on", "label": "GEF\n(GDP out, GTP in)", "bend": 0.5},
               {"from": "on", "to": "off", "label": "GAP / RGS\n(GTP → GDP + P_i)", "bend": 0.5, "color": RED}],
     "tags": ["signal_transduction", "gpcr"],
     "notes": "GTPase 本身水解 GTP 很慢；GEF 促進 GDP 釋出（活化），GAP 加速 GTP 水解（關閉）。"
              "例：G 蛋白（GEF = GPCR；GAP = RGS）、Ras（GEF = Sos；GAP = RasGAP/NF1）、EF-Tu（GEF = EF-Ts）、Ran、Rho、Rab、Arf。"
              "Ras 致癌突變（G12V）使 GTPase 失效 → 持續 on；NF1 缺失 → GAP 失去 → neurofibromatosis。"
              "脂質修飾：Gα（myristoyl、palmitoyl）、Gγ（farnesyl/geranylgeranyl）、Ras（farnesyl）把它們錨在膜內側。"},
    {"id": "p53_rb_checkpoint", "title": "DNA damage → p53 → G1 arrest", "layout": "free", "cell_size": 4.4,
     "nodes": [
         t("dmg", "DNA double-strand break", 0, 0, "name"), t("mrn", "MRN complex", 0, 0.9, "protein"),
         t("atm", "ATM / ATR (Ser/Thr kinases)", 0, 1.8, "protein"), t("p53", "p53 (stabilized)", 0, 2.7, "protein"),
         t("p21", "p21 (CKI)", 0, 3.6, "protein"), t("cdk", "CDK2 · cyclin E", 0, 4.6, "protein"),
         t("rb", "Rb · E2F (E2F held)", 0, 5.6, "protein"), t("g1", "G1 arrest → repair\n(or apoptosis via Bax)", -2.4, 2.7, "name"),
         t("rbp", "If no p21: Rb-P releases E2F\n→ S-phase genes", 2.9, 5.6, "name"),
     ],
     "edges": [{"from": "dmg", "to": "mrn"}, {"from": "mrn", "to": "atm"}, {"from": "atm", "to": "p53", "label": "P"},
               {"from": "p53", "to": "p21", "label": "transcription"}, {"from": "p21", "to": "cdk", "inhibit": True},
               {"from": "cdk", "to": "rb", "inhibit": True, "label": "cannot phosphorylate Rb"},
               {"from": "p53", "to": "g1"}, {"from": "rb", "to": "rbp", "color": GRAY}],
     "tags": ["signal_transduction", "cancer"],
     "notes": "DNA 雙股斷裂 → MRN → ATM/ATR 磷酸化 p53（不再被 MDM2 送去分解）→ p53 轉錄 p21（CDK inhibitor）→ CDK2/cyclin E 失活 → Rb 保持去磷酸、抓住 E2F → 停在 G1 修復；無法修復則 p53 → Bax → 凋亡。"
              "p53、Rb、BRCA1/2 為腫瘤抑制基因（loss-of-function、隱性、常為生殖細胞遺傳）；Ras、Myc、ErbB2 為原致癌基因（gain-of-function、顯性、體細胞突變）。"
              "原致癌基因活化方式：點突變、基因擴增、染色體轉位、promoter/enhancer 插入。c-onc 與 v-onc（病毒來源）。"},
    {"id": "two_component_system", "title": "Two-component system (bacteria, plants)", "layout": "linear",
     "nodes": [
         {"id": "sig", "text": "Signal", "kind": "name"},
         {"id": "sk", "text": "Sensor kinase\nHis → His-P\n(autophosphorylation)", "kind": "protein"},
         {"id": "rr", "text": "Response regulator\nAsp → Asp-P", "kind": "protein"},
         {"id": "out", "text": "Transcription /\nflagellar motor", "kind": "name"},
     ],
     "edges": [{"from": "sig", "to": "sk"}, {"from": "sk", "to": "rr", "label": "phosphorelay"}, {"from": "rr", "to": "out"}],
     "tags": ["signal_transduction"],
     "notes": "細菌與植物的二元系統：感應器 His kinase 自我磷酸化（His-P，phosphoramidate），再把磷酸轉給 response regulator 的 Asp（acyl-phosphate）→ 調控轉錄或鞭毛運動（趨化）。動物主要用 Ser/Thr、Tyr kinase。"},
    {"id": "bacterial_toxins", "title": "ADP-ribosylating bacterial toxins", "layout": "free", "cell_size": 5.6,
     "nodes": [
         t("ct", "Cholera toxin", 0, 0, "drug"), t("pt", "Pertussis toxin", 0, 1.3, "drug"), t("dt", "Diphtheria toxin", 0, 2.6, "drug"),
         t("gs", "G_sα (Arg)\nGTPase lost → always on", 2.2, 0, "protein"), t("gi", "G_iα (Cys)\ncannot inhibit AC", 2.2, 1.3, "protein"),
         t("ef", "eEF2 (diphthamide,\nmodified His)", 2.2, 2.6, "protein"),
         t("r1", "cAMP ↑ → CFTR Cl^- secretion\n→ watery diarrhea", 4.6, 0, "name"),
         t("r2", "cAMP ↑ (whooping cough)", 4.6, 1.3, "name"), t("r3", "Translation stops\n→ cell death", 4.6, 2.6, "name"),
     ],
     "edges": [{"from": "ct", "to": "gs", "label": "ADP-ribose\n(from NAD^+)"}, {"from": "pt", "to": "gi"}, {"from": "dt", "to": "ef"},
               {"from": "gs", "to": "r1"}, {"from": "gi", "to": "r2"}, {"from": "ef", "to": "r3"}],
     "tags": ["signal_transduction", "gpcr"],
     "notes": "三種毒素都是 ADP-ribosyltransferase：把 NAD⁺ 的 ADP-ribose 接到目標蛋白（釋出 nicotinamide）。"
              "Cholera：Gsα 的 Arg → GTPase 失活 → AC 持續 → 腸道 cAMP ↑ → CFTR 分泌 Cl⁻ → 水瀉。"
              "Pertussis：Giα 的 Cys → Gi 無法被受體活化 → 無法抑制 AC → cAMP ↑。"
              "Diphtheria（Pseudomonas exotoxin A 同）：eEF2 的 diphthamide（修飾的 His）→ 轉位停止 → 蛋白質合成停止。"},
    {"id": "second_messengers", "title": "Second messengers: source → target", "layout": "free", "cell_size": 5.0,
     "nodes": [
         t("h1", "Source", 0, -0.7, "heading"), t("h2", "Messenger", 2.4, -0.7, "heading"), t("h3", "Target", 4.6, -0.7, "heading"),
         t("s1", "Adenylyl cyclase", 0, 0, "protein"), n("m1", "camp", 2.4, 0), t("e1", "PKA, Epac", 4.6, 0, "protein"),
         t("s2", "Guanylyl cyclase", 0, 0.8, "protein"), n("m2", "cgmp", 2.4, 0.8), t("e2", "PKG, CNG channel", 4.6, 0.8, "protein"),
         t("s3", "PLC on PIP_2", 0, 1.6, "protein"), t("m3", "IP_3", 2.4, 1.6, "name"), t("e3", "IP_3 receptor (ER)", 4.6, 1.6, "protein"),
         t("m4", "DAG", 2.4, 2.4, "name"), t("e4", "PKC", 4.6, 2.4, "protein"),
         t("s5", "ER / channels", 0, 3.2, "protein"), t("m5", "Ca^{2+}", 2.4, 3.2, "ion"), t("e5", "Calmodulin (4 Ca^{2+}, EF hand),\nPKC, CaM kinase", 4.6, 3.2, "protein"),
         t("s6", "NO synthase", 0, 4.0, "protein"), n("m6", "nitric_oxide", 2.4, 4.0), t("e6", "Soluble guanylyl cyclase", 4.6, 4.0, "protein"),
         t("s7", "Sphingomyelinase", 0, 4.8, "protein"), n("m7", "ceramide", 2.4, 4.8), t("e7", "Protein kinases / phosphatases", 4.6, 4.8, "protein"),
     ],
     "edges": [
         {"from": "s1", "to": "m1", "step": "adenylyl_cyclase"}, {"from": "m1", "to": "e1"},
         {"from": "s2", "to": "m2", "step": "guanylyl_cyclase"}, {"from": "m2", "to": "e2"},
         {"from": "s3", "to": "m3"}, {"from": "s3", "to": "m4"}, {"from": "m3", "to": "e3"}, {"from": "m4", "to": "e4"},
         {"from": "e3", "to": "m5", "color": GRAY}, {"from": "s5", "to": "m5"}, {"from": "m5", "to": "e5"},
         {"from": "s6", "to": "m6", "step": "nitric_oxide_synthase"}, {"from": "m6", "to": "e6"},
         {"from": "s7", "to": "m7", "step": "sphingomyelinase"}, {"from": "m7", "to": "e7"},
     ],
     "tags": ["signal_transduction", "second_messenger"],
     "notes": "第二傳訊者：cAMP、cGMP、IP₃、DAG、Ca²⁺、NO（可兼第一與第二傳訊者）、ceramide、cADP-ribose；GTP、ATP、insulin 不是。"
              "Ser/Thr kinase 家族與活化物：PKA（cAMP）、PKB/Akt（PIP₃）、PKC（DAG、PS、Ca²⁺；phorbol ester TPA 活化 → 腫瘤促進劑）、PKG（cGMP）、CaMK（Ca²⁺–calmodulin）。"
              "PKC 三類：conventional（PS、DAG、Ca²⁺）、novel（PS、DAG）、atypical（都不需）。"
              "Calmodulin：4 個 Ca²⁺ 結合 EF hand（helix-loop-helix）；也是 phosphorylase kinase 的 δ 次單元。"},
    # ------------------------------------------------------------ membrane
    {"id": "membrane_transporters", "title": "Membrane transport proteins", "layout": "free", "cell_size": 5.0,
     "nodes": [
         t("out", "Outside", -3.3, -1.4, "heading"), t("in", "Cytosol", -3.3, 1.4, "heading"),
         t("d1", "Simple diffusion:\nO_2, CO_2, NO, ethanol", -1.0, -1.4, "ion"), t("d2", "O_2, CO_2 …", -1.0, 1.4, "ion"),
         t("glut", "GLUT\n(uniport, facilitated)", 0.9, 0, "protein", h=3.2),
         t("sglt", "SGLT1\n(symport,\n2° active)", 2.6, 0, "protein", h=3.2),
         t("nak", "Na^+/K^+-ATPase\n(antiport, 1° active,\nP-type)", 4.4, 0, "protein", h=3.2),
         t("ch", "Channel\n(K^+, AQP)", 6.1, 0, "protein", h=3.2),
         t("abc", "ABC transporter\n(MDR1, CFTR)", 7.8, 0, "protein", h=3.2),
         t("g1", "Glucose", 0.9, -1.4, "ion"), t("g2", "Glucose", 0.9, 1.4, "ion"),
         t("s1", "2 Na^+ + glucose", 2.6, -1.4, "ion"), t("s2", "2 Na^+ + glucose", 2.6, 1.4, "ion"),
         t("n1", "3 Na^+ out", 4.0, -1.4, "ion"), t("n2", "3 Na^+", 4.0, 1.4, "ion"),
         t("k1", "2 K^+", 4.9, -1.4, "ion"), t("k2", "2 K^+ in", 4.9, 1.4, "ion"),
         t("atp", "ATP → ADP + P_i", 4.4, 2.2, "ion"),
         t("c1", "K^+, H_2O", 6.1, 1.4, "ion"), t("c2", "K^+, H_2O", 6.1, -1.4, "ion"),
         t("a1", "Drug", 7.8, 1.4, "ion"), t("a2", "Drug (out)", 7.8, -1.4, "ion"),
     ],
     "edges": [
         {"from": "d1", "to": "d2", "color": GRAY}, {"from": "g1", "to": "g2"}, {"from": "s1", "to": "s2", "color": M}, {"from": "n2", "to": "n1", "color": M},
         {"from": "k1", "to": "k2", "color": C}, {"from": "c1", "to": "c2"}, {"from": "a1", "to": "a2", "color": M},
     ],
     "regions": [{"kind": "membrane", "y": 0, "x0": -3.6, "x1": 8.9, "label": ""}],
     "tags": ["membrane", "membrane_transport"],
     "notes": "Simple diffusion：氣體、乙醇、尿素（H₂O 略可）；離子、胺基酸、glucose、ATP、G6P 不能自由通過。"
              "Transporter（carrier）：GLUT（uniport，促進性擴散，不耗能）、SGLT（Na⁺ 同向 symport，二級主動運輸，靠 Na⁺ 梯度）、ABC（ATP，MDR1 = P-glycoprotein 把藥物打出；CFTR 為 Cl⁻ 通道；ABCA1 膽固醇外送、ABCD1 VLCFA 入 peroxisome）。"
              "Pump：Na⁺/K⁺-ATPase（3 Na⁺ 出、2 K⁺ 進，P-type：Asp 磷酸化中間物，ouabain／digoxin 抑制）。ATPase 類型：P（Na/K、Ca²⁺、flippase）、V（溶體、內體 H⁺）、F（F₀F₁ ATP synthase）、ABC（floppase、MDR）。"
              "Channel：單一閘門，voltage-gated 或 ligand-gated（nAChR），aquaporin；pore（porin）一直開。"
              "Ionophore：DNP（H⁺）、valinomycin（K⁺ mobile carrier）、gramicidin（channel former）。"},
    {"id": "membrane_asymmetry", "title": "Lipid asymmetry of the plasma membrane", "layout": "free", "cell_size": 4.0,
     "nodes": [
         t("out", "Outer leaflet: PC, sphingomyelin,\nglycolipids, GPI-anchored proteins", 1.8, -1.2, "name"),
         t("in", "Inner leaflet: PE, PS, PI (PIP_2);\nmyristoyl / palmitoyl / farnesyl-anchored proteins", 1.8, 1.2, "name"),
         t("flip", "Flippase\n(P-type, ATP)", -0.6, 0, "protein", h=3.0),
         t("flop", "Floppase\n(ABC, ATP)", 1.8, 0, "protein", h=3.0),
         t("scr", "Scramblase\n(Ca^{2+}, no ATP)", 4.2, 0, "protein", h=3.0),
         t("f1", "PS, PE", -1.4, -0.75, "ion"), t("f2", "PS, PE", -1.4, 0.75, "ion"),
         t("o1", "PL", 1.0, 0.75, "ion"), t("o2", "PL", 1.0, -0.75, "ion"),
         t("x1", "PS", 5.6, 0.75, "ion"), t("x2", "PS (apoptosis,\nplatelets)", 5.6, -0.8, "ion"),
     ],
     "edges": [{"from": "f1", "to": "f2", "color": M}, {"from": "o1", "to": "o2", "color": C},
               {"from": "x1", "to": "x2", "color": GRAY, "reversible": True}],
     "regions": [{"kind": "membrane", "y": 0, "x0": -2.4, "x1": 6.6, "label": ""}],
     "tags": ["membrane", "membrane_lipid"],
     "notes": "外層：PC、SM、醣脂、GPI-anchored 蛋白；內層：PE、PS、PI（PIP₂）、脂質錨定（myristoyl、palmitoyl、farnesyl）蛋白；膽固醇兩層約相等。"
              "Flip-flop 自發極慢；flippase（P-type ATPase）把 PS、PE 由外送入內；floppase（ABC）由內送外；scramblase（Ca²⁺ 活化、不耗 ATP）雙向打亂 → 凋亡細胞與活化血小板 PS 外露（吞噬訊號、凝血表面）。"
              "側向擴散快（~μm/s）、轉動快；膜蛋白分類：integral（transmembrane）、peripheral、lipid-anchored（Lehninger 另有 amphitropic）。"
              "Lipid raft：富含膽固醇、sphingolipid、GPI 蛋白；caveolin 形成 caveolae。分離膜蛋白：SDS（離子型）、Triton X-100（非離子型，保留活性）。"},
    {"id": "membrane_fluidity", "title": "Membrane fluidity and phase transition", "layout": "plot", "plot": FLUIDITY,
     "tags": ["membrane", "membrane_lipid"],
     "notes": "加熱使膜由有序的 gel 相轉成流動的 liquid-disordered 相，轉換溫度 T_m。"
              "飽和、長鏈脂肪酸 → T_m 高；不飽和（cis 彎折）、短鏈 → T_m 低、較流動。"
              "膽固醇是流動性緩衝劑（fluidity buffer）：T_m 不變但轉換變寬；高溫時限制運動（降低流動性），低溫時阻止鏈緊密排列（維持流動性）。"
              "冷環境下動物細胞靠增加不飽和脂肪酸與膽固醇維持流動性（homeoviscous adaptation）。"},
    {"id": "ldl_receptor_endocytosis", "title": "Receptor-mediated endocytosis of LDL", "layout": "cycle", "radius": 2.2,
     "cell_size": 5.0, "start_angle": -90, "center_text": "LDL\nreceptor\ncycle",
     "nodes": [
         {"id": "bind", "text": "LDL (apoB-100) binds\nLDL receptor", "kind": "protein"},
         {"id": "pit", "text": "Clathrin-coated pit\n(AP2 adaptor)", "kind": "protein"},
         {"id": "ves", "text": "Coated vesicle\n(dynamin pinches off)", "kind": "name"},
         {"id": "endo", "text": "Uncoated → endosome\n(low pH: LDL released)", "kind": "name"},
         {"id": "lys", "text": "Lysosome: cholesterol esters\n→ free cholesterol", "kind": "name"},
         {"id": "rec", "text": "Receptor recycles\nto surface", "kind": "protein"},
     ],
     "edges": [{"from": "bind", "to": "pit"}, {"from": "pit", "to": "ves"}, {"from": "ves", "to": "endo"},
               {"from": "endo", "to": "lys"}, {"from": "endo", "to": "rec", "bend": 0.3}, {"from": "rec", "to": "bind"}],
     "tags": ["membrane", "membrane_transport", "lipoprotein"],
     "notes": "Clathrin 為三腳形 triskelion（3 重鏈 + 3 輕鏈）；AP2 連接受體與 clathrin；dynamin（GTPase）切下小泡。"
              "內體酸化使 LDL 與受體分離，受體回到膜上；LDL 進溶體，cholesteryl ester 被 acid lipase 水解。"
              "游離膽固醇：抑制 HMG-CoA reductase（SREBP 停止）、減少 LDL receptor 合成、活化 ACAT 酯化儲存。"
              "家族性高膽固醇血症：LDL receptor 缺陷。胞吞種類：phagocytosis（>200 nm）、pinocytosis（clathrin、caveolin、macropinocytosis）。"},
]


# ================================================================ illustrated versions (variant `illus`)

def _lbl(x, text, y=14.8):
    return {"g": "text", "at": [x, y], "text": text, "size": 0.85}


def _io(x, y0, y1, text, hue, col, ion_at_end=True):
    """ion / solute crossing: sphere at both ends + arrow."""
    d = 1 if y1 > y0 else -1
    return [{"g": "ion", "at": [x, y0], "text": text, "hue": hue},
            {"g": "arrow", "from": [x, y0 + d * 0.65], "to": [x, y1 - d * 0.65], "color": col},
            {"g": "ion", "at": [x, y1], "text": text, "hue": hue}]


TRANSPORT_SCENE = {"items": [
    {"g": "text", "at": [-1.2, 1.4], "text": "Outside", "align": "left", "size": 0.9, "color": "gray"},
    {"g": "text", "at": [-1.2, 12.6], "text": "Cytosol", "align": "left", "size": 0.9, "color": "gray"},
    {"g": "membrane", "y": 7, "x0": -1.2, "x1": 44},
    *_io(2, 3.4, 10.6, "O_2", "gray", "gray"), _lbl(2, "Simple\ndiffusion"),
    {"id": "glut", "g": "transporter", "at": [8, 7], "hue": "blue", "open": "out"},
    *_io(8, 3.0, 11.0, "Glc", "orange", "black"), _lbl(8, "GLUT\n(uniport,\nfacilitated)"),
    {"id": "sglt", "g": "transporter", "at": [15, 7], "hue": "teal", "open": "out"},
    *_io(14.2, 3.0, 11.0, "Na^+", "purple", "magenta"), *_io(15.9, 3.0, 11.0, "Glc", "orange", "black"),
    _lbl(15, "SGLT1\n(symport,\n2° active)"),
    {"id": "pump", "g": "pump", "at": [23, 7], "hue": "orange", "accent": True},
    *_io(22.0, 11.0, 3.0, "Na^+", "purple", "magenta"), {"g": "text", "at": [20.6, 3.0], "text": "3×", "size": 0.8},
    *_io(24.2, 3.0, 11.0, "K^+", "green", "cyan"), {"g": "text", "at": [25.6, 3.0], "text": "2×", "size": 0.8},
    {"g": "text", "at": [26.4, 12.7], "text": "ATP → ADP + P_i", "font": "serif", "size": 0.8},
    _lbl(23, "Na^+/K^+-ATPase\n(antiport, 1° active,\nP-type)", 15.3),
    {"id": "channel", "g": "channel", "at": [31.5, 7], "hue": "teal", "open": True},
    {"g": "ion", "at": [31.5, 8.3], "text": "K^+", "hue": "green", "r": 0.36, "text_size": 0.45},
    {"g": "ion", "at": [31.5, 5.7], "text": "K^+", "hue": "green", "r": 0.36, "text_size": 0.45},
    {"g": "arrow", "from": [33.3, 10.4], "to": [33.3, 3.6], "color": "cyan"},
    _lbl(31.5, "K^+ channel\n(gated pore)"),
    {"id": "abc", "g": "transporter", "at": [39, 7], "hue": "pink", "open": "in"},
    {"g": "blob", "at": [38.2, 10.6], "w": 1.4, "h": 1.0, "shape": "kidney", "hue": "pink", "text": "ATP", "text_size": 0.5},
    {"g": "blob", "at": [39.8, 10.6], "w": 1.4, "h": 1.0, "shape": "kidney", "hue": "pink", "text": "ATP", "text_size": 0.5},
    {"g": "ligand", "at": [41.6, 9.6], "shape": "hexagon", "hue": "red", "r": 0.45},
    {"g": "arrow", "from": [41.6, 8.9], "to": [41.6, 3.8], "color": "magenta"},
    {"g": "ligand", "at": [41.6, 3.1], "shape": "hexagon", "hue": "red", "r": 0.45},
    _lbl(39.5, "ABC transporter\n(MDR1: drug efflux)"),
]}

CAMP_SCENE = {"items": [
    {"g": "text", "at": [-1.2, 0.8], "text": "Outside", "align": "left", "size": 0.9, "color": "gray"},
    {"g": "text", "at": [-1.2, 7.9], "text": "Cytosol", "align": "left", "size": 0.9, "color": "gray"},
    {"g": "membrane", "y": 5, "x0": -1.2, "x1": 31},
    {"id": "epi", "g": "ligand", "at": [5, 1.4], "shape": "triangle", "hue": "red", "label": "Epinephrine", "label_at": [8.6, 0.6]},
    {"id": "gpcr", "g": "gpcr", "at": [5, 5], "hue": "purple", "label": "β-adrenergic receptor\n(GPCR, GEF)", "label_at": [3.4, 10.4],
     "point": [5.6, 7.3]},
    {"id": "gs", "g": "gprotein", "at": [11.2, 7.25], "nucleotide": "GTP", "accent": True,
     "label": "G_s: α–GTP", "label_at": [10.8, 10.4], "point": [10.4, 8.4]},
    {"id": "ac", "g": "blob", "at": [17.5, 5], "w": 2.4, "h": 3.4, "shape": "box", "hue": "teal"},
    {"id": "ac_cat", "g": "blob", "at": [17.5, 7.9], "w": 2.6, "h": 1.5, "shape": "kidney", "hue": "teal", "accent": True,
     "label": "Adenylyl cyclase", "label_at": [21.6, 1.4], "point": [18.2, 3.6]},
    {"g": "arrow", "from": [12.6, 6.6], "to": [16.0, 7.4], "color": "green", "curve": -0.15},
    {"g": "ion", "at": [15.6, 10.6], "text": "ATP", "hue": "gray", "r": 0.55, "text_size": 0.5},
    {"g": "arrow", "from": [16.3, 10.6], "to": [18.6, 10.6]},
    {"id": "camp", "g": "ion", "at": [19.3, 10.6], "text": "cAMP", "hue": "yellow", "r": 0.62, "text_size": 0.45},
    {"g": "arrow", "from": [20.0, 10.6], "to": [23.2, 10.6], "color": "gray"},
    {"g": "text", "at": [21.6, 9.9], "text": "PDE", "size": 0.75, "color": "cyan"},
    {"g": "ion", "at": [23.9, 10.6], "text": "AMP", "hue": "gray", "r": 0.55, "text_size": 0.45},
    # PKA: R2C2 → R2(cAMP)4 + 2 C
    {"g": "blob", "at": [3.0, 13.2], "w": 1.6, "h": 1.0, "hue": "green"},
    {"g": "blob", "at": [3.0, 14.3], "w": 1.6, "h": 1.0, "hue": "green"},
    {"g": "blob", "at": [4.9, 13.2], "w": 2.0, "h": 1.0, "shape": "box", "hue": "gray", "text": "R", "text_size": 0.6},
    {"g": "blob", "at": [4.9, 14.3], "w": 2.0, "h": 1.0, "shape": "box", "hue": "gray", "text": "R", "text_size": 0.6},
    {"g": "text", "at": [4.0, 16.0], "text": "PKA (R_2C_2): inactive", "size": 0.8},
    {"g": "arrow", "from": [18.8, 11.3], "to": [9.6, 13.0], "curve": -0.12, "color": "black"},
    {"g": "text", "at": [13.6, 11.9], "text": "4 cAMP", "size": 0.8},
    {"g": "blob", "at": [12.5, 13.6], "w": 2.0, "h": 1.0, "shape": "box", "hue": "gray", "text": "R", "text_size": 0.6},
    {"g": "blob", "at": [12.5, 14.7], "w": 2.0, "h": 1.0, "shape": "box", "hue": "gray", "text": "R", "text_size": 0.6},
    *[{"g": "ion", "at": [13.8, y], "r": 0.27, "hue": "yellow"} for y in (13.35, 13.85, 14.45, 14.95)],
    {"g": "arrow", "from": [7.0, 13.8], "to": [10.6, 13.8], "style": "fat"},
    {"id": "c1", "g": "blob", "at": [17.4, 13.3], "w": 1.6, "h": 1.0, "hue": "green", "accent": True, "text": "C", "text_size": 0.6},
    {"id": "c2", "g": "blob", "at": [17.4, 14.8], "w": 1.6, "h": 1.0, "hue": "green", "accent": True, "text": "C", "text_size": 0.6},
    {"g": "text", "at": [17.4, 16.3], "text": "2 C: active", "size": 0.8, "color": "green"},
    {"g": "arrow", "from": [18.6, 14.0], "to": [22.4, 14.0]},
    {"g": "text", "at": [22.8, 14.0], "text": "Ser/Thr-Ⓟ:\nphosphorylase kinase ↑\nglycogen synthase ↓\nCREB → transcription",
     "align": "left", "size": 0.8},
]}

ILLUS = {"membrane_transporters": TRANSPORT_SCENE, "camp_pka_cascade": CAMP_SCENE}
for _p in PATHWAYS:
    if _p["id"] in ILLUS:
        _p["variants"] = {"illus": {"note_zh": "圖片版（插畫）", "scene": ILLUS[_p["id"]]}}


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


PARENT = {"gpcr": "signal_transduction", "rtk": "signal_transduction", "second_messenger": "signal_transduction",
          "apoptosis": "signal_transduction", "cancer": "signal_transduction",
          "membrane_transport": "membrane", "membrane_lipid": "membrane", "lipoprotein": "lipid_metabolism"}


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
        for k in ("chain", "free", "haworth", "ring", "decorate"):
            if m.get(k):
                d[k] = m[k]
        d["tags"] = with_parent(m["tags"])
        d["notes_zh"] = m["notes"]
        d["status"] = "pending"
        dump(data / "molecules" / f"{m['id']}.yaml", d)
    for r in REACTIONS:
        d = {"id": r["id"], "title": r["title"], "substrates": r["substrates"], "products": r["products"]}
        for k in ("cofactors_in", "cofactors_out", "extra_in", "extra_out", "prosthetic"):
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
    print(f"{len(MOLECULES)} molecules, {len(REACTIONS)} reactions, {len(PATHWAYS)} pathways")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data")
