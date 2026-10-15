"""Round 16: gaps found by comparing with the Book 2 metabolism lecture notes (生化_筆_代謝.pdf).

    python tools/gen_met.py data

Plots (energy charge, hexokinase vs glucokinase, glycolysis free-energy profile, phosphoryl-transfer
potential), reactions (pyruvate decarboxylase, lactose synthase, UDP-glucose dehydrogenase), and
overviews (fermentation, glucose carbon fate, PPP non-oxidative phase, hexose interconversion hub,
TCA / β-oxidation analogy, glycogen enzyme phosphorylation, xenobiotic detoxification, futile cycle,
ATP cleavage modes).
"""
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent))
from reaction_colors import COLORS, NADH_M, chain, frag, hw, whole  # noqa: E402

M, C = "magenta", "cyan"
RED, GREEN, ORANGE, GRAY = "red", "green", "orange", "gray"


def n(id_, ref, x, y, **kw):
    return {"id": id_, "ref": ref, "pos": [x, y], "view": "name", **kw}


def t(id_, text, x, y, kind="concept", **kw):
    return {"id": id_, "text": text, "pos": [x, y], "kind": kind, **kw}


# ================================================================ plots

ENERGY_CHARGE = {
    "size": [15, 9.5],
    "x": {"min": 0, "max": 1, "label": "Energy charge", "ticks": [0, 0.25, 0.5, 0.75, 1]},
    "y": {"min": 0, "max": 1.15, "label": "Relative rate", "label_pos": "center", "tick_labels": False, "ticks": []},
    "bands": [{"x0": 0.8, "x1": 0.95, "label": "most cells", "color": GREEN}],
    "curves": [
        {"id": "cat", "y": "1/(1+exp((x-0.78)*18))", "color": M, "width": 1.2},
        {"id": "ana", "y": "1/(1+exp(-(x-0.72)*18))", "color": C, "width": 1.2},
    ],
    "texts": [{"at": [0.05, 0.86], "text": "ATP-generating\n(catabolic) pathways", "align": "left", "color": M},
              {"at": [0.05, 0.16], "text": "ATP-utilizing\n(anabolic) pathways", "align": "left", "color": C}],
}

HK_GK = {
    "size": [15, 9.5],
    "params": {"KmHK": 0.1, "S05": 10, "nGK": 1.7},
    "x": {"min": 0, "max": 20, "label": "[Glucose] (mM)", "ticks": [0, 5, 10, 15, 20]},
    "y": {"min": 0, "max": 1.15, "label": "v_0 / V_{max}", "label_pos": "center", "ticks": [0, 0.5, 1]},
    "bands": [{"x0": 4, "x1": 6, "label": "fasting", "color": C},
              {"x0": 12, "x1": 16, "label": "portal vein (fed)", "color": M}],
    "curves": [
        {"id": "hk", "y": "mm(x, 1, KmHK)", "width": 1.2},
        {"id": "gk", "y": "hill(x, S05, nGK)", "color": M, "width": 1.2},
    ],
    "texts": [{"at": [6.6, 0.87], "text": "Hexokinase I (muscle, brain):\nK_m ~ 0.1 mM, inhibited by G6P", "align": "left", "size": 0.85},
              {"at": [10.6, 0.2], "text": "Glucokinase (hexokinase IV;\nliver, β cell): S_{0.5} ~ 10 mM,\nsigmoid, not inhibited by G6P",
               "align": "left", "size": 0.85, "color": M}],
    "points": [{"on": "hk", "x": "KmHK", "y": 0.5, "drop": "x"}, {"on": "gk", "x": "S05", "y": 0.5, "drop": "x", "color": M}],
}

# ΔG (actual, erythrocyte; Lehninger Table 14-2) — steps 6–10 occur twice per glucose
_DG = [("Glucose", 0), ("G6P", -33.4), ("F6P", -2.5), ("F1,6BP", -22.2), ("DHAP + GAP", -1.3), ("2 GAP", 2.5),
       ("2 1,3-BPG", -3.4), ("2 3-PG", 2.6), ("2 2-PG", 1.6), ("2 PEP", -6.6), ("2 Pyruvate", -33.4)]


def _glycolysis_profile():
    g, curves, texts, total = 0.0, [], [], 0.0
    for i, (name, d) in enumerate(_DG):
        g += d
        curves.append({"id": f"s{i}", "y": f"{g:.1f}", "from": i - 0.32, "to": i + 0.32, "width": 1.4,
                       "color": RED if i in (1, 3, 10) else None})
        texts.append({"at": [i, g + (5 if i % 2 == 0 else -6)], "text": name, "size": 0.75})
    curves = [{k: v for k, v in c.items() if v is not None} for c in curves]
    return {
        "size": [22, 11],
        "x": {"min": -0.9, "max": 10.6, "label": "Glycolytic step", "ticks": list(range(1, 11)), "at": -105},
        "y": {"min": -105, "max": 10, "label": "Cumulative ΔG (kJ/mol)", "label_pos": "center", "ticks": [0, -25, -50, -75, -100]},
        "curves": curves,
        "texts": texts + [{"at": [0.0, -88], "text": "Large drops (red) = steps 1, 3, 10:\nhexokinase, PFK-1, pyruvate kinase\n→ irreversible, regulated, bypassed\nin gluconeogenesis",
                           "align": "left", "size": 0.85, "color": RED},
                          {"at": [6.0, -20], "text": "steps 2, 4–9: ΔG ≈ 0\n(near equilibrium, reversible)", "size": 0.85, "color": GRAY}],
    }


PHOSPHO_POTENTIAL = {
    "size": [15, 11],
    "x": {"min": 0, "max": 10, "label": "", "ticks": [], "at": -16},
    "y": {"min": -16, "max": 0, "label": "ΔG°' of hydrolysis (kcal/mol)", "label_pos": "center",
          "ticks": [0, -4, -8, -12, -16], "at": 0},
    "curves": [
        {"id": "pep", "y": -14.8, "from": 0.4, "to": 3.2, "color": M, "width": 1.4},
        {"id": "bpg", "y": -11.8, "from": 0.4, "to": 3.2, "color": M, "width": 1.4},
        {"id": "pcr", "y": -10.3, "from": 0.4, "to": 3.2, "color": M, "width": 1.4},
        {"id": "atp", "y": -7.3, "from": 3.6, "to": 6.2, "width": 1.6},
        {"id": "g1p", "y": -5.0, "from": 6.6, "to": 7.6, "color": C, "width": 1.4},
        {"id": "g6p", "y": -3.3, "from": 6.6, "to": 7.6, "color": C, "width": 1.4},
        {"id": "g3p", "y": -2.2, "from": 6.6, "to": 7.6, "color": C, "width": 1.4},
    ],
    "texts": [
        {"at": [1.8, -14.2], "text": "PEP  −14.8", "size": 0.85}, {"at": [1.8, -11.2], "text": "1,3-BPG  −11.8", "size": 0.85},
        {"at": [1.8, -9.7], "text": "Creatine-P  −10.3", "size": 0.85},
        {"at": [5.7, -8.2], "text": "ATP → ADP  −7.3", "size": 0.9},
        {"at": [7.8, -5.0], "text": "Glucose 1-P  −5.0", "size": 0.85, "align": "left"},
        {"at": [7.8, -3.3], "text": "Glucose 6-P  −3.3", "size": 0.85, "align": "left"},
        {"at": [7.8, -2.2], "text": "Glycerol 3-P  −2.2", "size": 0.85, "align": "left"},
        {"at": [2.0, -2.5], "text": "High-energy\n(Ⓟ donors)", "size": 0.9, "color": M},
        {"at": [8.0, -12.5], "text": "Low-energy\n(Ⓟ acceptors\nform these)", "size": 0.9, "color": C},
    ],
    "arrows": [{"from": [3.25, -14.6], "to": [3.9, -7.5], "color": M}, {"from": [3.25, -11.6], "to": [4.2, -7.5], "color": M},
               {"from": [6.25, -7.1], "to": [6.6, -5.3], "color": C}],
}


# ================================================================ molecules

UDP_TAIL = "OP(=O)([O-])OP(=O)([O-])OC[C@H]5O[C@@H](n6ccc(=O)[nH]c6=O)[C@H](O)[C@@H]5O"
MOLECULES = [
    {"id": "udp_glucuronate", "name": "UDP-glucuronate", "label": "UDP-glucuronate", "aliases": ["UDP-GlcA", "UDP-葡萄糖醛酸"],
     "smiles": f"[O-]C(=O)[C@@H]1[C@H]([C@@H]([C@H]([C@H](O1){UDP_TAIL})O)O)O",
     "haworth": {"ring": "pyranose", "numbers": True,
                 "subs": {1: {"up": "H", "down": "O–UDP"}, 2: {"up": "H", "down": "OH"}, 3: {"up": "OH", "down": "H"},
                          4: {"up": "H", "down": "HO"}, 5: {"up": "COO^-", "down": "H"}}},
     "tags": ["carb_metabolism", "sugar_derivative", "detoxification"],
     "notes": "UDP-glucose 的 C6 被 UDP-glucose dehydrogenase 氧化成羧酸。活化的 glucuronate 供體：bilirubin、藥物、類固醇的 phase II 結合"
              "（UDP-glucuronosyltransferase），以及 GAG（hyaluronate、chondroitin、heparin）的 GlcA 單元。植物與多數動物再經 uronic acid pathway 合成 vitamin C；人缺 L-gulonolactone oxidase。"},
]

# ================================================================ reactions

REACTIONS = [
    {"id": "pyruvate_decarboxylase", "title": "Pyruvate → Acetaldehyde (pyruvate decarboxylase)",
     "substrates": ["pyruvate"], "products": ["acetaldehyde", "co2"], "cofactors_in": ["H+"],
     "prosthetic": ["TPP,", "Mg^{2+}"], "enzyme": "pyruvate decarboxylase", "reversible": False,
     "card": {"location": "酵母菌細胞質（人體沒有此酵素）",
              "type": "非氧化性脫羧：COO⁻ 以 CO₂ 離開（洋紅），原 C2 變成醛",
              "cofactors": "TPP（thiamine pyrophosphate）、Mg²⁺；不需 NAD⁺（與 PDH 的氧化性脫羧比較）",
              "mechanism": "TPP 的 thiazolium C2 carbanion 攻擊 pyruvate 的羰基 → 脫羧 → hydroxyethyl-TPP → 釋出 acetaldehyde",
              "exam_points": ["酒精發酵：pyruvate → acetaldehyde（PDC）→ ethanol（ADH，NADH 再生 NAD⁺）",
                              "C1-¹⁴C-glucose → ¹⁴CH₃–CH₂OH；CO₂ 來自 glucose C3、C4",
                              "TPP 依賴脫羧：PDH、α-KGDH、BCKDH（氧化性）vs PDC（非氧化性）；transketolase 也用 TPP"]},
     "tags": ["glycolysis", "ethanol"]},
    {"id": "lactose_synthase", "title": "UDP-galactose + Glucose → Lactose (lactose synthase)",
     "substrates": ["udp_galactose", "b_glucopyranose"], "products": ["lactose", {"id": "udp", "view": "formula"}],
     "extra_out": ["H+"], "enzyme": "lactose synthase\n(galactosyltransferase + α-lactalbumin)", "reversible": False,
     "card": {"location": "乳腺細胞 Golgi",
              "type": "醣基轉移：Gal 的 C1 接到 glucose 的 C4–O（洋紅）形成 β1→4 醣苷鍵；UDP 帶著原本 C1 的 O 離開（青）",
              "regulation": "α-lactalbumin（調節次單元）使 galactosyltransferase 的 K_m(glucose) 大幅下降；prolactin 誘導 α-lactalbumin，懷孕時 progesterone 抑制",
              "exam_points": ["平時 galactosyltransferase 把 Gal 接到 GlcNAc（醣蛋白）；哺乳期加上 α-lactalbumin → 改接 glucose 生成乳糖",
                              "Gal 來源：UDP-Glc ⇌ UDP-Gal（4-epimerase），所以不吃半乳糖也能合成乳糖"]},
     "tags": ["galactose_metabolism", "disaccharide"]},
    {"id": "udp_glucose_dehydrogenase", "title": "UDP-glucose → UDP-glucuronate (UDP-glucose dehydrogenase)",
     "substrates": ["udp_glucose"], "products": ["udp_glucuronate"],
     "cofactors_in": ["H2O", "2NAD+"], "cofactors_out": ["2NADH+2H+"], "extra_out": ["H+"],
     "enzyme": "UDP-glucose dehydrogenase", "reversible": False,
     "card": {"location": "細胞質（肝）",
              "type": "C6 的 –CH₂OH 經兩次 NAD⁺ 依賴氧化（醇 → 醛 → 酸）成為 –COO⁻；第二步加入 H₂O",
              "exam_points": ["產物 UDP-GlcA：bilirubin、藥物、類固醇的 glucuronidation（phase II，增加水溶性）",
                              "GAG 的 uronic acid 來源；L-iduronate 由 GlcA 的 C5 epimerization 而來",
                              "uronic acid pathway：glucuronate → L-gulonate → （人缺 gulonolactone oxidase）無法合成 vitamin C；L-xylulose 缺乏 → essential pentosuria"]},
     "tags": ["carb_metabolism", "sugar_derivative", "detoxification"]},
]

COLORS.update({
    "pyruvate_decarboxylase": {"pyruvate": chain(r0=M), "co2": whole(), "H+": whole(C),
                               "acetaldehyde": chain(r0={"fork_r": C})},
    "lactose_synthase": {"udp_galactose": hw(d1=C), "udp": whole(C), "b_glucopyranose": hw(d4=[M, "O"]),
                         "lactose": {"colors": {"link0": M}}},
    "udp_glucose_dehydrogenase": {"udp_glucose": hw(u5=[M, "H_2"]), "2NADH+2H+": NADH_M,
                                  "H2O": whole(C), "udp_glucuronate": hw(u5=[C, "O^-"])},
})

NEW_COFACTORS = {
    "2NAD+": {"label": "2 NAD^+", "smiles": "NC(=O)c1ccc[n+]([2*])c1.NC(=O)c1ccc[n+]([2*])c1", "molecule": "nad"},
    "2NADH+2H+": {"label": "2 NADH + 2 H^+", "smiles": "NC(=O)C1=CN([2*])C=CC1.[H+].NC(=O)C1=CN([2*])C=CC1.[H+]", "molecule": "nadh"},
}


# ================================================================ pathways

def mk(*nums):
    """marks for one chain row: mk("1") on the centre, mk(("l", "1"), ("c", "2")) on parts."""
    if len(nums) == 1 and isinstance(nums[0], str):
        return {"c": nums[0]}
    return {k: v for k, v in nums}


def glucose_fate():
    glc = {"marks": {i: mk(str(i + 1)) for i in range(6)}}
    return {
        "id": "glucose_carbon_fate", "title": "Fate of glucose carbons in glycolysis",
        "layout": "free", "cell_size": 5.2,
        "nodes": [
            {"id": "glc", "ref": "glucose", "pos": [0, 0], "decorate": glc},
            {"id": "dhap", "ref": "dhap", "pos": [-1.4, 1.85], "decorate": {"marks": {0: mk("1"), 1: mk("2"), 2: mk("3")}}},
            {"id": "gap2", "ref": "glyceraldehyde_3p", "pos": [1.4, 1.85], "decorate": {"marks": {0: mk("4"), 1: mk("5"), 2: mk("6")}}},
            {"id": "gap1", "ref": "glyceraldehyde_3p", "pos": [-1.4, 3.95], "decorate": {"marks": {0: mk("3"), 1: mk("2"), 2: mk("1")}}},
            {"id": "pyr1", "ref": "pyruvate", "pos": [-1.4, 5.95], "decorate": {"marks": {0: mk("3"), 1: mk("2"), 2: mk("1")}}},
            {"id": "pyr2", "ref": "pyruvate", "pos": [1.4, 5.95], "decorate": {"marks": {0: mk("4"), 1: mk("5"), 2: mk("6")}}},
            {"id": "ac1", "ref": "acetyl_coa", "pos": [-1.4, 7.4], "decorate": {"marks": {0: mk(("l", "1"), ("c", "2"))}}},
            {"id": "ac2", "ref": "acetyl_coa", "pos": [1.4, 7.4], "decorate": {"marks": {0: mk(("l", "6"), ("c", "5"))}}},
            t("co2", "CO_2 (PDH) from C3 and C4", 0, 8.05, "name"),
        ],
        "edges": [
            {"from": "glc", "to": "dhap", "label": "steps 1–4\n(aldolase)", "label_side": "right"},
            {"from": "glc", "to": "gap2", "label": "aldolase"},
            {"from": "dhap", "to": "gap1", "step": "triose_phosphate_isomerase", "label": "TPI", "reversible": True},
            {"from": "gap1", "to": "pyr1", "label": "steps 6–10"},
            {"from": "gap2", "to": "pyr2", "label": "steps 6–10", "label_side": "right"},
            {"from": "pyr1", "to": "ac1", "step": "pdh", "label": "PDH"},
            {"from": "pyr2", "to": "ac2", "step": "pdh", "label": "PDH", "label_side": "right"},
        ],
        "tags": ["glycolysis", "TCA"],
        "notes": "葡萄糖碳原子去向（¹⁴C 追蹤）：aldolase 把 C1–C3 分到 DHAP、C4–C6 分到 GAP；TPI 把 DHAP 的 C1（CH₂OⓅ）變成 GAP 的 C3。"
                 "因此 pyruvate：CH₃ = glucose C1／C6，C=O = C2／C5，COO⁻ = C3／C4。"
                 "PDH（或酵母的 pyruvate decarboxylase）放出的 CO₂ 來自 glucose C3、C4；acetyl-CoA 的 CH₃ 來自 C1／C6。"
                 "例：1-¹⁴C-glucose → ¹⁴CH₃–CH₂OH（酒精發酵）、¹⁴CH₃–CO–SCoA。"
                 "糖質新生反向：pyruvate C3（CH₃）→ glucose C1 與 C6；C2 → C2、C5；C1（COO⁻）→ C3、C4（PC 加入的 CO₂ 會在 PEPCK 離開）。",
    }


PATHWAYS = [
    {"id": "energy_charge", "title": "Energy charge regulates catabolism and anabolism", "layout": "plot", "plot": ENERGY_CHARGE,
     "caption": {"text": "EC = ([ATP] + 1/2 [ADP]) / ([ATP] + [ADP] + [AMP])"},
     "tags": ["bioenergetics"],
     "notes": "Energy charge（Atkinson）＝([ATP] + ½[ADP]) / ([ATP] + [ADP] + [AMP])，0–1；多數細胞維持 0.8–0.95。"
              "EC 高 → 產 ATP 的分解路徑（糖解、TCA）被抑制、耗 ATP 的合成路徑被促進；EC 低則相反。"
              "高能訊號：ATP、acetyl-CoA、citrate、NADH（抑制 PFK-1、PDH、IDH…）；低能訊號：ADP、AMP、NAD⁺（AMP 活化 PFK-1、glycogen phosphorylase b、AMPK）。"
              "Adenylate kinase：2 ADP ⇌ ATP + AMP，使 AMP 成為敏感的低能訊號。"},
    {"id": "hexokinase_glucokinase", "title": "Hexokinase vs glucokinase", "layout": "plot", "plot": HK_GK,
     "tags": ["glycolysis", "enzyme_kinetics"],
     "notes": "Hexokinase I–III（肌肉、腦等）：K_m ≈ 0.1 mM，血糖正常時即飽和；產物 G6P 抑制。"
              "Glucokinase（hexokinase IV，肝、胰 β 細胞）：S₀.₅ ≈ 10 mM、cooperative（S 形）、不受 G6P 抑制，飯後門脈血糖高時才大量作用 → 肝把多餘葡萄糖存成肝醣；β 細胞以它當血糖感測器（MODY2 為 GK 突變）。"
              "調控：GKRP（glucokinase regulatory protein）在 F6P 高時把 GK 拉進細胞核（失活），glucose 高時釋出；insulin 誘導 GK 合成。肝與 β 細胞都用 GLUT2（高 K_m）。"},
    {"id": "glycolysis_free_energy", "title": "Free-energy profile of glycolysis", "layout": "plot", "plot": _glycolysis_profile(),
     "tags": ["glycolysis", "bioenergetics"],
     "notes": "紅血球中實際 ΔG（Lehninger）：只有 hexokinase（−33.4）、PFK-1（−22.2）、pyruvate kinase（−16.7 × 2）三步大幅負值 → 不可逆、為調控點，"
              "糖質新生必須以 G6Pase、FBPase-1、PC + PEPCK 繞過。其餘七步 ΔG ≈ 0，接近平衡，兩方向共用。"
              "ΔG°′ 為正的步驟（aldolase +23.8、TPI +7.5 kJ/mol）在細胞內因產物迅速移除而 ΔG ≈ 0。第 6–10 步每分子葡萄糖做兩次。"},
    {"id": "phosphoryl_transfer_potential", "title": "Phosphoryl-transfer potential", "layout": "plot", "plot": PHOSPHO_POTENTIAL,
     "tags": ["bioenergetics", "glycolysis"],
     "notes": "水解 ΔG°′ 比 ATP 更負的化合物（PEP −14.8、1,3-BPG −11.8、phosphocreatine −10.3 kcal/mol）可把 Ⓟ 交給 ADP（受質層次磷酸化：PK、PGK、creatine kinase）；"
              "ATP（−7.3）再把 Ⓟ 交給 glucose 等形成低能磷酸酯（G6P −3.3）。ATP 位居中間 → 能量的「通貨」。"
              "PEP 特別高：水解後 enol → keto 互變（tautomerization）額外釋能。1,3-BPG：acyl-phosphate（混合酸酐）。"
              "ATP 高能的原因：① 磷酸間負電排斥 ② 產物的共振穩定較大 ③ 產物離子化／溶劑化（ADP³⁻ + Pi²⁻ + H⁺）。"},
    {"id": "fermentation", "title": "Fermentation: regenerating NAD+", "layout": "free", "cell_size": 4.6,
     "nodes": [
         n("glc", "glucose", 0, 0), n("pyr", "pyruvate", 0, 1.3),
         n("lac", "lactate", -1.4, 2.5), n("acal", "acetaldehyde", 1.4, 2.5), n("etoh", "ethanol", 1.4, 3.7),
         t("hom", "Muscle (anaerobic), RBC,\nlactic acid bacteria", -1.4, 3.2, "name"),
         t("yeast", "Yeast (alcoholic fermentation)", 1.4, 4.3, "name"),
         t("aer", "O_2 present: NADH → ETC,\npyruvate → acetyl-CoA (PDH)", -2.2, 1.3, "name"),
     ],
     "edges": [
         {"from": "glc", "to": "pyr", "label": "glycolysis\n(GAPDH)", "side_text": "2 NAD^+ → 2 NADH", "label_side": "left"},
         {"from": "pyr", "to": "lac", "step": "ldh", "label": "lactate DH", "side_text": "NADH → NAD^+"},
         {"from": "pyr", "to": "acal", "step": "pyruvate_decarboxylase", "label": "pyruvate\ndecarboxylase (TPP)",
          "side_text": "→ CO_2", "label_side": "left"},
         {"from": "acal", "to": "etoh", "step": "alcohol_dehydrogenase", "label": "alcohol DH",
          "side_text": "NADH → NAD^+", "label_side": "right"},
     ],
     "tags": ["glycolysis", "ethanol"],
     "notes": "無氧時 NADH 無法經 ETC 氧化，GAPDH 需要的 NAD⁺ 必須在細胞質再生：乳酸發酵（LDH）或酒精發酵（PDC + ADH）。"
              "兩者都淨得 2 ATP／glucose，沒有淨氧化還原（GAPDH 產生的 NADH 在最後一步用掉）。"
              "Pasteur effect：有氧時葡萄糖消耗速率下降（每分子產 ~30–32 ATP，ATP、citrate 抑制 PFK-1）。"
              "人體沒有 pyruvate decarboxylase；肝的 ADH 方向相反（ethanol → acetaldehyde）。"},
    glucose_fate(),
    {"id": "ppp_nonoxidative", "title": "Pentose phosphate pathway: non-oxidative phase", "layout": "free", "cell_size": 4.4,
     "nodes": [
         n("x5p", "xylulose_5p", 0, 0), n("r5p", "ribose_5p", 0, 1.2),
         n("s7p", "sedoheptulose_7p", 2.2, 0), n("gap", "glyceraldehyde_3p", 2.2, 1.2),
         n("f6p", "fructose_6p", 4.4, 0), n("e4p", "erythrose_4p", 4.4, 1.2),
         n("x5p2", "xylulose_5p", 4.4, 2.4), n("f6p2", "fructose_6p", 6.6, 1.2), n("gap2", "glyceraldehyde_3p", 6.6, 2.4),
         t("c1", "C_5 + C_5", 0, -0.75, "plain"), t("c2", "C_7 + C_3", 2.2, -0.75, "plain"),
         t("c3", "C_6 + C_4", 4.4, -0.75, "plain"), t("c4", "C_6 + C_3", 6.6, 0.45, "plain"),
         t("net", "Net: 3 C_5 → 2 C_6 + C_3\n(3 Ru5P → 2 F6P + GAP)", 1.1, 3.0, "concept"),
     ],
     "edges": [
         {"from": "r5p", "to": "gap", "step": "transketolase_1", "label": "transketolase\n(TPP, 2C)", "label_side": "right"},
         {"from": "x5p", "to": "s7p", "step": "transketolase_1"},
         {"from": "s7p", "to": "f6p", "step": "transaldolase", "label": "transaldolase\n(3C)"},
         {"from": "gap", "to": "e4p", "step": "transaldolase"},
         {"from": "e4p", "to": "f6p2", "step": "transketolase_2", "label": "transketolase\n(2C)"},
         {"from": "x5p2", "to": "gap2", "step": "transketolase_2"},
     ],
     "tags": ["PPP"],
     "notes": "非氧化期：Ru5P 經 isomerase → R5P（核苷酸）或 3-epimerase → Xu5P。"
              "Transketolase（TPP）把酮糖供體（Xu5P）的 C1–C2 轉給醛糖受體；transaldolase 轉 C1–C3（不需 TPP，Schiff base）。"
              "口訣：酮糖（上盤）供、醛糖（底座）受；TK「左 2」、TA「右 3」。"
              "全部可逆：需 NADPH 多於 R5P 時 → 轉成 F6P、GAP 回糖解；需 R5P 多於 NADPH（分裂中細胞）→ 由 F6P、GAP 逆向生成 R5P，不經氧化期。"
              "¹⁴C：glucose C1 在 6PG DH 以 CO₂ 離開。紅血球 transketolase 活性可評估 thiamine（B₁）狀態（Wernicke）。"},
    {"id": "hexose_interconversion", "title": "Hexose interconversions around UDP-glucose", "layout": "free", "cell_size": 4.4,
     "nodes": [
         n("glc", "glucose", 0, 0), n("g6p", "a_glucose_6p", 0, 1.2), n("g1p", "a_glucose_1p", 0, 2.4),
         n("udpg", "udp_glucose", 0, 3.6), n("gly", "glycogen", -2.0, 3.6),
         n("udpgal", "udp_galactose", 2.2, 3.6), n("lac", "lactose", 4.2, 3.6),
         t("gp", "Glycoproteins,\nglycolipids", 2.2, 4.9, "name"),
         n("udpga", "udp_glucuronate", 0, 4.9), t("conj", "Glucuronides (bilirubin, drugs),\nGAGs, vitamin C (not humans)", 0, 6.0, "name"),
         n("gal", "a_galactopyranose", 2.4, 1.2), n("gal1p", "a_galactose_1p", 2.4, 2.4),
         n("fru", "d_fructose", -2.2, 0), n("f1p", "fructose_1p", -2.2, 1.2),
         t("tri", "DHAP + glyceraldehyde\n→ glycolysis", -2.2, 2.4, "name"),
     ],
     "edges": [
         {"from": "glc", "to": "g6p", "step": "hexokinase", "label": "hexokinase"},
         {"from": "g6p", "to": "g1p", "step": "phosphoglucomutase", "label": "PGM", "reversible": True},
         {"from": "g1p", "to": "udpg", "step": "udp_glucose_pyrophosphorylase", "label": "UDP-Glc\npyrophosphorylase"},
         {"from": "udpg", "to": "gly", "step": "glycogen_synthase", "label": "glycogen\nsynthase"},
         {"from": "udpg", "to": "udpgal", "step": "udp_galactose_epimerase", "label": "4-epimerase", "reversible": True, "label_side": "right"},
         {"from": "udpgal", "to": "lac", "step": "lactose_synthase", "label": "lactose\nsynthase"},
         {"from": "udpgal", "to": "gp", "label": "galactosyl-\ntransferase", "label_side": "right"},
         {"from": "udpg", "to": "udpga", "step": "udp_glucose_dehydrogenase", "label": "UDP-Glc DH\n(2 NAD^+)"},
         {"from": "udpga", "to": "conj", "label": "UGT"},
         {"from": "gal", "to": "gal1p", "step": "galactokinase", "label": "galactokinase"},
         {"from": "gal1p", "to": "g1p", "step": "galt", "label": "GALT\n(+ UDP-Glc)", "label_side": "right"},
         {"from": "fru", "to": "f1p", "step": "fructokinase", "label": "fructokinase"},
         {"from": "f1p", "to": "tri", "step": "aldolase_b", "label": "aldolase B"},
     ],
     "tags": ["carb_metabolism", "galactose_metabolism", "fructose_metabolism", "glycogen"],
     "notes": "UDP-glucose 是醣類互換的樞紐（高能活化醣）：→ 肝醣；→ UDP-galactose（4-epimerase）→ 乳糖、醣蛋白；→ UDP-glucuronate（2 NAD⁺）→ 結合解毒、GAG。"
              "Galactose 經 GALT 換成 G1P 進入主線；fructose 在肝經 fructokinase、aldolase B 成為三碳糖（繞過 PFK-1，不受 insulin／PFK-1 調控 → 代謝比 glucose 快）。"
              "缺乏：GALT → classic galactosemia；aldolase B → hereditary fructose intolerance；fructokinase → essential fructosuria（良性）。"},
    {"id": "tca_beta_oxidation_analogy", "title": "Same chemistry: TCA steps 6–8 and β-oxidation", "layout": "free", "cell_size": 4.4,
     "nodes": [
         {"id": "a1", "ref": "acyl_coa", "pos": [0, 0]}, {"id": "a2", "ref": "enoyl_coa", "pos": [2.4, 0]},
         {"id": "a3", "ref": "hydroxyacyl_coa", "pos": [4.8, 0]}, {"id": "a4", "ref": "ketoacyl_coa", "pos": [7.2, 0]},
         {"id": "b1", "ref": "succinate", "pos": [0, 1.9]}, {"id": "b2", "ref": "fumarate", "pos": [2.4, 1.9]},
         {"id": "b3", "ref": "malate", "pos": [4.8, 1.9]}, {"id": "b4", "ref": "oxaloacetate", "pos": [7.2, 1.9]},
         t("h1", "β-Oxidation", -1.9, 0, "heading"), t("h2", "TCA cycle", -1.9, 1.9, "heading"),
         t("h3", "Ile\ndegradation", -1.9, 3.2, "heading"),
         t("c1", "2-Methyl-\nbutyryl-CoA", 0, 3.2, "name"), t("c2", "Tiglyl-CoA", 2.4, 3.2, "name"),
         t("c3", "2-Methyl-3-hydroxy-\nbutyryl-CoA", 4.8, 3.2, "name"), t("c4", "2-Methyl-\nacetoacetyl-CoA", 7.2, 3.2, "name"),
     ],
     "edges": [
         {"from": "a1", "to": "a2", "step": "acyl_coa_dehydrogenase", "label": "FAD"},
         {"from": "a2", "to": "a3", "step": "enoyl_coa_hydratase", "label": "+ H_2O"},
         {"from": "a3", "to": "a4", "step": "hydroxyacyl_coa_dehydrogenase", "label": "NAD^+"},
         {"from": "b1", "to": "b2", "step": "succinate_dehydrogenase", "label": "FAD"},
         {"from": "b2", "to": "b3", "step": "fumarase", "label": "+ H_2O"},
         {"from": "b3", "to": "b4", "step": "malate_dehydrogenase", "label": "NAD^+"},
         {"from": "c1", "to": "c2", "label": "FAD"}, {"from": "c2", "to": "c3", "label": "+ H_2O"},
         {"from": "c3", "to": "c4", "label": "NAD^+"},
     ],
     "tags": ["TCA", "beta_oxidation"],
     "notes": "同一組三步化學：–CH₂–CH₂– 經 FAD 脫氫成 trans 雙鍵 → 加水成 β-OH（L 型）→ NAD⁺ 脫氫成 β-酮。"
              "TCA 第 6–8 步（SDH、fumarase、MDH）＝β 氧化的 acyl-CoA DH、enoyl-CoA hydratase、3-hydroxyacyl-CoA DH ＝支鏈胺基酸（Ile、Val、Leu）分解中段。"
              "為什麼 C–C 脫氫用 FAD？烷→烯的 ΔE°′ 不足以還原 NAD⁺；醇→酮才用 NAD⁺。SDH 是 ETC complex II（唯一嵌在內膜的 TCA 酵素）。"},
    {"id": "glycogen_enzyme_phosphorylation", "title": "Glycogen phosphorylase and synthase: reciprocal phosphorylation",
     "layout": "free", "cell_size": 4.6,
     "nodes": [
         t("gpb", "Phosphorylase b\n(less active)", 0, 0, "protein"), t("gpa", "Phosphorylase a–Ⓟ\n(active)", 2.8, 0, "protein"),
         t("gsa", "Glycogen synthase a\n(active)", 0, 2.6, "protein"), t("gsb", "Glycogen synthase b–Ⓟ\n(inactive)", 2.8, 2.6, "protein"),
         t("glu", "Glucagon / epinephrine → cAMP → PKA;   insulin → PP1", 1.4, 4.0, "concept"),
         t("ampb", "+ AMP (muscle)", -1.6, 0, "name"), t("glca", "− glucose (liver),\n− ATP, G6P", 4.5, 0, "name"),
         t("g6pb", "+ G6P", 4.4, 2.6, "name"),
     ],
     "edges": [
         {"from": "gpb", "to": "gpa", "bend": 0.4, "label": "phosphorylase kinase\n(PKA; Ca^{2+} in muscle)", "side_text": "ATP → ADP"},
         {"from": "gpa", "to": "gpb", "bend": 0.4, "label": "PP1", "side_text": "→ P_i", "color": GREEN},
         {"from": "gsa", "to": "gsb", "bend": 0.4, "label": "PKA, GSK3", "side_text": "ATP → ADP"},
         {"from": "gsb", "to": "gsa", "bend": 0.4, "label": "PP1", "side_text": "→ P_i", "color": GREEN},
         {"from": "ampb", "to": "gpb"}, {"from": "glca", "to": "gpa", "inhibit": True},
         {"from": "g6pb", "to": "gsb"},
     ],
     "tags": ["glycogen", "integration"],
     "notes": "同一個訊號相反作用：glucagon（肝）／epinephrine（肝、肌）→ cAMP → PKA → phosphorylase kinase → 磷酸化 phosphorylase（b → a，活化）；"
              "PKA 與 GSK3 同時磷酸化 glycogen synthase（a → b，失活）→ 分解開、合成關。"
              "Insulin → PP1 去磷酸化兩者 → 合成開、分解關（insulin 也經 Akt 抑制 GSK3）。"
              "別構：AMP 活化肌肉 phosphorylase b（運動時不需激素）；肝 phosphorylase a 被 glucose 抑制（血糖感測）；G6P 活化 synthase b。"
              "記法：磷酸化 → 分解酵素活化、合成酵素失活（ACC、PFK-2 亦同）。"},
    {"id": "xenobiotic_detoxification", "title": "Xenobiotic metabolism: phase I and II", "layout": "free", "cell_size": 4.6,
     "nodes": [
         t("rh", "R–H\n(nonpolar drug,\nsteroid, toxin)", 0, 0, "name"), t("roh", "R–OH\n(hydroxylated)", 2.4, 0, "name"),
         t("conj", "R–O–glucuronide\n(polar)", 4.8, 0, "name"), t("urine", "Urine / bile", 6.8, 0, "name"),
         n("udpga", "udp_glucuronate", 3.6, 1.2),
         t("p450", "Cytochrome P450 monooxygenase (liver smooth ER)\nR–H + O_2 + NADPH + H^+ → R–OH + H_2O + NADP^+",
           2.4, -1.4, "concept"),
         t("other", "Other conjugates: sulfate (PAPS), glutathione,\nglycine, acetyl (NAT), methyl (SAM)", 5.0, 2.4, "name"),
         t("ox", "Oxidase:  RH + O_2 → R + H_2O\nMonooxygenase:  one O → R–OH, one O → H_2O\nDioxygenase:  both O atoms → product",
           0.4, 2.7, "plain"),
     ],
     "edges": [
         {"from": "rh", "to": "roh", "label": "Phase I\n(CYP450)", "side_text": "O_2, NADPH"},
         {"from": "roh", "to": "conj", "label": "Phase II\n(UGT)"},
         {"from": "udpga", "to": "edge:roh->conj"},
         {"from": "conj", "to": "urine"},
     ],
     "tags": ["detoxification"],
     "notes": "Phase I（官能基化）：主要是 cytochrome P450 單氧酶（mixed-function oxidase，肝 smooth ER）：O₂ 的一個 O 進入受質成 –OH、另一個成 H₂O，電子由 NADPH 經 cytochrome P450 reductase（FAD、FMN）供應。"
              "Phase II（結合）：UDP-glucuronosyltransferase（UGT，UDP-glucuronate）、sulfotransferase（PAPS）、GST（glutathione）、N-acetyltransferase、methyltransferase → 水溶性增加 → 尿或膽汁排出。"
              "同類單氧酶：類固醇合成的 P450（CYP11A1、21、11B1、17）、7α-hydroxylase（膽酸）、desaturase（配 cyt b₅）、NOS。"
              "CYP 誘導：酒精（CYP2E1）、phenobarbital、rifampin；抑制：葡萄柚汁、ketoconazole。新生兒 UGT 不足 → 生理性黃疸；Crigler–Najjar、Gilbert 為 UGT1A1 缺陷。"},
    {"id": "futile_cycle", "title": "Futile (substrate) cycles", "layout": "free", "cell_size": 4.6,
     "nodes": [
         n("glc", "glucose", 0, 0), n("g6p", "glucose_6p", 2.6, 0),
         n("f6p", "fructose_6p", 0, 1.8), n("fbp", "fructose_16bp", 2.6, 1.8),
         t("net", "Net: ATP + H_2O → ADP + P_i + heat", 1.3, 2.95, "concept"),
     ],
     "edges": [
         {"from": "glc", "to": "g6p", "step": "hexokinase", "bend": 0.35, "label": "hexokinase / GK", "side_text": "ATP → ADP"},
         {"from": "g6p", "to": "glc", "step": "g6pase", "bend": 0.35, "label": "G6Pase", "side_text": "H_2O → P_i"},
         {"from": "f6p", "to": "fbp", "step": "pfk1", "bend": 0.35, "label": "PFK-1", "side_text": "ATP → ADP"},
         {"from": "fbp", "to": "f6p", "step": "fbpase1", "bend": 0.35, "label": "FBPase-1", "side_text": "H_2O → P_i"},
     ],
     "tags": ["glycolysis", "gluconeogenesis", "bioenergetics"],
     "notes": "同一對受質由兩個不同的不可逆酵素正反向轉換，若同時進行，淨結果只是水解 ATP 產熱（substrate cycle）。"
              "用途：① 產熱（熊蜂飛行肌 PFK-1／FBPase-1；malignant hyperthermia 為另一種 Ca²⁺ 循環）"
              "② 放大訊號：兩方向同時被反向調控（F2,6BP 活化 PFK-1、抑制 FBPase-1），淨流量的變化比單一酵素變化大很多（筆記例：+20 % 與 −30 % → 淨流量變化 98 %）。"
              "正常情況下由 insulin／glucagon 互惠調控，避免同時進行。"},
    {"id": "atp_cleavage_modes", "title": "Two ways to cleave ATP", "layout": "free", "cell_size": 4.6,
     "nodes": [
         {"id": "atp", "ref": "atp_n", "pos": [1.2, 0]},
         t("bg", "β–γ cleavage\nATP → ADP + P_i\n(1 ~Ⓟ, ΔG°' −7.3 kcal/mol)", -0.6, 2.5, "concept"),
         t("ab", "α–β cleavage\nATP → AMP + PP_i\n→ PP_i + H_2O → 2 P_i\n(2 ~Ⓟ, ΔG°' ≈ −15 kcal/mol)", 3.0, 2.5, "concept"),
         t("bg1", "Kinases (hexokinase, PFK-1),\nATPases, active transport", -0.6, 3.8, "name"),
         t("ab1", "Aminoacyl-tRNA synthetase,\nacyl-CoA synthetase,\nUDP-Glc pyrophosphorylase (UTP),\nargininosuccinate synthetase,\nDNA / RNA polymerase", 3.0, 4.1, "name"),
     ],
     "edges": [
         {"from": "atp", "to": "bg"}, {"from": "atp", "to": "ab"},
         {"from": "bg", "to": "bg1", "head": False}, {"from": "ab", "to": "ab1", "head": False},
     ],
     "tags": ["bioenergetics"],
     "notes": "β–γ 斷裂（→ ADP + Pi）用掉 1 個高能鍵：激酶、ATPase、主動運輸。"
              "α–β 斷裂（→ AMP + PPi）再由 pyrophosphatase 水解 PPi → 2 Pi，相當於消耗 2 個高能鍵，讓反應強力往右（不可逆）："
              "胺基酸活化（aminoacyl-tRNA）、脂肪酸活化（acyl-CoA，故 β 氧化計算要扣 2 ATP）、UDP-glucose 合成（UTP → PPi）、尿素循環 ASS、核酸聚合。"
              "其他 NTP：CTP → 磷脂（CDP-choline、CDP-ethanolamine）；GTP → 蛋白質合成（IF2、EF-Tu、EF-G、RF3）、糖質新生 PEPCK；UTP → 醣（UDP-Glc）。"
              "高能可識別記號：UDP（醣）、CoA（醯基）、SAM（甲基）。"},
]


# ================================================================ writer

def dump(path, d):
    """Write the YAML, keeping a review status already given in the web UI (re-running a generator must not reset it)."""
    if path.exists() and "status" in d:
        old = (yaml.safe_load(path.read_text(encoding="utf-8")) or {}).get("status")
        if old in ("approved", "changes_requested"):
            d["status"] = old
    path.write_text(yaml.safe_dump(d, allow_unicode=True, sort_keys=False, width=200), encoding="utf-8")


PARENT = {"glycolysis": "carb_metabolism", "gluconeogenesis": "carb_metabolism", "glycogen": "carb_metabolism",
          "PPP": "carb_metabolism", "galactose_metabolism": "carb_metabolism", "fructose_metabolism": "carb_metabolism",
          "beta_oxidation": "lipid_metabolism", "bioenergetics": "energy_metabolism", "detoxification": "energy_metabolism",
          "integration": "energy_metabolism", "ethanol": "energy_metabolism", "enzyme_kinetics": "enzyme",
          "sugar_derivative": "carbohydrate", "disaccharide": "carbohydrate"}


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
        lines = ["", "# --- 第十六輪（代謝筆記補齊）：2 NAD⁺／2 NADH（由 tools/gen_met.py 產生）"]
        for k in new:
            v = NEW_COFACTORS[k]
            extra = f', molecule: {v["molecule"]}' if v.get("molecule") else ""
            lines.append(f'{k}: {{label: "{v["label"]}", smiles: "{v["smiles"]}"{extra}}}')
        p.write_text(text.rstrip("\n") + "\n" + "\n".join(lines) + "\n", encoding="utf-8")


def main(data):
    data = Path(data)
    add_cofactors(data)
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
