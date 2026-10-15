"""Energy metabolism (round 11): electron transport chain and oxidative phosphorylation,
NADH shuttles, ATP yield, ethanol / methanol metabolism and the integration of metabolism
(fed / fasting state, insulin vs glucagon, metabolic crossroads).

    .venv/bin/python tools/gen_energy.py data

Complexes I–V are written as their overall chemical equations (protons pumped across the
membrane are not part of the balance; they are given as `conditions`). Cytochrome c is a
cofactor "cyt c (Fe³⁺/Fe²⁺)" whose [15*] stands for the protein. Ubiquinone is drawn as a
schematic `free` drawing with R = the isoprenoid tail. Re-running overwrites the files listed
here, adds the cofactors below to data/cofactors.yaml and a `chain` to acetaldehyde.
"""
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent))
from reaction_colors import COLORS  # noqa: E402

ROOT_TAG = "energy_metabolism"
KIDS = ["etc_oxphos", "shuttle", "integration", "ethanol"]

# ---- ubiquinone / ubiquinol (Lehninger: methoxy groups left, CH3 and the tail right)
HEX = {"c1": (0, -0.82), "c2": (0.71, -0.41), "c3": (0.71, 0.41), "c4": (0, 0.82), "c5": (-0.71, 0.41), "c6": (-0.71, -0.41)}


def coq(reduced):
    atoms = {k: ["", x, y] for k, (x, y) in HEX.items()}
    o = "OH" if reduced else "O"
    atoms.update({"o1": [o, 0, -1.62], "o4": [o, 0, 1.62],
                  "m6": ["H_3CO", -1.5, -0.86, 3], "m5": ["H_3CO", -1.5, 0.86, 3],
                  "me": ["CH_3", 1.5, -0.86], "r": ["R", 1.5, 0.86]})
    ring = [("c1", "c2"), ("c2", "c3"), ("c3", "c4"), ("c4", "c5"), ("c5", "c6"), ("c6", "c1")]
    dbl = {("c2", "c3"), ("c5", "c6")} if not reduced else {("c1", "c2"), ("c3", "c4"), ("c5", "c6")}
    bonds = [[a, b, 2] if (a, b) in dbl else [a, b] for a, b in ring]
    bonds += [["c1", "o1"] + ([] if reduced else [2]), ["c4", "o4"] + ([] if reduced else [2]),
              ["c6", "m6"], ["c5", "m5"], ["c2", "me"], ["c3", "r"]]
    return {"atoms": atoms, "bonds": bonds, "schematic": True,
            "notes": [["R = (CH_2CH=C(CH_3)CH_2)_{10}H", 0, 2.55, "gray"]]}


def molecules():
    return [
        dict(id="ubiquinone", name="Ubiquinone (Q)", label="Ubiquinone (Q)", aliases=["coenzyme Q", "CoQ10", "Q", "輔酶 Q"],
             smiles="COC1=C(OC)C(=O)C(C[10*])=C(C)C1=O", free=coq(False), tags=["etc_oxphos"],
             notes="脂溶性的電子載體，在內膜中自由擴散。收集 Complex I、II、mitochondrial G3P DH、ETF:Q oxidoreductase（脂肪酸 β 氧化）、DHODH 的電子，交給 Complex III。可一次收 1 個電子（semiquinone •Q⁻）或 2 個電子＋2 H⁺（QH₂）。statin 也會減少 CoQ 合成（共用 mevalonate 路徑）。圖中 R 是類異戊二烯尾巴（人類 10 個單位），以示意圖繪製。"),
        dict(id="ubiquinol", name="Ubiquinol (QH2)", label="Ubiquinol (QH_2)", aliases=["QH2", "reduced coenzyme Q", "還原態輔酶 Q"],
             smiles="COc1c(O)c(C)c(C[10*])c(O)c1OC", free=coq(True), tags=["etc_oxphos"],
             notes="Q 接受 2 e⁻ + 2 H⁺ 後的還原態。Complex III 以 Q cycle 把 QH₂ 的兩個電子分開：一個經 Rieske Fe-S → cyt c₁ → cyt c，另一個經 cyt b 回到 Q。"),
        dict(id="adp_n", name="ADP", aliases=["adenosine diphosphate"],
             smiles="Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]1O",
             haworth={"ring": "furanose", "subs": {1: {"up": "base", "down": "H"}, 2: {"up": "H", "down": "OH"},
                                                   3: {"up": "H", "down": "OH"}, 4: {"up": "{P}–{P}–OCH_2", "down": "H"}},
                      "base": "[*]n1cnc2c(N)ncnc21", "base_right": 5, "numbers": True},
             tags=["etc_oxphos", "nucleotide", "purine"],
             notes="ATP synthase 的受質（與 Pi）。ADP/ATP 經 adenine nucleotide translocase（ANT）一進一出（atractyloside、bongkrekic acid 抑制）。"),
        dict(id="oxygen", name="Oxygen", aliases=["O2", "dioxygen", "氧"], smiles="O=O", formula="O_2", tags=["etc_oxphos"],
             notes="電子傳遞鏈最終的電子接受者（Complex IV：O₂ + 4 e⁻ + 4 H⁺ → 2 H₂O）。"),
        # ---- ethanol / methanol (explicit H on the carbon that is oxidised, Mathews style)
        dict(id="ethanol", name="Ethanol", aliases=["alcohol", "EtOH", "乙醇"], smiles="CCO",
             chain=[{"c": "C", "l": "H", "r": "OH", "u": "H"}, {"c": "CH_3"}], tags=["ethanol"],
             notes="主要在肝代謝：ADH（細胞質，主要）、MEOS／CYP2E1（微粒體，慢性飲酒誘導）、catalase（過氧化體，少量）。"),
        dict(id="acetate", name="Acetate", aliases=["ethanoate", "醋酸根"], smiles="CC(=O)[O-]",
             chain=[{"c": "C", "fork": {"l": "O", "lb": 2, "r": "O^-"}}, {"c": "CH_3"}], tags=["ethanol"],
             notes="ALDH 的產物；大多離開肝臟，在肌肉、心臟等由 acetyl-CoA synthetase 活化成 acetyl-CoA。"),
        dict(id="methanol", name="Methanol", aliases=["MeOH", "甲醇"], smiles="CO",
             chain=[{"c": "C", "l": "H", "r": "OH", "u": "H", "d": "H"}], tags=["ethanol"],
             notes="經 ADH → formaldehyde → formate；formate 抑制 cytochrome c oxidase → 視神經損傷（失明）、代謝性酸中毒（高陰離子間隙）。治療：fomepizole（ADH 抑制劑）或乙醇（競爭 ADH）、血液透析、folate。"),
        dict(id="formaldehyde", name="Formaldehyde", aliases=["methanal", "甲醛"], smiles="C=O",
             chain=[{"c": "C", "fork": {"l": "O", "lb": 2, "r": "H"}, "d": "H"}], tags=["ethanol"],
             notes="甲醇中毒的中間物，很快被 aldehyde dehydrogenase 氧化成 formate。"),
        dict(id="formate", name="Formate", aliases=["methanoate", "甲酸根"], smiles="O=C[O-]",
             chain=[{"c": "C", "fork": {"l": "O", "lb": 2, "r": "O^-"}, "d": "H"}], tags=["ethanol"],
             notes="甲醇中毒的毒性產物：抑制 Complex IV（cytochrome c oxidase），視網膜與視神經最敏感。由 folate（N¹⁰-formyl-THF）途徑緩慢清除。"),
    ]


COFACTORS = {      # name: (label, SMILES, linked molecule)
    "Q": ("Q", "O=C1C=CC(=O)C([10*])=C1", "ubiquinone"),
    "QH2": ("QH_2", "Oc1ccc(O)c([10*])c1", "ubiquinol"),
}
# [15*] = the cytochrome c protein around the heme iron
MULTI = {"2cytc_ox": ("2 cyt c (Fe^{3+})", "[15*][Fe+3]", 2), "2cytc_red": ("2 cyt c (Fe^{2+})", "[15*][Fe+2]", 2),
         "4cytc_ox": ("4 cyt c (Fe^{3+})", "[15*][Fe+3]", 4), "4cytc_red": ("4 cyt c (Fe^{2+})", "[15*][Fe+2]", 4),
         "2H+": ("2 H^+", "[H+]", 2), "4H+": ("4 H^+", "[H+]", 4)}
MARK = "# --- 能量代謝（第十一輪）"

REACTIONS = [
    dict(id="complex_i", title="Complex I: NADH + Q → NAD⁺ + QH₂ (NADH:ubiquinone oxidoreductase)",
         substrates=["ubiquinone"], products=["ubiquinol"], cofactors_in=["NADH+H+"], cofactors_out=["NAD+"],
         prosthetic=["FMN,", "Fe-S"], conditions="4 H^+ pumped (matrix → IMS)",
         enzyme="Complex I\n(NADH dehydrogenase)",
         card={"location": "粒線體內膜", "type": "氧化還原＋質子幫浦（每 2 e⁻ 打出 4 H⁺）",
               "cofactors": "FMN → 一串 Fe-S cluster → Q",
               "deficiency": ["rotenone（殺蟲劑、魚藤酮）、amytal（amobarbital）、piericidin A 抑制 → NADH 無法氧化，但 succinate 仍可經 Complex II 進來",
                              "Leber hereditary optic neuropathy（LHON）：粒線體 DNA 的 Complex I 基因突變，母系遺傳"],
               "exam_points": ["NADH 的 P/O ≈ 2.5（10 H⁺ / 4 H⁺ per ATP）", "Complex I、III、IV 是三個質子幫浦，II 不是"]},
         tags=["etc_oxphos"]),
    dict(id="complex_ii", title="Complex II: Succinate + Q → Fumarate + QH₂ (succinate dehydrogenase)",
         substrates=[{"id": "succinate", "chain": "explicit"}, "ubiquinone"], products=["fumarate", "ubiquinol"], prosthetic=["FAD,", "Fe-S"],
         conditions="no H^+ pumped", enzyme="Complex II\n(succinate dehydrogenase)",
         card={"location": "粒線體內膜（TCA 唯一的膜結合酵素）", "type": "氧化（FAD → Fe-S → Q）",
               "deficiency": ["malonate：競爭性抑制（結構類似 succinate）", "carboxin、TTFA 抑制 Q 結合位"],
               "exam_points": ["不打質子 → FADH₂ 的 P/O ≈ 1.5（6 H⁺）", "也是 TCA 的第 6 步"]},
         tags=["etc_oxphos"]),
    dict(id="complex_iii", title="Complex III: QH₂ + 2 cyt c (Fe³⁺) → Q + 2 cyt c (Fe²⁺) + 2 H⁺",
         substrates=["ubiquinol"], products=["ubiquinone"], cofactors_in=["2cytc_ox"], cofactors_out=["2cytc_red", "2H+"],
         prosthetic=["cyt b, Fe-S,", "cyt c_1"], conditions="Q cycle: 4 H^+ to IMS per 2 e^-",
         enzyme="Complex III\n(cytochrome bc_1)",
         card={"location": "粒線體內膜；cyt c 在膜間腔側", "type": "氧化還原＋質子轉移（Q cycle）",
               "mechanism": "Q cycle：QH₂ 在外側（Qo）放出 2 H⁺ 到膜間腔；一個 e⁻ 經 Rieske Fe-S → cyt c₁ → cyt c，另一個經 cyt b_L → b_H 到內側（Qi）把 Q 還原。兩輪合計：2 QH₂ 氧化、1 Q 還原，4 H⁺ 到膜間腔、2 個 cyt c 被還原",
               "deficiency": ["antimycin A 抑制 Qi 位；myxothiazol 抑制 Qo 位"],
               "exam_points": ["把 2 e⁻ 載體（QH₂）轉成 1 e⁻ 載體（cyt c）", "cyt c 是水溶性周邊蛋白；釋放到細胞質 → 活化 caspase-9（內在性凋亡）"]},
         tags=["etc_oxphos"]),
    dict(id="complex_iv", title="Complex IV: 4 cyt c (Fe²⁺) + O₂ + 4 H⁺ → 4 cyt c (Fe³⁺) + 2 H₂O",
         substrates=["oxygen"], products=[{"id": "water", "n": 2}], cofactors_in=["4cytc_red", "4H+"], cofactors_out=["4cytc_ox"],
         prosthetic=["Cu_A, heme a,", "heme a_3–Cu_B"], conditions="4 H^+ pumped per O_2",
         enzyme="Complex IV\n(cytochrome c oxidase)",
         card={"location": "粒線體內膜", "type": "O₂ 還原成水＋質子幫浦",
               "mechanism": "4 個 e⁻ 依序經 Cu_A → heme a → heme a₃–Cu_B 雙核中心，O₂ 在 a₃–Cu_B 一次接受 4 e⁻，不放出 ROS。另有 4 個「化學」H⁺ 從基質端進來生成水",
               "deficiency": ["CN⁻、N₃⁻（azide）、H₂S 結合 Fe³⁺（a₃）；CO 結合 Fe²⁺ → 細胞缺氧、乳酸中毒",
                              "氰化物中毒治療：hydroxocobalamin；或 nitrite（形成 methemoglobin 抓 CN⁻）＋thiosulfate（rhodanese → SCN⁻）"],
               "exam_points": ["含 heme 與 Cu", "基質端消耗 H⁺ 也讓梯度變大"]},
         tags=["etc_oxphos"]),
    dict(id="atp_synthase", title="Complex V: ADP + Pi + H⁺ → ATP + H₂O (ATP synthase)",
         substrates=["adp_n"], products=["atp_n"], cofactors_in=["Pi"], cofactors_out=["H2O"], extra_in=["H+"],
         conditions="≈ 3 H^+ flow IMS → matrix (F_o)", enzyme="ATP synthase\n(F_oF_1, Complex V)",
         card={"location": "粒線體內膜（F₁ 頭部朝基質）", "type": "化學滲透偶聯（chemiosmotic coupling）",
               "mechanism": "H⁺ 經 F_o（a 次單元與 c ring）流回基質 → c ring 與 γ 軸旋轉 → F₁ 的三個 β 次單元輪流經歷 open → loose → tight（binding change, Boyer）；ATP 的釋放才是耗能步驟",
               "deficiency": ["oligomycin 結合 F_o 阻斷 H⁺ 通道 → ATP 合成與電子傳遞都停止（呼吸控制）",
                              "uncoupler（2,4-DNP、大劑量 aspirin、FCCP、UCP1／thermogenin 在棕色脂肪）：H⁺ 不經 ATP synthase 回流 → O₂ 消耗 ↑、ATP ↓、產熱 ↑"],
               "exam_points": ["每個 ATP 約需 4 H⁺（3 個經 F_o ＋ 1 個隨 Pi 進基質）→ P/O：NADH 2.5、FADH₂ 1.5",
                               "呼吸控制：ADP 的供應決定電子傳遞的速率"]},
         tags=["etc_oxphos"]),
    dict(id="mito_glycerol_3p_dehydrogenase", title="Glycerol 3-phosphate + Q → DHAP + QH₂ (mitochondrial G3P dehydrogenase)",
         substrates=["glycerol_3p"], products=["dhap"], cofactors_in=["Q"], cofactors_out=["QH2"], prosthetic=["FAD"],
         enzyme="mitochondrial glycerol 3-phosphate\ndehydrogenase (GPD2)",
         card={"location": "粒線體內膜外側（面向膜間腔）", "type": "氧化（FAD → Q）",
               "exam_points": ["glycerol 3-phosphate shuttle 的第二步：細胞質 NADH 的電子以 FADH₂ 的形式進入 → 只得約 1.5 ATP",
                               "骨骼肌、腦常用；比 malate–aspartate shuttle 快但少 1 ATP"]},
         tags=["shuttle"]),
    dict(id="alcohol_dehydrogenase", title="Ethanol → Acetaldehyde (alcohol dehydrogenase)",
         substrates=["ethanol"], products=["acetaldehyde"], cofactors_in=["NAD+"], cofactors_out=["NADH+H+"],
         enzyme="alcohol dehydrogenase\n(ADH)", reversible=True,
         card={"location": "肝細胞細胞質（也在胃）", "type": "氧化（NAD⁺）",
               "regulation": "零級動力學（Km 很低，飲酒時已飽和）→ 每小時約代謝 7–10 g 乙醇",
               "deficiency": ["fomepizole 抑制 ADH：治療甲醇、乙二醇中毒", "女性胃 ADH 較少 → 同劑量血中濃度較高"],
               "exam_points": ["NADH/NAD⁺ ↑ 是酒精代謝異常的根源（見總圖）"]},
         tags=["ethanol"]),
    dict(id="aldehyde_dehydrogenase", title="Acetaldehyde → Acetate (aldehyde dehydrogenase, ALDH2)",
         substrates=["acetaldehyde"], products=["acetate"], cofactors_in=["NAD+", "H2O"], cofactors_out=["NADH+H+"],
         extra_out=["H+"], enzyme="aldehyde dehydrogenase\n(ALDH2, mitochondrial)",
         card={"location": "肝細胞粒線體（ALDH2）", "type": "氧化（NAD⁺）",
               "deficiency": ["ALDH2*2（Glu504Lys，東亞人常見）：乙醛堆積 → 臉紅、心悸、噁心（亞洲紅臉）；也增加食道癌風險",
                              "disulfiram 抑制 ALDH → 喝酒即出現類似反應（戒酒藥）；metronidazole、某些 cephalosporin、sulfonylurea 也會"],
               "exam_points": ["乙醛有毒：與蛋白質形成 adduct → 肝損傷"]},
         tags=["ethanol"]),
    dict(id="meos_cyp2e1", title="Ethanol → Acetaldehyde (MEOS, CYP2E1)",
         substrates=["ethanol"], products=["acetaldehyde"], cofactors_in=["NADPH+H+", "O2"], cofactors_out=["NADP+", "2H2O"],
         enzyme="CYP2E1 (MEOS)",
         card={"location": "肝細胞平滑內質網（微粒體）", "type": "單氧酶（cytochrome P450）",
               "regulation": "慢性飲酒誘導 CYP2E1 → 酒量變好（代謝耐受）",
               "deficiency": ["CYP2E1 也把 acetaminophen 變成 NAPQI → 酗酒者 acetaminophen 肝毒性 ↑", "產生 ROS → 肝損傷"],
               "exam_points": ["高濃度乙醇時才重要（Km 高）"]},
         tags=["ethanol"]),
    dict(id="acetyl_coa_synthetase", title="Acetate → Acetyl-CoA (acetyl-CoA synthetase)",
         substrates=["acetate"], products=["acetyl_coa"], cofactors_in=["CoA-SH", "ATP"], cofactors_out=["AMP", "PPi"],
         enzyme="acetyl-CoA synthetase\n(ACSS)",
         card={"location": "肝外組織（肌肉、心臟、腦）的細胞質與粒線體", "type": "活化（經 acetyl-AMP 中間物，耗 2 個高能鍵）",
               "exam_points": ["與脂肪酸的 acyl-CoA synthetase 同樣產生 AMP + PPi"]},
         tags=["ethanol"]),
    dict(id="methanol_dehydrogenase", title="Methanol → Formaldehyde (alcohol dehydrogenase)",
         substrates=["methanol"], products=["formaldehyde"], cofactors_in=["NAD+"], cofactors_out=["NADH+H+"],
         enzyme="alcohol dehydrogenase\n(ADH)",
         card={"location": "肝", "type": "氧化（NAD⁺）",
               "deficiency": ["治療：fomepizole（或乙醇）競爭 ADH，讓甲醇由腎與透析排除"],
               "exam_points": ["乙二醇同樣經 ADH → glycolaldehyde → glycolate → oxalate（草酸鈣結晶 → 急性腎損傷）"]},
         tags=["ethanol"]),
    dict(id="formaldehyde_dehydrogenase", title="Formaldehyde → Formate (aldehyde dehydrogenase)",
         substrates=["formaldehyde"], products=["formate"], cofactors_in=["NAD+", "H2O"], cofactors_out=["NADH+H+"],
         extra_out=["H+"], enzyme="aldehyde dehydrogenase",
         card={"location": "肝", "type": "氧化（NAD⁺）",
               "deficiency": ["formate 抑制 cytochrome c oxidase → 視神經病變、失明；陰離子間隙代謝性酸中毒"]},
         tags=["ethanol"]),
]


def n(id_, ref, x, y, **kw):
    return {"id": id_, "ref": ref, "pos": [x, y], **kw}


def t(id_, text, x, y, kind="concept", **kw):
    return {"id": id_, "text": text, "kind": kind, "pos": [x, y], **kw}


def redox_ladder():
    """Carriers on a staircase (E°' rises to the lower right); dashed boxes group the complexes."""
    seq = [("nadh", "NAD^+/NADH\n−0.32 V", None), ("fmn", "FMN", "I"), ("fes1", "Fe-S", "I"),
           ("q", "Q/QH_2\n+0.05 V", None), ("cytb", "cyt b", "III"), ("fes3", "Fe-S (Rieske)", "III"),
           ("c1", "cyt c_1\n+0.22 V", "III"), ("c", "cyt c\n+0.25 V", None), ("cua", "Cu_A", "IV"),
           ("a", "cyt a\n+0.29 V", "IV"), ("a3", "cyt a_3–Cu_B\n+0.35 V", "IV"), ("o2", "O_2/H_2O\n+0.82 V", None)]
    DX, DY = 2.0, 1.0
    nodes, edges, groups = [], [], {}
    for k, (id_, text, cx) in enumerate(seq):
        x, y = k * DX, k * DY
        nodes.append(t(id_, text, x, y, "carrier"))
        if cx:
            groups.setdefault(cx, []).append((x, y))
        if k:
            edges.append({"from": seq[k - 1][0], "to": id_, "color": E})
    regions = [{"kind": "box", "x0": min(x for x, _ in pts) - 0.6, "y0": min(y for _, y in pts) - 0.75,
                "x1": max(x for x, _ in pts) + 1.1, "y1": max(y for _, y in pts) + 0.45, "label": f"Complex {cx}"}
               for cx, pts in groups.items()]
    nodes += [t("succ", "Succinate → FAD\n(Complex II) +0.03 V", 2.6, 5.4, "carrier"),
              t("note", "electrons flow toward higher (more positive) E°'\nΔG°' = −nFΔE°' ≈ −220 kJ/mol per NADH", 15.0, 3.2, "plain")]
    edges.append({"from": "succ", "to": "q", "color": E})
    return dict(id="etc_redox_carriers", title="Redox carriers in order of E°′（電子載體與還原電位）", layout="free", cell_size=3.6,
                regions=regions, nodes=nodes, edges=edges, tags=["etc_oxphos"],
                notes="E°′ 由負到正排列：NADH −0.32 → Q +0.05 → cyt c₁ +0.22 → cyt c +0.25 → cyt a +0.29 → a₃ +0.35 → O₂ +0.82 V（Lehninger Table 19-2 的近似值）。cytochrome 的順序 b → c₁ → c → a → a₃ 是常考點。Q 與 cyt c 是可移動載體；Fe-S 與 cytochrome 一次只傳 1 個電子，FMN 與 Q 可以 1 或 2 個。")


NAME = {"view": "name"}
E, H = "orange", "magenta"              # electron flow, proton flow
PROT = {"h": 4.2}                       # protein complexes span the membrane


PATHWAYS = [
    dict(id="electron_transport_chain", title="Electron transport chain & oxidative phosphorylation（電子傳遞鏈與氧化磷酸化）",
         layout="free", cell_size=3.4,
         regions=[{"kind": "membrane", "y": 0, "x0": -7.4, "x1": 10.4, "label": "Inner\nmitochondrial\nmembrane"}],
         nodes=[t("ims", "Intermembrane space (high [H^+], +)", -7.4, -3.0, "heading"),
                t("mat", "Matrix (low [H^+], −)", -7.4, 4.0, "heading"),
                t("c1", "I", -5.2, 0, "protein", w=2.6, **PROT), t("c2", "II", -2.6, 0.55, "protein", w=2.0, h=2.4),
                t("q", "Q", -1.1, -0.25, "carrier"), t("c3", "III", 0.6, 0, "protein", w=2.4, **PROT),
                t("cytc", "cyt c", 2.0, -1.75, "carrier"), t("c4", "IV", 3.5, 0, "protein", w=2.4, **PROT),
                t("c5", "ATP\nsynthase\n(F_oF_1)", 7.0, 0, "protein", w=3.0, **PROT),
                t("nadh", "NADH + H^+", -6.4, 2.4, "ion"), t("nad", "NAD^+", -4.1, 2.4, "ion"),
                t("succ", "Succinate", -3.5, 3.2, "name"), t("fum", "Fumarate", -1.5, 3.2, "name"),
                t("o2", "½O_2 + 2H^+", 2.4, 2.4, "ion"), t("h2o", "H_2O", 4.5, 2.4, "ion"),
                t("adp", "ADP + P_i", 6.0, 2.6, "ion"), t("atp", "ATP", 8.1, 2.6, "ion"),
                t("h1", "4H^+", -5.2, -2.4, "ion"), t("h3", "4H^+", 0.6, -2.4, "ion"), t("h4", "2H^+", 3.5, -2.4, "ion"),
                t("h5", "H^+", 7.0, -2.4, "ion"), t("h6", "H^+", 9.4, -2.4, "ion"), t("h7", "H^+", 9.4, 2.6, "ion"),
                t("rot", "Rotenone,\namytal", -7.2, -1.6, "drug"), t("mal", "Malonate", -2.6, -1.9, "drug"),
                t("anti", "Antimycin A", 1.4, 3.4, "drug"), t("cn", "CN^-, CO, N_3^-", 4.8, -1.8, "drug"),
                t("oligo", "Oligomycin", 5.6, -3.3, "drug"),
                t("unc", "Uncouplers\n(2,4-DNP, UCP1)", 10.9, 2.0, "drug")],
         edges=[{"from": "nadh", "to": "c1", "step": "complex_i"}, {"from": "c1", "to": "nad", "step": "complex_i"},
                {"from": "succ", "to": "c2", "step": "complex_ii"}, {"from": "c2", "to": "fum", "step": "complex_ii"},
                {"from": "c1", "to": "q", "color": E, "bend": -0.2}, {"from": "c2", "to": "q", "color": E},
                {"from": "q", "to": "c3", "color": E}, {"from": "c3", "to": "cytc", "color": E},
                {"from": "cytc", "to": "c4", "color": E},
                {"from": "o2", "to": "c4", "step": "complex_iv"}, {"from": "c4", "to": "h2o", "step": "complex_iv"},
                {"from": "c1", "to": "h1", "color": H}, {"from": "c3", "to": "h3", "color": H}, {"from": "c4", "to": "h4", "color": H},
                {"from": "h5", "to": "c5", "color": H},
                {"from": "adp", "to": "c5", "step": "atp_synthase"}, {"from": "c5", "to": "atp", "step": "atp_synthase"},
                {"from": "rot", "to": "c1", "inhibit": True}, {"from": "mal", "to": "c2", "inhibit": True},
                {"from": "anti", "to": "c3", "inhibit": True}, {"from": "cn", "to": "c4", "inhibit": True},
                {"from": "oligo", "to": "c5", "inhibit": True},
                {"from": "h6", "to": "h7", "color": H}],
         tags=["etc_oxphos"],
         notes="橙色箭頭＝電子流（NADH → I → Q → III → cyt c → IV → O₂；succinate → II → Q）；洋紅箭頭＝質子流。每對電子：NADH 打出 10 H⁺（I 4、III 4、IV 2），FADH₂ 6 H⁺。ATP synthase 約 4 H⁺/ATP（含 Pi 運輸）→ P/O：NADH 2.5、FADH₂ 1.5。抑制劑：I rotenone、amytal；II malonate；III antimycin A；IV CN⁻、CO、N₃⁻、H₂S；V oligomycin；ANT atractyloside。uncoupler 讓 H⁺ 繞過 ATP synthase：O₂ 消耗 ↑、ATP ↓、體溫 ↑（UCP1 = thermogenin，新生兒棕色脂肪）。"),
    redox_ladder(),
    dict(id="atp_synthase_fof1", title="ATP synthase F_oF_1: structure and binding change（ATP 合成酶 F_oF_1）",
         layout="free", cell_size=3.4,
         regions=[{"kind": "membrane", "y": 0, "x0": -6.4, "x1": 11.6, "label": "Inner\nmitochondrial\nmembrane"},
                  {"kind": "box", "x0": -2.8, "y0": -1.1, "x1": 2.4, "y1": 1.1, "label": "F_o"},
                  {"kind": "box", "x0": -0.7, "y0": 1.5, "x1": 2.5, "y1": 6.5, "label": "F_1"}],
         nodes=[t("ims", "Intermembrane space (P side, high [H^+])", -6.4, -3.0, "heading"),
                t("mat", "Matrix (N side)", -6.4, 6.6, "heading"),
                t("hin", "H^+", -1.6, -2.4, "ion"),
                t("a", "a\n(H^+ half-\nchannels)", -1.6, 0, "protein", w=2.4, h=4.4),
                t("c", "c ring\n(8 c subunits)", 0.9, 0, "protein", w=3.2, h=3.4),
                t("hout", "H^+", 3.0, 1.8, "ion"),
                t("g", "γ ε\n(central stalk)", 0.9, 2.4, "protein", w=2.6, h=1.8),
                t("f1", "α_3β_3\n(F_1 head, 3 catalytic β)", 0.9, 4.3, "protein", w=5.6, h=2.4),
                t("d", "δ (OSCP)", 0.9, 5.9, "protein", w=2.8, h=1.2),
                t("b", "b_2\n(peripheral\nstalk)", -3.6, 3.2, "protein", w=2.0, h=4.0),
                t("adp", "ADP + P_i", 4.6, 3.6, "ion"), t("atp", "ATP", 4.6, 5.0, "ion"),
                t("oligo", "Oligomycin\n(binds F_o)", -4.4, -1.9, "drug"),
                t("rotor", "Rotor: c ring + γ + ε\nStator: a + b_2 + δ + α_3β_3", -3.6, 8.0, "concept"),
                t("ant", "ANT\n(ATP^{4-} out, ADP^{3-} in)", 7.0, 0, "protein", w=5.0, h=3.0),
                t("pic", "PiC\n(H^+ + P_i in)", 10.2, 0, "protein", w=3.0, h=3.0),
                t("atr", "Atractyloside,\nbongkrekic acid", 7.0, -2.3, "drug"),
                t("bc", "Binding change (Boyer): each β cycles O → L → T", 6.2, 2.6, "heading"),
                t("o", "O (open)\nATP leaves", 9.0, 3.6, "concept"),
                t("l", "L (loose)\nADP + P_i bind", 10.6, 5.8, "concept"),
                t("tt", "T (tight)\nATP forms (ΔG ≈ 0)", 7.4, 5.8, "concept"),
                t("note", "Energy is used to release ATP, not to make it", 9.0, 7.4, "plain")],
         edges=[{"from": "hin", "to": "a", "color": H}, {"from": "a", "to": "c", "color": H},
                {"from": "c", "to": "hout", "color": H},
                {"from": "c", "to": "g", "head": False}, {"from": "g", "to": "f1", "head": False},
                {"from": "a", "to": "b", "head": False}, {"from": "b", "to": "d", "head": False, "bend": -0.25},
                {"from": "d", "to": "f1", "head": False},
                {"from": "adp", "to": "f1"}, {"from": "f1", "to": "atp"},
                {"from": "oligo", "to": "a", "inhibit": True}, {"from": "atr", "to": "ant", "inhibit": True},
                {"from": "o", "to": "l", "label": "γ turns 120°"}, {"from": "l", "to": "tt", "label": "γ turns 120°"},
                {"from": "tt", "to": "o", "label": "γ turns 120°"}],
         tags=["etc_oxphos"],
         notes="F_o（o = oligomycin-sensitive）在膜內：H⁺ 由 a 次單元的膜間腔側半通道進入，與 c 次單元的 Asp/Glu 結合，c ring 轉約一圈後由基質側半通道放出。c ring、γ、ε 一起轉動（rotor）；a、b₂、δ、α₃β₃ 固定（stator），b₂ 與 δ 防止 α₃β₃ 跟著轉。γ 每轉 120°，三個 β 各自改變構形：O（空）→ L（鬆，ADP＋Pi 結合）→ T（緊，ATP 在活性位自發生成）→ O（釋放 ATP）。耗能的是 ATP 的釋放。哺乳類 c ring 有 8 個 c → 每轉一圈 8 H⁺ 產生 3 ATP（約 2.7 H⁺/ATP），加上 PiC 的 1 H⁺ → 約 3.7 ≈ 4 H⁺/ATP。ANT 由膜電位驅動（ATP⁴⁻ 出、ADP³⁻ 進），atractyloside 與 bongkrekic acid 抑制。"),
    dict(id="malate_aspartate_shuttle", title="Malate–aspartate shuttle（蘋果酸–天門冬胺酸穿梭）", layout="free", cell_size=3.6,
         regions=[{"kind": "membrane", "y": 0, "x0": -4.6, "x1": 4.6, "label": "Inner\nmitochondrial\nmembrane"}],
         nodes=[t("cy", "Cytosol", -4.6, -3.1, "heading"), t("mx", "Matrix", -4.6, 3.1, "heading"),
                n("asp_c", "aspartate", -3.2, -2.3, **NAME), n("oaa_c", "oxaloacetate", 0, -2.3, **NAME),
                n("mal_c", "malate", 3.2, -2.3, **NAME),
                t("t1", "malate–α-KG\nantiporter", 3.2, 0, "protein", w=3.0, h=2.2),
                t("t2", "glutamate–aspartate\nantiporter", -3.2, 0, "protein", w=3.6, h=2.2),
                n("mal_m", "malate", 3.2, 2.3, **NAME), n("oaa_m", "oxaloacetate", 0, 2.3, **NAME),
                n("asp_m", "aspartate", -3.2, 2.3, **NAME)],
         edges=[{"from": "oaa_c", "to": "mal_c", "step": "malate_dehydrogenase", "label": "malate DH",
                 "side_text": "NADH + H^+ → NAD^+", "reversible": True},
                {"from": "mal_c", "to": "t1"}, {"from": "t1", "to": "mal_m"},
                {"from": "mal_m", "to": "oaa_m", "step": "malate_dehydrogenase", "label": "malate DH",
                 "side_text": "NAD^+ → NADH + H^+", "reversible": True},
                {"from": "oaa_m", "to": "asp_m", "step": "ast", "label": "AST", "side_text": "Glu → α-KG", "reversible": True},
                {"from": "asp_m", "to": "t2"}, {"from": "t2", "to": "asp_c"},
                {"from": "asp_c", "to": "oaa_c", "step": "ast", "label": "AST", "side_text": "α-KG → Glu", "reversible": True}],
         tags=["shuttle"],
         notes="NADH 本身不能通過粒線體內膜。細胞質的 NADH 把 OAA 還原成 malate → malate–α-KG antiporter 進入基質 → 再氧化成 OAA，在基質內重新產生 NADH（→ Complex I，約 2.5 ATP）。OAA 也不能出膜，所以先轉胺成 Asp，經 Glu–Asp antiporter（aralar／citrin）出去。主要在肝、心、腎。citrin 缺乏 → 成人型 citrullinemia type II。可逆，依 NADH/NAD⁺ 比例決定方向。"),
    dict(id="glycerol_phosphate_shuttle", title="Glycerol 3-phosphate shuttle（甘油磷酸穿梭）", layout="free", cell_size=3.4,
         regions=[{"kind": "membrane", "y": 0, "x0": -3.8, "x1": 4.6, "label": "Inner\nmitochondrial\nmembrane"}],
         nodes=[t("cy", "Cytosol / intermembrane space", -3.8, -4.2, "heading"), t("mx", "Matrix", -3.8, 2.4, "heading"),
                n("dhap", "dhap", -2.2, -2.2, **NAME), n("g3p", "glycerol_3p", 2.2, -2.2, **NAME),
                t("gpd2", "mitochondrial\nG3P DH (FAD)", 0, -0.3, "protein", w=3.0, h=2.0),
                t("q", "Q → QH_2", 2.3, 0.2, "carrier"), t("c3", "III", 4.0, 0, "protein", w=1.6, h=3.0)],
         edges=[{"from": "dhap", "to": "g3p", "step": "glycerol_3p_dehydrogenase", "label": "cytosolic G3P DH (GPD1)",
                 "side_text": "NADH + H^+ → NAD^+", "bend": 0.35},
                {"from": "g3p", "to": "gpd2", "step": "mito_glycerol_3p_dehydrogenase"},
                {"from": "gpd2", "to": "dhap"},
                {"from": "gpd2", "to": "q", "color": E}, {"from": "q", "to": "c3", "color": E}],
         tags=["shuttle"],
         notes="細胞質 NADH 把 DHAP 還原成 glycerol 3-phosphate；位在內膜外側的 mitochondrial G3P DH（FAD）把它氧化回 DHAP，電子直接給 Q，跳過 Complex I → 每個細胞質 NADH 只得約 1.5 ATP。骨骼肌（白肌）、腦、棕色脂肪常用；不可逆（單向）。所以葡萄糖完全氧化：malate–aspartate shuttle 32 ATP，glycerol 3-P shuttle 30 ATP。"),
    dict(id="atp_yield_glucose", title="ATP yield from one glucose（一個葡萄糖完全氧化的 ATP）", layout="free", cell_size=3.2,
         nodes=[n("glc", "glucose", 0, 0, **NAME), n("pyr", "pyruvate", 0, 1.6, **NAME),
                n("acoa", "acetyl_coa", 0, 3.2, **NAME), t("tca", "Citric acid cycle", 0, 4.8, "concept"),
                t("co2", "CO_2", 3.2, 4.8, "ion"),
                t("y1", "2 ATP (net) + 2 NADH (cytosol)\n→ 2 + 3 or 5 = 5–7 ATP", 4.6, 0.8, "plain"),
                t("y2", "2 NADH → 5 ATP", 4.6, 2.4, "plain"),
                t("y3", "2 GTP + 6 NADH + 2 FADH_2\n→ 2 + 15 + 3 = 20 ATP", 4.6, 6.0, "plain"),
                t("tot", "Total ≈ 30–32 ATP\n(P/O: NADH 2.5, FADH_2 1.5;\ncytosolic NADH: 1.5 via G3P shuttle, 2.5 via Mal–Asp)", -4.4, 3.2, "concept")],
         edges=[{"from": "glc", "to": "pyr", "label": "glycolysis", "side_text": "×2 pyruvate"},
                {"from": "pyr", "to": "acoa", "step": "pdh", "label": "PDH", "side_text": "×2, 2 CO_2"},
                {"from": "acoa", "to": "tca"}, {"from": "tca", "to": "co2", "side_text": "4 CO_2"}],
         tags=["etc_oxphos", "integration"],
         notes="糖解淨 2 ATP＋2 NADH（細胞質）；PDH 2 NADH；TCA 2 GTP（ATP）＋6 NADH＋2 FADH₂。以 NADH 2.5、FADH₂ 1.5 計算：糖解 2＋(3 或 5)，PDH 5，TCA 2＋15＋3 → 30（G3P shuttle）或 32（Mal–Asp shuttle）。舊教科書的 36／38 是用 P/O 3 與 2 計算。無氧時只有糖解的 2 ATP（NAD⁺ 由 LDH 再生）。"),
    dict(id="ethanol_metabolism", title="Ethanol metabolism and its consequences（酒精代謝與後果）", layout="free", cell_size=3.4,
         nodes=[n("etoh", "ethanol", 0, 0, **NAME), n("acal", "acetaldehyde", 0, 1.8, **NAME),
                n("ac", "acetate", 0, 3.6, **NAME), n("acoa", "acetyl_coa", 0, 5.4, **NAME),
                t("fom", "Fomepizole", -2.6, 0.9, "drug"), t("dis", "Disulfiram", -2.6, 2.7, "drug"),
                t("flush", "ALDH2*2: flushing", 2.9, 2.7, "concept"),
                t("ratio", "NADH/NAD^+ ↑ (liver)", 3.9, 0.9, "concept"),
                t("lac", "Pyruvate → lactate:\nlactic acidosis, gout", 8.6, -0.7, "name"),
                t("gng", "OAA → malate:\ngluconeogenesis ↓ → hypoglycemia", 8.6, 0.9, "name"),
                t("fat", "β-oxidation ↓, G3P ↑ → TAG:\nfatty liver", 8.6, 2.5, "name"),
                t("keto", "acetyl-CoA → ketone bodies:\nalcoholic ketoacidosis", 3.9, 5.4, "name"),
                t("meoh", "Methanol → formaldehyde → formate\n(optic nerve, blindness)", -4.2, -1.0, "name"),
                t("eg", "Ethylene glycol → oxalate\n(kidney)", -4.2, -2.4, "name")],
         edges=[{"from": "etoh", "to": "acal", "step": "alcohol_dehydrogenase", "label": "ADH (cytosol);\nCYP2E1 (MEOS)",
                 "side_text": "NAD^+ → NADH"},
                {"from": "acal", "to": "ac", "step": "aldehyde_dehydrogenase", "label": "ALDH2 (mitochondria)",
                 "side_text": "NAD^+ → NADH"},
                {"from": "ac", "to": "acoa", "step": "acetyl_coa_synthetase", "label": "acetyl-CoA synthetase",
                 "side_text": "ATP → AMP + PP_i"},
                {"from": "fom", "to": "edge:etoh->acal", "inhibit": True}, {"from": "dis", "to": "edge:acal->ac", "inhibit": True},
                {"from": "flush", "to": "edge:acal->ac", "head": False},
                {"from": "ratio", "to": "lac"}, {"from": "ratio", "to": "gng"}, {"from": "ratio", "to": "fat"},
                {"from": "acoa", "to": "keto"}, {"from": "fom", "to": "meoh", "inhibit": True},
                {"from": "fom", "to": "eg", "inhibit": True}],
         tags=["ethanol"],
         notes="兩次脫氫都產生 NADH → 肝細胞 NADH/NAD⁺ 大增：pyruvate → lactate（乳酸中毒、乳酸與尿酸競爭排泄 → 痛風）；OAA → malate、pyruvate 被耗 → 糖質新生 ↓ → 空腹低血糖（尤其營養不良、肝醣耗竭時）；β 氧化與 TCA 被抑制、DHAP → G3P ↑ → TAG 堆積（脂肪肝）；acetyl-CoA 轉向酮體（酒精性酮酸中毒，以 β-hydroxybutyrate 為主）。甲醇（失明）與乙二醇（草酸 → 腎衰竭）都經 ADH 活化，以 fomepizole 阻斷。"),
    dict(id="fed_state", title="Fed (absorptive) state — insulin（進食狀態）", layout="free", cell_size=3.4,
         regions=[{"kind": "box", "x0": -6.6, "y0": -3.4, "x1": -1.2, "y1": 2.2, "label": "Liver"},
                  {"kind": "box", "x0": 2.6, "y0": -3.4, "x1": 7.8, "y1": -0.4, "label": "Muscle"},
                  {"kind": "box", "x0": 2.6, "y0": 0.2, "x1": 7.8, "y1": 3.6, "label": "Adipose tissue"},
                  {"kind": "box", "x0": -6.6, "y0": 2.8, "x1": -1.2, "y1": 4.6, "label": "Brain, RBC"}],
         nodes=[t("ins", "Insulin ↑, glucagon ↓", -6.6, -4.4, "heading"),
                t("bg", "Blood glucose\n(GLUT2 → liver;\nGLUT4 → muscle, fat)", 0.7, 0.4, "concept"),
                t("l1", "Glucose → glycogen", -3.9, -2.2, "name"), t("l2", "Glucose → acetyl-CoA\n→ fatty acids → TAG", -3.9, -0.8, "name"),
                t("l3", "TAG → VLDL", -3.9, 0.9, "name"),
                t("m1", "Glucose → glycogen", 5.2, -2.4, "name"), t("m2", "Amino acids → protein", 5.2, -1.3, "name"),
                t("a1", "Glucose → glycerol 3-P", 5.2, 1.2, "name"), t("a2", "Chylomicron / VLDL TAG\n(LPL) → FA → TAG", 5.2, 2.6, "name"),
                t("b1", "Glucose → CO_2 (brain); lactate (RBC)", -3.9, 3.9, "name")],
         edges=[{"from": "bg", "to": "l1"}, {"from": "bg", "to": "m1"}, {"from": "bg", "to": "a1"}, {"from": "bg", "to": "b1"},
                {"from": "l3", "to": "a2", "label": "VLDL", "bend": -0.15}],
         tags=["integration"],
         notes="胰島素（去磷酸化）：肝 glucokinase、PFK-2、pyruvate kinase、PDH、ACC、glycogen synthase、HMG-CoA reductase 活性 ↑；肌肉與脂肪 GLUT4 上膜；脂肪 LPL ↑、HSL ↓。合成作用為主：肝醣、脂肪酸、TAG、蛋白質。腦與 RBC 不依賴胰島素攝取葡萄糖（GLUT1/3）。"),
    dict(id="fasting_state", title="Fasting state — glucagon, epinephrine, cortisol（空腹／飢餓狀態）", layout="free", cell_size=3.4,
         regions=[{"kind": "box", "x0": -7.0, "y0": -3.5, "x1": -0.9, "y1": 2.4, "label": "Liver"},
                  {"kind": "box", "x0": 2.4, "y0": -3.5, "x1": 7.4, "y1": -0.5, "label": "Muscle"},
                  {"kind": "box", "x0": 2.4, "y0": 0.3, "x1": 7.4, "y1": 2.6, "label": "Adipose tissue"},
                  {"kind": "box", "x0": -7.0, "y0": 3.2, "x1": 7.4, "y1": 5.2, "label": "Brain, RBC"}],
         nodes=[t("glu", "Glucagon ↑, insulin ↓ (epinephrine, cortisol)", -7.0, -4.4, "heading"),
                t("l2", "Gluconeogenesis\n(lactate, alanine, glycerol)", -3.9, -2.4, "name"),
                t("l1", "Glycogenolysis (first ~24 h)", -3.9, -0.8, "name"),
                t("l3", "β-oxidation → ketone bodies", -3.9, 1.3, "name"),
                t("bg", "Blood glucose", 0.7, -0.3, "concept"),
                t("m1", "Proteolysis → Ala, Gln\n(glucose–alanine cycle)", 4.9, -2.6, "name"), t("m2", "Oxidises FA, ketone bodies", 4.9, -1.4, "name"),
                t("a1", "Lipolysis (HSL, ATGL)\n→ FA (on albumin) + glycerol", 4.9, 1.5, "name"),
                t("rbc", "RBC: glucose → lactate\n(Cori cycle)", -3.6, 4.3, "name"),
                t("br", "Brain: glucose; after days\nketone bodies (≈ 2/3)", 3.8, 4.3, "name")],
         edges=[{"from": "l1", "to": "bg"}, {"from": "l2", "to": "bg"},
                {"from": "m1", "to": "l2", "bend": 0.12},
                {"from": "a1", "to": "l2", "bend": -0.22},
                {"from": "a1", "to": "l3"},
                {"from": "l3", "to": "br"}, {"from": "bg", "to": "br"},
                {"from": "bg", "to": "rbc"}, {"from": "rbc", "to": "l2", "label": "lactate", "bend": 0.6}],
         tags=["integration"],
         notes="升糖素（cAMP/PKA 磷酸化）：肝 glycogen phosphorylase、FBPase-2（F2,6BP ↓ → 糖解 ↓、糖質新生 ↑）、PEPCK（轉錄）↑；pyruvate kinase L、ACC、glycogen synthase ↓。脂肪 HSL 由 epinephrine／低胰島素活化。肝不能用酮體（缺 SCOT）。長期飢餓：腦改用酮體 → 減少蛋白質分解（保留肌肉）；RBC 只能用葡萄糖。"),
    dict(id="insulin_glucagon_regulation", title="Insulin vs glucagon: key enzymes（胰島素與升糖素調控的酵素）", layout="free", cell_size=3.2,
         nodes=[t("ih", "Insulin (dephosphorylation)", -3.0, 0, "heading"), t("gh", "Glucagon / epinephrine (cAMP–PKA)", 3.0, 0, "heading"),
                t("i1", "Glycogen synthase ↑\nGlycogen phosphorylase ↓", -1.4, 1.3, "name"),
                t("g1", "Glycogen phosphorylase ↑\nGlycogen synthase ↓", 4.6, 1.3, "name"),
                t("i2", "PFK-2 ↑ (F2,6BP ↑) → glycolysis ↑\nPyruvate kinase (L) ↑", -1.4, 2.6, "name"),
                t("g2", "FBPase-2 ↑ (F2,6BP ↓) → gluconeogenesis ↑\nPyruvate kinase (L) ↓, PEPCK ↑", 4.6, 2.6, "name"),
                t("i3", "PDH ↑ (PDH phosphatase)\nACC ↑ → malonyl-CoA ↑", -1.4, 3.9, "name"),
                t("g3", "ACC ↓ (also AMPK)\nCPT-I ↑ → β-oxidation, ketogenesis", 4.6, 3.9, "name"),
                t("i4", "HMG-CoA reductase ↑\nLPL ↑, HSL ↓", -1.4, 5.2, "name"),
                t("g4", "HSL ↑ (lipolysis)\nHMG-CoA reductase ↓", 4.6, 5.2, "name")],
         edges=[], tags=["integration"],
         notes="原則：升糖素／腎上腺素經 cAMP–PKA 把酵素磷酸化；胰島素活化 protein phosphatase-1 去磷酸化。被磷酸化而「活化」的只有分解方向的酵素（glycogen phosphorylase、phosphorylase kinase、HSL、FBPase-2）；合成方向（glycogen synthase、ACC、PFK-2、PK-L、HMG-CoA reductase）被磷酸化則失活。"),
    dict(id="metabolic_crossroads", title="Metabolic crossroads: G6P, pyruvate, acetyl-CoA（代謝交叉點）", layout="free", cell_size=3.4,
         nodes=[n("g6p", "glucose_6p", 0, 0, **NAME), n("pyr", "pyruvate", 0, 2.6, **NAME), n("acoa", "acetyl_coa", 0, 5.2, **NAME),
                t("glyc", "Glycogen", -3.2, -1.0, "name"), t("glc", "Glucose (liver: G6Pase)", 3.4, -1.0, "name"),
                t("ppp", "Ribose 5-P + NADPH (PPP)", 3.6, 0.8, "name"),
                n("lac", "lactate", -3.2, 2.0, **NAME), n("ala", "alanine", -3.2, 3.3, **NAME),
                n("oaa", "oxaloacetate", 3.4, 2.6, **NAME),
                t("tca", "TCA cycle → CO_2", 3.4, 4.4, "name"), t("fa", "Fatty acids", -3.4, 4.6, "name"),
                t("chol", "Cholesterol", -3.4, 6.0, "name"), t("kb", "Ketone bodies", 3.4, 6.0, "name")],
         edges=[{"from": "glyc", "to": "g6p", "reversible": True}, {"from": "g6p", "to": "glc"}, {"from": "g6p", "to": "ppp"},
                {"from": "g6p", "to": "pyr", "label": "glycolysis", "reversible": True},
                {"from": "pyr", "to": "lac", "step": "ldh", "label": "LDH", "reversible": True},
                {"from": "pyr", "to": "ala", "step": "alt", "label": "ALT", "reversible": True},
                {"from": "pyr", "to": "oaa", "step": "pyruvate_carboxylase", "label": "pyruvate carboxylase"},
                {"from": "pyr", "to": "acoa", "step": "pdh", "label": "PDH (irreversible)"},
                {"from": "acoa", "to": "tca"}, {"from": "acoa", "to": "fa", "reversible": True},
                {"from": "acoa", "to": "chol"}, {"from": "acoa", "to": "kb"}],
         tags=["integration"],
         notes="G6P：肝醣、糖解、PPP、（肝、腎）葡萄糖輸出。Pyruvate：lactate（無氧、RBC）、alanine（肌肉 → 肝）、OAA（補充 TCA、糖質新生；acetyl-CoA 活化 pyruvate carboxylase）、acetyl-CoA（PDH，不可逆 → 脂肪酸不能變成葡萄糖）。Acetyl-CoA：TCA、脂肪酸、膽固醇、酮體。"),
]

ACETALDEHYDE_CHAIN = [{"c": "C", "fork": {"l": "O", "lb": 2, "r": "H"}}, {"c": "CH_3"}]


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


def with_parent(tags):
    out = list(tags)
    if any(t in KIDS for t in out) and ROOT_TAG not in out:
        out.insert(0, ROOT_TAG)
    return out


def add_cofactors(data):
    p = Path(data) / "cofactors.yaml"
    keys = list(COFACTORS) + list(MULTI)
    lines = [ln for ln in p.read_text(encoding="utf-8").splitlines()
             if not ln.startswith(MARK) and not any(ln.startswith(f"{k}:") for k in keys)]
    lines.append(MARK + "：Q／QH2 連到結構；cyt c（[15*] = 蛋白質）；多份 H⁺（由 tools/gen_energy.py 產生）")
    for k, (label, smi, mol) in COFACTORS.items():
        link = f", molecule: {mol}" if mol else ""
        lines.append(f'{k}: {{label: "{label}", smiles: "{smi}"{link}}}')
    for k, (label, frag, n_) in MULTI.items():
        lines.append(f'{k}: {{label: "{label}", smiles: "{".".join([frag] * n_)}"}}')
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(data):
    data = Path(data)
    add_cofactors(data)
    mols = molecules()
    for m in mols:
        d = {"id": m["id"], "name": m["name"]}
        if m.get("label"):
            d["label"] = m["label"]
        d["aliases"] = m.get("aliases", [])
        d["smiles"] = m["smiles"]
        for k in ("chain", "free", "haworth", "formula", "decorate"):
            if m.get(k):
                d[k] = m[k]
        d["tags"] = with_parent(m["tags"])
        d["notes_zh"] = m["notes"]
        d["status"] = "pending"
        dump(data / "molecules" / f"{m['id']}.yaml", d)
    # acetaldehyde (from the aldol round) gets a chain for ADH / ALDH; the aldol reaction keeps skeletal
    p = data / "molecules" / "acetaldehyde.yaml"
    acal = yaml.safe_load(p.read_text(encoding="utf-8"))
    acal["chain"] = ACETALDEHYDE_CHAIN
    acal["tags"] = list(dict.fromkeys(acal.get("tags", []) + [ROOT_TAG, "ethanol"]))
    acal.setdefault("notes_zh", "乙醇經 ADH（或 CYP2E1）氧化的產物，再由 ALDH2 氧化成 acetate；毒性與宿醉、臉紅有關。")
    dump(p, {k: acal[k] for k in ("id", "name", "aliases", "smiles", "chain", "tags", "notes_zh", "status") if k in acal})
    for r in REACTIONS:
        d = {"id": r["id"], "title": r["title"], "substrates": r["substrates"], "products": r["products"]}
        for k in ("cofactors_in", "cofactors_out", "extra_in", "extra_out", "prosthetic"):
            if r.get(k):
                d[k] = r[k]
        d["enzyme"] = r["enzyme"]
        if r.get("conditions"):
            d["conditions"] = r["conditions"]
        if r.get("layout"):
            d["layout"] = r["layout"]
        d["reversible"] = r.get("reversible", False)
        deco = COLORS.get(r["id"], r.get("decorate"))
        if deco:
            d["decorate"] = deco
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
