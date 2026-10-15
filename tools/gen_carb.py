"""Carbohydrate metabolism (round 8): gluconeogenesis, glycogen, PPP, fructose, galactose,
polyol pathway, Cori cycle.

    .venv/bin/python tools/gen_carb.py data

New SMILES are never typed stereocentre by stereocentre: ring sugars come from tools/sugar_lib.py
(PubChem) and open chains from molecules already validated in data/, and only atoms *after* an
existing O (phosphate, UDP, glycogen R) are appended, so every @/@@ keeps its meaning.
Re-running overwrites the files listed here.
"""
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent))
from reaction_colors import COLORS  # noqa: E402
from sugar_lib import A_GAL, A_GLC, B_FRU, B_GLC  # noqa: E402

P = "P(=O)([O-])[O-]"
URIDINE_TAIL = "P(=O)([O-])OP(=O)([O-])OC[C@H]5O[C@@H](n6ccc(=O)[nH]c6=O)[C@H](O)[C@@H]5O"   # 1-O … UDP
ANOMERIC_A = "[C@H](O1)O)"        # α-D-glucose / galactose (sugar_lib strings)
ANOMERIC_B = "[C@@H](O1)O)"


def last_oh(smi, add):
    """Append to the final -OH (C6 of the sugar_lib pyranoses, C1 of the furanose)."""
    assert smi.endswith(")O")
    return smi + add


def anomeric(smi, token, add):
    assert smi.count(token) == 1, (smi, token)
    return smi.replace(token, token[:-1] + add + ")")


def load(data, mid):
    return yaml.safe_load((Path(data) / "molecules" / f"{mid}.yaml").read_text(encoding="utf-8"))


GLC_SUBS = {2: {"up": "H", "down": "OH"}, 3: {"up": "OH", "down": "H"}, 4: {"up": "H", "down": "HO"}}
GAL_SUBS = {2: {"up": "H", "down": "OH"}, 3: {"up": "OH", "down": "H"}, 4: {"up": "HO", "down": "H"}}


def pyranose(c1, subs, c5="CH_2OH"):
    return {"ring": "pyranose", "numbers": True,
            "subs": {1: c1, **subs, 5: {"up": c5, "down": "H"}}}


def fru(c2_up, c5_up):
    return {"ring": "furanose", "numbers": True, "start": 2,
            "subs": {2: {"up": c2_up, "down": "CH_2OH"}, 3: {"up": "HO", "down": "H"},
                     4: {"up": "H", "down": "OH"}, 5: {"up": c5_up, "down": "H"}}}


def molecules(data):
    G6P = load(data, "glucose_6p")["smiles"]
    FRU = load(data, "d_fructose")["smiles"]
    RU = load(data, "d_ribulose")["smiles"]
    UTP = load(data, "utp")
    assert FRU.startswith("O=C(CO)") and RU.endswith("CO") and G6P.startswith("O=C[C@H]")
    phos_row = {"c": "CH_2O{P}"}
    M = [
        # ---- gluconeogenesis substrates
        dict(id="lactate", name="L-Lactate", aliases=["lactate", "lactic acid", "乳酸"],
             smiles="C[C@H](O)C(=O)[O-]",
             chain=[{"c": "C", "fork": {"l": "O", "lb": 2, "r": "O^-"}}, {"c": "C", "l": "HO", "r": "H"}, {"c": "CH_3"}],
             tags=["glycolysis", "gluconeogenesis"],
             notes="L-乳酸（Fischer：HO 在左）。無氧糖解的終產物（LDH 再生 NAD⁺，讓 GAPDH 能繼續）；紅血球沒有粒線體，只能做到 lactate。經 Cori cycle 回肝臟做糖質新生。"),
        dict(id="glycerol", name="Glycerol", aliases=["甘油", "glycerin"], smiles="OCC(O)CO",
             chain=[{"c": "CH_2OH"}, {"c": "C", "l": "HO", "r": "H"}, {"c": "CH_2OH"}],
             tags=["gluconeogenesis", "lipid"],
             notes="TAG 水解後的骨架。只有肝（與腎）有 glycerol kinase → 脂肪組織不能直接再利用 glycerol，要送回肝臟做糖質新生。"),
        dict(id="glycerol_3p", name="sn-Glycerol 3-phosphate", label="Glycerol\n3-phosphate",
             aliases=["G3P (glycerol)", "glycerol phosphate", "α-glycerophosphate", "甘油-3-磷酸"],
             smiles="OC[C@@H](O)COP(=O)([O-])[O-]",
             chain=[{"c": "CH_2OH"}, {"c": "C", "l": "HO", "r": "H"}, phos_row],
             tags=["gluconeogenesis", "lipid"],
             notes="TAG 與磷脂合成的骨架（脂肪組織由 DHAP 經 glycerol-3-P DH 還原而來）。也是 glycerol-phosphate shuttle 的載體（細胞質 NADH → 粒線體 FAD → 約 1.5 ATP）。"),
        # ---- PPP oxidative phase
        dict(id="phosphogluconolactone", name="6-Phosphoglucono-δ-lactone", label="6-Phospho-\nglucono-δ-lactone",
             aliases=["6-phosphogluconolactone", "6-磷酸葡萄糖酸內酯"],
             smiles=last_oh(B_GLC.replace(ANOMERIC_B, "C(O1)=O)"), P),
             haworth={"ring": "pyranose", "numbers": True,
                      "subs": {1: {"out": "O", "out_b": 2}, **GLC_SUBS, 5: {"up": "CH_2O{P}", "down": "H"}}},
             tags=["PPP"], notes="G6PD 把 β-G6P 的 C1（半縮醛）氧化成內酯（C1=O，洋紅）；lactonase 水解開環。"),
        dict(id="phosphogluconate", name="6-Phosphogluconate", label="6-Phospho-\ngluconate",
             aliases=["6PG", "6-磷酸葡萄糖酸"], smiles="[O-]C(=O)" + G6P[3:],
             chain=[{"c": "C", "fork": {"l": "O", "lb": 2, "r": "O^-"}}, {"c": "C", "l": "H", "r": "OH"},
                    {"c": "C", "l": "HO", "r": "H"}, {"c": "C", "l": "H", "r": "OH"}, {"c": "C", "l": "H", "r": "OH"},
                    phos_row],
             tags=["PPP"], notes="醛糖酸（C1 為 COO⁻）。6-phosphogluconate DH 把 C3 氧化並脫去 C1 → Ru5P + CO₂ + 第二個 NADPH。"),
        dict(id="ribulose_5p", name="Ribulose 5-phosphate", label="Ribulose\n5-phosphate", aliases=["Ru5P", "核酮糖-5-磷酸"],
             smiles=RU + P,
             chain=[{"c": "CH_2OH"}, {"c": "C", "r": "O", "rb": 2}, {"c": "C", "l": "H", "r": "OH"},
                    {"c": "C", "l": "H", "r": "OH"}, phos_row],
             tags=["PPP"], notes="氧化期的終產物。isomerase → R5P（核苷酸）；epimerase → X5P（非氧化期，transketolase 的供體）。"),
        # ---- ring sugar phosphates (Haworth)
        dict(id="a_glucose_6p", name="α-D-Glucose 6-phosphate", label="Glucose\n6-phosphate", aliases=["α-G6P", "glucose 6-phosphate (ring)"],
             smiles=last_oh(A_GLC, P), haworth=pyranose({"up": "H", "down": "OH"}, GLC_SUBS, "CH_2O{P}"),
             tags=["glycogen"], notes="環狀 G6P（α），phosphoglucomutase 與 G1P 互換。開鏈畫法見 Glucose 6-phosphate。"),
        dict(id="b_glucose_6p", name="β-D-Glucose 6-phosphate", label="Glucose\n6-phosphate", aliases=["β-G6P"],
             smiles=last_oh(B_GLC, P), haworth=pyranose({"up": "OH", "down": "H"}, GLC_SUBS, "CH_2O{P}"),
             tags=["PPP"], notes="G6PD 的受質是 β 異構物（C1-OH 朝上）。"),
        dict(id="a_glucose_1p", name="α-D-Glucose 1-phosphate", label="Glucose\n1-phosphate", aliases=["G1P", "glucose 1-phosphate", "葡萄糖-1-磷酸"],
             smiles=anomeric(A_GLC, ANOMERIC_A, P), haworth=pyranose({"up": "H", "down": "O{P}"}, GLC_SUBS),
             tags=["glycogen", "galactose_metabolism"],
             notes="glycogen phosphorylase 的產物（磷酸解，不耗 ATP）；也是肝醣合成的起點（→ UDP-glucose）。PGM ⇌ G6P。"),
        dict(id="udp_glucose", name="UDP-glucose", label="UDP-glucose", aliases=["UDP-Glc", "UDPG", "尿苷二磷酸葡萄糖"],
             smiles=anomeric(A_GLC, ANOMERIC_A, URIDINE_TAIL), haworth=pyranose({"up": "H", "down": "O–UDP"}, GLC_SUBS),
             tags=["glycogen", "galactose_metabolism"],
             notes="活化的葡萄糖供體（肝醣、乳糖、glucuronate、醣蛋白）。UDP 接在 α-C1 上；PPi 被 pyrophosphatase 水解，使合成不可逆。"),
        dict(id="udp", name="UDP", aliases=["uridine diphosphate", "尿苷二磷酸"],
             smiles=UTP["smiles"].replace("OP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-]", "OP(=O)([O-])OP(=O)([O-])[O-]"),
             formula="UDP", haworth={**UTP["haworth"], "subs": {**UTP["haworth"]["subs"], 4: {"up": "{P}–{P}–OCH_2", "down": "H"}}},
             tags=["glycogen", "nucleotide"], notes="glycogen synthase 放出的 UDP，由 nucleoside diphosphate kinase 用 ATP 回補成 UTP。"),
        dict(id="a_galactopyranose", name="α-D-Galactopyranose", label="α-D-Galactose", aliases=["α-D-galactose", "α-galactose"],
             smiles=A_GAL, haworth=pyranose({"up": "H", "down": "OH"}, GAL_SUBS),
             tags=["carbohydrate", "monosaccharide", "aldose", "galactose_metabolism"],
             notes="與 glucose 只差 C4（C4 epimer：Haworth 中 C4-OH 朝上）。galactokinase 作用於 α 異構物。"),
        dict(id="a_galactose_1p", name="α-D-Galactose 1-phosphate", label="Galactose\n1-phosphate", aliases=["Gal1P", "galactose 1-phosphate", "半乳糖-1-磷酸"],
             smiles=anomeric(A_GAL, ANOMERIC_A, P), haworth=pyranose({"up": "H", "down": "O{P}"}, GAL_SUBS),
             tags=["galactose_metabolism"],
             notes="GALT 缺乏（典型 galactosemia）時累積 → 肝損傷、黃疸、智能障礙、E. coli 敗血症；同時 galactitol 累積 → 白內障。"),
        dict(id="udp_galactose", name="UDP-galactose", label="UDP-galactose", aliases=["UDP-Gal", "尿苷二磷酸半乳糖"],
             smiles=anomeric(A_GAL, ANOMERIC_A, URIDINE_TAIL), haworth=pyranose({"up": "H", "down": "O–UDP"}, GAL_SUBS),
             tags=["galactose_metabolism"],
             notes="乳腺合成 lactose（lactose synthase）、醣脂與醣蛋白的 galactose 供體；UDP-Gal 4-epimerase 與 UDP-Glc 互換（所以不吃半乳糖也能合成）。"),
        # ---- fructose
        dict(id="fructose_1p", name="Fructose 1-phosphate", label="Fructose\n1-phosphate", aliases=["F1P", "果糖-1-磷酸"],
             smiles="O=C(COP(=O)([O-])[O-])" + FRU[len("O=C(CO)"):],
             chain=[phos_row, {"c": "C", "r": "O", "rb": 2}, {"c": "C", "l": "HO", "r": "H"},
                    {"c": "C", "l": "H", "r": "OH"}, {"c": "C", "l": "H", "r": "OH"}, {"c": "CH_2OH"}],
             tags=["fructose_metabolism"],
             notes="肝臟 fructokinase 的產物。aldolase B 缺乏（遺傳性果糖不耐症，HFI）→ F1P 累積、Pi 被困住 → ATP ↓、肝醣分解與糖質新生受抑制 → 低血糖、嘔吐、肝損傷。"),
        dict(id="b_fructose_6p", name="β-D-Fructose 6-phosphate", label="Fructose\n6-phosphate", aliases=["β-F6P", "fructose 6-phosphate (ring)"],
             smiles=last_oh(B_FRU, P), haworth=fru("OH", "{P}OCH_2"),
             tags=["glycolysis", "gluconeogenesis"], notes="環狀 F6P（β-furanose），PFK-2 在 C2 加磷酸 → F2,6BP。開鏈畫法見 Fructose 6-phosphate。"),
        dict(id="fructose_26bp", name="Fructose 2,6-bisphosphate", label="Fructose\n2,6-bisphosphate", aliases=["F2,6BP", "F-2,6-BP", "果糖-2,6-二磷酸"],
             smiles=last_oh(anomeric(B_FRU, "(CO)O)", P), P), haworth=fru("O{P}", "{P}OCH_2"),
             tags=["glycolysis", "gluconeogenesis"],
             notes="不是代謝中間物，而是訊號分子：活化 PFK-1、抑制 FBPase-1。胰島素 → F2,6BP ↑ → 糖解；glucagon（PKA 磷酸化 PFK-2/FBPase-2）→ F2,6BP ↓ → 糖質新生（肝）。"),
        # ---- glycogen chain ends (R = rest of the glycogen particle)
        dict(id="glycogen_n", name="Glycogen (n residues)", label="Glycogen\n(n residues)", aliases=["glycogen chain"],
             smiles=anomeric(A_GLC, ANOMERIC_A, "[9*]"),
             haworth=pyranose({"up": "H", "down": "OR"}, GLC_SUBS),
             tags=["glycogen"], notes="肝醣鏈的非還原端（C4-OH 自由）；R = 肝醣其餘部分（[9*]，平衡檢查視為不變的 R 基）。"),
        dict(id="glycogen_n1", name="Glycogen (n+1 residues)", label="Glycogen\n(n+1 residues)", aliases=["glycogen chain (n+1)"],
             smiles=None,         # filled below (glucosyl α1→4 on glycogen_n)
             haworth={"units": [
                 {"ring": "pyranose", "subs": {1: {"up": "H", "down": "link"}, **GLC_SUBS, 5: {"up": "CH_2OH", "down": "H"}}},
                 {"ring": "pyranose", "subs": {1: {"up": "H", "down": "OR"}, 2: {"up": "H", "down": "OH"},
                                               3: {"up": "OH", "down": "H"}, 4: {"up": "H", "down": "link"},
                                               5: {"up": "CH_2OH", "down": "H"}}, "at": [4.9, 0]}],
                 "links": [{"from": [0, 1, "down"], "to": [1, 4, "down"]}]},
             tags=["glycogen"], notes="多一個以 α1→4 接上的葡萄糖（左）。"),
    ]
    maltose = load(data, "maltose")["smiles"]
    # maltose = β-glucose (reducing end, written first) with its O4 glycosylated: make that end α and cap it with R
    M[-1]["smiles"] = maltose.replace("[C@@H](O1)O)", "[C@H](O1)O[9*])", 1)
    assert M[-1]["smiles"] != maltose
    return M


def ph(i):
    """decorate: colour the Ⓟ of chain row i magenta."""
    return {"colors": {i: {"c": ["magenta", "{P}"]}}}


def hw(key, frag):
    return {"colors": {key: ["magenta", frag]}}


ATP_ADP = dict(cofactors_in=["ATP"], cofactors_out=["ADP"], extra_out=["H+"])
CM = "carb_metabolism"

REACTIONS = [
    # ---------------- gluconeogenesis
    dict(id="pyruvate_carboxylase", title="Gluconeogenesis ① Pyruvate → Oxaloacetate (pyruvate carboxylase)",
         substrates=["pyruvate", {"id": "bicarbonate", "view": "formula"}], products=["oxaloacetate"],
         cofactors_in=["ATP"], cofactors_out=["ADP", "Pi"], extra_out=["H+"], prosthetic=["biotin"],
         enzyme="pyruvate carboxylase", decorate={"oxaloacetate": {"colors": {3: "magenta"}}},
         card={"location": "粒線體基質", "type": "羧化（消耗 1 ATP）", "cofactors": "biotin（B7）、Mg²⁺",
               "regulation": "acetyl-CoA 活化（β 氧化旺盛 → 把 pyruvate 導向糖質新生而非 PDH）",
               "mechanism": "ATP + HCO₃⁻ 先形成 carboxybiotin，再把 CO₂ 轉到 pyruvate 的 CH₃",
               "deficiency": "pyruvate carboxylase 缺乏 → 乳酸中毒、低血糖、發展遲緩",
               "exam_points": ["糖質新生第一個繞道：PC（粒線體）+ PEPCK", "也是 TCA 的 anaplerotic 反應（補 OAA）"]},
         tags=["gluconeogenesis", "TCA"]),
    dict(id="pepck", title="Gluconeogenesis ② Oxaloacetate → PEP (PEP carboxykinase)",
         substrates=["oxaloacetate"], products=["pep", {"id": "co2", "view": "formula"}],
         cofactors_in=["GTP"], cofactors_out=["GDP"], enzyme="PEP carboxykinase\n(PEPCK)",
         decorate={"oxaloacetate": {"colors": {3: "magenta"}}},
         card={"location": "細胞質與粒線體（人類兩種都有）", "type": "脫羧＋磷酸化（消耗 1 GTP）",
               "regulation": "轉錄調控：glucagon、cortisol → ↑；insulin → ↓",
               "exam_points": ["剛加上去的 CO₂ 又被拿掉：羧化只是為了活化", "OAA 不能穿過粒線體內膜：先變成 malate（或 Asp）再出來"]},
         tags=["gluconeogenesis"]),
    dict(id="fbpase1", title="Gluconeogenesis ③ Fructose 1,6-bisphosphate → Fructose 6-phosphate (FBPase-1)",
         substrates=["fructose_16bp"], products=["fructose_6p"], cofactors_in=["H2O"], cofactors_out=["Pi"],
         enzyme="fructose 1,6-\nbisphosphatase", decorate={"fructose_16bp": ph(0)},
         card={"location": "細胞質", "type": "水解（不產生 ATP）",
               "regulation": "F2,6BP 與 AMP 抑制（與 PFK-1 相反）；限速步驟",
               "exam_points": ["與 PFK-1 互為「無效循環」：F2,6BP 同時開一個、關一個"]},
         tags=["gluconeogenesis"]),
    dict(id="g6pase", title="Gluconeogenesis ④ Glucose 6-phosphate → Glucose (glucose 6-phosphatase)",
         substrates=["glucose_6p"], products=["glucose"], cofactors_in=["H2O"], cofactors_out=["Pi"],
         enzyme="glucose 6-phosphatase", decorate={"glucose_6p": ph(5)},
         card={"location": "內質網腔（G6P 經 G6P translocase 進入）；只有肝、腎（腸）有", "type": "水解",
               "deficiency": ["von Gierke（GSD Ia：G6Pase；Ib：translocase）→ 嚴重空腹低血糖、乳酸 ↑、尿酸 ↑、TG ↑、肝腫大"],
               "exam_points": ["肌肉沒有 G6Pase → 肌肝醣只能自己用，不能補血糖"]},
         tags=["gluconeogenesis", "glycogen"]),
    dict(id="ldh", title="Pyruvate ⇌ Lactate (lactate dehydrogenase)",
         substrates=["pyruvate"], products=["lactate"], cofactors_in=["NADH+H+"], cofactors_out=["NAD+"],
         reversible=True, enzyme="lactate\ndehydrogenase",
         decorate={"pyruvate": {"colors": {1: {"r": "magenta"}}}, "lactate": {"colors": {1: {"l": "magenta", "r": "magenta"}}}},
         card={"location": "細胞質", "type": "氧化還原（可逆）",
               "exam_points": ["無氧時再生 NAD⁺ 給 GAPDH", "肝臟反向：lactate → pyruvate → 糖質新生（Cori cycle）",
                               "血中 LDH ↑：組織損傷（溶血、心肌梗塞、腫瘤）"]},
         tags=["glycolysis", "gluconeogenesis"]),
    dict(id="glycerol_kinase", title="Glycerol → Glycerol 3-phosphate (glycerol kinase)",
         substrates=["glycerol"], products=["glycerol_3p"], **ATP_ADP, enzyme="glycerol kinase",
         decorate={"glycerol_3p": ph(2)},
         card={"location": "肝（脂肪組織沒有）", "type": "磷酸化（1 ATP）"}, tags=["gluconeogenesis", "lipid"]),
    dict(id="glycerol_3p_dehydrogenase", title="Glycerol 3-phosphate ⇌ DHAP (glycerol-3-phosphate dehydrogenase)",
         substrates=["glycerol_3p"], products=["dhap"], cofactors_in=["NAD+"], cofactors_out=["NADH+H+"], reversible=True,
         enzyme="glycerol 3-phosphate\ndehydrogenase",
         card={"location": "細胞質（另有粒線體內膜的 FAD 型，構成 glycerol-phosphate shuttle）", "type": "氧化還原",
               "exam_points": ["glycerol 由此進入糖質新生（DHAP）", "反方向：脂肪組織由 DHAP 做出 TAG 的骨架"]},
         tags=["gluconeogenesis", "lipid"]),
    # ---------------- F2,6BP
    dict(id="pfk2", title="Fructose 6-phosphate → Fructose 2,6-bisphosphate (PFK-2)",
         substrates=["b_fructose_6p"], products=["fructose_26bp"], **ATP_ADP, enzyme="PFK-2\n(bifunctional enzyme)",
         decorate={"fructose_26bp": hw("u1", "{P}")},
         card={"location": "細胞質", "type": "磷酸化（調控用，非代謝主幹）",
               "regulation": "同一條多肽有 PFK-2 與 FBPase-2 兩個活性。肝：glucagon → PKA 磷酸化 → PFK-2 關、FBPase-2 開 → F2,6BP ↓；insulin 去磷酸化 → 相反",
               "exam_points": ["F2,6BP：PFK-1 最強的活化劑、FBPase-1 抑制劑"]},
         tags=["glycolysis", "gluconeogenesis"]),
    dict(id="fbpase2", title="Fructose 2,6-bisphosphate → Fructose 6-phosphate (FBPase-2)",
         substrates=["fructose_26bp"], products=["b_fructose_6p"], cofactors_in=["H2O"], cofactors_out=["Pi"],
         enzyme="FBPase-2\n(bifunctional enzyme)", decorate={"fructose_26bp": hw("u1", "{P}")},
         card={"type": "水解", "regulation": "肝中被 PKA 磷酸化後活化（空腹、glucagon）"}, tags=["glycolysis", "gluconeogenesis"]),
    # ---------------- glycogen
    dict(id="phosphoglucomutase", title="Glucose 1-phosphate ⇌ Glucose 6-phosphate (phosphoglucomutase)",
         substrates=["a_glucose_1p"], products=["a_glucose_6p"], reversible=True, enzyme="phosphoglucomutase",
         decorate={"a_glucose_1p": hw("d1", "{P}"), "a_glucose_6p": hw("u5", "{P}")},
         card={"location": "細胞質", "type": "異構化（磷酸從 C1 移到 C6）",
               "mechanism": "酵素的 Ser-Ⓟ 先把磷酸給 C6 → glucose 1,6-bisphosphate → 再從 C1 收回",
               "exam_points": ["肝醣分解的 G1P 由此進入糖解（或在肝變成血糖）"]},
         tags=["glycogen"]),
    dict(id="udp_glucose_pyrophosphorylase", title="Glucose 1-phosphate + UTP → UDP-glucose (UDP-glucose pyrophosphorylase)",
         substrates=["a_glucose_1p", {"id": "utp", "view": "formula"}], products=["udp_glucose"], cofactors_out=["PPi"],
         extra_in=["H+"], enzyme="UDP-glucose\npyrophosphorylase",
         card={"location": "細胞質", "type": "核苷酸轉移（活化葡萄糖）",
               "exam_points": ["PPi → 2 Pi（inorganic pyrophosphatase）把反應拉向合成", "肝醣每加一個葡萄糖耗 1 UTP（≈ 1 ATP）"]},
         tags=["glycogen"]),
    dict(id="glycogen_synthase", title="UDP-glucose + Glycogen(n) → Glycogen(n+1) (glycogen synthase)",
         substrates=["udp_glucose", "glycogen_n"], products=["glycogen_n1", {"id": "udp", "view": "formula"}],
         extra_out=["H+"], enzyme="glycogen synthase",
         card={"location": "細胞質（肝、肌肉）", "type": "醣基轉移：形成 α1→4 醣苷鍵",
               "regulation": "insulin → PP1 去磷酸化 → 活化；glucagon／epinephrine → PKA、GSK3 磷酸化 → 失活；G6P 別構活化",
               "mechanism": "需要 primer：glycogenin 自我醣化起頭；分支由 branching enzyme（α1→6，每 8–12 個殘基）完成",
               "deficiency": ["GSD 0：glycogen synthase 缺乏（空腹低血糖）", "Andersen（GSD IV）：branching enzyme 缺乏 → 長而少分支的肝醣、肝硬化"],
               "exam_points": ["肝醣合成的限速酵素"]},
         tags=["glycogen"]),
    dict(id="glycogen_phosphorylase", title="Glycogen(n+1) + Pi → Glycogen(n) + Glucose 1-phosphate (glycogen phosphorylase)",
         substrates=["glycogen_n1"], products=["glycogen_n", "a_glucose_1p"], cofactors_in=["Pi"],
         prosthetic=["PLP"], enzyme="glycogen phosphorylase", decorate={"a_glucose_1p": hw("d1", "{P}")},
         card={"location": "細胞質（肝、肌肉）", "type": "磷酸解（phosphorolysis）：斷 α1→4，不耗 ATP", "cofactors": "PLP（B6，用它的磷酸當酸鹼觸媒）",
               "regulation": "glucagon（肝）／epinephrine（肝、肌）→ PKA → phosphorylase kinase → 磷酸化活化（a 型）；肌肉另有 AMP、Ca²⁺ 活化；肝中 glucose 抑制",
               "mechanism": "停在分支點前 4 個殘基；debranching enzyme（transferase 搬 3 個 + α1→6 glucosidase 放出 1 個游離 glucose）",
               "deficiency": ["McArdle（GSD V）：肌肉 phosphorylase → 運動不耐、肌肉痙攣、肌紅素尿，乳酸不上升", "Hers（GSD VI）：肝 phosphorylase",
                              "Cori（GSD III）：debranching enzyme → 外層短分支（limit dextrin）"],
               "exam_points": ["肝醣分解的限速酵素", "磷酸解得到 G1P（已磷酸化），比水解多省 1 ATP"]},
         tags=["glycogen"]),
    # ---------------- PPP
    dict(id="g6pd", title="PPP ① Glucose 6-phosphate → 6-Phosphoglucono-δ-lactone (G6PD)",
         substrates=["b_glucose_6p"], products=["phosphogluconolactone"], cofactors_in=["NADP+"], cofactors_out=["NADPH+H+"],
         enzyme="glucose 6-phosphate\ndehydrogenase (G6PD)",
         decorate={"b_glucose_6p": hw("u1", "OH"), "phosphogluconolactone": hw("o1", "O")},
         card={"location": "細胞質", "type": "氧化（第一個 NADPH）",
               "regulation": "限速步驟；NADP⁺ 活化、NADPH 抑制；insulin 誘導",
               "deficiency": ["G6PD 缺乏（X 性聯）：紅血球 NADPH ↓ → glutathione 無法還原 → 氧化壓力下溶血（蠶豆、primaquine、sulfa、感染）；Heinz bodies、bite cells",
                              "對瘧疾有保護作用"],
               "exam_points": ["PPP 的限速酵素", "紅血球 NADPH 的唯一來源"]},
         tags=["PPP"]),
    dict(id="lactonase", title="PPP ② 6-Phosphoglucono-δ-lactone → 6-Phosphogluconate (lactonase)",
         substrates=["phosphogluconolactone"], products=["phosphogluconate"], cofactors_in=["H2O"], extra_out=["H+"],
         enzyme="lactonase", card={"type": "水解（開環）"}, tags=["PPP"]),
    dict(id="phosphogluconate_dehydrogenase", title="PPP ③ 6-Phosphogluconate → Ribulose 5-phosphate (6-phosphogluconate DH)",
         substrates=["phosphogluconate"], products=["ribulose_5p", {"id": "co2", "view": "formula"}],
         cofactors_in=["NADP+"], cofactors_out=["NADPH"], enzyme="6-phosphogluconate\ndehydrogenase",
         decorate={"phosphogluconate": {"colors": {0: "magenta"}}},
         card={"type": "氧化脫羧（第二個 NADPH，C1 → CO₂）", "exam_points": ["氧化期：1 G6P → 2 NADPH + CO₂ + Ru5P，不可逆"]},
         tags=["PPP"]),
    dict(id="ribose_5p_isomerase", title="PPP ④ Ribulose 5-phosphate ⇌ Ribose 5-phosphate (isomerase)",
         substrates=["ribulose_5p"], products=["ribose_5p"], reversible=True, enzyme="ribose 5-phosphate\nisomerase",
         card={"type": "異構化（酮糖 ⇌ 醛糖，經 enediol）", "exam_points": ["R5P → PRPP → 核苷酸"]}, tags=["PPP"]),
    dict(id="ribulose_5p_epimerase", title="PPP ⑤ Ribulose 5-phosphate ⇌ Xylulose 5-phosphate (epimerase)",
         substrates=["ribulose_5p"], products=["xylulose_5p"], reversible=True, enzyme="ribulose 5-phosphate\n3-epimerase",
         decorate={"ribulose_5p": {"colors": {2: {"l": "magenta", "r": "magenta"}}}, "xylulose_5p": {"colors": {2: {"l": "magenta", "r": "magenta"}}}},
         card={"type": "C3 差向異構化", "exam_points": ["X5P 是 transketolase 的 2C 供體；X5P 也活化 PP2A → 肝中 ChREBP、F2,6BP ↑"]},
         tags=["PPP"]),
    dict(id="transaldolase", title="Transaldolase S7P + G3P ⇌ E4P + F6P",
         substrates=["sedoheptulose_7p", "glyceraldehyde_3p"], products=["erythrose_4p", "fructose_6p"], reversible=True,
         enzyme="transaldolase",
         decorate={"sedoheptulose_7p": {"colors": {0: "magenta", 1: "magenta", 2: "magenta"}},
                   "fructose_6p": {"colors": {0: "magenta", 1: "magenta", 2: "magenta"}}},
         card={"location": "細胞質（PPP 非氧化期）", "type": "轉移 3 碳（dihydroxyacetone 單位）",
               "mechanism": "酵素 Lys 與酮糖 C2 形成 Schiff base（不需 TPP）",
               "exam_points": ["transketolase 轉 2C（TPP）；transaldolase 轉 3C"]},
         tags=["PPP"]),
    # ---------------- fructose
    dict(id="fructokinase", title="Fructose → Fructose 1-phosphate (fructokinase)",
         substrates=["d_fructose"], products=["fructose_1p"], **ATP_ADP, enzyme="fructokinase",
         decorate={"fructose_1p": ph(0)},
         card={"location": "肝（以及腎、小腸）", "type": "磷酸化",
               "deficiency": "essential fructosuria：fructokinase 缺乏 → 尿中有果糖，無症狀（良性）",
               "exam_points": ["繞過 PFK-1 → 果糖進入糖解比葡萄糖快、不受調控（大量果糖 → 脂肪生成、尿酸 ↑）"]},
         tags=["fructose_metabolism"]),
    dict(id="aldolase_b", title="Fructose 1-phosphate → DHAP + Glyceraldehyde (aldolase B)",
         substrates=["fructose_1p"], products=["dhap", "d_glyceraldehyde"], reversible=True, enzyme="aldolase B",
         card={"location": "肝", "type": "醛醇裂解",
               "deficiency": "遺傳性果糖不耐症（HFI）：F1P 累積、Pi 耗竭 → 吃果糖／蔗糖後低血糖、嘔吐、黃疸、肝衰竭；治療：避免果糖、蔗糖、sorbitol",
               "exam_points": ["aldolase A（肌肉）主要切 F1,6BP；aldolase B（肝）兩者都能切"]},
         tags=["fructose_metabolism"]),
    dict(id="triose_kinase", title="Glyceraldehyde → Glyceraldehyde 3-phosphate (triose kinase)",
         substrates=["d_glyceraldehyde"], products=["glyceraldehyde_3p"], **ATP_ADP, enzyme="triose kinase",
         card={"location": "肝", "type": "磷酸化"}, tags=["fructose_metabolism"]),
    # ---------------- galactose
    dict(id="galactokinase", title="Galactose → Galactose 1-phosphate (galactokinase)",
         substrates=["a_galactopyranose"], products=["a_galactose_1p"], **ATP_ADP, enzyme="galactokinase",
         decorate={"a_galactose_1p": hw("d1", "{P}")},
         card={"type": "磷酸化",
               "deficiency": "galactokinase 缺乏：galactose → galactitol 累積 → 白內障（嬰兒），其餘症狀輕",
               "exam_points": ["比 GALT 缺乏輕微：沒有 Gal1P 累積"]},
         tags=["galactose_metabolism"]),
    dict(id="galt", title="Galactose 1-phosphate + UDP-glucose → UDP-galactose + Glucose 1-phosphate (GALT)",
         substrates=["a_galactose_1p", "udp_glucose"], products=["udp_galactose", "a_glucose_1p"], enzyme="galactose 1-phosphate\nuridylyltransferase (GALT)",
         card={"type": "核苷酸交換（UMP 從 glucose 換到 galactose）",
               "deficiency": "典型 galactosemia（GALT 缺乏，體染色體隱性）：餵奶後嘔吐、黃疸、肝腫大、白內障、智能障礙、E. coli 敗血症；治療：去除乳糖與半乳糖",
               "exam_points": ["新生兒篩檢項目"]},
         tags=["galactose_metabolism"]),
    dict(id="udp_galactose_epimerase", title="UDP-galactose ⇌ UDP-glucose (UDP-galactose 4-epimerase)",
         substrates=["udp_galactose"], products=["udp_glucose"], reversible=True, prosthetic=["NAD^+"],
         enzyme="UDP-galactose\n4-epimerase", decorate={"udp_galactose": hw("u4", "HO"), "udp_glucose": hw("d4", "HO")},
         card={"type": "C4 差向異構化", "mechanism": "緊密結合的 NAD⁺ 把 C4-OH 氧化成酮，再從另一面還原回來",
               "exam_points": ["可逆 → 不吃半乳糖也能做出 UDP-Gal（乳糖、醣脂）"]},
         tags=["galactose_metabolism"]),
    # ---------------- polyol
    dict(id="aldose_reductase", title="Glucose → Sorbitol (aldose reductase)",
         substrates=["glucose"], products=["sorbitol"], cofactors_in=["NADPH+H+"], cofactors_out=["NADP+"],
         enzyme="aldose reductase", decorate={"glucose": {"colors": {0: "magenta"}}, "sorbitol": {"colors": {0: "magenta"}}},
         card={"location": "水晶體、視網膜、神經（Schwann 細胞）、腎", "type": "還原（醛 → 醇）",
               "deficiency": "高血糖時 sorbitol 累積（缺 sorbitol DH 的組織）→ 滲透壓損傷：糖尿病白內障、神經病變、視網膜病變；也耗用 NADPH",
               "exam_points": ["polyol pathway 第一步"]},
         tags=["polyol"]),
    dict(id="aldose_reductase_galactose", title="Galactose → Galactitol (aldose reductase)",
         substrates=["d_galactose"], products=["galactitol"], cofactors_in=["NADPH+H+"], cofactors_out=["NADP+"],
         enzyme="aldose reductase", decorate={"d_galactose": {"colors": {0: "magenta"}}, "galactitol": {"colors": {0: "magenta"}}},
         card={"type": "還原", "exam_points": ["galactitol 不能再被代謝 → galactosemia 的白內障"]},
         tags=["polyol", "galactose_metabolism"]),
    dict(id="sorbitol_dehydrogenase", title="Sorbitol → Fructose (sorbitol dehydrogenase)",
         substrates=["sorbitol"], products=["d_fructose"], cofactors_in=["NAD+"], cofactors_out=["NADH+H+"],
         enzyme="sorbitol\ndehydrogenase", decorate={"sorbitol": {"colors": {1: "magenta"}}, "d_fructose": {"colors": {1: "magenta"}}},
         card={"location": "肝、卵巢、精囊（精液的果糖）；水晶體、視網膜、神經很少", "type": "氧化（C2 醇 → 酮）"},
         tags=["polyol"]),
]


def n(id_, ref, x, y, **kw):
    return {"id": id_, "ref": ref, "pos": [x, y], **kw}


def t(id_, text, x, y, kind="concept"):
    return {"id": id_, "text": text, "kind": kind, "pos": [x, y]}


NAME = {"view": "name"}
PATHWAYS = [
    dict(id="gluconeogenesis_overview", title="Gluconeogenesis 總覽（名稱版，三個繞道）", layout="free", cell=0.62,
         nodes=[n("glc", "glucose", 0, 0, **NAME), n("g6p", "glucose_6p", 0, 1, **NAME), n("f6p", "fructose_6p", 0, 2, **NAME),
                n("fbp", "fructose_16bp", 0, 3, **NAME), n("dhap", "dhap", -2.1, 4, **NAME), n("g3p", "glyceraldehyde_3p", 0, 4, **NAME),
                n("bpg", "bpg_13", 0, 5, **NAME), n("pg3", "pg_3", 0, 6, **NAME), n("pg2", "pg_2", 0, 7, **NAME),
                n("pep", "pep", 0, 8, **NAME), n("oaa", "oxaloacetate", 2.4, 8.6, **NAME), n("pyr", "pyruvate", 0, 9.2, **NAME),
                n("lac", "lactate", -2.1, 9.2, **NAME), n("ala", "alanine", 2.1, 9.9, **NAME),
                n("gol", "glycerol", -2.1, 2.4, **NAME), n("gol3p", "glycerol_3p", -2.1, 3.2, **NAME)],
         edges=[{"from": "pyr", "to": "oaa", "step": "pyruvate_carboxylase", "label": "pyruvate carboxylase\n(mitochondria)", "side_text": "ATP, CO_2", "label_side": "right"},
                {"from": "oaa", "to": "pep", "step": "pepck", "label": "PEPCK", "side_text": "GTP → CO_2"},
                {"from": "pep", "to": "pg2", "step": "enolase", "reversible": True},
                {"from": "pg2", "to": "pg3", "step": "phosphoglycerate_mutase", "reversible": True},
                {"from": "pg3", "to": "bpg", "step": "phosphoglycerate_kinase", "reversible": True, "side_text": "ATP"},
                {"from": "bpg", "to": "g3p", "step": "gapdh", "reversible": True, "side_text": "NADH"},
                {"from": "dhap", "to": "g3p", "step": "triose_phosphate_isomerase", "reversible": True},
                {"from": "g3p", "to": "fbp", "step": "aldolase", "reversible": True},
                {"from": "fbp", "to": "f6p", "step": "fbpase1", "label": "FBPase-1", "side_text": "H_2O → P_i"},
                {"from": "f6p", "to": "g6p", "step": "phosphoglucose_isomerase", "reversible": True},
                {"from": "g6p", "to": "glc", "step": "g6pase", "label": "glucose 6-phosphatase\n(ER; liver, kidney)", "side_text": "H_2O → P_i"},
                {"from": "lac", "to": "pyr", "step": "ldh", "label": "LDH", "reversible": True},
                {"from": "ala", "to": "pyr", "step": "alt", "label": "ALT"},
                {"from": "gol", "to": "gol3p", "step": "glycerol_kinase", "label": "glycerol kinase"},
                {"from": "gol3p", "to": "dhap", "step": "glycerol_3p_dehydrogenase"}],
         tags=["gluconeogenesis"],
         notes="糖解的三個不可逆步驟（hexokinase、PFK-1、pyruvate kinase）以四個酵素繞過：PC + PEPCK、FBPase-1、G6Pase。2 pyruvate → glucose 共耗 4 ATP + 2 GTP + 2 NADH。原料：lactate（Cori）、Ala（glucose–alanine cycle）、glycerol、生糖胺基酸、propionyl-CoA（奇數碳脂肪酸）；偶數碳脂肪酸（acetyl-CoA）不能淨生成葡萄糖。"),
    dict(id="glycogen_metabolism", title="Glycogen synthesis & degradation（肝醣合成與分解）", layout="free", cell=0.62,
         nodes=[n("glc", "a_glucopyranose", 0, 0), n("g6p", "a_glucose_6p", 0, 1.15), n("g1p", "a_glucose_1p", 1.25, 1.15),
                n("udpg", "udp_glucose", 2.5, 0.3), n("gly", "glycogen", 2.5, 2.0),
                t("gk", "Glycolysis\n(muscle)", 0, 2.1, "plain")],
         edges=[{"from": "glc", "to": "g6p", "step": "hexokinase", "label": "hexokinase /\nglucokinase", "side_text": "ATP", "bend": 0.25},
                {"from": "g6p", "to": "glc", "step": "g6pase", "label": "G6Pase\n(liver)", "side_text": "P_i", "bend": 0.25},
                {"from": "g6p", "to": "g1p", "step": "phosphoglucomutase", "reversible": True, "label": "PGM"},
                {"from": "g1p", "to": "udpg", "step": "udp_glucose_pyrophosphorylase", "label": "UDP-glucose\npyrophosphorylase", "side_text": "UTP → PP_i"},
                {"from": "udpg", "to": "gly", "step": "glycogen_synthase", "label": "glycogen synthase\n+ branching enzyme", "side_text": "→ UDP"},
                {"from": "gly", "to": "g1p", "step": "glycogen_phosphorylase", "label": "glycogen phosphorylase\n+ debranching enzyme", "side_text": "P_i", "label_side": "right"},
                {"from": "g6p", "to": "gk"}],
         tags=["glycogen"],
         notes="合成：G6P → G1P → UDP-Glc → 肝醣（synthase 限速，insulin 活化）。分解：phosphorylase（PLP，限速，glucagon／epinephrine 活化）→ G1P → G6P；肝有 G6Pase 可釋出葡萄糖，肌肉沒有。GSD：0 synthase、Ia von Gierke（G6Pase）、II Pompe（溶小體 acid α-glucosidase）、III Cori（debranching）、IV Andersen（branching）、V McArdle（肌 phosphorylase）、VI Hers（肝 phosphorylase）。"),
    dict(id="ppp_oxidative", title="Pentose phosphate pathway: oxidative phase（PPP 氧化期）", layout="free", cell=0.8,
         nodes=[n("g6p", "b_glucose_6p", 0, 0), n("lac", "phosphogluconolactone", 1.15, 0), n("pg", "phosphogluconate", 2.3, 0),
                n("ru5p", "ribulose_5p", 2.3, 1.1), n("r5p", "ribose_5p", 1.5, 2.1), n("x5p", "xylulose_5p", 3.1, 2.1),
                t("non", "Non-oxidative phase (reversible):\ntransketolase (TPP, 2C), transaldolase (3C)\n→ F6P + G3P", 3.6, 1.1, "plain")],
         edges=[{"from": "g6p", "to": "lac", "step": "g6pd", "label": "G6PD", "side_text": "NADP^+ → NADPH"},
                {"from": "lac", "to": "pg", "step": "lactonase", "label": "lactonase", "side_text": "H_2O"},
                {"from": "pg", "to": "ru5p", "step": "phosphogluconate_dehydrogenase", "label": "6PG DH", "side_text": "NADP^+ → NADPH, CO_2"},
                {"from": "ru5p", "to": "r5p", "step": "ribose_5p_isomerase", "label": "isomerase", "reversible": True},
                {"from": "ru5p", "to": "x5p", "step": "ribulose_5p_epimerase", "label": "epimerase", "reversible": True, "label_side": "right"}],
         tags=["PPP"],
         notes="氧化期不可逆：G6P → 2 NADPH + CO₂ + Ru5P；G6PD 限速（NADP⁺ 活化）。非氧化期可逆：依需求在 R5P（核苷酸）與 F6P／G3P（糖解）之間轉換。NADPH 用於：脂肪酸、膽固醇合成，glutathione 還原，CYP450，呼吸爆。"),
    dict(id="fructose_metabolism", title="Fructose metabolism（果糖代謝，肝）", layout="free", cell=0.9,
         nodes=[n("fru", "d_fructose", 0.5, 0), n("f1p", "fructose_1p", 0.5, 1.4), n("dhap", "dhap", 0, 2.6),
                n("gald", "d_glyceraldehyde", 1.1, 2.6), n("g3p", "glyceraldehyde_3p", 1.1, 3.6),
                t("gly", "Glycolysis /\ngluconeogenesis", 0.0, 3.6, "plain")],
         edges=[{"from": "fru", "to": "f1p", "step": "fructokinase", "label": "fructokinase", "side_text": "ATP → ADP"},
                {"from": "f1p", "to": "dhap", "step": "aldolase_b", "label": "aldolase B"},
                {"from": "f1p", "to": "gald", "step": "aldolase_b", "label_side": "right"},
                {"from": "gald", "to": "g3p", "step": "triose_kinase", "label": "triose kinase", "side_text": "ATP → ADP", "label_side": "right"},
                {"from": "dhap", "to": "g3p", "step": "triose_phosphate_isomerase", "reversible": True},
                {"from": "dhap", "to": "gly"}],
         tags=["fructose_metabolism"],
         notes="果糖在肝經 F1P 進入糖解，繞過 PFK-1（不受調控）。essential fructosuria（fructokinase，良性）vs 遺傳性果糖不耐症（aldolase B，F1P 累積、Pi 耗竭 → 低血糖、肝損傷）。肌肉、脂肪則由 hexokinase 直接做成 F6P。"),
    dict(id="galactose_metabolism", title="Galactose metabolism（半乳糖代謝，Leloir pathway）", layout="free", cell=0.9,
         nodes=[n("gal", "a_galactopyranose", 0, 0), n("gal1p", "a_galactose_1p", 0, 1.15), n("g1p", "a_glucose_1p", 0, 2.3),
                n("udpglc", "udp_glucose", 1.35, 1.15), n("udpgal", "udp_galactose", 1.35, 2.3),
                n("galol", "galactitol", -1.2, 0.2)],
         edges=[{"from": "gal", "to": "gal1p", "step": "galactokinase", "label": "galactokinase", "side_text": "ATP → ADP"},
                {"from": "gal1p", "to": "g1p", "step": "galt", "label": "GALT"},
                {"from": "udpglc", "to": "udpgal", "step": "galt", "label": "GALT", "label_side": "right"},
                {"from": "udpgal", "to": "udpglc", "step": "udp_galactose_epimerase", "label": "epimerase", "bend": 0.5},
                {"from": "gal", "to": "galol", "step": "aldose_reductase_galactose", "label": "aldose\nreductase"}],
         tags=["galactose_metabolism"],
         notes="Galactose → Gal1P →（GALT 與 UDP-Glc 交換）→ G1P → G6P。GALT 缺乏（典型 galactosemia）最嚴重；galactokinase 缺乏只有白內障。兩者都有 galactitol 累積（aldose reductase）→ 白內障。"),
    dict(id="polyol_pathway", title="Polyol pathway（多元醇路徑）", layout="linear",
         nodes=[n("glc", "glucose", 0, 0), n("sor", "sorbitol", 1, 0), n("fru", "d_fructose", 2, 0)],
         edges=[{"from": "glc", "to": "sor", "step": "aldose_reductase", "label": "aldose reductase", "side_text": "NADPH → NADP^+"},
                {"from": "sor", "to": "fru", "step": "sorbitol_dehydrogenase", "label": "sorbitol DH", "side_text": "NAD^+ → NADH"}],
         tags=["polyol"],
         notes="水晶體、視網膜、Schwann 細胞缺 sorbitol DH → 高血糖時 sorbitol 堆積 → 滲透壓損傷（白內障、神經病變）。精囊用這條路徑把血糖變成精液中的果糖。"),
    dict(id="cori_cycle", title="Cori cycle（乳酸循環）", layout="free", cell_size=6.0, text=False,
         nodes=[t("hm", "Muscle / RBC", 0, -0.7, "heading"), t("hl", "Liver", 2.6, -0.7, "heading"),
                n("glc_m", "glucose", 0, 0, **NAME), n("pyr_m", "pyruvate", 0, 1, **NAME), n("lac_m", "lactate", 0, 2, **NAME),
                n("lac_l", "lactate", 2.6, 2, **NAME), n("pyr_l", "pyruvate", 2.6, 1, **NAME), n("glc_l", "glucose", 2.6, 0, **NAME)],
         edges=[{"from": "glc_m", "to": "pyr_m", "label": "glycolysis", "side_text": "+2 ATP"},
                {"from": "pyr_m", "to": "lac_m", "step": "ldh", "label": "LDH", "side_text": "NADH → NAD^+"},
                {"from": "lac_m", "to": "lac_l", "label": "blood"},
                {"from": "lac_l", "to": "pyr_l", "step": "ldh", "label": "LDH", "label_side": "right"},
                {"from": "pyr_l", "to": "glc_l", "label": "gluconeogenesis", "side_text": "−6 ATP", "label_side": "right"},
                {"from": "glc_l", "to": "glc_m", "label": "blood"}],
         tags=["gluconeogenesis", "glycolysis"],
         notes="肌肉（劇烈運動）與紅血球產生的乳酸送到肝臟做回葡萄糖：肝臟淨付出 4 ATP（6 − 2），代價由肝的脂肪酸氧化支付。類似的 glucose–alanine cycle 則同時把肌肉的胺基送去肝臟做尿素。"),
]

TAGS_CHILDREN = ["glycolysis", "gluconeogenesis", "glycogen", "PPP", "fructose_metabolism", "galactose_metabolism", "polyol"]


def dump(path, d):
    path.write_text(yaml.safe_dump(d, allow_unicode=True, sort_keys=False, width=200), encoding="utf-8")


def with_parent(tags):
    out = list(tags)
    if any(t in TAGS_CHILDREN for t in out) and CM not in out:
        out.insert(0, CM)
    return out


def main(data):
    data = Path(data)
    for m in molecules(data):
        d = {"id": m["id"], "name": m["name"]}
        if m.get("label"):
            d["label"] = m["label"]
        d["aliases"] = m.get("aliases", [])
        d["smiles"] = m["smiles"]
        for k in ("chain", "haworth", "formula"):
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
        d["reversible"] = r.get("reversible", False)
        deco = COLORS.get(r["id"], r.get("decorate"))       # tools/reaction_colors.py wins
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
    print(f"{len(molecules(data))} molecules, {len(REACTIONS)} reactions, {len(PATHWAYS)} pathways")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data")
