"""Lipid metabolism (round 9): fatty-acid activation and β-oxidation, odd-chain fatty acids,
fatty-acid synthesis, ketone bodies, cholesterol synthesis, TAG metabolism, lipoproteins.

    .venv/bin/python tools/gen_lipmet.py data

Generic acyl chains are written R–CH₂–CH₂–C(=O)–S-CoA with R = [8*] (the rest of the fatty acyl
chain; a chain row may be written "R"). Stereocentres next to R are taken from the real
3-hydroxybutyryl-CoA with its CH₃ swapped for [8*], so the spatial arrangement is unchanged.
Also appends multi-copy cofactors (8 CoA-SH, 14 NADPH + H⁺ …) to data/cofactors.yaml.
"""
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent))
from reaction_colors import COLORS  # noqa: E402

CM = "lipid_metabolism"
KIDS = ["beta_oxidation", "fa_synthesis", "ketone_bodies", "cholesterol_synthesis", "tag_metabolism", "lipoprotein"]

COO = {"c": "C", "fork": {"l": "O", "lb": 2, "r": "O^-"}}
TE = {"c": "C", "l": "O", "lb": 2, "r": "S-CoA"}           # thioester carbon, last row
OPP = "OP(=O)([O-])OP(=O)([O-])[O-]"

MOLS = [
    # ---- activation / transport
    dict(id="palmitoyl_coa", name="Palmitoyl-CoA", aliases=["palmitoyl coenzyme A", "C16:0-CoA", "棕櫚醯輔酶A"],
         smiles="[1*]SC(=O)CCCCCCCCCCCCCCC", lipid={"type": "fatty_acid", "chains": [{"n": 16}], "end": "coa"},
         tags=["beta_oxidation", "fa_synthesis"],
         notes="脂肪酸活化後的形式（高能硫酯）。細胞質中由 acyl-CoA synthetase（耗 2 個高能磷酸鍵：ATP → AMP + PPi）生成；要進粒線體必須先換成 acylcarnitine。"),
    dict(id="carnitine", name="L-Carnitine", aliases=["carnitine", "肉鹼"],
         smiles="[O-]C(=O)C[C@@H](O)C[N+](C)(C)C",
         chain=[COO, {"c": "CH_2"}, {"c": "C", "l": "HO", "r": "H"}, {"c": "CH_2"}, {"c": "N^+(CH_3)_3", "anchor": 0}],
         tags=["beta_oxidation"],
         notes="由 Lys + Met 合成（需 vitamin C）。長鏈 acyl-CoA 以 acylcarnitine 的形式穿過粒線體內膜。缺乏（原發或 CPT-II 缺乏）→ 空腹低酮性低血糖、肌肉無力、肌紅素尿。"),
    dict(id="palmitoylcarnitine", name="Palmitoylcarnitine", aliases=["acylcarnitine", "棕櫚醯肉鹼"],
         smiles="C[N+](C)(C)C[C@@H](CC(=O)[O-])OC(=O)CCCCCCCCCCCCCCC",
         lipid={"type": "fatty_acid", "chains": [{"n": 16}], "end": "carnitine"}, tags=["beta_oxidation"],
         notes="CPT-I（外膜）生成 → carnitine–acylcarnitine translocase（內膜）→ CPT-II（基質側）換回 acyl-CoA。CPT-I 被 malonyl-CoA 抑制：合成與分解不會同時進行。"),
    # ---- β-oxidation (generic chain)
    dict(id="acyl_coa", name="Fatty acyl-CoA", label="Fatty acyl-CoA", aliases=["acyl-CoA", "脂醯輔酶A"],
         smiles="[8*]CCC(=O)S[1*]", chain=[{"c": "R"}, {"c": "CH_2"}, {"c": "CH_2"}, TE], tags=["beta_oxidation", "fa_synthesis"],
         notes="β 氧化的受質（R = 其餘的脂肪酸碳鏈）。由上而下：β 碳（C3）、α 碳（C2）、硫酯羰基（C1）。"),
    dict(id="enoyl_coa", name="trans-Δ2-Enoyl-CoA", label="trans-Δ^2-Enoyl-CoA", aliases=["enoyl-CoA", "烯醯輔酶A"],
         smiles="[8*]/C=C/C(=O)S[1*]", chain=[{"c": "R"}, {"c": "CH", "down": 2}, {"c": "CH"}, TE], tags=["beta_oxidation"],
         notes="acyl-CoA DH（FAD）在 α、β 碳之間形成 trans 雙鍵（此直式畫法看不出 trans，名稱標明）。天然不飽和脂肪酸是 cis，需要 isomerase 才能進入此步。"),
    dict(id="hydroxyacyl_coa", name="L-3-Hydroxyacyl-CoA", label="L-3-Hydroxyacyl-CoA", aliases=["β-hydroxyacyl-CoA", "3-羥醯輔酶A"],
         smiles="[8*][C@H](O)CC(=O)S[1*]", chain=[{"c": "R"}, {"c": "C", "l": "H", "r": "OH"}, {"c": "CH_2"}, TE], tags=["beta_oxidation"],
         notes="enoyl-CoA hydratase 加水，OH 在 β 碳（L 型；合成時的中間物則是 D 型，接在 ACP 上）。"),
    dict(id="ketoacyl_coa", name="3-Ketoacyl-CoA", label="β-Ketoacyl-CoA", aliases=["β-ketoacyl-CoA", "3-oxoacyl-CoA"],
         smiles="[8*]C(=O)CC(=O)S[1*]", chain=[{"c": "R"}, {"c": "C", "r": "O", "rb": 2}, {"c": "CH_2"}, TE], tags=["beta_oxidation"],
         notes="3-hydroxyacyl-CoA DH（NAD⁺）把 β 碳氧化成酮；thiolase 在 α–β 之間以 CoA-SH 切開。"),
    dict(id="acyl_coa_n2", name="Fatty acyl-CoA (n−2)", label="Fatty acyl-CoA\n(2 C shorter)", aliases=["shortened acyl-CoA"],
         smiles="[8*]C(=O)S[1*]", chain=[{"c": "R"}, TE], tags=["beta_oxidation"],
         notes="原本的 β 碳變成新的硫酯羰基；少兩個碳，回到第一步再循環。"),
    # ---- odd chain
    dict(id="propionyl_coa", name="Propionyl-CoA", aliases=["丙醯輔酶A"], smiles="CCC(=O)S[1*]",
         chain=[{"c": "CH_3"}, {"c": "CH_2"}, TE], tags=["beta_oxidation"],
         notes="奇數碳脂肪酸 β 氧化的最後 3 碳；也來自 Val、Ile、Met、Thr 與膽固醇側鏈。→ succinyl-CoA（可淨生成葡萄糖）。"),
    dict(id="s_methylmalonyl_coa", name="(S)-Methylmalonyl-CoA", label="(S)-Methylmalonyl-CoA", aliases=["D-methylmalonyl-CoA"],
         smiles="[O-]C(=O)[C@H](C)C(=O)S[1*]", chain=[COO, {"c": "C", "l": "H", "r": "CH_3"}, TE], tags=["beta_oxidation"],
         notes="propionyl-CoA carboxylase（biotin）的產物。"),
    dict(id="r_methylmalonyl_coa", name="(R)-Methylmalonyl-CoA", label="(R)-Methylmalonyl-CoA", aliases=["L-methylmalonyl-CoA"],
         smiles="[O-]C(=O)[C@@H](C)C(=O)S[1*]", chain=[COO, {"c": "C", "l": "H_3C", "r": "H"}, TE], tags=["beta_oxidation"],
         notes="epimerase 轉成 (R)，才是 methylmalonyl-CoA mutase（adenosylcobalamin）的受質。B12 缺乏 → methylmalonic acid ↑（與 folate 缺乏的鑑別）。"),
    # ---- synthesis
    dict(id="malonyl_coa", name="Malonyl-CoA", aliases=["丙二醯輔酶A"], smiles="[O-]C(=O)CC(=O)S[1*]",
         chain=[COO, {"c": "CH_2"}, TE], tags=["fa_synthesis"],
         notes="ACC（biotin）的產物、脂肪酸合成每輪的 2C 供體（加入時放出 CO₂ 推動縮合）。同時抑制 CPT-I → 合成時不會同時 β 氧化。"),
    # ---- ketone bodies
    dict(id="acetoacetyl_coa", name="Acetoacetyl-CoA", aliases=["乙醯乙醯輔酶A"], smiles="CC(=O)CC(=O)S[1*]",
         chain=[{"c": "CH_3"}, {"c": "C", "r": "O", "rb": 2}, {"c": "CH_2"}, TE], tags=["ketone_bodies", "cholesterol_synthesis"],
         notes="2 acetyl-CoA 經 thiolase 縮合；酮體合成與膽固醇合成的共同起點（前者在粒線體，後者在細胞質）。"),
    dict(id="hmg_coa", name="HMG-CoA", label="HMG-CoA", aliases=["3-hydroxy-3-methylglutaryl-CoA", "β-hydroxy-β-methylglutaryl-CoA"],
         smiles="[O-]C(=O)C[C@](C)(O)CC(=O)S[1*]",
         chain=[COO, {"c": "CH_2"}, {"c": "C", "l": "H_3C", "r": "OH"}, {"c": "CH_2"}, TE], tags=["ketone_bodies", "cholesterol_synthesis"],
         notes="(S)-HMG-CoA。粒線體：HMG-CoA lyase → acetoacetate（酮體）；細胞質（ER）：HMG-CoA reductase → mevalonate（膽固醇）。"),
    dict(id="acetoacetate", name="Acetoacetate", aliases=["乙醯乙酸"], smiles="[O-]C(=O)CC(C)=O",
         chain=[COO, {"c": "CH_2"}, {"c": "C", "r": "O", "rb": 2}, {"c": "CH_3"}], tags=["ketone_bodies"],
         notes="酮體。肝臟製造但不能使用（缺 SCOT）；腦（長期飢餓）、心、肌肉、腎皮質用它。自發脫羧 → acetone（呼氣水果味）。"),
    dict(id="beta_hydroxybutyrate", name="D-β-Hydroxybutyrate", label="D-β-Hydroxybutyrate", aliases=["3-hydroxybutyrate", "BHB", "β-羥丁酸"],
         smiles="[O-]C(=O)C[C@@H](C)O", chain=[COO, {"c": "CH_2"}, {"c": "C", "l": "H", "r": "OH"}, {"c": "CH_3"}], tags=["ketone_bodies"],
         notes="血中最多的酮體（DKA 時 β-HB : AcAc 可達 10 : 1）。尿液試紙（nitroprusside）只測 acetoacetate，所以會低估。"),
    dict(id="acetone", name="Acetone", aliases=["丙酮"], smiles="CC(C)=O",
         chain=[{"c": "CH_3"}, {"c": "C", "r": "O", "rb": 2}, {"c": "CH_3"}], tags=["ketone_bodies"],
         notes="acetoacetate 自發脫羧的產物，不能再利用，由肺排出（DKA 的水果味口氣）。"),
    # ---- cholesterol
    dict(id="mevalonate", name="Mevalonate", aliases=["mevalonic acid", "甲羥戊酸"], smiles="[O-]C(=O)C[C@](C)(O)CCO",
         chain=[COO, {"c": "CH_2"}, {"c": "C", "l": "H_3C", "r": "OH"}, {"c": "CH_2"}, {"c": "CH_2OH"}], tags=["cholesterol_synthesis"],
         notes="(R)-mevalonate。HMG-CoA reductase（限速、statin 標的）的產物；硫酯被 2 NADPH 還原成一級醇。"),
    dict(id="ipp", name="Isopentenyl pyrophosphate", label="Isopentenyl\npyrophosphate (IPP)", aliases=["IPP", "isopentenyl diphosphate"],
         smiles="C=C(C)CC" + OPP, chain=[{"c": "CH_2", "down": 2}, {"c": "C", "r": "CH_3"}, {"c": "CH_2"}, {"c": "CH_2O{P}–{P}"}],
         tags=["cholesterol_synthesis"], notes="活化的異戊二烯（C5）單位。mevalonate + 3 ATP → IPP + CO₂。"),
    dict(id="dmapp", name="Dimethylallyl pyrophosphate", label="Dimethylallyl\npyrophosphate (DMAPP)", aliases=["DMAPP", "dimethylallyl diphosphate"],
         smiles="CC(C)=CC" + OPP, chain=[{"c": "C", "l": "H_3C", "r": "CH_3", "down": 2}, {"c": "CH"}, {"c": "CH_2O{P}–{P}"}],
         tags=["cholesterol_synthesis"], notes="IPP 的異構物（C5）。DMAPP + IPP → GPP（C10）→ + IPP → FPP（C15）；2 FPP → squalene（C30）。"),
    dict(id="squalene", name="Squalene", aliases=["鯊烯", "角鯊烯"],
         smiles="CC(C)=CCC/C(C)=C/CC/C(C)=C/CC/C=C(\\C)CC/C=C(\\C)CCC=C(C)C", ring={"methyls": False}, tags=["cholesterol_synthesis"],
         notes="30 碳直鏈三萜（6 個異戊二烯單位，中間尾對尾）。squalene epoxidase → 環化成 lanosterol（第一個固醇）→（約 19 步）→ cholesterol。"),
]

# multi-copy cofactors for overall equations (label, fragment, copies)
MULTI = {"2NADPH+H+": ("2 NADPH + 2 H^+", "NC(=O)C1=CN([6*])C=CC1.[H+]", 2),
         "2NADP+": ("2 NADP^+", "NC(=O)c1ccc[n+]([6*])c1", 2),
         "14NADPH+H+": ("14 NADPH + 14 H^+", "NC(=O)C1=CN([6*])C=CC1.[H+]", 14),
         "14NADP+": ("14 NADP^+", "NC(=O)c1ccc[n+]([6*])c1", 14),
         "7CO2": ("7 CO_2", "O=C=O", 7), "6H2O": ("6 H_2O", "O", 6), "8CoA-SH": ("8 CoA-SH", "S[1*]", 8),
         "3ATP": ("3 ATP", "[4*]OP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-]", 3),
         "3ADP": ("3 ADP", "[4*]OP(=O)([O-])OP(=O)([O-])[O-]", 3)}


def ph(i, frag="{P}"):
    return {"colors": {i: {"c": ["magenta", frag]}}}


REACTIONS = [
    dict(id="acyl_coa_synthetase", title="Palmitate → Palmitoyl-CoA (acyl-CoA synthetase)",
         substrates=["palmitate"], products=["palmitoyl_coa"], cofactors_in=["CoA-SH", "ATP"], cofactors_out=["AMP", "PPi"],
         enzyme="acyl-CoA synthetase\n(thiokinase)",
         card={"location": "粒線體外膜、ER（細胞質側）", "type": "活化：先形成 acyl-adenylate，再接上 CoA",
               "exam_points": ["耗 2 個高能磷酸鍵（ATP → AMP + PPi，PPi 再水解）", "算 ATP 收支時，活化記為 −2 ATP"]},
         tags=["beta_oxidation"]),
    dict(id="cpt1", title="Palmitoyl-CoA + Carnitine → Palmitoylcarnitine (CPT-I)",
         substrates=["palmitoyl_coa", "carnitine"], products=["palmitoylcarnitine"], cofactors_out=["CoA-SH"],
         enzyme="carnitine palmitoyl-\ntransferase I (CPT-I)",
         card={"location": "粒線體外膜", "type": "醯基轉移（硫酯 → 酯）",
               "regulation": "malonyl-CoA 抑制（β 氧化的限速步驟）",
               "deficiency": ["CPT-I 缺乏（肝）：空腹低酮性低血糖", "CPT-II 缺乏（肌肉）：運動後肌肉痛、肌紅素尿"],
               "exam_points": ["中、短鏈脂肪酸不需要 carnitine 就能進粒線體"]},
         tags=["beta_oxidation"]),
    dict(id="acyl_coa_dehydrogenase", title="β-Oxidation ① Acyl-CoA → trans-Δ²-Enoyl-CoA (acyl-CoA dehydrogenase)",
         substrates=["acyl_coa"], products=["enoyl_coa"], cofactors_in=["FAD"], cofactors_out=["FADH2"],
         enzyme="acyl-CoA dehydrogenase", decorate={"acyl_coa": {"colors": {1: "magenta", 2: "magenta"}}, "enoyl_coa": {"colors": {1: "magenta", 2: "magenta"}}},
         card={"location": "粒線體基質", "type": "氧化（FAD）：α、β 碳間形成 trans 雙鍵",
               "cofactors": "FAD → 電子經 ETF → Q（約 1.5 ATP）",
               "deficiency": "MCAD 缺乏（最常見）：空腹低酮性低血糖、dicarboxylic aciduria，嬰兒猝死；避免空腹",
               "exam_points": ["依鏈長有 VLCAD、LCAD、MCAD、SCAD"]},
         tags=["beta_oxidation"]),
    dict(id="enoyl_coa_hydratase", title="β-Oxidation ② trans-Δ²-Enoyl-CoA → L-3-Hydroxyacyl-CoA (enoyl-CoA hydratase)",
         substrates=["enoyl_coa"], products=["hydroxyacyl_coa"], cofactors_in=["H2O"], reversible=True,
         enzyme="enoyl-CoA hydratase", decorate={"hydroxyacyl_coa": {"colors": {1: {"r": "magenta"}}}},
         card={"location": "粒線體基質", "type": "加水（只產生 L 型）"}, tags=["beta_oxidation"]),
    dict(id="hydroxyacyl_coa_dehydrogenase", title="β-Oxidation ③ L-3-Hydroxyacyl-CoA → 3-Ketoacyl-CoA (3-hydroxyacyl-CoA DH)",
         substrates=["hydroxyacyl_coa"], products=["ketoacyl_coa"], cofactors_in=["NAD+"], cofactors_out=["NADH+H+"],
         enzyme="L-3-hydroxyacyl-CoA\ndehydrogenase",
         decorate={"hydroxyacyl_coa": {"colors": {1: {"l": "cyan", "r": "cyan"}}}, "ketoacyl_coa": {"colors": {1: {"r": "magenta"}}}},
         card={"location": "粒線體基質", "type": "氧化（NAD⁺）：β-OH → β-酮",
               "exam_points": ["長鏈的 ②③④ 由 trifunctional protein 執行（LCHAD 缺乏與孕婦 AFLP、HELLP 有關）"]},
         tags=["beta_oxidation"]),
    dict(id="beta_ketothiolase", title="β-Oxidation ④ 3-Ketoacyl-CoA → Acyl-CoA (n−2) + Acetyl-CoA (thiolase)",
         substrates=["ketoacyl_coa"], products=["acyl_coa_n2", "acetyl_coa"], cofactors_in=["CoA-SH"],
         enzyme="β-ketothiolase",
         decorate={"ketoacyl_coa": {"colors": {2: "magenta", 3: "magenta"}}, "acetyl_coa": {"colors": {0: "magenta"}}},
         card={"location": "粒線體基質", "type": "硫解（thiolysis）：CoA-SH 在 α–β 碳之間切開",
               "exam_points": ["palmitate（16C）：7 輪 → 8 acetyl-CoA + 7 FADH₂ + 7 NADH", "淨 ATP：8×10 + 7×1.5 + 7×2.5 − 2 = 106"]},
         tags=["beta_oxidation"]),
    dict(id="propionyl_coa_carboxylase", title="Propionyl-CoA → (S)-Methylmalonyl-CoA (propionyl-CoA carboxylase)",
         substrates=["propionyl_coa", {"id": "bicarbonate", "view": "formula"}], products=["s_methylmalonyl_coa"],
         cofactors_in=["ATP"], cofactors_out=["ADP", "Pi"], extra_out=["H+"], prosthetic=["biotin"],
         enzyme="propionyl-CoA\ncarboxylase", decorate={"s_methylmalonyl_coa": {"colors": {0: "magenta"}}},
         card={"location": "粒線體基質", "type": "羧化（biotin、ATP）",
               "deficiency": "propionic acidemia：propionyl-CoA carboxylase 缺乏 → 代謝性酸中毒、高血氨"},
         tags=["beta_oxidation"]),
    dict(id="methylmalonyl_coa_epimerase", title="(S)-Methylmalonyl-CoA ⇌ (R)-Methylmalonyl-CoA (epimerase)",
         substrates=["s_methylmalonyl_coa"], products=["r_methylmalonyl_coa"], reversible=True, enzyme="methylmalonyl-CoA\nepimerase",
         card={"type": "差向異構化"}, tags=["beta_oxidation"]),
    dict(id="methylmalonyl_coa_mutase", title="(R)-Methylmalonyl-CoA → Succinyl-CoA (methylmalonyl-CoA mutase)",
         substrates=["r_methylmalonyl_coa"], products=["succinyl_coa"], prosthetic=["B_{12}"], enzyme="methylmalonyl-CoA\nmutase",
         decorate={"r_methylmalonyl_coa": {"colors": {2: {"r": "magenta"}}}, "succinyl_coa": {"colors": {3: "magenta"}}},
         card={"location": "粒線體基質", "type": "分子內重排（自由基機制）", "cofactors": "adenosylcobalamin（B12）",
               "deficiency": ["B12 缺乏 → methylmalonic acid ↑（尿、血），神經症狀", "mutase 缺乏 → methylmalonic acidemia"],
               "exam_points": ["succinyl-CoA 進入 TCA → OAA → 葡萄糖：奇數碳脂肪酸的末端可生糖"]},
         tags=["beta_oxidation"]),
    dict(id="citrate_lyase", title="Citrate → Acetyl-CoA + Oxaloacetate (ATP-citrate lyase)",
         substrates=["citrate"], products=["acetyl_coa", "oxaloacetate"], cofactors_in=["CoA-SH", "ATP"], cofactors_out=["ADP", "Pi"],
         enzyme="ATP-citrate lyase",
         card={"location": "細胞質", "type": "裂解（耗 ATP）",
               "exam_points": ["acetyl-CoA 不能穿過粒線體內膜：以 citrate 形式出來（citrate shuttle）", "insulin 誘導"]},
         tags=["fa_synthesis"]),
    dict(id="malic_enzyme", title="Malate → Pyruvate (malic enzyme)",
         substrates=["malate"], products=["pyruvate", {"id": "co2", "view": "formula"}], cofactors_in=["NADP+"], cofactors_out=["NADPH"],
         enzyme="malic enzyme", decorate={"malate": {"colors": {0: "magenta"}}},
         card={"location": "細胞質", "type": "氧化脫羧（產生 NADPH）",
               "exam_points": ["脂肪酸合成的 NADPH 來源：PPP（主要）＋ malic enzyme"]},
         tags=["fa_synthesis"]),
    dict(id="acetyl_coa_carboxylase", title="Acetyl-CoA → Malonyl-CoA (acetyl-CoA carboxylase)",
         substrates=["acetyl_coa", {"id": "bicarbonate", "view": "formula"}], products=["malonyl_coa"],
         cofactors_in=["ATP"], cofactors_out=["ADP", "Pi"], extra_out=["H+"], prosthetic=["biotin"],
         enzyme="acetyl-CoA carboxylase\n(ACC)", decorate={"malonyl_coa": {"colors": {0: "magenta"}}},
         card={"location": "細胞質", "type": "羧化（biotin、ATP）",
               "regulation": "限速。citrate 活化（聚合）；palmitoyl-CoA 抑制；AMPK 與 glucagon／epinephrine（PKA）磷酸化 → 失活；insulin 去磷酸化 → 活化",
               "exam_points": ["脂肪酸合成的限速酵素", "malonyl-CoA 同時抑制 CPT-I"]},
         tags=["fa_synthesis"]),
    dict(id="fatty_acid_synthase", title="Fatty acid synthase (overall): Acetyl-CoA + 7 Malonyl-CoA → Palmitate",
         substrates=["acetyl_coa", {"id": "malonyl_coa", "n": 7}], products=["palmitate"],
         cofactors_in=["14NADPH+H+"], cofactors_out=["14NADP+", "7CO2", "8CoA-SH", "6H2O"], extra_in=["H+"] * 6,
         enzyme="fatty acid synthase (7 cycles)",
         card={"location": "細胞質（肝、脂肪、乳腺）", "type": "多功能酵素（同源二聚體），醯基接在 ACP 的 phosphopantetheine 上",
               "mechanism": "每輪：縮合（KS，放出 CO₂）→ 還原（KR，NADPH）→ 脫水（DH）→ 還原（ER，NADPH）；與 β 氧化方向相反，中間物為 D-3-hydroxyacyl",
               "exam_points": ["每加 2 碳耗 2 NADPH", "產物是 palmitate（16:0）；延長在 ER、去飽和（Δ9 等）也在 ER", "人類不能在 Δ9 之後加雙鍵 → linoleate、α-linolenate 為必需脂肪酸"]},
         tags=["fa_synthesis"]),
    dict(id="acetoacetyl_coa_thiolase", title="2 Acetyl-CoA ⇌ Acetoacetyl-CoA (thiolase)",
         substrates=[{"id": "acetyl_coa", "n": 2}], products=["acetoacetyl_coa"], cofactors_out=["CoA-SH"], reversible=True,
         enzyme="thiolase", card={"location": "粒線體（酮體）／細胞質（膽固醇）", "type": "Claisen 縮合（β 氧化最後一步的逆反應）"},
         tags=["ketone_bodies", "cholesterol_synthesis"]),
    dict(id="hmg_coa_synthase", title="Acetoacetyl-CoA + Acetyl-CoA → HMG-CoA (HMG-CoA synthase)",
         substrates=["acetoacetyl_coa", "acetyl_coa"], products=["hmg_coa"], cofactors_in=["H2O"], cofactors_out=["CoA-SH"], extra_out=["H+"],
         enzyme="HMG-CoA synthase", decorate={"hmg_coa": {"colors": {0: "magenta", 1: "magenta"}}},
         card={"location": "粒線體（肝，酮體）；細胞質同功酶（膽固醇）", "type": "縮合＋硫酯水解",
               "regulation": "粒線體型是酮體合成的限速酵素（glucagon、飢餓誘導）"},
         tags=["ketone_bodies", "cholesterol_synthesis"]),
    dict(id="hmg_coa_lyase", title="HMG-CoA → Acetoacetate + Acetyl-CoA (HMG-CoA lyase)",
         substrates=["hmg_coa"], products=["acetoacetate", "acetyl_coa"], enzyme="HMG-CoA lyase",
         decorate={"hmg_coa": {"colors": {3: "magenta", 4: "magenta"}}, "acetyl_coa": {"colors": {0: "magenta"}}},
         card={"location": "肝粒線體", "type": "裂解", "deficiency": "HMG-CoA lyase 缺乏 → 低酮性低血糖、代謝性酸中毒（也影響 Leu 分解）"},
         tags=["ketone_bodies"]),
    dict(id="hydroxybutyrate_dehydrogenase", title="Acetoacetate ⇌ D-β-Hydroxybutyrate (β-hydroxybutyrate DH)",
         substrates=["acetoacetate"], products=["beta_hydroxybutyrate"], cofactors_in=["NADH+H+"], cofactors_out=["NAD+"], reversible=True,
         enzyme="D-β-hydroxybutyrate\ndehydrogenase", decorate={"acetoacetate": {"colors": {2: {"r": "magenta"}}}, "beta_hydroxybutyrate": {"colors": {2: {"l": "magenta", "r": "magenta"}}}},
         card={"location": "粒線體", "type": "氧化還原（可逆）", "exam_points": ["NADH/NAD⁺ 越高（飢餓、酒精、DKA）→ β-HB 比例越高"]},
         tags=["ketone_bodies"]),
    dict(id="acetoacetate_decarboxylation", title="Acetoacetate → Acetone + CO₂ (spontaneous)",
         substrates=["acetoacetate"], products=["acetone", {"id": "co2", "view": "formula"}], extra_in=["H+"], enzyme="", conditions="spontaneous",
         decorate={"acetoacetate": {"colors": {0: "magenta"}}}, card={"type": "非酵素性脫羧"}, tags=["ketone_bodies"]),
    dict(id="scot", title="Acetoacetate + Succinyl-CoA → Acetoacetyl-CoA + Succinate (SCOT)",
         substrates=["acetoacetate", "succinyl_coa"], products=["acetoacetyl_coa", "succinate"],
         enzyme="succinyl-CoA:3-ketoacid\nCoA transferase (thiophorase)",
         card={"location": "肝外組織粒線體（心、腦、肌肉、腎）", "type": "CoA 轉移",
               "exam_points": ["肝臟沒有 SCOT → 肝製造酮體但不能利用", "接著 thiolase → 2 acetyl-CoA → TCA"]},
         tags=["ketone_bodies"]),
    dict(id="hmg_coa_reductase", title="HMG-CoA → Mevalonate (HMG-CoA reductase)",
         substrates=["hmg_coa"], products=["mevalonate"], cofactors_in=["2NADPH+H+"], cofactors_out=["2NADP+", "CoA-SH"],
         enzyme="HMG-CoA reductase", decorate={"mevalonate": {"colors": {4: "magenta"}}},
         card={"location": "ER 膜（細胞質側）", "type": "還原（2 NADPH）：硫酯 → 一級醇",
               "regulation": "限速。膽固醇 → SREBP 不被活化（轉錄 ↓）並促進酵素降解；insulin 去磷酸化 → 活化；glucagon、AMPK 磷酸化 → 失活",
               "exam_points": ["statin 競爭性抑制 → 肝細胞膽固醇 ↓ → LDL receptor ↑ → 血 LDL ↓", "副作用：肌病變、肝酵素上升"]},
         tags=["cholesterol_synthesis"]),
    dict(id="mevalonate_to_ipp", title="Mevalonate → Isopentenyl pyrophosphate (3 steps)",
         substrates=["mevalonate"], products=["ipp", {"id": "co2", "view": "formula"}], cofactors_in=["3ATP"], cofactors_out=["3ADP", "Pi"], extra_out=["H+"],
         enzyme="kinases ×2 +\ndecarboxylase", card={"location": "細胞質", "type": "兩次磷酸化＋脫羧（共 3 ATP）"},
         tags=["cholesterol_synthesis"]),
    dict(id="ipp_isomerase", title="IPP ⇌ DMAPP (isopentenyl-diphosphate isomerase)",
         substrates=["ipp"], products=["dmapp"], reversible=True, enzyme="IPP isomerase",
         decorate={"ipp": {"colors": {0: "magenta", 1: "magenta"}}, "dmapp": {"colors": {0: "magenta", 1: "magenta"}}},
         card={"type": "異構化（雙鍵移位）"}, tags=["cholesterol_synthesis"]),
]


def n(id_, ref, x, y, **kw):
    return {"id": id_, "ref": ref, "pos": [x, y], **kw}


def t(id_, text, x, y, kind="concept"):
    return {"id": id_, "text": text, "kind": kind, "pos": [x, y]}


NAME = {"view": "name"}
PATHWAYS = [
    dict(id="fatty_acid_activation", title="Fatty acid activation & carnitine shuttle（脂肪酸活化與肉鹼穿梭）", layout="free", cell=0.7,
         nodes=[n("pal", "palmitate", 0, 0), n("pcoa", "palmitoyl_coa", 0, 1), n("pcar", "palmitoylcarnitine", 0, 2),
                t("mat", "Mitochondrial matrix:\nCPT-II → acyl-CoA → β-oxidation", 0, 2.85, "plain"),
                t("mal", "Malonyl-CoA", 1.25, 1.5, "drug")],
         edges=[{"from": "pal", "to": "pcoa", "step": "acyl_coa_synthetase", "label": "acyl-CoA synthetase", "side_text": "CoA, ATP → AMP + PP_i"},
                {"from": "pcoa", "to": "pcar", "step": "cpt1", "label": "CPT-I (outer membrane)", "side_text": "carnitine → CoA"},
                {"from": "pcar", "to": "mat", "label": "translocase"},
                {"from": "mal", "to": "edge:pcoa->pcar", "inhibit": True}],
         tags=["beta_oxidation"], notes="細胞質活化（−2 ATP）→ CPT-I（限速，malonyl-CoA 抑制）→ translocase → CPT-II → 基質中的 acyl-CoA。"),
    dict(id="beta_oxidation", title="β-Oxidation（β 氧化，一輪四步）", layout="free", cell_size=9.5,
         nodes=[n("acyl", "acyl_coa", 0, 0), n("enoyl", "enoyl_coa", 1.5, 0), n("oh", "hydroxyacyl_coa", 1.5, 1.25),
                n("keto", "ketoacyl_coa", 0, 1.25), n("n2", "acyl_coa_n2", -0.7, 2.45), n("ac", "acetyl_coa", 0.7, 2.45)],
         edges=[{"from": "acyl", "to": "enoyl", "step": "acyl_coa_dehydrogenase", "label": "acyl-CoA DH", "side_text": "FAD → FADH_2"},
                {"from": "enoyl", "to": "oh", "step": "enoyl_coa_hydratase", "label": "enoyl-CoA\nhydratase", "side_text": "H_2O", "label_side": "right"},
                {"from": "oh", "to": "keto", "step": "hydroxyacyl_coa_dehydrogenase", "label": "3-hydroxyacyl-CoA DH", "side_text": "NAD^+ → NADH", "label_side": "right"},
                {"from": "keto", "to": "n2", "step": "beta_ketothiolase", "label": "thiolase", "side_text": "CoA-SH"},
                {"from": "keto", "to": "ac", "step": "beta_ketothiolase"},
                {"from": "n2", "to": "acyl", "label": "next cycle", "bend": 0.45, "label_side": "right"}],
         tags=["beta_oxidation"],
         notes="氧化（FAD）→ 加水 → 氧化（NAD⁺）→ 硫解。每輪少 2 碳，產生 1 FADH₂、1 NADH、1 acetyl-CoA。palmitate：7 輪 → 8 acetyl-CoA，淨 106 ATP。不飽和脂肪酸另需 enoyl-CoA isomerase（與 2,4-dienoyl-CoA reductase）；極長鏈先在 peroxisome 縮短（Zellweger、X-ALD）。"),
    dict(id="odd_chain_fatty_acids", title="Odd-chain fatty acids: Propionyl-CoA → Succinyl-CoA", layout="linear",
         nodes=[n("pro", "propionyl_coa", 0, 0), n("smm", "s_methylmalonyl_coa", 1, 0), n("rmm", "r_methylmalonyl_coa", 2, 0), n("suc", "succinyl_coa", 3, 0)],
         edges=[{"from": "pro", "to": "smm", "step": "propionyl_coa_carboxylase", "label": "carboxylase\n(biotin)", "side_text": "ATP, HCO_3^-"},
                {"from": "smm", "to": "rmm", "step": "methylmalonyl_coa_epimerase", "label": "epimerase", "reversible": True},
                {"from": "rmm", "to": "suc", "step": "methylmalonyl_coa_mutase", "label": "mutase\n(B_{12})"}],
         tags=["beta_oxidation"], notes="propionyl-CoA（奇數碳脂肪酸、Val、Ile、Met、Thr）→ succinyl-CoA → TCA。biotin 與 B12 各用一次。B12 缺乏 → MMA ↑。"),
    dict(id="fatty_acid_synthesis", title="Fatty acid synthesis（脂肪酸合成，名稱版）", layout="free", cell_size=6.2, text=False,
         nodes=[t("hm", "Mitochondrion", 0, -0.7, "heading"), t("hc", "Cytosol", 1.6, -0.7, "heading"),
                n("pyr", "pyruvate", 0, 0, **NAME), n("ac_m", "acetyl_coa", 0, 1, **NAME), n("cit_m", "citrate", 0, 2, **NAME),
                n("oaa_m", "oxaloacetate", -1.1, 1.5, **NAME),
                n("cit_c", "citrate", 1.6, 2, **NAME), n("ac_c", "acetyl_coa", 1.6, 3, **NAME), n("mal", "malonyl_coa", 1.6, 4, **NAME),
                n("pal", "palmitate", 1.6, 5, **NAME), n("oaa_c", "oxaloacetate", 3.2, 2.5, **NAME), n("mal_c", "malate", 3.2, 1.5, **NAME),
                n("pyr_c", "pyruvate", 3.2, 0.5, **NAME)],
         edges=[{"from": "pyr", "to": "ac_m", "step": "pdh", "label": "PDH"},
                {"from": "ac_m", "to": "cit_m", "step": "citrate_synthase", "label": "citrate synthase"},
                {"from": "oaa_m", "to": "cit_m"},
                {"from": "cit_m", "to": "cit_c", "label": "citrate\ntransporter"},
                {"from": "cit_c", "to": "ac_c", "step": "citrate_lyase", "label": "ATP-citrate lyase", "label_side": "right"},
                {"from": "cit_c", "to": "oaa_c"},
                {"from": "ac_c", "to": "mal", "step": "acetyl_coa_carboxylase", "label": "ACC (biotin)\nrate-limiting", "side_text": "ATP, HCO_3^-"},
                {"from": "mal", "to": "pal", "step": "fatty_acid_synthase", "label": "fatty acid synthase", "side_text": "14 NADPH"},
                {"from": "oaa_c", "to": "mal_c", "label": "MDH", "side_text": "NADH"},
                {"from": "mal_c", "to": "pyr_c", "step": "malic_enzyme", "label": "malic enzyme", "side_text": "NADP^+ → NADPH"}],
         tags=["fa_synthesis"],
         notes="餐後（insulin）：粒線體的 acetyl-CoA 以 citrate 形式送出 → ATP-citrate lyase → ACC（限速；citrate 活化、palmitoyl-CoA 抑制；AMPK／PKA 磷酸化抑制）→ FAS。NADPH 來自 PPP 與 malic enzyme。"),
    dict(id="ketone_bodies", title="Ketone body synthesis（酮體合成，肝粒線體）", layout="free", cell_size=9.0,
         nodes=[n("ac", "acetyl_coa", 0, 0), n("aac", "acetoacetyl_coa", 0, 1.1), n("hmg", "hmg_coa", 0, 2.4),
                n("aa", "acetoacetate", 1.6, 2.4), n("bhb", "beta_hydroxybutyrate", 3.2, 2.4), n("ace", "acetone", 1.6, 3.7)],
         edges=[{"from": "ac", "to": "aac", "step": "acetoacetyl_coa_thiolase", "label": "thiolase", "side_text": "2 acetyl-CoA → CoA"},
                {"from": "aac", "to": "hmg", "step": "hmg_coa_synthase", "label": "HMG-CoA synthase\n(rate-limiting)", "side_text": "acetyl-CoA, H_2O"},
                {"from": "hmg", "to": "aa", "step": "hmg_coa_lyase", "label": "HMG-CoA lyase", "side_text": "→ acetyl-CoA"},
                {"from": "aa", "to": "bhb", "step": "hydroxybutyrate_dehydrogenase", "label": "β-HB DH", "side_text": "NADH → NAD^+", "reversible": True},
                {"from": "aa", "to": "ace", "step": "acetoacetate_decarboxylation", "label": "spontaneous", "side_text": "CO_2"}],
         tags=["ketone_bodies"],
         notes="飢餓、糖尿病（insulin ↓、glucagon ↑）時，β 氧化產生的 acetyl-CoA 超過 TCA 的處理量（OAA 被拿去做糖質新生）→ 酮體。肝製造、肝外利用（SCOT：succinyl-CoA → succinate；肝缺 SCOT）。腦在長期飢餓時可用酮體提供約 2/3 能量；紅血球不能用（無粒線體）。"),
    dict(id="cholesterol_synthesis", title="Cholesterol synthesis（膽固醇合成，名稱版）", layout="free", cell=0.62,
         nodes=[n("ac", "acetyl_coa", 0, 0, **NAME), n("aac", "acetoacetyl_coa", 0, 1, **NAME), n("hmg", "hmg_coa", 0, 2, **NAME),
                n("mev", "mevalonate", 0, 3, **NAME), n("ipp", "ipp", 0, 4, **NAME), n("dmapp", "dmapp", 2.2, 4, **NAME),
                t("gpp", "Geranyl-PP (C10)", 0, 5, "name"), t("fpp", "Farnesyl-PP (C15)", 0, 6, "name"),
                n("sq", "squalene", 0, 7, **NAME), t("lan", "Lanosterol", 0, 8, "name"), n("cho", "cholesterol", 0, 9, **NAME),
                t("statin", "Statins", -1.7, 2.5, "drug"),
                t("fpp_out", "Dolichol, CoQ,\nprenylated proteins", 2.2, 6, "plain")],
         edges=[{"from": "ac", "to": "aac", "step": "acetoacetyl_coa_thiolase", "label": "thiolase"},
                {"from": "aac", "to": "hmg", "step": "hmg_coa_synthase", "label": "HMG-CoA synthase"},
                {"from": "hmg", "to": "mev", "step": "hmg_coa_reductase", "label": "HMG-CoA reductase\n(rate-limiting)", "side_text": "2 NADPH", "label_side": "right"},
                {"from": "mev", "to": "ipp", "step": "mevalonate_to_ipp", "side_text": "3 ATP, −CO_2"},
                {"from": "ipp", "to": "dmapp", "step": "ipp_isomerase", "reversible": True},
                {"from": "ipp", "to": "gpp", "side_text": "+ DMAPP"},
                {"from": "gpp", "to": "fpp", "side_text": "+ IPP"},
                {"from": "fpp", "to": "sq", "label": "squalene synthase", "side_text": "2 FPP, NADPH"},
                {"from": "sq", "to": "lan", "label": "epoxidase, cyclase", "side_text": "O_2, NADPH"},
                {"from": "lan", "to": "cho", "label": "~19 steps"},
                {"from": "statin", "to": "edge:hmg->mev", "inhibit": True},
                {"from": "fpp", "to": "fpp_out", "head": True}],
         tags=["cholesterol_synthesis"],
         notes="全部在細胞質／ER，由 acetyl-CoA 開始（27 碳全部來自 acetyl-CoA），需要 NADPH 與 ATP。限速：HMG-CoA reductase（膽固醇經 SREBP–SCAP 回饋抑制；insulin 活化；statin 抑制）。Smith–Lemli–Opitz：7-DHC reductase 缺乏（最後一步）。"),
    dict(id="tag_metabolism", title="TAG synthesis & lipolysis（三酸甘油酯合成與分解，名稱版）", layout="free", cell_size=6.2, text=False,
         nodes=[t("hs", "Synthesis (liver, adipose)", 0, -0.7, "heading"), t("hl", "Lipolysis (adipose)", 2.6, -0.7, "heading"),
                n("g3p", "glycerol_3p", 0, 0, **NAME), t("lpa", "Lysophosphatidate", 0, 1, "name"),
                n("pa", "phosphatidate", 0, 2, **NAME), t("dag", "Diacylglycerol", 0, 3, "name"), n("tag", "triacylglycerol", 0, 4, **NAME),
                n("tag2", "triacylglycerol", 2.6, 0, **NAME), t("dag2", "Diacylglycerol", 2.6, 1, "name"), t("mag", "Monoacylglycerol", 2.6, 2, "name"),
                n("gol", "glycerol", 2.6, 3, **NAME), t("fa", "3 Fatty acids\n(→ albumin)", 3.9, 1.5, "plain"),
                n("dhap", "dhap", -1.3, -0.6, **NAME), n("gol0", "glycerol", -1.3, 0.6, **NAME),
                t("pl", "Phospholipids\n(PC, PE, PS, PI…)", -1.4, 2.5, "plain")],
         edges=[{"from": "dhap", "to": "g3p", "step": "glycerol_3p_dehydrogenase", "label": "G3P DH", "label_side": "right"},
                {"from": "gol0", "to": "g3p", "step": "glycerol_kinase", "label": "glycerol kinase\n(liver only)"},
                {"from": "g3p", "to": "lpa", "label": "acyltransferase", "side_text": "acyl-CoA"},
                {"from": "lpa", "to": "pa", "label": "acyltransferase", "side_text": "acyl-CoA"},
                {"from": "pa", "to": "dag", "label": "PA phosphatase", "side_text": "P_i"},
                {"from": "dag", "to": "tag", "label": "DGAT", "side_text": "acyl-CoA"},
                {"from": "pa", "to": "pl"}, {"from": "dag", "to": "pl"},
                {"from": "tag2", "to": "dag2", "label": "ATGL"},
                {"from": "dag2", "to": "mag", "label": "HSL\n(PKA ↑, insulin ↓)"},
                {"from": "mag", "to": "gol", "label": "MGL"},
                {"from": "dag2", "to": "fa"}],
         tags=["tag_metabolism"],
         notes="合成：glycerol-3-P（脂肪組織只能由 DHAP 來；肝也可用 glycerol kinase）＋ 3 acyl-CoA。PA 與 DAG 也是磷脂的前驅物。分解：ATGL → HSL（epinephrine／glucagon → PKA 活化；insulin 抑制）→ MGL；脂肪酸接 albumin 運送，glycerol 回肝做糖質新生。"),
    dict(id="lipoproteins", title="Lipoproteins（脂蛋白代謝，概念圖）", layout="free", cell_size=6.4, text=False,
         nodes=[t("gut", "Intestine", 0, 0, "heading"), t("liv", "Liver", 0, 2.2, "heading"),
                t("cm", "Chylomicron\n(apoB-48, apoC-II, apoE)", 0.6, 0.6, "name"),
                t("cmr", "Chylomicron remnant", 2.6, 0.6, "name"),
                t("vldl", "VLDL\n(apoB-100, apoC-II, apoE)", 0.6, 2.8, "name"), t("idl", "IDL", 2.6, 2.8, "name"),
                t("ldl", "LDL\n(apoB-100)", 4.2, 2.8, "name"), t("rec", "LDL receptor\n(liver, tissues)", 4.2, 4.0, "name"),
                t("hdl", "HDL (apoA-I)\nLCAT, CETP", 2.6, 4.6, "name"), t("tis", "Peripheral cells\n(ABCA1)", 0.6, 4.6, "name"),
                t("lpl", "LPL (capillary wall)\nTAG → FA + glycerol", 1.6, 1.7, "concept")],
         edges=[{"from": "cm", "to": "cmr", "label": "LPL (apoC-II)"},
                {"from": "liv", "to": "vldl"},
                {"from": "cmr", "to": "liv", "label": "apoE receptor", "bend": -0.2},
                {"from": "vldl", "to": "idl", "label": "LPL"},
                {"from": "idl", "to": "ldl", "label": "hepatic lipase"},
                {"from": "ldl", "to": "rec", "label": "apoB-100"},
                {"from": "tis", "to": "hdl", "label": "cholesterol efflux"},
                {"from": "hdl", "to": "liv", "label": "SR-BI (reverse transport)", "bend": -0.45}],
         tags=["lipoprotein"],
         notes="外源：乳糜微粒（腸，apoB-48）→ LPL（需 apoC-II）→ 殘餘物由肝經 apoE 回收。內源：VLDL（肝，apoB-100）→ IDL → LDL → LDL receptor（家族性高膽固醇血症：受體缺陷）。HDL：apoA-I、LCAT 酯化，逆向膽固醇運輸；CETP 把膽固醇酯換到 VLDL／LDL。Ⅰ型（LPL 或 apoC-II 缺乏）→ 乳糜微粒血症、胰臟炎；Ⅲ型 apoE2；Tangier（ABCA1）HDL 極低；abetalipoproteinemia（MTP）。"),
]


def dump(path, d):
    path.write_text(yaml.safe_dump(d, allow_unicode=True, sort_keys=False, width=200), encoding="utf-8")


def with_parent(tags):
    out = list(tags)
    if any(t in KIDS for t in out) and CM not in out:
        out.insert(0, CM)
    return out


def add_cofactors(data):
    p = Path(data) / "cofactors.yaml"
    text = p.read_text(encoding="utf-8")
    lines = [ln for ln in text.splitlines() if not ln.startswith("# --- 脂質代謝") and
             not any(ln.startswith(f"{k}:") for k in MULTI)]
    lines.append("# --- 脂質代謝（第九輪）：總反應式用的多份輔因子（由 tools/gen_lipmet.py 產生）")
    for k, (label, frag, n_) in MULTI.items():
        lines.append(f'{k}: {{label: "{label}", smiles: "{".".join([frag] * n_)}"}}')
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(data):
    data = Path(data)
    add_cofactors(data)
    for m in MOLS:
        d = {"id": m["id"], "name": m["name"]}
        if m.get("label"):
            d["label"] = m["label"]
        d["aliases"] = m.get("aliases", [])
        d["smiles"] = m["smiles"]
        for k in ("lipid", "chain", "ring"):
            if m.get(k):
                d[k] = m[k]
        d["tags"] = with_parent(m["tags"])
        d["notes_zh"] = m["notes"]
        d["status"] = "pending"
        dump(data / "molecules" / f"{m['id']}.yaml", d)
    for r in REACTIONS:
        d = {"id": r["id"], "title": r["title"], "substrates": r["substrates"], "products": r["products"]}
        for k in ("products_layout", "cofactors_in", "cofactors_out", "extra_in", "extra_out", "prosthetic"):
            if r.get(k):
                d[k] = r[k]
        d["enzyme"] = r["enzyme"]
        if r.get("conditions"):
            d["conditions"] = r["conditions"]
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
    print(f"{len(MOLS)} molecules, {len(REACTIONS)} reactions, {len(PATHWAYS)} pathways")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data")
