"""Amino-acid metabolism, completed (round 10): catecholamines and tyrosine degradation,
tryptophan / histidine / glutamate / arginine derivatives, the methionine cycle and one-carbon
metabolism, carbon skeletons into the TCA cycle, BCAA, heme synthesis and degradation.

    .venv/bin/python tools/gen_aamet.py data

Amino-acid derivatives keep the Lehninger layout of their parent amino acid (COO⁻ on top,
H₃N⁺ on the left, R group below) so that a decarboxylation simply removes the top row.
Heme is built with RDKit from protoporphyrin IX (Fe bonded to the four N). Re-running
overwrites the files listed here and appends cofactors (BH4/BH2 …) to data/cofactors.yaml.
"""
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent))
from reaction_colors import COLORS  # noqa: E402
from rdkit import Chem

CM = "amino_acid_metabolism"
KIDS = ["aa_derivatives", "aa_degradation", "one_carbon", "heme"]

COO = {"c": "COO^-"}
ALPHA = {"c": "C", "l": "^+H_3N", "r": "H"}
CATECHOL = "[*]c1ccc(O)c(O)c1"
PHENOL = "[*]c1ccc(O)cc1"
HYDROXYINDOLE = "[*]c1c[nH]c2ccc(O)cc12"
IMIDAZOLE = "[*]c1c[nH]cn1"
PROTOPORPHYRIN = "C=CC1=C(C)C2=CC3=NC(=CC4=C(C)C(CCC(=O)[O-])=C(N4)C=C4N=C(C=C1N2)C(C)=C4C=C)C(CCC(=O)[O-])=C3C"


def heme_smiles():
    """Fe(II) protoporphyrin IX: the two pyrrole N lose H, the two imine N become N⁺ and Fe
    carries −2 with four bonds — overall the neutral heme core plus two propionates."""
    m = Chem.RWMol(Chem.MolFromSmiles(PROTOPORPHYRIN))
    Chem.Kekulize(m, clearAromaticFlags=True)
    fe = m.AddAtom(Chem.Atom(26))
    for a in list(m.GetAtoms()):
        if a.GetSymbol() != "N":
            continue
        if a.GetTotalNumHs() == 1:
            a.SetNumExplicitHs(0)
            a.SetNoImplicit(True)
        else:
            a.SetFormalCharge(1)
        m.AddBond(a.GetIdx(), fe, Chem.BondType.SINGLE)
    m.GetAtomWithIdx(fe).SetFormalCharge(-2)
    mol = m.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


def load(data, mid):
    return yaml.safe_load((Path(data) / "molecules" / f"{mid}.yaml").read_text(encoding="utf-8"))


def molecules(data):
    sam = load(data, "sam")["smiles"]
    thf = load(data, "thf")["smiles"]
    assert "[S+](C)" in sam and ")Nc2c(=O)" in thf
    return [
        # ---- catecholamines
        dict(id="l_dopa", name="L-DOPA", aliases=["levodopa", "3,4-dihydroxyphenylalanine", "左旋多巴"],
             smiles="[NH3+][C@@H](Cc1ccc(O)c(O)c1)C(=O)[O-]", chain=[COO, ALPHA, {"c": "CH_2"}, {"ring": CATECHOL}],
             tags=["aa_derivatives"],
             notes="tyrosine hydroxylase（BH4，兒茶酚胺合成的限速步驟）的產物。Parkinson 治療：L-DOPA 能過血腦屏障，dopamine 不能；併用 carbidopa（周邊 AADC 抑制劑）。也是黑色素的前驅物（tyrosinase）。"),
        dict(id="dopamine", name="Dopamine", aliases=["多巴胺"], smiles="[NH3+]CCc1ccc(O)c(O)c1",
             chain=[{"c": "CH_2", "l": "^+H_3N"}, {"c": "CH_2"}, {"ring": CATECHOL}], tags=["aa_derivatives"],
             notes="DOPA decarboxylase（AADC，PLP）去掉 COO⁻。分解：MAO、COMT → homovanillic acid（HVA）。"),
        dict(id="norepinephrine", name="Norepinephrine", aliases=["noradrenaline", "正腎上腺素"],
             smiles="[NH3+]C[C@H](O)c1ccc(O)c(O)c1",
             chain=[{"c": "CH_2", "l": "^+H_3N"}, {"c": "C", "l": "H", "r": "OH"}, {"ring": CATECHOL}], tags=["aa_derivatives"],
             notes="(R)-norepinephrine。dopamine β-hydroxylase（Cu²⁺、vitamin C、O₂）在 β 碳加 OH（在儲存囊泡內）。分解 → vanillylmandelic acid（VMA；嗜鉻細胞瘤時 ↑）。"),
        dict(id="epinephrine", name="Epinephrine", aliases=["adrenaline", "腎上腺素"],
             smiles="C[NH2+]C[C@H](O)c1ccc(O)c(O)c1",
             chain=[{"c": "CH_2", "l": "CH_3NH_2^+"}, {"c": "C", "l": "H", "r": "OH"}, {"ring": CATECHOL}], tags=["aa_derivatives"],
             notes="PNMT（SAM 提供 CH₃）在腎上腺髓質；cortisol 誘導 PNMT（皮質的血液流經髓質）。"),
        # ---- tyrosine degradation
        dict(id="hydroxyphenylpyruvate", name="p-Hydroxyphenylpyruvate", label="p-Hydroxyphenyl-\npyruvate",
             aliases=["4-hydroxyphenylpyruvate", "HPP"], smiles="[O-]C(=O)C(=O)Cc1ccc(O)cc1",
             chain=[COO, {"c": "C", "l": "O", "lb": 2}, {"c": "CH_2"}, {"ring": PHENOL}], tags=["aa_degradation"],
             notes="tyrosine aminotransferase 的產物（Tyr 的 α-酮酸）。缺乏 TAT → tyrosinemia II。"),
        dict(id="homogentisate", name="Homogentisate", aliases=["homogentisic acid", "尿黑酸"], smiles="[O-]C(=O)Cc1cc(O)ccc1O",
             chain=[COO, {"c": "CH_2"}, {"ring": "[*]c1cc(O)ccc1O"}], tags=["aa_degradation"],
             notes="HPP dioxygenase（Fe²⁺、vitamin C）同時做脫羧、羥化、側鏈移位。homogentisate oxidase 缺乏 → alkaptonuria：尿放置變黑、褐黃病（ochronosis，結締組織發黑）、關節炎。"),
        # ---- tryptophan, histidine, glutamate
        dict(id="five_htp", name="5-Hydroxytryptophan", label="5-Hydroxy-\ntryptophan", aliases=["5-HTP"],
             smiles="[NH3+][C@@H](Cc1c[nH]c2ccc(O)cc12)C(=O)[O-]", chain=[COO, ALPHA, {"c": "CH_2"}, {"ring": HYDROXYINDOLE}],
             tags=["aa_derivatives"], notes="tryptophan hydroxylase（BH4）的產物，serotonin 合成的限速步驟。"),
        dict(id="serotonin", name="Serotonin", aliases=["5-hydroxytryptamine", "5-HT", "血清素"],
             smiles="[NH3+]CCc1c[nH]c2ccc(O)cc12", chain=[{"c": "CH_2", "l": "^+H_3N"}, {"c": "CH_2"}, {"ring": HYDROXYINDOLE}],
             tags=["aa_derivatives"],
             notes="AADC（PLP）脫羧。腸嗜鉻細胞、血小板、腦。分解：MAO → 5-HIAA（carcinoid 時尿 5-HIAA ↑；Trp 被耗用 → pellagra）。松果體 → N-acetylserotonin → melatonin（SAM）。"),
        dict(id="histamine", name="Histamine", aliases=["組織胺"], smiles="[NH3+]CCc1c[nH]cn1",
             chain=[{"c": "CH_2", "l": "^+H_3N"}, {"c": "CH_2"}, {"ring": IMIDAZOLE}], tags=["aa_derivatives"],
             notes="histidine decarboxylase（PLP）。肥大細胞、嗜鹼性球（過敏）、胃 ECL 細胞（H2 → 胃酸）、腦。"),
        dict(id="gaba", name="GABA", aliases=["γ-aminobutyrate", "4-aminobutanoate", "γ-胺基丁酸"], smiles="[NH3+]CCCC(=O)[O-]",
             chain=[{"c": "CH_2", "l": "^+H_3N"}, {"c": "CH_2"}, {"c": "CH_2"}, COO], tags=["aa_derivatives"],
             notes="glutamate decarboxylase（PLP）去掉 Glu 的 α-COO⁻。中樞主要抑制性神經傳導物質；B6 缺乏 → GABA ↓ → 抽搐。"),
        dict(id="nitric_oxide", name="Nitric oxide", aliases=["NO", "一氧化氮"], smiles="[N]=O", formula="NO",
             tags=["aa_derivatives"], notes="NOS（NADPH、FAD、FMN、BH4、heme）由 Arg 的胍基 N 生成，副產物 citrulline。內皮 NO → 平滑肌 guanylyl cyclase → cGMP ↑ → 血管擴張（nitroglycerin、sildenafil 的機制相關）。"),
        # ---- methionine cycle / one carbon
        dict(id="sah", name="S-Adenosylhomocysteine (SAH)", label="S-Adenosyl-\nhomocysteine (SAH)", aliases=["SAH", "AdoHcy"],
             smiles=sam.replace("[S+](C)", "S"), ring={"coordgen": True}, tags=["one_carbon"],
             notes="SAM 把 CH₃ 給出去之後的產物；SAH hydrolase → homocysteine + adenosine。SAH 本身是甲基轉移酶的抑制劑。"),
        dict(id="homocysteine", name="Homocysteine", aliases=["Hcy", "同半胱胺酸"], smiles="SCC[C@H]([NH3+])C(=O)[O-]",
             chain=[COO, ALPHA, {"c": "CH_2"}, {"c": "CH_2"}, {"c": "SH"}], tags=["one_carbon"],
             notes="Met 少一個 CH₃。兩條出路：再甲基化回 Met（methionine synthase：B12 + N⁵-CH₃-THF；或肝的 betaine）；transsulfuration → Cys（CBS，B6）。Hcy ↑：B12、folate、B6 缺乏或 CBS 缺乏 → 血栓、動脈硬化。"),
        dict(id="cystathionine", name="Cystathionine", aliases=["胱硫醚"], smiles="[O-]C(=O)[C@@H]([NH3+])CCSC[C@H]([NH3+])C(=O)[O-]",
             chain=[COO, ALPHA, {"c": "CH_2"}, {"c": "CH_2"}, {"c": "S"}, {"c": "CH_2"}, {"c": "C", "l": "H", "r": "NH_3^+"}, COO],
             tags=["one_carbon"],
             notes="上半來自 homocysteine，下半（CH₂–CH–COO⁻）來自 serine。cystathionine β-synthase（PLP）缺乏 → homocystinuria：水晶體向下脫位、類 Marfan 體型、智能障礙、血栓（治療：B6、限 Met、補 Cys、betaine）。"),
        dict(id="alpha_ketobutyrate", name="α-Ketobutyrate", aliases=["2-oxobutanoate", "α-酮丁酸"], smiles="CCC(=O)C(=O)[O-]",
             chain=[COO, {"c": "C", "l": "O", "lb": 2}, {"c": "CH_2"}, {"c": "CH_3"}], tags=["one_carbon", "aa_degradation"],
             notes="cystathionase 的另一產物（homocysteine 的碳骨架）→ propionyl-CoA → succinyl-CoA。Thr 分解也會產生。"),
        dict(id="methyl_thf", name="N5-Methyl-THF", label="N^5-Methyl-THF", aliases=["5-methyltetrahydrofolate", "CH3-THF", "N5-甲基四氫葉酸"],
             smiles=thf.replace(")Nc2c(=O)", ")N(C)c2c(=O)"), formula="N^5-CH_3-THF",
             ring={"rotate": 0}, decorate={"highlight": [{"smarts": "[CH3:1][N:2]", "color": "magenta"}]},
             tags=["one_carbon", "vitamin", "coenzyme"],
             notes="MTHFR（NADPH）由 N⁵,N¹⁰-methylene-THF 不可逆還原而來；唯一的出路是 methionine synthase（B12）。B12 缺乏 → folate 全部卡在 CH₃-THF（methyl-folate trap）→ 功能性 folate 缺乏 → 巨球性貧血。"),
        # ---- heme
        dict(id="aminolevulinate", name="δ-Aminolevulinate (ALA)", label="δ-Aminolevulinate\n(ALA)", aliases=["ALA", "5-aminolevulinate", "δ-胺基乙醯丙酸"],
             smiles="[O-]C(=O)CCC(=O)C[NH3+]",
             chain=[COO, {"c": "CH_2"}, {"c": "CH_2"}, {"c": "C", "r": "O", "rb": 2}, {"c": "CH_2NH_3^+"}], tags=["heme"],
             notes="ALA synthase（粒線體，PLP）由 succinyl-CoA + Gly 縮合並脫羧 → heme 合成的限速步驟（heme 回饋抑制；藥物、酒精誘導 ALAS1）。"),
        dict(id="porphobilinogen", name="Porphobilinogen (PBG)", label="Porphobilinogen\n(PBG)", aliases=["PBG", "膽色素原"],
             smiles="[NH3+]Cc1[nH]cc(CCC(=O)[O-])c1CC(=O)[O-]", ring={"rotate": 0}, tags=["heme"],
             notes="2 ALA 經 ALA dehydratase（Zn²⁺，鉛抑制）縮合成的吡咯（pyrrole）。4 個 PBG → hydroxymethylbilane（PBG deaminase；缺乏 → 急性間歇性紫質症 AIP：腹痛、神經精神症狀、尿 PBG/ALA ↑、無光敏感）。"),
        dict(id="protoporphyrin", name="Protoporphyrin IX", aliases=["protoporphyrin", "原紫質"], smiles=PROTOPORPHYRIN,
             ring={"methyls": False}, tags=["heme"],
             notes="4 個吡咯以 =CH– 相連成共軛大環；側鏈：4 CH₃、2 乙烯基、2 丙酸。ferrochelatase（粒線體，鉛抑制）放入 Fe²⁺ → heme。"),
        dict(id="heme", name="Heme", aliases=["heme b", "Fe-protoporphyrin IX", "血基質"], smiles=heme_smiles(),
             ring={"methyls": False, "hide_charges": ["N", "Fe"]}, tags=["heme"],
             notes="Fe²⁺–protoporphyrin IX（圖中略去形式電荷）。Hb、Mb、cytochrome、catalase、NOS 的輔基。Fe 另外兩個配位：近端 His（F8）與 O₂。"),
        dict(id="bilirubin", name="Bilirubin", aliases=["unconjugated bilirubin", "膽紅素"],
             smiles="CC1=C(C=C)C(=O)N/C1=C\\c1[nH]c(Cc2[nH]c(/C=C3\\NC(=O)C(C=C)=C3C)c(C)c2CCC(=O)[O-])c(CCC(=O)[O-])c1C",
             ring={"methyls": False}, tags=["heme"],
             notes="heme oxygenase 打開環（放出 CO 與 Fe）→ biliverdin（綠）→ biliverdin reductase（NADPH）→ bilirubin（黃，(4Z,15Z)）。分子內氫鍵使其不溶於水 → 接 albumin 運送；肝 UGT1A1 接 2 個 glucuronate 後才能排入膽汁。"),
    ]


ATP_ADP = dict(cofactors_in=["ATP"], cofactors_out=["ADP"])
REACTIONS = [
    dict(id="phenylalanine_hydroxylase", title="Phenylalanine → Tyrosine (phenylalanine hydroxylase)",
         substrates=["phenylalanine"], products=["tyrosine"], cofactors_in=["O2", "BH4"], cofactors_out=["H2O", "BH2"],
         enzyme="phenylalanine\nhydroxylase", card={"location": "肝（細胞質）", "type": "羥化（單氧酶：一個 O 進 Tyr、一個變成 H₂O）",
               "cofactors": "BH4（tetrahydrobiopterin）→ BH2，由 dihydropteridine reductase（NADH）再生",
               "deficiency": ["PKU（PAH 缺乏）：Phe ↑ → phenylpyruvate、phenyllactate、phenylacetate（鼠臭味）；智能障礙、淺色皮膚（Tyr ↓ → 黑色素 ↓）；治療：低 Phe 飲食、補 Tyr、避免 aspartame",
                              "BH4 缺乏（DHPR 等）：惡性 PKU，也影響 dopamine 與 serotonin 合成",
                              "母體 PKU：胎兒小頭、心臟缺陷"],
               "exam_points": ["Tyr 在 PKU 變成必需胺基酸", "新生兒篩檢項目"]},
         tags=["aa_derivatives", "aa_degradation"]),
    dict(id="tyrosine_hydroxylase", title="Tyrosine → L-DOPA (tyrosine hydroxylase)",
         substrates=["tyrosine"], products=["l_dopa"], cofactors_in=["O2", "BH4"], cofactors_out=["H2O", "BH2"],
         enzyme="tyrosine hydroxylase", card={"location": "腎上腺髓質、交感神經、多巴胺神經元", "type": "羥化（BH4、Fe²⁺）",
               "regulation": "兒茶酚胺合成的限速步驟；兒茶酚胺回饋抑制，PKA 磷酸化活化"},
         tags=["aa_derivatives"]),
    dict(id="dopa_decarboxylase", title="L-DOPA → Dopamine (aromatic L-amino acid decarboxylase)",
         substrates=["l_dopa"], products=["dopamine", {"id": "co2", "view": "formula"}], extra_in=["H+"], prosthetic=["PLP"],
         enzyme="DOPA decarboxylase\n(AADC)", decorate={"l_dopa": {"colors": {0: "magenta"}}},
         card={"type": "脫羧（PLP）", "exam_points": ["同一個酵素也把 5-HTP 變成 serotonin", "carbidopa 抑制周邊 AADC（不過 BBB）"]},
         tags=["aa_derivatives"]),
    dict(id="dopamine_beta_hydroxylase", title="Dopamine → Norepinephrine (dopamine β-hydroxylase)",
         substrates=["dopamine"], products=["norepinephrine"], cofactors_in=["O2", "ascorbate"], cofactors_out=["H2O", "dehydroascorbate"],
         enzyme="dopamine\nβ-hydroxylase", decorate={"norepinephrine": {"colors": {1: {"r": "magenta"}}}},
         card={"location": "儲存囊泡內", "type": "羥化（Cu²⁺、vitamin C）", "exam_points": ["含銅酵素（Menkes 病會受影響）"]},
         tags=["aa_derivatives"]),
    dict(id="pnmt", title="Norepinephrine → Epinephrine (PNMT)",
         substrates=["norepinephrine"], products=["epinephrine"], cofactors_in=["SAM"], cofactors_out=["SAH"], extra_out=["H+"],
         enzyme="phenylethanolamine\nN-methyltransferase", decorate={"epinephrine": {"colors": {0: {"l": ["magenta", "CH_3"]}}}},
         card={"location": "腎上腺髓質", "type": "N-甲基化（SAM）", "regulation": "cortisol 誘導 PNMT 的表現"},
         tags=["aa_derivatives"]),
    dict(id="tyrosine_aminotransferase", title="Tyrosine + α-KG ⇌ p-Hydroxyphenylpyruvate + Glutamate (TAT)",
         substrates=["tyrosine", "alpha_kg"], products=["hydroxyphenylpyruvate", "glutamate"], reversible=True, prosthetic=["PLP"],
         enzyme="tyrosine\naminotransferase", card={"location": "肝", "type": "轉胺（PLP）", "deficiency": "tyrosinemia II（眼、皮膚病灶）"},
         tags=["aa_degradation"]),
    dict(id="hppd", title="p-Hydroxyphenylpyruvate → Homogentisate (HPP dioxygenase)",
         substrates=["hydroxyphenylpyruvate"], products=["homogentisate", {"id": "co2", "view": "formula"}], cofactors_in=["O2"],
         enzyme="p-hydroxyphenylpyruvate\ndioxygenase", card={"type": "雙氧酶（Fe²⁺、vitamin C）：脫羧＋羥化＋側鏈移位",
               "exam_points": ["nitisinone 抑制此酶，用於治療 tyrosinemia I（fumarylacetoacetate hydrolase 缺乏）"]},
         tags=["aa_degradation"]),
    dict(id="tryptophan_hydroxylase", title="Tryptophan → 5-Hydroxytryptophan (tryptophan hydroxylase)",
         substrates=["tryptophan"], products=["five_htp"], cofactors_in=["O2", "BH4"], cofactors_out=["H2O", "BH2"],
         enzyme="tryptophan\nhydroxylase", card={"type": "羥化（BH4）", "regulation": "serotonin 合成的限速步驟"},
         tags=["aa_derivatives"]),
    dict(id="htp_decarboxylase", title="5-Hydroxytryptophan → Serotonin (AADC)",
         substrates=["five_htp"], products=["serotonin", {"id": "co2", "view": "formula"}], extra_in=["H+"], prosthetic=["PLP"],
         enzyme="aromatic L-amino acid\ndecarboxylase", decorate={"five_htp": {"colors": {0: "magenta"}}},
         card={"type": "脫羧（PLP）"}, tags=["aa_derivatives"]),
    dict(id="histidine_decarboxylase", title="Histidine → Histamine (histidine decarboxylase)",
         substrates=["histidine"], products=["histamine", {"id": "co2", "view": "formula"}], extra_in=["H+"], prosthetic=["PLP"],
         enzyme="histidine\ndecarboxylase", decorate={"histidine": {"colors": {0: "magenta"}}},
         card={"location": "肥大細胞、嗜鹼性球、胃 ECL 細胞", "type": "脫羧（PLP）"}, tags=["aa_derivatives"]),
    dict(id="glutamate_decarboxylase", title="Glutamate → GABA (glutamate decarboxylase)",
         substrates=["glutamate"], products=["gaba", {"id": "co2", "view": "formula"}], extra_in=["H+"], prosthetic=["PLP"],
         enzyme="glutamate\ndecarboxylase (GAD)", decorate={"glutamate": {"colors": {0: "magenta"}}},
         card={"location": "GABA 神經元、胰島 β 細胞", "type": "脫羧（PLP）",
               "exam_points": ["anti-GAD 抗體：第 1 型糖尿病、stiff-person syndrome"]},
         tags=["aa_derivatives"]),
    dict(id="nitric_oxide_synthase", title="2 Arginine → 2 Citrulline + 2 NO (nitric oxide synthase)",
         substrates=[{"id": "arginine", "n": 2}], products=[{"id": "citrulline", "n": 2}, {"id": "nitric_oxide", "n": 2}],
         cofactors_in=["3NADPH", "4O2"], cofactors_out=["3NADP+", "4H2O"], extra_in=["H+"],
         enzyme="nitric oxide synthase", card={"type": "五電子氧化（FAD、FMN、BH4、heme；NADPH、O₂）",
               "regulation": "eNOS、nNOS：Ca²⁺–calmodulin 活化；iNOS：細胞激素誘導（巨噬細胞，殺菌）",
               "exam_points": ["NO 的 N 來自 Arg 的胍基（guanidino）", "NO → guanylyl cyclase → cGMP → 血管擴張"]},
         tags=["aa_derivatives"]),
    dict(id="methionine_adenosyltransferase", title="Methionine + ATP → SAM (methionine adenosyltransferase)",
         substrates=["methionine", "atp_n"], products=["sam"], cofactors_out=["PPi", "Pi"], cofactors_in=["H2O"],
         enzyme="methionine\nadenosyltransferase",
         card={"type": "ATP 的腺苷接到 S 上，三個磷酸全部離開（PPi + Pi）", "exam_points": ["SAM 是主要的甲基供體"]},
         tags=["one_carbon"]),
    dict(id="sah_hydrolase", title="SAH → Homocysteine + Adenosine (SAH hydrolase)",
         substrates=["sah"], products=["homocysteine", "adenosine"], cofactors_in=["H2O"], reversible=True,
         enzyme="S-adenosylhomocysteine\nhydrolase", card={"type": "水解"}, tags=["one_carbon"]),
    dict(id="methionine_synthase", title="Homocysteine + N5-Methyl-THF → Methionine + THF (methionine synthase)",
         substrates=["homocysteine", {"id": "methyl_thf", "view": "formula"}], products=["methionine", {"id": "thf", "view": "formula"}],
         prosthetic=["B_{12}"], enzyme="methionine synthase", decorate={"methionine": {"colors": {5: "magenta"}}},
         card={"location": "細胞質", "type": "甲基轉移（methylcobalamin）",
               "deficiency": ["B12 缺乏 → Hcy ↑ 且 CH₃-THF 卡住（methyl-folate trap）→ 巨球性貧血；另有 MMA ↑ 與神經症狀",
                              "folate 缺乏 → Hcy ↑，MMA 正常"],
               "exam_points": ["連接 folate 循環與甲硫胺酸循環的唯一反應"]},
         tags=["one_carbon"]),
    dict(id="mthfr", title="N5,N10-Methylene-THF → N5-Methyl-THF (MTHFR)",
         substrates=[{"id": "methylene_thf", "view": "formula"}], products=[{"id": "methyl_thf", "view": "formula"}],
         cofactors_in=["NADPH+H+"], cofactors_out=["NADP+"], enzyme="methylene-THF reductase\n(MTHFR)",
         card={"type": "還原（不可逆）", "deficiency": "MTHFR C677T 多型 → 輕度高同半胱胺酸血症"}, tags=["one_carbon"]),
    dict(id="cystathionine_beta_synthase", title="Homocysteine + Serine → Cystathionine (cystathionine β-synthase)",
         substrates=["homocysteine", "serine"], products=["cystathionine"], cofactors_out=["H2O"], prosthetic=["PLP"],
         enzyme="cystathionine\nβ-synthase",
         card={"location": "肝", "type": "縮合（PLP）", "regulation": "SAM 活化（Met 多時把 Hcy 導向分解）",
               "deficiency": "homocystinuria（最常見的型）：水晶體向下脫位、類 Marfan、智能障礙、血栓；部分對 B6 有反應"},
         tags=["one_carbon"]),
    dict(id="cystathionase", title="Cystathionine → Cysteine + α-Ketobutyrate + NH₄⁺ (cystathionine γ-lyase)",
         substrates=["cystathionine"], products=["cysteine", "alpha_ketobutyrate", {"id": "ammonium", "view": "formula"}],
         cofactors_in=["H2O"], prosthetic=["PLP"], enzyme="cystathionine\nγ-lyase",
         card={"type": "裂解（PLP）", "exam_points": ["Cys 的 S 來自 Met，碳骨架來自 Ser → Cys 是條件必需胺基酸", "cystathioninuria 為良性"]},
         tags=["one_carbon"]),
    dict(id="serine_hydroxymethyltransferase", title="Serine + THF ⇌ Glycine + N5,N10-Methylene-THF (SHMT)",
         substrates=["serine", {"id": "thf", "view": "formula"}], products=["glycine", {"id": "methylene_thf", "view": "formula"}],
         cofactors_out=["H2O"], reversible=True, prosthetic=["PLP"], enzyme="serine hydroxymethyl-\ntransferase (SHMT)",
         decorate={"serine": {"colors": {2: "magenta"}}, "methylene_thf": {"colors": {0: {"c": ["magenta", "CH_2"]}}}},
         card={"type": "一碳單位轉移（PLP）", "exam_points": ["THF 一碳單位的主要來源（Ser 的 β 碳）"]},
         tags=["one_carbon"]),
    dict(id="ala_synthase", title="Glycine + Succinyl-CoA → δ-Aminolevulinate (ALA synthase)",
         substrates=["succinyl_coa", "glycine"], products=["aminolevulinate", {"id": "co2", "view": "formula"}],
         cofactors_out=["CoA-SH"], extra_in=["H+"], prosthetic=["PLP"], enzyme="ALA synthase",
         decorate={"glycine": {"colors": {1: "magenta"}}, "aminolevulinate": {"colors": {4: "magenta"}}},
         card={"location": "粒線體基質", "type": "縮合＋脫羧（PLP）",
               "regulation": "限速。heme（與 hemin）抑制合成與運入粒線體；ALAS1（肝）受藥物、酒精、低醣誘導；ALAS2（紅血球前驅）受鐵調控",
               "deficiency": ["B6 缺乏、isoniazid → 鐵粒母細胞性貧血", "ALAS2 突變 → X 性聯鐵粒母細胞性貧血"],
               "exam_points": ["治療急性紫質症：給 hemin 與葡萄糖（抑制 ALAS1）"]},
         tags=["heme"]),
    dict(id="ala_dehydratase", title="2 δ-Aminolevulinate → Porphobilinogen (ALA dehydratase)",
         substrates=[{"id": "aminolevulinate", "n": 2}], products=["porphobilinogen"], cofactors_out=["2H2O"], extra_out=["H+"],
         enzyme="ALA dehydratase\n(PBG synthase)",
         card={"location": "細胞質", "type": "縮合（Zn²⁺）", "deficiency": ["鉛中毒：抑制 ALAD 與 ferrochelatase → ALA ↑、紅血球 Zn-protoporphyrin ↑、小球性貧血、嗜鹼性點彩"]},
         tags=["heme"]),
]


def n(id_, ref, x, y, **kw):
    return {"id": id_, "ref": ref, "pos": [x, y], **kw}


def t(id_, text, x, y, kind="concept"):
    return {"id": id_, "text": text, "kind": kind, "pos": [x, y]}


NAME = {"view": "name"}
PATHWAYS = [
    dict(id="catecholamine_synthesis", title="Catecholamine synthesis（兒茶酚胺合成）", layout="free", cell=0.85,
         nodes=[n("phe", "phenylalanine", 0, 0), n("tyr", "tyrosine", 1.15, 0), n("dopa", "l_dopa", 2.3, 0),
                n("da", "dopamine", 2.3, 1.4), n("ne", "norepinephrine", 1.15, 1.4), n("epi", "epinephrine", 0, 1.4)],
         edges=[{"from": "phe", "to": "tyr", "step": "phenylalanine_hydroxylase", "label": "PAH", "side_text": "BH_4"},
                {"from": "tyr", "to": "dopa", "step": "tyrosine_hydroxylase", "label": "tyrosine hydroxylase\n(rate-limiting)", "side_text": "BH_4"},
                {"from": "dopa", "to": "da", "step": "dopa_decarboxylase", "label": "AADC", "side_text": "PLP, −CO_2", "label_side": "right"},
                {"from": "da", "to": "ne", "step": "dopamine_beta_hydroxylase", "label": "DBH", "side_text": "vit C, Cu"},
                {"from": "ne", "to": "epi", "step": "pnmt", "label": "PNMT", "side_text": "SAM"}],
         tags=["aa_derivatives"],
         notes="Phe → Tyr（PAH）→ DOPA（TH，限速，BH4）→ dopamine（AADC，B6）→ NE（DBH，Cu、vit C）→ Epi（PNMT，SAM；cortisol 誘導）。分解：MAO、COMT → HVA（dopamine）、VMA（NE、Epi）。"),
    dict(id="tyrosine_metabolism", title="Phenylalanine & tyrosine metabolism（名稱版，含缺陷疾病）", layout="free", cell=0.62,
         nodes=[n("phe", "phenylalanine", 0, 0, **NAME), n("tyr", "tyrosine", 0, 1, **NAME),
                n("hpp", "hydroxyphenylpyruvate", 0, 2, **NAME), n("hga", "homogentisate", 0, 3, **NAME),
                t("maa", "Maleylacetoacetate", 0, 4, "name"), t("faa", "Fumarylacetoacetate", 0, 5, "name"),
                n("fum", "fumarate", -1.0, 6, **NAME), n("aa", "acetoacetate", 1.0, 6, **NAME),
                n("dopa", "l_dopa", 2.4, 1, **NAME), t("mel", "Melanin", 4.0, 1, "name"),
                t("cat", "Catecholamines", 2.4, 2, "name"), t("thy", "Thyroid hormones\n(T_3, T_4)", -2.3, 1, "name"),
                t("pku", "PKU", -1.6, 0.5, "drug"), t("alk", "Alkaptonuria", -1.6, 3.5, "drug"),
                t("tyr1", "Tyrosinemia I", 1.8, 5.5, "drug"), t("alb", "Albinism", 3.3, 1.5, "drug")],
         edges=[{"from": "phe", "to": "tyr", "step": "phenylalanine_hydroxylase", "label": "PAH (BH_4)"},
                {"from": "tyr", "to": "hpp", "step": "tyrosine_aminotransferase", "label": "TAT (PLP)"},
                {"from": "hpp", "to": "hga", "step": "hppd", "label": "HPP dioxygenase (vit C)"},
                {"from": "hga", "to": "maa", "label": "homogentisate oxidase"},
                {"from": "maa", "to": "faa", "label": "isomerase"},
                {"from": "faa", "to": "fum", "label": "FAH"}, {"from": "faa", "to": "aa"},
                {"from": "tyr", "to": "dopa", "step": "tyrosine_hydroxylase", "label": "TH"},
                {"from": "dopa", "to": "mel", "label": "tyrosinase (Cu)"},
                {"from": "dopa", "to": "cat"}, {"from": "tyr", "to": "thy"},
                {"from": "pku", "to": "edge:phe->tyr", "inhibit": True},
                {"from": "alk", "to": "edge:hga->maa", "inhibit": True},
                {"from": "tyr1", "to": "edge:faa->fum", "inhibit": True},
                {"from": "alb", "to": "edge:dopa->mel", "inhibit": True}],
         tags=["aa_degradation", "aa_derivatives"],
         notes="Phe、Tyr 同時生糖（fumarate）也生酮（acetoacetate）。缺陷：PKU（PAH 或 BH4）、alkaptonuria（homogentisate oxidase：黑尿、褐黃病）、tyrosinemia I（FAH：肝衰竭、腎小管病變、煮包心菜味）、albinism（tyrosinase）。"),
    dict(id="tryptophan_metabolism", title="Tryptophan metabolism（名稱版）", layout="free", cell=0.62,
         nodes=[n("trp", "tryptophan", 0, 0, **NAME), n("htp", "five_htp", -1.4, 1, **NAME), n("ser", "serotonin", -1.4, 2, **NAME),
                t("nas", "N-Acetylserotonin", -1.4, 3, "name"), t("mel", "Melatonin\n(pineal)", -1.4, 4, "name"),
                t("hiaa", "5-HIAA (urine)", -3.2, 2, "name"),
                t("kyn", "N-Formylkynurenine →\nkynurenine", 1.4, 1, "name"), t("nad", "NAD^+ (niacin)", 1.4, 2.2, "name"),
                t("ac", "Alanine + acetyl-CoA\n(gluco- & ketogenic)", 3.4, 1.6, "name")],
         edges=[{"from": "trp", "to": "htp", "step": "tryptophan_hydroxylase", "label": "Trp hydroxylase\n(BH_4)"},
                {"from": "htp", "to": "ser", "step": "htp_decarboxylase", "label": "AADC (PLP)"},
                {"from": "ser", "to": "nas", "label": "N-acetyltransferase"},
                {"from": "nas", "to": "mel", "label": "O-methyltransferase", "side_text": "SAM"},
                {"from": "ser", "to": "hiaa", "label": "MAO"},
                {"from": "trp", "to": "kyn", "label": "Trp 2,3-dioxygenase\n(IDO, TDO)"},
                {"from": "kyn", "to": "nad", "label": "B_6, B_2"},
                {"from": "kyn", "to": "ac"}],
         tags=["aa_derivatives", "aa_degradation"],
         notes="Trp → serotonin → melatonin；Trp → kynurenine → NAD⁺（60 mg Trp ≈ 1 mg niacin，需 B6、B2）。Hartnup（中性胺基酸吸收障礙）與 carcinoid（Trp 被大量拿去做 serotonin）→ pellagra。"),
    dict(id="amino_acid_derivatives", title="Specialized products of amino acids（胺基酸衍生物，名稱版）", layout="free", cell_size=6.0, text=False,
         nodes=[n("his", "histidine", 0, 0, **NAME), n("hist", "histamine", 1.6, 0, **NAME),
                n("glu", "glutamate", 0, 1, **NAME), n("gaba", "gaba", 1.6, 1, **NAME),
                n("arg", "arginine", 0, 2, **NAME), n("no", "nitric_oxide", 1.6, 2, **NAME),
                t("cre", "Creatine (Gly + Arg + SAM)\n→ creatine-P, creatinine", 3.6, 2.5, "name"),
                n("trp", "tryptophan", 0, 3, **NAME), n("sero", "serotonin", 1.6, 3, **NAME),
                n("tyr", "tyrosine", 0, 4, **NAME), t("cats", "Catecholamines, melanin,\nT_3/T_4", 1.6, 4, "name"),
                n("gly", "glycine", 0, 5, **NAME), n("heme", "heme", 1.6, 5, **NAME),
                t("gsh", "Glutathione (γ-Glu–Cys–Gly)", 3.6, 5.5, "name"),
                t("ser", "Serine → sphingosine,\nphosphatidylserine", 3.6, 0.5, "name")],
         edges=[{"from": "his", "to": "hist", "step": "histidine_decarboxylase", "label": "HDC (PLP)"},
                {"from": "glu", "to": "gaba", "step": "glutamate_decarboxylase", "label": "GAD (PLP)"},
                {"from": "arg", "to": "no", "step": "nitric_oxide_synthase", "label": "NOS"},
                {"from": "arg", "to": "cre", "bend": -0.15},
                {"from": "trp", "to": "sero", "label": "BH_4, PLP"},
                {"from": "tyr", "to": "cats", "label": "BH_4, PLP, Cu"},
                {"from": "gly", "to": "heme", "label": "ALA synthase"},
                {"from": "gly", "to": "gsh", "bend": 0.15}],
         tags=["aa_derivatives"],
         notes="胺基酸脫羧多半需要 PLP（B6）：histamine、GABA、serotonin、dopamine。Gly 是 heme、creatine、glutathione、purine 的原料；Arg → NO、creatine、urea、polyamines（經 ornithine）。"),
    dict(id="methionine_cycle", title="Methionine cycle & one-carbon metabolism（甲硫胺酸循環與一碳代謝）", layout="free", cell=0.5,
         nodes=[n("met", "methionine", 0, 0, **NAME), n("sam", "sam", 2.0, 0, **NAME), n("sah", "sah", 2.0, 1.4, **NAME),
                n("hcy", "homocysteine", 0, 1.4, **NAME), n("cth", "cystathionine", 0, 2.6, **NAME),
                n("cys", "cysteine", -0.8, 3.7, **NAME), n("akb", "alpha_ketobutyrate", 0.8, 3.7, **NAME),
                n("mthf", "methyl_thf", -1.6, 0.7, view="formula"), n("thf", "thf", -1.6, -0.5, view="formula"),
                n("ch2", "methylene_thf", -3.1, 0.1, view="formula"),
                t("acc", "Acceptor → CH_3–acceptor\n(DNA, NE → Epi, creatine…)", 3.6, 0.7, "plain")],
         edges=[{"from": "met", "to": "sam", "step": "methionine_adenosyltransferase", "label": "MAT", "side_text": "ATP → PP_i + P_i"},
                {"from": "sam", "to": "sah", "label": "methyltransferases", "label_side": "right"},
                {"from": "sah", "to": "hcy", "step": "sah_hydrolase", "label": "SAH hydrolase", "side_text": "→ adenosine"},
                {"from": "hcy", "to": "met", "step": "methionine_synthase", "label": "methionine synthase\n(B_{12})"},
                {"from": "mthf", "to": "edge:hcy->met"}, {"from": "edge:hcy->met", "to": "thf"},
                {"from": "thf", "to": "ch2", "step": "serine_hydroxymethyltransferase", "label": "SHMT", "side_text": "Ser → Gly"},
                {"from": "ch2", "to": "mthf", "step": "mthfr", "label": "MTHFR", "side_text": "NADPH"},
                {"from": "hcy", "to": "cth", "step": "cystathionine_beta_synthase", "label": "CBS (B_6)", "side_text": "Ser"},
                {"from": "cth", "to": "cys", "step": "cystathionase", "label": "γ-lyase (B_6)"},
                {"from": "cth", "to": "akb", "step": "cystathionase"},
                {"from": "sam", "to": "acc", "head": False}],
         tags=["one_carbon"],
         notes="Met → SAM（甲基供體）→ SAH → Hcy。Hcy 再甲基化（methionine synthase：B12 + CH₃-THF）或經 transsulfuration（CBS、γ-lyase：B6）變成 Cys。Hcy ↑：B12、folate、B6 缺乏、CBS 或 MTHFR 缺陷。B12 缺乏 → methyl-folate trap。"),
    dict(id="amino_acid_carbon_skeletons", title="Amino acid carbon skeletons → TCA（胺基酸碳骨架的去路）", layout="free", cell_size=5.6, text=False,
         nodes=[n("pyr", "pyruvate", 0, 0, **NAME), n("acc", "acetyl_coa", 1.6, 0, **NAME), t("aacoa", "Acetoacetyl-CoA", 3.4, 0, "name"),
                n("oaa", "oxaloacetate", 0, 1.6, **NAME), n("cit", "citrate", 1.6, 1.0, **NAME), n("akg", "alpha_kg", 3.0, 2.0, **NAME),
                n("suc", "succinyl_coa", 2.4, 3.3, **NAME), n("fum", "fumarate", 0.8, 3.3, **NAME),
                t("tca", "TCA\ncycle", 1.6, 2.1, "plain"),
                t("a_pyr", "Ala, Cys, Gly,\nSer, Thr, Trp", -1.8, 0, "plain"),
                t("a_ac", "Ile, Leu, Thr, Trp", 1.6, -0.9, "plain"),
                t("a_aac", "Leu, Lys, Phe,\nTrp, Tyr", 3.4, -0.9, "plain"),
                t("a_oaa", "Asn, Asp", -1.6, 1.6, "plain"),
                t("a_akg", "Arg, Glu, Gln,\nHis, Pro", 4.8, 2.0, "plain"),
                t("a_suc", "Ile, Met,\nThr, Val", 3.6, 4.2, "plain"),
                t("a_fum", "Phe, Tyr", -0.6, 4.2, "plain"),
                t("ket", "Ketogenic only: Leu, Lys", 3.4, 5.2, "concept"),
                t("both", "Both: Ile, Phe, Thr, Trp, Tyr", 0.0, 5.2, "concept")],
         edges=[{"from": "a_pyr", "to": "pyr"}, {"from": "a_ac", "to": "acc"}, {"from": "a_aac", "to": "aacoa"},
                {"from": "aacoa", "to": "acc"}, {"from": "pyr", "to": "acc", "step": "pdh"},
                {"from": "pyr", "to": "oaa", "step": "pyruvate_carboxylase"},
                {"from": "a_oaa", "to": "oaa"}, {"from": "a_akg", "to": "akg"}, {"from": "a_suc", "to": "suc"}, {"from": "a_fum", "to": "fum"},
                {"from": "oaa", "to": "cit", "step": "citrate_synthase"}, {"from": "acc", "to": "cit"},
                {"from": "cit", "to": "akg"}, {"from": "akg", "to": "suc", "step": "alpha_kg_dehydrogenase"},
                {"from": "suc", "to": "fum"}, {"from": "fum", "to": "oaa"}],
         tags=["aa_degradation"],
         notes="生糖（glucogenic）：進入 pyruvate 或 TCA 中間物 → OAA → 葡萄糖。生酮（ketogenic）：只進 acetyl-CoA／acetoacetyl-CoA。只生酮：Leu、Lys。兩者皆是：Ile、Phe、Thr、Trp、Tyr。其餘只生糖。"),
    dict(id="bcaa_catabolism", title="Branched-chain amino acid catabolism（支鏈胺基酸分解，名稱版）", layout="free", cell_size=6.0, text=False,
         nodes=[n("leu", "leucine", 0, 0, **NAME), n("ile", "isoleucine", 1.5, 0, **NAME), n("val", "valine", 3.0, 0, **NAME),
                t("kic", "α-Ketoisocaproate", 0, 1, "name"), t("kmv", "α-Keto-β-methylvalerate", 1.5, 1, "name"), t("kiv", "α-Ketoisovalerate", 3.0, 1, "name"),
                t("ivc", "Isovaleryl-CoA", 0, 2.2, "name"), t("mbc", "α-Methylbutyryl-CoA", 1.5, 2.2, "name"), t("ibc", "Isobutyryl-CoA", 3.0, 2.2, "name"),
                n("ac", "acetyl_coa", 0.6, 3.4, **NAME), n("aa", "acetoacetate", -0.9, 3.4, **NAME),
                n("pro", "propionyl_coa", 2.4, 3.4, **NAME), n("suc", "succinyl_coa", 2.4, 4.4, **NAME),
                t("msud", "MSUD\n(BCKDH: all three)", 4.6, 1.6, "drug")],
         edges=[{"from": "leu", "to": "kic", "label": "BCAT (PLP)"}, {"from": "ile", "to": "kmv"}, {"from": "val", "to": "kiv"},
                {"from": "kic", "to": "ivc", "label": "BCKDH"}, {"from": "kmv", "to": "mbc", "label": "BCKDH"}, {"from": "kiv", "to": "ibc", "label": "BCKDH"},
                {"from": "ivc", "to": "ac"}, {"from": "ivc", "to": "aa"}, {"from": "mbc", "to": "ac"}, {"from": "mbc", "to": "pro"},
                {"from": "ibc", "to": "pro"}, {"from": "pro", "to": "suc", "label": "B_7, B_{12}"},
                {"from": "msud", "to": "edge:kiv->ibc", "inhibit": True}],
         tags=["aa_degradation"],
         notes="主要在肌肉轉胺（BCAT），α-酮酸由 BCKDH（輔因子同 PDH：TPP、lipoate、CoA、FAD、NAD⁺）氧化脫羧。缺乏 → 楓糖尿症（MSUD）：尿有楓糖漿味、餵食困難、神經退化；治療：限制 BCAA、部分對 thiamine 有反應。Leu 只生酮、Val 只生糖、Ile 兩者皆是。"),
    dict(id="heme_synthesis", title="Heme synthesis（血基質合成，含紫質症）", layout="free", cell=0.62,
         nodes=[t("hm", "Mitochondrion", -2.6, -0.6, "heading"), t("hc", "Cytosol", 2.4, -0.6, "heading"),
                n("gly", "glycine", -2.6, 0.2, **NAME), n("suc", "succinyl_coa", -0.9, 0.2, **NAME),
                n("ala", "aminolevulinate", -1.6, 1.3, **NAME), n("ala_c", "aminolevulinate", 2.4, 1.3, **NAME),
                n("pbg", "porphobilinogen", 2.4, 2.5, **NAME), t("hmb", "Hydroxymethylbilane", 2.4, 3.7, "name"),
                t("uro", "Uroporphyrinogen III", 2.4, 4.9, "name"), t("copro", "Coproporphyrinogen III", 2.4, 6.1, "name"),
                t("ppgen", "Protoporphyrinogen IX", -1.6, 6.1, "name"), n("pp", "protoporphyrin", -1.6, 4.9, **NAME),
                n("heme", "heme", -1.6, 3.7, **NAME),
                t("lead", "Lead", 0.2, 2.0, "drug"), t("lead2", "Lead", 0.2, 4.3, "drug"),
                t("aip", "AIP", 4.2, 3.1, "drug"), t("pct", "PCT", 4.2, 5.5, "drug"), t("cep", "CEP", 4.2, 4.3, "drug"),
                t("epp", "EPP", -3.4, 4.3, "drug")],
         edges=[{"from": "gly", "to": "ala", "step": "ala_synthase", "label": "ALA synthase\n(PLP, rate-limiting)"},
                {"from": "suc", "to": "ala", "step": "ala_synthase"},
                {"from": "ala", "to": "ala_c"},
                {"from": "ala_c", "to": "pbg", "step": "ala_dehydratase", "label": "ALA dehydratase (×2)", "label_side": "right"},
                {"from": "pbg", "to": "hmb", "label": "PBG deaminase (×4)", "label_side": "right"},
                {"from": "hmb", "to": "uro", "label": "UROS", "label_side": "right"},
                {"from": "uro", "to": "copro", "label": "UROD", "label_side": "right"},
                {"from": "copro", "to": "ppgen", "label": "CPOX"},
                {"from": "ppgen", "to": "pp", "label": "PPOX"},
                {"from": "pp", "to": "heme", "label": "ferrochelatase", "side_text": "Fe^{2+}"},
                {"from": "heme", "to": "edge:gly->ala", "inhibit": True},
                {"from": "lead", "to": "edge:ala_c->pbg", "inhibit": True}, {"from": "lead2", "to": "edge:pp->heme", "inhibit": True},
                {"from": "aip", "to": "edge:pbg->hmb", "inhibit": True}, {"from": "cep", "to": "edge:hmb->uro", "inhibit": True},
                {"from": "pct", "to": "edge:uro->copro", "inhibit": True}, {"from": "epp", "to": "edge:pp->heme", "inhibit": True}],
         tags=["heme"],
         notes="頭尾兩步在粒線體，中間在細胞質；骨髓（約 85%）與肝。紫質症：AIP（PBG deaminase；急性腹痛、神經精神、紅尿，無光敏感，藥物／酒精／飢餓誘發，治療 hemin + 葡萄糖）、PCT（UROD；最常見，光敏感水泡，與 HCV、酒精、鐵過多有關）、CEP（UROS）、EPP（ferrochelatase）。鉛：抑制 ALAD 與 ferrochelatase。"),
    dict(id="heme_degradation", title="Heme degradation & bilirubin（血基質分解與膽紅素）", layout="free", cell=0.62,
         nodes=[t("hs", "Spleen / macrophage", -0.2, -0.7, "heading"), n("heme", "heme", 0, 0, **NAME),
                t("bv", "Biliverdin (green)", 0, 1, "name"), n("br", "bilirubin", 0, 2, **NAME),
                t("alb", "Blood: bilirubin–albumin\n(unconjugated, indirect)", 0, 3.1, "name"),
                t("hl", "Liver", 2.6, 2.4, "heading"), t("dg", "Bilirubin diglucuronide\n(conjugated, direct)", 2.6, 3.1, "name"),
                t("hi", "Intestine", 2.6, 4.3, "heading"), t("ubg", "Urobilinogen", 2.6, 5.0, "name"),
                t("ster", "Stercobilin (feces, brown)", 4.6, 5.9, "name"), t("urob", "Urobilin (urine, yellow)", 0.6, 5.9, "name"),
                t("gil", "Gilbert, Crigler–Najjar", 4.8, 2.4, "drug")],
         edges=[{"from": "heme", "to": "bv", "label": "heme oxygenase", "side_text": "O_2, NADPH → CO + Fe^{2+}"},
                {"from": "bv", "to": "br", "label": "biliverdin reductase", "side_text": "NADPH"},
                {"from": "br", "to": "alb"}, {"from": "alb", "to": "dg", "label": "UGT1A1", "side_text": "2 UDP-glucuronate"},
                {"from": "dg", "to": "ubg", "label": "bile → bacteria", "label_side": "right"},
                {"from": "ubg", "to": "ster"}, {"from": "ubg", "to": "urob", "label": "kidney"},
                {"from": "gil", "to": "edge:alb->dg", "inhibit": True}],
         tags=["heme"],
         notes="間接（未結合）膽紅素 ↑：溶血、新生兒、Gilbert（UGT1A1 ↓，壓力或空腹時輕度黃疸）、Crigler–Najjar（I 型完全缺乏 → 核黃疸）。直接（結合）膽紅素 ↑：Dubin–Johnson（MRP2，肝色黑）、Rotor、膽道阻塞（灰白便、茶色尿、尿 urobilinogen ↓）。heme oxygenase 是體內 CO 的唯一來源。"),
]

COFACTORS = {      # name: (label, SMILES, linked molecule)
    "BH4": ('BH_4', "Nc1nc2NCC([13*])Nc2c(=O)[nH]1", None),
    "BH2": ('BH_2', "Nc1nc2NCC([13*])=Nc2c(=O)[nH]1", None),
    "SAM": ("SAM", "C[S+]([14*])[14*]", "sam"),                 # [14*] = adenosyl / homocysteinyl
    "SAH": ("SAH", "[14*]S[14*]", "sah"),
    "ascorbate": ("ascorbate", "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O", "ascorbate"),
    "dehydroascorbate": ("dehydroascorbate", "OC[C@H](O)[C@H]1OC(=O)C(=O)C1=O", "dehydroascorbate"),
}
MULTI = {"3NADPH": ("3 NADPH", "NC(=O)C1=CN([6*])C=CC1", 3),
         "3NADP+": ("3 NADP^+", "NC(=O)c1ccc[n+]([6*])c1", 3),
         "4O2": ("4 O_2", "O=O", 4), "4H2O": ("4 H_2O", "O", 4), "2H2O": ("2 H_2O", "O", 2)}
MARK = "# --- 胺基酸代謝補完（第十輪）"


def dump(path, d):
    path.write_text(yaml.safe_dump(d, allow_unicode=True, sort_keys=False, width=200), encoding="utf-8")


def with_parent(tags):
    out = list(tags)
    if any(t in KIDS for t in out) and CM not in out:
        out.insert(0, CM)
    return out


def add_cofactors(data):
    p = Path(data) / "cofactors.yaml"
    keys = list(COFACTORS) + list(MULTI)
    lines = [ln for ln in p.read_text(encoding="utf-8").splitlines()
             if not ln.startswith(MARK) and not any(ln.startswith(f"{k}:") for k in keys)]
    lines.append(MARK + "：BH4 = tetrahydrobiopterin（[13*] = 側鏈）；SAM/SAH（[14*] = 其餘部分）；多份輔因子（由 tools/gen_aamet.py 產生）")
    for k, (label, smi, mol) in COFACTORS.items():
        link = f", molecule: {mol}" if mol else ""
        lines.append(f'{k}: {{label: "{label}", smiles: "{smi}"{link}}}')
    for k, (label, frag, n_) in MULTI.items():
        lines.append(f'{k}: {{label: "{label}", smiles: "{".".join([frag] * n_)}"}}')
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(data):
    data = Path(data)
    add_cofactors(data)
    mols = molecules(data)
    for m in mols:
        d = {"id": m["id"], "name": m["name"]}
        if m.get("label"):
            d["label"] = m["label"]
        d["aliases"] = m.get("aliases", [])
        d["smiles"] = m["smiles"]
        for k in ("chain", "ring", "formula", "decorate"):
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
    print(f"{len(mols)} molecules, {len(REACTIONS)} reactions, {len(PATHWAYS)} pathways")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data")
