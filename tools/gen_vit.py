"""Vitamins & coenzymes (round 7): molecules, activation reactions and overview sheets.

    .venv/bin/python tools/gen_vit.py data

All SMILES are built from the same templates as the nucleotides (adenosine written base-first
like AMP, so the ribose stereo matches the validated Haworth drawings). Re-running overwrites
the files listed here, including the THF family (thf, dhf, methylene_thf).
"""
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent))
from reaction_colors import COLORS  # noqa: E402

# ---------------------------------------------------------------- SMILES building blocks
ADO_TAIL = "OC[C@H]5O[C@@H](n6cnc7c(N)ncnc76)[C@H](%s)[C@@H]5O"   # 5'-O … adenosine (2'-substituent)
PPI_LINK = "P(=O)([O-])OP(=O)([O-])"
PABA_GLU = "c8ccc(cc8)C(=O)N[C@@H](CCC(=O)[O-])C(=O)[O-]"          # p-aminobenzoyl-L-glutamate
RIBITYL = "C[C@H](O)[C@H](O)[C@H](O)C"                                # D-ribityl (2S,3S,4R), C5' open
FLAVIN_OX = "Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(%s)c2cc1C"
FLAVIN_RED = "Cc1cc2Nc3c(=O)[nH]c(=O)[nH]c3N(%s)c2cc1C"
NIC_OX = "NC(=O)c1ccc[n+](c1)"
NIC_RED = "NC(=O)C1=CN(C=CC1)"
PHYTYL = "C/C=C(\\C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C"


def nucleotide_sugar(base, o5):
    """Ribose written base-first (AMP template): base atom, C1′ … C5′–O–o5."""
    return f"{base}[C@@H]4O[C@H](CO{o5})[C@@H](O)[C@H]4O"


SECO = "C[C@H](CCCC(C)(C)%s)[C@H]1CC[C@@H]\\2[C@@]1(CCC/C2=C\\C=C/3\\C[C@H](%s)O)C"

# ---------------------------------------------------------------- highlight SMARTS
HL_FLAVIN = [{"smarts": "[#7:1]1~[#6:2]2~[#6:3](~[#6]~[#6]~[#6]~[#6]~2)~[#7:4]~[#6:5]~[#6:6]~1~[#7:7]",
              "color": "magenta"}]
HL_NIC_OX = [{"smarts": "NC(=O)c1ccc[n+]c1", "color": "magenta"}]
HL_NIC_RED = [{"smarts": "NC(=O)C1=CNC=CC1", "color": "magenta"}]
HL_2P = [{"smarts": "n[CH]1O[CH][CH](O)[CH]1[O:1][P:2](=[O:3])([O-:4])[O-:5]", "color": "cyan"}]

WS, FS, CO = "water_soluble", "fat_soluble", "coenzyme"

MOLS = [
    # ---- B1
    dict(id="thiamine", name="Thiamine (vitamin B1)", label="Thiamine (B_1)", aliases=["vitamin B1", "維生素B1", "硫胺"],
         smiles="Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1", ring={"methyls": True}, tags=[WS],
         notes="嘧啶環＋噻唑（thiazolium）環以 CH₂ 相連。活化：thiamine pyrophosphokinase + ATP → TPP。缺乏：beriberi（乾性：周邊神經病變；濕性：高輸出心衰竭）、Wernicke–Korsakoff（酒精；給葡萄糖前先補 thiamine）。"),
    dict(id="tpp", name="Thiamine pyrophosphate (TPP)", label="Thiamine\npyrophosphate (TPP)",
         aliases=["TPP", "ThDP", "thiamine diphosphate", "硫胺焦磷酸"],
         smiles="Cc1ncc(C[n+]2csc(CCOP(=O)([O-])OP(=O)([O-])[O-])c2C)c(N)n1",
         ring={"highlight": [{"smarts": "[n+]:[cH]:s", "color": "magenta"}]}, tags=[CO],
         notes="洋紅：噻唑環 C2（夾在 N⁺ 與 S 之間）。C2–H 解離成碳陰離子（ylide），攻擊 α-酮酸的羰基 → 脫羧後以「活性醛」形式攜帶。酵素：PDH（E1）、α-KGDH、BCKDH（MSUD）、transketolase（PPP；紅血球 transketolase 活性可評估 B1 狀態）。"),
    # ---- B2
    dict(id="riboflavin", name="Riboflavin (vitamin B2)", label="Riboflavin (B_2)", aliases=["vitamin B2", "維生素B2", "核黃素"],
         smiles=FLAVIN_OX % (RIBITYL + "O"), ring={"methyls": True}, tags=[WS],
         notes="isoalloxazine 環（黃色、可吸光）＋ ribityl（糖醇，不是 ribose）。缺乏：口角炎（angular cheilitis）、舌炎、角膜血管增生。"),
    dict(id="fmn", name="Flavin mononucleotide (FMN)", label="FMN", aliases=["flavin mononucleotide", "黃素單核苷酸"],
         smiles=FLAVIN_OX % (RIBITYL + "OP(=O)([O-])[O-]"), ring={"highlight": HL_FLAVIN}, tags=[CO],
         notes="riboflavin kinase（flavokinase）+ ATP → FMN。洋紅：中央環與 N1、N5，接受 2 個 H（2 e⁻ + 2 H⁺）的位置；可一次接 1 個 e⁻（semiquinone），所以能銜接 2e⁻ 與 1e⁻ 載體。Complex I 的輔基。"),
    dict(id="fad", name="Flavin adenine dinucleotide (FAD)", label="FAD", aliases=["flavin adenine dinucleotide", "黃素腺嘌呤二核苷酸"],
         smiles=FLAVIN_OX % (RIBITYL + "O" + PPI_LINK + ADO_TAIL % "O"), ring={"highlight": HL_FLAVIN}, tags=[CO],
         notes="FMN + ATP → FAD + PPi（FAD synthetase）。洋紅：氧化還原中心（N1、N5）。常為緊密結合的輔基（prosthetic group）：succinate DH（Complex II，與膜共價結合）、acyl-CoA DH、PDH／α-KGDH 的 E3（dihydrolipoyl DH）。"),
    dict(id="fadh2", name="FADH2", label="FADH_2", aliases=["reduced FAD", "還原態FAD"],
         smiles=FLAVIN_RED % (RIBITYL + "O" + PPI_LINK + ADO_TAIL % "O"), ring={"highlight": HL_FLAVIN}, tags=[CO],
         notes="還原態：N1、N5 各接一個 H（洋紅），中央環不再共軛 → 失去黃色。電子經 Q 進入 ETC，每對電子約 1.5 ATP（NADH 約 2.5）。"),
    # ---- B3
    dict(id="nicotinate", name="Nicotinate (niacin, vitamin B3)", label="Nicotinate (niacin)", aliases=["niacin", "nicotinic acid", "vitamin B3", "菸鹼酸"],
         smiles="[O-]C(=O)c1cccnc1", ring={"rotate": 0}, tags=[WS],
         notes="可由 Trp 合成（約 60 mg Trp → 1 mg niacin，需 B6、B2）。缺乏：pellagra（3D：dermatitis、diarrhea、dementia，可再加 death）；Hartnup 病（中性胺基酸吸收障礙）、carcinoid（Trp 被拿去做 serotonin）、isoniazid（耗 B6）也會造成。大劑量 niacin 降血脂（副作用：潮紅，可用 aspirin 預防）。"),
    dict(id="nicotinamide", name="Nicotinamide", label="Nicotinamide", aliases=["niacinamide", "菸鹼醯胺"],
         smiles="NC(=O)c1cccnc1", ring={"rotate": 0}, tags=[WS],
         notes="NAD⁺ 分解（sirtuin、PARP、ADP-ribosylation）後的產物，可經 salvage 重新合成 NAD⁺。不會造成潮紅。"),
    dict(id="nad", name="NAD+", label="NAD^+", aliases=["nicotinamide adenine dinucleotide", "NAD", "菸鹼醯胺腺嘌呤二核苷酸"],
         smiles=nucleotide_sugar(NIC_OX, PPI_LINK + ADO_TAIL % "O"), ring={"flip": True, "highlight": HL_NIC_OX}, tags=[CO],
         notes="洋紅：nicotinamide 環，C4（醯胺的對位）接受 hydride（H⁻ = 2 e⁻ + H⁺），另一個 H⁺ 釋入溶液 → NADH + H⁺。可溶性輔酶（不是輔基），主要用在分解代謝的氧化反應（NAD⁺/NADH 比值高）。"),
    dict(id="nadh", name="NADH", label="NADH", aliases=["reduced NAD"],
         smiles=nucleotide_sugar(NIC_RED, PPI_LINK + ADO_TAIL % "O"), ring={"flip": True, "highlight": HL_NIC_RED}, tags=[CO],
         notes="還原態：nicotinamide 環的 C4 變成 CH₂（多一個 H），環不再芳香 → 340 nm 吸光（酵素活性測定的依據）。電子經 Complex I 進入 ETC。"),
    dict(id="nadp", name="NADP+", label="NADP^+", aliases=["nicotinamide adenine dinucleotide phosphate", "NADP"],
         smiles=nucleotide_sugar(NIC_OX, PPI_LINK + ADO_TAIL % "OP(=O)([O-])[O-]"),
         ring={"coordgen": False, "seed": 14, "flip": True, "highlight": HL_NIC_OX + HL_2P}, tags=[CO],
         notes="與 NAD⁺ 的差別只有腺苷 2′-OH 上多一個磷酸（青藍）。NAD kinase 由 NAD⁺ 生成。酵素靠這個磷酸分辨 NAD 與 NADP。"),
    dict(id="nadph", name="NADPH", label="NADPH", aliases=["reduced NADP"],
         smiles=nucleotide_sugar(NIC_RED, PPI_LINK + ADO_TAIL % "OP(=O)([O-])[O-]"),
         ring={"coordgen": False, "seed": 14, "flip": True, "highlight": HL_NIC_RED + HL_2P}, tags=[CO],
         notes="合成代謝的還原力（NADPH/NADP⁺ 比值高）：脂肪酸、膽固醇合成，glutathione 還原（紅血球，G6PD 缺乏 → 溶血），CYP450，NADPH oxidase（呼吸爆）。主要來源：PPP 氧化相（G6PD）、malic enzyme。"),
    # ---- B5
    dict(id="pantothenate", name="Pantothenate (vitamin B5)", label="Pantothenate (B_5)", aliases=["pantothenic acid", "vitamin B5", "泛酸"],
         smiles="CC(C)(CO)[C@@H](O)C(=O)NCCC(=O)[O-]", ring={"methyls": False}, tags=[WS],
         notes="pantoate + β-alanine（醯胺鍵）。是 CoA 與 ACP（fatty acid synthase 的醯基載體）的組成。缺乏罕見（burning feet）。"),
    dict(id="coa", name="Coenzyme A (CoA-SH)", label="CoA-SH", aliases=["CoA", "coenzyme A", "輔酶A"],
         smiles="SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]2O[C@@H](n3cnc4c(N)ncnc43)[C@H](O)[C@@H]2OP(=O)([O-])[O-]",
         ring={"methyls": False, "highlight": [{"smarts": "[SH1:1]", "color": "magenta"}]}, tags=[CO],
         notes="β-mercaptoethylamine – pantothenate – 3′-phospho-ADP。洋紅：末端 –SH，與醯基形成高能硫酯（acetyl-CoA、succinyl-CoA、acyl-CoA）。由 pantothenate、Cys、ATP 合成。"),
    # ---- B6
    dict(id="pyridoxine", name="Pyridoxine (vitamin B6)", label="Pyridoxine (B_6)", aliases=["vitamin B6", "維生素B6", "吡哆醇"],
         smiles="Cc1ncc(CO)c(CO)c1O", ring={"rotate": 0}, tags=[WS],
         notes="B6 三種形式：pyridoxine（醇，C4′ = CH₂OH）、pyridoxal（醛）、pyridoxamine（胺）。isoniazid 會耗竭 B6（周邊神經病變，需同時補充）；缺乏 → 鐵粒母細胞性貧血（ALA synthase 需要 PLP）、抽搐。"),
    dict(id="pyridoxal", name="Pyridoxal", label="Pyridoxal", aliases=["吡哆醛"],
         smiles="Cc1ncc(CO)c(C=O)c1O", ring={"rotate": 0}, tags=[WS],
         notes="C4′ 為醛基。pyridoxal kinase + ATP → PLP。"),
    dict(id="pyridoxamine", name="Pyridoxamine", label="Pyridoxamine", aliases=["吡哆胺"],
         smiles="Cc1ncc(CO)c(C[NH3+])c1O", ring={"rotate": 0}, tags=[WS], notes="C4′ 為 CH₂NH₃⁺。"),
    dict(id="plp", name="Pyridoxal phosphate (PLP)", label="Pyridoxal\nphosphate (PLP)", aliases=["PLP", "pyridoxal 5'-phosphate", "磷酸吡哆醛"],
         smiles="Cc1ncc(COP(=O)([O-])[O-])c(C=O)c1O", ring={"highlight": [{"smarts": "[CH1]=O", "color": "magenta"}]}, tags=[CO],
         notes="洋紅：C4′ 醛基，平時與酵素的 Lys 形成 Schiff base（internal aldimine），受質胺基酸進來後換成 external aldimine；吡啶環當「電子槽」。反應：轉胺（ALT、AST）、胺基酸脫羧（GABA、serotonin、dopamine、histamine）、ALA synthase（heme）、glycogen phosphorylase、cystathionine β-synthase。"),
    dict(id="pmp", name="Pyridoxamine phosphate (PMP)", label="Pyridoxamine\nphosphate (PMP)", aliases=["PMP", "pyridoxamine 5'-phosphate"],
         smiles="Cc1ncc(COP(=O)([O-])[O-])c(C[NH3+])c1O", ring={"highlight": [{"smarts": "[CH2:1][NH3+:2]", "color": "magenta"}]}, tags=[CO],
         notes="轉胺反應的中間型：PLP 從胺基酸拿走 –NH₂ 變成 PMP（胺基酸 → α-酮酸），再把 –NH₂ 交給另一個 α-酮酸而回到 PLP（乒乓機制）。"),
    # ---- B7
    dict(id="biotin", name="Biotin (vitamin B7)", label="Biotin (B_7)", aliases=["vitamin B7", "vitamin H", "生物素"],
         smiles="[O-]C(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12",
         ring={"highlight": [{"smarts": "[NH:1]([CH][CH2]S)C(=O)", "color": "magenta"}]}, tags=[WS, CO],
         notes="羧基以醯胺鍵接在酵素的 Lys 上（biocytin，holocarboxylase synthetase）。洋紅：N1′，ATP + HCO₃⁻ 先在這裡形成 carboxybiotin，再把 CO₂ 轉給受質。酵素：pyruvate carboxylase、acetyl-CoA carboxylase、propionyl-CoA carboxylase（ABC：ATP、Biotin、CO₂）。生蛋白的 avidin 會結合 biotin。"),
    # ---- lipoate
    dict(id="lipoate", name="Lipoate", label="Lipoate", aliases=["lipoic acid", "α-lipoic acid", "硫辛酸"],
         smiles="[O-]C(=O)CCCC[C@@H]1CCSS1", ring={"highlight": [{"smarts": "SS", "color": "magenta"}]}, tags=[CO],
         notes="不是維生素（人體可合成）。以醯胺鍵接在 E2 的 Lys 上（lipoamide，長擺臂）。洋紅：雙硫鍵，接受醯基並被還原成二硫醇，再由 E3（FAD）氧化回來。PDH、α-KGDH、BCKDH 的 E2。砷（arsenite）會與還原態的二硫醇結合。"),
    # ---- B9
    dict(id="folate", name="Folate (vitamin B9)", label="Folate (B_9)", aliases=["folic acid", "vitamin B9", "pteroylglutamate", "葉酸"],
         smiles="Nc1nc2ncc(CN" + PABA_GLU + ")nc2c(=O)[nH]1", ring={"rotate": 0}, tags=[WS],
         notes="pteridine + p-aminobenzoate（PABA）+ Glu。人類不能合成 PABA → 必須攝取（細菌自己合成：sulfonamide 抑制）。folate → DHF → THF 兩步都由 DHFR（NADPH）還原。缺乏：巨球性貧血（無神經症狀，與 B12 區分）、懷孕早期缺乏 → 神經管缺陷；homocysteine ↑，MMA 正常。"),
    # ---- B12 (schematic)
    dict(id="cobalamin", name="Cobalamin (vitamin B12)", label="Cobalamin (B_{12})",
         aliases=["vitamin B12", "cobalamin", "methylcobalamin", "adenosylcobalamin", "維生素B12", "鈷胺素"],
         smiles="[Co]",
         free={"schematic": True,
               "atoms": {"co": ["Co", 0, 0],
                         "na": ["N", -1.9, 0.25], "nb": ["N", -0.75, -0.6], "nc": ["N", 1.9, -0.25], "nd": ["N", 0.75, 0.6],
                         "p1": ["", -3.0, 0.7], "p2": ["", -1.2, -1.0], "p3": ["", 3.0, -0.7], "p4": ["", 1.2, 1.0],
                         "r": ["R", 0, -2.0], "dmb": ["DMB", 0, 2.0]},
               "bonds": [["p1", "p2"], ["p2", "p3"], ["p3", "p4"], ["p4", "p1"],
                         ["co", "na"], ["co", "nb"], ["co", "nc"], ["co", "nd"], ["co", "r"], ["co", "dmb"]],
               "notes": [["Corrin ring", 4.1, 0, "gray"],
                         ["R = CH_3 (methylcobalamin)\n      5'-deoxyadenosyl (adenosylcobalamin)\n      CN (cyanocobalamin), OH", 0, -3.3, "magenta"],
                         ["5,6-dimethylbenzimidazole", 0, 2.7, "gray"]]},
         decorate={"colors": {"r": "magenta"}}, tags=[WS, CO],
         notes="示意圖（只驗證 SMILES 格式，不檢查原子數）：corrin 環中央一個 Co，下方配位 DMB，上方配位 R（洋紅）決定功能。methylcobalamin：methionine synthase（Hcy → Met，同時 N⁵-CH₃-THF → THF；缺乏 → methyl-folate trap）。adenosylcobalamin：methylmalonyl-CoA mutase（→ succinyl-CoA；缺乏 → MMA ↑、神經症狀）。吸收需胃壁細胞的 intrinsic factor、迴腸末端。缺乏：巨球性貧血 + 亞急性聯合變性；Hcy ↑ 且 MMA ↑。"),
    # ---- methyl donor
    dict(id="sam", name="S-Adenosylmethionine (SAM)", label="S-Adenosyl-\nmethionine (SAM)", aliases=["SAM", "AdoMet", "S-腺苷甲硫胺酸"],
         smiles=nucleotide_sugar("Nc8ncnc9c8ncn9", "").replace("CO)", "C[S+](C)CC[C@H]([NH3+])C(=O)[O-])"),
         ring={"highlight": [{"smarts": "[CH3][S+]", "color": "magenta"}]}, tags=[CO],
         notes="Met + ATP → SAM（三個磷酸全部離開）。洋紅：鋶離子（S⁺）上的 CH₃，是體內主要的甲基供體 → SAH → homocysteine。例：norepinephrine → epinephrine（PNMT）、guanidinoacetate → creatine、DNA／組蛋白甲基化、PC 合成。"),
    # ---- C
    dict(id="ascorbate", name="Ascorbic acid (vitamin C)", label="Ascorbic acid (C)", aliases=["vitamin C", "ascorbate", "維生素C", "抗壞血酸"],
         smiles="OC[C@H](O)[C@H]1OC(=O)C(O)=C1O", ring={"highlight": [{"smarts": "[OH:1]C=C[OH:2]", "color": "magenta"}]}, tags=[WS],
         notes="洋紅：烯二醇（enediol），失去 2 H（2 e⁻ + 2 H⁺）→ dehydroascorbate，所以是還原劑／抗氧化劑。功能：prolyl／lysyl hydroxylase 的輔因子（膠原蛋白，維持 Fe²⁺）、dopamine β-hydroxylase、促進非血基質鐵（Fe³⁺ → Fe²⁺）吸收。缺乏：壞血病（牙齦出血、毛囊周圍出血、傷口不癒）。人類缺 L-gulonolactone oxidase，不能自行合成。"),
    dict(id="dehydroascorbate", name="Dehydroascorbate", label="Dehydroascorbate", aliases=["dehydroascorbic acid", "去氫抗壞血酸"],
         smiles="OC[C@H](O)[C@H]1OC(=O)C(=O)C1=O", ring={"rotate": 0}, tags=[WS],
         notes="ascorbate 的氧化態（C2、C3 變成酮）。可被 glutathione 還原回 ascorbate；以 GLUT 進入細胞。"),
    # ---- A
    dict(id="beta_carotene", name="β-Carotene", label="β-Carotene", aliases=["beta-carotene", "provitamin A", "β-胡蘿蔔素"],
         smiles="CC1=C(C(C)(C)CCC1)/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C2=C(C)CCCC2(C)C",
         ring={"methyls": False, "rotate": 0}, tags=[FS],
         notes="40 碳（類胡蘿蔔素，8 個異戊二烯單位）。腸道的 β-carotene 15,15′-dioxygenase 從正中間的雙鍵切開 → 2 分子 retinal。"),
    dict(id="retinol", name="Retinol (vitamin A)", label="Retinol (A)", aliases=["vitamin A", "維生素A", "視黃醇"],
         smiles="CC1=C(C(C)(C)CCC1)/C=C/C(C)=C/C=C/C(C)=C/CO", ring={"methyls": False}, tags=[FS],
         notes="20 碳，β-ionone 環＋全反式多烯側鏈，末端 CH₂OH。以 retinyl ester 存於肝 stellate（Ito）細胞，血中由 retinol-binding protein 運送。缺乏：夜盲、乾眼症（xerophthalmia）、Bitot spots、免疫力下降；過量：顱內壓上升、肝毒性、致畸胎（isotretinoin）。"),
    dict(id="retinal", name="all-trans-Retinal", label="all-trans-Retinal", aliases=["retinaldehyde", "視黃醛"],
         smiles="CC1=C(C(C)(C)CCC1)/C=C/C(C)=C/C=C/C(C)=C/C=O", ring={"methyls": False}, tags=[FS],
         notes="retinol 的醛（可逆）。光照使 rhodopsin 中的 11-cis-retinal 異構化成 all-trans-retinal → 視覺訊號；再於 RPE 異構化回 11-cis。"),
    dict(id="cis_retinal", name="11-cis-Retinal", label="11-cis-Retinal", aliases=["11-cis-retinaldehyde"],
         smiles="CC1=C(C(C)(C)CCC1)/C=C/C(C)=C/C=C\\C(C)=C\\C=O",
         ring={"methyls": False, "highlight": [{"atoms": [14, 15], "color": "magenta"}]}, tags=[FS],
         notes="洋紅：C11=C12 的 cis 雙鍵。醛基與 opsin 的 Lys 形成 Schiff base → rhodopsin；吸收光子後變成 all-trans（唯一的光化學步驟）。"),
    dict(id="retinoic_acid", name="Retinoic acid", label="Retinoic acid", aliases=["tretinoin", "all-trans-retinoic acid", "ATRA", "視黃酸"],
         smiles="CC1=C(C(C)(C)CCC1)/C=C/C(C)=C/C=C/C(C)=C/C(=O)O", ring={"methyls": False}, tags=[FS],
         notes="retinal 不可逆氧化而成（所以 retinoic acid 不能恢復視覺功能）。作用於核受體 RAR／RXR，調控上皮分化與胚胎發育。ATRA 治療 APL（t(15;17)，PML-RARA）。"),
    # ---- D
    dict(id="dehydrocholesterol", name="7-Dehydrocholesterol", label="7-Dehydrocholesterol",
         aliases=["7-DHC", "provitamin D3", "7-去氫膽固醇"],
         smiles="CC(C)CCC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2=CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
         ring={"orient": "steroid", "letters": True, "methyls": False,
               "highlight": [{"smarts": "[C;R]=[C;R][C;R]=[C;R]", "color": "magenta"}]}, tags=[FS, "lipid", "steroid"],
         notes="膽固醇合成的最後中間物（B 環多一個 C7=C8 雙鍵）。洋紅：B 環的共軛雙烯，皮膚中吸收 UVB 後 B 環在 C9–C10 斷開 → previtamin D₃ → （熱異構化）vitamin D₃。"),
    dict(id="cholecalciferol", name="Cholecalciferol (vitamin D3)", label="Cholecalciferol (D_3)",
         aliases=["vitamin D3", "維生素D3", "膽鈣化醇"],
         smiles=(SECO % ("[H]", "CCC3=C")).replace("C(C)(C)[H]", "C(C)C"), ring={"methyls": False}, tags=[FS],
         notes="secosteroid（B 環已斷開）。皮膚由 7-DHC 經 UVB 生成，或由食物（D₃ 動物、D₂ ergocalciferol 植物／菇類）取得。本身沒有活性：肝 25-hydroxylase → 25-OH-D₃ → 腎 1α-hydroxylase → calcitriol。"),
    dict(id="calcidiol", name="25-Hydroxycholecalciferol (calcidiol)", label="25-OH-D_3 (calcidiol)",
         aliases=["calcidiol", "25-hydroxyvitamin D3", "25(OH)D", "骨化二醇"],
         smiles=SECO % ("O", "CCC3=C"), ring={"methyls": False, "highlight": [{"smarts": "[CX4](C)(C)(C)[OH:1]", "color": "magenta"}]}, tags=[FS],
         notes="肝的 25-hydroxylase（CYP2R1）。洋紅：新加的 25-OH。血中主要的循環形式、半衰期長 → 臨床用來評估 vitamin D 狀態。"),
    dict(id="calcitriol", name="1,25-Dihydroxycholecalciferol (calcitriol)", label="1,25-(OH)_2-D_3 (calcitriol)",
         aliases=["calcitriol", "1,25-dihydroxyvitamin D3", "活性維生素D", "骨化三醇"],
         smiles=SECO % ("O", "C[C@@H](C3=C)O"),
         ring={"methyls": False, "highlight": [{"smarts": "C(=[CH2])[CH][OH:1]", "color": "magenta"},
                                               {"smarts": "[CX4](C)(C)(C)[OH:1]", "color": "magenta"}]}, tags=[FS],
         notes="活性形式。腎近曲小管的 1α-hydroxylase（CYP27B1）：PTH ↑、低磷 → 活化；calcitriol 本身與 FGF23 → 抑制。洋紅：1α-OH 與 25-OH。作用於核受體 VDR：腸道吸收 Ca²⁺ 與磷 ↑。缺乏：佝僂病（兒童）、骨軟化症（成人）；肉芽腫（sarcoidosis）的巨噬細胞也有 1α-hydroxylase → 高血鈣。"),
    # ---- E
    dict(id="tocopherol", name="α-Tocopherol (vitamin E)", label="α-Tocopherol (E)", aliases=["vitamin E", "alpha-tocopherol", "維生素E", "生育醇"],
         smiles="Cc1c(C)c2O[C@](C)(CCC[C@H](C)CCC[C@H](C)CCCC(C)C)CCc2c(C)c1O",
         ring={"methyls": False, "highlight": [{"smarts": "c[OH:1]", "color": "magenta"}]}, tags=[FS],
         notes="chromanol 環＋飽和的 phytyl 尾巴（嵌在膜中）。洋紅：酚 OH，把 H· 給脂質過氧自由基 → 中斷脂質過氧化連鎖反應；生成的 tocopheryl 自由基由 vitamin C 還原回來。缺乏（脂肪吸收不良）：溶血性貧血、棘狀紅血球、脊髓小腦退化（類似 B12 缺乏但無巨球性貧血）。"),
    # ---- K
    dict(id="phylloquinone", name="Phylloquinone (vitamin K1)", label="Phylloquinone (K_1)", aliases=["vitamin K1", "vitamin K", "維生素K1", "葉綠醌"],
         smiles="CC1=C(" + PHYTYL + ")C(=O)c2ccccc2C1=O", ring={"methyls": False}, tags=[FS],
         notes="2-methyl-1,4-naphthoquinone（醌）＋ phytyl 側鏈。植物（綠葉）來源；K₂（menaquinone）由腸道細菌合成 → 新生兒無菌腸道，出生時常規注射 vitamin K。"),
    dict(id="phylloquinol", name="Vitamin K hydroquinone (KH2)", label="Vitamin K\nhydroquinone (KH_2)", aliases=["KH2", "phylloquinol", "reduced vitamin K"],
         smiles="Cc1c(" + PHYTYL + ")c(O)c2ccccc2c1O", ring={"methyls": False, "highlight": [{"smarts": "c[OH:1]", "color": "magenta"}]}, tags=[FS],
         notes="還原態（活性形式）。洋紅：兩個酚 OH。γ-glutamyl carboxylase 用它與 O₂、CO₂ 把凝血因子 II、VII、IX、X 與 protein C、S 的 Glu 變成 γ-carboxyglutamate（Gla，可螯合 Ca²⁺），KH₂ 同時被氧化成 epoxide。"),
    dict(id="vitamin_k_epoxide", name="Vitamin K 2,3-epoxide", label="Vitamin K\n2,3-epoxide", aliases=["vitamin K epoxide", "KO"],
         smiles="CC12OC1(" + PHYTYL + ")C(=O)c1ccccc1C2=O",
         ring={"methyls": False, "highlight": [{"smarts": "[#6;r3:1]1[O;r3:2][#6;r3:3]1", "color": "magenta"}]}, tags=[FS],
         notes="洋紅：C2–C3 環氧（epoxide）。由 VKOR（vitamin K epoxide reductase）還原回醌；warfarin 抑制 VKOR → 無法再生 KH₂ → 凝血因子 γ-羧化不足（PT/INR ↑）。解藥：vitamin K、急性出血給 FFP／PCC。"),
    dict(id="menadione", name="Menadione (vitamin K3)", label="Menadione (K_3)", aliases=["vitamin K3", "menadione", "甲萘醌"],
         smiles="CC1=CC(=O)c2ccccc2C1=O", ring={"rotate": 0}, tags=[FS],
         notes="合成的 vitamin K（沒有側鏈），體內可接上側鏈。新生兒大量使用可能造成溶血與核黃疸，臨床已少用。"),
]

# THF family: full structures replace the old [7*] stand-in; reactions keep the short formula view
THF = [
    dict(id="thf", name="Tetrahydrofolate", label="Tetrahydrofolate (THF)", aliases=["THF", "H4folate", "四氫葉酸"],
         smiles="Nc1nc2NC[C@H](CN" + PABA_GLU + ")Nc2c(=O)[nH]1", formula="THF",
         ring={"rotate": 0, "highlight": [{"smarts": "[NH;R:1][CH;R]C[NH:2]c", "color": "magenta"}]},
         tags=["nucleotide", CO],
         notes="一碳單位載體（folate → DHF → THF，皆由 DHFR 還原）。洋紅：N⁵、N¹⁰，一碳單位（甲基、亞甲基、甲醯基…）接在這裡。天然為 (6S)。"),
    dict(id="dhf", name="Dihydrofolate", label="Dihydrofolate (DHF)", aliases=["DHF", "H2folate", "二氫葉酸"],
         smiles="Nc1nc2NCC(CN" + PABA_GLU + ")=Nc2c(=O)[nH]1", formula="DHF", ring={"rotate": 0},
         tags=["nucleotide", CO],
         notes="thymidylate synthase 唯一同時「給碳又給電子」的 THF 反應，所以產生 DHF，必須由 DHFR 再生成 THF（methotrexate、trimethoprim 抑制 DHFR）。"),
    dict(id="methylene_thf", name="N5,N10-Methylene-THF", label="N^5,N^{10}-Methylene-THF", aliases=["CH2-THF", "5,10-methylene-tetrahydrofolate"],
         smiles="Nc1nc2NC[C@@H]3CN(" + PABA_GLU + ")CN3c2c(=O)[nH]1", formula="N^5,N^{10}-CH_2-THF",
         ring={"coordgen": False, "highlight": [{"smarts": "[#7:1][CH2;r5:2][#7:3]", "color": "magenta"}]},
         tags=["nucleotide", CO],
         notes="洋紅：N⁵–CH₂–N¹⁰ 橋。thymidylate synthase 的甲基（＋還原力）來源；由 serine hydroxymethyltransferase（Ser → Gly）生成。"),
]

H = ["H+"]
REACTIONS = [
    dict(id="thiamine_pyrophosphokinase", title="Thiamine → TPP (thiamine pyrophosphokinase)",
         substrates=["thiamine"], products=["tpp"], cofactors_in=["ATP"], cofactors_out=["AMP"], extra_out=H,
         enzyme="thiamine\npyrophosphokinase", decorate={"tpp": {"highlight": [{"smarts": "P(=O)([O-])OP(=O)([O-])[O-]", "color": "magenta"}]}}, card={"type": "活化：轉移焦磷酸（ATP 的 β,γ-磷酸，AMP 離開）",
                                                    "exam_points": ["TPP 是 thiamine 的活性形式；酒精抑制 thiamine 吸收與活化"]}, tags=[CO]),
    dict(id="riboflavin_kinase", title="Riboflavin → FMN (riboflavin kinase)",
         substrates=["riboflavin"], products=["fmn"], cofactors_in=["ATP"], cofactors_out=["ADP"], extra_out=H,
         enzyme="riboflavin kinase\n(flavokinase)", decorate={"fmn": {"highlight": [{"smarts": "C[O:1][P:2](=[O:3])([O-:4])[O-:5]", "color": "magenta"}]}}, card={"type": "磷酸化（ribityl 的 5′-OH）"}, tags=[CO]),
    dict(id="fad_synthetase", title="FMN → FAD (FAD synthetase)",
         substrates=["fmn", "atp_n"], products=["fad"], cofactors_out=["PPi"], extra_in=H,
         enzyme="FAD synthetase\n(FMN adenylyltransferase)",
         card={"type": "腺苷酸化：ATP 的 AMP 部分接到 FMN 的磷酸上，放出 PPi",
               "exam_points": ["FAD = FMN + AMP（兩個核苷酸以焦磷酸相連 → dinucleotide）"]}, tags=[CO]),
    dict(id="pyridoxal_kinase", title="Pyridoxal → PLP (pyridoxal kinase)",
         substrates=["pyridoxal"], products=["plp"], cofactors_in=["ATP"], cofactors_out=["ADP"], extra_out=H,
         enzyme="pyridoxal kinase", decorate={"plp": {"highlight": [{"smarts": "C[O:1][P:2](=[O:3])([O-:4])[O-:5]", "color": "magenta"}]}}, card={"type": "磷酸化（C5′-CH₂OH）"}, tags=[CO]),
    dict(id="folate_reduction", title="Folate → DHF (dihydrofolate reductase)",
         substrates=["folate"], products=["dhf"], cofactors_in=["NADPH+H+"], cofactors_out=["NADP+"],
         enzyme="dihydrofolate\nreductase",
         card={"type": "還原（NADPH）：pteridine 環 C7=N8 加氫",
               "exam_points": ["DHFR 做兩次：folate → DHF → THF", "methotrexate 也擋這一步"]}, tags=[CO]),
    dict(id="vitamin_d_photolysis", title="7-Dehydrocholesterol → Vitamin D3 (UVB, skin)",
         substrates=["dehydrocholesterol"], products=["cholecalciferol"], enzyme="", conditions="UVB (skin)",
         card={"location": "皮膚（表皮）", "type": "非酵素：光解（B 環 C9–C10 斷開）→ previtamin D₃ → 熱異構化",
               "exam_points": ["日照不足、深膚色、高齡、防曬 → 合成減少"]}, tags=[FS]),
    dict(id="vitamin_d_25_hydroxylase", title="Vitamin D3 → Calcidiol (25-hydroxylase, liver)",
         substrates=["cholecalciferol"], products=["calcidiol"], cofactors_in=["O2", "NADPH+H+"], cofactors_out=["NADP+", "H2O"],
         enzyme="25-hydroxylase\n(CYP2R1, liver)", decorate={"calcidiol": {"highlight": [{"smarts": "[CX4](C)(C)(C)[OH:1]", "color": "magenta"}]}},
         card={"location": "肝（微粒體）", "type": "羥化（cytochrome P450 單氧酶）", "regulation": "幾乎不受調控（受質多就做多）",
               "exam_points": ["25-OH-D₃ 是評估 vitamin D 存量的指標"]}, tags=[FS]),
    dict(id="vitamin_d_1a_hydroxylase", title="Calcidiol → Calcitriol (1α-hydroxylase, kidney)",
         substrates=["calcidiol"], products=["calcitriol"], cofactors_in=["O2", "NADPH+H+"], cofactors_out=["NADP+", "H2O"],
         enzyme="1α-hydroxylase\n(CYP27B1, kidney)", decorate={"calcitriol": {"highlight": [{"smarts": "C(=[CH2])[CH][OH:1]", "color": "magenta"}]}},
         card={"location": "腎近曲小管（粒線體）", "type": "羥化（cytochrome P450 單氧酶）",
               "regulation": "PTH、低血磷 → 活化；calcitriol、FGF23 → 抑制",
               "deficiency": ["慢性腎病 → calcitriol ↓ → 低血鈣、次發性副甲狀腺亢進（腎性骨病變）",
                              "sarcoidosis 等肉芽腫的巨噬細胞也有 1α-hydroxylase（不受 PTH 調控）→ 高血鈣"],
               "exam_points": ["限速、受調控的一步"]}, tags=[FS]),
    dict(id="beta_carotene_cleavage", title="β-Carotene → 2 Retinal (β-carotene dioxygenase)",
         substrates=["beta_carotene"], products=[{"id": "retinal", "n": 2}], cofactors_in=["O2"],
         enzyme="β-carotene\n15,15'-dioxygenase", products_layout="stack",
         card={"location": "小腸黏膜", "type": "氧化斷裂（從正中間的雙鍵切開）"}, tags=[FS]),
    dict(id="retinol_dehydrogenase", title="Retinol ⇌ Retinal (retinol dehydrogenase)",
         substrates=["retinol"], products=["retinal"], cofactors_in=["NAD+"], cofactors_out=["NADH+H+"], reversible=True,
         enzyme="retinol\ndehydrogenase", decorate={"retinol": {"highlight": [{"smarts": "[CH2:1][OH:2]", "color": "magenta"}]},
                                              "retinal": {"highlight": [{"smarts": "[CH:1]=[O:2]", "color": "magenta"}]}}, card={"type": "氧化還原（可逆）：醇 ⇌ 醛"}, tags=[FS]),
    dict(id="retinal_dehydrogenase", title="Retinal → Retinoic acid (retinal dehydrogenase)",
         substrates=["retinal"], products=["retinoic_acid"], cofactors_in=["NAD+", "H2O"], cofactors_out=["NADH+H+"],
         enzyme="retinal\ndehydrogenase", decorate={"retinal": {"highlight": [{"smarts": "[CH:1]=[O:2]", "color": "magenta"}]},
                                              "retinoic_acid": {"highlight": [{"smarts": "[C:1](=[O:2])[OH:3]", "color": "magenta"}]}},
         card={"type": "氧化（不可逆）：醛 → 羧酸", "exam_points": ["retinoic acid 無法回到 retinal → 不能支持視覺"]}, tags=[FS]),
]


PATHWAYS = [
    dict(id="vitamin_d_activation", title="Vitamin D activation（維生素 D 的活化）", layout="free", cell=1.0,
         nodes=[{"id": "dhc", "ref": "dehydrocholesterol", "pos": [0, 0]},
                {"id": "d3", "ref": "cholecalciferol", "pos": [1.45, 0]},
                {"id": "c25", "ref": "calcidiol", "pos": [1.45, 1.15]},
                {"id": "c1", "ref": "calcitriol", "pos": [0, 1.15]}],
         edges=[{"from": "dhc", "to": "d3", "step": "vitamin_d_photolysis", "label": "UVB (skin)"},
                {"from": "d3", "to": "c25", "step": "vitamin_d_25_hydroxylase", "label": "25-hydroxylase\n(liver)", "label_side": "right"},
                {"from": "c25", "to": "c1", "step": "vitamin_d_1a_hydroxylase", "label": "1α-hydroxylase\n(kidney)"}],
         tags=[FS], notes="皮膚（UVB）→ 肝（25-OH）→ 腎（1α-OH）。慢性腎病缺 1α-hydroxylase → 需直接給 calcitriol；肝病影響較小。"),
    dict(id="vitamin_a_forms", title="Vitamin A: carotene, retinol, retinal, retinoic acid（維生素 A 的各種形式）",
         layout="free", cell=0.62,
         nodes=[{"id": "car", "ref": "beta_carotene", "pos": [0.9, 0]},
                {"id": "ral", "ref": "retinal", "pos": [0, 1]},
                {"id": "rol", "ref": "retinol", "pos": [1.8, 1]},
                {"id": "ra", "ref": "retinoic_acid", "pos": [0, 2]},
                {"id": "cis", "ref": "cis_retinal", "pos": [1.8, 2]},
                {"id": "rho", "text": "Rhodopsin\n(+ opsin)", "kind": "concept", "pos": [1.8, 2.75]},
                {"id": "rar", "text": "RAR / RXR\n(gene expression)", "kind": "concept", "pos": [0, 2.75]}],
         edges=[{"from": "car", "to": "ral", "step": "beta_carotene_cleavage", "label": "dioxygenase", "side_text": "O_2"},
                {"from": "ral", "to": "rol", "step": "retinol_dehydrogenase", "reversible": True,
                 "label": "retinol DH"},
                {"from": "ral", "to": "ra", "step": "retinal_dehydrogenase", "label": "retinal DH", "label_side": "right"},
                {"from": "ral", "to": "cis", "label": "isomerase (RPE)", "reversible": True},
                {"from": "cis", "to": "rho", "head": False},
                {"from": "ra", "to": "rar", "head": False}],
         tags=[FS], notes="retinol ⇌ retinal 可逆；retinal → retinoic acid 不可逆。視覺用 retinal（11-cis），細胞分化用 retinoic acid。"),
    dict(id="vitamin_k_cycle", title="Vitamin K cycle（維生素 K 循環）", layout="free", cell=0.85,
         nodes=[{"id": "k", "ref": "phylloquinone", "pos": [0, 0]},
                {"id": "kh2", "ref": "phylloquinol", "pos": [1.25, 1]},
                {"id": "ko", "ref": "vitamin_k_epoxide", "pos": [0, 2]},
                {"id": "war", "text": "Warfarin", "kind": "drug", "pos": [-0.75, 1]}],
         edges=[{"from": "k", "to": "kh2", "label": "VKOR /\nquinone reductase", "side_text": "NAD(P)H"},
                {"from": "kh2", "to": "ko", "label": "γ-glutamyl\ncarboxylase", "side_text": "Glu → Gla\n(CO_2, O_2)"},
                {"from": "ko", "to": "k", "label": "VKOR"},
                {"from": "war", "to": "edge:ko->k", "inhibit": True}],
         tags=[FS], notes="KH₂ 是 γ-carboxylase 真正的輔因子（因子 II、VII、IX、X、protein C、S 的 Glu → Gla）。warfarin 擋 VKOR（epoxide → K），使 KH₂ 無法再生；大量 vitamin K 可經另一個不受 warfarin 抑制的還原酶（NAD(P)H 依賴）恢復凝血。"),
    dict(id="vitamin_coenzymes", title="Water-soluble vitamins → coenzymes（水溶性維生素與其輔酶）", layout="free",
         cell_size=6.6, text=False,
         nodes=[
             {"id": "h1", "text": "Vitamin", "kind": "heading", "pos": [0, -0.7]},
             {"id": "h2", "text": "Coenzyme", "kind": "heading", "pos": [2.3, -0.7]},
             {"id": "h3", "text": "Reaction type / enzymes", "kind": "heading", "pos": [4.7, -0.7]},
             {"id": "b1", "ref": "thiamine", "view": "name", "pos": [0, 0]},
             {"id": "tpp", "ref": "tpp", "view": "name", "pos": [2.3, 0]},
             {"id": "f1", "text": "Oxidative decarboxylation (PDH, α-KGDH, BCKDH),\ntransketolase", "kind": "plain", "pos": [4.9, 0]},
             {"id": "b2", "ref": "riboflavin", "view": "name", "pos": [0, 1]},
             {"id": "fmn", "ref": "fmn", "view": "name", "pos": [1.3, 1]},
             {"id": "fad", "ref": "fad", "view": "name", "pos": [2.3, 1]},
             {"id": "f2", "text": "Redox, 1 or 2 e^- (Complex I, II,\nacyl-CoA DH, E3 of PDH)", "kind": "plain", "pos": [4.9, 1]},
             {"id": "b3", "ref": "nicotinate", "view": "name", "pos": [0, 2]},
             {"id": "nad", "ref": "nad", "view": "name", "pos": [2.3, 2]},
             {"id": "nadp", "ref": "nadp", "view": "name", "pos": [2.3, 2.6]},
             {"id": "f3", "text": "Redox, hydride (2 e^-): NAD^+ catabolic,\nNADPH anabolic", "kind": "plain", "pos": [4.9, 2.3]},
             {"id": "b5", "ref": "pantothenate", "view": "name", "pos": [0, 3.4]},
             {"id": "coa", "ref": "coa", "view": "name", "pos": [2.3, 3.4]},
             {"id": "f5", "text": "Acyl carrier (thioester): acetyl-CoA,\nfatty acid synthase (ACP)", "kind": "plain", "pos": [4.9, 3.4]},
             {"id": "b6", "ref": "pyridoxal", "view": "name", "pos": [0, 4.4]},
             {"id": "plp", "ref": "plp", "view": "name", "pos": [2.3, 4.4]},
             {"id": "f6", "text": "Amino acids: transamination, decarboxylation;\nALA synthase, glycogen phosphorylase", "kind": "plain", "pos": [4.9, 4.4]},
             {"id": "b7", "ref": "biotin", "view": "name", "pos": [0, 5.4]},
             {"id": "bl", "text": "Biotin–Lys\n(carboxylases)", "kind": "name", "pos": [2.3, 5.4]},
             {"id": "f7", "text": "Carboxylation (+CO_2): pyruvate carboxylase,\nacetyl-CoA C, propionyl-CoA C", "kind": "plain", "pos": [4.9, 5.4]},
             {"id": "b9", "ref": "folate", "view": "name", "pos": [0, 6.4]},
                          {"id": "thf", "ref": "thf", "view": "name", "pos": [2.3, 6.4]},
             {"id": "f9", "text": "One-carbon transfer: purines, dTMP,\nSer → Gly, Met", "kind": "plain", "pos": [4.9, 6.4]},
             {"id": "b12", "ref": "cobalamin", "view": "name", "pos": [0, 7.4]},
             {"id": "mecbl", "text": "Methyl-Cbl\nAdenosyl-Cbl", "kind": "name", "pos": [2.3, 7.4]},
             {"id": "f12", "text": "Methionine synthase (Hcy → Met),\nmethylmalonyl-CoA mutase", "kind": "plain", "pos": [4.9, 7.4]},
             {"id": "c", "ref": "ascorbate", "view": "name", "pos": [0, 8.4]},
             {"id": "fc", "text": "Hydroxylation (Pro, Lys: collagen),\nantioxidant, Fe absorption", "kind": "plain", "pos": [4.9, 8.4]},
         ],
         edges=[{"from": "b1", "to": "tpp", "step": "thiamine_pyrophosphokinase", "label": "ATP → AMP"},
                {"from": "b2", "to": "fmn", "step": "riboflavin_kinase", "label": "ATP"},
                {"from": "fmn", "to": "fad", "step": "fad_synthetase", "label": "ATP"},
                {"from": "b3", "to": "nad", "label": "Preiss–Handler"},
                {"from": "nad", "to": "nadp", "label": "NAD kinase", "label_side": "right"},
                {"from": "b5", "to": "coa", "label": "+ Cys, ATP"},
                {"from": "b6", "to": "plp", "step": "pyridoxal_kinase", "label": "ATP"},
                {"from": "b7", "to": "bl", "label": "holocarboxylase\nsynthetase"},
                {"from": "b9", "to": "thf", "label": "DHFR × 2 (via DHF)"},
                {"from": "b12", "to": "mecbl"}],
         tags=[WS, CO], notes="B1 TPP、B2 FAD/FMN、B3 NAD(P)、B5 CoA、B6 PLP、B7 biotin、B9 THF、B12 cobalamin。口訣：PDH 需要 B1、B2、B3、B5、lipoate（「Tender Loving Care For Nancy」：TPP、lipoate、CoA、FAD、NAD⁺）。"),
]


def dump(path, d):
    path.write_text(yaml.safe_dump(d, allow_unicode=True, sort_keys=False, width=200), encoding="utf-8")


def mol_doc(m):
    d = {"id": m["id"], "name": m["name"]}
    if m.get("label"):
        d["label"] = m["label"]
    d["aliases"] = m.get("aliases", [])
    d["smiles"] = m["smiles"]
    ring = dict(m["ring"]) if m.get("ring") is not None else None
    deco = dict(m.get("decorate") or {})
    if ring and ring.get("highlight"):            # reactive part: standalone figure and stickers only
        deco["highlight"] = ring.pop("highlight")
        ring = ring or {"coordgen": True}           # an empty ring spec would drop the view
    for k, v in (("ring", ring or None), ("free", m.get("free")), ("formula", m.get("formula")), ("decorate", deco or None)):
        if v is not None:
            d[k] = v
    d["tags"] = ["vitamin"] + m["tags"] if "vitamin" not in m["tags"] else m["tags"]
    d["notes_zh"] = m["notes"]
    d["status"] = "pending"
    return d


def main(data):
    data = Path(data)
    for m in MOLS + THF:
        dump(data / "molecules" / f"{m['id']}.yaml", mol_doc(m))
    for r in REACTIONS:
        d = {"id": r["id"], "title": r["title"], "substrates": r["substrates"], "products": r["products"]}
        for k in ("products_layout", "cofactors_in", "cofactors_out", "extra_in", "extra_out"):
            if r.get(k):
                d[k] = r[k]
        d["enzyme"] = r.get("enzyme", "")
        deco = COLORS.get(r["id"], r.get("decorate"))       # tools/reaction_colors.py wins
        if deco:
            d["decorate"] = deco
        if r.get("conditions"):
            d["conditions"] = r["conditions"]
        d["reversible"] = r.get("reversible", False)
        d["card_zh"] = r["card"]
        d["tags"] = ["vitamin"] + r["tags"]
        d["status"] = "pending"
        dump(data / "reactions" / f"{r['id']}.yaml", d)
    for p in PATHWAYS:
        d = {k: v for k, v in p.items() if k not in ("tags", "notes")}
        d["tags"] = ["vitamin"] + p["tags"]
        d["notes_zh"] = p["notes"]
        d["status"] = "pending"
        dump(data / "pathways" / f"{p['id']}.yaml", d)
    print(f"{len(MOLS) + len(THF)} molecules, {len(REACTIONS)} reactions, {len(PATHWAYS)} pathways")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data")
