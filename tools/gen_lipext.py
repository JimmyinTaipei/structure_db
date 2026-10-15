"""Lipid metabolism, extended (round 12): eicosanoids, steroid hormones, bile acids,
sphingolipid synthesis and the lysosomal sphingolipidoses.

    .venv/bin/python tools/gen_lipext.py data

SMILES come from PubChem (isomeric) and are converted to the physiological ionisation used in
this database (COO⁻, SO₃⁻, R–NH₃⁺) by `physio()`. Ring drawings do not show wedges, so the
stereocentres in the SMILES are kept for the formula / CIP information only; C=C geometry is
drawn. Glycosphingolipids are SNFG cartoons (`free`, schematic) on top of a real SMILES built
from ceramide + sugars (no stereo), so the hydrolase reactions still balance.
Re-running overwrites the files listed here and appends cofactors (GSH, GSSG …).
"""
import sys
from pathlib import Path

import yaml
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

sys.path.insert(0, str(Path(__file__).parent))
from reaction_colors import COLORS  # noqa: E402

CM = "lipid_metabolism"
KIDS = ["eicosanoid_synthesis", "steroid_hormone", "bile_acid", "sphingolipid_metabolism"]
STEROID = {"letters": True, "orient": "steroid", "methyls": False}


def physio(smiles):
    """Neutral PubChem SMILES → pH 7.4 form: COOH/SO3H deprotonated, aliphatic amines protonated."""
    m = Chem.RWMol(Chem.MolFromSmiles(smiles))
    for a in m.GetAtoms():
        if a.GetSymbol() == "O" and a.GetTotalNumHs() == 1 and a.GetDegree() == 1:
            c = a.GetNeighbors()[0]
            acid = (c.GetSymbol() == "C" and any(b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(c).GetSymbol() == "O"
                                                 for b in c.GetBonds())) or c.GetSymbol() == "S"
            if acid:
                a.SetFormalCharge(-1)
                a.SetNumExplicitHs(0)
                a.SetNoImplicit(True)
        if a.GetSymbol() == "N" and not a.GetIsAromatic() and a.GetFormalCharge() == 0 and a.GetTotalNumHs() == 2:
            c = a.GetNeighbors()[0]
            amide = any(b.GetBondType() == Chem.BondType.DOUBLE for b in c.GetBonds())
            if not amide:
                a.SetFormalCharge(1)
                a.SetNumExplicitHs(3)
    mol = m.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


PC = {  # PubChem isomeric SMILES (neutral)
    "pgg2": "CCCCC[C@@H](/C=C/[C@H]1[C@H]2C[C@@H]([C@@H]1C/C=C\\CCCC(=O)O)OO2)OO",
    "pgh2": "CCCCC[C@@H](/C=C/[C@H]1[C@H]2C[C@@H]([C@@H]1C/C=C\\CCCC(=O)O)OO2)O",
    "pgf2a": "CCCCC[C@@H](/C=C/[C@H]1[C@@H](C[C@@H]([C@@H]1C/C=C\\CCCC(=O)O)O)O)O",
    "pgi2": "CCCCC[C@@H](/C=C/[C@H]1[C@@H](C[C@H]2[C@@H]1C/C(=C/CCCC(=O)O)/O2)O)O",
    "txa2": "CCCCC[C@@H](/C=C/[C@@H]1[C@H]([C@@H]2C[C@@H](O2)O1)C/C=C\\CCCC(=O)O)O",
    "lta4": "CCCCC/C=C\\C/C=C\\C=C\\C=C\\[C@H]1[C@@H](O1)CCCC(=O)O",
    "ltb4": "CCCCC/C=C\\C[C@H](/C=C/C=C/C=C\\[C@H](CCCC(=O)O)O)O",
    "ltc4": "CCCCC/C=C\\C/C=C\\C=C\\C=C\\[C@H]([C@H](CCCC(=O)O)O)SC[C@@H](C(=O)NCC(=O)O)NC(=O)CC[C@@H](C(=O)O)N",
    "hpete5": "CCCCC/C=C\\C/C=C\\C/C=C\\C=C\\[C@H](CCCC(=O)O)OO",
    "pregnenolone": "CC(=O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
    "oh_pregnenolone": "CC(=O)[C@]1(CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)O",
    "oh_progesterone": "CC(=O)[C@]1(CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4=CC(=O)CC[C@]34C)C)O",
    "doc": "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2C(=O)CO)CCC4=CC(=O)CC[C@]34C",
    "corticosterone": "C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2[C@H](C[C@]4([C@H]3CC[C@@H]4C(=O)CO)C)O",
    "deoxycortisol": "C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2CC[C@]4([C@H]3CC[C@@]4(C(=O)CO)O)C",
    "dhea": "C[C@]12CC[C@H]3[C@H]([C@@H]1CCC2=O)CC=C4[C@@]3(CC[C@@H](C4)O)C",
    "androstenedione": "C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2CC[C@]4([C@H]3CCC4=O)C",
    "dht": "C[C@]12CCC(=O)C[C@@H]1CC[C@@H]3[C@@H]2CC[C@]4([C@H]3CC[C@@H]4O)C",
    "estrone": "C[C@]12CC[C@H]3[C@H]([C@@H]1CCC2=O)CCC4=C3C=CC(=C4)O",
    "oh7_cholesterol": "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2[C@@H](C=C4[C@@]3(CC[C@@H](C4)O)C)O)C",
    "chenodeoxycholate": "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)C",
    "deoxycholate": "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2CC[C@H]4[C@@]3(CC[C@H](C4)O)C)O)C",
    "lithocholate": "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@H]4[C@@]3(CC[C@H](C4)O)C)C",
    "glycocholate": "C[C@H](CCC(=O)NCC(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)O)C",
    "taurocholate": "C[C@H](CCC(=O)NCCS(=O)(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)O)C",
    "sphinganine": "CCCCCCCCCCCCCCC[C@H]([C@H](CO)N)O",
    "ketosphinganine": "CCCCCCCCCCCCCCCC(=O)[C@H](CO)N",
    "sphingosine": "CCCCCCCCCCCCC/C=C/[C@H]([C@H](CO)N)O",
    "isocaproaldehyde": "CC(C)CCC=O",
    "gsh": "C(CC(=O)N[C@@H](CS)C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
    "gssg": "C(CC(=O)N[C@@H](CSSC[C@@H](C(=O)NCC(=O)O)NC(=O)CC[C@@H](C(=O)O)N)C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
}

# ---- glycosphingolipids: ceramide (d18:1/18:0, as in data/molecules/ceramide.yaml) + sugars, no stereo
CER = "C(O{})C(NC(=O)CCCCCCCCCCCCCCCCC)C(O)/C=C/CCCCCCCCCCCCC"
NEUAC = "C9(C(=O)[O-])CC(O)C(NC(C)=O)C(O9)C(O)C(O)CO"           # α2→3 linked through C2


def glc(o4=""):
    return f"C1OC(CO)C(O{o4})C(O)C1O" if o4 else "C1OC(CO)C(O)C(O)C1O"


def gal(o4="", o3="", ring="2"):
    a = f"O{o4}" if o4 else "O"
    b = f"O{o3}" if o3 else "O"
    return f"C{ring}OC(CO)C({a})C({b})C{ring}O"


GALNAC_GM2 = "C4OC(CO)C(O)C(O)C4NC(C)=O"
GALNAC_GM1 = "C4OC(CO)C(O)C(OC5OC(CO)C(O)C(O)C5O)C4NC(C)=O"
GSL = {
    "lactosylceramide": CER.format(glc(gal())),
    "globotriaosylceramide": CER.format(glc(gal(o4="C3OC(CO)C(O)C(O)C3O"))),
    "gm3": CER.format(glc(gal(o3=NEUAC))),
    "gm2": CER.format(glc(gal(o4=GALNAC_GM2, o3=NEUAC))),
    "gm1": CER.format(glc(gal(o4=GALNAC_GM1, o3=NEUAC))),
}


def snfg_cartoon(chain, branch=None, extra=None):
    """Ceramide on the left, sugars to the right (SNFG symbols); branch = NeuAc under the Gal."""
    atoms = {"cer": ["Ceramide", 0, 0, 7]}                  # text runs to the left of its last letter
    bonds = []
    prev = "cer"
    for k, s in enumerate(chain):
        aid = f"s{k}"
        atoms[aid] = [f"@{s}", 0.95 + 0.8 * k, 0]
        bonds.append([prev, aid])
        prev = aid
    if branch:
        k, s = branch
        atoms["b"] = [f"@{s}", 0.95 + 0.8 * k, 0.8]
        bonds.append([f"s{k}", "b"])
    if extra:
        k, s = extra
        atoms["e"] = [f"@{s}", 0.95 + 0.8 * (k + 1), 0]
        bonds.append([f"s{k}", "e"])
    return {"atoms": atoms, "bonds": bonds, "schematic": True}


SNFG_NOTE = "（SNFG 示意圖：藍圓 Glc、黃圓 Gal、黃方 GalNAc、紫菱 NeuAc）"


def molecules():
    E = ["eicosanoid_synthesis"]
    S = ["steroid_hormone"]
    B = ["bile_acid"]
    SP = ["sphingolipid_metabolism"]
    ring = {"rotate": 0}
    out = [
        # ---- eicosanoids
        dict(id="pgg2", name="Prostaglandin G2", label="PGG_2", aliases=["PGG₂"], smiles=physio(PC["pgg2"]), ring=ring, tags=E,
             notes="COX 的環氧酶（cyclooxygenase）活性：arachidonate + 2 O₂ → 環內過氧化物（9,11-endoperoxide）＋ 15-OOH。"),
        dict(id="pgh2", name="Prostaglandin H2", label="PGH_2", aliases=["PGH₂"], smiles=physio(PC["pgh2"]), ring=ring, tags=E,
             notes="COX 的過氧化酶活性把 15-OOH 還原成 15-OH。所有 prostanoid（PGE₂、PGF₂α、PGI₂、TXA₂）的共同前驅物，由各組織的 synthase 決定去向。"),
        dict(id="pgf2a", name="Prostaglandin F2α", label="PGF_2α", aliases=["PGF₂α", "dinoprost"], smiles=physio(PC["pgf2a"]), ring=ring, tags=E,
             notes="子宮收縮（分娩、月經痛）、支氣管收縮；carboprost、latanoprost（青光眼）為類似物。"),
        dict(id="pgi2", name="Prostacyclin (PGI2)", label="Prostacyclin (PGI_2)", aliases=["PGI₂", "epoprostenol"], smiles=physio(PC["pgi2"]), ring=ring, tags=E,
             notes="血管內皮產生：血管擴張、抑制血小板凝集（與 TXA₂ 相反）。epoprostenol 治療肺動脈高壓。"),
        dict(id="txa2", name="Thromboxane A2", label="Thromboxane A_2", aliases=["TXA₂"], smiles=physio(PC["txa2"]), ring=ring, tags=E,
             notes="血小板（COX-1）產生：血管收縮、促進血小板凝集。極不穩定（半衰期約 30 秒）→ TXB₂。低劑量 aspirin 不可逆抑制血小板 COX-1（血小板無核，無法再合成）→ 抗血栓。"),
        dict(id="hpete5", name="5-HPETE", aliases=["5-hydroperoxyeicosatetraenoate"], smiles=physio(PC["hpete5"]), ring=ring, tags=E,
             notes="5-lipoxygenase（需 FLAP）在 C5 加入 O₂。"),
        dict(id="lta4", name="Leukotriene A4", label="LTA_4", aliases=["LTA₄"], smiles=physio(PC["lta4"]), ring=ring, tags=E,
             notes="5-LOX 的第二個活性把 5-HPETE 脫水成環氧化物（5,6-epoxide），共軛三烯。"),
        dict(id="ltb4", name="Leukotriene B4", label="LTB_4", aliases=["LTB₄"], smiles=physio(PC["ltb4"]), ring=ring, tags=E,
             notes="LTA₄ hydrolase 產生：嗜中性球趨化（chemotaxis）。"),
        dict(id="ltc4", name="Leukotriene C4", label="LTC_4", aliases=["LTC₄", "cysteinyl leukotriene"],
             smiles="CCCCC/C=C\\C/C=C\\C=C\\C=C\\[C@H]([C@H](CCCC(=O)[O-])O)S[16*]", ring=ring, tags=E,
             notes="S–G：glutathione（γ-Glu-Cys-Gly）以 Cys 的 S 接在 C6。LTC₄ synthase 接上 glutathione → LTD₄ → LTE₄（cysteinyl LT，舊稱 SRS-A）：支氣管收縮、血管通透性 ↑。montelukast、zafirlukast 阻斷受體；zileuton 抑制 5-LOX。aspirin-exacerbated respiratory disease：COX 被抑制 → 更多 arachidonate 走 LOX。"),
        # ---- steroid hormones
        dict(id="pregnenolone", name="Pregnenolone", aliases=["孕烯醇酮"], smiles=physio(PC["pregnenolone"]), ring=STEROID, tags=S,
             notes="所有固醇類荷爾蒙的共同前驅物（21 碳）。CYP11A1（side-chain cleavage, desmolase）在粒線體切掉膽固醇 C22–C27；ACTH（腎上腺）、LH（性腺）調控，StAR 蛋白把膽固醇送進粒線體內膜是限速步驟。"),
        dict(id="oh_pregnenolone", name="17α-Hydroxypregnenolone", label="17α-Hydroxy-\npregnenolone", aliases=["17-OH pregnenolone"],
             smiles=physio(PC["oh_pregnenolone"]), ring=STEROID, tags=S, notes="CYP17A1 的 17α-hydroxylase 活性（束狀帶、網狀帶、性腺；球狀帶沒有）。"),
        dict(id="oh_progesterone", name="17α-Hydroxyprogesterone", label="17α-Hydroxy-\nprogesterone", aliases=["17-OHP"],
             smiles=physio(PC["oh_progesterone"]), ring=STEROID, tags=S,
             notes="21-hydroxylase 缺乏時累積 → 新生兒篩檢 CAH 的指標；轉向雄性素 → 女嬰外生殖器男性化。"),
        dict(id="doc", name="11-Deoxycorticosterone", label="11-Deoxy-\ncorticosterone", aliases=["DOC", "deoxycorticosterone"],
             smiles=physio(PC["doc"]), ring=STEROID, tags=S,
             notes="本身有礦物皮質素活性：11β-hydroxylase 缺乏時 DOC ↑ → 高血壓、低血鉀（同時男性化）。"),
        dict(id="corticosterone", name="Corticosterone", aliases=["皮質酮"], smiles=physio(PC["corticosterone"]), ring=STEROID, tags=S,
             notes="球狀帶由 aldosterone synthase（CYP11B2）再經 18-hydroxylation、18-oxidation → aldosterone。"),
        dict(id="deoxycortisol", name="11-Deoxycortisol", label="11-Deoxycortisol", aliases=["cortodoxone", "compound S"],
             smiles=physio(PC["deoxycortisol"]), ring=STEROID, tags=S, notes="CYP11B1（11β-hydroxylase）→ cortisol。metyrapone 抑制此步（ACTH 儲備測試）。"),
        dict(id="dhea", name="Dehydroepiandrosterone (DHEA)", label="DHEA", aliases=["dehydroepiandrosterone", "DHEA-S"],
             smiles=physio(PC["dhea"]), ring=STEROID, tags=S, notes="CYP17A1 的 17,20-lyase 活性切掉 C20–C21；網狀帶主要產物（多為 DHEA-S），腎上腺雄性素。"),
        dict(id="androstenedione", name="Androstenedione", aliases=["雄烯二酮"], smiles=physio(PC["androstenedione"]), ring=STEROID, tags=S,
             notes="3β-HSD 作用於 DHEA 或 17,20-lyase 作用於 17-OHP 而來；→ testosterone（17β-HSD）或 estrone（aromatase）。"),
        dict(id="dht", name="Dihydrotestosterone (DHT)", label="DHT", aliases=["5α-dihydrotestosterone"], smiles=physio(PC["dht"]), ring=STEROID, tags=S,
             notes="5α-reductase（NADPH）還原 C4=C5；比 testosterone 強的雄性素：外生殖器、前列腺、毛囊。finasteride、dutasteride 抑制（BPH、雄性禿）。5α-reductase 缺乏 → 出生時外生殖器似女性，青春期男性化。"),
        dict(id="estrone", name="Estrone", aliases=["E1", "雌酮"], smiles=physio(PC["estrone"]), ring=STEROID, tags=S,
             notes="aromatase 作用於 androstenedione；停經後主要的 estrogen（脂肪組織）。"),
        dict(id="isocaproaldehyde", name="Isocaproaldehyde", aliases=["4-methylpentanal"], smiles=physio(PC["isocaproaldehyde"]), tags=S,
             notes="CYP11A1 切下的膽固醇側鏈（C22–C27）。"),
        # ---- bile acids
        dict(id="oh7_cholesterol", name="7α-Hydroxycholesterol", label="7α-Hydroxy-\ncholesterol", aliases=["7α-OH-cholesterol"],
             smiles=physio(PC["oh7_cholesterol"]), ring=STEROID, tags=B,
             notes="CYP7A1（cholesterol 7α-hydroxylase）是膽酸合成的限速步驟：膽酸經 FXR 回饋抑制；cholestyramine 打斷腸肝循環 → CYP7A1 ↑ → 肝膽固醇 ↓ → LDL 受體 ↑。"),
        dict(id="chenodeoxycholate", name="Chenodeoxycholate", aliases=["chenodeoxycholic acid", "CDCA", "鵝去氧膽酸"],
             smiles=physio(PC["chenodeoxycholate"]), ring=STEROID, tags=B, notes="初級膽酸（3α,7α-OH）。"),
        dict(id="deoxycholate", name="Deoxycholate", aliases=["deoxycholic acid", "去氧膽酸"], smiles=physio(PC["deoxycholate"]), ring=STEROID, tags=B,
             notes="次級膽酸：腸道細菌去掉 cholate 的 7α-OH（3α,12α-OH）。"),
        dict(id="lithocholate", name="Lithocholate", aliases=["lithocholic acid", "石膽酸"], smiles=physio(PC["lithocholate"]), ring=STEROID, tags=B,
             notes="次級膽酸：chenodeoxycholate 去 7α-OH（只剩 3α-OH）；最不溶、具肝毒性，大多隨糞便排出。"),
        dict(id="glycocholate", name="Glycocholate", aliases=["glycocholic acid"], smiles=physio(PC["glycocholate"]), ring=STEROID, tags=B,
             notes="膽鹽：cholyl-CoA + glycine（BAAT）→ 醯胺鍵；pKa 降低 → 在小腸 pH 下完全解離、兩性 → 乳化脂質、形成 micelle。"),
        dict(id="taurocholate", name="Taurocholate", aliases=["taurocholic acid"], smiles=physio(PC["taurocholate"]), ring=STEROID, tags=B,
             notes="與 taurine 結合的膽鹽（磺酸根，pKa 更低）。"),
        # ---- sphingolipids
        dict(id="ketosphinganine", name="3-Ketosphinganine", label="3-Keto-\nsphinganine", aliases=["3-dehydrosphinganine"],
             smiles=physio(PC["ketosphinganine"]), ring={"coordgen": False, "rotate": 180}, tags=SP,
             notes="serine palmitoyltransferase（PLP，ER，限速）：palmitoyl-CoA + Ser → 縮合並脫羧。"),
        dict(id="sphinganine", name="Sphinganine", aliases=["dihydrosphingosine"], smiles=physio(PC["sphinganine"]), ring=ring, tags=SP,
             notes="3-ketosphinganine reductase（NADPH）。接著 ceramide synthase 接上脂肪醯基（N-acyl）→ dihydroceramide → desaturase 產生 trans-Δ4 → ceramide。"),
        dict(id="sphingosine", name="Sphingosine", aliases=["鞘胺醇"], smiles=physio(PC["sphingosine"]), ring=ring, tags=SP,
             notes="ceramidase 分解 ceramide 的產物（Farber disease：acid ceramidase 缺乏）。磷酸化 → S1P（淋巴球移出、血管新生；fingolimod 作用點）。"),
        dict(id="glucocerebroside", name="Glucocerebroside", aliases=["glucosylceramide", "GlcCer"],
             smiles=CER.format(glc()), lipid={"type": "sphingolipid", "chains": [{"n": 18}], "head": "glucose"}, tags=SP,
             notes="ceramide + β-glucose（UDP-glucose 提供）。大多數鞘醣脂與神經節苷脂的起點。Gaucher disease：glucocerebrosidase 缺乏（最常見的溶小體儲積症；肝脾腫大、骨痛、Gaucher cell 似皺紙；enzyme replacement imiglucerase）。"),
        dict(id="sulfatide", name="Sulfatide", aliases=["galactocerebroside 3-sulfate", "sulfogalactosylceramide"],
             smiles=CER.format("C1OC(CO)C(O)C(OS(=O)(=O)[O-])C1O"),
             lipid={"type": "sphingolipid", "chains": [{"n": 18}], "head": "sulfogalactose"}, tags=SP,
             notes="galactocerebroside 的 Gal C3 接硫酸（PAPS 提供）。髓鞘成分。metachromatic leukodystrophy：arylsulfatase A 缺乏 → 中樞與周邊脫髓鞘。"),
        dict(id="lactosylceramide", name="Lactosylceramide", aliases=["LacCer"], smiles=GSL["lactosylceramide"],
             free=snfg_cartoon(["glc", "gal"]), tags=SP, notes="ceramide–Glc–Gal（β1→4）。GM3 與 Gb3 的共同前驅物。{SNFG}"),
        dict(id="globotriaosylceramide", name="Globotriaosylceramide (Gb3)", label="Globotriaosyl-\nceramide (Gb3)", aliases=["Gb3", "ceramide trihexoside"],
             smiles=GSL["globotriaosylceramide"], free=snfg_cartoon(["glc", "gal", "gal"]), tags=SP,
             notes="ceramide–Glc–Gal–Gal（末端 α1→4）。Fabry disease：α-galactosidase A 缺乏（X-linked）→ 末梢神經痛、angiokeratoma、少汗、腎衰竭、心肌病。也是 Shiga toxin 受體。{SNFG}"),
        dict(id="gm3", name="Ganglioside GM3", label="GM3", aliases=["GM3"], smiles=GSL["gm3"],
             free=snfg_cartoon(["glc", "gal"], branch=(1, "neu5ac")), tags=SP, notes="ceramide–Glc–Gal(–NeuAc)。最簡單的 ganglioside。{SNFG}"),
        dict(id="gm2", name="Ganglioside GM2", label="GM2", aliases=["GM2"], smiles=GSL["gm2"],
             free=snfg_cartoon(["glc", "gal", "galnac"], branch=(1, "neu5ac")), tags=SP,
             notes="GM3 + GalNAc。Tay-Sachs：hexosaminidase A 缺乏 → GM2 堆積（神經退化、黃斑櫻桃紅點、無肝脾腫大，Ashkenazi Jews）；Sandhoff：hexosaminidase A + B（β 次單元），另有 globoside 堆積、可有肝脾腫大。{SNFG}"),
        dict(id="gm1", name="Ganglioside GM1", label="GM1", aliases=["GM1"], smiles=GSL["gm1"],
             free=snfg_cartoon(["glc", "gal", "galnac", "gal"], branch=(1, "neu5ac")), tags=SP,
             notes="GM2 的 GalNAc 再接 Gal。cholera toxin 的受體。GM1 gangliosidosis：β-galactosidase 缺乏。{SNFG}"),
        dict(id="galnac", name="N-Acetylgalactosamine", label="GalNAc", aliases=["GalNAc", "β-GalNAc"],
             smiles="OC1OC(CO)C(O)C(O)C1NC(C)=O", free={"atoms": {"s": ["@galnac", 0, 0]}, "bonds": [], "schematic": True}, tags=SP,
             notes="hexosaminidase 從 GM2 切下的末端糖。{SNFG}"),
        dict(id="phosphocholine", name="Phosphocholine", aliases=["phosphorylcholine"], smiles="C[N+](C)(C)CCOP(=O)([O-])[O-]", tags=SP,
             notes="sphingomyelinase 的產物；合成 sphingomyelin 時則由 PC 提供（ceramide + PC → SM + DAG）。"),
        dict(id="sulfate", name="Sulfate", aliases=["SO4 2-"], smiles="[O-]S(=O)(=O)[O-]", formula="SO_4^{2-}", tags=SP,
             notes="arylsulfatase A 從 sulfatide 切下。"),
    ]
    for m in out:
        Chem.MolFromSmiles(m["smiles"])
        m["notes"] = m["notes"].replace("{SNFG}", SNFG_NOTE)
    return out


COFACTORS = {
    "GSH": ("GSH", "S[16*]", None),                  # [16*] = the rest of glutathione (γ-Glu-Cys-Gly)
    "GSSG": ("GSSG", "[16*]SS[16*]", None),
}
MULTI = {"2GSH": ("2 GSH", "S[16*]", 2), "2O2": ("2 O_2", "O=O", 2), "3O2": ("3 O_2", "O=O", 3),
         "3H+": ("3 H^+", "[H+]", 3), "5H2O": ("5 H_2O", "O", 5), "3H2O": ("3 H_2O", "O", 3)}
MARK = "# --- 脂質代謝延伸（第十二輪）"

CYP = dict(cofactors_in=["O2", "NADPH+H+"], cofactors_out=["NADP+", "H2O"])

REACTIONS = [
    # ---- eicosanoids
    dict(id="cyclooxygenase", title="Arachidonate → PGG₂ (COX, cyclooxygenase activity)",
         substrates=["arachidonate"], products=["pgg2"], cofactors_in=["2O2"], enzyme="cyclooxygenase\n(COX-1, COX-2)",
         card={"location": "ER 膜（PGH synthase）", "type": "雙氧化（環化＋15-OOH）",
               "regulation": "COX-1 常駐（胃黏膜、血小板、腎）；COX-2 誘導型（發炎、cytokine）",
               "deficiency": ["aspirin：不可逆乙醯化 COX Ser530（低劑量 → 血小板 TXA₂ ↓）",
                              "NSAID（ibuprofen、naproxen）：可逆競爭抑制；celecoxib：選擇性 COX-2（胃副作用少、心血管風險 ↑）",
                              "corticosteroid：抑制 PLA₂（annexin/lipocortin）並減少 COX-2 表現"],
               "exam_points": ["arachidonate 先由 cPLA₂ 從膜磷脂 sn-2 釋出（限速）"]},
         tags=["eicosanoid_synthesis"]),
    dict(id="cox_peroxidase", title="PGG₂ → PGH₂ (COX, peroxidase activity)",
         substrates=["pgg2"], products=["pgh2"], cofactors_in=["2GSH"], cofactors_out=["GSSG", "H2O"],
         enzyme="cyclooxygenase\n(peroxidase site)",
         card={"location": "ER 膜（同一個 PGH synthase 的另一個活性位）", "type": "還原 15-OOH → 15-OH（電子由還原劑提供，圖中以 GSH 表示）",
               "exam_points": ["aspirin 只抑制環氧酶活性，不抑制過氧化酶活性"]},
         tags=["eicosanoid_synthesis"]),
    dict(id="prostacyclin_synthase", title="PGH₂ → PGI₂ (prostacyclin synthase)", substrates=["pgh2"], products=["pgi2"],
         enzyme="prostacyclin synthase", card={"location": "血管內皮", "type": "異構化（CYP8A1）",
                                                "exam_points": ["PGI₂ 與 TXA₂ 的平衡決定血小板凝集；ω-3（EPA）讓 TXA₃ 較弱 → 抗血栓"]},
         tags=["eicosanoid_synthesis"]),
    dict(id="thromboxane_synthase", title="PGH₂ → TXA₂ (thromboxane synthase)", substrates=["pgh2"], products=["txa2"],
         enzyme="thromboxane synthase", card={"location": "血小板", "type": "異構化（CYP5A1）",
                                               "exam_points": ["TXA₂：血管收縮＋血小板凝集；PGI₂：相反"]},
         tags=["eicosanoid_synthesis"]),
    dict(id="pge_synthase", title="PGH₂ → PGE₂ (PGE synthase)", substrates=["pgh2"], products=["pge2"],
         enzyme="PGE synthase", card={"type": "異構化", "exam_points": ["PGE₂：發燒（下視丘）、疼痛、血管擴張、胃黏液 ↑ 胃酸 ↓、維持 PDA 開放（indomethacin 關閉 PDA）"]},
         tags=["eicosanoid_synthesis"]),
    dict(id="lipoxygenase_5", title="Arachidonate → 5-HPETE (5-lipoxygenase)", substrates=["arachidonate"], products=["hpete5"],
         cofactors_in=["O2"], enzyme="5-lipoxygenase\n(5-LOX + FLAP)",
         card={"location": "白血球（嗜中性球、肥大細胞、嗜酸性球）", "type": "雙氧化",
               "deficiency": ["zileuton 抑制 5-LOX（氣喘）"]},
         tags=["eicosanoid_synthesis"]),
    dict(id="lta4_synthase", title="5-HPETE → LTA₄ (5-lipoxygenase, dehydrase activity)", substrates=["hpete5"], products=["lta4"],
         cofactors_out=["H2O"], enzyme="5-lipoxygenase", card={"type": "脫水形成 5,6-環氧化物"}, tags=["eicosanoid_synthesis"]),
    dict(id="lta4_hydrolase", title="LTA₄ → LTB₄ (LTA₄ hydrolase)", substrates=["lta4"], products=["ltb4"], cofactors_in=["H2O"],
         enzyme="LTA_4 hydrolase", card={"location": "嗜中性球", "type": "環氧化物水解", "exam_points": ["LTB₄：嗜中性球趨化、黏附"]},
         tags=["eicosanoid_synthesis"]),
    dict(id="ltc4_synthase", title="LTA₄ + GSH → LTC₄ (LTC₄ synthase)", substrates=["lta4"], products=["ltc4"], cofactors_in=["GSH"],
         enzyme="LTC_4 synthase",
         card={"location": "肥大細胞、嗜酸性球", "type": "與 glutathione 結合（開環）",
               "exam_points": ["LTC₄ → LTD₄ → LTE₄：支氣管收縮（比 histamine 強且久）；montelukast 阻斷 CysLT1 受體"]},
         tags=["eicosanoid_synthesis"]),
    # ---- steroid hormones
    dict(id="cyp11a1", title="Cholesterol → Pregnenolone (CYP11A1, side-chain cleavage)",
         substrates=["cholesterol"], products=["pregnenolone", "isocaproaldehyde"], cofactors_in=["3O2", "3NADPH"],
         cofactors_out=["3NADP+", "4H2O"], extra_in=["3H+"], enzyme="CYP11A1 (desmolase,\nside-chain cleavage)",
         card={"location": "粒線體內膜（adrenodoxin 傳電子）", "type": "三次單氧化：22R-OH、20,22-diOH，再切斷 C20–C22",
               "regulation": "ACTH（腎上腺）、LH（Leydig、黃體）經 cAMP；真正的限速是 StAR 把膽固醇送到粒線體內膜",
               "deficiency": ["StAR 或 CYP11A1 缺乏 → lipoid CAH：所有固醇類荷爾蒙 ↓，腎上腺充滿膽固醇"]},
         tags=["steroid_hormone"]),
    dict(id="hsd3b", title="Pregnenolone → Progesterone (3β-HSD)", substrates=["pregnenolone"], products=["progesterone"],
         cofactors_in=["NAD+"], cofactors_out=["NADH+H+"], enzyme="3β-hydroxysteroid dehydrogenase\n/Δ5,4-isomerase",
         card={"location": "ER", "type": "3β-OH 氧化成酮＋雙鍵由 Δ5 移到 Δ4",
               "deficiency": ["3β-HSD 缺乏：所有活性荷爾蒙 ↓、DHEA ↑；男嬰外生殖器模糊、女嬰輕度男性化、失鹽"]},
         tags=["steroid_hormone"]),
    dict(id="cyp17_hydroxylase", title="Progesterone → 17α-Hydroxyprogesterone (CYP17A1, 17α-hydroxylase)",
         substrates=["progesterone"], products=["oh_progesterone"], **CYP, enzyme="CYP17A1\n(17α-hydroxylase)",
         card={"location": "ER（束狀帶、網狀帶、性腺）", "type": "單氧化",
               "deficiency": ["17α-hydroxylase 缺乏：cortisol 與性荷爾蒙 ↓，DOC ↑ → 高血壓、低血鉀；男性外生殖器女性化、女性青春期不發育"],
               "exam_points": ["ketoconazole、abiraterone（前列腺癌）抑制 CYP17"]},
         tags=["steroid_hormone"]),
    dict(id="cyp17_lyase", title="17α-Hydroxypregnenolone → DHEA (CYP17A1, 17,20-lyase)",
         substrates=["oh_pregnenolone"], products=["dhea", "acetate"], **CYP, extra_out=["H+"], enzyme="CYP17A1\n(17,20-lyase)",
         card={"location": "ER（網狀帶、性腺）", "type": "切斷 C17–C20（需 cytochrome b5）",
               "exam_points": ["人類主要走 Δ5 路徑：17-OH pregnenolone → DHEA"]},
         tags=["steroid_hormone"]),
    dict(id="cyp21_doc", title="Progesterone → 11-Deoxycorticosterone (CYP21A2, 21-hydroxylase)",
         substrates=["progesterone"], products=["doc"], **CYP, enzyme="CYP21A2\n(21-hydroxylase)",
         card={"location": "ER（球狀帶、束狀帶）", "type": "單氧化（C21）",
               "deficiency": ["21-hydroxylase 缺乏（CAH 最常見，約 90%）：cortisol 與 aldosterone ↓、17-OHP ↑、雄性素 ↑ → 女嬰男性化、失鹽型（低血鈉、高血鉀、低血壓）"]},
         tags=["steroid_hormone"]),
    dict(id="cyp21_deoxycortisol", title="17α-Hydroxyprogesterone → 11-Deoxycortisol (CYP21A2)",
         substrates=["oh_progesterone"], products=["deoxycortisol"], **CYP, enzyme="CYP21A2\n(21-hydroxylase)",
         card={"location": "ER（束狀帶）", "type": "單氧化（C21）", "exam_points": ["新生兒篩檢以 17-OHP 偵測 21-hydroxylase 缺乏"]},
         tags=["steroid_hormone"]),
    dict(id="cyp11b1", title="11-Deoxycortisol → Cortisol (CYP11B1, 11β-hydroxylase)",
         substrates=["deoxycortisol"], products=["cortisol"], **CYP, enzyme="CYP11B1\n(11β-hydroxylase)",
         card={"location": "粒線體（束狀帶）", "type": "單氧化（C11）",
               "deficiency": ["11β-hydroxylase 缺乏：cortisol ↓、11-deoxycortisol 與 DOC ↑ → 高血壓、低血鉀＋男性化"],
               "exam_points": ["metyrapone 抑制"]},
         tags=["steroid_hormone"]),
    dict(id="cyp11b_corticosterone", title="11-Deoxycorticosterone → Corticosterone (11β-hydroxylase)",
         substrates=["doc"], products=["corticosterone"], **CYP, enzyme="CYP11B1 / CYP11B2\n(11β-hydroxylase)",
         card={"location": "粒線體", "type": "單氧化（C11）"}, tags=["steroid_hormone"]),
    dict(id="aldosterone_synthase", title="Corticosterone → Aldosterone (CYP11B2, aldosterone synthase)",
         substrates=["corticosterone"], products=["aldosterone"], cofactors_in=["2O2", "2NADPH+H+"], cofactors_out=["2NADP+", "3H2O"],
         enzyme="CYP11B2\n(aldosterone synthase)",
         card={"location": "粒線體（只在球狀帶）", "type": "18-hydroxylation → 18-oxidation（醛）",
               "regulation": "angiotensin II、高血鉀（ACTH 只有短暫作用）"},
         tags=["steroid_hormone"]),
    dict(id="hsd17b", title="Androstenedione → Testosterone (17β-HSD)", substrates=["androstenedione"], products=["testosterone"],
         cofactors_in=["NADPH+H+"], cofactors_out=["NADP+"], enzyme="17β-hydroxysteroid\ndehydrogenase",
         card={"location": "Leydig 細胞（17β-HSD3）、周邊組織", "type": "17-酮還原成 17β-OH",
               "deficiency": ["17β-HSD3 缺乏：男性 46,XY 外生殖器女性化，青春期部分男性化"]},
         tags=["steroid_hormone"]),
    dict(id="aromatase", title="Testosterone → Estradiol (aromatase, CYP19A1)", substrates=["testosterone"], products=["estradiol", "formate"],
         cofactors_in=["3O2", "3NADPH"], cofactors_out=["3NADP+", "4H2O"], extra_in=["2H+"], enzyme="aromatase\n(CYP19A1)",
         card={"location": "顆粒細胞（FSH 誘導）、脂肪、胎盤、腦", "type": "三次單氧化：切掉 C19（甲酸）並使 A 環芳香化",
               "deficiency": ["aromatase 抑制劑（anastrozole、letrozole）：停經後乳癌", "aromatase 缺乏：女胎男性化（雄性素累積）、高個子（骨骺不閉合）"]},
         tags=["steroid_hormone"]),
    dict(id="srd5a", title="Testosterone → DHT (5α-reductase)", substrates=["testosterone"], products=["dht"],
         cofactors_in=["NADPH+H+"], cofactors_out=["NADP+"], enzyme="5α-reductase",
         card={"location": "前列腺、皮膚、毛囊、外生殖器", "type": "還原 C4=C5",
               "deficiency": ["finasteride（5α-R2）、dutasteride（1 + 2）：BPH、雄性禿"]},
         tags=["steroid_hormone"]),
    # ---- bile acids
    dict(id="cyp7a1", title="Cholesterol → 7α-Hydroxycholesterol (CYP7A1)", substrates=["cholesterol"], products=["oh7_cholesterol"],
         **CYP, enzyme="cholesterol 7α-hydroxylase\n(CYP7A1)",
         card={"location": "肝 ER", "type": "單氧化（膽酸合成的限速步驟）",
               "regulation": "膽酸 → FXR → SHP 抑制 CYP7A1；膽固醇（LXR）在囓齒類誘導",
               "exam_points": ["每天約 0.5 g 膽固醇轉成膽酸 → 排出膽固醇的主要途徑", "膽酸結合劑（cholestyramine）→ CYP7A1 ↑ → LDL ↓"]},
         tags=["bile_acid"]),
    dict(id="bile_7a_dehydroxylation", title="Cholate → Deoxycholate (bacterial 7α-dehydroxylation)",
         substrates=["cholate"], products=["deoxycholate"], cofactors_in=["NADH+H+"], cofactors_out=["NAD+", "H2O"],
         enzyme="gut bacteria\n(7α-dehydroxylase)",
         card={"location": "大腸（腸道菌）", "type": "去羥（先去結合，再脫去 7α-OH）",
               "exam_points": ["初級膽酸（肝）：cholate、chenodeoxycholate；次級（細菌）：deoxycholate、lithocholate"]},
         tags=["bile_acid"]),
    dict(id="bile_7a_dehydroxylation_cdca", title="Chenodeoxycholate → Lithocholate (bacterial 7α-dehydroxylation)",
         substrates=["chenodeoxycholate"], products=["lithocholate"], cofactors_in=["NADH+H+"], cofactors_out=["NAD+", "H2O"],
         enzyme="gut bacteria\n(7α-dehydroxylase)", card={"location": "大腸（腸道菌）", "type": "去羥",
                                                          "exam_points": ["lithocholate 有肝毒性，大部分隨糞便排出"]},
         tags=["bile_acid"]),
    # ---- sphingolipids
    dict(id="serine_palmitoyltransferase", title="Palmitoyl-CoA + Serine → 3-Ketosphinganine (serine palmitoyltransferase)",
         substrates=["palmitoyl_coa", "serine"], products=["ketosphinganine", {"id": "co2", "view": "formula"}],
         cofactors_out=["CoA-SH"], extra_in=["H+"], prosthetic=["PLP"], enzyme="serine\npalmitoyltransferase",
         card={"location": "ER 細胞質面", "type": "縮合＋脫羧（PLP）；鞘脂合成的限速步驟",
               "deficiency": ["SPTLC1 突變 → 遺傳性感覺神經病變（HSAN1）"]},
         tags=["sphingolipid_metabolism"]),
    dict(id="ketosphinganine_reductase", title="3-Ketosphinganine → Sphinganine (3-ketosphinganine reductase)",
         substrates=["ketosphinganine"], products=["sphinganine"], cofactors_in=["NADPH+H+"], cofactors_out=["NADP+"],
         enzyme="3-ketosphinganine\nreductase", card={"location": "ER", "type": "還原（NADPH）"}, tags=["sphingolipid_metabolism"]),
    dict(id="glucocerebrosidase", title="Glucocerebroside → Ceramide + Glucose (glucocerebrosidase)",
         substrates=["glucocerebroside"], products=["ceramide", "b_glucopyranose"], cofactors_in=["H2O"],
         enzyme="glucocerebrosidase\n(β-glucosidase, GBA)",
         card={"location": "溶小體", "type": "醣苷鍵水解",
               "deficiency": ["Gaucher disease（最常見的溶小體儲積症，AR）：肝脾腫大、全血球低下、骨痛與骨壞死（Erlenmeyer flask）、Gaucher cell（皺紙狀胞質）；治療 imiglucerase、miglustat",
                              "GBA 帶因者 → Parkinson 風險 ↑"]},
         tags=["sphingolipid_metabolism"]),
    dict(id="sphingomyelinase", title="Sphingomyelin → Ceramide + Phosphocholine (sphingomyelinase)",
         substrates=["sphingomyelin"], products=["ceramide", "phosphocholine"], cofactors_in=["H2O"], extra_out=["H+"],
         enzyme="acid sphingomyelinase", card={"location": "溶小體", "type": "磷酸二酯鍵水解",
                                               "deficiency": ["Niemann-Pick A/B（AR）：肝脾腫大、黃斑櫻桃紅點、foam cell；A 型嚴重神經退化（Ashkenazi Jews）",
                                                              "Niemann-Pick C：NPC1/2 膽固醇運輸缺陷（不是此酵素）"]},
         tags=["sphingolipid_metabolism"]),
    dict(id="galactocerebrosidase", title="Galactocerebroside → Ceramide + Galactose (galactocerebrosidase)",
         substrates=["galactocerebroside"], products=["ceramide", "b_galactopyranose"], cofactors_in=["H2O"],
         enzyme="galactocerebrosidase\n(GALC)", card={"location": "溶小體", "type": "醣苷鍵水解",
                                                    "deficiency": ["Krabbe disease（globoid cell leukodystrophy）：psychosine 累積毒殺寡突膠細胞 → 脫髓鞘、視神經萎縮、globoid cell"]},
         tags=["sphingolipid_metabolism"]),
    dict(id="arylsulfatase_a", title="Sulfatide → Galactocerebroside + Sulfate (arylsulfatase A)",
         substrates=["sulfatide"], products=["galactocerebroside", "sulfate"], cofactors_in=["H2O"], extra_out=["H+"],
         enzyme="arylsulfatase A", card={"location": "溶小體", "type": "硫酸酯水解",
                                         "deficiency": ["metachromatic leukodystrophy（AR）：中樞與周邊脫髓鞘、運動失調、失智；尿中 sulfatide ↑"]},
         tags=["sphingolipid_metabolism"]),
    dict(id="hexosaminidase_a", title="GM2 → GM3 + GalNAc (β-hexosaminidase A)", substrates=["gm2"], products=["gm3", "galnac"],
         cofactors_in=["H2O"], enzyme="β-hexosaminidase A\n(+ GM2 activator)",
         card={"location": "溶小體", "type": "醣苷鍵水解",
               "deficiency": ["Tay-Sachs（HEXA，α 次單元）：進行性神經退化、黃斑櫻桃紅點、對聲音過度驚嚇；無肝脾腫大；Ashkenazi Jews",
                              "Sandhoff（HEXB，β 次單元 → Hex A + B 皆缺）：另有 globoside 累積，可有肝脾腫大",
                              "GM2 activator 缺乏（AB variant）"]},
         tags=["sphingolipid_metabolism"]),
    dict(id="alpha_galactosidase_a", title="Gb3 → Lactosylceramide + Galactose (α-galactosidase A)",
         substrates=["globotriaosylceramide"], products=["lactosylceramide", "a_galactopyranose"], cofactors_in=["H2O"],
         enzyme="α-galactosidase A", card={"location": "溶小體", "type": "醣苷鍵水解（α1→4）",
                                           "deficiency": ["Fabry disease（X-linked 隱性，唯一 X-linked 的鞘脂儲積症＋Hunter）：肢端疼痛、angiokeratoma、少汗、角膜混濁、蛋白尿 → 腎衰竭、心肌病、中風；治療 agalsidase"]},
         tags=["sphingolipid_metabolism"]),
    dict(id="beta_galactosidase_gm1", title="GM1 → GM2 + Galactose (β-galactosidase)", substrates=["gm1"], products=["gm2", "b_galactopyranose"],
         cofactors_in=["H2O"], enzyme="β-galactosidase (GLB1)",
         card={"location": "溶小體", "type": "醣苷鍵水解", "deficiency": ["GM1 gangliosidosis：神經退化、櫻桃紅點、肝脾腫大、骨骼異常"]},
         tags=["sphingolipid_metabolism"]),
]


def n(id_, ref, x, y, **kw):
    return {"id": id_, "ref": ref, "pos": [x, y], **kw}


def t(id_, text, x, y, kind="concept", **kw):
    return {"id": id_, "text": text, "kind": kind, "pos": [x, y], **kw}


NAME = {"view": "name"}

PATHWAYS = [
    dict(id="eicosanoid_synthesis", title="Eicosanoid synthesis（類二十烷酸合成：COX 與 LOX）", layout="free", cell_size=3.4,
         nodes=[t("pl", "Membrane phospholipid\n(arachidonate at sn-2)", 0, -1.6, "name"),
                n("aa", "arachidonate", 0, 0, **NAME),
                t("ster", "Corticosteroids", -3.6, -0.8, "drug"),
                n("pgg2", "pgg2", -3.0, 1.5, **NAME), n("pgh2", "pgh2", -3.0, 3.0, **NAME),
                n("pgi2", "pgi2", -6.2, 4.6, **NAME), n("txa2", "txa2", -3.9, 4.6, **NAME),
                n("pge2", "pge2", -1.8, 4.6, **NAME), n("pgf2a", "pgf2a", 0.1, 4.6, **NAME),
                t("pgi2e", "endothelium: vasodilation,\nplatelet aggregation ↓", -6.2, 5.7, "name"),
                t("txa2e", "platelets: vasoconstriction,\naggregation ↑", -3.5, 6.6, "name"),
                t("pge2e", "fever, pain, PDA open,\ngastric mucus ↑", -0.9, 5.7, "name"),
                t("nsaid", "Aspirin (irreversible),\nNSAIDs, celecoxib (COX-2)", -6.0, 1.0, "drug"),
                n("hpete", "hpete5", 3.2, 1.5, **NAME), n("lta4", "lta4", 3.2, 3.0, **NAME),
                n("ltb4", "ltb4", 1.9, 4.6, **NAME), n("ltc4", "ltc4", 4.6, 4.6, **NAME),
                t("ltb4e", "neutrophil chemotaxis", 1.9, 5.6, "name"),
                t("ltc4e", "LTC_4 → LTD_4 → LTE_4:\nbronchoconstriction", 4.6, 5.7, "name"),
                t("zil", "Zileuton", 5.8, 1.0, "drug"), t("mon", "Montelukast,\nzafirlukast (receptor)", 6.9, 4.6, "drug")],
         edges=[{"from": "pl", "to": "aa", "label": "phospholipase A_2"},
                {"from": "ster", "to": "edge:pl->aa", "inhibit": True},
                {"from": "aa", "to": "pgg2", "step": "cyclooxygenase", "label": "COX", "side_text": "2 O_2"},
                {"from": "pgg2", "to": "pgh2", "step": "cox_peroxidase", "label": "peroxidase"},
                {"from": "pgh2", "to": "pgi2", "step": "prostacyclin_synthase"}, {"from": "pgh2", "to": "txa2", "step": "thromboxane_synthase"},
                {"from": "pgh2", "to": "pge2", "step": "pge_synthase"}, {"from": "pgh2", "to": "pgf2a"},
                {"from": "nsaid", "to": "edge:aa->pgg2", "inhibit": True},
                {"from": "aa", "to": "hpete", "step": "lipoxygenase_5", "label": "5-LOX", "side_text": "O_2"},
                {"from": "hpete", "to": "lta4", "step": "lta4_synthase"},
                {"from": "lta4", "to": "ltb4", "step": "lta4_hydrolase"}, {"from": "lta4", "to": "ltc4", "step": "ltc4_synthase", "side_text": "GSH"},
                {"from": "zil", "to": "edge:aa->hpete", "inhibit": True}, {"from": "mon", "to": "ltc4", "inhibit": True}],
         tags=["eicosanoid_synthesis"],
         notes="cPLA₂ 從膜磷脂 sn-2 釋出 arachidonate（20:4 ω-6）。COX 路徑 → prostanoid（PG、TX）；LOX 路徑 → leukotriene、lipoxin。PGI₂（內皮）與 TXA₂（血小板）作用相反。aspirin 不可逆乙醯化 COX；低劑量時只有無核的血小板無法再生 COX-1 → 選擇性 TXA₂ ↓。aspirin-exacerbated respiratory disease：COX 被抑制 → 更多 arachidonate 走 LOX → cysteinyl LT ↑。"),
    dict(id="steroidogenesis", title="Adrenal steroidogenesis & CAH（腎上腺固醇類荷爾蒙合成與 CAH）", layout="free", cell_size=3.6,
         nodes=[t("hg", "Zona glomerulosa\n(mineralocorticoid)", -3.0, -2.0, "name"),
                t("hf", "Zona fasciculata\n(glucocorticoid)", 0.5, -2.0, "name"),
                t("hr", "Zona reticularis / gonads\n(androgens, estrogens)", 4.6, -2.0, "name"),
                n("chol", "cholesterol", -3.0, -0.9, **NAME), n("preg", "pregnenolone", -3.0, 0.8, **NAME),
                n("ohpreg", "oh_pregnenolone", 0.5, 0.8, **NAME), n("dhea", "dhea", 4.0, 0.8, **NAME),
                n("prog", "progesterone", -3.0, 2.3, **NAME), n("ohprog", "oh_progesterone", 0.5, 2.3, **NAME),
                n("andro", "androstenedione", 4.0, 2.3, **NAME),
                n("doc", "doc", -3.0, 3.8, **NAME), n("dcort", "deoxycortisol", 0.5, 3.8, **NAME), n("testo", "testosterone", 4.0, 3.8, **NAME),
                n("cortico", "corticosterone", -3.0, 5.3, **NAME), n("cort", "cortisol", 0.5, 5.3, **NAME), n("e2", "estradiol", 4.0, 5.3, **NAME),
                n("aldo", "aldosterone", -3.0, 6.8, **NAME), n("dht", "dht", 6.4, 3.8, **NAME), n("estrone", "estrone", 6.4, 2.3, **NAME),
                t("d21", "21-OH deficiency (~90% CAH)", -5.3, 3.0, "drug"), t("d11", "11β-OH deficiency", -5.3, 4.5, "drug"),
                t("d17", "17α-OH deficiency", -1.2, -0.2, "drug"),
                t("acth", "ACTH ↑ in CAH (cortisol ↓)", 0.5, 7.0, "concept")],
         edges=[{"from": "chol", "to": "preg", "step": "cyp11a1", "label": "CYP11A1 (desmolase)"},
                {"from": "preg", "to": "ohpreg", "label": "17α-hydroxylase", "label_side": "right"}, {"from": "ohpreg", "to": "dhea", "step": "cyp17_lyase", "label": "17,20-lyase"},
                {"from": "preg", "to": "prog", "step": "hsd3b", "label": "3β-HSD"}, {"from": "ohpreg", "to": "ohprog", "label": "3β-HSD"},
                {"from": "dhea", "to": "andro", "label": "3β-HSD"},
                {"from": "prog", "to": "ohprog", "step": "cyp17_hydroxylase", "label": "17α-hydroxylase"},
                {"from": "ohprog", "to": "andro", "label": "17,20-lyase"},
                {"from": "prog", "to": "doc", "step": "cyp21_doc", "label": "21-hydroxylase"},
                {"from": "ohprog", "to": "dcort", "step": "cyp21_deoxycortisol", "label": "21-hydroxylase"},
                {"from": "andro", "to": "testo", "step": "hsd17b", "label": "17β-HSD"},
                {"from": "doc", "to": "cortico", "step": "cyp11b_corticosterone", "label": "11β-hydroxylase"},
                {"from": "dcort", "to": "cort", "step": "cyp11b1", "label": "11β-hydroxylase"},
                {"from": "testo", "to": "e2", "step": "aromatase", "label": "aromatase"},
                {"from": "andro", "to": "estrone", "label": "aromatase"},
                {"from": "testo", "to": "dht", "step": "srd5a", "label": "5α-reductase"},
                {"from": "cortico", "to": "aldo", "step": "aldosterone_synthase", "label": "aldosterone synthase"},
                {"from": "d21", "to": "edge:prog->doc", "inhibit": True}, {"from": "d11", "to": "edge:doc->cortico", "inhibit": True},
                {"from": "d17", "to": "edge:preg->ohpreg", "inhibit": True}],
         tags=["steroid_hormone"],
         notes="CAH 三型比較（新生兒、皆 cortisol ↓ → ACTH ↑ → 腎上腺增生）：21-hydroxylase 缺乏（90%）：aldosterone ↓（失鹽、低血壓、高血鉀）、17-OHP ↑、雄性素 ↑（女嬰男性化）；11β-hydroxylase 缺乏：DOC ↑ → 高血壓、低血鉀＋雄性素 ↑；17α-hydroxylase 缺乏：性荷爾蒙與 cortisol ↓、DOC ↑ → 高血壓、低血鉀、性徵不發育。口訣：數字以 1 開頭（11β、17α）→ 高血壓；以 1 結尾（21、11）→ 男性化。"),
    dict(id="bile_acid_synthesis", title="Bile acid synthesis & enterohepatic circulation（膽酸合成與腸肝循環）", layout="free", cell_size=3.4,
         regions=[{"kind": "box", "x0": -4.2, "y0": -0.9, "x1": 4.4, "y1": 5.2, "label": "Liver"},
                  {"kind": "box", "x0": -4.2, "y0": 6.0, "x1": 4.4, "y1": 9.6, "label": "Intestine (bacteria)"}],
         nodes=[n("chol", "cholesterol", 0, 0.2, **NAME), n("oh7", "oh7_cholesterol", 0, 1.6, **NAME),
                n("ca", "cholate", -2.2, 3.0, **NAME), n("cdca", "chenodeoxycholate", 2.2, 3.0, **NAME),
                n("gca", "glycocholate", -3.0, 4.4, **NAME), n("tca", "taurocholate", -0.6, 4.4, **NAME),
                t("cdcac", "glyco-/tauro-CDCA", 2.4, 4.4, "name"),
                n("dca", "deoxycholate", -2.2, 7.6, **NAME), n("lca", "lithocholate", 2.2, 7.6, **NAME),
                t("ileum", "~95% reabsorbed (terminal ileum, ASBT)\n→ portal vein → liver", 6.8, 4.0, "concept"),
                t("feces", "~5% (≈ 0.5 g/day) lost in feces", 0, 9.0, "name"),
                t("fxr", "FXR (bile acids)", -3.6, 0.9, "drug"), t("chol_res", "Cholestyramine\n(sequestrant)", 6.8, 7.4, "drug")],
         edges=[{"from": "chol", "to": "oh7", "step": "cyp7a1", "label": "CYP7A1 (rate-limiting)"},
                {"from": "oh7", "to": "ca", "label": "CYP8B1"}, {"from": "oh7", "to": "cdca"},
                {"from": "ca", "to": "gca"}, {"from": "ca", "to": "tca"},
                {"from": "cdca", "to": "cdcac"},
                {"from": "ca", "to": "dca", "step": "bile_7a_dehydroxylation", "bend": 0.12},
                {"from": "cdca", "to": "lca", "step": "bile_7a_dehydroxylation_cdca", "bend": -0.12,
                 "label": "gut bacteria: deconjugate,\nremove 7α-OH", "label_side": "left"},
                {"from": "fxr", "to": "edge:chol->oh7", "inhibit": True},
                {"from": "ileum", "to": "chol", "bend": 0.15},
                {"from": "chol_res", "to": "ileum", "inhibit": True}],
         tags=["bile_acid"],
         notes="初級膽酸（肝）：cholate（3α,7α,12α）、chenodeoxycholate（3α,7α）；與 glycine 或 taurine 結合 → 膽鹽（兩性、乳化脂質、形成 micelle、吸收脂溶性維生素）。腸道細菌去結合並去掉 7α-OH → 次級膽酸 deoxycholate、lithocholate。約 95% 在迴腸末端回收（腸肝循環，每天 6–10 次）。迴腸切除或 Crohn → 膽酸流失 → 脂肪瀉、膽結石（膽固醇過飽和）。ursodeoxycholate 治療膽固醇結石、PBC。"),
    dict(id="sphingolipid_synthesis", title="Sphingolipid synthesis（鞘脂合成）", layout="free", cell_size=3.4,
         nodes=[n("pal", "palmitoyl_coa", -2.2, 0, **NAME), n("ser", "serine", 2.2, 0, **NAME),
                n("ks", "ketosphinganine", 0, 1.4, **NAME), n("sa", "sphinganine", 0, 2.8, **NAME),
                t("dhc", "Dihydroceramide", 0, 4.2, "name"), n("cer", "ceramide", 0, 5.6, **NAME),
                n("sm", "sphingomyelin", -4.4, 7.4, **NAME), n("glccer", "glucocerebroside", 0, 7.4, **NAME),
                n("galcer", "galactocerebroside", 4.4, 7.4, **NAME), n("sulf", "sulfatide", 4.4, 9.0, **NAME),
                n("lac", "lactosylceramide", 0, 9.0, **NAME),
                t("gang", "Gangliosides (GM3 → GM2 → GM1 …),\nglobosides (Gb3)", 0, 10.6, "name")],
         edges=[{"from": "pal", "to": "ks", "step": "serine_palmitoyltransferase", "label": "serine palmitoyltransferase\n(PLP, rate-limiting)"},
                {"from": "ser", "to": "ks", "step": "serine_palmitoyltransferase", "side_text": "CO_2 + CoA"},
                {"from": "ks", "to": "sa", "step": "ketosphinganine_reductase", "side_text": "NADPH"},
                {"from": "sa", "to": "dhc", "label": "ceramide synthase", "side_text": "acyl-CoA"},
                {"from": "dhc", "to": "cer", "label": "desaturase (trans-Δ4)"},
                {"from": "cer", "to": "sm", "label": "SM synthase\n(PC → DAG)", "label_side": "right"},
                {"from": "cer", "to": "glccer", "label": "+ UDP-Glc", "label_side": "right"}, {"from": "cer", "to": "galcer", "label": "+ UDP-Gal", "label_side": "left"},
                {"from": "galcer", "to": "sulf", "label": "+ PAPS"}, {"from": "glccer", "to": "lac", "label": "+ UDP-Gal"},
                {"from": "lac", "to": "gang", "label": "+ CMP-NeuAc, UDP-sugars"}],
         tags=["sphingolipid_metabolism"],
         notes="在 ER 合成 ceramide，到 Golgi 加頭基。sphingomyelin 的 phosphocholine 來自 PC（不是 CDP-choline）。醣基由 UDP-sugar 提供，sialic acid 由 CMP-NeuAc 提供；sulfatide 的硫酸來自 PAPS。"),
    dict(id="sphingolipidoses", title="Sphingolipidoses: lysosomal degradation（鞘脂儲積症）", layout="free", cell_size=3.6,
         nodes=[n("gm1", "gm1", 0, 0, **NAME), n("gm2", "gm2", 0, 1.7, **NAME), n("gm3", "gm3", 0, 3.4, **NAME),
                n("gb3", "globotriaosylceramide", 4.6, 3.4, **NAME), n("lac", "lactosylceramide", 2.3, 5.1, **NAME),
                n("glc", "glucocerebroside", 2.3, 6.8, **NAME), n("sulf", "sulfatide", -4.6, 5.1, **NAME),
                n("gal", "galactocerebroside", -4.6, 6.8, **NAME), n("sm", "sphingomyelin", -1.2, 6.8, **NAME),
                n("cer", "ceramide", 0, 8.6, **NAME), n("sph", "sphingosine", 0, 10.3, **NAME),
                t("d1", "GM1 gangliosidosis", -3.4, 0.85, "drug"), t("d2", "Tay-Sachs (Hex A),\nSandhoff (Hex A + B)", -3.6, 2.55, "drug"),
                t("d3", "Fabry (X-linked)", 6.0, 4.6, "drug"), t("d4", "Gaucher", 5.0, 7.8, "drug"),
                t("d5", "Metachromatic\nleukodystrophy", -7.4, 5.95, "drug"), t("d6", "Krabbe", -6.4, 8.0, "drug"),
                t("d7", "Niemann-Pick A/B\n(sphingomyelinase)", -2.8, 9.9, "drug"), t("d8", "Farber", 2.6, 9.45, "drug")],
         edges=[{"from": "gm1", "to": "gm2", "step": "beta_galactosidase_gm1", "label": "β-galactosidase", "label_side": "right"},
                {"from": "gm2", "to": "gm3", "step": "hexosaminidase_a", "label": "β-hexosaminidase A", "label_side": "right"},
                {"from": "gm3", "to": "lac", "label": "neuraminidase", "label_side": "left"},
                {"from": "gb3", "to": "lac", "step": "alpha_galactosidase_a", "label": "α-galactosidase A", "label_side": "left"},
                {"from": "lac", "to": "glc", "label": "β-galactosidase", "label_side": "right"},
                {"from": "glc", "to": "cer", "step": "glucocerebrosidase", "label": "glucocerebrosidase", "label_side": "left"},
                {"from": "sulf", "to": "gal", "step": "arylsulfatase_a", "label": "arylsulfatase A", "label_side": "left"},
                {"from": "gal", "to": "cer", "step": "galactocerebrosidase", "label": "galactocerebrosidase", "label_side": "left"},
                {"from": "sm", "to": "cer", "step": "sphingomyelinase"},
                {"from": "cer", "to": "sph", "label": "ceramidase", "label_side": "right"},
                {"from": "d1", "to": "edge:gm1->gm2", "inhibit": True}, {"from": "d2", "to": "edge:gm2->gm3", "inhibit": True},
                {"from": "d3", "to": "edge:gb3->lac", "inhibit": True}, {"from": "d4", "to": "edge:glc->cer", "inhibit": True},
                {"from": "d5", "to": "edge:sulf->gal", "inhibit": True}, {"from": "d6", "to": "edge:gal->cer", "inhibit": True},
                {"from": "d7", "to": "edge:sm->cer", "inhibit": True}, {"from": "d8", "to": "edge:cer->sph", "inhibit": True}],
         tags=["sphingolipid_metabolism"],
         notes="溶小體由非還原端逐一切下糖。全部體染色體隱性，只有 Fabry 是 X-linked。櫻桃紅點：Tay-Sachs、Niemann-Pick（＋GM1、Sandhoff）。肝脾腫大：Gaucher、Niemann-Pick（Tay-Sachs 沒有）。Ashkenazi Jews：Tay-Sachs、Niemann-Pick A、Gaucher。脫髓鞘：Krabbe（globoid cell）、metachromatic leukodystrophy（sulfatide）。Gaucher 治療：酵素替代（imiglucerase）或受質減少（miglustat、eliglustat）。"),
]


def dump(path, d):
    path.write_text(yaml.safe_dump(d, allow_unicode=True, sort_keys=False, width=200), encoding="utf-8")


def with_parent(tags):
    out = list(tags)
    if any(t_ in KIDS for t_ in out) and CM not in out:
        out.insert(0, CM)
    return out


def add_cofactors(data):
    p = Path(data) / "cofactors.yaml"
    keys = list(COFACTORS) + list(MULTI)
    lines = [ln for ln in p.read_text(encoding="utf-8").splitlines()
             if not ln.startswith(MARK) and not any(ln.startswith(f"{k}:") for k in keys)]
    lines.append(MARK + "：GSH/GSSG（[16*] = glutathione 其餘部分）、多份 O₂、H₂O、H⁺（由 tools/gen_lipext.py 產生）")
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
        d["formula_check"] = CalcMolFormula(Chem.MolFromSmiles(m["smiles"]))
        for k in ("ring", "lipid", "free", "formula", "decorate"):
            if m.get(k):
                d[k] = m[k]
        d["tags"] = with_parent(m["tags"])
        d["notes_zh"] = m["notes"]
        d["status"] = "pending"
        d.pop("formula_check")
        dump(data / "molecules" / f"{m['id']}.yaml", d)
    for r in REACTIONS:
        d = {"id": r["id"], "title": r["title"], "substrates": r["substrates"], "products": r["products"]}
        for k in ("cofactors_in", "cofactors_out", "extra_in", "extra_out", "prosthetic"):
            if r.get(k):
                d[k] = r[k]
        d["enzyme"] = r["enzyme"]
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
