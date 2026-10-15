"""Generator: lipid structure library (fatty acids, TAG, membrane lipids, steroids, eicosanoids).

Glycerolipids / sphingolipids / fatty acids: the `lipid:` spec builds both the drawing and the
SMILES (sdb.draw.lipid.build), so the YAML SMILES is written from the spec and the validator
re-checks it. Steroid and eicosanoid SMILES are PubChem isomeric SMILES (ionised where relevant);
their formulas are printed so they can be compared with PubChem.

    .venv/bin/python tools/gen_lipid.py data && .venv/bin/python -m sdb validate
"""
import sys
from pathlib import Path

sys.path.insert(0, ".")
from rdkit import Chem                                   # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula   # noqa: E402
from sdb.draw.lipid import build                         # noqa: E402

DATA = Path(sys.argv[1])
MOL, PW = DATA / "molecules", DATA / "pathways"


def q(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n") + '"'


def flow(d):
    if isinstance(d, dict):
        return "{" + ", ".join(f"{k}: {flow(v)}" for k, v in d.items()) + "}"
    if isinstance(d, (list, tuple)):
        return "[" + ", ".join(flow(v) for v in d) + "]"
    if isinstance(d, bool):
        return "true" if d else "false"
    if isinstance(d, (int, float)):
        return str(d)
    return q(d)


def write(id_, name, smiles, tags, notes, aliases=(), label=None, lipid=None, ring=None):
    assert Chem.MolFromSmiles(smiles) is not None, smiles
    L = [f"id: {id_}", f"name: {q(name)}"]
    if label:
        L.append(f"label: {q(label)}")
    L += [f"aliases: [{', '.join(aliases)}]", f"smiles: {q(smiles)}"]
    if lipid:
        L.append(f"lipid: {flow(lipid)}")
    if ring:
        L.append(f"ring: {flow(ring)}")
    L += [f"tags: [{', '.join(tags)}]", f"notes_zh: {q(notes)}", "status: pending"]
    (MOL / f"{id_}.yaml").write_text("\n".join(L) + "\n", encoding="utf-8")


def lipid(id_, name, spec, tags, notes, aliases=(), label=None):
    _f, smi, problems, _n = build(spec)
    assert not problems, (id_, problems)
    write(id_, name, smi, tags, notes, aliases, label, lipid=spec)


FA = ["lipid", "fatty_acid"]
# ---------------------------------------------------------------- fatty acids
for id_, name, n, cis, trans, sym, note, al in [
    ("palmitate", "Palmitate", 16, [], [], "16:0", "飽和脂肪酸。fatty acid synthase 的產物（7 次循環，消耗 7 ATP + 14 NADPH）；人體最常見的飽和脂肪酸之一。", ["palmitic acid", "棕櫚酸"]),
    ("stearate", "Stearate", 18, [], [], "18:0", "飽和脂肪酸。可被 Δ9 desaturase 轉成 oleate。", ["stearic acid", "硬脂酸"]),
    ("palmitoleate", "Palmitoleate", 16, [9], [], "16:1(Δ9)", "單元不飽和（ω-7）。", ["palmitoleic acid"]),
    ("oleate", "Oleate", 18, [9], [], "18:1(Δ9)", "單元不飽和、cis-Δ9（ω-9）。cis 雙鍵使鏈彎折 → 熔點降低、膜流動性增加。人體可自行合成（Δ9 desaturase）。", ["oleic acid", "油酸"]),
    ("elaidate", "Elaidate", 18, [], [9], "18:1(trans-Δ9)", "oleate 的 trans 異構物（反式脂肪）：鏈幾乎是直的，物性接近飽和脂肪；增加 LDL、心血管風險。", ["elaidic acid", "trans fat"]),
    ("linoleate", "Linoleate", 18, [9, 12], [], "18:2(Δ9,12)", "必需脂肪酸（ω-6）：人體缺 Δ12、Δ15 desaturase，不能在 C9 之後（往甲基端）加雙鍵。arachidonate 的前驅物。", ["linoleic acid", "LA", "亞麻油酸"]),
    ("alpha_linolenate", "α-Linolenate", 18, [9, 12, 15], [], "18:3(Δ9,12,15)", "必需脂肪酸（ω-3）。可轉成 EPA、DHA（效率低）。", ["ALA", "α-linolenic acid", "次亞麻油酸"]),
    ("arachidonate", "Arachidonate", 20, [5, 8, 11, 14], [], "20:4(Δ5,8,11,14)", "ω-6，由 linoleate 延長與去飽和而來。phospholipase A₂ 從膜磷脂 sn-2 釋出 → COX（prostaglandin、thromboxane）、LOX（leukotriene）→ eicosanoids。", ["arachidonic acid", "AA", "花生四烯酸"]),
    ("epa", "Eicosapentaenoate", 20, [5, 8, 11, 14, 17], [], "20:5(Δ5,8,11,14,17)", "EPA，ω-3（魚油）。", ["EPA", "eicosapentaenoic acid"]),
    ("dha", "Docosahexaenoate", 22, [4, 7, 10, 13, 16, 19], [], "22:6(Δ4,7,10,13,16,19)", "DHA，ω-3（魚油），視網膜與腦的膜脂大量含有。", ["DHA", "docosahexaenoic acid"]),
]:
    omega = n - max(cis + trans) if (cis or trans) else None
    extra = f"；ω-{omega}" if omega else ""
    lipid(id_, name, {"type": "fatty_acid", "chains": [{"n": n, "cis": cis, "trans": trans}], "numbers": True},
          FA, f"{sym}{extra}。{note}", al, label=f"{name}\n{sym}")

# ---------------------------------------------------------------- triacylglycerol / phospholipids
lipid("triacylglycerol", "Triacylglycerol", {"type": "tag", "numbers": False, "chains": [
    {"n": 16}, {"n": 18}, {"n": 18, "cis": [9]}]},
      ["lipid", "triacylglycerol"],
      "例：1-palmitoyl-2-stearoyl-3-oleoyl-glycerol（不飽和鏈畫在 sn-3，彎折才不會穿過其他鏈；天然 TAG 的 sn-2 多為不飽和）。三個脂肪酸以酯鍵接在 glycerol 上，中性、疏水 → 能量儲存（脂肪組織）。hormone-sensitive lipase／ATGL 水解；腸道由 pancreatic lipase 分解。",
      ["TAG", "triglyceride", "三酸甘油酯"], label="Triacylglycerol")
PL = ["lipid", "phospholipid"]
PO = [{"n": 16}, {"n": 18, "cis": [9]}]                    # palmitoyl / oleoyl (sn-1 / sn-2)
for id_, name, head, note, al in [
    ("phosphatidate", "Phosphatidate", "phosphate", "PA：最簡單的甘油磷脂，TAG 與其他磷脂合成的分歧點（→ DAG → TAG、PC、PE；→ CDP-DAG → PI、PG、cardiolipin）。", ["PA", "phosphatidic acid"]),
    ("phosphatidylcholine", "Phosphatidylcholine", "choline", "PC（lecithin）：膜外葉主要磷脂；肺表面活性劑主成分為 dipalmitoyl-PC（L/S ratio 判斷胎兒肺成熟度）。", ["PC", "lecithin", "卵磷脂"]),
    ("phosphatidylethanolamine", "Phosphatidylethanolamine", "ethanolamine", "PE：膜內葉為主；可甲基化成 PC（肝，PEMT）。", ["PE", "cephalin"]),
    ("phosphatidylserine", "Phosphatidylserine", "serine", "PS：帶負電，正常位於膜內葉；細胞凋亡時翻到外葉 → 巨噬細胞辨識（annexin V 染色）；血小板活化時外翻促進凝血。", ["PS"]),
    ("phosphatidylinositol", "Phosphatidylinositol", "inositol", "PI：磷酸化成 PIP₂ → phospholipase C 切成 IP₃（釋放 Ca²⁺）+ DAG（活化 PKC）；GPI anchor。", ["PI"]),
    ("phosphatidylglycerol", "Phosphatidylglycerol", "glycerol", "PG：肺表面活性劑成分；cardiolipin（粒線體內膜）的前驅物。", ["PG"]),
]:
    lipid(id_, name, {"type": "phospholipid", "chains": PO, "head": head}, PL,
          f"{note}（圖例：sn-1 palmitoyl、sn-2 oleoyl）", al)
SL = ["lipid", "sphingolipid"]
for id_, name, head, note, al in [
    ("ceramide", "Ceramide", "oh", "sphingosine + 脂肪酸（醯胺鍵，接在 C2 的 N）。所有 sphingolipid 的母體；也是凋亡訊息分子。", ["N-acylsphingosine"]),
    ("sphingomyelin", "Sphingomyelin", "phosphocholine", "ceramide + phosphocholine（C1）。髓鞘含量高。Niemann-Pick disease：sphingomyelinase 缺乏。", ["SM"]),
    ("galactocerebroside", "Galactocerebroside", "galactose", "ceramide + 一個 β-galactose（C1）。神經髓鞘。Krabbe disease：galactocerebrosidase 缺乏。（glucocerebroside 累積 → Gaucher disease）", ["GalCer", "cerebroside"]),
]:
    lipid(id_, name, {"type": "sphingolipid", "chains": [{"n": 18}], "head": head}, SL,
          f"{note}（圖例：N-stearoyl；sphingosine 有 trans-Δ4 雙鍵）", al)

# ---------------------------------------------------------------- steroids & eicosanoids (PubChem)
ST = ["lipid", "steroid"]
STEROIDS = [
    ("cholesterol", "Cholesterol", "CC(C)CCC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
     "C27H46O", "27 碳。A–D 四環固醇核，C3-OH（兩性：OH 極性、其餘疏水）、C5=C6 雙鍵、C17 側鏈。膜的流動性緩衝；固醇類荷爾蒙、膽酸、維生素 D 的前驅物。限速酵素 HMG-CoA reductase（statin）。", ["膽固醇"]),
    ("cortisol", "Cortisol", "C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2[C@H](C[C@]4([C@H]3CC[C@@]4(C(=O)CO)O)C)O",
     "C21H30O5", "glucocorticoid（腎上腺皮質束狀帶）。21-hydroxylase 缺乏 → cortisol、aldosterone ↓ → CAH。", ["hydrocortisone", "皮質醇"]),
    ("aldosterone", "Aldosterone", "C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2[C@H](C[C@]4([C@H]3CC[C@@H]4C(=O)CO)C=O)O",
     "C21H28O5", "mineralocorticoid（腎上腺皮質球狀帶），C18 為醛基；作用於集尿管保鈉排鉀。", ["醛固酮"]),
    ("testosterone", "Testosterone", "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C",
     "C19H28O2", "androgen（19 碳）。5α-reductase → DHT；aromatase → estradiol。", ["睪固酮"]),
    ("estradiol", "Estradiol", "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O",
     "C18H24O2", "estrogen（18 碳），A 環為芳香環（aromatase 由 testosterone 轉來），沒有 C19 甲基。", ["E2", "雌二醇"]),
    ("progesterone", "Progesterone", "CC(=O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4=CC(=O)CC[C@]34C)C",
     "C21H30O2", "progestin（21 碳）；pregnenolone → progesterone，是 corticosteroid 的前驅物。", ["黃體素"]),
    ("cholate", "Cholate", "C[C@H](CCC(=O)[O-])[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)O)C",
     "C24H39O5-", "初級膽酸（3α,7α,12α-OH）。限速酵素 7α-hydroxylase；與 glycine／taurine 結合成膽鹽，乳化脂質。", ["cholic acid", "膽酸"]),
]
for id_, name, smi, formula, note, al in STEROIDS:
    got = CalcMolFormula(Chem.MolFromSmiles(smi))
    assert got.replace("-", "") == formula.replace("-", ""), (id_, got, formula)
    write(id_, name, smi, ST, note, al, ring={"letters": True, "orient": "steroid", "methyls": False})
EIC = ["lipid", "eicosanoid"]
for id_, name, smi, formula, note, al in [
    ("pge2", "Prostaglandin E2", "CCCCC[C@@H](/C=C/[C@H]1[C@@H](CC(=O)[C@@H]1C/C=C\\CCCC(=O)[O-])O)O", "C20H31O5-",
     "PGE₂：arachidonate 經 COX（cyclooxygenase）→ PGH₂ → PGE₂。發燒、疼痛、發炎；保護胃黏膜。NSAID（aspirin 不可逆乙醯化 COX）抑制。", ["PGE₂", "dinoprostone"]),
]:
    got = CalcMolFormula(Chem.MolFromSmiles(smi))
    assert got.replace("-", "") == formula.replace("-", ""), (id_, got, formula)
    write(id_, name, smi, EIC, note, al, label="PGE_2", ring={"rotate": 0})


# ---------------------------------------------------------------- overview sheets

def write_pw(id_, title, body, tags, notes):
    L = [f"id: {id_}", f"title: {q(title)}"] + body + [f"tags: [{', '.join(tags)}]", f"notes_zh: {q(notes)}",
                                                     "status: pending"]
    (PW / f"{id_}.yaml").write_text("\n".join(L) + "\n", encoding="utf-8")


write_pw("fatty_acids", "Fatty acids（脂肪酸：飽和、cis／trans、ω-3／ω-6）", [
    "layout: free", "cell_size: 4.6", "node_anchor: top", "text: false", "nodes:",
    "  - {id: h1, text: \"Saturated\", kind: heading, pos: [-0.2, -0.75]}",
    "  - {id: pal, ref: palmitate, pos: [2.3, 0]}",
    "  - {id: ste, ref: stearate, pos: [7.4, 0]}",
    "  - {id: h2, text: \"cis / trans\", kind: heading, pos: [-0.2, 1.65]}",
    "  - {id: ole, ref: oleate, pos: [2.3, 2.4]}",
    "  - {id: ela, ref: elaidate, pos: [7.4, 2.4]}",
    "  - {id: h3, text: \"Essential (ω-6, ω-3)\", kind: heading, pos: [-0.2, 5.25]}",
    "  - {id: lin, ref: linoleate, pos: [2.3, 6.0]}",
    "  - {id: aln, ref: alpha_linolenate, pos: [7.4, 6.0]}",
    "  - {id: h4, text: \"Long-chain polyunsaturated\", kind: heading, pos: [-0.2, 8.9]}",
    "  - {id: ara, ref: arachidonate, pos: [1.6, 9.65]}",
    "  - {id: epa, ref: epa, pos: [5.0, 9.65]}",
    "  - {id: dha, ref: dha, pos: [8.6, 9.65]}",
    "edges: []"],
    FA, "Δ 編號從羧基 C1 數；ω 編號從甲基端數（ω-3：最後一個雙鍵距甲基端 3 個碳）。cis 雙鍵讓鏈彎折、熔點降低；trans 接近直鏈。人類不能合成 linoleate（ω-6）與 α-linolenate（ω-3）→ 必需脂肪酸。")
GP = [("pa", "phosphatidate"), ("pc", "phosphatidylcholine"), ("pe", "phosphatidylethanolamine"),
      ("ps", "phosphatidylserine"), ("pi", "phosphatidylinositol"), ("pg", "phosphatidylglycerol")]
SPH = [("cer", "ceramide"), ("sm", "sphingomyelin"), ("gc", "galactocerebroside")]
body = ["layout: free", "cell_size: 5.0", "node_anchor: top", "text: false", "nodes:",
        "  - {id: h1, text: \"Glycerophospholipids (head group on the sn-3 phosphate)\", kind: heading, pos: [-0.3, -0.8]}"]
body += [f"  - {{id: {i}, ref: {r}, pos: [{1.4 + 3.2 * k}, 0]}}" for k, (i, r) in enumerate(GP)]
body += ["  - {id: h2, text: \"Sphingolipids (sphingosine backbone)\", kind: heading, pos: [-0.3, 7.0]}"]
body += [f"  - {{id: {i}, ref: {r}, pos: [{1.4 + 3.2 * k}, 7.8]}}" for k, (i, r) in enumerate(SPH)]
body += ["  - {id: h3, text: \"Neutral fat (storage, not in membranes)\", kind: heading, pos: [9.4, 7.0]}",
         "  - {id: tag, ref: triacylglycerol, pos: [11.8, 7.8]}", "edges: []"]
write_pw("membrane_lipids", "Membrane lipids（甘油磷脂與鞘脂）", body,
    ["lipid", "phospholipid", "sphingolipid"],
    "甘油磷脂：glycerol 的 sn-1、sn-2 接脂肪酸（sn-2 常為不飽和），sn-3 接磷酸＋頭基。鞘脂：sphingosine（C1 頭基、C2 醯胺接脂肪酸、C3-OH、trans-Δ4）。TAG 放在右下角作對照（中性、不在膜上）。")
write_pw("steroids", "Steroids（固醇類：膽固醇與固醇類荷爾蒙）", [
    "layout: free", "cell: 1.0", "text: false", "nodes:",
    "  - {id: cho, ref: cholesterol, pos: [0.5, 0]}",
    "  - {id: pro, ref: progesterone, pos: [0, 1.05]}",
    "  - {id: cor, ref: cortisol, pos: [1.05, 1.05]}",
    "  - {id: ald, ref: aldosterone, pos: [2.1, 1.05]}",
    "  - {id: tes, ref: testosterone, pos: [0, 2.1]}",
    "  - {id: est, ref: estradiol, pos: [1.05, 2.1]}",
    "  - {id: cha, ref: cholate, pos: [2.1, 2.1]}",
    "edges: []"],
    ST, "共同骨架：A–D 四環（cyclopentanoperhydrophenanthrene）。碳數：cholesterol 27 → progestin／corticoid 21 → androgen 19 → estrogen 18（A 環芳香化）；膽酸 24。立體化學（楔形鍵）目前未畫出。")
print("ok")
