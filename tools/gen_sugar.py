"""One-off generator: carbohydrate structure library (molecules + overview pathways).

Fischer SMILES are computed from the drawing (fischer_smiles), cyclic SMILES come from PubChem
(sugar_lib) or are spliced from them with stereo-safe glycosyl fragments; every drawing is then
checked by `python -m sdb validate`.
"""
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, ".")
from fischer_smiles import fischer_smiles          # noqa: E402
from sugar_lib import (A_FRU, A_GAL, A_GLC, A_MAN, B_FRU, B_GAL, B_GLC, glycosyl, check)  # noqa: E402,F401

DATA = Path(sys.argv[1])
MOL, PW = DATA / "molecules", DATA / "pathways"


def q(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n") + '"'


def flow(d):
    """Tiny YAML flow writer for nested dicts / lists of str / numbers / bools."""
    if isinstance(d, dict):
        return "{" + ", ".join(f"{k}: {flow(v)}" for k, v in d.items()) + "}"
    if isinstance(d, (list, tuple)):
        return "[" + ", ".join(flow(v) for v in d) + "]"
    if isinstance(d, bool):
        return "true" if d else "false"
    if isinstance(d, (int, float)):
        return str(d)
    return q(d)


def write_mol(id_, name, smiles, tags, notes, aliases=(), label=None, chain=None, haworth=None, chair=None,
              free=None, default_view=None, decorate=None):
    L = [f"id: {id_}", f"name: {q(name)}"]
    if label:
        L.append(f"label: {q(label)}")
    L.append(f"aliases: [{', '.join(aliases)}]")
    L.append(f"smiles: {q(check(smiles))}")
    if default_view:
        L.append(f"default_view: {default_view}")
    if chain:
        L.append("chain:")
        L += [f"  - {flow(r)}" for r in chain]
    for key, spec in (("haworth", haworth), ("chair", chair)):
        if not spec:
            continue
        L.append(f"{key}:")
        for k, v in spec.items():
            if k == "subs":
                L.append("  subs:")
                L += [f"    {c}: {flow(s)}" for c, s in v.items()]
            elif k == "units":
                L.append("  units:")
                for u in v:
                    L.append(f"    - ring: {u['ring']}")
                    for kk, vv in u.items():
                        if kk == "ring":
                            continue
                        if kk == "subs":
                            L.append("      subs:")
                            L += [f"        {c}: {flow(s)}" for c, s in vv.items()]
                        else:
                            L.append(f"      {kk}: {flow(vv)}")
            else:
                L.append(f"  {k}: {flow(v)}")
    if free:
        L.append("free:")
        for k, v in free.items():
            if k == "atoms":
                L.append("  atoms:")
                L += [f"    {a}: {flow(t)}" for a, t in v.items()]
            else:
                L.append(f"  {k}: {flow(v)}")
    if decorate:
        L.append(f"decorate: {flow(decorate)}")
    L.append(f"tags: [{', '.join(tags)}]")
    L.append(f"notes_zh: {q(notes)}")
    L.append("status: pending")
    (MOL / f"{id_}.yaml").write_text("\n".join(L) + "\n", encoding="utf-8")


# ---------------------------------------------------------------- Fischer (open chain)

CHO = {"c": "C", "fork": {"l": "O", "lb": 2, "r": "H"}}
CH2OH = {"c": "CH_2OH"}
KETO = {"c": "C", "r": "O", "rb": 2}


def row(side):
    return {"c": "C", "l": "H", "r": "OH"} if side == "R" else {"c": "C", "l": "HO", "r": "H"}


def aldose_chain(sides):
    return [CHO] + [row(s) for s in sides] + [CH2OH]


def ketose_chain(sides):
    return [CH2OH, KETO] + [row(s) for s in sides] + [CH2OH]


MONO = ["carbohydrate", "monosaccharide"]
ALDOSES = [  # id, name, sides (C2 → last stereocentre), carbons, note
    ("d_glyceraldehyde", "D-Glyceraldehyde", "R", "醛糖（aldotriose）。D/L 的參考分子：OH 在右 = D。"),
    ("d_erythrose", "D-Erythrose", "RR", "醛四碳糖。erythrose 4-phosphate 是 PPP 與芳香族胺基酸合成（shikimate，植物／細菌）的原料。"),
    ("d_threose", "D-Threose", "LR", "醛四碳糖，erythrose 的 C2 epimer。"),
    ("d_ribose", "D-Ribose", "RRR", "醛五碳糖。RNA、ATP、NAD⁺、FAD、CoA 的五碳糖；由 PPP 產生（ribose 5-phosphate）。"),
    ("d_arabinose", "D-Arabinose", "LRR", "ribose 的 C2 epimer。"),
    ("d_xylose", "D-Xylose", "RLR", "木糖；蛋白聚醣（proteoglycan）連接區的第一個糖。D-xylose 吸收試驗：檢查小腸吸收功能。"),
    ("d_lyxose", "D-Lyxose", "LLR", "醛五碳糖。"),
    ("d_allose", "D-Allose", "RRRR", "醛六碳糖（四個 OH 都在右）。"),
    ("d_altrose", "D-Altrose", "LRRR", "allose 的 C2 epimer。"),
    ("d_mannose", "D-Mannose", "LLRR", "glucose 的 C2 epimer。醣蛋白 N-連結寡醣的主要成分；mannose 6-phosphate 是溶酶體酵素的標記（I-cell disease 缺此標記）。"),
    ("d_gulose", "D-Gulose", "RRLR", "醛六碳糖。"),
    ("d_idose", "D-Idose", "LRLR", "L-iduronate（idose 衍生物）是 heparan sulfate、dermatan sulfate 的成分。"),
    ("d_galactose", "D-Galactose", "RLLR", "glucose 的 C4 epimer。乳糖的組成；galactosemia（GALT 缺乏）→ galactitol 堆積 → 白內障。"),
    ("d_talose", "D-Talose", "LLLR", "galactose 的 C2 epimer。"),
]
KETOSES = [
    ("d_erythrulose", "D-Erythrulose", "R", "酮四碳糖。"),
    ("d_ribulose", "D-Ribulose", "RR", "酮五碳糖。ribulose 5-phosphate 是 PPP 氧化相的產物；RuBisCO 的受質為 ribulose 1,5-bisphosphate。"),
    ("d_xylulose", "D-Xylulose", "LR", "ribulose 的 C3 epimer；xylulose 5-phosphate 是 transketolase 的 2C 供體。"),
    ("d_psicose", "D-Psicose", "RRR", "酮六碳糖（allulose）。"),
    ("d_fructose", "D-Fructose", "LRR", "果糖。與 glucose 同分子式（C₆H₁₂O₆）的酮糖；肝中經 fructokinase → F1P → aldolase B。遺傳性果糖不耐（aldolase B 缺乏）。"),
    ("d_sorbose", "D-Sorbose", "RLR", "酮六碳糖。"),
    ("d_tagatose", "D-Tagatose", "LLR", "酮六碳糖（galactose 的酮糖對應物）。"),
]

for id_, name, sides, note in ALDOSES:
    n = len(sides) + 2
    sub = {3: "aldotriose", 4: "aldotetrose", 5: "aldopentose", 6: "aldohexose"}[n]
    write_mol(id_, name, fischer_smiles("aldose", sides), MONO + ["aldose"], f"{note}（{sub}）",
              aliases=[name.split("-", 1)[1]], chain=aldose_chain(sides))
write_mol("l_glyceraldehyde", "L-Glyceraldehyde", fischer_smiles("aldose", "L"), MONO + ["aldose"],
          "D-glyceraldehyde 的鏡像異構物（enantiomer）：OH 在左 = L。", aliases=["Glyceraldehyde"],
          chain=aldose_chain("L"))
write_mol("l_glucose", "L-Glucose", fischer_smiles("aldose", "LRLL"), MONO + ["aldose"],
          "D-glucose 的鏡像：所有立體中心都相反。判斷 D/L 只看離羰基最遠的立體中心（C5）。",
          aliases=["Glucose"], chain=aldose_chain("LRLL"))
write_mol("dihydroxyacetone", "Dihydroxyacetone", "OCC(=O)CO", MONO + ["ketose"],
          "最簡單的酮糖（ketotriose），沒有立體中心。DHAP 是 aldolase 的產物之一。",
          aliases=["DHA"], chain=[CH2OH, KETO, CH2OH])
for id_, name, sides, note in KETOSES:
    n = len(sides) + 3
    sub = {4: "ketotetrose", 5: "ketopentose", 6: "ketohexose"}[n]
    write_mol(id_, name, fischer_smiles("ketose", sides), MONO + ["ketose"], f"{note}（{sub}）",
              aliases=[name.split("-", 1)[1]], chain=ketose_chain(sides))

# alditols (polyol pathway)
write_mol("sorbitol", "D-Sorbitol", fischer_smiles("alditol", "RLRR"), ["carbohydrate", "sugar_derivative"],
          "D-glucitol。aldose reductase（NADPH）還原 glucose → sorbitol；高血糖時累積於水晶體、神經、腎 → 糖尿病白內障、神經病變。sorbitol dehydrogenase → fructose（polyol pathway）。",
          aliases=["glucitol", "山梨醇"], chain=[CH2OH] + [row(s) for s in "RLRR"] + [CH2OH])
write_mol("galactitol", "Galactitol", fischer_smiles("alditol", "RLLR"), ["carbohydrate", "sugar_derivative"],
          "dulcitol。galactosemia 時 galactose 經 aldose reductase 還原成 galactitol → 白內障。分子有對稱面（meso），無光學活性。",
          aliases=["dulcitol", "半乳糖醇"], chain=[CH2OH] + [row(s) for s in "RLLR"] + [CH2OH])

# ---------------------------------------------------------------- Haworth / chair (cyclic)

GLC_B = {1: {"up": "OH", "down": "H"}, 2: {"up": "H", "down": "OH"}, 3: {"up": "OH", "down": "H"},
         4: {"up": "H", "down": "HO"}, 5: {"up": "CH_2OH", "down": "H"}}


def subs(base, **over):
    d = {k: dict(v) for k, v in base.items()}
    for k, v in over.items():
        d[int(k[1:])] = v
    return d


GLC_A = subs(GLC_B, c1={"up": "H", "down": "OH"})
GAL_B = subs(GLC_B, c4={"up": "HO", "down": "H"})
GAL_A = subs(GAL_B, c1={"up": "H", "down": "OH"})
MAN_A = subs(GLC_A, c2={"up": "HO", "down": "H"})      # inner groups point into the ring
CYC = ["carbohydrate", "monosaccharide"]
PYR = {"ring": "pyranose", "numbers": True}

CH_B = {1: {"ax": "H", "eq": "OH"}, 2: {"ax": "H", "eq": "OH"}, 3: {"ax": "H", "eq": "HO"},
        4: {"ax": "H", "eq": "HO"}, 5: {"ax": "H", "eq": "CH_2OH"}}
CH_A = dict(CH_B)
CH_A[1] = {"ax": "OH", "eq": "H"}

write_mol("b_glucopyranose", "β-D-Glucopyranose", B_GLC, CYC + ["aldose"],
          "β：C1-OH 與 C5 的 CH₂OH 同側（Haworth 都在上）。椅式中所有大基團都在 equatorial → 最穩定，水溶液中約 64%。纖維素的單元。",
          aliases=["β-D-glucose", "β-glucose"], label="β-D-Glucose",
          haworth={**PYR, "subs": GLC_B}, chair={"ring": "pyranose", "subs": CH_B})
write_mol("a_glucopyranose", "α-D-Glucopyranose", A_GLC, CYC + ["aldose"],
          "α：C1-OH 在 Haworth 下方（與 CH₂OH 反側），椅式中為 axial。水溶液中約 36%。肝醣、澱粉的單元（α1→4、α1→6）。",
          aliases=["α-D-glucose", "α-glucose"], label="α-D-Glucose",
          haworth={**PYR, "subs": GLC_A}, chair={"ring": "pyranose", "subs": CH_A})
write_mol("b_galactopyranose", "β-D-Galactopyranose", B_GAL, CYC + ["aldose"],
          "glucose 的 C4 epimer：Haworth 中 C4-OH 在上。乳糖、醣脂（cerebroside）、血型抗原的成分。",
          aliases=["β-D-galactose"], label="β-D-Galactose", haworth={**PYR, "subs": GAL_B})
write_mol("a_mannopyranose", "α-D-Mannopyranose", A_MAN, CYC + ["aldose"],
          "glucose 的 C2 epimer：Haworth 中 C2-OH 在上。N-連結醣蛋白核心（Man₃GlcNAc₂）。",
          aliases=["α-D-mannose"], label="α-D-Mannose", haworth={**PYR, "subs": MAN_A})
FRU_B = {2: {"up": "OH", "down": "CH_2OH"}, 3: {"up": "HO", "down": "H"}, 4: {"up": "H", "down": "OH"},
         5: {"up": "HOCH_2", "down": "H"}}
FRU_A = dict(FRU_B)
FRU_A[2] = {"up": "CH_2OH", "down": "OH"}
FUR = {"ring": "furanose", "numbers": True, "start": 2}
write_mol("b_fructofuranose", "β-D-Fructofuranose", B_FRU, CYC + ["ketose"],
          "果糖以五員環（furanose）存在於蔗糖中；C2 是酮糖的 anomeric carbon。β：C2-OH 與 C6 的 CH₂OH 同側。",
          aliases=["β-D-fructose"], label="β-D-Fructofuranose", haworth={**FUR, "subs": FRU_B})
write_mol("a_fructofuranose", "α-D-Fructofuranose", A_FRU, CYC + ["ketose"],
          "α：C2-OH 在下（與 C6 的 CH₂OH 反側）。",
          aliases=["α-D-fructose"], label="α-D-Fructofuranose", haworth={**FUR, "subs": FRU_A})
RIB_B = {1: {"up": "OH", "down": "H"}, 2: {"up": "H", "down": "OH"}, 3: {"up": "H", "down": "OH"},
         4: {"up": "HOCH_2", "down": "H"}}
write_mol("b_ribofuranose", "β-D-Ribofuranose", "C([C@@H]1[C@H]([C@H]([C@@H](O1)O)O)O)O", CYC + ["aldose"],
          "RNA 與核苷酸中的核糖（β-furanose）。",
          aliases=["β-D-ribose"], label="β-D-Ribose", haworth={"ring": "furanose", "numbers": True, "subs": RIB_B})
write_mol("b_deoxyribofuranose", "2-Deoxy-β-D-ribofuranose", "C1[C@@H]([C@H](O[C@H]1O)CO)O",
          CYC + ["aldose", "sugar_derivative"],
          "DNA 的去氧核糖：C2 沒有 OH（只有兩個 H）。由 ribonucleotide reductase 在 NDP 層級還原產生。",
          aliases=["β-D-deoxyribose", "deoxyribose"], label="2-Deoxy-β-D-ribose",
          haworth={"ring": "furanose", "numbers": True, "subs": subs(RIB_B, c2={"up": "H", "down": "H"})})

# derivatives
DER = ["carbohydrate", "sugar_derivative"]
GLCA = "[O-]C(=O)[C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O"
write_mol("b_glucuronate", "β-D-Glucuronate", GLCA, DER,
          "C6 氧化成羧酸（COO⁻）。UDP-glucuronate 用於膽紅素、藥物的 glucuronidation（肝，增加水溶性）；也是玻尿酸、heparin 等 GAG 的組成。",
          aliases=["glucuronic acid", "GlcA", "葡萄糖醛酸"], label="β-D-Glucuronate",
          haworth={**PYR, "subs": subs(GLC_B, c5={"up": "COO^-", "down": "H"})})
GLCN = B_GLC.replace("[C@@H](O1)O)O)", "[C@@H](O1)O)[NH3+])")
write_mol("b_glucosamine", "β-D-Glucosamine", GLCN, DER,
          "胺基糖：C2-OH 換成胺基（生理 pH 為 –NH₃⁺）。由 fructose 6-P + Gln 合成（GFAT），是 GAG、幾丁質的前驅物。",
          aliases=["GlcN", "葡萄糖胺"], label="β-D-Glucosamine",
          haworth={**PYR, "subs": subs(GLC_B, c2={"up": "H", "down": "NH_3^+"})})
GLCNAC = B_GLC.replace("[C@@H](O1)O)O)", "[C@@H](O1)O)NC(C)=O)")
write_mol("b_glcnac", "N-Acetyl-β-D-glucosamine", GLCNAC, DER,
          "GlcNAc：C2 為乙醯胺基（NHAc = –NH–CO–CH₃）。玻尿酸、幾丁質、N-連結醣蛋白核心的成分；O-GlcNAc 修飾（營養感測）。",
          aliases=["GlcNAc", "N-乙醯葡萄糖胺"], label="N-Acetyl-β-D-\nglucosamine",
          haworth={**PYR, "subs": subs(GLC_B, c2={"up": "H", "down": "NHAc"})})

# ---------------------------------------------------------------- disaccharides


def at_c4(reducing, glycosyl_frag):
    """Replace the C4 hydroxyl of a PubChem-style glucose SMILES (C6 first, …C1(O1)O)O)O)O)O)
    with a glycosidic O carrying glycosyl_frag; the O keeps its place, so stereo tags hold."""
    for c1 in ("[C@@H](O1)O", "[C@H](O1)O"):
        tail = c1 + ")O)O)O)O"                 # C1-OH, then the OH of C2, C3, C4 and C6
        if reducing.endswith(tail):
            return reducing[:-2] + glycosyl_frag + ")O"
    raise ValueError(reducing)


DI = ["carbohydrate", "disaccharide"]
G1_A = subs(GLC_A, c1={"up": "H", "down": "link"})
G1_B = subs(GLC_B, c1={"up": "link", "down": "H"})
GAL1_B = subs(GAL_B, c1={"up": "link", "down": "H"})
RED_A4 = subs(GLC_B, c4={"up": "H", "down": "link"})          # reducing glucose, normal orientation
# reducing glucose turned over (flip): every up ↔ down
RED_FLIP = {1: {"up": "H", "down": "OH"}, 2: {"up": "OH", "down": "H"}, 3: {"up": "H", "down": "OH"},
            4: {"up": "link", "down": "H"}, 5: {"up": "H", "down": "CH_2OH"}}

write_mol("maltose", "Maltose", at_c4(B_GLC, glycosyl(A_GLC)), DI,
          "Glc(α1→4)Glc。澱粉、肝醣經 α-amylase 分解的產物；還原糖（右側 glucose 的 C1 自由）。小腸 maltase 分解。",
          aliases=["麥芽糖"], haworth={"units": [
              {"ring": "pyranose", "subs": G1_A},
              {"ring": "pyranose", "subs": RED_A4, "at": [4.9, 0]}],
              "links": [{"from": [0, 1, "down"], "to": [1, 4, "down"]}]})
write_mol("lactose", "Lactose", at_c4(B_GLC, glycosyl(B_GAL)), DI,
          "Gal(β1→4)Glc。乳糖；還原糖。小腸 lactase（β-galactosidase）分解；成人 lactase 不足 → 乳糖不耐（腹脹、滲透性腹瀉）。右側 glucose 依慣例翻面畫，讓 β 鍵兩端都朝上。",
          aliases=["乳糖"], haworth={"units": [
              {"ring": "pyranose", "subs": GAL1_B},
              {"ring": "pyranose", "subs": RED_FLIP, "at": [4.9, -0.2], "flip": True}],
              "links": [{"from": [0, 1, "up"], "to": [1, 4, "up"]}]})
write_mol("cellobiose", "Cellobiose", at_c4(B_GLC, glycosyl(B_GLC)), DI,
          "Glc(β1→4)Glc。纖維素的重複單元；人類沒有 cellulase，無法消化 β1→4 鍵（膳食纖維）。右側 glucose 翻面畫。",
          aliases=["纖維二糖"], haworth={"units": [
              {"ring": "pyranose", "subs": G1_B},
              {"ring": "pyranose", "subs": RED_FLIP, "at": [4.9, -0.2], "flip": True}],
              "links": [{"from": [0, 1, "up"], "to": [1, 4, "up"]}]})
SUC = A_GLC.replace("[C@H](O1)O)", "[C@H](O1)O" + glycosyl(B_FRU) + ")")
FRU_IN_SUC = {2: {"up": "HOCH_2", "down": "link"}, 3: {"up": "H", "down": "OH"}, 4: {"up": "OH", "down": "H"},
              5: {"up": "H", "down": "CH_2OH"}}
write_mol("sucrose", "Sucrose", SUC, DI,
          "Glc(α1↔β2)Fru。蔗糖；兩個 anomeric carbon 互相連結 → 非還原糖。小腸 sucrase 分解。右側 fructose 依 Lehninger 慣例在平面內轉 180°（C2 在左）。",
          aliases=["蔗糖", "table sugar"], haworth={"units": [
              {"ring": "pyranose", "subs": G1_A},
              {"ring": "furanose", "start": 2, "mirror": True, "subs": FRU_IN_SUC, "at": [4.7, 0.1]}],
              "links": [{"from": [0, 1, "down"], "to": [1, 2, "down"]}]})

# ---------------------------------------------------------------- polysaccharide cartoons (SNFG)


def chain_atoms(prefix, x0, y, n, dx=0.62):
    return {f"{prefix}{i}": ["@glc", round(x0 + i * dx, 3), y] for i in range(n)}


def glycogen():
    """Main chain left → right (reducing end at glycogenin), α1→6 branches above and below,
    each branch starting one step to the right of its branch point; no two chains overlap."""
    D = 0.62
    atoms, bonds = {}, []

    def chain(prefix, x0, y, n):
        atoms.update(chain_atoms(prefix, x0, y, n, D))
        bonds.extend([f"{prefix}{i}", f"{prefix}{i + 1}"] for i in range(n - 1))

    chain("m", 0, 0, 13)
    branches = [  # (name, parent residue, length, y)
        ("a", "m2", 4, -1.0), ("b", "m7", 4, -1.0), ("c", "m4", 5, 1.0), ("d", "m10", 3, 1.0),
        ("e", "a2", 3, -2.0), ("f", "c1", 3, 2.0)]
    for name, parent, n, y in branches:
        px = atoms[parent][1]
        chain(f"{name}", px + D, y, n)
        bonds.append([parent, f"{name}0"])
    notes = [["Glycogenin", -1.55, 0, "black"], ["reducing end", -1.55, 0.62, "gray"],
             ["α1→4", 8.3, 0, "gray"], ["α1→6\nbranch", 0.55, -1.0, "magenta"],
             ["nonreducing ends", 7.5, -2.0, "gray"]]
    return {"schematic": True, "atoms": atoms, "bonds": bonds, "notes": notes}


def linear(n, notes):
    atoms = chain_atoms("m", 0, 0, n)
    return {"schematic": True, "atoms": atoms, "bonds": [[f"m{i}", f"m{i + 1}"] for i in range(n - 1)],
            "notes": notes}


POLY = ["carbohydrate", "polysaccharide"]
write_mol("glycogen", "Glycogen", A_GLC, POLY,
          "示意圖（SNFG：藍色圓 = glucose）。α1→4 直鏈，每 8–12 個殘基一個 α1→6 分支；分支多 → 非還原端多，合成（glycogen synthase + branching enzyme）與分解（glycogen phosphorylase + debranching enzyme）都在非還原端進行。核心為 glycogenin（Tyr-OH 接第一個 glucose）。SMILES 僅以 α-glucose 代表單元。",
          aliases=["肝醣"], free=glycogen())
write_mol("amylose", "Amylose", A_GLC, POLY,
          "示意圖。澱粉的直鏈部分：α1→4，形成螺旋（遇碘呈藍色）。另一成分 amylopectin 結構似肝醣，但每 24–30 個殘基才一個分支。",
          aliases=["直鏈澱粉"], free=linear(9, [["α1→4 (helix)", 2.5, -0.75, "gray"], ["reducing end", 0, 0.75, "gray"]]))
write_mol("cellulose", "Cellulose", B_GLC, POLY,
          "示意圖。β1→4 直鏈（重複單元 cellobiose），鏈間大量氫鍵 → 纖維。人類缺 cellulase → 膳食纖維。",
          aliases=["纖維素"], free=linear(9, [["β1→4 (straight chains)", 2.5, -0.75, "gray"], ["reducing end", 0, 0.75, "gray"]]))

# ---------------------------------------------------------------- overview sheets


def write_pw(id_, title, body, tags, notes):
    L = [f"id: {id_}", f"title: {q(title)}"] + body + [f"tags: [{', '.join(tags)}]", f"notes_zh: {q(notes)}",
                                                     "status: pending"]
    (PW / f"{id_}.yaml").write_text("\n".join(L) + "\n", encoding="utf-8")


def tree(ids_rows, edges, newrow, cell_size, extra_nodes=()):
    """Fischer family tree: rows of molecule ids, node tops aligned (node_anchor: top)."""
    body = ["layout: free", f"cell_size: {cell_size}", "node_anchor: top", "nodes:"]
    for (y, xs) in ids_rows:
        for mid, x in xs:
            body.append(f"  - {{id: {mid}, ref: {mid}, pos: [{x}, {y}], decorate: {{colors: {{{newrow}: magenta}}}}}}")
    body += [f"  - {n}" for n in extra_nodes]
    body.append("edges:")
    body += [f"  - {{from: {a}, to: {b}, head: false, bend: 0}}" for a, b in edges]
    return body


H = 1.9                                   # hexose spacing in cells; parents sit midway
ALD_ROWS = [(0, [("d_glyceraldehyde", 3.5 * H)]),
            (2.6, [("d_erythrose", 1.5 * H), ("d_threose", 5.5 * H)]),
            (5.8, [("d_ribose", 0.5 * H), ("d_arabinose", 2.5 * H), ("d_xylose", 4.5 * H), ("d_lyxose", 6.5 * H)]),
            (9.6, [("d_allose", 0), ("d_altrose", H), ("glucose", 2 * H), ("d_mannose", 3 * H),
                   ("d_gulose", 4 * H), ("d_idose", 5 * H), ("d_galactose", 6 * H), ("d_talose", 7 * H)])]
ALD_EDGES = [("d_glyceraldehyde", "d_erythrose"), ("d_glyceraldehyde", "d_threose"),
             ("d_erythrose", "d_ribose"), ("d_erythrose", "d_arabinose"), ("d_threose", "d_xylose"),
             ("d_threose", "d_lyxose"), ("d_ribose", "d_allose"), ("d_ribose", "d_altrose"),
             ("d_arabinose", "glucose"), ("d_arabinose", "d_mannose"), ("d_xylose", "d_gulose"),
             ("d_xylose", "d_idose"), ("d_lyxose", "d_galactose"), ("d_lyxose", "d_talose")]
write_pw("d_aldose_family", "D-Aldoses（D-醛糖家族樹）", tree(ALD_ROWS, ALD_EDGES, 1, 3.6), MONO + ["aldose"],
         "每往下一層，在羰基下方多一個立體中心（洋紅），OH 左右各一 → 2ⁿ 個異構物。全部 D 型：離羰基最遠的立體中心 OH 在右。必背：glucose、mannose（C2 epimer）、galactose（C4 epimer）、ribose。")
KET_ROWS = [(0, [("dihydroxyacetone", 1.5 * H)]),
            (2.6, [("d_erythrulose", 1.5 * H)]),
            (5.8, [("d_ribulose", 0.5 * H), ("d_xylulose", 2.5 * H)]),
            (9.6, [("d_psicose", 0), ("d_fructose", H), ("d_sorbose", 2 * H), ("d_tagatose", 3 * H)])]
KET_EDGES = [("dihydroxyacetone", "d_erythrulose"), ("d_erythrulose", "d_ribulose"),
             ("d_erythrulose", "d_xylulose"), ("d_ribulose", "d_psicose"), ("d_ribulose", "d_fructose"),
             ("d_xylulose", "d_sorbose"), ("d_xylulose", "d_tagatose")]
write_pw("d_ketose_family", "D-Ketoses（D-酮糖家族樹）", tree(KET_ROWS, KET_EDGES, 2, 3.6), MONO + ["ketose"],
         "酮基在 C2；新增的立體中心在 C3（洋紅）。必背：fructose、ribulose、xylulose、dihydroxyacetone。")

write_pw("epimers", "Epimers: Mannose · Glucose · Galactose", [
    "layout: free", "cell_size: 4.2", "node_anchor: top", "nodes:",
    "  - {id: man, ref: d_mannose, pos: [0, 0], decorate: {colors: {1: magenta}}}",
    "  - {id: glc, ref: glucose, pos: [3, 0], decorate: {colors: {1: magenta, 3: magenta}}}",
    "  - {id: gal, ref: d_galactose, pos: [6, 0], decorate: {colors: {3: magenta}}}",
    "edges:",
    "  - {from: man, to: glc, head: false, bend: 0, label: \"C-2 epimers\"}",
    "  - {from: glc, to: gal, head: false, bend: 0, label: \"C-4 epimers\"}"],
    MONO + ["aldose"],
    "epimer：只差一個立體中心。mannose 與 glucose 差在 C2；galactose 與 glucose 差在 C4。mannose 與 galactose 差兩個，不是 epimer。")
write_pw("d_l_sugars", "D- and L-sugars（D/L 構型）", [
    "layout: free", "cell_size: 4.2", "node_anchor: top", "nodes:",
    "  - {id: lga, ref: l_glyceraldehyde, pos: [0, 0], decorate: {colors: {1: magenta}}}",
    "  - {id: dga, ref: d_glyceraldehyde, pos: [2.4, 0], decorate: {colors: {1: magenta}}}",
    "  - {id: lgl, ref: l_glucose, pos: [4.8, 0], decorate: {colors: {4: magenta}}}",
    "  - {id: dgl, ref: glucose, pos: [7.2, 0], decorate: {colors: {4: magenta}}}",
    "edges:",
    "  - {from: lga, to: dga, head: false, bend: 0, label: \"mirror images\"}",
    "  - {from: lgl, to: dgl, head: false, bend: 0, label: \"mirror images\"}"],
    MONO,
    "D/L 由離羰基最遠的立體中心決定（洋紅）：OH 在右 = D。天然糖幾乎都是 D 型（胺基酸是 L 型）。L-glucose 是 D-glucose 的鏡像（每個中心都反過來），不是 epimer。")
write_pw("glucose_anomers", "Mutarotation of D-glucose（α/β 互變）", [
    "layout: free", "cell: 1.0", "nodes:",
    "  - {id: a, ref: a_glucopyranose, pos: [0, 0], decorate: {colors: {d1: magenta}}}",
    "  - {id: open, ref: glucose, pos: [1.25, 0], decorate: {colors: {0: magenta}}}",
    "  - {id: b, ref: b_glucopyranose, pos: [2.5, 0], decorate: {colors: {u1: magenta}}}",
    "  - {id: pa, text: \"≈ 36 %\", kind: plain, pos: [0, 0.45]}",
    "  - {id: pb, text: \"≈ 64 %\", kind: plain, pos: [2.5, 0.45]}",
    "edges:",
    "  - {from: a, to: open, reversible: true, bend: 0}",
    "  - {from: open, to: b, reversible: true, bend: 0}"],
    CYC,
    "C5-OH 攻擊 C1 醛基形成半縮醛（hemiacetal）→ 產生新的立體中心 C1（anomeric carbon）。α、β 互為 anomer，經開鏈互變（mutarotation）。β 較多：椅式中所有大基團 equatorial。")
write_pw("disaccharides", "Disaccharides（雙醣）", [
    "layout: free", "cell: 1.0", "nodes:",
    "  - {id: malt, ref: maltose, pos: [0, 0]}",
    "  - {id: lac, ref: lactose, pos: [1.05, 0]}",
    "  - {id: suc, ref: sucrose, pos: [0, 0.82]}",
    "  - {id: cel, ref: cellobiose, pos: [1.05, 0.82]}",
    "edges: []"],
    DI,
    "還原糖：maltose、lactose、cellobiose（右側單元的 anomeric C 自由）。非還原糖：sucrose（α1↔β2，兩個 anomeric C 都用掉）、trehalose。消化酵素：maltase、lactase、sucrase-isomaltase（都在小腸刷狀緣）。")
print("ok")
