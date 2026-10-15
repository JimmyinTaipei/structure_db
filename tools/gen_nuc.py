"""One-off generator for nucleotide molecule YAMLs (keeps the sugar stereo identical everywhere)."""
import re
import sys
from pathlib import Path
from rdkit import Chem

OUT = Path(sys.argv[1])
P1 = "COP(=O)([O-])[O-]"
P2 = "COP(=O)([O-])OP(=O)([O-])[O-]"
P3 = "COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-]"
OH5 = "CO"
LAB = {P1: "{P}OCH_2", P2: "{P}–{P}–OCH_2", P3: "{P}–{P}–{P}–OCH_2", OH5: "HOCH_2"}


def sugar(d, c5, deoxy=False):
    c2 = f"C{d}" if deoxy else f"[C@H]{d}O"
    return f"[C@@H]{d}O[C@H]({c5})[C@@H](O){c2}"


# base: (prefix written before C1' for purines | wrapper with {S} for pyrimidines, [*]-SMILES for drawing)
BASES = {
    "A": ("Nc1ncnc2c1ncn2{S}", "[*]n1cnc2[c:2](N)ncnc21", True),
    "G": ("O=c1[nH]c(N)nc2c1ncn2{S}", "[*]n1cnc2[c:2](=O)[nH]c(N)nc21", True),
    "I": ("O=c1[nH]cnc2c1ncn2{S}", "[*]n1cnc2[c:2](=O)[nH]cnc21", True),
    "X": ("O=c1[nH]c(=O)[nH]c2c1ncn2{S}", "[*]n1cnc2[c:2](=O)[nH]c(=O)[nH]c21", True),
    "U": ("O=c1ccn({S})c(=O)[nH]1", "[*]n1ccc(=O)[nH][c:2]1=O", False),
    "C": ("Nc1ccn({S})c(=O)n1", "[*]n1ccc(N)n[c:2]1=O", False),
    "T": ("Cc1cn({S})c(=O)[nH]c1=O", "[*]n1cc(C)c(=O)[nH][c:2]1=O", False),
    "O": ("O=c1cc(C(=O)[O-])n({S})c(=O)[nH]1", "[*]n1[c:2](C(=O)[O-])cc(=O)[nH]c1=O", False),
    "SAMP": ("[O-]C(=O)C[C@H](Nc1ncnc2c1ncn2{S})C(=O)[O-]", "[*]n1cnc2[c:2](N[C@@H](CC(=O)[O-])C(=O)[O-])ncnc21", True),
}


def nuc(base, c5, deoxy=False):
    tpl, draw, purine = BASES[base]
    d = "1" if purine else "2"
    smi = tpl.replace("{S}", sugar(d, c5, deoxy))
    Chem.CanonSmiles(smi)
    subs = {1: {"up": "base", "down": "H"},
            2: {"up": "H", "down": "H" if deoxy else "OH"},
            3: {"up": "H", "down": "OH"},
            4: {"up": LAB[c5], "down": "H"}}
    draw, maps = strip(draw)
    h = {"ring": "furanose", "subs": subs, "base": draw}
    if 2 in maps:
        h["base_right"] = maps[2]
    if base == "O":                       # C2=O would sit on the 5′ label
        h["base_len"] = 1.45
    return smi, h


def q(s):
    return '"' + s.replace('\\', '\\\\').replace('"', '\\"').replace('\n', '\\n') + '"'


def write(id_, name, smiles, view, tags, notes, aliases=(), label=None):
    lines = [f"id: {id_}", f"name: {q(name)}"]
    if label:
        lines.append(f"label: {q(label)}")
    lines.append(f"aliases: [{', '.join(aliases)}]")
    lines.append(f"smiles: {q(smiles)}")
    if "haworth" in view:
        h = view["haworth"]
        lines.append("haworth:")
        lines.append(f"  ring: {h['ring']}")
        lines.append("  subs:")
        for k, v in h["subs"].items():
            lines.append(f"    {k}: {{up: {q(v['up'])}, down: {q(v['down'])}}}")
        if h.get("base"):
            lines.append(f"  base: {q(h['base'])}")
            if h.get("base_flip"):
                lines.append("  base_flip: true")
            if "base_len" in h:
                lines.append(f"  base_len: {h['base_len']}")
            if "base_right" in h:
                lines.append(f"  base_right: {h['base_right']}     # 這個原子放右邊（嘧啶 C2；OMP 為 C6，讓 COO⁻ 避開糖）")
        lines.append("  numbers: true")
    else:
        r = view["ring"]
        lines.append("ring: {" + ", ".join(f"{k}: {v}" for k, v in r.items()) + "}")
    lines.append(f"tags: [{', '.join(tags)}]")
    lines.append(f"notes_zh: {q(notes)}")
    lines.append("status: pending")
    (OUT / f"{id_}.yaml").write_text("\n".join(lines) + "\n", encoding="utf-8")


def strip(smiles_mapped):
    """'…[nH:9]…' → (smiles without atom maps, {map number: atom index})."""
    m = Chem.MolFromSmiles(smiles_mapped)
    maps = {a.GetAtomMapNum(): a.GetIdx() for a in m.GetAtoms() if a.GetAtomMapNum()}
    # drop the map text only, so atom order (= indices) stays exactly as written
    smi = re.sub(r":\d+\]", "]", smiles_mapped)
    smi = re.sub(r"\[([cnoCNO])\]", r"\1", smi)        # [c] → c: keep implicit H
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    assert Chem.CanonSmiles(smi) == Chem.MolToSmiles(m), smi
    return smi, maps


def free_base(smiles_mapped):
    smi, maps = strip(smiles_mapped)
    r = {"down": maps.get(9, maps.get(1))}
    if 2 in maps:
        r["right"] = maps[2]
    return smi, r


def base_full(mapped, purine):
    """Atom maps = ring numbers. Vertical (N9 / N1 straight down), C6 / C2 on the right,
    skeletal ring with ring N–H drawn as a bond, plus a mirrored twin (textbook base cards)."""
    smi, maps = strip(mapped)
    if purine:      # fused C4–C5 bond vertical, six-membered ring pointing up, N9–H straight down
        r = {"axis": f"[{maps[4]}, {maps[5]}]", "right": maps[6], "h_down": maps[9]}
    else:           # N1 at the bottom, C4 on top
        r = {"down": maps[1], "right": maps[2]}
    r.update({"nh_bond": "true", "mirror": "true"})
    return smi, r


PUR, PYR = ["nucleotide", "purine"], ["nucleotide", "pyrimidine"]
FULL = {
    "adenine": ("N[c:6]1[n:1][cH:2][n:3][c:4]2[nH:9][cH:8][n:7][c:5]12", True),
    "guanine": ("N[c:2]1[n:3][c:4]2[nH:9][cH:8][n:7][c:5]2[c:6](=O)[nH:1]1", True),
    "hypoxanthine": ("O=[c:6]1[nH:1][cH:2][n:3][c:4]2[nH:9][cH:8][n:7][c:5]12", True),
    "xanthine": ("O=[c:6]1[nH:1][c:2](=O)[nH:3][c:4]2[nH:9][cH:8][n:7][c:5]12", True),
    "uric_acid": ("O=[c:6]1[nH:1][c:2](=O)[nH:3][c:4]2[nH:9][c:8](=O)[nH:7][c:5]12", True),
    "cytosine": ("N[c:4]1[cH:5][cH:6][nH:1][c:2](=O)[n:3]1", False),
    "uracil": ("O=[c:4]1[cH:5][cH:6][nH:1][c:2](=O)[nH:3]1", False),
    "thymine": ("C[c:5]1[cH:6][nH:1][c:2](=O)[nH:3][c:4]1=O", False),
}
# ---------------------------------------------------------------- free bases
for id_, name, mapped, tags, flip, notes, al in [
    ("adenine", "Adenine", "N[c:2]1ncnc2[nH:9]cnc12", PUR, True, "嘌呤鹼基（6-胺基嘌呤）。A 與 T/U 以 2 個氫鍵配對。", ["A", "腺嘌呤"]),
    ("guanine", "Guanine", "Nc1nc2[nH:9]cnc2[c:2](=O)[nH]1", PUR, True, "嘌呤鹼基（2-胺基-6-氧嘌呤）。G≡C 3 個氫鍵。分解：guanine deaminase（guanase）→ xanthine。", ["G", "鳥嘌呤"]),
    ("hypoxanthine", "Hypoxanthine", "O=[c:2]1[nH]cnc2[nH:9]cnc12", PUR, True, "IMP 的鹼基。救援：HGPRT（+ PRPP → IMP）；分解：xanthine oxidase → xanthine。", ["次黃嘌呤"]),
    ("xanthine", "Xanthine", "O=[c:2]1[nH]c(=O)[nH]c2[nH:9]cnc12", PUR, True, "XMP 的鹼基。xanthine oxidase：hypoxanthine → xanthine → uric acid（allopurinol 抑制）。", ["黃嘌呤"]),
    ("uric_acid", "Uric acid", "O=[c:2]1[nH]c(=O)[nH]c2[nH:9]c(=O)[nH]c12", PUR, True, "人類嘌呤分解終產物（缺 urate oxidase）。過多 → 痛風（gout）、Lesch-Nyhan。生理 pH 下主要為 urate。", ["尿酸", "urate"]),
    ("cytosine", "Cytosine", "Nc1cc[nH:1][c:2](=O)n1", PYR, False, "嘧啶鹼基。脫胺 → uracil（DNA 用 T 而非 U 的理由之一）。", ["C", "胞嘧啶"]),
    ("uracil", "Uracil", "O=c1cc[nH:1][c:2](=O)[nH]1", PYR, False, "RNA 的嘧啶鹼基。5-FU 為其類似物（→ FdUMP 抑制 thymidylate synthase）。", ["U", "尿嘧啶"]),
    ("thymine", "Thymine", "Cc1c[nH:1][c:2](=O)[nH]c1=O", PYR, False, "5-methyluracil，DNA 專用。甲基來自 N⁵,N¹⁰-methylene-THF（thymidylate synthase）。", ["T", "胸腺嘧啶"]),
    ("orotate", "Orotate", "O=C([O-])c1cc(=O)[nH][c:2](=O)[nH:1]1", PYR, False, "嘧啶從頭合成中第一個完整的嘧啶環（DHODH 產物）；OPRT 接上 PRPP → OMP。UMP synthase 缺乏 → orotic aciduria。", ["乳清酸"]),
]:
    smi, r = base_full(*FULL[id_]) if id_ in FULL else free_base(mapped)
    write(id_, name, smi, {"ring": r}, tags, notes, al)

smi, r = free_base("O=C([O-])[C@@H]1CC(=O)N[C:2](=O)[NH:1]1")
write("dihydroorotate", "Dihydroorotate", smi, {"ring": r}, ["nucleotide", "pyrimidine"],
      "L-dihydroorotate。dihydroorotase 環化產物；DHODH（粒線體內膜，接 CoQ）氧化成 orotate。", ["DHO"])

# ---------------------------------------------------------------- nucleosides / nucleotides
for id_, name, base, c5, deoxy, tags, notes, al in [
    ("amp", "AMP", "A", P1, False, PUR, "Adenosine 5′-monophosphate。IMP + Asp（GTP）→ adenylosuccinate → AMP + fumarate。", ["adenylate", "腺苷酸"]),
    ("gmp", "GMP", "G", P1, False, PUR, "Guanosine 5′-monophosphate。IMP → XMP（IMP DH，NAD⁺）→ GMP（GMP synthetase，Gln + ATP）。", ["guanylate"]),
    ("imp", "IMP", "I", P1, False, PUR, "Inosine monophosphate：嘌呤從頭合成的第一個完整嘌呤核苷酸，AMP 與 GMP 的分歧點。", ["inosinate", "肌苷酸"]),
    ("xmp", "XMP", "X", P1, False, PUR, "Xanthosine monophosphate：IMP → GMP 的中間物。", ["xanthylate"]),
    ("adenylosuccinate", "Adenylosuccinate", "SAMP", P1, False, PUR, "IMP + Asp（adenylosuccinate synthetase，耗 GTP）的產物；lyase 切出 fumarate → AMP（與尿素循環 ASL 類似）。", ["S-AMP"]),
    ("ump", "UMP", "U", P1, False, PYR, "Uridine monophosphate：嘧啶從頭合成的產物（OMP decarboxylase）。", ["uridylate"]),
    ("cmp", "CMP", "C", P1, False, PYR, "Cytidine monophosphate。CTP 由 UTP 胺化而來（CTP synthetase，Gln）。", ["cytidylate"]),
    ("omp", "OMP", "O", P1, False, PYR, "Orotidine monophosphate（orotidylate）：OPRT 產物，OMP decarboxylase → UMP。", ["orotidylate"]),
    ("dump", "dUMP", "U", P1, True, PYR, "Deoxyuridine monophosphate：thymidylate synthase 的受質（→ dTMP）。", []),
    ("dtmp", "dTMP", "T", P1, True, PYR, "Deoxythymidine monophosphate。5-FU（FdUMP）抑制 thymidylate synthase；methotrexate 抑制 DHFR 使 THF 耗竭。", ["thymidylate"]),
    ("utp", "UTP", "U", P3, False, PYR, "Uridine triphosphate：CTP synthetase 的受質。", []),
    ("ctp", "CTP", "C", P3, False, PYR, "Cytidine triphosphate。CTP 回饋抑制 ATCase（細菌）／CPS II。", []),
    ("atp_n", "ATP", "A", P3, False, PUR, "三磷酸核苷酸的結構範例：α、β、γ 磷酸。（輔因子表中的 ATP 是另一個簡寫項目）", ["adenosine triphosphate"]),
    ("adenosine", "Adenosine", "A", OH5, False, PUR, "核苷（無磷酸）。ADA（adenosine deaminase）→ inosine；ADA 缺乏 → dATP 累積抑制 RNR → SCID。", ["腺苷"]),
    ("inosine", "Inosine", "I", OH5, False, PUR, "PNP（purine nucleoside phosphorylase）→ hypoxanthine + ribose 1-P。", ["肌苷"]),
    ("guanosine", "Guanosine", "G", OH5, False, PUR, "PNP → guanine + ribose 1-P。", ["鳥苷"]),
]:
    smi, h = nuc(base, c5, deoxy)
    write(id_, name, smi, {"haworth": h}, tags, notes, al)

# generic NDP / dNDP for ribonucleotide reductase
GEN = "[*]" + sugar("1", P2)
GEN_D = "[*]" + sugar("1", P2, True)
for id_, name, smi, deoxy, notes in [
    ("ndp", "Ribonucleoside diphosphate", GEN, False, "NDP（ADP、GDP、CDP、UDP）：ribonucleotide reductase 的受質。"),
    ("dndp", "Deoxyribonucleoside diphosphate", GEN_D, True, "dNDP：RNR 產物，2′-OH 被還原成 H。"),
]:
    h = {"ring": "furanose", "subs": {1: {"up": "Base", "down": "H"}, 2: {"up": "H", "down": "H" if deoxy else "OH"},
                                      3: {"up": "H", "down": "OH"}, 4: {"up": LAB[P2], "down": "H"}}}
    write(id_, name, smi, {"haworth": h}, ["nucleotide"], notes, ["dNDP" if deoxy else "NDP"],
          label="dNDP" if deoxy else "NDP")

# PRPP (α: pyrophosphate down) and 5-phosphoribosylamine (β: amine up)
prpp = "[O-]P(=O)([O-])OP(=O)([O-])O[C@H]1O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]1O"
pra = "N[C@@H]1O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]1O"
write("prpp", "PRPP", prpp, {"haworth": {"ring": "furanose", "subs": {
    1: {"up": "H", "down": "O–{P}–{P}"}, 2: {"up": "H", "down": "OH"}, 3: {"up": "H", "down": "OH"},
    4: {"up": "{P}OCH_2", "down": "H"}}}}, ["nucleotide"],
    "5-Phosphoribosyl-1-pyrophosphate（α）。R5P + ATP（PRPP synthetase）。嘌呤、嘧啶從頭合成與救援途徑都需要；PRPP 過多 → 嘌呤過度合成 → 痛風。",
    ["5-phosphoribosyl 1-pyrophosphate"])
write("pra", "5-Phosphoribosylamine", pra, {"haworth": {"ring": "furanose", "subs": {
    1: {"up": "NH_2", "down": "H"}, 2: {"up": "H", "down": "OH"}, 3: {"up": "H", "down": "OH"},
    4: {"up": "{P}OCH_2", "down": "H"}}}}, ["nucleotide", "purine"],
    "PRA：Gln-PRPP amidotransferase 的產物（嘌呤合成限速步驟）。注意 C1′ 由 α 轉成 β（反轉）。",
    ["PRA"], label="5-Phosphoribosyl-\namine")
print("ok")
