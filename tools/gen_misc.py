"""Round 14: glycosaminoglycans & structural polysaccharides, glyoxylate cycle, Calvin cycle,
photorespiration, mucopolysaccharidoses.

    python tools/gen_misc.py data

Repeating disaccharides are assembled with RDKit from the monosaccharide SMILES already in the
database (glycosidic bond made by re-pointing the donor's C1, keeping its chirality), so the
validator can check every Haworth up/down against real stereochemistry.
"""
import sys
from pathlib import Path

import yaml
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).parent))
from reaction_colors import COLORS  # noqa: E402

M, C = "magenta", "cyan"

# ---------------------------------------------------------------- monosaccharide building blocks

GLC_B = "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O"
GAL_B = "C([C@@H]1[C@@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O"
GLCA_B = "[O-]C(=O)[C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O"
GLCNAC_B = "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)NC(C)=O)O)O)O"
GALNAC_B = "C([C@@H]1[C@@H]([C@@H]([C@H]([C@@H](O1)O)NC(C)=O)O)O)O"
GLCN_A = "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)[NH3+])O)O)O"        # α-D-glucosamine


def positions(mol, offset=0):
    """{1..6: (carbon idx, exocyclic heteroatom idx)} for a pyranose (C6 hetero only if CH2–X)."""
    ring = next(r for r in mol.GetRingInfo().AtomRings() if len(r) == 6)
    o = next(i for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == "O")
    nb = [n.GetIdx() for n in mol.GetAtomWithIdx(o).GetNeighbors()]

    def exo_het(c):
        return [n.GetIdx() for n in mol.GetAtomWithIdx(c).GetNeighbors()
                if n.GetIdx() not in ring and n.GetSymbol() in "ON"]
    c1 = next(c for c in nb if exo_het(c))
    seq = [c1]
    prev = o
    while len(seq) < 5:
        cur = seq[-1]
        nxt = next(n.GetIdx() for n in mol.GetAtomWithIdx(cur).GetNeighbors()
                   if n.GetIdx() in ring and n.GetIdx() not in (prev, o) and n.GetIdx() not in seq)
        prev = cur
        seq.append(nxt)
    out = {}
    for k, c in enumerate(seq, start=1):
        h = exo_het(c) if k < 5 else []
        out[k] = (c + offset, (h[0] + offset) if h else None)
    c5 = seq[4]
    c6 = next((n.GetIdx() for n in mol.GetAtomWithIdx(c5).GetNeighbors()
               if n.GetIdx() not in ring and n.GetSymbol() == "C"), None)
    if c6 is not None:
        h = [n.GetIdx() for n in mol.GetAtomWithIdx(c6).GetNeighbors() if n.GetSymbol() == "O"]
        out[6] = (c6 + offset, (h[0] + offset) if len(h) == 1 else None)
    return out


def _parity(a, b):
    """+1 if list b is an even permutation of list a, else −1."""
    b = list(b)
    sign = 1
    for i, x in enumerate(a):
        j = b.index(x)
        if j != i:
            b[i], b[j] = b[j], b[i]
            sign = -sign
    return sign


def _repoint(rw, atom, old, new):
    """Swap neighbour `old` of `atom` for `new` without changing its configuration."""
    before = [b.GetOtherAtomIdx(atom) for b in rw.GetAtomWithIdx(atom).GetBonds()]
    rw.RemoveBond(atom, old)
    rw.AddBond(atom, new, Chem.BondType.SINGLE)
    after = [b.GetOtherAtomIdx(atom) for b in rw.GetAtomWithIdx(atom).GetBonds()]
    mapped = [new if x == old else x for x in before]
    if rw.GetAtomWithIdx(atom).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED and _parity(mapped, after) < 0:
        rw.GetAtomWithIdx(atom).InvertChirality()


def _free(atom):
    atom.SetNoImplicit(False)
    atom.SetNumExplicitHs(0)


def _sulfate(rw, het):
    """X–H → X–SO₃⁻ (O-sulfate or N-sulfate)."""
    a = rw.GetAtomWithIdx(het)
    a.SetFormalCharge(0)
    _free(a)
    s = rw.AddAtom(Chem.Atom(16))
    rw.AddBond(het, s, Chem.BondType.SINGLE)
    for k in range(3):
        o = rw.AddAtom(Chem.Atom(8))
        rw.AddBond(s, o, Chem.BondType.DOUBLE if k < 2 else Chem.BondType.SINGLE)
        if k == 2:
            rw.GetAtomWithIdx(o).SetFormalCharge(-1)


def disaccharide(donor, acceptor, at, sulf_d=(), sulf_a=(), nsulf_a=False, invert_c5_d=False, lactyl_a=False,
                 invert_c5_a=False):
    """donor C1 –O– acceptor C<at>. sulf_*: positions carrying O-sulfate; nsulf_a: N-sulfate on the
    acceptor's C2 amine; invert_c5_d: C5 epimer of the donor (GlcA → L-IdoA)."""
    d, a = Chem.MolFromSmiles(donor), Chem.MolFromSmiles(acceptor)
    pd, pa = positions(d), positions(a, d.GetNumAtoms())
    rw = Chem.RWMol(Chem.CombineMols(d, a))
    if invert_c5_d:
        rw.GetAtomWithIdx(pd[5][0]).InvertChirality()
    if invert_c5_a:
        rw.GetAtomWithIdx(pa[5][0]).InvertChirality()
    for k in sulf_d:
        _sulfate(rw, pd[k][1])
    for k in sulf_a:
        _sulfate(rw, pa[k][1])
    if nsulf_a:
        _sulfate(rw, pa[2][1])
    if lactyl_a:                                   # 3-O-(R)-lactyl ether: NAM
        o3 = pa[3][1]
        lac = Chem.MolFromSmiles("C[C@H](C(=O)[O-])Cl")
        off = rw.GetNumAtoms()
        rw2 = Chem.RWMol(Chem.CombineMols(rw, lac))
        cl = next(x.GetIdx() for x in rw2.GetAtoms() if x.GetSymbol() == "Cl")
        cc = next(x.GetIdx() for x in rw2.GetAtomWithIdx(cl).GetNeighbors())
        _repoint(rw2, cc, cl, o3)
        _free(rw2.GetAtomWithIdx(o3))
        rw2.RemoveAtom(cl)
        rw = rw2
        del off
    c1, o1 = pd[1]
    oa = pa[at][1]
    _repoint(rw, c1, o1, oa)
    _free(rw.GetAtomWithIdx(oa))
    rw.RemoveAtom(o1)
    m = rw.GetMol()
    Chem.SanitizeMol(m)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    return Chem.MolToSmiles(m)


# Haworth subs. Unit 0 = donor (normal orientation, β: C1 up). Unit 1 = acceptor; flipped when the
# link goes to its C4 (like lactose), normal when it goes to C3.

def glc_like(c1=("OH", "H"), c2=("H", "OH"), c3=("OH", "H"), c4=("H", "HO"), c5=("CH_2OH", "H")):
    return {1: {"up": c1[0], "down": c1[1]}, 2: {"up": c2[0], "down": c2[1]}, 3: {"up": c3[0], "down": c3[1]},
            4: {"up": c4[0], "down": c4[1]}, 5: {"up": c5[0], "down": c5[1]}}


def flipped(subs):
    return {k: {"up": v["down"], "down": v["up"]} for k, v in subs.items()}


def two(u0, u1, link_to, flip1, at=(4.9, -0.2), label=None):
    units = [{"ring": "pyranose", "subs": u0},
             {"ring": "pyranose", "subs": u1, "at": list(at), **({"flip": True} if flip1 else {})}]
    side = "up"
    ln = {"from": [0, 1, "up"], "to": [1, link_to, side]}
    if label:
        ln["label"] = label
    return {"units": units, "links": [ln]}


SULF = "OSO_3^-"
SULF_L = ["^-O_3SO", 4]                    # on the left of a vertical bond: bond meets the last O
NS = "NHSO_3^-"
NHAC = "NHAc"

IDO_C5 = {"up": "H", "down": "COO^-", "down_len": 0.6, "down_dx": 0.35}               # L-IdoA: COO⁻ below C5
IDOA2S = {**glc_like(c1=("link", "H"), c2=("H", SULF)), 5: IDO_C5}     # α-L: C1–O up, COO⁻ down
IDOA2S[2].pop("up")                       # ring-interior H left out (textbook style) to make room for COO⁻
IDOA2S[3]["up_len"] = 0.4
GAL_D = glc_like(c1=("link", "H"), c4=("HO", "H"))
GLCNAC_D = glc_like(c1=("link", "H"), c2=("H", NHAC))
GALNAC4S_D = glc_like(c1=("link", "H"), c2=("H", NHAC), c4=(SULF_L, "H"))
GLCA_A = flipped(glc_like(c4=("H", "link"), c5=("COO^-", "H")))                 # acceptor, turned over
IDOA_A = flipped({**glc_like(c4=("H", "link")), 5: {"up": "H", "down": "COO^-"}})
IDOA_A[5].update(up_len=0.6, up_dx=0.35)                                                       # (flipped: COO⁻ now "up")
IDOA_A[2].pop("down")                     # ring-interior H's left out (textbook style) to make room for COO⁻
IDOA_A[3]["down_len"] = 0.4

# Each repeating disaccharide is drawn as its β1→4 (or α1→4) half, which a Haworth pair shows
# cleanly (second ring turned over, like lactose); the other linkage of the polymer is in the notes.
GAG = [
    {
        "id": "hyaluronate_unit", "name": "Hyaluronate (repeating unit)", "label": "Hyaluronate",
        "aliases": ["hyaluronic acid", "hyaluronan", "玻尿酸", "透明質酸"],
        "smiles": disaccharide(GLCNAC_B, GLCA_B, 4),
        "haworth": two(GLCNAC_D, GLCA_A, 4, True),
        "notes": "重複雙醣：[–4)GlcA(β1→3)GlcNAc(β1–]ₙ，圖為 GlcNAc(β1→4)GlcA 這一半。唯一「不含硫酸」、也不以共價鍵接 core protein 的 GAG；"
                 "分子極大、高度水合 → 關節液潤滑、玻璃體、軟骨中作為 aggrecan 的骨架（經 link protein 非共價結合）。"
                 "嚴格說玻尿酸是 GAG，不是 proteoglycan。",
    },
    {
        "id": "chondroitin_4_sulfate_unit", "name": "Chondroitin 4-sulfate (repeating unit)",
        "label": "Chondroitin 4-sulfate",
        "aliases": ["chondroitin sulfate", "軟骨素硫酸"],
        "smiles": disaccharide(GALNAC_B, GLCA_B, 4, sulf_d=(4,)),
        "haworth": two(GALNAC4S_D, GLCA_A, 4, True),
        "notes": "重複雙醣：[–4)GlcA(β1→3)GalNAc-4-sulfate(β1–]ₙ，圖為 GalNAc4S(β1→4)GlcA。軟骨中最多的 GAG（aggrecan），抗壓。"
                 "GalNAc 的 C4–O 朝上（galacto），硫酸接在這裡。",
    },
    {
        "id": "dermatan_sulfate_unit", "name": "Dermatan sulfate (repeating unit)", "label": "Dermatan sulfate",
        "aliases": ["dermatan sulfate", "chondroitin sulfate B", "皮膚素硫酸"],
        "smiles": disaccharide(GALNAC_B, GLCA_B, 4, sulf_d=(4,), invert_c5_a=True),
        "haworth": two(GALNAC4S_D, IDOA_A, 4, True),
        "notes": "重複雙醣：[–4)L-IdoA(α1→3)GalNAc-4-sulfate(β1–]ₙ，圖為 GalNAc4S(β1→4)IdoA。與 chondroitin 4-sulfate 只差在 GlcA 的 C5 被翻成 L-iduronate（C5 epimerase；COO⁻ 改朝下）。"
                 "皮膚、血管、心瓣膜。Hurler／Hunter 症候群會累積 dermatan sulfate 與 heparan sulfate。",
    },
    {
        "id": "heparin_unit", "name": "Heparin (repeating unit)", "label": "Heparin",
        "aliases": ["heparin", "heparan sulfate", "肝素"],
        "smiles": disaccharide(GLCA_B, GLCN_A, 4, sulf_d=(2,), sulf_a=(6,), nsulf_a=True, invert_c5_d=True),
        "haworth": two(IDOA2S, flipped(glc_like(c1=("H", "OH"), c2=("H", NS), c4=("H", "link"),
                                                c5=("CH_2OSO_3^-", "H"))), 4, True),
        "notes": "IdoA-2-sulfate(α1→4)GlcNS-6-sulfate，單元間 α1→4。帶負電密度最高的生物分子。"
                 "由 mast cell 分泌；結合 antithrombin III 使其抑制 thrombin、factor Xa 的速度加快約 1000 倍 → 抗凝血。"
                 "Heparan sulfate 結構相近但硫酸較少、GlcA 較多，分布在細胞表面（syndecan）。",
    },
    {
        "id": "keratan_sulfate_unit", "name": "Keratan sulfate (repeating unit)", "label": "Keratan sulfate",
        "aliases": ["keratan sulfate", "角質素硫酸"],
        "smiles": disaccharide(GAL_B, GLCNAC_B, 4, sulf_a=(6,)),
        "haworth": two(GAL_D, flipped(glc_like(c2=("H", NHAC), c4=("H", "link"), c5=("CH_2OSO_3^-", "H"))), 4, True),
        "notes": "Gal(β1→4)GlcNAc-6-sulfate，單元間 β1→3。唯一不含 uronic acid 的 GAG（以 galactose 代替）。"
                 "角膜（透明度）、軟骨、骨。Morquio 症候群（MPS IV）累積 keratan sulfate。",
    },
    {
        "id": "chitin_unit", "name": "Chitin (repeating unit)", "label": "Chitin",
        "aliases": ["chitin", "chitobiose", "幾丁質"],
        "smiles": disaccharide(GLCNAC_B, GLCNAC_B, 4),
        "haworth": two(GLCNAC_D, flipped(glc_like(c2=("H", NHAC), c4=("H", "link"))), 4, True),
        "notes": "GlcNAc(β1→4)GlcNAc 的同元多醣（homopolysaccharide）；與 cellulose 只差在 C2 的 OH 換成 NHAc。"
                 "節肢動物外骨骼、真菌細胞壁。",
    },
    {
        "id": "peptidoglycan_unit", "name": "Peptidoglycan (glycan unit)", "label": "Peptidoglycan (NAG–NAM)",
        "aliases": ["peptidoglycan", "murein", "NAG-NAM", "肽聚醣"],
        "smiles": disaccharide(GLCNAC_B, GLCNAC_B, 4, lactyl_a=True),
        "haworth": two(GLCNAC_D, flipped(glc_like(c2=("H", NHAC), c3=("OR", "H"), c4=("H", "link"))), 4, True),
        "notes": "細菌細胞壁：NAG(β1→4)NAM 交替（雜多醣）。NAM = GlcNAc 的 C3–O 接 D-lactyl（圖中 OR，R = –CH(CH₃)COO⁻），"
                 "lactyl 的 COO⁻ 再接四胜肽（L-Ala–D-Glu–L-Lys/DAP–D-Ala），胜肽之間交聯。"
                 "Lysozyme 切 NAM–NAG 的 β1→4 鍵；penicillin 抑制交聯的 transpeptidase。",
        "lactyl_note": True,
    },
]


# ---------------------------------------------------------------- small molecules

MOLECULES = [
    {"id": "glyoxylate", "name": "Glyoxylate", "aliases": ["乙醛酸", "glyoxylic acid"],
     "smiles": "O=CC(=O)[O-]",
     "chain": [{"c": "COO^-"}, {"c": "C", "l": "H", "r": "O", "rb": 2}],
     "tags": ["glyoxylate_cycle"],
     "notes": "乙醛酸：isocitrate lyase 把 isocitrate 切成 succinate + glyoxylate；malate synthase 再把它和 acetyl-CoA 接成 malate。"
              "人體內 glyoxylate 由 glycine／hydroxyproline 代謝產生，經 AGT 轉成 glycine；AGT 缺乏（primary hyperoxaluria type 1）→ oxalate 腎結石。"},
    {"id": "ribulose_15bp", "name": "Ribulose 1,5-bisphosphate", "label": "Ribulose\n1,5-bisphosphate",
     "aliases": ["RuBP", "核酮糖-1,5-二磷酸"],
     "smiles": "O=C(COP(=O)([O-])[O-])[C@H](O)[C@H](O)COP(=O)([O-])[O-]",
     "chain": [{"c": "CH_2O{P}"}, {"c": "C", "r": "O", "rb": 2}, {"c": "C", "l": "H", "r": "OH"},
               {"c": "C", "l": "H", "r": "OH"}, {"c": "CH_2O{P}"}],
     "tags": ["calvin_cycle"],
     "notes": "RuBP：Rubisco 的受質。接上 CO₂ 後斷成 2 分子 3-phosphoglycerate（Calvin cycle）；接上 O₂ 則為光呼吸。"},
    {"id": "phosphoglycolate_2", "name": "2-Phosphoglycolate", "label": "2-Phospho-\nglycolate",
     "aliases": ["2-PG", "phosphoglycolate", "2-磷酸乙醇酸"],
     "smiles": "[O-]C(=O)COP(=O)([O-])[O-]",
     "chain": [{"c": "C", "fork": {"l": "O", "lb": 2, "r": "O^-"}}, {"c": "CH_2O{P}"}],
     "tags": ["calvin_cycle"],
     "notes": "Rubisco 的 oxygenase 產物（光呼吸）：2 碳、無法直接用，要經過 peroxisome 與粒線體（photorespiration）回收成 3-PG，損失 CO₂ 與 NH₃。"},
]


REACTIONS = [
    {
        "id": "isocitrate_lyase", "title": "Isocitrate lyase", "substrates": ["isocitrate"],
        "products": ["succinate", "glyoxylate"], "enzyme": "Isocitrate lyase", "reversible": False,
        "card": {"location": "Glyoxysome（植物種子）、細菌、真菌；動物沒有",
                 "type": "醛醇裂解（retro-aldol）：C2–C3 斷開",
                 "exam_points": ["乙醛酸循環的兩個獨特酵素之一（另一個是 malate synthase）",
                                 "跳過 TCA 中兩個脫羧步驟（IDH、α-KGDH），碳不以 CO₂ 流失"]},
        "tags": ["glyoxylate_cycle"],
    },
    {
        "id": "malate_synthase", "title": "Malate synthase", "substrates": ["glyoxylate", "acetyl_coa"],
        "products": ["malate"], "cofactors_in": ["H2O"], "cofactors_out": ["CoA-SH", "H+"],
        "enzyme": "Malate synthase", "reversible": False,
        "card": {"location": "Glyoxysome；動物沒有",
                 "type": "Claisen 縮合（同 citrate synthase）：acetyl-CoA 的甲基碳接到 glyoxylate 的醛基",
                 "exam_points": ["2 acetyl-CoA → 1 succinate（4 碳）淨合成 → 可走糖質新生",
                                 "動物缺此酵素與 isocitrate lyase → 偶數碳脂肪酸無法淨合成葡萄糖"]},
        "tags": ["glyoxylate_cycle"],
    },
    {
        "id": "rubisco_carboxylase", "title": "Rubisco (carboxylation)",
        "substrates": ["ribulose_15bp"], "products": [{"id": "pg_3", "n": 2}],
        "cofactors_in": ["CO2", "H2O"], "cofactors_out": ["2H+"],
        "enzyme": "Rubisco\n(ribulose-1,5-bisphosphate carboxylase/oxygenase)", "reversible": False,
        "card": {"location": "葉綠體 stroma",
                 "type": "羧化：CO₂ 接到 RuBP 的 C2（經 enediolate 中間物），6 碳中間物水解成 2 × 3-PG；洋紅 CO₂ 成為其中一個 3-PG 的 COO⁻",
                 "regulation": ["光照 → stroma pH ↑、Mg²⁺ ↑ → 活化；需 carbamoylation（Lys + CO₂ + Mg²⁺）",
                                "Rubisco activase 移除抑制性的糖磷酸"],
                 "exam_points": ["地球上最多的蛋白質；催化很慢（約 3 次／秒）",
                                 "同時有 oxygenase 活性 → 光呼吸；C4、CAM 植物先把 CO₂ 濃縮以減少光呼吸"]},
        "tags": ["calvin_cycle"],
    },
    {
        "id": "rubisco_oxygenase", "title": "Rubisco (oxygenation, photorespiration)",
        "substrates": ["ribulose_15bp"], "products": ["pg_3", "phosphoglycolate_2"],
        "cofactors_in": ["O2"], "cofactors_out": ["2H+"],
        "enzyme": "Rubisco (oxygenase)", "reversible": False,
        "card": {"location": "葉綠體 stroma",
                 "type": "O₂ 取代 CO₂ 加到 RuBP：產物 1 個 3-PG + 1 個 2-phosphoglycolate；洋紅 O 進入 2-phosphoglycolate 的羧基",
                 "exam_points": ["高溫、低 CO₂/O₂ 比時增加",
                                 "2-phosphoglycolate 的回收（photorespiration）跨葉綠體、peroxisome、粒線體，放出 CO₂ 與 NH₃"]},
        "tags": ["calvin_cycle"],
    },
]


def n(id_, ref, x, y, **kw):
    return {"id": id_, "ref": ref, "pos": [x, y], "view": "name", **kw}


def t(id_, text, x, y, kind="concept", **kw):
    return {"id": id_, "text": text, "pos": [x, y], "kind": kind, **kw}


PATHWAYS = [
    {
        "id": "glyoxylate_cycle", "title": "Glyoxylate cycle", "layout": "free", "cell_size": 4.2,
        "nodes": [
            n("ac1", "acetyl_coa", -1.6, -2.2), n("oaa", "oxaloacetate", -1.6, 0), n("cit", "citrate", 0, -1.2),
            n("iso", "isocitrate", 1.6, 0), n("gly", "glyoxylate", 0, 1.3), n("suc", "succinate", 3.2, 1.3),
            n("mal", "malate", -1.6, 1.9), n("ac2", "acetyl_coa", 1.4, 2.8),
            t("gng", "→ mitochondria → OAA\n→ gluconeogenesis", 3.2, 2.6, "name"),
            t("note", "Plants (glyoxysome), bacteria, fungi\nNOT animals", 1.6, -2.2, "heading"),
        ],
        "edges": [
            {"from": "oaa", "to": "cit", "step": "citrate_synthase", "label": "Citrate\nsynthase"},
            {"from": "ac1", "to": "edge:oaa->cit", "bend": 0},
            {"from": "cit", "to": "iso", "label": "Aconitase"},
            {"from": "iso", "to": "gly", "step": "isocitrate_lyase", "label": "Isocitrate\nlyase", "color": "magenta"},
            {"from": "iso", "to": "suc", "step": "isocitrate_lyase", "color": "magenta"},
            {"from": "gly", "to": "mal", "step": "malate_synthase", "label": "Malate\nsynthase", "color": "magenta"},
            {"from": "ac2", "to": "edge:gly->mal", "bend": 0},
            {"from": "mal", "to": "oaa", "step": "malate_dehydrogenase", "label": "Malate DH",
             "side_text": "NAD^+ → NADH"},
            {"from": "suc", "to": "gng", "head": True},
        ],
        "tags": ["glyoxylate_cycle"],
        "notes": "乙醛酸循環：2 acetyl-CoA → 1 succinate（淨合成 4 碳），跳過 TCA 的兩個脫羧步驟。"
                 "洋紅箭頭為兩個獨特酵素：isocitrate lyase、malate synthase。"
                 "植物種子把脂肪轉成糖靠它；動物沒有這兩個酵素 → 偶數碳脂肪酸（acetyl-CoA）無法淨合成葡萄糖"
                 "（PDH 不可逆，且 TCA 中 2 個碳以 CO₂ 離開）。奇數碳脂肪酸的 propionyl-CoA → succinyl-CoA 則可以。",
    },
    {
        "id": "calvin_cycle", "title": "Calvin cycle (overview)", "layout": "cycle", "radius": 1.85,
        "start_angle": -90, "cell_size": 4.4, "label_side": "left",
        "center_text": "Calvin\ncycle",
        "nodes": [
            {"id": "rubp", "ref": "ribulose_15bp", "view": "name"},
            {"id": "pg3", "ref": "pg_3", "view": "name"},
            {"id": "bpg", "ref": "bpg_13", "view": "name"},
            {"id": "g3p", "ref": "glyceraldehyde_3p", "view": "name"},
            {"id": "ru5p", "ref": "ribulose_5p", "view": "name"},
            {"id": "out", "text": "1 G3P exported\n→ sucrose, starch", "kind": "name", "pos": [-0.2, 2.75]},
            {"id": "sum", "text": "3 CO_2 + 9 ATP + 6 NADPH → 1 G3P", "kind": "heading", "pos": [-2.9, -2.55]},
        ],
        "edges": [
            {"from": "rubp", "to": "pg3", "step": "rubisco_carboxylase", "label": "Rubisco",
             "side_text": "3 CO_2 + 3 H_2O"},
            {"from": "pg3", "to": "bpg", "label": "PG kinase", "side_text": "6 ATP\n→ 6 ADP"},
            {"from": "bpg", "to": "g3p", "label": "G3P DH (NADP^+)", "side_text": "6 NADPH\n→ 6 NADP^+"},
            {"from": "g3p", "to": "ru5p", "label": "Regeneration\n(aldolase,\ntransketolase …)", "side_text": "5 G3P"},
            {"from": "ru5p", "to": "rubp", "label": "Phosphoribulo-\nkinase", "side_text": "3 ATP\n→ 3 ADP"},
            {"from": "g3p", "to": "out", "bend": 0},
        ],
        "tags": ["calvin_cycle"],
        "notes": "Calvin cycle（暗反應，stroma）：3 CO₂ + 9 ATP + 6 NADPH → 1 G3P。三階段：固定（Rubisco：3 RuBP → 6 3-PG）、還原（類似糖質新生的 PGK、GAPDH，但用 NADPH）、"
                 "再生（與 PPP 非氧化相同的 transketolase／aldolase 重排：5 G3P → 3 Ru5P）。ATP 與 NADPH 來自光反應。",
    },
    {
        "id": "photorespiration", "title": "Photorespiration", "layout": "free", "cell_size": 4.4,
        "regions": [
            {"kind": "box", "x0": -0.75, "y0": -0.7, "x1": 2.4, "y1": 3.55, "label": "Chloroplast"},
            {"kind": "box", "x0": 2.65, "y0": -0.7, "x1": 4.9, "y1": 3.55, "label": "Peroxisome"},
            {"kind": "box", "x0": 5.1, "y0": -0.7, "x1": 7.45, "y1": 3.55, "label": "Mitochondrion"},
        ],
        "nodes": [
            n("rubp", "ribulose_15bp", 0, 0.3), n("pgly", "phosphoglycolate_2", 1.75, 0.3),
            t("glyc", "Glycolate", 1.75, 1.3, "name"), n("pg3", "pg_3", 0, 1.9),
            t("glyc2", "Glycolate", 3.7, 0.3, "name"), n("glyox", "glyoxylate", 3.7, 1.3), n("gly", "glycine", 3.7, 2.2),
            t("gly2", "2 Glycine", 6.05, 2.2, "name"), n("ser", "serine", 6.05, 3.05),
            t("hp", "Hydroxypyruvate\n→ glycerate", 3.7, 3.05, "name"), t("gla", "Glycerate", 1.75, 3.05, "name"),
        ],
        "edges": [
            {"from": "rubp", "to": "pgly", "step": "rubisco_oxygenase", "label": "Rubisco + O_2", "color": "red"},
            {"from": "rubp", "to": "pg3", "step": "rubisco_oxygenase"},
            {"from": "pgly", "to": "glyc", "side_text": "P_i"},
            {"from": "glyc", "to": "glyc2"},
            {"from": "glyc2", "to": "glyox", "label": "Glycolate\noxidase", "side_text": "O_2 → H_2O_2"},
            {"from": "glyox", "to": "gly", "label": "Trans-\namination"},
            {"from": "gly", "to": "gly2"},
            {"from": "gly2", "to": "ser", "label": "Gly cleavage\n+ SHMT", "side_text": "CO_2 + NH_3\nlost"},
            {"from": "ser", "to": "hp"},
            {"from": "hp", "to": "gla"},
            {"from": "gla", "to": "pg3", "side_text": "ATP"},
        ],
        "tags": ["calvin_cycle"],
        "notes": "光呼吸：Rubisco 把 O₂ 加到 RuBP，產生 2-phosphoglycolate；回收路徑跨三個胞器，最後 2 個 2-PG（4 C）只換回 1 個 3-PG（3 C），"
                 "損失 1 個 CO₂ 與 NH₃，消耗 ATP。高溫、乾燥（氣孔關閉 → CO₂ 低）時更嚴重；C4 與 CAM 植物用 PEP carboxylase 先固定 CO₂ 來避免。",
    },
    {
        "id": "glycosaminoglycans", "title": "Glycosaminoglycans: repeating units", "layout": "free", "cell_size": 17,
        "nodes": [
            {"id": "ha", "ref": "hyaluronate_unit", "pos": [0, 0]},
            {"id": "cs", "ref": "chondroitin_4_sulfate_unit", "pos": [1.12, 0]},
            {"id": "ds", "ref": "dermatan_sulfate_unit", "pos": [0, 0.62]},
            {"id": "hep", "ref": "heparin_unit", "pos": [1.12, 0.62]},
            {"id": "ks", "ref": "keratan_sulfate_unit", "pos": [0, 1.24]},
            {"id": "note", "text": "Drawn: the β1→4 (heparin α1→4) half\nof each repeating disaccharide", "kind": "name",
             "pos": [1, 1.24]},
        ],
        "edges": [],
        "text": False,
        "tags": ["polysaccharide", "gag"],
        "notes": "GAG：由「胺基糖（GlcNAc／GalNAc／GlcN）＋ uronic acid（GlcA／IdoA）」的雙醣重複組成，多數帶硫酸 → 高負電、高度水合、潤滑抗壓。"
                 "例外：hyaluronate 無硫酸、不接蛋白；keratan sulfate 以 galactose 代替 uronic acid。"
                 "除 hyaluronate 外，GAG 以 Ser–Xyl–Gal–Gal–GlcA 連接區共價接在 core protein 上 → proteoglycan。",
    },
    {
        "id": "polysaccharide_classes", "title": "Homo- vs heteropolysaccharides", "layout": "free", "cell_size": 4.5,
        "regions": [{"kind": "box", "x0": -1.45, "y0": 1.55, "x1": 1.45, "y1": 3.2},
                    {"kind": "box", "x0": 1.6, "y0": 1.55, "x1": 6.15, "y1": 3.2}],
        "nodes": [
            t("root", "Polysaccharides", 2.05, -0.1, "concept"),
            t("homo", "Homopolysaccharides\n(one kind of monomer)", 0, 0.9, "concept"),
            t("hetero", "Heteropolysaccharides\n(two or more kinds)", 3.85, 0.9, "concept"),
            n("amy", "amylose", -0.7, 2.0), n("glyg", "glycogen", 0.7, 2.0), n("cel", "cellulose", -0.7, 2.75),
            n("chi", "chitin_unit", 0.7, 2.75),
            n("ha", "hyaluronate_unit", 2.65, 2.0), n("hep", "heparin_unit", 4.95, 2.0),
            n("cs", "chondroitin_4_sulfate_unit", 2.65, 2.75), n("pg", "peptidoglycan_unit", 4.95, 2.75),
            t("hnote", "Glc α1→4 / α1→6: energy storage\nGlc β1→4, GlcNAc β1→4: structure", 0, 3.7, "name"),
            t("tnote", "GAGs: amino sugar + uronic acid (± sulfate)\nPeptidoglycan: NAG–NAM + peptide", 3.65, 3.7, "name"),
        ],
        "edges": [{"from": "root", "to": "homo", "head": False}, {"from": "root", "to": "hetero", "head": False}],
        "tags": ["polysaccharide"],
        "notes": "同元多醣：starch（amylose、amylopectin）、glycogen、cellulose、chitin、dextran、inulin（fructose）。"
                 "雜多醣：GAG（hyaluronate、chondroitin、dermatan、heparin、keratan）、peptidoglycan、agarose。"
                 "考點：玻尿酸是 GAG（雜多醣），本身不是 proteoglycan；chitin 是同元多醣（只有 GlcNAc）。",
    },
    {
        "id": "mucopolysaccharidoses", "title": "Mucopolysaccharidoses (lysosomal GAG degradation)",
        "layout": "free", "cell_size": 4.6,
        "nodes": [
            t("h", "Lysosomal storage of GAGs (urine GAGs ↑)", -1.4, -0.75, "heading"),
            t("hurler", "Hurler (MPS I)\nα-L-iduronidase\nDS + HS · corneal clouding\nautosomal recessive", 0, 0.4),
            t("hunter", "Hunter (MPS II)\niduronate-2-sulfatase\nDS + HS · NO corneal clouding\nX-linked recessive", 3.05, 0.4),
            t("sanfilippo", "Sanfilippo (MPS III)\nheparan sulfate–\ndegrading enzymes\nHS · severe CNS", 6.1, 0.4),
            t("morquio", "Morquio (MPS IV)\ngalactose-6-sulfatase /\nβ-galactosidase\nKS · skeletal dysplasia", 0, 1.75),
            t("sly", "Sly (MPS VII)\nβ-glucuronidase\nDS + HS + CS", 3.05, 1.75),
            t("ml", "I-cell disease\nno mannose-6-P tag →\nenzymes secreted, not\ndelivered to lysosomes", 6.1, 1.75, "drug"),
        ],
        "edges": [],
        "tags": ["polysaccharide", "gag"],
        "notes": "黏多醣症（MPS）：溶體分解 GAG 的酵素缺乏，GAG 堆積、尿中 GAG 上升；共同表現為粗糙臉孔、骨骼變形、肝脾腫大、智能退化。"
                 "口訣：Hunter「看得清楚獵物」→ 無角膜混濁、X-linked；Hurler 有角膜混濁、體染色體隱性。"
                 "DS = dermatan sulfate，HS = heparan sulfate，KS = keratan sulfate，CS = chondroitin sulfate。"
                 "I-cell disease（mucolipidosis II）：GlcNAc phosphotransferase 缺乏，無法加上 mannose-6-phosphate 標籤。",
    },
]


def dump(path, d):
    path.write_text(yaml.safe_dump(d, allow_unicode=True, sort_keys=False, width=200), encoding="utf-8")


PARENT = {"glyoxylate_cycle": "carb_metabolism", "calvin_cycle": "carb_metabolism",
          "polysaccharide": "carbohydrate", "gag": "carbohydrate"}


def with_parent(tags):
    out = []
    for x in tags:
        if x in PARENT:
            out.append(PARENT[x])
        out.append(x)
    return list(dict.fromkeys(out))


def main(data):
    data = Path(data)
    for g in GAG:
        d = {"id": g["id"], "name": g["name"], "label": g["label"], "aliases": g["aliases"], "smiles": g["smiles"],
             "haworth": g["haworth"], "tags": with_parent(["polysaccharide", "gag"] if g["id"] not in
                                                          ("chitin_unit", "peptidoglycan_unit") else ["polysaccharide"]),
             "notes_zh": g["notes"], "status": "pending"}
        dump(data / "molecules" / f"{g['id']}.yaml", d)
    for m in MOLECULES:
        d = {"id": m["id"], "name": m["name"]}
        if m.get("label"):
            d["label"] = m["label"]
        d.update(aliases=m["aliases"], smiles=m["smiles"], chain=m["chain"], tags=with_parent(m["tags"]),
                 notes_zh=m["notes"], status="pending")
        dump(data / "molecules" / f"{m['id']}.yaml", d)
    for r in REACTIONS:
        d = {"id": r["id"], "title": r["title"], "substrates": r["substrates"], "products": r["products"]}
        for k in ("cofactors_in", "cofactors_out"):
            if r.get(k):
                d[k] = r[k]
        d["enzyme"] = r["enzyme"]
        d["reversible"] = r["reversible"]
        deco = COLORS.get(r["id"])
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
    print(f"{len(GAG) + len(MOLECULES)} molecules, {len(REACTIONS)} reactions, {len(PATHWAYS)} pathways")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data")
