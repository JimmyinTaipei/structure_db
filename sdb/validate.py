"""Chemistry checks that catch the mistakes an AI author typically makes.

1. every SMILES parses
2. reactions balance atoms and charge (cofactors and extra_in/extra_out included)
3. chain drawings agree with the SMILES: carbon count, and for Fischer rows the
   left/right placement of H vs OH matches the SMILES stereochemistry
4. reaction colours (decorate) point at something that is drawn: a participant, an existing
   row / part / text fragment, a SMARTS that matches, atom indices inside the molecule
"""
import re
from collections import Counter

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Geometry import Point3D

from .text import anchor_index, plain

RDLogger.DisableLog("rdApp.*")


def mol_from(smiles):
    m = Chem.MolFromSmiles(smiles) if smiles else None
    if m is None:
        raise ValueError(f"SMILES 無法解析：{smiles!r}")
    return m


def composition(smiles):
    """Atom counts (H explicit; dummy atoms [n*] count as 'R<n>') and net charge."""
    m = Chem.AddHs(mol_from(smiles))
    c = Counter()
    for a in m.GetAtoms():
        c[f"R{a.GetIsotope()}" if a.GetAtomicNum() == 0 else a.GetSymbol()] += 1
    return c, Chem.GetFormalCharge(m)


def formula(smiles):
    c, q = composition(smiles)
    order = ["C", "H"] + sorted(k for k in c if k not in ("C", "H"))
    s = "".join(f"{k}{c[k] if c[k] > 1 else ''}" for k in order if c.get(k))
    if q == 0:
        return s
    sign = "+" if q > 0 else "−"
    return s + (sign if abs(q) == 1 else f"{abs(q)}{sign}")


def mol_weight(smiles):
    return round(Descriptors.MolWt(mol_from(smiles)), 2)


def _smiles_of(db, name):
    if name in db.cofactors:
        return db.cofactors[name]["smiles"]
    return db.molecules[name]["smiles"]


def side(db, parts, cofactors, extra):
    total, charge = Counter(), 0
    items = [(p["id"], p["n"]) for p in parts] + [(c, 1) for c in cofactors] + [(e, 1) for e in extra]
    for name, n in items:
        c, q = composition(_smiles_of(db, name))
        for k, v in c.items():
            total[k] += v * n
        charge += q * n
    return total, charge


def check_balance(db, r):
    left, ql = side(db, r["substrates"], r["cofactors_in"], r.get("extra_in", []))
    right, qr = side(db, r["products"], r["cofactors_out"], r.get("extra_out", []))
    problems = []
    for k in sorted(set(left) | set(right)):
        if left[k] != right[k]:
            problems.append(f"{k}: 左 {left[k]} ≠ 右 {right[k]}")
    if ql != qr:
        problems.append(f"電荷: 左 {ql:+d} ≠ 右 {qr:+d}")
    return problems


# ---------------------------------------------------------------- chain ↔ SMILES

def _element(markup, last):
    """'HO' → 'O' (last=True), 'OH' → 'O', 'H' → 'H', 'H_3N^+' → 'N'."""
    p = plain(markup)
    syms = []
    i = 0
    while i < len(p):
        if p[i].isascii() and p[i].isupper():          # skips Ⓟ
            s = p[i] + (p[i + 1] if i + 1 < len(p) and p[i + 1].islower() else "")
            syms.append(s)
            i += len(s)
        else:
            i += 1
    if not syms:
        return None
    heavy = [s for s in syms if s != "H"]
    if heavy:
        return heavy[-1] if last else heavy[0]
    return "H"


def _row_element(row):
    p = plain(row["c"])
    a = row.get("anchor", anchor_index(row["c"]))
    return p[a] + (p[a + 1] if a + 1 < len(p) and p[a + 1].islower() else "")


def _drawn_h(row):
    """H atoms drawn on the row's chain atom: 'CH_2OH' → 2, '^+NH_3' → 3, 'HC' → 1, 'H_3N' → 3,
    plus any substituent that is just 'H'."""
    p = plain(row["c"])
    a = row.get("anchor", anchor_index(row["c"]))
    n = 0
    rest = p[a + 1:]
    if rest[:1].islower():                              # two-letter element: 'SeH' → H after 'Se'
        rest = rest[1:]
    after = re.match(r"H(\d?)", rest)
    if after:
        n += int(after.group(1) or 1)
    before = re.search(r"(?:^|[^A-Za-z])H(\d?)$", p[:a])
    if before:
        n += int(before.group(1) or 1)
    fork = row.get("fork") or row.get("dfork") or {}
    sides = [row.get(k) for k in ("l", "r", "u", "d")] + [fork.get(k) for k in ("l", "r")]
    return n + sum(1 for t in sides if t == "H")


def check_chain(m):
    """Check the default `chain` and every alternative layout in `chains`."""
    layouts = ([("chain", m["chain"])] if m.get("chain") else []) + \
              [(f"chains.{k}", v) for k, v in (m.get("chains") or {}).items()]
    out = []
    for label, rows in layouts:
        out += [p if label == "chain" else f"{label}: {p}" for p in check_rows(m["smiles"], rows)]
    return out


def check_rows(smiles, rows):
    mol = Chem.AddHs(mol_from(smiles))
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    ring = rows[-1].get("ring") if rows else None
    if ring:
        rows = rows[:-1]
    if any(r.get("ring") for r in rows):
        return ["chain: ring 只能放在最後一列"]

    def n_h(idx):
        return sum(1 for n in mol.GetAtomWithIdx(idx).GetNeighbors() if n.GetSymbol() == "H")

    want = [_drawn_h(r) for r in rows]
    els = [_row_element(r) for r in rows]

    def sym(atom):
        """A row written 'R' is the rest of the molecule: a dummy atom ([8*]) in the SMILES."""
        return "R" if atom.GetAtomicNum() == 0 else atom.GetSymbol()

    def phospho(idx):
        """Does this chain atom carry an O–P (a phosphate ester)?"""
        return any(o.GetSymbol() == "O" and any(x.GetSymbol() == "P" for x in o.GetNeighbors())
                   for o in mol.GetAtomWithIdx(idx).GetNeighbors())

    # a row written 'CH_2O{P}' must sit on the phosphorylated carbon and 'CH_2OH' on a free one;
    # otherwise a symmetric chain (glycerol 3-phosphate) also matches read upside down
    has_p = [("{P}" in r["c"]) if "O" in plain(r["c"]) else None for r in rows]

    # the chain = a simple path of heavy atoms whose element and H count match the drawing row by
    # row (branches such as the COO⁻ of citrate are substituents, not rows); any atom may start
    # it, so one SMILES serves several layouts of the same molecule
    def walk(path):
        if len(path) == len(rows):
            yield path
            return
        k = len(path)
        for nb in mol.GetAtomWithIdx(path[-1]).GetNeighbors():
            j = nb.GetIdx()
            if sym(nb) == els[k] and j not in path and n_h(j) == want[k] and \
                    (has_p[k] is None or phospho(j) == has_p[k]):
                yield from walk(path + [j])

    starts = [a.GetIdx() for a in mol.GetAtoms() if sym(a) == els[0] and n_h(a.GetIdx()) == want[0]
              and (has_p[0] is None or phospho(a.GetIdx()) == has_p[0])]
    paths = [p for c in starts for p in walk([c])]
    if not paths:
        return [f"chain 的 {len(rows)} 列（主原子 {els}、各列 H 數 {want}）在 SMILES 中找不到對應的鏈（結構有誤）"]
    stereo = {i for i, _ in Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)}
    results = []
    for path in paths:
        tail = None
        if ring:                                   # ring attachment atom = the neighbour below the last row
            tail = next((n.GetIdx() for n in mol.GetAtomWithIdx(path[-1]).GetNeighbors()
                         if n.IsInRing() and n.GetIdx() not in path), None)
        r = _check_path(mol, rows, path, stereo, tail)
        if ring and not r:
            r = _check_ring(mol, ring, path)
        results.append(r)
    return next((r for r in results if not r), results[0])


def _check_ring(mol, ring_smiles, path):
    """The ring fragment must exist in the molecule, attached ([*]) to the last chain atom."""
    q = Chem.MolFromSmarts(ring_smiles)
    if q is None:
        return [f"ring：無法解析 {ring_smiles}"]
    dummy = next(a.GetIdx() for a in q.GetAtoms() if a.GetAtomicNum() == 0)
    plain_mol = Chem.RemoveHs(mol)
    # map the H-free match back: atom indices are preserved for heavy atoms by RemoveHs
    for match in plain_mol.GetSubstructMatches(q):
        if match[dummy] == path[-1]:
            return []
    return [f"ring：{ring_smiles} 沒有接在主鏈最後一個原子上（或環的結構與 SMILES 不符）"]


def _substituents(row):
    """Heavy atoms the drawing bonds to this row's carbon ('HO' → O, 'COO^-' → C, 'S-CoA' → S),
    and the subset drawn with a double bond."""
    out, double = Counter(), Counter()
    fork = row.get("fork") or {}
    dfork = row.get("dfork") or {}
    for text, last, order in ((row.get("l"), True, row.get("lb", 1)), (row.get("r"), False, row.get("rb", 1)),
                              (row.get("u"), False, row.get("ub", 1)), (row.get("d"), False, row.get("db", 1)),
                              (fork.get("l"), True, fork.get("lb", 1)), (fork.get("r"), False, fork.get("rb", 1)),
                              (dfork.get("l"), True, dfork.get("lb", 1)), (dfork.get("r"), False, dfork.get("rb", 1))):
        el = _element(text, last) if text else None
        if el and el != "H":
            out[el] += 1
            if order == 2:
                double[el] += 1
    return out, double


def _check_path(mol, rows, path, stereo, tail=None):
    problems = []
    on_path = set(path)
    for i, row in enumerate(rows):
        atom = mol.GetAtomWithIdx(path[i])
        # every substituent drawn must really be bonded there (catches O= written as HO–, etc.)
        side = [b for b in atom.GetBonds() if b.GetOtherAtomIdx(path[i]) not in on_path
                and b.GetOtherAtom(atom).GetSymbol() != "H"]
        have = Counter(b.GetOtherAtom(atom).GetSymbol() for b in side)
        have2 = Counter(b.GetOtherAtom(atom).GetSymbol() for b in side if b.GetBondType() == Chem.BondType.DOUBLE)
        want, want2 = _substituents(row)
        missing = want - have
        if missing:
            problems.append(f"chain 第 {i} 列：圖上畫了 {dict(missing)}，但 SMILES 中此碳沒有接這些原子")
            continue
        # a centre written as COO⁻ / CH_2OH already contains its heteroatoms (and their double
        # bonds) inside the text, so only rows whose centre is plain C/CH_n get the double-bond check
        el = _row_element(row)
        condensed = any(ch.isascii() and ch.isupper() and ch not in ("H", el[0]) for ch in plain(row["c"])) \
            or "\u24C5" in plain(row["c"])
        if not condensed and want2 != have2:
            problems.append(f"chain 第 {i} 列：圖上的雙鍵 {dict(want2) or '無'} 與 SMILES 的 {dict(have2) or '無'} 不符")
            continue
        if path[i] not in stereo:
            continue
        cip = atom.GetPropsAsDict().get("_CIPCode")
        # the atom below: the next row, or for the last row the ring hung under it (norepinephrine)
        below = path[i + 1] if i < len(rows) - 1 else tail
        if not (row.get("l") and row.get("r")) or i == 0 or below is None:
            continue
        ex = {path[i - 1], below}
        others = [n for n in atom.GetNeighbors() if n.GetIdx() not in ex]
        if len(others) != 2:
            continue
        el_l, el_r = _element(row["l"], True), _element(row["r"], False)
        if el_l == el_r:
            continue
        try:
            left = next(n for n in others if n.GetSymbol() == el_l)
            right = next(n for n in others if n.GetSymbol() == el_r and n is not left)
        except StopIteration:
            problems.append(f"chain 第 {i} 列：左右取代基 {row['l']}/{row['r']} 與 SMILES 不符")
            continue
        if cip is None:
            problems.append(f"chain 第 {i} 列：此碳是立體中心，但 SMILES 沒有寫立體化學（@/@@）")
            continue
        # place the four neighbours in Fischer geometry and let RDKit read the configuration
        test = Chem.Mol(mol)
        conf = Chem.Conformer(test.GetNumAtoms())
        for a in range(test.GetNumAtoms()):
            conf.SetAtomPosition(a, Point3D(5 + a * 0.37, 7 + a * 0.11, 3 + a * 0.53))
        pos = {path[i]: (0, 0, 0), path[i - 1]: (0, 1, -0.7), below: (0, -1, -0.7),
               left.GetIdx(): (-1, 0, 0.7), right.GetIdx(): (1, 0, 0.7)}
        for a, (x, y, z) in pos.items():
            conf.SetAtomPosition(a, Point3D(x, y, z))
        test.RemoveAllConformers()
        test.AddConformer(conf)
        Chem.AssignStereochemistryFrom3D(test)
        drawn = test.GetAtomWithIdx(path[i]).GetPropsAsDict().get("_CIPCode")
        if drawn != cip:
            problems.append(f"chain 第 {i} 列：Fischer 圖是 {drawn}，但 SMILES 是 {cip}（左右放反了，或 SMILES 的 @/@@ 寫錯）")
    return problems + _ez_problems(mol, rows, path)


def _ez_problems(mol, rows, path):
    """A vertical C=C (`down: 2`) with one side group on each carbon: groups drawn on the same
    side (both l or both r) are cis. Compare with the E/Z of the SMILES (cis-aconitate)."""
    out = []
    for i in range(len(rows) - 1):
        if rows[i].get("down") != 2:
            continue
        a, b = path[i], path[i + 1]
        bond = mol.GetBondBetweenAtoms(a, b)
        if bond is None or bond.GetStereo() == Chem.BondStereo.STEREONONE:
            continue
        sides = []
        for k, idx, other in ((i, a, b), (i + 1, b, a)):
            row = rows[k]
            if bool(row.get("l")) == bool(row.get("r")):
                break
            text = row.get("l") or row.get("r")
            el = _element(text, bool(row.get("l")))
            on = {path[k - 1]} if k > 0 else set()
            if k + 1 < len(path):
                on.add(path[k + 1])
            cand = [n.GetIdx() for n in mol.GetAtomWithIdx(idx).GetNeighbors()
                    if n.GetIdx() not in on and n.GetIdx() != other and n.GetSymbol() == el]
            if len(cand) != 1:
                break
            sides.append((cand[0], "l" if row.get("l") else "r"))
        if len(sides) != 2:
            continue
        from rdkit.Chem import AllChem
        m2 = Chem.Mol(mol)
        AllChem.Compute2DCoords(m2)
        c = m2.GetConformer()
        pa, pb = c.GetAtomPosition(a), c.GetAtomPosition(b)

        def side_of(j):
            q = c.GetAtomPosition(j)
            return (pb.x - pa.x) * (q.y - pa.y) - (pb.y - pa.y) * (q.x - pa.x) > 0
        real_cis = side_of(sides[0][0]) == side_of(sides[1][0])
        drawn_cis = sides[0][1] == sides[1][1]
        if real_cis != drawn_cis:
            out.append(f"chain 第 {i}–{i + 1} 列的雙鍵：圖上 {rows[i].get('l') or rows[i].get('r')} 與 "
                       f"{rows[i + 1].get('l') or rows[i + 1].get('r')} 是{'順式' if drawn_cis else '反式'}，"
                       f"但 SMILES 是{'順式' if real_cis else '反式'}（把其中一個換到另一側）")
    return out


def check_free(m):
    from .draw.free import free_composition
    if not m.get("free") or m["free"].get("schematic"):
        return []
    drawn, q = free_composition(m["free"])
    real, rq = composition(m["smiles"])
    real = Counter({k: v for k, v in real.items() if not k.startswith("R")})
    problems = []
    if drawn != real:
        diff = {k: (drawn[k], real[k]) for k in set(drawn) | set(real) if drawn[k] != real[k]}
        problems.append(f"free：圖上的原子數與 SMILES 不符（圖, SMILES）= {diff}")
    if q != rq:
        problems.append(f"free：圖上的電荷 {q:+d} 與 SMILES 的 {rq:+d} 不符")
    return problems


def _label_element(text):
    if isinstance(text, list):                          # [text, anchor index] (haworth subs)
        return plain(text[0])[text[1]]
    p = plain(text)
    if p in ("Base", "B", "R"):                         # generic base → dummy atom in the SMILES
        return "*"
    return p[anchor_index(text)] if p else None


def _cip_from(mol, center, placed):
    """CIP code of `center` when the atoms in `placed` sit at the given 3D points."""
    test = Chem.Mol(mol)
    conf = Chem.Conformer(test.GetNumAtoms())
    for a in range(test.GetNumAtoms()):
        conf.SetAtomPosition(a, Point3D(5 + a * 0.37, 7 + a * 0.11, 3 + a * 0.53))
    for a, (x, y, z) in placed.items():
        conf.SetAtomPosition(a, Point3D(x, y, z))
    test.RemoveAllConformers()
    test.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(test)
    return test.GetAtomWithIdx(center).GetPropsAsDict().get("_CIPCode")


def _unit_problems(mol, u, carbons, ring, ring3, subs3, sides):
    """Problems when the ring carbons `carbons` (position 1 …) carry the drawn substituents.
    ring3 / subs3: 3D points of the drawing (ring keys "O", 1…n; subs keyed (pos, side));
    sides: the two side names ("up", "down") or ("ax", "eq")."""
    problems = []
    base_el = None
    if u.get("base"):
        q = Chem.MolFromSmiles(u["base"])
        d = next(a for a in q.GetAtoms() if a.GetAtomicNum() == 0)
        base_el = d.GetNeighbors()[0].GetSymbol()
    start = u.get("start", 1)
    ring_o = next(i for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == "O")
    pts = {ring_o: ring3["O"]}
    for k, idx in enumerate(carbons, start=1):
        pts[idx] = ring3[k]
    subs = {int(c) - start + 1: v for c, v in (u.get("subs") or {}).items()}
    for k, idx in enumerate(carbons, start=1):
        sub = subs.get(k) or {}
        name = f"C{k + start - 1}"
        atom = mol.GetAtomWithIdx(idx)
        exo = [n for n in atom.GetNeighbors() if n.GetIdx() not in ring]
        drawn = {}
        for side in sides + ("out",):
            t = sub.get(side)
            if t:
                drawn[side] = base_el if (t == "base" and base_el) else ("O" if t == "link" else _label_element(t))
        have = Counter(n.GetSymbol() for n in exo if n.GetSymbol() != "H")
        want = Counter(e for e in drawn.values() if e and e != "H")
        if have != want and (len(drawn) >= 2 or "out" in drawn or not want <= have):
            problems.append(f"{name}：圖上取代基 {dict(want) or '只有 H'}，SMILES 是 {dict(have) or '只有 H'}")
            continue
        cip = atom.GetPropsAsDict().get("_CIPCode")
        if not cip:
            continue
        el_a, el_b = drawn.get(sides[0], "H"), drawn.get(sides[1], "H")
        if el_a == el_b:
            problems.append(f"{name}：此碳是立體中心，但兩側畫了相同的取代基")
            continue
        try:
            na = next(n for n in exo if n.GetSymbol() == el_a)
            nb = next(n for n in exo if n.GetSymbol() == el_b and n.GetIdx() != na.GetIdx())
        except StopIteration:
            problems.append(f"{name}：取代基 {el_a}/{el_b} 與 SMILES 不符")
            continue
        placed = {idx: pts[idx], na.GetIdx(): subs3[(k, sides[0])], nb.GetIdx(): subs3[(k, sides[1])]}
        for n in atom.GetNeighbors():
            if n.GetIdx() in pts:
                placed[n.GetIdx()] = pts[n.GetIdx()]
        drawn_cip = _cip_from(mol, idx, placed)
        if drawn_cip != cip:
            problems.append(f"{name}：圖是 {drawn_cip}，但 SMILES 是 {cip}（上下／軸向放反了，或 α/β、D/L 弄錯）")
    if u.get("base"):
        pos = next((int(c) - start + 1 for c, v in (u.get("subs") or {}).items() if "base" in v.values()), None)
        q = Chem.MolFromSmarts(u["base"])
        d = next(a.GetIdx() for a in q.GetAtoms() if a.GetAtomicNum() == 0)
        if pos is None or not any(mt[d] == carbons[pos - 1] for mt in mol.GetSubstructMatches(q)):
            problems.append(f"鹼基 {u['base']} 沒有接在 C{pos}′ 上（或鹼基結構與 SMILES 不符）")
    return problems


def _candidates(mol, n):
    """Every (ring atom set, carbons in order) for rings of size n with exactly one O."""
    out = []
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != n or sum(mol.GetAtomWithIdx(i).GetSymbol() == "O" for i in ring) != 1:
            continue
        o = next(i for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == "O")
        k = ring.index(o)
        seq = list(ring[k:] + ring[:k])
        out += [(frozenset(ring), seq[1:]), (frozenset(ring), seq[1:][::-1])]
    return out


def _prepared(m):
    mol = Chem.AddHs(Chem.RemoveHs(mol_from(m["smiles"])))
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol


def check_haworth(m):
    spec = m.get("haworth")
    if not spec:
        return []
    from .draw.haworth import RINGS, haworth_3d, units_of
    mol = _prepared(m)
    units = units_of(spec)
    tag = "haworth" + ("" if len(units) == 1 else " 第 {} 環")
    fits, problems = [], []
    for ui, u in enumerate(units):
        n = len(RINGS[u.get("ring", "furanose")])
        ring3, subs3 = haworth_3d(u)
        good, best = [], None
        for ring, carbons in _candidates(mol, n):
            probs = _unit_problems(mol, u, carbons, ring, ring3, subs3, ("up", "down"))
            if not probs:
                good.append((ring, carbons))
            elif best is None or len(probs) < len(best):
                best = probs
        fits.append(good)
        if not good:
            pre = tag.format(ui + 1)
            problems += [f"{pre}：{p}" for p in (best or [f"SMILES 中找不到含一個 O 的 {n} 員環"])]
    if problems:
        return problems
    # each unit needs its own ring, and every glycosidic O must really join the two carbons
    links = spec.get("links") or []

    def linked(chosen):
        for ln in links:
            (ua, ca, _), (ub, cb, _) = ln["from"], ln["to"]
            a = chosen[ua][int(ca) - units[ua].get("start", 1)]
            b = chosen[ub][int(cb) - units[ub].get("start", 1)]
            oa = {n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetSymbol() == "O"}
            ob = {n.GetIdx() for n in mol.GetAtomWithIdx(b).GetNeighbors() if n.GetSymbol() == "O"}
            if not oa & ob:
                return ln
        return None

    bad = None
    def assign(i, used, chosen):
        nonlocal bad
        if i == len(units):
            miss = linked(chosen)
            if miss is None:
                return True
            bad = miss
            return False
        return any(assign(i + 1, used | {ring}, chosen + [carbons])
                   for ring, carbons in fits[i] if ring not in used)

    if assign(0, frozenset(), []):
        return []
    if bad is None:
        return ["haworth：兩個環對到 SMILES 中的同一個環"]
    return [f"haworth：連結 {bad['from']} → {bad['to']} 在 SMILES 中沒有共用的 O（醣苷鍵位置錯）"]


def check_chair(m):
    spec = m.get("chair")
    if not spec:
        return []
    from .draw.haworth import chair_validation_3d
    mol = _prepared(m)
    ring3, subs3 = chair_validation_3d(spec)
    best = None
    for ring, carbons in _candidates(mol, 6):
        probs = _unit_problems(mol, spec, carbons, ring, ring3, subs3, ("ax", "eq"))
        if not probs:
            return []
        if best is None or len(probs) < len(best):
            best = probs
    return [f"chair：{p}" for p in (best or ["SMILES 中找不到吡喃環"])]


def _flat(smiles):
    """Canonical SMILES without tetrahedral stereo (cis/trans double bonds are kept)."""
    mol = mol_from(smiles)
    for a in mol.GetAtoms():
        a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    return Chem.MolToSmiles(mol)


def check_lipid(m):
    """The lipid spec builds its own SMILES; it must be the molecule (every cis/trans bond
    included). Tetrahedral centres (glycerol C2, sphingosine) are not drawn, so not compared."""
    if not m.get("lipid"):
        return []
    from .draw.lipid import build
    _free, smi, problems, _n = build(m["lipid"])
    out = [f"lipid：{p}" for p in problems]
    if _flat(smi) != _flat(m["smiles"]):
        out.append(f"lipid：圖（{smi}）與 SMILES 不是同一個分子（鏈長、雙鍵位置或 cis/trans 不符）")
    return out


def check_ring_molecule(m):
    if not m.get("ring"):
        return []
    spec = m["ring"] if isinstance(m["ring"], dict) else {}
    if "smiles" in spec and Chem.CanonSmiles(spec["smiles"]) != Chem.CanonSmiles(m["smiles"]):
        return [f"ring：畫圖用的 SMILES 與分子的 SMILES 不同"]
    mol = Chem.MolFromSmiles(spec.get("smiles", m["smiles"]))
    probs = []
    for h in (m.get("decorate") or {}).get("highlight") or []:
        if h.get("smarts"):
            patt = Chem.MolFromSmarts(h["smarts"])
            if patt is None or not mol.HasSubstructMatch(patt):
                probs.append(f"decorate.highlight：SMARTS {h['smarts']} 在分子中找不到")
        if any(int(i) >= mol.GetNumAtoms() for i in h.get("atoms") or []):
            probs.append("decorate.highlight：atoms 的索引超出範圍")
    return probs


def check_cofactor_link(db, c):
    """`molecule: nad` in cofactors.yaml links the label to a structure page; the cofactor's
    changing part (dummy atoms = the rest) must really be part of that molecule."""
    m = db.molecules.get(c["molecule"])
    if m is None:
        return [f"molecule: 找不到分子 {c['molecule']}"]
    frag = max(c["smiles"].split("."), key=len)
    patt = Chem.MolFromSmarts(re.sub(r"\[\d*\*\]", "*", frag))
    if patt is None or not Chem.MolFromSmiles(m["smiles"]).HasSubstructMatch(patt):
        return [f"molecule: {c['molecule']} 不含輔因子片段 {frag}"]
    return []


_PARTS = {"c", "l", "r", "u", "d", "fork_l", "fork_r", "dfork_l", "dfork_r", "ring"}


def _frags(spec):
    """colour spec → list of text fragments it needs ([] for a whole-part colour)."""
    if isinstance(spec, (list, tuple)):
        pairs = spec if spec and isinstance(spec[0], (list, tuple)) else [spec]
        return [f for _, f in pairs]
    return []


def _text_has(text, frag):
    from .text import TextLine
    return TextLine(text).span(frag) is not None


def _check_highlight(smiles, specs, what):
    out = []
    mol = Chem.MolFromSmiles(smiles)
    for h in specs or []:
        if h.get("smarts"):
            patt = Chem.MolFromSmarts(h["smarts"])
            if patt is None or not mol.GetSubstructMatches(patt):
                out.append(f"{what}: SMARTS {h['smarts']} 沒有比對到")
        for i in h.get("atoms") or []:
            if int(i) >= mol.GetNumAtoms():
                out.append(f"{what}: 原子編號 {i} 超出範圍（共 {mol.GetNumAtoms()} 個）")
    return out


def check_decorate(db, r):
    """Colours that silently do nothing are the typical mistake: wrong key, row or fragment."""
    out = []
    decos = [("", r.get("decorate") or {})] + [(f"variants.{k}.", v.get("decorate") or {})
                                                for k, v in (r.get("variants") or {}).items()]
    parts = {p["id"]: p for p in r["substrates"] + r["products"]}
    cofs = set(r.get("cofactors_in") or []) | set(r.get("cofactors_out") or [])
    for pre, deco in decos:
        for key, d in deco.items():
            what = f"decorate.{pre}{key}"
            if key in cofs and key not in parts:
                from .compose.reaction import _cof_label
                label = _cof_label(db, key)
                for row, spec in ((d or {}).get("colors") or {}).items():
                    for f in _frags(spec.get("c") if isinstance(spec, dict) else spec):
                        if not _text_has(label, f):
                            out.append(f"{what}: 「{label}」裡找不到 {f}")
                continue
            if key not in parts:
                out.append(f"{what}: 不是這個反應的反應物、產物或輔因子")
                continue
            m = db.molecules[key]
            view = parts[key].get("view") or m["default_view"]
            colors = (d or {}).get("colors") or {}
            if view == "chain":
                rows = (m.get("chains") or {}).get(parts[key]["chain"]) if parts[key].get("chain") else m["chain"]
                for i, spec in colors.items():
                    i = int(i)
                    if i >= len(rows):
                        out.append(f"{what}: 第 {i} 列不存在（共 {len(rows)} 列）")
                        continue
                    row = rows[i]
                    for part, ps in (spec.items() if isinstance(spec, dict) else []):
                        if part not in _PARTS:
                            out.append(f"{what}: 第 {i} 列沒有「{part}」這種部位")
                            continue
                        if part == "ring":
                            if not row.get("ring"):
                                out.append(f"{what}: 第 {i} 列不是環")
                            elif isinstance(ps, list):
                                out += _check_highlight(row["ring"], ps, f"{what} 第 {i} 列")
                            continue
                        src = row.get("fork" if part.startswith("fork") else "dfork" if part.startswith("dfork") else None) or row
                        text = src.get(part.split("_")[-1] if "fork" in part else part)
                        if text is None:
                            out.append(f"{what}: 第 {i} 列沒有 {part}")
                            continue
                        for f in _frags(ps):
                            if not _text_has(text, f):
                                out.append(f"{what}: 第 {i} 列 {part}「{text}」裡找不到 {f}")
            elif view == "formula":
                for spec in colors.values():
                    for f in _frags(spec.get("c") if isinstance(spec, dict) else spec):
                        if not _text_has(m["formula"], f):
                            out.append(f"{what}: 「{m['formula']}」裡找不到 {f}")
            elif view in ("ring", "ring_mirror"):
                ring = m["ring"] if isinstance(m["ring"], dict) else {}
                out += _check_highlight(ring.get("smiles", m["smiles"]), d.get("highlight"), what)
            elif view == "skeletal" and d.get("smarts"):
                out += _check_highlight(m["smiles"], [{"smarts": d["smarts"]}], what)
            elif view == "haworth":
                units = m["haworth"].get("units") or [m["haworth"]]
                for k, spec in colors.items():
                    mk = re.fullmatch(r"(?:u(\d+)_)?([udo])(\d+)", str(k))
                    if not mk:
                        if not re.fullmatch(r"(?:u\d+_)?(c\d+|link\d+|o)", str(k)):
                            out.append(f"{what}: 不認得的 Haworth 鍵名 {k}")
                        continue
                    u = units[int(mk.group(1) or 0)]
                    pos = int(mk.group(3)) + u.get("start", 1) - 1
                    sub = (u.get("subs") or {}).get(pos) or (u.get("subs") or {}).get(str(pos)) or {}
                    text = sub.get({"u": "up", "d": "down", "o": "out"}[mk.group(2)])
                    if text is None:
                        out.append(f"{what}: Haworth {k} 位置沒有取代基")
                        continue
                    for f in _frags(spec):
                        if not _text_has(text, f):
                            out.append(f"{what}: Haworth {k}「{text}」裡找不到 {f}")
                for u in units:
                    if u.get("base") and (d.get("base_colors") or d.get("base_highlight")):
                        out += _check_highlight(u["base"], [{"atoms": list((d.get("base_colors") or {}).keys())}]
                                                + list(d.get("base_highlight") or []), what)
            elif view in ("lipid", "free"):
                pass                                   # atom ids are checked by the drawer itself
    return out


def validate(db):
    """Return (errors, warnings, info) where info[(kind,id)] has formula/MW."""
    errors, warnings, info = list(db.errors), [], {}
    for m in db.molecules.values():
        try:
            info[("molecule", m["id"])] = {"formula": formula(m["smiles"]), "mw": mol_weight(m["smiles"])}
            errors += [f"molecule {m['id']}: {p}" for p in check_chain(m) + check_free(m)
                                                               + check_haworth(m) + check_chair(m) + check_ring_molecule(m)
                                                               + check_lipid(m)]
        except Exception as e:
            errors.append(f"molecule {m['id']}: {e}")
        if m.get("ring") and isinstance(m["ring"], dict) and not m["ring"].get("explicit"):
            from .draw.ring import label_overlaps
            hits = label_overlaps(m["ring"].get("smiles", m["smiles"]), m["ring"])
            if hits:
                warnings.append(f"molecule {m['id']}: ring 圖上有 {len(hits)} 組原子標籤重疊 {hits}"
                                "（試試 coordgen: false、seed: 1…、rotate）")
    for name, c in db.cofactors.items():
        try:
            mol_from(c["smiles"])
            if c.get("molecule"):
                errors += [f"cofactor {name}: {p}" for p in check_cofactor_link(db, c)]
        except Exception as e:
            errors.append(f"cofactor {name}: {e}")
    for r in db.reactions.values():
        try:
            errors += [f"reaction {r['id']}: {p}" for p in check_decorate(db, r)]
        except Exception as e:
            errors.append(f"reaction {r['id']}: decorate 檢查失敗 {e!r}")
        if r.get("balance_check", True) is False:
            warnings.append(f"reaction {r['id']}: 已關閉平衡檢查")
            continue
        try:
            for p in check_balance(db, r):
                errors.append(f"reaction {r['id']}: 不平衡 → {p}")
        except Exception as e:
            errors.append(f"reaction {r['id']}: {e}")
    for p in db.pathways.values():
        scenes = ([p["scene"]] if p["layout"] == "scene" else []) + \
                 [v["scene"] for v in (p.get("variants") or {}).values() if v.get("scene")]
        for sc in scenes:
            from .draw.scene import check_scene
            e, w = check_scene(sc)
            errors += [f"pathway {p['id']}: {q}" for q in e]
            warnings += [f"pathway {p['id']}: {q}" for q in w]
        for n in p.get("nodes") or []:                     # node highlights: SMARTS must parse and hit the molecule
            m = db.molecules.get(n.get("ref"))
            for h in ((n.get("decorate") or {}).get("highlight") or []):
                if not m or not h.get("smarts"):
                    continue
                patt = Chem.MolFromSmarts(h["smarts"])
                if patt is None:
                    errors.append(f"pathway {p['id']}: 節點 {n['id']} 的 SMARTS 無法解析 {h['smarts']!r}")
                elif not Chem.MolFromSmiles(m["smiles"]).HasSubstructMatch(patt):
                    errors.append(f"pathway {p['id']}: 節點 {n['id']} 的 SMARTS {h['smarts']!r} 在 {m['id']} 找不到")
        if p["layout"] != "plot":
            continue
        from .draw.plot import check_plot
        specs = list(p.get("panels") or [p.get("plot")])
        for v in (p.get("variants") or {}).values():
            specs += v.get("panels") or ([{**specs[0], **v["plot"]}] if v.get("plot") else [])
        for spec in specs:
            try:
                errors += [f"pathway {p['id']}: {q}" for q in check_plot(spec)]
            except Exception as e:
                errors.append(f"pathway {p['id']}: 作圖設定錯誤 {e!r}")
    return errors, warnings, info
