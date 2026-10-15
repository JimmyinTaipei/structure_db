import copy
import re
import shutil
import struct
import zlib

import pytest

from sdb import export, model, review, style
from sdb.compose.reaction import compose_reaction
from sdb.draw.chain import PITCH, draw_chain
from sdb.validate import check_balance, check_chain, validate


@pytest.fixture()
def sandbox(tmp_path, monkeypatch):
    """Copy data/ into a temp dir and redirect every output path there."""
    shutil.copytree(model.DATA, tmp_path / "data")
    monkeypatch.setattr(model, "DATA", tmp_path / "data")
    monkeypatch.setattr(export, "OUT", tmp_path / "out")
    monkeypatch.setattr(review, "FILE", tmp_path / "review" / "comments.jsonl")
    return tmp_path


def test_pilot_data_is_valid():
    errors, _, _ = validate(model.load())
    assert errors == []


def test_balance_catches_missing_cofactor():
    db = model.load()
    r = copy.deepcopy(db.reactions["pdh"])
    r["cofactors_in"] = ["CoA-SH"]            # forgot NAD+
    assert check_balance(db, r)
    r["cofactors_out"] = []                   # dropping NADH too still leaves the hydride unaccounted
    assert any(p.startswith("H:") for p in check_balance(db, r))


def test_balance_catches_wrong_stoichiometry():
    db = model.load()
    r = copy.deepcopy(db.reactions["aldol_condensation"])
    r["substrates"][0]["n"] = 1
    assert any(p.startswith("C:") for p in check_balance(db, r))


def test_fischer_flip_is_detected():
    m = copy.deepcopy(model.load().molecules["xylulose_5p"])
    assert check_chain(m) == []
    m["chain"][2] = {"c": "C", "l": "H", "r": "OH"}     # C3 OH moved to the right
    assert "第 2 列" in check_chain(m)[0]


def test_smiles_order_is_free_but_structure_must_match():
    m = copy.deepcopy(model.load().molecules["pyruvate"])
    m["smiles"] = "CC(=O)C(=O)[O-]"                     # written from the CH3 end: still fine
    assert check_chain(m) == []
    m["smiles"] = "CC(O)C(=O)[O-]"                      # lactate: C2 has an H the drawing lacks
    assert check_chain(m)


@pytest.mark.parametrize("mid,rows", [("glyceraldehyde_3p", 3), ("xylulose_5p", 5), ("sedoheptulose_7p", 7)])
def test_fischer_row_count_sets_height(mid, rows):
    m = model.load().molecules[mid]
    assert len(m["chain"]) == rows
    g = draw_chain(m["chain"])
    # vertical bonds: one per gap between rows
    n_vertical = len(re.findall(r'<line x1="0.00" y1="[\d.]+" x2="0.00"', g.body))
    assert n_vertical == rows - 1
    assert g.h == pytest.approx((rows - 1) * PITCH, abs=3 * style.S)


def test_arrow_points_right():
    db = model.load()
    body = compose_reaction(db, db.reactions["pdh"]).body
    arrow = body[body.index('class="sdb-arrow"'):]
    x2 = max(float(x) for x in re.findall(r'<line x1="[\d.]+" y1="0.00" x2="([\d.]+)" y2="0.00"', arrow))
    tips = [float(m) for m in re.findall(r'L([\d.]+),0\.00 L', arrow)]
    assert x2 in tips


def _png_info(path):
    data = path.read_bytes()
    w, h = struct.unpack(">II", data[16:24])
    color_type = data[25]
    idat, i = b"", 8
    while i < len(data):
        n = struct.unpack(">I", data[i:i + 4])[0]
        if data[i + 4:i + 8] == b"IDAT":
            idat += data[i + 8:i + 8 + n]
        i += 12 + n
    raw = zlib.decompress(idat)
    first_px_alpha = raw[4]          # row 0: filter byte, then RGBA of pixel 0
    return w, h, color_type, first_px_alpha


def test_build_png_is_3x_and_transparent(sandbox):
    errors, _ = export.build(model.load())
    assert errors == []
    svg = (sandbox / "out/svg/reactions/pdh.svg").read_text()
    sw, sh = (float(v) for v in re.search(r'width="([\d.]+)" height="([\d.]+)"', svg).groups())
    w, h, ct, alpha = _png_info(sandbox / "out/png/reactions/pdh.png")
    assert (w, h) == (pytest.approx(sw * 3, abs=3), pytest.approx(sh * 3, abs=3))
    assert ct == 6 and alpha == 0                        # RGBA, corner pixel transparent
    assert (sandbox / "out/md/TCA.md").exists()
    assert (sandbox / "out/index.json").exists()


def test_review_flow(sandbox):
    db = model.load()
    path = db.files[("reaction", "pdh")]
    before = path.read_text(encoding="utf-8")
    model.set_status(db, "reaction", "pdh", "changes_requested")
    after = path.read_text(encoding="utf-8")
    assert "status: changes_requested" in after
    # only the status line changed (comments & formatting preserved)
    assert [l for l in before.splitlines() if not l.startswith("status:")] == \
           [l for l in after.splitlines() if not l.startswith("status:")]
    c = review.add("reaction", "pdh", "comment", "NADH 的箭頭太短")
    assert review.open_comments()[0]["cid"] == c["cid"]
    review.resolve(c["cid"], "已加長")
    assert review.open_comments() == []


def test_parts_split(sandbox):
    from sdb.compose.pathway import compose_pathway
    db = model.load()
    parts = export.molecule_parts(db.molecules["xylulose_5p"])
    assert {"structure", "name", "row0", "row4"} <= set(parts)
    rparts = compose_reaction(db, db.reactions["pdh"]).parts
    assert {"arrow", "enzyme", "mol_pyruvate", "mol_acetyl_coa", "name_pyruvate"} <= set(rparts)
    assert not any(k.startswith("highlight") for k in rparts)       # colour, not boxes, by default
    # the arrow part must not contain the enzyme text path
    assert rparts["enzyme"].body not in rparts["arrow"].body
    g = compose_pathway(db, db.pathways["tca_cycle"])
    assert g.body.count('data-ref="reaction:') == 9


def test_branched_chain_and_isocitrate_stereo():
    m = copy.deepcopy(model.load().molecules["isocitrate"])
    assert check_chain(m) == []                         # COO⁻ branches are substituents, not rows
    m["smiles"] = m["smiles"].replace("[C@H]", "[C@@H]", 1)
    assert any("第 2 列" in p and "SMILES 是" in p for p in check_chain(m))


MAG = style.COLORS["magenta"]


def test_part_colouring_and_bond_rule():
    rows = [{"c": "C", "l": "CH_3", "r": "S-CoA", "u": "O", "ub": 2}]
    g = draw_chain(rows, {"colors": {0: {"c": "magenta", "l": "magenta", "u": "magenta"}}})
    coloured_lines = re.findall(rf'<line [^>]*stroke="{MAG}"', g.body)
    black_lines = re.findall(rf'<line [^>]*stroke="{style.INK}"', g.body)
    assert len(coloured_lines) == 3                    # CH3–C and the two lines of C=O
    assert len(black_lines) == 1                       # C–S stays black (S-CoA is not coloured)


def test_fragment_colouring_and_phosphoryl_symbol():
    g = draw_chain([{"c": "CH_2O{P}"}], {"colors": {0: {"c": ["magenta", "{P}"]}}})
    assert 'fill="#F3C388"' in g.body                  # the Ⓟ disc keeps its orange fill
    assert f'stroke="{MAG}"' in g.body                 # …its ring turns magenta
    assert f'fill="{style.INK}"' in g.body             # …while CH2O stays black


def test_isotope_marks_drawn():
    plain_g = draw_chain([{"c": "CH_2"}, {"c": "C"}])
    marked = draw_chain([{"c": "CH_2"}, {"c": "C"}], {"marks": {0: {"c": "*"}, 1: {"c": ["tri", "2"]}}})
    orange = style.COLORS["orange"]
    assert orange not in plain_g.body
    assert marked.body.count(orange) >= 3              # *, △ and the ring + digit of ②


def test_alternative_layouts_are_validated():
    m = copy.deepcopy(model.load().molecules["oxaloacetate"])
    assert "cs" in m["chains"] and check_chain(m) == []
    m["chains"]["cs"][0] = {"c": "C", "l": "HO", "r": "COO^-"}     # wrong: that C has no OH
    assert any(p.startswith("chains.cs") for p in check_chain(m))


def test_variants_exported(sandbox):
    errors, _ = export.build(model.load())
    assert errors == []
    out = sandbox / "out"
    assert (out / "svg/reactions/citrate_synthase.isotope.svg").exists()
    assert (out / "png/pathways/tca_cycle.isotope.png").exists()
    assert (out / "parts/reactions/citrate_synthase/isotope/mol_citrate.svg").exists()
    assert "isotope" not in (out / "svg/reactions/pdh.svg").read_text()


# ------------------------------------------------------------ round 3: amino acids & tracking

def test_chain_through_heteroatoms():
    db = model.load()
    for mid in ("methionine", "arginine", "lysine", "argininosuccinate", "citrulline"):
        assert check_chain(db.molecules[mid]) == [], mid
    m = copy.deepcopy(db.molecules["methionine"])
    m["chain"][4] = {"c": "O"}                          # S drawn as O
    assert check_chain(m)


def test_ring_rows_are_checked():
    db = model.load()
    assert check_chain(db.molecules["tryptophan"]) == []
    m = copy.deepcopy(db.molecules["phenylalanine"])
    m["chain"][-1] = {"ring": "[*]c1ccc(O)cc1"}         # Tyr ring on Phe
    assert any("ring" in p for p in check_chain(m))


def test_free_drawing_composition():
    from sdb.validate import check_free
    m = copy.deepcopy(model.load().molecules["proline"])
    assert check_free(m) == []
    m["free"]["atoms"]["cg"][0] = "CH"
    assert check_free(m)


def test_dfork_substituents_checked():
    m = copy.deepcopy(model.load().molecules["valine"])
    assert check_chain(m) == []
    m["chain"][2]["dfork"]["r"] = "OH"
    assert check_chain(m)


def test_background_layers_are_separate_parts():
    parts = export.molecule_parts(model.load().molecules["methionine"])
    assert {"structure", "box_R", "box_S", "name"} <= set(parts)
    assert style.COLORS["yellow"] not in parts["structure"].body
    assert style.COLORS["yellow"] in parts["box_S"].body


def test_inline_cofactors_and_multi_span():
    db = model.load()
    body = compose_reaction(db, db.reactions["malate_dehydrogenase"]).body
    assert 'data-ref="molecule:NAD+"' in body and 'data-ref="molecule:NADH+H+"' in body
    from sdb.text import TextLine
    from sdb.draw.chain import Paint
    t = TextLine("NADH + H^+")
    _, spans = Paint({0: {"c": [["magenta", "H"], ["magenta", "H^+"]]}}).text(0, "c", t)
    assert [(a, b) for a, b, _ in spans] == [(3, 4), (7, 9)]


def test_edges_branching_from_an_edge():
    from sdb.compose.pathway import compose_pathway
    db = model.load()
    body = compose_pathway(db, db.pathways["urea_cycle"]).body
    assert body.count('data-ref="reaction:') == 5


# ---------------------------------------------------------------- round 4: Haworth, bases, text versions

def _haworth_problems(mid, pos, sub):
    from sdb.validate import check_haworth
    m = copy.deepcopy(model.load().molecules[mid])
    m["haworth"]["subs"][pos] = sub
    return check_haworth(m)


def test_haworth_stereo_alpha_beta_and_sugar():
    from sdb.validate import check_haworth
    db = model.load()
    assert check_haworth(db.molecules["amp"]) == []
    assert any("C1" in p for p in _haworth_problems("amp", 1, {"up": "H", "down": "base"}))       # α-anomer
    assert any("C2" in p for p in _haworth_problems("amp", 2, {"up": "OH", "down": "H"}))         # arabinose
    assert any("C2" in p for p in _haworth_problems("amp", 2, {"up": "H", "down": "H"}))          # drawn deoxy
    assert _haworth_problems("dtmp", 2, {"up": "H", "down": "OH"})                                 # dTMP drawn ribo
    assert check_haworth(db.molecules["prpp"]) == []                                               # α-PRPP


def test_haworth_base_must_match():
    from sdb.validate import check_haworth
    m = copy.deepcopy(model.load().molecules["amp"])
    m["haworth"]["base"] = "[*]n1cnc2c(=O)[nH]cnc21"                                               # hypoxanthine on AMP
    assert any("鹼基" in p for p in check_haworth(m))


def test_text_variants_and_pathway_layers():
    from sdb.compose.pathway import compose_pathway
    db = model.load()
    p = db.pathways["tca_cycle"]
    assert {"text", "text_full"} <= set(p["variants"])
    assert "text_full" not in db.pathways["urea_cycle"]["variants"]          # already name-only
    assert "variants" not in db.pathways["amino_acids"] or "text" not in db.pathways["amino_acids"]["variants"]
    full = compose_pathway(db, p, "text_full")
    assert {"path", "enzymes", "cofactors"} <= set(full.parts)
    plain = compose_pathway(db, p, "text")
    assert "enzymes" not in plain.parts and "cofactors" not in plain.parts
    assert plain.w < compose_pathway(db, p).w                                  # names are smaller than structures


def test_tag_tree_in_index(sandbox):
    db = model.load()
    export.build(db)
    import json
    idx = json.loads((sandbox / "out" / "index.json").read_text(encoding="utf-8"))
    assert "nonpolar" in idx["tag_tree"]["tree"]["amino_acid"]
    assert idx["tag_tree"]["labels"]["purine"] == "嘌呤"


# ---------------------------------------------------------------- round 5: sugars, bases, reaction text

def test_reaction_text_variants():
    from sdb.compose.reaction import compose_reaction
    db = model.load()
    r = db.reactions["malate_dehydrogenase"]
    assert {"text", "text_full"} <= set(r["variants"])
    full = compose_reaction(db, r, "text_full")
    assert "enzyme" in full.parts and "Malate" in model.load().molecules["malate"]["name"]
    plain = compose_reaction(db, r, "text")
    assert "enzyme" not in plain.parts
    assert plain.h < compose_reaction(db, r).h                     # names are flatter than structures


def test_chair_axial_equatorial_checked():
    from sdb.validate import check_chair
    m = copy.deepcopy(model.load().molecules["b_glucopyranose"])
    assert check_chair(m) == []
    m["chair"]["subs"][1] = {"ax": "OH", "eq": "H"}                # α drawing, β SMILES
    assert any("C1" in p for p in check_chair(m))


def test_disaccharide_links_checked():
    from sdb.validate import check_haworth
    db = model.load()
    for mid in ("maltose", "lactose", "sucrose", "cellobiose"):
        assert check_haworth(db.molecules[mid]) == [], mid
    m = copy.deepcopy(db.molecules["maltose"])
    m["haworth"]["links"][0]["to"] = [1, 3, "down"]                # α1→3 instead of α1→4
    m["haworth"]["units"][1]["subs"][3] = {"up": "OH", "down": "link"}
    assert check_haworth(m)


def test_fischer_family_and_bases():
    from sdb.compose.molecule import structure, views_of
    db = model.load()
    for mid in ("d_ribose", "d_galactose", "d_fructose", "sorbitol", "l_glucose"):
        assert check_chain(db.molecules[mid]) == [], mid
    ade = db.molecules["adenine"]
    assert {"ring", "ring_mirror"} <= set(views_of(ade))
    a, b = structure(ade, "ring"), structure(ade, "ring_mirror")
    assert abs(a.w - b.w) < style.S and a.body != b.body                 # same card, mirrored


# ---------------------------------------------------------------- round 6: lipids, base cards, topics

def test_lipid_spec_matches_smiles():
    from sdb.validate import check_lipid
    db = model.load()
    for mid in ("oleate", "elaidate", "arachidonate", "phosphatidylcholine", "sphingomyelin", "triacylglycerol"):
        assert check_lipid(db.molecules[mid]) == [], mid
    m = copy.deepcopy(db.molecules["oleate"])
    m["lipid"]["chains"][0] = {"n": 18, "trans": [9]}                 # drawn trans, SMILES cis
    assert check_lipid(m)
    m = copy.deepcopy(db.molecules["phosphatidylcholine"])
    m["lipid"]["head"] = "ethanolamine"                                 # wrong head group
    assert check_lipid(m)


def test_zigzag_cis_is_drawn_cis():
    from sdb.draw.lipid import check_geometry, zigzag
    for cis in ([9], [9, 12], [5, 8, 11, 14], [4, 7, 10, 13, 16, 19]):
        dbl = {k: "cis" for k in cis}
        pos, _ = zigzag(0, 0, 1, 22, dbl)
        assert check_geometry(pos, dbl) == []
        pos, _ = zigzag(0, 0, 1, 22, dbl, heading=60.0, theta=90.0)
        assert check_geometry(pos, dbl) == []


def test_base_cards_and_steroid_pose():
    from sdb.compose.molecule import structure, views_of
    db = model.load()
    ade = db.molecules["adenine"]
    assert views_of(ade)[:2] == ["ring", "ring_mirror"] and not ade["ring"].get("explicit")
    g = structure(db.molecules["cholesterol"], "ring")
    assert g.w > g.h * 0.8                                              # lying steroid, not standing


def test_topics_follow_reactions_and_pathways():
    db = model.load()
    t = export.topics_of(db)
    assert {"glycolysis", "TCA"} <= t[("molecule", "pyruvate")]
    assert "TCA" in t[("molecule", "malate")]
    assert export.topic_folder("acidic", export.tag_tree()) == "胺基酸/酸性"


def test_coenzyme_highlight_and_cofactor_links():
    from sdb.compose.molecule import structure
    from sdb.validate import check_cofactor_link, check_ring_molecule
    db = model.load()
    nad = db.molecules["nad"]
    # the reactive part is coloured on the standalone figure only
    assert style.COLORS["magenta"] in structure(nad, "ring", "default").body
    assert style.COLORS["magenta"] not in structure(nad, "ring").body
    bad = copy.deepcopy(nad)
    bad["decorate"]["highlight"] = [{"smarts": "[Fe]"}]
    assert check_ring_molecule(bad)
    # cofactor labels link to structures that really contain them
    assert check_cofactor_link(db, db.cofactors["NADH"]) == []
    assert check_cofactor_link(db, {**db.cofactors["NADH"], "molecule": "nad"})   # oxidised ≠ NADH


def test_coenzyme_stereo_templates():
    from rdkit import Chem
    db = model.load()
    ribose = Chem.MolFromSmarts("[#7][C@@H]1O[C@H](C[O,S])[C@@H](O)[C@H]1O")     # β-D-ribose as in AMP
    for mid, n in [("nad", 2), ("nadph", 2), ("fad", 1), ("coa", 1), ("sam", 1)]:
        mol = Chem.MolFromSmiles(db.molecules[mid]["smiles"])
        assert len(mol.GetSubstructMatches(ribose, useChirality=True)) == n, mid
    from rdkit.Chem import rdCIPLabeler
    thf = Chem.MolFromSmiles(db.molecules["thf"]["smiles"])
    rdCIPLabeler.AssignCIPLabels(thf)
    assert sorted(a.GetProp("_CIPCode") for a in thf.GetAtoms() if a.HasProp("_CIPCode")) == ["S", "S"]   # (6S), L-Glu


def test_phosphate_bonds_stay_visible_and_overlaps_warn():
    from sdb.draw.ring import _stretch_label_bonds, label_overlaps, ring_geometry, GAP_H, GAP_V, MIN_VISIBLE
    from sdb.text import TextLine
    db = model.load()
    smi = db.molecules["coa"]["smiles"]
    atoms, bonds, pos = ring_geometry(smi, 0.82 * 1.75 * style.S)
    pos = _stretch_label_bonds(atoms, bonds, pos)
    sym = {a[0]: a[1] for a in atoms}
    half = {}
    for i, s in sym.items():
        t = TextLine(s)
        x0, x1 = t.char_box(0)
        half[i] = ((x1 - x0) / 2 + GAP_H, t.cap / 2 + GAP_V)
    for i, j, _o, centre in bonds:
        if centre is None and "P" in (sym[i], sym[j]):
            (x1, y1), (x2, y2) = pos[i], pos[j]
            L = ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5
            ux, uy = (x2 - x1) / L, (y2 - y1) / L
            cut = sum(min(hw / abs(ux) if ux else 1e9, hh / abs(uy) if uy else 1e9) for hw, hh in (half[i], half[j]))
            assert L - cut >= MIN_VISIBLE - 1e-6
    _, warnings, _ = validate(db)
    assert not [w for w in warnings if "重疊" in w]
    assert label_overlaps(db.molecules["nadp"]["smiles"], {"coordgen": False})      # the unfixed layout is caught


def test_round8_carb_stereo_checks():
    from sdb.validate import check_chain, check_haworth
    db = model.load()
    g3p = copy.deepcopy(db.molecules["glycerol_3p"])
    row = g3p["chain"][1]
    row["l"], row["r"] = row["r"], row["l"]
    assert check_chain(g3p)                     # symmetric CH₂ ends: the Ⓟ fixes the reading direction
    for mid in ["a_glucose_1p", "udp_galactose", "fructose_26bp", "phosphogluconolactone", "glycogen_n1"]:
        assert check_haworth(db.molecules[mid]) == [], mid
    for rid in ["glycogen_synthase", "glycogen_phosphorylase", "galt", "pyruvate_carboxylase", "fad_synthetase"]:
        assert check_balance(db, db.reactions[rid]) == [], rid


def test_round9_generic_chains_and_overall_equations():
    from sdb.validate import check_chain, check_lipid
    db = model.load()
    oh = copy.deepcopy(db.molecules["hydroxyacyl_coa"])
    assert check_chain(oh) == []
    row = oh["chain"][1]
    row["l"], row["r"] = row["r"], row["l"]
    assert check_chain(oh)                                   # D/L next to an R row is still checked
    for mid in ["palmitoyl_coa", "palmitoylcarnitine"]:
        assert check_lipid(db.molecules[mid]) == [], mid
    for rid in ["fatty_acid_synthase", "beta_ketothiolase", "hmg_coa_reductase", "scot", "methylmalonyl_coa_mutase"]:
        assert check_balance(db, db.reactions[rid]) == [], rid


def test_round10_heme_catecholamines_and_one_carbon():
    from rdkit import Chem
    from sdb.draw.ring import label_overlaps
    from sdb.validate import check_chain
    db = model.load()
    heme = Chem.MolFromSmiles(db.molecules["heme"]["smiles"])
    fe = next(a for a in heme.GetAtoms() if a.GetSymbol() == "Fe")
    assert fe.GetDegree() == 4 and Chem.GetFormalCharge(heme) == -2          # Fe bonded to the four N, 2 propionates
    assert not label_overlaps(db.molecules["heme"]["smiles"], db.molecules["heme"]["ring"])
    ne = copy.deepcopy(db.molecules["norepinephrine"])
    row = ne["chain"][1]
    row["l"], row["r"] = row["r"], row["l"]
    assert check_chain(ne)                                       # (R)-norepinephrine is enforced
    for rid in ["nitric_oxide_synthase", "pnmt", "methionine_synthase", "ala_synthase", "ala_dehydratase", "cystathionase"]:
        assert check_balance(db, db.reactions[rid]) == [], rid


def test_colours_are_paired_and_point_at_something():
    """Every decarboxylation colours its CO₂; every colour setting hits a drawn atom or text."""
    from sdb.validate import check_decorate
    db = model.load()
    for r in db.reactions.values():
        assert check_decorate(db, r) == [], r["id"]
        deco = r.get("decorate") or {}
        co2 = [p["id"] for p in r["products"] if p["id"] == "co2"] + [c for c in r.get("cofactors_out") or [] if "CO2" in c]
        if co2 and any(k in deco for k in [p["id"] for p in r["substrates"]]):
            assert any(c in deco for c in co2), f"{r['id']}: CO₂ 沒有上色"
    hppd = copy.deepcopy(db.reactions["hppd"])
    hppd["decorate"]["homogentisate"]["colors"][2]["ring"][0]["atoms"] = [42]
    hppd["decorate"]["co3"] = {"colors": {0: "magenta"}}
    hppd["decorate"]["O2"] = {"colors": {0: {"c": ["magenta", "N"]}}}
    assert len(check_decorate(db, hppd)) == 3
    svg = compose_reaction(db, db.reactions["tyrosine_hydroxylase"]).body
    assert svg.count(style.color("magenta")) >= 3          # new ring OH + O₂ + H₂O


def test_double_bond_cis_trans_is_checked():
    db = model.load()
    assert check_chain(db.molecules["cis_aconitate"]) == []
    assert check_chain(db.molecules["fumarate"]) == []
    wrong = copy.deepcopy(db.molecules["cis_aconitate"])
    wrong["chain"][3] = {"c": "C", "r": "H"}                # H on the same side as C2's COO⁻ = trans
    assert any("雙鍵" in p for p in check_chain(wrong))


def test_round11_etc_and_ethanol():
    from sdb.compose.pathway import compose_pathway
    db = model.load()
    for rid in ["complex_i", "complex_ii", "complex_iii", "complex_iv", "atp_synthase", "mito_glycerol_3p_dehydrogenase",
                "alcohol_dehydrogenase", "aldehyde_dehydrogenase", "meos_cyp2e1", "acetyl_coa_synthetase",
                "methanol_dehydrogenase", "formaldehyde_dehydrogenase"]:
        assert check_balance(db, db.reactions[rid]) == [], rid
    g = compose_pathway(db, db.pathways["electron_transport_chain"])
    assert style.color("membrane") in g.body and style.color("protein") in g.body      # bilayer + complexes
    assert style.color("orange") in g.body                                             # electron-flow arrows
    assert db.cofactors["Q"].get("molecule") == "ubiquinone"
    assert db.molecules["acetaldehyde"]["default_view"] == "chain"
    aldol = db.reactions["aldol_condensation"]
    assert aldol["substrates"][0].get("view") == "skeletal"                            # organic round keeps skeletal


def test_round12_lipids_and_finer_parts():
    from sdb.compose.pathway import compose_pathway
    db = model.load()
    for rid in ["cyp11a1", "aromatase", "aldosterone_synthase", "cyp17_lyase", "ltc4_synthase", "cox_peroxidase",
                "serine_palmitoyltransferase", "hexosaminidase_a", "arylsulfatase_a", "sphingomyelinase"]:
        assert check_balance(db, db.reactions[rid]) == [], rid
    from sdb.draw.ring import ring_geometry
    atoms, _, _ = ring_geometry(db.molecules["ltc4"]["smiles"], 20)
    assert "G" in [a[1] for a in atoms]                                                   # [16*] drawn as "G"
    r = compose_reaction(db, db.reactions["complex_i"])
    assert {"arrow_line", "cofactors", "prosthetic"} <= set(r.parts)                       # finer reaction layers
    p = compose_pathway(db, db.pathways["electron_transport_chain"])
    assert {"regions", "proteins", "electrons", "protons", "inhibitions", "drugs"} <= set(p.parts)
    assert "atp_synthase_fof1" in db.pathways


def test_round13_plots_are_checked():
    """Plot points on curves are verified (K_m at V_max/2, LB intercepts, P50, pK at half-equivalence)."""
    from copy import deepcopy
    from sdb.compose.pathway import compose_pathway
    from sdb.draw.plot import check_plot
    db = model.load()
    for pid in ["michaelis_menten", "lineweaver_burk", "inhibition_competitive", "inhibition_noncompetitive",
                "inhibition_uncompetitive", "o2_binding_mb_hb", "hill_plot", "titration_glycine", "ramachandran_plot"]:
        p = db.pathways[pid]
        assert p["layout"] == "plot"
        for spec in p.get("panels") or [p["plot"]]:
            assert check_plot(spec) == [], pid
        g = compose_pathway(db, p)
        assert {"axes", "curves"} <= set(g.parts), pid
    assert {"panel_1", "panel_2"} <= set(compose_pathway(db, db.pathways["inhibition_competitive"]).parts)
    # a wrong mark is caught: K_m does not sit at 0.4·V_max, and −1/K_m is not the LB x intercept at −1
    mm = deepcopy(db.pathways["michaelis_menten"]["plot"])
    mm["points"] = [{"on": "mm", "x": "Km", "y": "0.4*Vmax"}]
    assert check_plot(mm)
    lb = deepcopy(db.pathways["lineweaver_burk"]["plot"])
    lb["points"] = [{"on": "lb", "x": -1, "y": 0}]
    assert check_plot(lb)
    # glycine titration: pI = (pK1 + pK2)/2 at one equivalent
    t = db.pathways["titration_glycine"]["plot"]
    assert abs(t["params"]["pK1"] + t["params"]["pK2"] - 2 * 5.97) < 0.01


def test_round13_peptides():
    db = model.load()
    assert check_balance(db, db.reactions["dipeptide_formation"]) == []
    errors, _, _ = validate(db)
    assert not [e for e in errors if "glycylalanine" in e or "glutathione" in e]
    edges = db.pathways["peptide_bond_resonance"]["edges"]
    assert edges[0]["resonance"]


def test_round14_polysaccharides_and_glyoxylate():
    """GAG repeating units are checked against SMILES assembled with RDKit (stereo kept), and a wrong
    Haworth (L-IdoA drawn with COO⁻ up, i.e. D-GlcA) is caught."""
    from sdb.validate import check_haworth
    db = model.load()
    for mid in ["hyaluronate_unit", "chondroitin_4_sulfate_unit", "dermatan_sulfate_unit", "heparin_unit",
                "keratan_sulfate_unit", "chitin_unit", "peptidoglycan_unit"]:
        assert check_haworth(db.molecules[mid]) == [], mid
    ds = copy.deepcopy(db.molecules["dermatan_sulfate_unit"])
    ido = ds["haworth"]["units"][1]["subs"][5]
    ido["up"], ido["down"] = ido["down"], ido["up"]                  # COO⁻ to the other face → D-GlcA
    assert check_haworth(ds)
    for rid in ["isocitrate_lyase", "malate_synthase", "rubisco_carboxylase", "rubisco_oxygenase"]:
        assert check_balance(db, db.reactions[rid]) == [], rid
    from sdb.compose.pathway import compose_pathway
    for pid in ["glyoxylate_cycle", "calvin_cycle", "photorespiration", "glycosaminoglycans",
                "polysaccharide_classes", "mucopolysaccharidoses"]:
        compose_pathway(db, db.pathways[pid])


def test_round15_notes_gaps():
    """Items added from the lecture-note comparison: stereo-checked sugars/PTMs, balanced reactions,
    plots whose marked points sit on their curves, and dashed H-bonds in base pairs."""
    from sdb.compose.molecule import structure
    from sdb.compose.pathway import compose_pathway
    from sdb.draw.plot import check_plot
    from sdb.validate import check_haworth
    db = model.load()
    assert check_haworth(db.molecules["trehalose"]) == []
    for mid in ["d_gluconate", "d_glucuronate_open", "enediol", "bpg_23", "hydroxylysine", "gla", "selenocysteine"]:
        assert check_chain(db.molecules[mid]) == [], mid
    for rid in ["bpg_mutase", "bpg_phosphatase", "disulfide_reduction", "prolyl_hydroxylase", "lysyl_hydroxylase",
                "gamma_glutamyl_carboxylase", "protein_kinase_ser", "histone_acetyltransferase",
                "histone_methyltransferase", "cytosine_deamination", "methylcytosine_deamination"]:
        assert check_balance(db, db.reactions[rid]) == [], rid
    for pid in ["inhibition_mixed", "substrate_inhibition", "burst_kinetics", "reaction_coordinate",
                "lysozyme_ph_profile", "hb_t_r_states", "titration_his_glu", "dna_melting_curve"]:
        p = db.pathways[pid]
        for spec in p.get("panels") or [p["plot"]]:
            assert check_plot(spec) == [], pid
        compose_pathway(db, p)
    g = structure(db.molecules["base_pair_gc"], "free", None)
    assert g.body.count('stroke-dasharray="3 3"') == 3                                   # G≡C: three H-bonds
    assert structure(db.molecules["base_pair_at"], "free", None).body.count('stroke-dasharray="3 3"') == 2


def test_round16_metabolism_notes():
    """Items added from the metabolism lecture notes: balanced reactions, a stereo-checked UDP-GlcA,
    plots that pass check_plot, and carbon-fate marks that render on every glucose carbon."""
    from sdb.compose.pathway import compose_pathway
    from sdb.draw.plot import check_plot
    from sdb.validate import check_haworth
    db = model.load()
    assert check_haworth(db.molecules["udp_glucuronate"]) == []
    for rid in ["pyruvate_decarboxylase", "lactose_synthase", "udp_glucose_dehydrogenase"]:
        assert check_balance(db, db.reactions[rid]) == [], rid
    for pid in ["energy_charge", "hexokinase_glucokinase", "glycolysis_free_energy", "phosphoryl_transfer_potential"]:
        assert check_plot(db.pathways[pid]["plot"]) == [], pid
    for pid in ["glucose_carbon_fate", "fermentation", "ppp_nonoxidative", "hexose_interconversion",
                "tca_beta_oxidation_analogy", "glycogen_enzyme_phosphorylation", "xenobiotic_detoxification",
                "futile_cycle", "atp_cleavage_modes"]:
        compose_pathway(db, db.pathways[pid])


def test_round17_signal_transduction():
    """Book 3 items: cyclic nucleotides keep the sugar stereo of AMP/GMP, the cyclase/PDE reactions
    balance, and every signalling / membrane figure composes."""
    from rdkit import Chem
    from sdb.compose.pathway import compose_pathway
    from sdb.draw.plot import check_plot
    from sdb.validate import check_haworth
    db = model.load()
    assert check_haworth(db.molecules["gtp"]) == []
    for cyc, lin in (("camp", "amp"), ("cgmp", "gmp")):
        a = Chem.FindMolChiralCenters(Chem.MolFromSmiles(db.molecules[cyc]["smiles"]))
        b = Chem.FindMolChiralCenters(Chem.MolFromSmiles(db.molecules[lin]["smiles"]))
        assert len(a) == len(b) == 4, cyc
    for rid in ["adenylyl_cyclase", "camp_phosphodiesterase", "guanylyl_cyclase", "cgmp_phosphodiesterase"]:
        assert check_balance(db, db.reactions[rid]) == [], rid
    assert check_plot(db.pathways["membrane_fluidity"]["plot"]) == []
    for pid in ["receptor_types", "gpcr_g_protein_classes", "camp_pka_cascade", "no_cgmp_relaxation",
                "insulin_rtk_signaling", "jak_stat_pathway", "nfkb_pathway", "apoptosis_pathways", "gtpase_switch",
                "p53_rb_checkpoint", "two_component_system", "bacterial_toxins", "second_messengers",
                "membrane_transporters", "membrane_asymmetry", "ldl_receptor_endocytosis"]:
        compose_pathway(db, db.pathways[pid])


def test_round18_scene_engine():
    """Illustrated scenes: every glyph renders on its own, the pilot scenes compose with one part per
    item, validation flags unknown glyphs, and `parts_png: false` / scene variants export SVG-only parts."""
    from sdb.compose.pathway import compose_pathway
    from sdb.draw.glyphs import GLYPHS
    from sdb.draw.scene import check_scene
    samples = {"membrane": {"y": 0, "x0": 0, "x1": 5}, "region": {"x0": 0, "y0": 0, "x1": 3, "y1": 2},
               "strand": {"points": [[0, 0], [4, 0]], "bases": "left", "arrow": True, "ends": ["5'", "3'"]},
               "dsdna": {"from": [0, 0], "to": [4, 0], "ends": True}, "arrow": {"from": [0, 0], "to": [3, 1]},
               "text": {"at": [0, 0], "text": "x"}, "helix": {"from": [0, 0], "to": [6, 0]}}
    for g, fn in GLYPHS.items():
        svg, bb, _ = fn({"at": [0, 0], **samples.get(g, {})})
        assert svg and bb[2] > bb[0] and bb[3] > bb[1], g
    db = model.load()
    p = db.pathways["replication_fork"]
    assert p.get("parts_png") is False
    g = compose_pathway(db, p)
    assert len(g.parts) == len(p["scene"]["items"])
    for pid in ["membrane_transporters", "camp_pka_cascade"]:
        assert compose_pathway(db, db.pathways[pid], "illus").parts
    assert check_scene({"items": [{"g": "nope"}]})[0]


def test_round18_dna_replication():
    """Round 18 content: chain terminators keep the deoxyribose stereo (Haworth check), new glyphs render,
    every replication scene composes without label overlaps, all new items are SVG-only stickers."""
    from sdb.compose.pathway import compose_pathway
    from sdb.draw.scene import check_scene
    from sdb.validate import check_haworth
    db = model.load()
    for mid in ["datp", "ddatp", "azt"]:
        assert check_haworth(db.molecules[mid]) == [], mid
    for pid in ["dna_double_helix", "chromatin_levels", "oric_initiation", "bacterial_chromosome_replication",
                "topoisomerase_types", "telomere_replication", "eukaryotic_origin_licensing", "dna_polymerase_mechanism",
                "dna_polymerases", "replication_fork"]:
        p = db.pathways[pid]
        assert p.get("parts_png") is False, pid
        assert check_scene(p["scene"]) == ([], []), pid
        assert compose_pathway(db, p).parts, pid
    compose_pathway(db, db.pathways["chain_terminators"])
    compose_pathway(db, db.pathways["dna_ligase_mechanism"])


def test_round19_repair_recombination():
    """Round 19: MGMT reaction balances, damage structures parse, every repair / recombination scene composes
    with no label overlaps and SVG-only stickers."""
    from sdb.compose.pathway import compose_pathway
    from sdb.draw.scene import check_scene
    db = model.load()
    assert check_balance(db, db.reactions["mgmt"]) == []
    for pid in ["mismatch_repair", "base_excision_repair", "nucleotide_excision_repair", "homologous_recombination", "nhej",
                "sos_response", "site_specific_recombination", "transposition", "vdj_recombination", "dna_repair_pathways",
                "mutation_types"]:
        p = db.pathways[pid]
        assert p.get("parts_png") is False, pid
        assert check_scene(p["scene"]) == ([], []), pid
        assert compose_pathway(db, p).parts, pid
    compose_pathway(db, db.pathways["dna_damage_structures"])


def test_round20_transcription():
    """Round 20: editing reactions balance, cap / Ψ structures parse, every transcription / RNA-processing scene
    composes with no label overlaps and SVG-only stickers; the tRNA glyph draws with all options."""
    from sdb.compose.pathway import compose_pathway
    from sdb.draw.glyphs import GLYPHS
    from sdb.draw.scene import check_scene
    db = model.load()
    for rid in ["apobec1_editing", "adar_editing"]:
        assert check_balance(db, db.reactions[rid]) == [], rid
    for pid in ["mrna_prok_vs_euk", "bacterial_transcription_initiation", "transcription_elongation", "transcription_termination",
                "pol2_preinitiation_complex", "rna_polymerases", "mrna_capping", "polyadenylation", "spliceosome_splicing",
                "splicing_mechanisms", "alternative_splicing", "rna_editing", "trna_processing", "rrna_processing"]:
        p = db.pathways[pid]
        assert p.get("parts_png") is False, pid
        assert check_scene(p["scene"]) == ([], []), pid
        assert compose_pathway(db, p).parts, pid
    for pid in ["mrna_cap_structure", "rna_modified_nucleosides"]:
        compose_pathway(db, db.pathways[pid])
    svg, bb, _ = GLYPHS["trna"]({"g": "trna", "at": [0, 0], "leader": 2, "trailer": 2, "cca": False, "intron": True})
    assert svg and bb[3] > bb[1]


def test_round21_translation():
    """Round 21: aminoacylation / peptidyl-transfer reactions balance, every translation / protein-fate scene composes
    with no label overlaps and SVG-only stickers; the ribosome glyph draws both ways up and the tRNA icon draws."""
    from sdb.compose.pathway import compose_pathway
    from sdb.draw.glyphs import GLYPHS
    from sdb.draw.scene import check_scene
    db = model.load()
    for rid in ["aars_activation", "aars_transfer", "peptidyl_transfer"]:
        assert check_balance(db, db.reactions[rid]) == [], rid
    for pid in ["genetic_code", "wobble_pairing", "aminoacyl_trna_synthetase", "ribosome_structure", "translation_initiation",
                "translation_elongation", "translation_termination", "translation_inhibitors", "srp_er_targeting",
                "ubiquitin_proteasome", "translation_factors", "translation_energy_cost", "protein_sorting_signals",
                "post_translational_modifications"]:
        p = db.pathways[pid]
        assert p.get("parts_png") is False, pid
        assert check_scene(p["scene"]) == ([], []), pid
        assert compose_pathway(db, p).parts, pid
    compose_pathway(db, db.pathways["puromycin_vs_trna"])
    up = GLYPHS["ribosome"]({"g": "ribosome", "at": [0, 0], "sites": True, "large": "50S"})[1]
    down = GLYPHS["ribosome"]({"g": "ribosome", "at": [0, 0], "flip": True})[1]
    assert up[1] < 0 < up[3] and -up[1] > up[3] and down[3] > -down[1]       # large subunit above / below the mRNA
    assert GLYPHS["trna"]({"g": "trna", "shape": "icon", "at": [0, 0], "anticodon": "UAC", "peptide": 3})[0]


def test_round22_gene_regulation():
    """Round 22: DNMT reaction balances, every gene-regulation scene composes with no label overlaps, SVG-only stickers
    and no long sentences on the picture (explanations belong in notes_zh)."""
    from sdb.compose.pathway import compose_pathway
    from sdb.draw.scene import check_scene
    db = model.load()
    assert check_balance(db, db.reactions["dna_methyltransferase"]) == []
    for pid in ["control_modes", "lac_operon", "trp_attenuation", "stringent_response", "bacterial_riboswitch_srna", "phase_variation",
                "eukaryotic_activation", "dna_binding_motifs", "chromatin_epigenetics", "rnai_mirna_sirna", "nuclear_receptor_signaling"]:
        p = db.pathways[pid]
        assert p.get("parts_png") is False, pid
        assert check_scene(p["scene"]) == ([], []), pid
        assert compose_pathway(db, p).parts, pid
        for it in p["scene"]["items"]:
            t = it.get("text") if it["g"] == "text" else it.get("label")
            assert not t or len(t) <= 45, (pid, t)
    for pid in ["small_regulatory_rnas", "lac_inducers", "epigenetic_marks"]:
        compose_pathway(db, db.pathways[pid])


def test_round23_techniques():
    """Round 23: X-gal hydrolysis balances; every technique scene composes with no label overlaps, SVG-only stickers and
    short on-figure text; the qPCR plot's Ct points sit on their curves."""
    from sdb.compose.pathway import compose_pathway
    from sdb.draw.plot import check_plot
    from sdb.draw.scene import check_scene
    db = model.load()
    assert check_balance(db, db.reactions["xgal_hydrolysis"]) == []
    for pid in ["pcr_cycle", "sanger_sequencing", "blotting_methods", "restriction_cloning", "blue_white_selection", "cdna_synthesis",
                "crispr_cas9", "dna_protein_assays", "yeast_two_hybrid", "rflp_sickle_cell"]:
        p = db.pathways[pid]
        assert p.get("parts_png") is False, pid
        assert check_scene(p["scene"]) == ([], []), pid
        assert compose_pathway(db, p).parts, pid
        for it in p["scene"]["items"]:
            t = it.get("text") if it["g"] == "text" else it.get("label")
            assert not t or len(t) <= 45, (pid, t)
    assert check_plot(db.pathways["qpcr_ct"]["plot"]) == []
    for pid in ["qpcr_ct", "cloning_vectors", "techniques_by_purpose", "xgal_blue_dye"]:
        compose_pathway(db, db.pathways[pid])


def test_round25_signalling_illustrations():
    """Round 25: the membrane / signalling / energy figures gained an `illus` variant (boxed version untouched) and two new
    membrane scenes; every scene composes with no label overlaps and short on-figure text."""
    from sdb.compose.pathway import compose_pathway
    from sdb.draw.scene import check_scene
    db = model.load()
    illus = ["receptor_types", "gpcr_g_protein_classes", "g_protein_cycle", "gtpase_switch", "phosphoinositide_signaling",
             "insulin_rtk_signaling", "jak_stat_pathway", "nfkb_pathway", "no_cgmp_relaxation", "two_component_system",
             "bacterial_toxins", "apoptosis_pathways", "ldl_receptor_endocytosis", "membrane_asymmetry",
             "electron_transport_chain", "atp_synthase_fof1", "malate_aspartate_shuttle", "glycerol_phosphate_shuttle"]
    scenes = [(pid, db.pathways[pid]["variants"]["illus"]["scene"]) for pid in illus]
    scenes += [(pid, db.pathways[pid]["scene"]) for pid in ("membrane_protein_classes", "membrane_permeability")]
    for pid, sc in scenes:
        assert check_scene(sc) == ([], []), pid
        for it in sc["items"]:
            t = it.get("text") if it["g"] == "text" else it.get("label")
            assert not t or len(t) <= 45, (pid, t)
    for pid in illus:
        p = db.pathways[pid]
        assert p.get("layout") != "scene", pid               # the boxed version is still the default
        assert compose_pathway(db, p, "illus").parts, pid
        compose_pathway(db, p)
    for pid in ("membrane_protein_classes", "membrane_permeability"):
        assert db.pathways[pid].get("parts_png") is False
        assert compose_pathway(db, db.pathways[pid]).parts


def test_round26_biochem_gaps():
    """Round 26: serine-protease reactions balance; the new scenes compose with no label overlaps and short text; the
    starvation plot and the THF / phospholipase figures compose."""
    from sdb.compose.pathway import compose_pathway
    from sdb.compose.reaction import compose_reaction
    from sdb.draw.plot import check_plot
    from sdb.draw.scene import check_scene
    db = model.load()
    for rid in ("chymotrypsin_acylation", "chymotrypsin_deacylation"):
        assert check_balance(db, db.reactions[rid]) == [], rid
        compose_reaction(db, db.reactions[rid])
    for pid in ("serine_protease_mechanism", "chaperonin_groel_cycle", "pdh_complex_mechanism", "rbc_glutathione_cycle"):
        p = db.pathways[pid]
        assert p.get("parts_png") is False
        assert check_scene(p["scene"]) == ([], []), pid
        for it in p["scene"]["items"]:
            t = it.get("text") if it["g"] == "text" else it.get("label")
            assert not t or len(t) <= 45, (pid, t)
        assert compose_pathway(db, p).parts
    assert check_plot(db.pathways["starvation_fuel_plot"]["plot"]) == []
    for pid in ("starvation_fuel_plot", "thf_one_carbon_forms", "phospholipase_sites", "protein_purification",
                "protein_sequencing_strategy"):
        compose_pathway(db, db.pathways[pid])
    for n in db.pathways["phospholipase_sites"]["nodes"]:
        assert not n.get("text") or len(n["text"]) <= 45, n["text"]
