"""Atom-tracking colours (reaction `decorate`) for the reactions written by the generators and
for the hand-written nucleotide reactions — one place to review them all.

    .venv/bin/python tools/reaction_colors.py data      # writes `decorate` into the YAML files

The generators (gen_vit / gen_carb / gen_lipmet / gen_aamet) also read COLORS, so re-running
them keeps these colours. Rules (Mathews / Lehninger, see AUTHORING.md):
  - decarboxylation: the leaving COO⁻ and the CO₂ share a colour
  - dehydrogenase (oxidation): the H atoms that leave and the H / H⁺ of NAD(P)H + H⁺ share a colour;
    reduction: the H atoms that arrive and NAD(P)H + H⁺ share a colour
  - group transfer / hydrolysis / carboxylation: source and destination share a colour
    (Ⓟ ↔ Pi, PPi; HCO₃⁻ ↔ new COO⁻; H₂O ↔ new OH; Gln ↔ new NH₂; SAM ↔ CH₃ …)
  - a second tracked fragment is cyan
"""
import sys
from pathlib import Path

import yaml

M, C = "magenta", "cyan"


def whole(col=M):
    """whole-label colour for a cofactor / formula molecule."""
    return {"colors": {0: col}}


def frag(*pairs):
    """colour fragments of a cofactor / formula label: frag((M, "H"), (M, "H^+"))."""
    return {"colors": {0: {"c": [list(p) for p in pairs]}}}


def rows(**kw):
    return {"colors": kw}


def hl(*smarts, col=M):
    """ring-view molecule: highlight every SMARTS match (mapped atoms only when the SMARTS has maps)."""
    return {"highlight": [{"smarts": sm, "color": col} for sm in smarts]}


def ring_atoms(*specs):
    """ring-view molecule: highlight([7, 8], M), …"""
    return {"highlight": [{"atoms": list(a), "color": col} for a, col in specs]}


NADH_M = frag((M, "H"), (M, "H^+"))        # NADH + H⁺ / NADPH + H⁺ / 2 NADPH + 2 H⁺
NADH_C = frag((C, "H"), (C, "H^+"))


def chain(**rows_):
    """chain colours keyed by row number: chain(r0=M, r1={"l": M})."""
    return {"colors": {int(k[1:]): v for k, v in rows_.items()}}


def ring_row(row, atoms, col=M):
    """chain row that is a ring (Tyr, Trp …): colour only these SMILES atoms of the ring."""
    return {row: {"ring": [{"atoms": list(atoms), "color": col}]}}


def hw(**keys):
    """Haworth substituent colours: hw(d1=[M, "{P}"], u4=M)."""
    return {"colors": keys}


COLORS = {
    # ------------------------------------------------------------------ round 10: amino acids
    "phenylalanine_hydroxylase": {"tyrosine": {"colors": ring_row(3, [5])}, "O2": whole(), "H2O": whole()},
    "tyrosine_hydroxylase": {"l_dopa": {"colors": ring_row(3, [7])}, "O2": whole(), "H2O": whole()},
    "tryptophan_hydroxylase": {"five_htp": {"colors": ring_row(3, [8])}, "O2": whole(), "H2O": whole()},
    "dopa_decarboxylase": {"l_dopa": chain(r0=M), "co2": whole()},
    "htp_decarboxylase": {"five_htp": chain(r0=M), "co2": whole()},
    "histidine_decarboxylase": {"histidine": chain(r0=M), "co2": whole()},
    "glutamate_decarboxylase": {"glutamate": chain(r0=M), "co2": whole()},
    "dopamine_beta_hydroxylase": {"norepinephrine": chain(r1={"r": M}), "O2": whole(), "H2O": whole()},
    "pnmt": {"SAM": whole(), "epinephrine": chain(r0={"l": [M, "CH_3"]})},
    "tyrosine_aminotransferase": {"tyrosine": chain(r1={"l": M}), "glutamate": chain(r1={"l": M})},
    # HPPD: the α-COO⁻ leaves as CO₂ (magenta); one O of O₂ hydroxylates the ring (cyan)
    "hppd": {"hydroxyphenylpyruvate": chain(r0=M), "co2": whole(), "O2": whole(C),
             "homogentisate": {"colors": ring_row(2, [8], C)}},
    # NOS: a guanidino N becomes NO (magenta); O from O₂ goes to NO and to the citrulline C=O (cyan)
    "nitric_oxide_synthase": {"arginine": chain(r6={"r": M}), "nitric_oxide": frag((M, "N"), (C, "O")),
                              "citrulline": chain(r6={"r": C}), "4O2": whole(C)},
    # MAT: the S of Met attacks C5′ of ATP (magenta); all three phosphates leave (cyan)
    "methionine_adenosyltransferase": {"methionine": chain(r4=M), "sam": ring_atoms(([13, 14], M)),
                                       "atp_n": hw(u4=[C, "{P}–{P}–{P}"]), "PPi": whole(C), "Pi": whole(C)},
    "sah_hydrolase": {"sah": ring_atoms(([13, 14], M)), "homocysteine": chain(r4=M),
                      "H2O": whole(C), "adenosine": hw(u4=[C, "HO"])},
    "methionine_synthase": {"methyl_thf": frag((M, "CH_3")), "methionine": chain(r5=M)},
    "mthfr": {"methylene_thf": frag((M, "CH_2")), "methyl_thf": frag((M, "CH_3"))},
    "cystathionine_beta_synthase": {"homocysteine": chain(r4=M), "cystathionine": chain(r4=M),
                                    "serine": chain(r3=C), "H2O": whole(C)},
    # cystathionase: the Hcy carbon skeleton (rows 0–3) → α-ketobutyrate; its NH₃⁺ → NH₄⁺
    "cystathionase": {"cystathionine": chain(r0=M, r1={"c": M, "l": C}, r2=M, r3=M),
                      "alpha_ketobutyrate": chain(r0=M, r1={"c": M}, r2=M, r3=M), "ammonium": whole(C)},
    "serine_hydroxymethyltransferase": {"serine": chain(r2=M, r3=C), "methylene_thf": frag((M, "CH_2")),
                                        "H2O": whole(C)},
    # ALAS: Gly's Cα–N stays (magenta), its COO⁻ leaves as CO₂ (cyan)
    "ala_synthase": {"glycine": chain(r0=C, r1=M), "aminolevulinate": chain(r4=M), "co2": whole(C)},
    "ala_dehydratase": {"aminolevulinate": chain(r3={"r": M}), "2H2O": whole()},

    # ------------------------------------------------------------------ round 9: lipids
    "acyl_coa_synthetase": {"palmitate": {"colors": {"o2": M}}, "AMP": whole(),
                            "CoA-SH": whole(C), "palmitoyl_coa": {"colors": {"o2": C, "rest": C}}},
    "cpt1": {"palmitoyl_coa": {"colors": {"o2": M, "rest": M}}, "CoA-SH": whole(),
             "carnitine": chain(r2={"l": C}), "palmitoylcarnitine": {"colors": {"o2": C, "rest": C}}},
    "acyl_coa_dehydrogenase": {"acyl_coa": chain(r1=M, r2=M), "enoyl_coa": chain(r1=M, r2=M),
                               "FADH2": frag((M, "H_2"))},
    "enoyl_coa_hydratase": {"H2O": whole(), "hydroxyacyl_coa": chain(r1={"r": M})},
    "hydroxyacyl_coa_dehydrogenase": {"hydroxyacyl_coa": chain(r1={"l": M, "r": [M, "H"]}), "NADH+H+": NADH_M},
    # thiolase: CoA-SH (cyan) attacks C3; C1–C2 with the old CoA (magenta) leave as acetyl-CoA
    "beta_ketothiolase": {"ketoacyl_coa": chain(r2=M, r3=M), "acetyl_coa": chain(r0=M),
                          "CoA-SH": whole(C), "acyl_coa_n2": chain(r1={"r": C})},
    "propionyl_coa_carboxylase": {"bicarbonate": whole(), "s_methylmalonyl_coa": chain(r0=M)},
    "methylmalonyl_coa_epimerase": {"s_methylmalonyl_coa": chain(r1={"l": M, "r": C}),
                                    "r_methylmalonyl_coa": chain(r1={"l": C, "r": M})},
    # ATP-citrate lyase = citrate synthase backwards (the acetyl part of citrate is rows 0–1)
    "citrate_lyase": {"citrate": chain(r0=M, r1=M), "CoA-SH": whole(C),
                      "acetyl_coa": {"colors": {0: {"c": M, "l": M, "u": M, "r": C}}}},
    # malic enzyme: like IDH — the H on C2 → NADPH (cyan), the C4 COO⁻ → CO₂ (magenta)
    "malic_enzyme": {"malate": chain(r0=M, r1={"r": C}), "co2": whole(), "NADPH": frag((C, "H"))},
    "acetyl_coa_carboxylase": {"bicarbonate": whole(), "malonyl_coa": chain(r0=M)},
    "fatty_acid_synthase": {"malonyl_coa": chain(r0=M), "7CO2": whole()},
    # the first acetyl-CoA acylates the enzyme (its CoA leaves); the second one keeps its CoA
    "acetoacetyl_coa_thiolase": {"acetoacetyl_coa": chain(r0=M, r1=M), "CoA-SH": whole()},
    "hmg_coa_synthase": {"acetyl_coa": {"colors": {0: {"c": M, "l": M, "u": M, "r": C}}}, "CoA-SH": whole(C),
                         "hmg_coa": chain(r0=M, r1=M)},
    "hmg_coa_lyase": {"hmg_coa": chain(r3=M, r4=M), "acetyl_coa": chain(r0=M)},
    "hydroxybutyrate_dehydrogenase": {"NADH+H+": NADH_M, "beta_hydroxybutyrate": chain(r2={"l": M, "r": [M, "H"]})},
    "acetoacetate_decarboxylation": {"acetoacetate": chain(r0=M), "co2": whole()},
    "scot": {"succinyl_coa": chain(r3={"r": M}), "acetoacetyl_coa": chain(r3={"r": M})},
    # HMG-CoA reductase: the thioester C becomes CH₂OH with 4 H from 2 NADPH; CoA leaves (cyan)
    "hmg_coa_reductase": {"hmg_coa": chain(r4={"c": M, "l": M, "r": C}), "2NADPH+H+": NADH_M,
                          "mevalonate": chain(r4=M), "CoA-SH": whole(C)},
    "mevalonate_to_ipp": {"mevalonate": chain(r0=M, r2={"r": C}), "co2": whole(), "Pi": whole(C)},

    # ------------------------------------------------------------------ round 8: carbohydrates
    "pyruvate_carboxylase": {"bicarbonate": whole(), "oxaloacetate": chain(r3=M)},
    "pepck": {"oxaloacetate": chain(r3=M), "co2": whole(), "GTP": whole(C), "pep": chain(r1={"r": [C, "{P}"]})},
    "fbpase1": {"fructose_16bp": chain(r0={"c": [M, "{P}"]}), "Pi": whole()},
    "g6pase": {"glucose_6p": chain(r5={"c": [M, "{P}"]}), "Pi": whole()},
    "fbpase2": {"fructose_26bp": hw(u1=[M, "{P}"]), "Pi": whole()},
    "ldh": {"NADH+H+": NADH_M, "lactate": chain(r1={"l": [M, "H"], "r": M})},
    "glycerol_kinase": {"glycerol_3p": chain(r2={"c": [M, "{P}"]})},
    "glycerol_3p_dehydrogenase": {"glycerol_3p": chain(r1={"l": [M, "H"], "r": M}), "NADH+H+": NADH_M},
    "pfk2": {"fructose_26bp": hw(u1=[M, "{P}"])},
    "phosphoglucomutase": {"a_glucose_1p": hw(d1=[M, "{P}"]), "a_glucose_6p": hw(u5=[M, "{P}"])},
    "udp_glucose_pyrophosphorylase": {"utp": whole(), "udp_glucose": hw(d1=[M, "UDP"])},
    # glycogen synthase: the 4-OH O of the chain end forms the new α1→4 link; UDP leaves with C1's O
    "glycogen_synthase": {"udp_glucose": hw(d1=C), "udp": whole(C), "glycogen_n": hw(d4=[M, "O"]),
                          "glycogen_n1": {"colors": {"link0": M}}},
    # phosphorylase: Pi attacks C1 (magenta); the link O stays on the shortened chain (cyan)
    "glycogen_phosphorylase": {"Pi": whole(), "a_glucose_1p": hw(d1=M),
                               "glycogen_n1": {"colors": {"link0": C}}, "glycogen_n": hw(d4=[C, "O"])},
    "g6pd": {"b_glucose_6p": hw(u1=[M, "H"], d1=M), "NADPH+H+": NADH_M},
    "lactonase": {"H2O": whole(), "phosphogluconate": chain(r0={"fork_r": M})},
    # 6PGD: like IDH — H from C3 → NADPH (cyan), C1 COO⁻ → CO₂ (magenta)
    "phosphogluconate_dehydrogenase": {"phosphogluconate": chain(r0=M, r2={"l": [C, "H"], "r": C}),
                                       "co2": whole(), "NADPH": frag((C, "H"))},
    "ribose_5p_isomerase": {"ribulose_5p": chain(r0=M, r1=M), "ribose_5p": chain(r0=M, r1=M)},
    "ribulose_5p_epimerase": {"ribulose_5p": chain(r2={"l": M, "r": M}), "xylulose_5p": chain(r2={"l": M, "r": M})},
    "transaldolase": {"sedoheptulose_7p": chain(r0=M, r1=M, r2=M), "fructose_6p": chain(r0=M, r1=M, r2=M)},
    "fructokinase": {"fructose_1p": chain(r0={"c": [M, "{P}"]})},
    "aldolase_b": {"fructose_1p": chain(r0=M, r1=M, r2=M, r3=C, r4=C, r5=C), "dhap": chain(r0=M, r1=M, r2=M),
                   "d_glyceraldehyde": chain(r0=C, r1=C, r2=C)},
    "triose_kinase": {"glyceraldehyde_3p": chain(r2={"c": [M, "{P}"]})},
    "galactokinase": {"a_galactose_1p": hw(d1=[M, "{P}"])},
    # GALT: the UMP part moves from glucose 1-P to galactose 1-P
    "galt": {"udp_glucose": hw(d1=[C, "UDP"]), "udp_galactose": hw(d1=[C, "UDP"])},
    "udp_galactose_epimerase": {"udp_galactose": hw(u4=[M, "HO"]), "udp_glucose": hw(d4=[M, "HO"])},
    "aldose_reductase": {"NADPH+H+": NADH_M, "sorbitol": chain(r0=M)},
    "aldose_reductase_galactose": {"NADPH+H+": NADH_M, "galactitol": chain(r0=M)},
    "sorbitol_dehydrogenase": {"sorbitol": chain(r1={"l": M, "r": [M, "H"]}), "NADH+H+": NADH_M},

    # ------------------------------------------------------------------ round 7: vitamins
    "thiamine_pyrophosphokinase": {"tpp": {"highlight": [{"smarts": "P(=O)([O-])OP(=O)([O-])[O-]", "color": M}]}},
    "riboflavin_kinase": {"fmn": {"highlight": [{"smarts": "C[O:1][P:2](=[O:3])([O-:4])[O-:5]", "color": M}]}},
    "pyridoxal_kinase": {"plp": {"highlight": [{"smarts": "C[O:1][P:2](=[O:3])([O-:4])[O-:5]", "color": M}]}},
    "fad_synthetase": {"atp_n": hw(u4=[M, "{P}–{P}"]), "PPi": whole()},
    "folate_reduction": {"NADPH+H+": NADH_M, "dhf": ring_atoms(([4, 5], M))},
    "vitamin_d_photolysis": {"dehydrocholesterol": ring_atoms(([15, 20], M))},
    "vitamin_d_25_hydroxylase": {"O2": whole(), "calcidiol": {"highlight": [{"smarts": "[CX4](C)(C)(C)[OH:1]", "color": M}]}},
    "vitamin_d_1a_hydroxylase": {"O2": whole(), "calcitriol": {"highlight": [{"smarts": "C(=[CH2])[CH][OH:1]", "color": M}]}},
    "beta_carotene_cleavage": {"beta_carotene": ring_atoms(([19, 20], M)), "O2": whole(),
                               "retinal": {"highlight": [{"smarts": "[CH:1]=[O:2]", "color": M}]}},
    "retinol_dehydrogenase": {"retinol": {"highlight": [{"smarts": "[CH2:1][OH:2]", "color": M}]}, "NADH+H+": NADH_M},
    # the aldehyde H (not drawn) goes to NADH; the new OH comes from water
    "retinal_dehydrogenase": {"H2O": whole(), "retinoic_acid": {"highlight": [{"smarts": "C(=O)[OH:1]", "color": M}]}},

    # ------------------------------------------------------------------ round 4: nucleotides (hand-written YAML)
    "adenylosuccinate_synthetase": {"aspartate": chain(r0=C, r1={"c": C, "r": C, "l": M}, r2=C, r3=C),
                                    "imp": {"base_colors": {6: M}}, "Pi": whole(),
                                    "adenylosuccinate": {"base_colors": {6: M, **{i: C for i in range(7, 15)}}}},
    "adenylosuccinate_lyase": {"adenylosuccinate": {"base_colors": {6: M, **{i: C for i in range(7, 15)}}},
                               "amp": {"base_colors": {6: M}}, "fumarate": chain(r0=C, r1=C, r2=C, r3=C)},
    "amidophosphoribosyltransferase": {"glutamine": whole(), "pra": hw(u1=M),
                                       "prpp": hw(d1=[C, "{P}–{P}"]), "PPi": whole(C)},
    "aprt": {"prpp": hw(d1=[M, "{P}–{P}"]), "PPi": whole()},
    "hgprt": {"prpp": hw(d1=[M, "{P}–{P}"]), "PPi": whole()},
    "oprt": {"prpp": hw(d1=[M, "{P}–{P}"]), "PPi": whole()},
    "prpp_synthetase": {"ATP": whole(), "prpp": hw(d1=[M, "{P}–{P}"])},
    "cps2": {"glutamine": whole(), "bicarbonate": whole(C),
             "carbamoyl_phosphate": {"colors": {0: {"l": M, "c": C, "u": C}}}},
    "atcase": {"carbamoyl_phosphate": {"colors": {0: {"l": M, "c": M, "u": M, "r": [C, "{P}"]}}}, "Pi": whole(C),
               "carbamoyl_aspartate": ring_atoms(([0, 1, 2], M))},
    # dihydroorotase: the carbamoyl NH₂ closes the ring on the side-chain C=O (cyan); O⁻ leaves as H₂O
    "dihydroorotase": {"carbamoyl_aspartate": ring_atoms(([0, 6, 7], C), ([8], M)),
                       "dihydroorotate": ring_atoms(([5, 6, 7], C)), "H2O": whole()},
    "dhodh": {"dihydroorotate": ring_atoms(([3, 4], M)), "orotate": ring_atoms(([3, 4], M)), "QH2": frag((M, "H_2"))},
    "ctp_synthetase": {"glutamine": whole(), "ctp": {"base_colors": {5: M}}},
    "gmp_synthetase": {"glutamine": whole(), "gmp": {"base_colors": {9: M}}},
    "nucleotidase_amp": {"amp": hw(u4=[M, "{P}"]), "Pi": whole()},
    "nucleotidase_gmp": {"gmp": hw(u4=[M, "{P}"]), "Pi": whole()},
    "pnp": {"Pi": whole(), "R1P": frag((M, "{P}"))},
    "pnp_guanosine": {"Pi": whole(), "R1P": frag((M, "{P}"))},
    "dhfr": {"NADPH+H+": NADH_M, "thf": ring_atoms(([6, 27], M))},
    # ---- glycolysis (round 1): substrate-level phosphorylation — the Ⓟ ends up in ATP
    "phosphoglycerate_kinase": {"bpg_13": chain(r0={"fork_r": [M, "{P}"]}), "pg_3": chain(r0={"fork_r": M}), "ATP": whole()},
    "pyruvate_kinase": {"pep": chain(r1={"r": [M, "{P}"]}), "pyruvate": chain(r1={"r": M}), "ATP": whole()},
    # GAPDH: the aldehyde H → NADH (magenta); Pi attacks the thioester and becomes the acyl Ⓟ (cyan)
    "gapdh": {"glyceraldehyde_3p": chain(r0={"fork_r": M}), "NADH+H+": frag((M, "H")), "Pi": whole(C),
              "bpg_13": chain(r0={"fork_r": C})},
    "aldol_condensation": {"crotonaldehyde": {"smarts": "[CH]=[CH]"}, "water": whole()},
    # ------------------------------------------------------------------ round 11: energy metabolism
    "complex_i": {"NADH+H+": NADH_M, "ubiquinol": {"colors": {"o1": [M, "H"], "o4": [M, "H"]}}},
    "complex_ii": {"succinate": chain(r1={"l": M}, r2={"r": M}), "ubiquinol": {"colors": {"o1": [M, "H"], "o4": [M, "H"]}}},
    # QH₂: its two H leave as 2 H⁺ (magenta); its two electrons reduce two cyt c (cyan)
    "complex_iii": {"ubiquinol": {"colors": {"o1": [M, "H"], "o4": [M, "H"]}}, "2H+": whole(),
                    "2cytc_red": frag((C, "Fe^{2+}"))},
    "complex_iv": {"oxygen": whole(), "4cytc_red": frag((C, "Fe^{2+}")), "4H+": whole(C),
                   "water": frag((C, "H_2"), (M, "O"))},
    "atp_synthase": {"Pi": whole(), "atp_n": hw(u4=[M, "{P}"])},
    "mito_glycerol_3p_dehydrogenase": {"glycerol_3p": chain(r1={"l": [M, "H"], "r": M}), "QH2": frag((M, "H_2"))},
    "alcohol_dehydrogenase": {"ethanol": chain(r0={"u": M, "r": [M, "H"]}), "NADH+H+": NADH_M},
    "aldehyde_dehydrogenase": {"acetaldehyde": chain(r0={"fork_r": M}), "NADH+H+": frag((M, "H")),
                               "H2O": whole(C), "acetate": chain(r0={"fork_r": C})},
    "meos_cyp2e1": {"ethanol": chain(r0={"u": M, "r": [M, "H"]}), "O2": whole(C)},
    "acetyl_coa_synthetase": {"acetate": chain(r0={"fork_r": M}), "AMP": whole(), "CoA-SH": whole(C),
                              "acetyl_coa": {"colors": {0: {"r": C}}}},
    "methanol_dehydrogenase": {"methanol": chain(r0={"u": M, "r": [M, "H"]}), "NADH+H+": NADH_M},
    "formaldehyde_dehydrogenase": {"formaldehyde": chain(r0={"fork_r": M}), "NADH+H+": frag((M, "H")),
                                   "H2O": whole(C), "formate": chain(r0={"fork_r": C})},
    # ------------------------------------------------------------------ round 12: lipid metabolism, extended
    "cyclooxygenase": {"2O2": whole(), "pgg2": hl("[OX2:1][OX2:2]")},
    # peroxidase: the outer O of 15-OOH leaves as water; GSH supplies the electrons
    "cox_peroxidase": {"pgg2": hl("[OX2H1:1][OX2]"), "H2O": whole(), "2GSH": whole(C), "GSSG": whole(C)},
    "prostacyclin_synthase": {"pgh2": hl("[OX2:1][OX2:2]"), "pgi2": hl("[O;R:1]", "[C;R][OH:1]")},
    "thromboxane_synthase": {"pgh2": hl("[OX2:1][OX2:2]"), "txa2": hl("[O;R:1]")},
    "pge_synthase": {"pgh2": hl("[OX2:1][OX2:2]"), "pge2": hl("[C;R]=[O:1]", "[C;R][OH:1]")},
    "lipoxygenase_5": {"O2": whole(), "hpete5": hl("[O:1][OH:2]")},
    # 5-HPETE → LTA4: the outer O leaves as water (magenta), the inner one becomes the epoxide O (cyan)
    "lta4_synthase": {"hpete5": {"highlight": [{"smarts": "[OH1:1]O", "color": M}, {"smarts": "C[O:1][OH]", "color": C}]},
                      "H2O": whole(), "lta4": hl("[O;R:1]", col=C)},
    # LTA4 hydrolase: water adds the 12-OH (magenta); the epoxide O ends up as the 5-OH (cyan)
    "lta4_hydrolase": {"lta4": hl("[O;R:1]", col=C), "H2O": whole(),
                       "ltb4": {"highlight": [{"smarts": "[OH:1][CH1]([CH2][CH1]=[CH1])[CH1]=[CH1]", "color": M},
                                              {"smarts": "[OH:1][CH1]([CH2][CH2])[CH1]=[CH1]", "color": C}]}},
    "ltc4_synthase": {"GSH": whole(), "lta4": hl("[O;R:1]", col=C),
                      "ltc4": {"highlight": [{"smarts": "[S:1][#0:2]", "color": M}, {"smarts": "[OH:1][CH1]([CH2])[CH1]S", "color": C}]}},
    # CYP11A1: C22–C27 of cholesterol leave as isocaproaldehyde (cyan); the new C20 ketone O comes from O2
    "cyp11a1": {"cholesterol": ring_atoms(([0, 1, 2, 3, 4, 5], C)), "isocaproaldehyde": {"smarts": "CC(C)CCC=O", "color": C},
                "3O2": whole(), "pregnenolone": hl("[CX3]=[O:1]")},
    "hsd3b": {"pregnenolone": hl("[C;R][OH:1]"), "NADH+H+": NADH_M, "progesterone": hl("[C;R](=[O:1])C=C")},
    "cyp17_hydroxylase": {"O2": whole(), "oh_progesterone": hl("[C;R](C(C)=O)[OH:1]")},
    "cyp17_lyase": {"O2": whole(), "oh_pregnenolone": {"highlight": [{"smarts": "[CH3:1][C:2](=[O:3])", "color": C},
                                                                      {"smarts": "[C;R](C(C)=O)[OH:1]", "color": M}]},
                    "dhea": hl("[C;R]=[O:1]"), "acetate": chain(r0=C, r1=C)},
    "cyp21_doc": {"O2": whole(), "doc": hl("C(=O)C[OH:1]")},
    "cyp21_deoxycortisol": {"O2": whole(), "deoxycortisol": hl("C(=O)C[OH:1]")},
    "cyp11b1": {"O2": whole(), "cortisol": hl("[CH1;R][OH:1]")},
    "cyp11b_corticosterone": {"O2": whole(), "corticosterone": hl("[CH1;R][OH:1]")},
    # aldosterone synthase: the C18 methyl (cyan) becomes the C18 aldehyde
    "aldosterone_synthase": {"corticosterone": hl("[CH3:1][C;R:2]~[C;R]C(=O)CO", col=C), "2O2": whole(C),
                             "aldosterone": hl("[CH1:1](=[O:2])[C;R:3]", col=C)},
    "hsd17b": {"NADPH+H+": NADH_M, "testosterone": hl("[CH1;R][OH:1]")},
    # aromatase: the C19 methyl (cyan) leaves as formate
    "aromatase": {"testosterone": hl("[CH3:1][C;R:2]C=C", col=C), "formate": chain(r0=C), "3O2": whole()},
    "srd5a": {"testosterone": hl("[C;R:1]=[CH1;R:2]C=O"), "NADPH+H+": NADH_M, "dht": hl("[CH1;R:2][CH2;R:1]C=O")},
    "cyp7a1": {"O2": whole(), "oh7_cholesterol": hl("[CH1;R]([OH:1])C=C")},
    "bile_7a_dehydroxylation": {"cholate": hl("[OH:1][CH1;R]([CH2;R])[CH1;R]"), "H2O": whole(), "NADH+H+": NADH_C},
    "bile_7a_dehydroxylation_cdca": {"chenodeoxycholate": hl("[OH:1][CH1;R]([CH2;R])[CH1;R]"), "H2O": whole(), "NADH+H+": NADH_C},
    # SPT: Ser's COO⁻ leaves as CO2 (magenta); the CoA leaves the palmitoyl group (cyan)
    "serine_palmitoyltransferase": {"serine": chain(r0=M), "co2": whole(), "palmitoyl_coa": {"colors": {"o2": C, "rest": C}},
                                    "CoA-SH": whole(C)},
    "ketosphinganine_reductase": {"NADPH+H+": NADH_M, "sphinganine": hl("[CH1][OH:1]")},
    # lysosomal hydrolases: water (magenta) ends up on the released sugar / phosphate / sulfate
    "glucocerebrosidase": {"H2O": whole(), "b_glucopyranose": hw(u1=M)},
    "galactocerebrosidase": {"H2O": whole(), "b_galactopyranose": hw(u1=M)},
    "beta_galactosidase_gm1": {"H2O": whole(), "b_galactopyranose": hw(u1=M)},
    "alpha_galactosidase_a": {"H2O": whole(), "a_galactopyranose": hw(d1=M)},
    "hexosaminidase_a": {"H2O": whole()},
    "sphingomyelinase": {"H2O": whole(), "phosphocholine": {"smarts": "P(=O)([O-])[O-]", "color": M}},
    # ------------------------------------------------------------------ round 13: peptides
    # leaving O⁻ of Gly + 2 H of Ala's NH₃⁺ (magenta) → H₂O; the new C(=O)–N peptide bond in cyan
    "dipeptide_formation": {"glycine": chain(r0={"c": [[C, "CO"], [M, "O^-"]]}),
                            "alanine": chain(r1={"l": [[M, "^+H_3"], [C, "N"]]}),
                            "glycylalanine": {"colors": {"c1": C, "o1": C, "n2": C}}, "H2O": whole()},
    # ------------------------------------------------------------------ round 14: glyoxylate / Calvin
    # isocitrate's lower two carbons (magenta) leave as glyoxylate
    "isocitrate_lyase": {"isocitrate": chain(r3=M, r4=M), "glyoxylate": chain(r0=M, r1=M)},
    # glyoxylate (magenta) → malate C1–C2; acetyl (cyan) → malate C3–C4
    "malate_synthase": {"glyoxylate": chain(r0=M, r1=M), "acetyl_coa": {"colors": {0: {"c": C, "l": C, "u": C}}},
                        "malate": chain(r0=M, r1=M, r2=C, r3=C)},
    # fixed CO₂ (magenta) becomes the carboxylate of one of the two 3-PG
    "rubisco_carboxylase": {"CO2": whole(), "pg_3": chain(r0=M)},
    # one O of O₂ (magenta) ends up in the carboxylate of 2-phosphoglycolate
    "rubisco_oxygenase": {"O2": whole(), "phosphoglycolate_2": chain(r0={"fork_r": M})},
    "arylsulfatase_a": {"H2O": whole(), "sulfate": whole(), "sulfatide": {"colors": {"head": [M, "Sulfo"]}}},
}


def apply(data):
    """Write COLORS into data/reactions/<id>.yaml, replacing the `decorate:` block and keeping the rest."""
    for rid, deco in COLORS.items():
        p = Path(data) / "reactions" / f"{rid}.yaml"
        lines = p.read_text(encoding="utf-8").splitlines()
        out, i = [], 0
        while i < len(lines):
            if lines[i].startswith("decorate:"):
                i += 1
                while i < len(lines) and (lines[i].startswith(" ") or not lines[i].strip()):
                    i += 1
                continue
            out.append(lines[i])
            i += 1
        block = yaml.safe_dump({"decorate": deco}, allow_unicode=True, sort_keys=False, default_flow_style=None, width=120)
        at = next((k for k, ln in enumerate(out) if ln.startswith(("variants:", "card_zh:", "cofactors_inline:", "balance_check:", "layout:", "tags:", "status:"))), len(out))
        out[at:at] = block.rstrip("\n").splitlines()
        p.write_text("\n".join(out) + "\n", encoding="utf-8")
    print(f"{len(COLORS)} reactions coloured")


if __name__ == "__main__":
    apply(sys.argv[1] if len(sys.argv) > 1 else "data")
