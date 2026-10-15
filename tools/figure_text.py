"""Short on-figure wording for the illustrated scenes of rounds 18–21 (user, 2026-10-13: less descriptive text on the
pictures). Each long original text / label maps to a short replacement, or to None = take it off the picture and
append it to notes_zh (【補充】). Applied by tidy() in every molecular-biology generator's writer.
"""

SHORT = {
    # alternative_splicing
    'Calcitonin / CGRP gene:\nalternative poly(A) site + splicing': 'Calcitonin / CGRP gene',
    # bacterial_transcription_initiation
    "1. Closed complex: holoenzyme (α_2ββ'ωσ)\n    σ reads −35 (TTGACA) and −10 (TATAAT)": '1. Closed complex',
    '2. Open complex: DNA melted from −10 to +2\n    (~13 bp bubble; no helicase needed)': '2. Open complex (−10 → +2)',
    '3. Abortive initiation: short RNAs (2–9 nt)\n    are made and released again and again': '3. Abortive initiation',
    '4. Promoter escape (RNA ≥ ~10 nt): σ leaves;\n    core (+ NusA) elongates; σ is reused': '4. Promoter escape; σ leaves',
    # base_excision_repair
    'DNA glycosylase (e.g. uracil-DNA glycosylase)\ncleaves the N-glycosidic bond': 'DNA glycosylase',
    "AP endonuclease (APE1)\nnicks 5' of the AP site": 'AP endonuclease (APE1)',
    'Pol β (eukaryotes) / Pol I (E. coli):\ninserts 1 nt, removes dRP': 'Pol β / Pol I',
    # chromatin_epigenetics
    'Genomic imprinting (one parental allele silenced)': 'Genomic imprinting',
    # chromatin_levels, dna_double_helix, dna_polymerase_mechanism, eukaryotic_origin_licensing: paragraphs → notes
    'Nucleosome core: histone octamer (H2A, H2B, H3, H4) × 2 + ~147 bp DNA wrapped 1.7 turns; histones rich in Lys/Arg (+).\nEuchromatin: open, acetylated, active.  Heterochromatin: condensed, H3K9me3, inactive.': None,
    'B-DNA: right-handed, antiparallel; bases stacked 0.34 nm apart, perpendicular to the axis;\nA–T 2 H-bonds, G≡C 3 H-bonds; proteins read the sequence mainly in the major groove.\nA-DNA: right-handed, wider, 11 bp/turn (dehydrated; RNA–DNA hybrids, dsRNA).   Z-DNA: left-handed zigzag, 12 bp/turn (alternating GC, methylated).': None,
    'Palm (catalytic Asp + 2 Mg^{2+}) · fingers (bind incoming dNTP) · thumb (holds DNA)': 'palm · fingers · thumb',
    "New strand grows 5'→3' only; template read 3'→5'. A primer with a free 3'-OH is required (no de novo start; primase/RNA pol can start).\nProofreading: mismatch → 3' end shifts to the 3'→5' exonuclease site (Pol I, III, δ, ε; absent in Pol α, β, RT, Taq). Error ~10^{-9} with MMR.": None,
    'Licensing only in G1 (low CDK); in S/G2 high CDK and geminin (blocks Cdt1) prevent reloading → each origin fires once per cycle.\nEukaryotes: many origins (ARS in yeast ≈ 11 bp core), replicons; pols α (primase), ε (leading), δ (lagging); PCNA clamp, RFC clamp loader, FEN1 + RNase H remove primers.': None,
    # homologous_recombination
    'RecBCD (helicase + nuclease)\nresects until chi (GCTGGTGG)': 'RecBCD resects to chi',
    "3' single-strand tail coated by RecA\n(eukaryotes: Rad51, loaded with BRCA2)": 'RecA (Rad51) filament',
    "strand invasion → D-loop;\nthe invading 3' end is extended": 'strand invasion → D-loop',
    'Holliday junction: RuvA (binds)\n+ RuvB (branch migration, ATP)': 'Holliday junction: RuvAB',
    'RuvC (resolvase) cuts the junction\n→ crossover (or non-crossover) products': 'RuvC resolves',
    # mismatch_repair
    'excised: UvrD (helicase II) + exonuclease, SSB': 'UvrD + exonuclease',
    'DNA Pol III fills the gap;\nDNA ligase (NAD^+) seals': 'Pol III + ligase',
    'Dam methylase (SAM) methylates\nthe new GATC → fully methylated': 'Dam methylates new GATC',
    # mrna_capping
    'RNA triphosphatase removes the γ-phosphate (→ P_i)': 'RNA triphosphatase (− P_i)',
    'guanylyltransferase: GTP → GMP joined 5′→5′ (→ PP_i)': 'guanylyltransferase (GTP → PP_i)',
    '2′-O-methyltransferase (SAM) on nucleotide 1 (→ cap 1; nt 2 too → cap 2)': '2′-O-methyltransferase (SAM)',
    # mrna_prok_vs_euk
    'Prokaryote (E. coli): polycistronic mRNA, translated while it is still being made': 'Prokaryote: polycistronic',
    'each ORF has its own Shine–Dalgarno + AUG; ribosomes (purple) load before transcription ends': None,
    'Eukaryote: monocistronic; processed in the nucleus, translated in the cytoplasm': 'Eukaryote: monocistronic',
    # nucleotide_excision_repair
    'UvrA_2B: scans, recognises\nthe distortion (ATP)': 'UvrA_2B recognises',
    "UvrC: dual incision\n(8 nt 5', 4–5 nt 3' of the lesion)": 'UvrC: dual incision',
    'UvrD (helicase II) releases the\n12–13-mer oligonucleotide': 'UvrD removes 12–13-mer',
    # oric_initiation
    'oriC (245 bp): 11 GATC sites — fully methylated (Dam) = ready; hemimethylated → SeqA blocks re-initiation': None,
    'IHF, HU, FIS bend DNA;\ngyrase supplies\nnegative supercoils': 'IHF, HU, FIS',
    # pol2_preinitiation_complex
    'Assembly: TFIID (TBP binds TATA and bends DNA) → TFIIA → TFIIB → Pol II·TFIIF → TFIIE → TFIIH = preinitiation complex (PIC)':
        'TFIID → IIA → IIB → Pol II·IIF → IIE → IIH',
    # polyadenylation
    'downstream RNA degraded by Xrn2 → Pol II terminates': 'Xrn2 → termination',
    'poly(A) polymerase: ATP → AMP + PP_i, no template': 'poly(A) polymerase',
    # replication_fork
    'Topoisomerase II\n(gyrase): relieves\nsupercoils': 'Gyrase',
    'Lagging strand: Okazaki fragments (1000–2000 nt)': 'Lagging strand (Okazaki fragments)',
    # ribosome_structure
    'P site: peptidyl-tRNA (initiator tRNA starts here)': 'P: peptidyl-tRNA',
    'large subunit: peptidyl transferase centre (23S / 28S rRNA = ribozyme), exit tunnel': None,
    'small subunit: binds mRNA first, decoding centre (16S / 18S rRNA)': None,
    '30S (16S rRNA + 21 proteins) + 50S (23S + 5S rRNA + 34 proteins)': '30S (16S) + 50S (23S, 5S)',
    '40S (18S rRNA + ~33 proteins) + 60S (28S, 5.8S, 5S rRNA + ~47 proteins)': '40S (18S) + 60S (28S, 5.8S, 5S)',
    # rna_editing
    'C → U: apolipoprotein B mRNA (APOBEC-1 cytidine deaminase)': 'C → U: apoB mRNA (APOBEC-1)',
    'A → I: ADAR (adenosine deaminase acting on RNA) in double-stranded RNA; I is read as G': 'A → I: ADAR (dsRNA)',
    'GluR-B (AMPA receptor) Q/R site: CAG (Gln) → CIG = CGG (Arg) → Ca^{2+}-impermeable channel; also 5-HT_{2C} receptor; tRNA wobble I (ADAT)': None,
    'U insertion / deletion: trypanosome mitochondria (kinetoplastid); a guide RNA (gRNA) is the template': 'U insertion: trypanosome (gRNA)',
    # rrna_processing
    'E. coli: rrn operon (7 copies) → one 30S pre-rRNA, one RNA polymerase': 'E. coli: 30S pre-rRNA',
    'Eukaryotes: Pol I makes 45S pre-rRNA in the nucleolus; the 5S rRNA gene is transcribed separately by Pol III': 'Eukaryotes: 45S pre-rRNA (Pol I) + 5S (Pol III)',
    # site_specific_recombination
    'Integration / excision (direct repeats) — phage λ': 'Integration / excision (λ)',
    'Inversion (inverted repeats) — e.g. Salmonella flagellar phase variation (Hin), Cre–loxP in opposite orientation': 'Inversion (Hin)',
    # sos_response
    'SOS genes OFF:\nrecA, lexA, uvrA/B, umuDC (Pol V), dinB (Pol IV), sulA': 'SOS genes OFF',
    'RecA* acts as a co-protease:\nLexA self-cleaves → SOS genes ON': 'LexA self-cleaves → ON',
    "Pol V (UmuD'_2C) / Pol IV inserts across\nthe lesion: no proofreading → mutations": 'Pol V / IV: error-prone',
    # spliceosome_splicing
    '1. U1 snRNP pairs with the 5′ splice site (GU);\n    U2 pairs at the branch point (A bulges out);\n    U2AF binds the Py tract and the 3′ AG': '1. U1, U2, U2AF bind',
    'spliceosome: U2, U5, U6 snRNPs + ~150 proteins': 'spliceosome',
    '2. U4/U6·U5 join; U1 and U4 leave; U6 pairs with\n    the 5′ splice site and with U2 → RNA catalytic core.\n    (1) the branch-A 2′-OH attacks the 5′ splice site': '2. branch-A 2′-OH attacks 5′ss',
    '3. Lariat intermediate (2′–5′ link at the branch A);\n    (2) the exon-1 3′-OH attacks the 3′ splice site': '3. exon-1 3′-OH attacks 3′ss',
    '4. Exons joined (mRNA); lariat released,\n    debranched and degraded': '4. exons joined + lariat',
    # splicing_mechanisms
    'linear intron (5′ G added)\nself-splicing ribozyme:\nTetrahymena pre-rRNA, organelles': 'linear intron',
    'lariat intron\nself-splicing ribozyme:\nfungal / plant organelles': 'lariat intron',
    'lariat intron\nsnRNPs; catalytic core = U2/U6 RNA:\nnuclear pre-mRNA': 'lariat intron',
    # srp_er_targeting
    'SRP·SR (both GTP) dock on Sec61;\nGTP hydrolysis releases SRP': 'docks on Sec61',
    'translation resumes; signal peptidase cuts; OST adds the N-glycan (Asn-X-Ser/Thr) from dolichol-PP; BiP folds': None,
    # telomere_replication
    'C-rich strand: lagging strand cannot be finished\n(last RNA primer removed) → ends shorten each division': 'end-replication problem',
    'Telomerase (TERT reverse\ntranscriptase + RNA template)': 'Telomerase',
    "extends the 3' end using its own RNA template (3'-AAUCCC-5')": "extends the 3' end",
    'Then Pol α–primase fills in the C-rich strand; the end folds into a T-loop (D-loop) or G-quadruplex, protected by shelterin (TRF1/2, POT1).\nTelomerase active in germ cells, stem cells and ~90 % of cancers; somatic cells lack it → senescence. Telomere 10–20 kb, overhang > 100 nt.': None,
    # topoisomerase_types
    'cuts BOTH strands (G segment);\nT segment passes through': 'cuts both strands',
    'Inhibitors: camptothecin, irinotecan, topotecan (topo I); etoposide, doxorubicin (topo II);\nquinolones / fluoroquinolones (bacterial gyrase, topo IV); novobiocin (gyrase ATPase).': None,
    # transcription_termination
    'pause at the hairpin + weak\nhybrid → RNA released': 'RNA released',
    'ρ (hexameric RNA helicase) moves\n5′→3′ along the RNA (ATP)': 'ρ moves 5′→3′ (ATP)',
    # translation_elongation
    'codon–anticodon OK →\nGTP hydrolysis;\nEF-Tu·GDP leaves': 'EF-Tu·GDP leaves',
    '2. peptidyl transfer:\nα-NH_2 of aa-tRNA (A)\nattacks the peptidyl\nester (P); 23S rRNA': '2. peptidyl\ntransfer',
    '4. tRNA leaves\nfrom E; A site\nfree for the\nnext aa-tRNA': '4. E-site\ntRNA leaves',
    # translation_inhibitors
    'Macrolides (erythromycin): 50S exit\ntunnel, block elongation / translocation': 'Macrolides (exit tunnel)',
    'Chloramphenicol, linezolid, clindamycin\n(50S peptidyl transferase centre)': 'Chloramphenicol (PTC)',
    'Puromycin: mimics Tyr-tRNA end, enters A,\naccepts the peptide → premature release': 'Puromycin (A site)',
    'Tetracyclines (30S A site):\nblock aa-tRNA binding': 'Tetracyclines (30S A site)',
    'Fusidic acid: EF-G·GDP stuck;\ncycloheximide (80S): translocation': 'Fusidic acid (EF-G)',
    'Streptomycin, aminoglycosides (30S):\nmisreading; high dose blocks initiation': 'Aminoglycosides (30S)',
    # translation_initiation
    '70S initiation complex: fMet-tRNA in P, A site empty': '70S initiation complex',
    # translation_termination
    'RF1 reads UAA, UAG; RF2 reads UAA, UGA (eRF1 reads all three). Its GGQ motif puts water into the peptidyl transferase centre. RF3·GTP (eRF3) removes RF1/2.': None,
    'Rescue of stalled ribosomes: trans-translation by tmRNA (bacteria)': 'tmRNA rescue (bacteria)',
    'ribosome switches to the tmRNA reading frame:\n11-aa ssrA tag (AANDENYALAA) + stop': 'ssrA tag added',
    'tagged protein → ClpXP protease\nfaulty mRNA → RNase R': 'protein + mRNA degraded',
    # transposition
    'Bacterial transposons (inverted repeats ▲ at both ends)': 'Bacterial transposons',
    'cut-and-paste (non-replicative): donor left broken': 'cut-and-paste',
    'copy-and-paste (replicative): donor keeps a copy': 'copy-and-paste',
    'yellow = target-site duplication (short direct repeats)': 'yellow = target-site duplication',
    # trna_processing
    '3′ trailer: RNase Z (endo, eukaryotes)\nor RNase D / PH / T (exo, E. coli)': 'RNase Z / RNase D',
    'intron in the anticodon loop:\nsplicing endonuclease + ligase': 'intron (eukaryotes)',
    '3′ CCA added by tRNA nucleotidyl-\ntransferase (no template; encoded in E. coli)': 'CCA added',
    'modified bases: D (dihydro-U),\nΨ, T (ribothymidine), I, m^1G …': 'modified bases (D, Ψ, T, I)',
    # ubiquitin_proteasome
    '19S: binds the Ub chain, removes Ub (DUB),\nunfolds the protein (ATP)': '19S: unfolds (ATP)',
    '20S core (α7β7β7α7): β1 caspase-like,\nβ2 trypsin-like, β5 chymotrypsin-like (N-terminal Thr)': '20S core: β1, β2, β5',
    # vdj_recombination
    'V–DJ joining; intervening DNA lost as a circle': 'V–DJ joining',
    'RSS = heptamer – spacer (12 or 23 bp) – nonamer; 12/23 rule: a 12-RSS pairs only with a 23-RSS. RAG1/2 make hairpin coding ends,\nArtemis opens them, TdT adds N nucleotides → junctional diversity. RAG or Artemis defect → SCID (no B and T cells).': None,
    # labels a little under the 45-character mark that are still sentences
    'pause at the hairpin + weak hybrid → RNA released': 'RNA released',
}


def tidy(d):
    """Shorten on-figure texts of a pathway dict (scene and illus variants); None entries move into notes_zh."""
    extra = []
    scenes = [d["scene"]] if d.get("scene") else []
    scenes += [v["scene"] for v in (d.get("variants") or {}).values() if v.get("scene")]
    for sc in scenes:
        keep = []
        for it in sc["items"]:
            key = "text" if it.get("g") == "text" else "label"
            t = it.get(key)
            if t in SHORT:
                if SHORT[t] is None:
                    extra.append(" ".join(t.split()))
                    if key == "text":
                        continue
                    it = {k: v for k, v in it.items() if k not in ("label", "label_at", "point", "label_align", "label_size")}
                else:
                    it = {**it, key: SHORT[t]}
            keep.append(it)
        sc["items"] = keep
    if extra:
        d["notes_zh"] = d.get("notes_zh", "") + "\n【補充】" + " ".join(extra)
    return d
