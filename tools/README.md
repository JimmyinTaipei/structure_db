# tools/ — 一次性的資料產生器

這些腳本會產生 `data/` 下的 YAML。產生後的 YAML 就是資料本身，平常不需要再跑；要大量修改同一類分子時，改腳本再重跑比逐一改檔安全。

從 `structure_db/` 執行：

```sh
.venv/bin/python tools/gen_nuc.py data/molecules     # 鹼基、核苷、核苷酸（同一個糖模板）
.venv/bin/python tools/gen_sugar.py data             # 醣類結構庫與總表
.venv/bin/python tools/gen_aamet.py data             # 胺基酸代謝補完與 heme
.venv/bin/python tools/gen_energy.py data            # ETC、氧化磷酸化、shuttle、代謝整合、酒精代謝
.venv/bin/python tools/gen_lipext.py data            # 類二十烷酸、固醇類荷爾蒙、膽酸、鞘脂與儲積症（SMILES 取自 PubChem）
.venv/bin/python tools/gen_plots.py data             # 作圖：酵素動力學、抑制、O₂ 結合曲線、Hill plot、滴定、Ramachandran；胜肽鍵、GSH
.venv/bin/python tools/gen_misc.py data              # GAG 重複雙醣、chitin、peptidoglycan、黏多醣症；乙醛酸循環；Calvin 循環與光呼吸
.venv/bin/python tools/gen_notes.py data             # 依上課筆記補齊：抑制／動力學作圖、Hb T/R、滴定、melting、PTM、鹼基配對、PUFA、類二十烷酸系列…
.venv/bin/python tools/gen_met.py data               # 依代謝筆記補齊：energy charge、HK vs GK、糖解 ΔG、磷酸轉移勢、發酵、葡萄糖碳去向、PPP 非氧化期、醣互換、解毒…
.venv/bin/python tools/gen_sig.py data               # 訊號傳遞與細胞膜（Book 3）：cAMP／cGMP 與其環化酶、PDE，GPCR、RTK、JAK-STAT、NF-κB、凋亡、p53、毒素、膜運輸、流動性…
.venv/bin/python tools/gen_mb_rep.py data            # 分生 I：DNA 結構與複製（插畫 scene：複製叉…）
.venv/bin/python tools/gen_mb_repair.py data         # 分生 II：DNA 損傷、修復（MMR、BER、NER、HR、NHEJ、SOS）、重組（site-specific、轉位、V(D)J）、MGMT
.venv/bin/python tools/gen_mb_tx.py data             # 分生 III：轉錄（σ cycle、延長、終止、Pol II PIC）、RNA 加工（cap、poly(A)、剪接、editing、tRNA／rRNA）
.venv/bin/python tools/gen_mb_tl.py data             # 分生 IV：轉譯（密碼表、wobble、aaRS、核糖體、起始／延長／終止、抗生素）、蛋白質命運（SRP、泛素、修飾）
.venv/bin/python tools/gen_mb_reg.py data            # 分生 V：基因調控（operon、attenuation、ppGpp、sRNA/riboswitch、Hin、真核活化、DNA 結合 motif、表觀遺傳、RNAi、核受體）
.venv/bin/python tools/gen_mb_tech.py data           # 分生 VI：技術（PCR、qPCR、Sanger、blotting、選殖與篩選、cDNA、CRISPR、EMSA/footprint/ChIP、Y2H、RFLP、X-gal）
.venv/bin/python tools/gen_illus_sig.py data        # 圖片版：細胞膜、訊號傳遞、能量代謝 18 張 + 膜蛋白分類、通透性（只加 variants.illus；gen_sig / gen_notes / gen_energy 重跑後要再跑一次）
.venv/bin/python tools/gen_gap.py data              # 生化筆記缺口：serine protease 兩步反應與三聯體、GroEL、PDH 複合體、紅血球 GSH、THF 一碳單位、phospholipase 切點、純化、定序、飢餓曲線
# tools/figure_text.py：第 18–21 輪圖上長字的縮短對照表（None = 搬進 notes_zh）
.venv/bin/python tools/reaction_colors.py data       # 反應的原子追蹤配色（產生器也會讀）
.venv/bin/python tools/gen_lipmet.py data            # 脂質代謝：β 氧化、合成、酮體、膽固醇、TAG、脂蛋白
.venv/bin/python tools/gen_carb.py data              # 醣類代謝：糖質新生、肝醣、PPP、果糖、半乳糖、多元醇
.venv/bin/python tools/gen_vit.py data               # 維生素與輔酶、活化反應、總表（含 THF 家族）
.venv/bin/python -m sdb validate                     # 一定要跑
```

- `fischer_smiles.py`：依 Fischer 圖（OH 在左或右）算出含立體化學的 SMILES，用 D-glucose、D-fructose、D-ribose 對照過 PubChem
- `sugar_lib.py`：PubChem 的環狀單醣 SMILES，以及 `glycosyl()`。`glycosyl()` 產生以 anomeric OH 為根的醣基片段，拼雙醣時立體標記不會跑掉
- 注意：重跑會覆寫對應的 YAML，包括人工改過的內容。審核意見請改在腳本裡
