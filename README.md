# structure_db：生化／有機結構式資料庫

結構式和反應都存成 YAML 資料，由程式自動畫出風格一致的圖。風格參考 Lehninger，參考圖在 `reference/`。

- **SVG** 放在 `out/svg/`：給 Typora、VS Code 筆記用，匯出 PDF 時是向量圖，不會糊
- **PNG** 放在 `out/png/`：3 倍解析度、透明背景，可以直接當 GoodNotes 貼紙
- **本機網頁**：可以搜尋、點圖跳轉、看中文資料卡，也可以在上面審核

## 打開網頁

側欄標籤分兩層：例如點「胺基酸」才會出現「非極性、極性、鹼性…」子標籤。上下層關係寫在 `data/tags.yaml`。
反應圖上可以跳轉的分子，滑鼠移上去會出現青色虛線框和名稱提示；圖下方也有一排「跳到：」按鈕。


在終端機輸入：

```bash
cd ~/Projects/structure_db
.venv/bin/python -m sdb serve
```

然後在瀏覽器打開 <http://localhost:8765>。要關掉時，在終端機按 Ctrl+C。

## 程式和圖片放在不同的地方

- **程式和資料**：放在 `~/Projects/structure_db`（git repo，不在 iCloud 裡，避免出現「xxx 2」這種同步衝突複本）
- **產出的圖**：還是寫到 iCloud 的 `82_考題整理/structure_db/out/`，所以 Typora 筆記的相對路徑和 iPad GoodNotes 都照常可以用。要改位置的話，設定環境變數 `SDB_OUT`（見 `sdb/style.py`）

## 外面看的網頁（GitHub Pages）

```bash
.venv/bin/python -m sdb build        # 先出圖（寫到 iCloud 的 out/）
tools/deploy_site.sh                 # 產出 site/ 並推到 gh-pages branch
```

`site/` 只放網頁、`index.json` 和 SVG（約 70 MB），**不放** PNG、拆件、審核紀錄，也不放 `reference/` 的教科書截圖。網頁是唯讀的：沒有審核佇列，也不能重新出圖，這些都要在家裡用 `sdb serve` 操作。

- **分子頁**：可以切換「縮寫 / Fischer」和「骨架式」兩種畫法，並列出這個分子參與的所有反應
- **反應頁**：大圖加上中文資料卡（位置、ΔG、輔因子、調控、缺乏症、考點）。點圖上的分子，會跳到該分子的頁面。TCA 的反應有「一般版／同位素版」可以切換：同位素版用 * ○ △ 追蹤 acetyl-CoA 的碳，IDH 另外標 ①②③
- **路徑頁**：點箭頭，會跳到那一步反應的細部圖
- **審核佇列**：每一項都可以按「✓ 核准」，或寫下意見再按「送出意見」

## 審核流程

1. 打開網頁，點「審核佇列」，檢查每一項的**化學正確性**和**畫法**
2. 正確就按「核准」；要修改就寫下意見，例如「NADH 的箭頭太短」或「C3 的 OH 應該在左邊」
3. 對 AI 說「**處理審核意見**」，AI 會修改內容、重新出圖，再把項目放回審核佇列

## 在 Typora 筆記裡引用圖片

在網頁上按「**複製 Markdown 引用**」，直接貼進筆記就可以。路徑是以 `typora_file/生化/`（或其他科目資料夾）為起點的相對路徑：

```html
<img src="../../structure_db/out/svg/reactions/pdh.svg" alt="Pyruvate → Acetyl-CoA">
```

想先看整理好的圖庫，可以用 Typora 打開 `out/md/<標籤>.md`，例如 `out/md/PPP.md`、`out/md/TCA.md`。

## 放進 GoodNotes

這個資料夾在 iCloud 裡，所以 iPad 的「檔案」App 也找得到：
iCloud Drive → 後醫 → 08_考古題與模擬考 → 82_考題整理 → structure_db → out → png。
把 PNG 拖進 GoodNotes，存成貼紙就可以用。

## 圖的畫法

- **上色**：反應中有變化的部分塗洋紅，其餘維持黑色，不用色塊
- **磷酸**：一律畫成 Ⓟ
- **TCA**：用老師的直線排法，同一個原子在反應前後留在同一個位置
- **檔名**：同位素版的檔名多一個 `.isotope`，例如 `out/png/reactions/citrate_synthase.isotope.png`

## 胺基酸

- **20 個胺基酸**：照 Lehninger 畫，R 基有鮭魚色底色。另外照 Mathews 替特殊元素上色：鹼性 N 是青藍、酸性 COO⁻ 是洋紅、S 有黃色底色塊
- **底色可以分開**：每個胺基酸的拆件裡，`structure` 是沒有底色的版本，`box_R`、`box_S` 各是一層獨立的底色
- **總表**：`pathways/amino_acids`，依 R 基分組
- **胺基酸代謝**：ALT、AST、GDH、glutamine synthetase、glutaminase、尿素循環 5 步，以及尿素循環、glucose–alanine cycle 的總圖

## 核苷酸

- **結構**：
  - 核苷酸用 Haworth 畫（Lehninger 風格），鹼基接在 C1′
  - 嘌呤的六員環在右，嘧啶的 C2=O 在右
  - 鹼基單獨也有一頁
- **原子來源圖**：`purine_atom_sources`、`pyrimidine_atom_sources`，每個原子依來源上色
- **反應**：
  - PRPP、嘌呤從頭合成（限速步驟與 IMP 之後）
  - 嘧啶從頭合成 6 步、CTP synthetase
  - RNR、thymidylate synthase、DHFR
  - 救援途徑（HGPRT、APRT）
  - 嘌呤分解
- **總圖**：
  - 嘌呤從頭合成（名稱版）
  - 嘧啶從頭合成
  - 去氧核苷酸與 dTMP（含藥物作用點）
  - 嘌呤分解與救援

## 醣類

- **單醣（Fischer）**：
  - 所有 D-醛糖（3–6 碳）、D-酮糖（3–6 碳）
  - L-glyceraldehyde、L-glucose
  - sorbitol、galactitol
- **環狀**：
  - α/β-D-glucose（Haworth 和椅式）
  - β-galactose、α-mannose、α/β-fructofuranose
  - β-ribose、2-deoxyribose
- **衍生物**：glucuronate、glucosamine、GlcNAc
- **雙醣**：maltose、lactose、sucrose、cellobiose，醣苷鍵的位置與 α/β 都經過驗證
- **多醣示意**：glycogen、amylose、cellulose，用 SNFG 符號（藍色圓 = glucose）
- **總表**：
  - D-醛糖家族樹、D-酮糖家族樹
  - epimer（mannose · glucose · galactose）
  - D/L 構型
  - 葡萄糖 α/β 互變
  - 雙醣

## 鹼基

- **畫法**：照課本卡片的骨架式。嘌呤的 C4–C5 鍵垂直、六員環在右、N9–H 朝正下方；嘧啶的 N1–H 在下。
- **兩個版本**：每個鹼基都有一般和左右相反兩版；網頁上可切換，拆件中各一張。

## 脂質

- **脂肪酸**：palmitate、stearate、palmitoleate、oleate、elaidate（trans）、linoleate、α-linolenate、arachidonate、EPA、DHA。附 Δ 與 ω 編號。
- **TAG、甘油磷脂**：PA、PC、PE、PS、PI、PG。
- **鞘脂**：ceramide、sphingomyelin、galactocerebroside。
- **固醇類**：cholesterol、cortisol、aldosterone、testosterone、estradiol、progesterone、cholate，標出 A–D 環。
- **類二十烷酸**：PGE₂。
- **總表**：脂肪酸（飽和、cis／trans、必需、長鏈多元不飽和）、膜脂、固醇類。

## 維生素與輔酶

- **水溶性**：B1 thiamine → TPP、B2 riboflavin → FMN → FAD／FADH₂、B3 niacin／nicotinamide → NAD⁺／NADH／NADP⁺／NADPH、B5 pantothenate → CoA、B6 pyridoxine／pyridoxal／pyridoxamine → PLP／PMP、B7 biotin、B9 folate → DHF → THF、B12（示意圖）、C ascorbate／dehydroascorbate；另有 lipoate、SAM。
- **脂溶性**：A（β-carotene、retinol、all-trans／11-cis-retinal、retinoic acid）、D（7-DHC、D₃、25-OH-D₃、calcitriol）、E（α-tocopherol）、K（K₁、KH₂、epoxide、menadione）。
- **洋紅**：單獨的分子圖標出反應部位（例如 NAD⁺ 的 nicotinamide 環、TPP 的噻唑 C2、CoA 的 –SH、PLP 的醛基）。反應圖則改標這一步變化的部分。
- **反應**：活化（TPP、FMN、FAD、PLP、folate → DHF）、維生素 D 的活化、β-carotene 切割、retinol ⇌ retinal → retinoic acid。
- **總表**：水溶性維生素 → 輔酶 → 功能對照、維生素 D 活化、維生素 A 的各種形式、維生素 K 循環（warfarin）。
- **輔因子可以點**：反應圖上的 NAD⁺、FAD、CoA-SH、ATP 等，點了會跳到完整結構頁。

## 醣類代謝

側欄的「醣類代謝」底下有：糖解、糖質新生、肝醣代謝、五碳糖磷酸路徑、果糖代謝、半乳糖代謝、多元醇路徑。
- **糖質新生**：pyruvate carboxylase、PEPCK、FBPase-1、G6Pase 四個繞道酵素，加上 LDH、glycerol kinase、glycerol-3-P DH；PFK-2／FBPase-2 與 F2,6BP。
- **肝醣**：phosphoglucomutase、UDP-glucose pyrophosphorylase、glycogen synthase、glycogen phosphorylase（肝醣鏈末端用 Haworth，R 代表其餘部分）。
- **PPP**：G6PD、lactonase、6PG DH、isomerase、epimerase、transaldolase（transketolase 原本就有）。
- **果糖、半乳糖、多元醇**：fructokinase、aldolase B、triose kinase；galactokinase、GALT、UDP-Gal 4-epimerase；aldose reductase、sorbitol DH。
- **總圖**：糖質新生總覽、肝醣合成與分解、PPP 氧化期、果糖代謝、半乳糖代謝（Leloir）、多元醇路徑、Cori cycle。
- 卡片寫了調控與缺乏症（GSD 各型、HFI、galactosemia、G6PD 缺乏…）。

## 脂質代謝

側欄「脂質代謝」底下有：β 氧化、脂肪酸合成、酮體、膽固醇合成、三酸甘油酯代謝、脂蛋白。
- **活化與運輸**：acyl-CoA synthetase、CPT-I（malonyl-CoA 抑制）；palmitoyl-CoA、palmitoylcarnitine 用鋸齒骨架式。
- **β 氧化四步**：用 R–CH₂–CH₂–CO–SCoA 表示任意長度；加上奇數碳（propionyl-CoA → succinyl-CoA，biotin、B12）。
- **脂肪酸合成**：ATP-citrate lyase、malic enzyme、ACC、fatty acid synthase 總反應式（原子平衡檢查通過）。
- **酮體**：thiolase、HMG-CoA synthase、HMG-CoA lyase、β-HB DH、acetone、SCOT（肝外利用）。
- **膽固醇**：HMG-CoA reductase（statin）、mevalonate → IPP、IPP ⇌ DMAPP；squalene 結構。
- **總圖**：活化與肉鹼穿梭、β 氧化、奇數碳、脂肪酸合成、酮體、膽固醇合成、TAG 合成與分解、脂蛋白概念圖。

## 胺基酸代謝（補完）

「胺基酸代謝」底下新增：胺基酸衍生物、碳骨架分解、一碳與甲硫胺酸循環、血基質。
- **衍生物**：Phe → Tyr → DOPA → dopamine → NE → Epi；Trp → 5-HTP → serotonin；His → histamine；Glu → GABA；Arg → NO。
- **分解**：Tyr → HPP → homogentisate（alkaptonuria）；碳骨架進入 TCA 的總圖（生糖／生酮）；BCAA（MSUD）。
- **甲硫胺酸循環**：MAT、SAH hydrolase、methionine synthase（B12）、MTHFR、CBS、cystathionase、SHMT。
- **Heme**：ALA synthase、ALA dehydratase；ALA、PBG、protoporphyrin IX、heme（Fe 配位）、bilirubin 的結構；合成總圖（紫質症、鉛）與分解總圖（Gilbert、Crigler–Najjar）。

## 依主題分類的圖檔 `out/分類/`

每次 build 都會把圖依主題再整理一份。結構是 `out/分類/<主題>/<分子|反應|路徑|拆件>/`，例如 `out/分類/TCA/分子/`、`out/分類/胺基酸/酸性/`。
- **怎麼分**：一個分子會出現在所有用到它的主題。例如 pyruvate 同時在「糖解」和「TCA」。
- **不占空間**：這份是硬連結，不另外占硬碟空間。
- **筆記引用**：Typora 的引用路徑仍是 `out/svg/`，維持不變。

## 路徑的文字版

每張總圖和每個反應都會自動多出兩個版本，網頁上用按鈕切換：
- **文字版**：只有物質名稱和粗箭頭，用來指出代謝路徑
- **文字版＋酵素**：不畫結構式，但保留酵素和輔因子，用來省空間

拆件區有 `enzymes`（酵素名稱層），可以先貼「文字版＋酵素」，再自行拿掉酵素層。

## 拆件貼紙

每個分子和反應都會自動拆成好幾張 PNG，放在 `out/parts/` 裡，所有貼紙都用同一個比例：
- **分子**：結構式、名稱、Fischer 的每一列
- **反應**：每個分子、每個名稱、箭頭加輔因子、酶素名稱、標示層
- **路徑**：節點＋箭頭、酵素名稱層、輔因子層

網頁上的分子頁和反應頁有「拆件貼紙」區塊，可以單張下載，也可以下載整張總表 `_sheet.png`。

## 請 AI 新增內容的範例指令

- 「幫我加 TCA 的 citrate synthase，畫法照 Lehninger」
- 「把 glycolysis 10 步都加進資料庫」
- 「pyruvate 加一個 LDH 反應，lactate 用 Fischer 畫」
- 「處理審核意見」

AI 會依照 `AUTHORING.md` 的規範新增內容。每次新增後都會自動做三種化學檢查：
- 反應的原子數和電荷是否平衡
- Fischer 圖上 H/OH 的左右位置是否和立體化學一致
- 每個碳上的 H 數是否正確

## 資料夾說明

| 路徑 | 內容 |
|---|---|
| `data/` | 分子、反應、路徑的 YAML（AI 撰寫） |
| `out/` | 自動產生的圖。**不要手動編輯**，每次 build 都會覆蓋 |
| `out/parts/` | 拆件貼紙 |
| `review/comments.jsonl` | 你在網頁上寫的審核意見 |
| `reference/` | 風格參考圖 |
| `fonts/` | 內建字型：TeX Gyre Schola 與 Source Sans 3，皆為開源授權 |
| `sdb/` | 繪圖引擎程式 |
