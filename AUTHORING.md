# AUTHORING — AI 撰寫規範

給負責新增／修改內容的 AI 看。使用者只負責審核；資料輸入、排版、出圖、修正都由 AI 完成。

## 每次改完都要跑

```bash
cd structure_db
.venv/bin/python -m sdb validate   # 必須 0 錯誤
.venv/bin/python -m sdb build      # 產出 out/
.venv/bin/python -m pytest -q      # 改到 sdb/ 程式時
```
接著**用眼睛看圖**：把新產生的 `out/svg/...` 轉成 PNG 讀進來檢查（有沒有重疊、截斷、顏色是否正確）。只有 validate 通過還不夠。

新增或修改的項目一律設成 `status: pending`，讓它出現在審核佇列。

## 命名

- 檔名＝`id`，用小寫底線：`xylulose_5p`、`acetyl_coa`、`transketolase_1`
- 一個檔案放一個分子／反應／路徑：`data/molecules/<id>.yaml`、`data/reactions/<id>.yaml`、`data/pathways/<id>.yaml`
- `tags` 用自訂標籤，例如 `TCA`、`glycolysis`、`PPP`、`pyruvate`、`organic`。每個標籤會自動產生 `out/md/<tag>.md` 圖庫

## 文字標記語法（所有圖上文字通用）

| 寫法 | 結果 |
|---|---|
| `CH_3`、`E_{12}` | 下標 |
| `NAD^+`、`O^-`、`^{2-}` | 上標（`-` 會自動換成真正的減號） |
| `^{2-}` 緊接在下標後面，例如 `PO_3^{2-}` | 下標後面緊接上標時，上下疊排 |
| `{P}`，例如 `CH_2O{P}` | 磷酸符號 Ⓟ（淺橙底、黑框的圓）。**所有磷酸基一律寫成 Ⓟ**，不寫 `OPO_3^{2-}` |
| `\n`（在 YAML 雙引號字串中） | 換行 |

圖上一律寫英文；中文只放在 `notes_zh` 和 `card_zh`。

## 分子 `data/molecules/<id>.yaml`

```yaml
id: xylulose_5p
name: Xylulose 5-phosphate
label: "Xylulose\n5-phosphate"   # 選填：圖上名稱要換行時才寫
aliases: [X5P, 木酮糖-5-磷酸]       # 搜尋用，中英文都可以
smiles: "C(O)C(=O)[C@@H](O)[C@H](O)COP(=O)([O-])[O-]"
chain:                             # 直立鏈：縮寫結構式或 Fischer 投影（兩者共用一種格式）
  - {c: CH_2OH}
  - {c: C, r: O, rb: 2}            # rb/lb: 2 = 雙鍵
  - {c: C, l: HO, r: H}
  - {c: C, l: H, r: OH}
  - {c: "CH_2O{P}"}
chains:                            # 選填：同一分子的其他排法，反應裡用 {id: ..., chain: cs} 指定
  cs: [...]
formula: CO_2                      # 選填：只用文字顯示的小分子（CO₂、H₂O）
tags: [PPP]
notes_zh: 一兩句中文重點
status: pending
```

`chain` 每一列可以用的欄位：
- `c`：中心基團
- `l`、`r`：左、右取代基
- `lb`、`rb`：左、右的鍵級
- 每一列的「主原子」：中心文字裡有 C 就是第一個 C，否則是第一個不是 H 的元素（`NH` → N、`^+NH_3` → N、`S` → S、`OH` → O）。所以主鏈可以經過 N、S、O，例如 Arg、Lys、Met、citrulline
- 一列寫 `R`：代表「其餘的碳鏈」，對應 SMILES 的虛擬原子 `[8*]`，例如 β 氧化的 `R / CH_2 / CH_2 / C(=O)S-CoA`。旁邊有立體中心時，SMILES 請由真實分子（例如 3-hydroxybutyryl-CoA）把 CH₃ 換成 `[8*]`，空間排列才不會變
- `dfork`：只用在最後一列，寫成 `{l, r, lb, rb}`，往左下、右下 45° 分出取代基，例如 Val 的 `CH(CH₃)₂`、Asn 的醯胺 `{l: H_2N, r: O, rb: 2}`
- `ring`：只能當最後一列，例如 `{ring: "[*]c1ccc(O)cc1"}`。`[*]` 是接到上一列的位置，環會畫在正下方。驗證會檢查這個環片段確實存在，而且接在主鏈最後一個原子上
- `u`、`d`：畫在中心原子正上方、正下方的取代基，鍵級用 `ub`、`db`。只能用在上下沒有接主鏈的列，例如橫式 acetyl-CoA 寫成 `{c: C, l: CH_3, r: S-CoA, u: O, ub: 2}`
- `fork`：只用在第一列，寫成 `{l, r, lb, rb}`，表示往左上、右上 45° 分出的取代基。醛基寫 `{l: O, lb: 2, r: H}`，羧酸根寫 `{l: O, lb: 2, r: "O^-"}`
- `anchor`：鍵要接在 `c` 的第幾個字元。預設接第一個 C
- `down`：和下一列之間的鍵級

**主鏈畫不出來的分子**（例如 Pro 的環含 Cα 和 N）：改用 `free`，直接給座標。

```yaml
free:
  atoms: {a: ["COO^-", 0, -0.15], ca: [C, 0, 1], n: ["^+H_2N", -1.15, 1], ...}   # [文字, x, y]，單位是主鏈列距
  bonds: [[a, ca], [ca, n], ...]                                                 # 第 3 項：1、2、bold，或 hb（灰色虛線氫鍵：鹼基配對、催化三元組）
```
atoms 可以有第 4 項：標籤的錨點字元索引。例如 `[CoA, x, y, 2]` 讓「CoA」接在 A 上、文字往左延伸。
`free` 只驗證分子式和電荷，不驗證立體化學，請在 `notes_zh` 註明。

**Haworth（核苷酸、PRPP、之後的醣類）**：用 `haworth`，預設檢視就是它。
```yaml
haworth:
  ring: furanose                    # 或 pyranose
  subs:                             # 環上位置（C1 = 接環 O 右邊的碳）→ 取代基在環平面上方或下方
    1: {up: base, down: H}          # base = 下面的 base 環；也可以寫 "Base"（通用鹼基，SMILES 用 [*]）
    2: {up: H, down: OH}            # 去氧核糖寫 {up: H, down: H}
    3: {up: H, down: OH}
    4: {up: "{P}OCH_2", down: H}    # 取代基文字的「主原子」（第一個 C，否則第一個非 H 元素）放在鍵的末端
  base: "[*]n1cnc2[c](N)ncnc21"     # 鹼基 SMILES，[*] 接 C1′
  base_right: 5                     # 選填：這個原子（base 的索引）放右邊；嘌呤放 C6、嘧啶放 C2
  base_len: 0.95                    # 選填：C1′–鹼基的鍵長（列距），鹼基撞到 5′ 文字時加長
  numbers: true                     # 灰色 1′–4′，另外輸出成 numbers 圖層
```
- 驗證會做三件事：
  1. 逐位置比對取代基元素
  2. 把環放進 3D（上 = +y），比對每個立體中心的 CIP。α/β、核糖與阿拉伯糖、D/L 畫錯都會報錯
  3. 鹼基片段必須接在 C1′
- 報錯時，**以教科書的 Haworth 圖為準去修 SMILES**。PRPP 就是這樣抓到 α/β 寫反的
- 核苷酸請用一致的糖 SMILES 模板，鹼基寫在前面：`<鹼基…n9>[C@@H]1O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]1O`

**環平面內的取代基**：`subs: {1: {out: O, out_b: 2}}` 從環碳往外（環平面內）畫一個取代基，`out_b: 2` 為雙鍵，例如 6-phosphoglucono-δ-lactone 的 C1=O。

**多個環（雙醣）**：`haworth: {units: [環, 環], links: [...]}`。
- 每個環可以用以下欄位：
  - `at: [dx, dy]`：位置（列距）
  - `start: 2`：第一個環碳的編號（果糖）
  - `mirror: true`：左右翻，用於蔗糖的果糖
  - `flip: true`：翻面，也就是繞水平軸轉 180°，環 O 在前、上下互換。用於乳糖、纖維二糖右側的 glucose
- 連結的那一側在 subs 寫 `link`，例如 `links: [{from: [0, 1, down], to: [1, 4, down]}]`
- 上下照「圖上看到的」寫，驗證器直接用畫出來的 3D 幾何比對
- 驗證器會：
  - 確認每個環對到 SMILES 中不同的環
  - 確認醣苷 O 真的連著兩端的碳（例如 α1→4 寫成 α1→3 會報錯）
- **SMILES 不要手打**：用 PubChem 的單醣 SMILES，加上「以 anomeric OH 為根的醣基片段」（見 `tools/sugar_lib.py` 的 `glycosyl`）。拼接時 O 留在原位，立體標記不會跑掉
- GAG、chitin、peptidoglycan 的雙醣用 `tools/gen_misc.py` 的 `disaccharide()` 以 RDKit 組裝：把 donor 的 C1 改接到 acceptor 的 O（保持原本的立體組態），可加 O-／N-硫酸、C5 翻轉（GlcA → L-IdoA）、3-O-lactyl（NAM）
- β1→3 的連結在 Haworth 圖上一定會穿過環，所以重複雙醣一律畫它 1→4 的那一半（第二環翻面，同乳糖），另一個連結寫在 notes

**擁擠的取代基**：
- `up: ["^-O_3SO", 4]`：標籤寫成 [文字, 第幾個字元接鍵]，讓鍵接在最後的 O（寫在左邊的硫酸）
- `down_len: 0.6`、`up_len`：這個取代基的鍵長（列距）
- `down_dx: 0.35`、`up_dx`：把標籤往右推一點，例如 L-iduronate 環內的 COO⁻
- 環內放不下時可以省略 H（不寫那一側），驗證時視為 H

**椅式**：`chair: {subs: {1: {ax: H, eq: OH}, …}}`。
- 圖是真實 3D 椅式（⁴C₁：C1 右端朝下、C4 左端朝上、環 O 在右後）的正投影
- 驗證器用同一個 3D 模型，axial／equatorial 寫錯會報錯
- 環內側的取代基（C2、C3 朝上，C5 朝下）會畫得比較短。標籤請朝環內寫：右半寫 `HO`、左半寫 `OH`，避免與相鄰的字重疊

**環狀或小分子（鹼基、orotate、carbamoyl aspartate）**：用 `ring: {down: 4, right: 5, rotate: 0}`，用 RDKit 座標配上本資料庫的字型和線寬。
- `explicit: true`：環上每個原子都寫出（包含 C），環上的 N–H 畫成獨立的鍵。這是課本的鹼基卡片畫法
- `numbers: {原子索引: 環編號}`：灰色小編號，另外輸出成 `numbers` 圖層
- `mirror: true`：多一個左右相反的檢視 `ring_mirror`，拆件中多一張 `structure_mirror`
- `down`：讓這個原子朝正下方，例如嘌呤 N9、嘧啶 N1
- `right`：讓這個原子落在右側
- 索引就是 SMILES 中的原子順序
- 驗證只比對 SMILES 本身

**脂質 `lipid:`**（規格同時產生圖與 SMILES，驗證會比對雙鍵位置與 cis/trans）：
```yaml
lipid: {type: fatty_acid, chains: [{n: 18, cis: [9]}], numbers: true}          # Δ 編號從羧基 C1 數
lipid: {type: phospholipid, chains: [{n: 16}, {n: 18, cis: [9]}], head: choline}  # sn-1, sn-2
lipid: {type: tag, chains: [{n: 16}, {n: 18}, {n: 18, cis: [9]}]}
lipid: {type: sphingolipid, chains: [{n: 18}], head: phosphocholine}          # oh | galactose | glucose
```
- 頭基可以是 choline、ethanolamine、serine、inositol、glycerol、phosphate
- **脂肪酸**：橫向鋸齒排列，cis 雙鍵會彎折 60°；三個以上 cis 雙鍵時改成較小的角度，排成 U 形。可加 1、Δ 碳、末端碳與 ω 的灰色編號
- **甘油脂與鞘脂**：主鏈橫放，頭基朝上，尾巴往下垂；第一個 cis 彎折一律往右
- 不飽和鏈請放在右側（sn-2 或 sn-3），彎折才不會穿過其他鏈
- 脂肪酸的末端可以改：`end: coa`（–S–CoA，palmitoyl-CoA）、`end: carnitine`（acylcarnitine 酯）

**固醇類**：`ring: {orient: steroid, letters: true, methyls: false}`。
- `orient: steroid`：擺成課本姿勢（A、B 環在左下且同高，C、D 在右上）
- `letters: true`：標 A–D 環
- `methyls: false`：角甲基畫成線，不寫 CH₃
- 立體化學（楔形鍵）尚未畫出

**鹼基卡片**：`ring: {axis: [C4, C5], right: C6, h_down: N9, nh_bond: true, mirror: true}`（嘌呤）、`ring: {down: N1, right: C2, nh_bond: true, mirror: true}`（嘧啶）。
- 骨架式，只標 N
- 環上的 N–H 畫成鍵
- 自動多一張左右相反版
- 由 `tools/gen_nuc.py` 產生

**維生素與輔酶（大分子骨架式）**：用 `ring:`，整個分子由 RDKit 排座標。由 `tools/gen_vit.py` 產生。
- `methyls: false`：末端甲基畫成線，不寫 CH₃。類異戊二烯（retinoid、carotene、vitamin D、E、K）都用這個，否則標籤會互相重疊
- `coordgen: false`：改用 RDKit 的另一種排法。預設排法有字疊在一起時試試看（NADP⁺、methylene-THF）
- **反應部位上色**寫在 `decorate.highlight`，只用在單獨的分子圖與拆件，反應圖不會帶進去：
  ```yaml
  decorate:
    highlight:
      - {smarts: "NC(=O)c1ccc[n+]c1", color: magenta}       # 整個符合的片段（原子標籤＋兩端同色的鍵）
      - {smarts: "[CX4](C)(C)[OH:1]", color: magenta}       # 有原子編號時只塗編號的原子
      - {atoms: [14, 15], color: magenta}                   # 直接寫 SMILES 原子索引
  ```
  驗證會確認每個 SMARTS 在分子裡找得到。
- 反應圖要標「這一步變化的部分」時，寫在反應的 `decorate`：`decorate: {tpp: {highlight: [{smarts: "P(=O)([O-])OP(=O)([O-])[O-]"}]}}`
- 含核糖的輔酶（NAD、FAD、CoA、SAM）一律用 AMP 的糖模板組 SMILES（見 `gen_vit.py` 的 `nucleotide_sugar`、`ADO_TAIL`），立體化學才會和驗證過的 Haworth 一致
- `hide_charges: [N, Fe]`：圖上不寫這些元素的形式電荷（heme 的 Fe–N 配位在 SMILES 中寫成形式電荷，課本只寫 Fe）
- Fe、Co、Mg 這類兩個字母的元素，標籤以兩個字母一起置中
- **B₁₂** 用 `free` 示意圖（`schematic: true`，SMILES 只寫 `[Co]`）
- `seed: 14`：用打亂原子順序的副本重新排座標（只對 `coordgen: false` 有效），原子索引不受影響。`validate` 遇到**標籤重疊**會發出警告，照提示換 `coordgen`、`seed` 或 `rotate`
- 兩端都有字的鏈狀鍵（P–O、O–P、S–S）太短時，程式會自動把較小的一側整塊往外推，磷酸基的鍵不會被字吃掉

**醣鏈示意（SNFG 符號）**：在 `free` 的 atoms 中把文字寫成 `"@glc"`、`"@gal"`、`"@man"`、`"@glcnac"`、`"@fru"` 等，就會畫成 SNFG 符號（例如藍色圓 = glucose），例如 `glycogen`。

**示意圖（原子來源）**：用 `free` 加 `schematic: true`（不檢查 H 數），`notes: [[文字, x, y, 顏色]]` 是不算原子的說明文字，例如 `purine_atom_sources`。

**分子的預設標示 `decorate`**：寫在分子 YAML 裡，用於單獨的分子圖、拆件，以及寫了 `molecule_deco: default` 的路徑。反應圖不會使用它。胺基酸的 R 基底色就放在這裡：
```yaml
decorate:
  boxes: [{rows: [2, 3, 4, 5], name: R}, {parts: ["4.c"], color: yellow, name: S}]   # 每個 name 輸出成獨立圖層 box_<name>
  colors: {6: cyan}
```

**同一個糖的環狀與開鏈**：開鏈（Fischer）與環狀（Haworth）SMILES 不同，所以是兩個分子，例如 `glucose_6p`（開鏈，糖解用）與 `a_glucose_6p`、`b_glucose_6p`（環狀，肝醣與 PPP 用）。由 `tools/gen_carb.py` 產生：環狀磷酸糖一律由 `sugar_lib.py` 的 PubChem 單醣「接上」磷酸或 UDP，不重寫立體標記。

**Fischer 鏈的 Ⓟ**：列的中心文字寫 `CH_2O{P}` 時，驗證器要求那個碳真的接著磷酸酯，`CH_2OH` 則不能。這樣兩端都是 CH₂ 的分子（glycerol 3-phosphate）不會被倒過來配對而漏掉左右寫反。

**排法原則**：反應前後，同一個原子盡量留在同一個位置（老師的畫法）。例如 citrate 寫成 `COO⁻ / CH₂ / HO–C–COO⁻ / CH₂ / COO⁻` 的直線，轉成 isocitrate 時只看到 OH 往下移一列。一個分子在不同反應需要不同排法時，用 `chains` 另外寫，不要硬改預設排法。

**圖上文字只能用英文**：字型沒有中文字，中文只放在 `notes_zh`、`card_zh`。

**SMILES 的規則**（validate 會強制檢查 `chain` 和 `chains` 的每一種排法）：
1. SMILES 的書寫順序不限。程式會在 SMILES 裡找出一條碳鏈，其各列 H 數和圖相符；`COO⁻` 這類分支算取代基，不算列
2. 每一列在圖上畫出的 H 數必須和 SMILES 一致
3. 圖上畫的取代基（HO、COO⁻、S-CoA、O…）必須真的接在該碳上；中心是單純 C／CHₙ 的列，畫出的雙鍵也要和 SMILES 一致
4. Fischer 列（左右都有取代基、而且是立體中心）的 H/OH 左右位置，必須和 SMILES 的 `@`/`@@` 一致。報錯時，**以教科書的 Fischer 圖為準去修 SMILES**
5. 糖類和胺基酸請寫成生理 pH 下的離子型：磷酸 `OP(=O)([O-])[O-]`，羧酸 `C(=O)[O-]`，胺基 `[NH3+]`

## 輔因子 `data/cofactors.yaml`

反應圖上用彎箭頭畫的小分子（NAD⁺、CoA-SH、ATP…）。為了做平衡檢查，SMILES 只寫**會變化的部分**，其餘部分用虛擬原子代表：

| 虛擬原子 | 代表 |
|---|---|
| `[1*]` | CoA 的其餘部分 |
| `[2*]` | NAD |
| `[3*]` | FAD |
| `[4*]` | adenosine |
| `[5*]` | guanosine |
| `[6*]` | NADP |

**連到結構頁**：輔因子可以加 `molecule: nad`。網頁上點反應圖中的 NAD⁺ 就會跳到 NAD⁺ 的完整結構。驗證會確認輔因子的片段（虛擬原子當成任意原子）真的在那個分子裡。

分子裡也可以用同樣的虛擬原子，例如 acetyl-CoA 寫成 `O=C(S[1*])C`。

**ring 畫法裡的代號原子**：`[*]`（沒有編號）是接點，不會畫出；有編號的 `[16*]` 會依 `sdb/draw/skeletal.py` 的 `R_LABELS` 畫成文字（例如 LTC₄ 的 `–S–G`），讓大分子的一部分可以縮寫，同時保留原子平衡檢查。

**輔因子也可以是片段＋連結**：例如 `SAM: {label: SAM, smiles: "C[S+]([14*])[14*]", molecule: sam}`，反應圖上只寫「SAM」，點了跳到完整結構。BH4／BH2 用 `[13*]` 代表側鏈。

**有係數的輔因子**：總反應式要 8 CoA-SH、14 NADPH + H⁺ 時，用 `8CoA-SH`、`14NADPH+H+`、`7CO2` 這類名稱（SMILES 是同一片段重複 n 次，由 `tools/gen_lipmet.py` 寫入 `cofactors.yaml`）。`[9*]` 是「肝醣其餘部分」：`glycogen_n`、`glycogen_n1` 的 C1 寫 `OR`，glycogen synthase／phosphorylase 才能做原子平衡檢查。需要新的代號時，要同時加到 `cofactors.yaml` 的註解和 `sdb/draw/skeletal.py` 的 `R_LABELS`。

## 反應 `data/reactions/<id>.yaml`

```yaml
id: pdh
title: Pyruvate → Acetyl-CoA
substrates: [pyruvate]                 # 需要係數時寫 {id: acetaldehyde, n: 2}
products: [co2, acetyl_coa]
products_layout: stack                 # row（預設，左右排）| stack（上下疊）
cofactors_in: [CoA-SH, NAD+]           # 上方彎箭頭進入，左邊的放在較高層
cofactors_out: [NADH]                  # 上方彎箭頭離開
extra_in: [H+]                         # 選填：不畫在圖上、只參與平衡計算
extra_out: []
prosthetic: ["TPP,", "lipoate,", "FAD"]  # 箭頭上方中央，一個元素一行
prosthetic_font: serif                 # serif（預設）| sans（transketolase 的 TPP 用 sans）
enzyme: "pyruvate dehydrogenase\ncomplex (E_1 + E_2 + E_3)"   # 箭頭下方，青藍色
conditions: heat                       # 選填：箭頭下方的黑字，給有機反應用
reversible: false                      # true → 畫雙半箭頭
decorate:                              # 只在這張圖上生效的標示（上色，不用色塊）
  pyruvate: {colors: {0: magenta}}                            # 整列上色
  acetyl_coa: {colors: {0: {c: magenta, l: magenta, u: magenta}}}   # 只塗部分：c l r u d fork_l fork_r
  glucose_6p: {colors: {5: {c: [magenta, "{P}"]}}}            # 只塗一段字（這裡是 Ⓟ）
  xylulose_5p: {colors: {0: magenta, 1: magenta, 2: cyan}}   # 片段追蹤（transketolase）
  crotonaldehyde: {smarts: "[CH]=[CH]"}                       # 骨架式：把符合的鍵和原子上色
variants:                              # 選填：另一個版本，輸出成 <id>.isotope.svg，網頁可切換
  isotope:
    note_zh: "…版本說明，顯示在網頁圖下方…"
    decorate:                          # 疊加在上面的 decorate 之上
      citrate: {marks: {0: {c: circle}, 1: {c: "*"}, 2: {c: [tri, "2"], r: "3"}}}
      # 符號：* | circle（○）| tri（△）| "1"…"9"（①…⑨）；寫成 list 時，第一個放左上、第二個放右上
cofactors_inline: true                 # 選填：輔因子寫進方程式（+ H₂O、+ NAD⁺），不畫彎箭頭
balance_check: true                    # 預設 true；只在確實無法平衡時才關掉，並在 card_zh 說明原因
layout: {arrow_min: 6, arrow_at: 0.36, name_dy: 0, nudge: {co2: [0, -0.5]}}   # 手動微調，單位是 S（字級）；也可寫在 data/overrides/<id>.yaml
card_zh:
  location: 粒線體基質
  type: ...
  deltaG: ...
  cofactors: ...
  mechanism: ...
  regulation: ...
  deficiency: ...
  exam_points: [..., ...]
tags: [TCA]
status: pending
```

**配色的意思**（照 Mathews／Stryer：讓物質本身有顏色，**不用色塊**）：
- 洋紅：這一步中有變化的部分，例如 acetyl-CoA 的 acetyl 基、新接上的 Ⓟ、離開的 COO⁻
- 青藍：第二個追蹤片段（例如 transketolase 中供體剩下的部分），也是酶素文字的顏色
- 黑色：其餘全部，包括 OAA 這類沒有被強調的分子，以及輔因子
- 橙色：只用在同位素版的符號（* ○ △ ①②③）
- 鍵的顏色：兩端同色時，鍵跟著上色；否則畫黑色
- 色塊 `box` 程式還支援，但預設不要用

**原子追蹤**（Mathews 的畫法，用在氧化還原、加水、轉胺、尿素循環）：
- 會移動的 H 要明確畫出：寫成 `H–C–H`，不寫 CH₂。需要的話，另外寫一種 `chains.explicit` 排法
- 來源和去向用同一個顏色：第一組原子洋紅，第二組青藍
- 輔因子也可以上色，例如 `"NADH+H+": {colors: {0: {c: [[cyan, "H"], [cyan, "H^+"]]}}}`。一個部位可以有好幾段上色，依序往後尋找
- 例子：MDH 中 HO 的 H 與 C–H 的 H → NADH 的 H 與 H⁺；IDH 中青藍的 H → NADH + H⁺、洋紅的 COO⁻ → CO₂

**每一種反應都要成對上色**（來源和去向用同一個顏色，不能只塗一邊）：

| 反應類型 | 要上色的兩端 |
|---|---|
| 脫羧 | 離開的 COO⁻ ↔ CO₂（CO₂ 不論是產物分子還是輔因子，都要上色） |
| 脫氫（氧化） | 離開的 H（`l: [magenta, "H"]`、`r: [magenta, "H"]` 只塗 H）↔ NAD(P)H + H⁺ 的 H 與 H⁺ |
| 還原 | NAD(P)H + H⁺ ↔ 產物上新增的 H |
| 激酶／水解／磷解 | 新接上或離開的 Ⓟ ↔ ATP（生成 ATP 時）、Pi、PPi |
| 羧化 | HCO₃⁻ ↔ 新的 COO⁻ |
| 加水／脫水 | H₂O ↔ 新的 OH（或離開的 OH） |
| 胺基轉移／醯胺轉移 | ⁺H₃N ↔ ⁺H₃N；Gln（輔因子）↔ 新的 NH₂ |
| 甲基轉移 | SAM／N⁵-CH₃-THF 的 CH₃ ↔ 新的 CH₃ |
| 羥化（單氧酶） | O₂ ↔ 新的 OH（與 H₂O） |

- chain 畫法中掛在下面的環，只塗幾個原子：`{3: {ring: [{atoms: [5], color: magenta}]}}`（編號是 ring SMILES 的原子索引）
- Haworth 的鹼基：`base_colors: {6: magenta}`，或 `base_highlight: [{smarts: …}]`；相鄰兩個同色原子之間的鍵也會上色
- 第 4、7–10 輪的反應配色集中寫在 `tools/reaction_colors.py`，執行 `.venv/bin/python tools/reaction_colors.py data` 寫回 YAML；產生器也會讀這份檔案
- `validate` 會檢查每一個上色設定都有作用：對象在反應中、列號與部位存在、文字片段找得到、SMARTS 比對得到、原子編號在範圍內

**同位素版**：一般版保持乾淨，不加符號。同位素的標記一律放進 `variants.isotope`，不要混進一般版。

平衡檢查會把 substrates + cofactors_in + extra_in 和 products + cofactors_out + extra_out 做比較，原子數和電荷都要平衡。缺 H⁺ 時，用 `extra_in` / `extra_out` 補上，**不要關掉檢查**。

## 路徑 `data/pathways/<id>.yaml`

有四種版面（`layout`）：

| 版面 | 用途 | 節點位置怎麼決定 |
|---|---|---|
| `linear` | 一排、從左到右 | 依照箭頭順序自動排 |
| `free` | 直欄、分支，例如糖解 | 每個節點都寫 `pos: [x, y]`，單位是「格」（1 格＝最大節點的尺寸 × `cell`） |
| `cycle` | 環，例如 TCA | 沒寫 `pos` 的節點會依照順序、順時針排在半徑 `radius` 的圓上；寫了 `pos` 的節點（例如 acetyl-CoA）放在指定位置 |
| `plot` | 座標圖：酵素動力學、O₂ 結合曲線、滴定曲線、Ramachandran | 不用節點，改寫 `plot:`（見下方「作圖」） |

```yaml
id: tca_cycle
layout: cycle
radius: 2.15
start_angle: -110          # 第一個環上節點的角度；-90 是正上方
center_text: "Citric acid\ncycle"
cell: 1.0                  # 選填：格子縮放；名稱版總覽可以用 0.6
nodes:
  - {id: oaa, ref: oxaloacetate}                 # view 可以用 chain | skeletal | formula | name
  - {id: accoa, ref: acetyl_coa, pos: [-0.15, -3.25]}
  - {id: tca, text: "TCA cycle", kind: concept}  # 文字節點；kind 可以用 concept | enzyme | drug | plain
edges:
  - {from: oaa, to: cit, step: citrate_synthase, label: "citrate\nsynthase", side_text: "H_2O → CoA-SH"}
  - {from: accoa, to: "edge:oaa->cit", bend: 0}  # 箭頭指向另一條箭頭的中點（側邊匯入）
  - {from: cit, to: cis, reversible: true, label_side: left, bend: 0.12}
```

各欄位的意思：
- `label`：箭頭旁的青藍色酶素名稱
- `side_text`：放在箭頭另一側的黑色輔因子文字
- 這兩段文字會自動避開節點、箭頭和其他文字
- `label_side`：酶素名稱放在前進方向的哪一側。`cycle` 預設是 `right`（環內），其他版面預設是 `left`
- `bend`：箭頭的彎曲程度；正值會往前進方向的左側彎
- `from: "edge:a->b"`：箭頭從另一條箭頭的中點分出來，例如 fumarate 從 ASL 那一步離開；`to: "edge:a->b"` 則是匯入另一條箭頭
- 路徑層級的選項：
  - `cell_size: 7.2`：固定格子大小，單位是 S
  - `node_anchor: top`：`pos` 改成節點的頂端，例如讓一列胺基酸的 COO⁻ 對齊
  - `molecule_deco: default`：每個分子節點顯示自己的預設標示
- 文字節點的 `kind` 還有：
  - `heading`：靠左的分組標題
  - `name`：沒有分子檔的中間物名稱，例如 GAR、AICAR
- `head: false`：只畫線、不畫箭頭（家族樹）
- `resonance: true`：兩端都有箭頭的直線（↔），用於共振式
- `inhibit: true`：畫成紅色的抑制線（⊣），用於藥物和回饋抑制。例如 `{from: mtx, to: "edge:dhf->thf", inhibit: true}`
- `color: orange`：改變箭頭顏色。ETC 用橙色畫電子流，洋紅畫質子流
- 文字節點的 `kind` 另有 `protein`（淡藍方塊，膜上的複合體、轉運蛋白）、`carrier`（橙色圓角，Q、cyt c 等可移動載體）、`ion`（襯線字，H⁺、NADH + H⁺）；`w`、`h`（單位 S）指定方塊的最小寬高，讓複合體跨過膜
- `regions:` 畫在所有節點後面的背景，座標單位同 `pos`（格）：
  ```yaml
  regions:
    - {kind: membrane, y: 0, x0: -7, x1: 10, label: "Inner\nmitochondrial\nmembrane"}   # 水平磷脂雙層；直的寫 x, y0, y1
    - {kind: box, x0: -7, y0: -3.5, x1: -1, y1: 2.4, label: Liver}                        # 器官或胞器的虛線框
  ```
- **文字版（自動產生）**：凡是有分子節點的路徑，`build` 都會多做兩個版本
  - `text`：只有名稱和粗箭頭，不放酵素和輔因子
  - `text_full`：名稱，加上酵素和輔因子
  - 如果路徑本來就是名稱版，只會多出 `text`
  - 寫 `text: false` 可以關掉
  - `text_layout: {radius: 1.8, cell: 0.9, pos: {accoa: [0, -2.5]}}`：只覆寫文字版的排法
- `notes_zh`：路徑頁的中文筆記
- **反應也有文字版（自動產生）**：
  - `text`：只有名稱與箭頭
  - `text_full`：名稱，加上輔因子與酵素
  - 小分子（CO₂、NH₄⁺、HCO₃⁻）用化學式
  - 寫 `text: false` 可以關掉
- 有寫 `step` 的箭頭，在網頁上點了會跳到該反應

圈號（①…⑩）不在內建字型裡，編號請寫 `1.`、`2.`。

### 作圖 `layout: plot`（`sdb/draw/plot.py`，由 `tools/gen_plots.py` 產生）

```yaml
layout: plot
plot:                                  # 多張並排改寫 panels: [{title: …, <同樣的鍵>}, …]
  size: [15, 9.5]                      # 繪圖區大小（S）
  params: {Vmax: 1, Km: 2}             # 下面所有數字都可以寫成這些參數的算式，例如 "-1/Km"
  x: {min: 0, max: 20, label: "[S]", ticks: [0, 10, 20], marks: [{at: Km, label: "K_m"}]}
  y: {min: 0, max: 1.2, label: "v_0", label_pos: center}   # center = 沿軸旋轉的標題
  curves:
    - {id: mm, y: "mm(x, Vmax, Km)", color: red, dash: true, label: "…", label_at: 16, label_off: [0.3, -0.6]}
    - {id: gly, t: [0.6, 12.6], x: "…", y: "t"}            # 參數式曲線（滴定曲線）
  points: [{on: mm, x: Km, y: "Vmax/2", drop: true, label: "…", leader: true}]
  hlines: [{y: Vmax, label: "V_{max}"}]                     # 虛線輔助線；vlines 同理
  bands: [{x0: 20, x1: 40, label: Tissues, color: cyan}]    # 淡色直條
  areas: [{points: [[x, y], …], color: magenta, label: "α_R"}]   # 多邊形色塊（Ramachandran）
  arrows: [{from: [x, y], to: [x, y], color: red, label: "…"}]
  texts: [{at: [x, y], text: "…", align: left}]             # 多行文字依 align 對齊
  legend: {at: [x, y], items: [mm, {curve: inh, text: "+ Inhibitor"}]}
caption: {text: "…", color: red}       # 整張圖下方的一行重點（路徑層級）
```
- 內建函數：`mm(s, Vmax, Km)`、`hill(p, P50, n)`、`mwc(p, KR, KT, L)`、`logit(Y)`、`bell(x, mu, sigma)`、`exp`、`log10` 等
- **驗證**：`points` 若同時寫了 `x` 和 `y`，`validate` 會檢查這個點真的在曲線上（K_m 對應 V_max/2、−1/K_m 與 1/V_max 截距、P₅₀ 對應 Y = 0.5、pK 對應半當量）。只寫其中一個時，另一個由曲線算出
- 拆件：`axes`（座標軸＋刻度）、`curves`（曲線＋點）、`guides`（虛線、色帶、色塊）、`labels`、`arrows`，多張並排時另有 `panel_1`、`panel_2`

### 插畫 `scene`（`sdb/draw/scene.py`、元件在 `sdb/draw/glyphs.py`）

Alberts／Becker 風格的圖片版：膜、蛋白質、DNA 都畫成圖形，依座標擺放（單位 S）。有兩種用法：
- 整張圖都是插畫：寫 `layout: scene` 加 `scene: {items: [...]}`，例如 `replication_fork`
- 替既有的文字框路徑加一個「圖片版」：寫 `variants: {illus: {note_zh: 圖片版（插畫）, scene: {...}}}`，文字框版保留

```yaml
scene:
  caption: "…"                        # 選填：圖下方一行
  items:                              # 依序繪製，後面的蓋在前面上
    - {g: membrane, y: 7, x0: 0, x1: 40}
    - {id: pol, g: blob, at: [18, 6.3], w: 2.6, h: 2.1, shape: hand, hue: blue, accent: true, text: "Pol III",
       label: "DNA polymerase III", label_at: [12.5, 2.4], point: [17.4, 5.5]}
```
- 每個 item 都可以加 `label`：
  - 有 `label_at` 時，畫灰色引線指到元件（`point` 改指定引線終點）
  - 沒有 `label_at` 時，標籤放在元件正下方，不畫引線
  - `leader: false` 關閉引線
- `id` 是拆件名稱；沒寫時用 `<g><序號>`
- 配色（混合風格）：
  - 淡色平塗＋同色系深外框
  - `accent: true` 加放射漸層，用在這張圖的主角（聚合酶、被活化的酵素）
  - `hue` 可以用 blue、green、orange、purple、pink、yellow、teal、gray、red、brown，或 `#rrggbb`

| 元件 `g` | 參數 |
|---|---|
| `membrane` | `y`、`x0`、`x1`、`t`（厚度，預設 2.6） |
| `region` | 胞器或區域底色：`x0, y0, x1, y1`、`hue`、`dashed` |
| `blob` | 蛋白質：`at`、`w`、`h`、`shape`（oval、round、box、kidney、hand、bilobe）、`rot`、`text` |
| `gpcr` | 7 根跨膜螺旋：`at` 是膜中心 |
| `channel` | `open: true/false` |
| `transporter` | `open: out/in` |
| `pump` | 跨膜體＋胞內 ATP 結合葉 |
| `gprotein` | `nucleotide: GDP/GTP` |
| `ring` | `kind: clamp`（側視環）或 `hexamer`（helicase） |
| `strand` | 單股核酸：`points`、`kind`（parent、new、rna、primer）、`bases: left/right`（朝前進方向的左或右畫鹼基）、`arrow`、`ends: ["5'", "3'"]`、`seq` |
| `dsdna` | 雙股：`from`、`to`、`width`、`kinds`、`ends: true`、`seq` |
| `ion` | 小球加符號：`text`、`r` |
| `ligand` | `shape`：triangle、diamond、hexagon |
| `arrow` | `from`、`to`、`curve`、`style: fat`、`color`、`inhibit`、`dashed` |
| `text` | `text`、`size`、`font`、`color`、`align` |
| `helix` | 側視 B-DNA：`from`、`to`、`pitch`（每圈 S）、`amp`、`groove`（兩股相位差，決定大小溝） |
| `plasmid` | 環狀 DNA：`kind`（relaxed、supercoiled 加 `twists`、catenane、theta 加 `bubble: [deg0, deg1]`）、`r` |
| `nucleosome` | 組蛋白八聚體＋纏繞的 DNA：`at`、`r` |
| `chromosome` | 中期染色體：`at`、`h` |
| `mark` | DNA 上的小記號：`kind`（lesion 紅星、cut 切口三角加 `side: up/down`、x 紅叉、methyl 加 CH₃）、`at`、`color` |
| `trna` | tRNA 三葉草（acceptor 朝上）：`at`（中心）、`s`（縮放）、`anticodon`（文字）、`cca`（預設 true）、`leader` / `trailer`（前驅物兩端長度，灰色）、`intron`（反密碼子環的紅色 intron） ；`shape: icon` = 簡化 tRNA（腳站在 mRNA 上）：`anticodon`、`aa`（一顆胺基酸球）或 `peptide`（n 顆胜肽鏈）、`hue` |
| `ribosome` | 核糖體：`at` = P 位中心（mRNA 線上）、`w`、`sites`（E/P/A 槽，字母在小次單元）、`large` / `small`（次單元名稱）、`flip`（大次單元在下：ER 上的核糖體）、`hue` |

- **圖上只放短標籤**（使用者要求，2026-10-13）：名稱、幾個字、步驟編號；一行約 6 個字以內，不寫句子、不放整段說明。
  機轉、例子、考點一律寫進 `notes_zh`。產生器裡的 `note()` 段落會自動搬進 `notes_zh` 的【補充】；舊圖的長字由 `tools/figure_text.py` 縮短或搬走。
- **文字不要壓在圖形上**：
  - 說明一律用 `label`（加引線）或另外放 `text`
  - 蛋白質內部只放極短的代號（R、C、α/β/γ）
  - `validate` 會檢查標籤或文字是否蓋到蛋白質、小球、DNA 股，DNA 股沿著線段實際判斷
  - 方塊裡的字請用 `blob` 的 `text` 參數，不要另外疊一個 `text`
- **鹼基朝向**：`bases: left/right` 是相對於 `points` 的前進方向（y 軸朝下：往東走時 right = 朝下）。兩股要互相朝對方

- 鹼基顏色：A 綠、T/U 紅、G 黃、C 藍
- 股的顏色：親代深藍灰、新股橙、RNA 引子洋紅
- `validate` 會檢查元件名稱與參數，標籤重疊時發出警告
- 拆件：每個 item 一張。從第 18 輪起，新項目寫 `parts_png: false`，拆件只出 SVG 和 `_sheet.png`；插畫版本一律只出 SVG

## 拆件貼紙（自動產生）

`build` 會把每個項目拆成好幾張貼紙，放在 `out/parts/<kind>s/<id>/`，所有貼紙的比例都相同：
- **分子**：`structure`（不含名稱）、`name`、`row0…rowN`（chain 的每一列，不含列與列之間的直鍵）
- **反應**：
  - `mol_<id>`：分子，含這張圖上的標示
  - `name_<id>`：分子名稱
  - `highlight_<id>`：只有鮭魚色方框的標示層
  - `arrow`：箭頭加輔因子，不含酶素名稱
  - `enzyme`、`conditions`：箭頭下方的文字
- **路徑**：`path`（節點＋箭頭）、`enzymes`（酵素名稱層）、`cofactors`（輔因子層）。每個版本各一組，文字版也有
- **總表**：每個資料夾都有一張 `_sheet`，把所有部件排在同一張圖上

新增繪圖元件時，請把子圖放進 `Graphic.parts`，拆件就會自動輸出。

## 處理審核意見（使用者說「處理審核意見」時）

1. 執行 `.venv/bin/python -m sdb review`，列出所有尚未處理的意見，每則都有 `cid` 和對應的項目
2. 逐則修改對應的 YAML。只有共用的排版需要改時，才去動 `sdb/`；只想調整某一張圖時，用 `layout:` 微調
3. 執行 `validate` 和 `build`，再用眼睛確認修改後的圖
4. 每則意見處理完，執行 `.venv/bin/python -m sdb resolve <cid> "簡述做了什麼"`。這會把意見標成已處理，並把該項目的狀態改回 `pending`
5. 回報時列出每則意見怎麼處理；若有無法處理或需要使用者決定的，要明確說出來

## 化學正確性

- 酶素名稱、ΔG、輔因子、調控、缺乏症，都以 Lehninger 或 Harper 為準；有疑慮時寫在回報裡，不要自己猜
- 使用者提供參考圖（放在 `reference/`）時，以參考圖的畫法為準
