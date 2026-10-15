// structure_db local UI — plain JS, no build step.
const $ = (s, el = document) => el.querySelector(s);
const esc = (s) => String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const KIND_ZH = { molecule: "分子", reaction: "反應", pathway: "路徑" };
const STATUS_ZH = { draft: "草稿", pending: "待審", changes_requested: "待修改", approved: "已核准" };
const CARD = [["location", "位置"], ["type", "類型"], ["deltaG", "ΔG / 可逆性"], ["cofactors", "輔因子"],
  ["mechanism", "機制"], ["regulation", "調控"], ["deficiency", "缺乏症 / 臨床"], ["exam_points", "考點"]];

let IDX = null;
const state = { tab: "all", tag: null, sub: null, q: "" };
try { state.base = localStorage.getItem("sdb-base") || "rel"; } catch { state.base = "rel"; }

async function load() {
  IDX = await (await fetch("out/index.json", { cache: "no-store" })).json();
  IDX.all = [...IDX.molecules, ...IDX.reactions, ...IDX.pathways];
  IDX.byKey = Object.fromEntries(IDX.all.map((x) => [`${x.kind}:${x.id}`, x]));
  const tt = IDX.tag_tree || {};
  IDX.tree = tt.tree || {};
  IDX.tagZh = tt.labels || {};
  IDX.parentOf = Object.fromEntries(Object.entries(IDX.tree).flatMap(([p, cs]) => cs.map((c) => [c, p])));
  if (IDX.static) document.body.classList.add("static");
}

const tagLabel = (t) => IDX.tagZh[t] || t;
// "胺基酸／鹼性, TCA" — a child tag replaces its parent
function tagText(x) {
  const kids = x.tags.filter((t) => IDX.parentOf[t]);
  const parents = new Set(kids.map((t) => IDX.parentOf[t]));
  return x.tags.filter((t) => !parents.has(t))
    .map((t) => IDX.parentOf[t] ? `${tagLabel(IDX.parentOf[t])}／${tagLabel(t)}` : tagLabel(t)).join(", ");
}

const title = (x) => x.title || x.name || x.id;
const fileOf = (x, view) => (x.files.find((f) => f.view === view) || x.files[0] || {}).file;
const svgUrl = (f) => `out/svg/${f}.svg?t=${IDX.generated}`;
const pngUrl = (f) => `out/png/${f}.png?t=${IDX.generated}`;

function toast(msg) {
  const t = $("#toast"); t.textContent = msg; t.classList.add("show");
  clearTimeout(t._h); t._h = setTimeout(() => t.classList.remove("show"), 1600);
}

function matches(x) {
  if (state.tag && !x.tags.includes(state.tag)) return false;
  if (state.sub && !x.tags.includes(state.sub)) return false;
  if (state.tab === "review") return x.status === "pending" || x.status === "changes_requested";
  if (state.tab !== "all" && x.kind !== state.tab) return false;
  if (!state.q) return true;
  const hay = [x.id, title(x), ...(x.aliases || []), x.enzyme, x.notes_zh, JSON.stringify(x.card_zh || {}), ...x.tags]
    .join(" ").toLowerCase();
  return state.q.toLowerCase().split(/\s+/).every((w) => hay.includes(w));
}

function renderChrome() {
  document.querySelectorAll("#tabs button").forEach((b) => b.classList.toggle("on", b.dataset.tab === state.tab));
  const used = new Set(IDX.all.flatMap((x) => x.tags));
  const chip = (t, attr, on) => `<span class="tag ${on ? "on" : ""} ${IDX.tree[t] ? "parent" : ""}" ${attr}="${esc(t)}" title="${esc(t)}">${esc(tagLabel(t))}</span>`;
  const top = [...used].filter((t) => !IDX.parentOf[t]).sort((a, b) => tagLabel(a).localeCompare(tagLabel(b)));
  const kids = (IDX.tree[state.tag] || []).filter((t) => used.has(t));
  $("#tags").innerHTML = `<div class="tagrow">${top.map((t) => chip(t, "data-tag", state.tag === t)).join("")}</div>` +
    (kids.length ? `<div class="tagrow sub">${kids.map((t) => chip(t, "data-sub", state.sub === t)).join("")}</div>` : "");
  const n = IDX.all.filter((x) => x.status === "pending" || x.status === "changes_requested").length;
  $("#qcount").textContent = n || "";
  $("#base").value = state.base;
}

function renderList(currentKey) {
  const items = IDX.all.filter(matches);
  $("#list").innerHTML = items.map((x) => {
    const key = `${x.kind}:${x.id}`;
    return `<a class="item ${key === currentKey ? "on" : ""}" href="#/${x.kind}/${x.id}">
      <div>${esc(title(x))}<span class="st ${x.status}">${STATUS_ZH[x.status]}</span></div>
      <div class="k">${KIND_ZH[x.kind]} · ${esc(tagText(x))}</div></a>`;
  }).join("") || `<div class="sub">沒有符合的項目</div>`;
}

async function inlineSvg(el, file, zoom = 1.35) {
  const txt = await (await fetch(svgUrl(file))).text();
  el.innerHTML = txt;
  el.querySelectorAll("svg > title").forEach((t) => t.remove());   // no native tooltip over the figure
  // the root <svg> itself carries data-ref of the page we are on; don't make it a link
  const root = el.querySelector("svg");
  const here = root && root.dataset.ref;
  if (root) {
    root.removeAttribute("data-ref");
    root.style.cursor = "default";
    for (const a of ["width", "height"]) root.setAttribute(a, Math.round(+root.getAttribute(a) * zoom));
  }
  el.querySelectorAll("[data-ref]").forEach((g) => {
    let [kind, id] = g.dataset.ref.split(/:(.*)/s);
    // a cofactor (NAD⁺, CoA-SH …) links to its structure page when cofactors.yaml names one
    const cof = kind === "cofactor";
    if (cof) { const c = (IDX.cofactors || {})[id]; if (!c || !c.molecule) return; kind = "molecule"; id = c.molecule; }
    const x = IDX.byKey[`${kind}:${id}`];
    if (!x || g.dataset.ref === here || `${kind}:${id}` === here) return;
    g.classList.add("link");
    g.addEventListener("click", (ev) => { ev.stopPropagation(); hideTip(); location.hash = `#/${kind}/${id}`; });
    g.addEventListener("mousemove", (ev) => showTip(ev, `→ ${KIND_ZH[kind]}：${title(x)}`));
    g.addEventListener("mouseleave", hideTip);
    // molecules get a padded hit area (the whole block is clickable, not only lines and glyphs);
    // arrows already carry their own hit rectangles and their bbox would cover neighbouring nodes
    if (!g.classList.contains("sdb-mol") && !cof) return;
    const b = g.getBBox();
    const pad = Math.max(4, 0.08 * Math.min(b.width, b.height) + 4);
    const r = document.createElementNS("http://www.w3.org/2000/svg", "rect");
    Object.entries({ x: b.x - pad, y: b.y - pad, width: b.width + 2 * pad, height: b.height + 2 * pad, rx: 6, class: "hit" })
      .forEach(([k, v]) => r.setAttribute(k, v));
    g.prepend(r);
  });
}

function showTip(ev, text) {
  let t = $("#hovertip");
  if (!t) { t = document.createElement("div"); t.id = "hovertip"; document.body.appendChild(t); }
  t.textContent = text;
  t.style.left = `${ev.clientX + 14}px`; t.style.top = `${ev.clientY + 16}px`;
  t.classList.add("show");
}
function hideTip() { const t = $("#hovertip"); if (t) t.classList.remove("show"); }

function mdRef(x, file) {
  const base = state.base === "abs" ? IDX.abs_out : IDX.notes_base;
  return `<img src="${base}/svg/${file}.svg" alt="${title(x)}">`;
}

function toolbar(x, file) {
  if (IDX.static) return `<div class="row">
    <a class="btn" href="${svgUrl(file)}" download="${file.split("/").pop()}.svg">下載 SVG</a></div>`;
  return `<div class="row">
    <button class="btn" data-copy="${esc(mdRef(x, file))}">複製 Markdown 引用</button>
    <a class="btn" href="${svgUrl(file)}" download="${file.split("/").pop()}.svg">下載 SVG</a>
    <a class="btn" href="${pngUrl(file)}" download="${file.split("/").pop()}.png">下載 PNG（3x 透明）</a>
    <span class="sub">PNG 路徑：<code>out/png/${esc(file)}.png</code></span></div>`;
}

function reviewBox(x) {
  if (IDX.static) return "";
  const comments = (x.comments || []).map((c) => `<div class="comment"><b>${esc(c.cid)}</b> · ${esc(c.ts)}<br>${esc(c.text)}</div>`).join("");
  return `<div class="review" data-key="${x.kind}:${x.id}">
    <div class="row"><b>審核</b><span class="st ${x.status}">${STATUS_ZH[x.status]}</span></div>
    ${comments ? `<div class="sub" style="margin-top:6px">未處理的意見（對 AI 說「處理審核意見」即可）：</div>${comments}` : ""}
    <textarea placeholder="寫下要修改的地方，例如：「pyruvate 的 C=O 應該畫在左邊」"></textarea>
    <div class="row">
      <button class="btn primary" data-act="approve">✓ 核准</button>
      <button class="btn warn" data-act="comment">送出意見</button>
      ${x.status === "approved" ? `<button class="btn" data-act="reopen">改回待審</button>` : ""}
    </div></div>`;
}

function thumbs(keys) {
  return `<div class="grid">${keys.map((k) => IDX.byKey[k]).filter(Boolean).map((x) =>
    `<a class="thumb" href="#/${x.kind}/${x.id}"><span class="im"><img src="${svgUrl(fileOf(x))}" alt=""></span><div>${esc(title(x))}</div></a>`).join("")}</div>`;
}

const PART_FINE = { arrow_line: "箭頭（不含輔因子）", cofactors: "輔因子（含彎箭頭）", prosthetic: "輔基／條件（箭頭上方）",
  regions: "膜／方框", molecules: "分子節點", proteins: "蛋白質／載體", labels: "文字節點", drugs: "藥物／疾病",
  arrows: "一般箭頭", electrons: "電子流（橙）", protons: "質子流（洋紅）", inhibitions: "抑制線",
  axes: "座標軸＋刻度", curves: "曲線＋標點", guides: "虛線／區域" };
const PART_ZH = (n) => PART_FINE[n] ? PART_FINE[n] : n === "structure" ? "結構式" : n === "structure_mirror" ? "結構式（左右相反）"
  : n === "numbers" ? "編號層" : n === "base" ? "鹼基" : n === "path" ? "路徑（節點＋箭頭）"
  : n === "enzymes" ? "酵素名稱層" : n === "cofactors" ? "輔因子層" : n === "name" ? "名稱" : n === "arrow" ? "箭頭＋輔因子（整組）"
  : n === "enzyme" ? "酶素名稱" : n === "conditions" ? "反應條件" : n.startsWith("row") ? `第 ${+n.slice(3) + 1} 列`
  : n.startsWith("panel_") ? `第 ${n.slice(6)} 張圖` : n.startsWith("mol_") ? `分子：${n.slice(4)}` : n.startsWith("name_") ? `名稱：${n.slice(5)}`
  : n.startsWith("highlight") ? `標示層${n.includes("_") ? "：" + n.split("_").slice(1).join("_") : ""}`
  : n.startsWith("box_") ? `底色層：${{ R: "R 基", S: "S（硫）" }[n.slice(4)] || n.slice(4)}` : n;

const VIEW_ZH = { default: "一般版", isotope: "同位素版", text: "文字版", text_full: "文字版＋酵素", illus: "圖片版",
  haworth: "Haworth", chair: "椅式", ring: "結構式", ring_mirror: "結構式（左右相反）", free: "結構式", chain: "縮寫 / Fischer", skeletal: "骨架式", formula: "分子式" };

function viewButtons(x, view) {
  if (x.files.length < 2) return "";
  return `<div class="row views" style="margin-top:12px">${x.files.map((f) =>
    `<button class="btn ${f.view === view ? "on" : ""}" data-view="${f.view}">${VIEW_ZH[f.view] || f.view}</button>`).join("")}</div>`;
}

function partsSection(x, view = "default") {
  const list = (x.parts || {})[view] || [];
  if (!list.length || IDX.static) return "";
  const dir = list[0].file.split("/").slice(0, -1).join("/");
  const v = `?t=${IDX.generated}`;
  return `<h2>拆件貼紙（${list.length}）</h2>
    <div class="sub">每個部件各一張 PNG（3x、透明背景、所有圖同一比例）。iPad 路徑：<code>out/parts/${esc(dir)}/</code>
      · <a href="out/parts/${dir}/_sheet.png${v}" download="${esc(x.id)}_sheet.png">下載總表 PNG</a></div>
    <div class="grid parts">${list.map((p) => `<a class="thumb" href="out/parts/${p.file}.png${v}" download="${esc(x.id)}_${esc(p.name)}.png">
      <span class="im"><img src="out/parts/${p.file}.svg${v}" alt=""></span><div>${esc(PART_ZH(p.name))}</div></a>`).join("")}</div>`;
}

function cardTable(card) {
  const rows = CARD.filter(([k]) => card && card[k]).map(([k, label]) => {
    const v = card[k];
    return `<tr><th>${label}</th><td>${Array.isArray(v) ? `<ul>${v.map((s) => `<li>${esc(s)}</li>`).join("")}</ul>` : esc(v)}</td></tr>`;
  }).join("");
  return rows ? `<table class="card">${rows}</table>` : "";
}

async function showMolecule(m, view) {
  view = view || m.default_view;
  const file = fileOf(m, view);
  $("#detail").innerHTML = `
    <h1>${esc(m.name)}</h1>
    <div class="sub">${esc(m.aliases.join("、"))}</div>
    <div class="sub">分子式 <code>${esc(m.formula || "")}</code> · MW ${esc(m.mw ?? "")} · SMILES <code>${esc(m.smiles)}</code></div>
    ${viewButtons(m, view)}
    <div class="fig" id="fig"></div>
    ${toolbar(m, file)}
    ${partsSection(m)}
    ${m.notes_zh ? `<h2>筆記</h2><p>${esc(m.notes_zh)}</p>` : ""}
    <h2>參與的反應（${m.reactions.length}）</h2>${m.reactions.length ? thumbs(m.reactions.map((r) => `reaction:${r}`)) : `<div class="sub">尚無</div>`}
    ${m.pathways.length ? `<h2>出現在路徑</h2>${thumbs(m.pathways.map((p) => `pathway:${p}`))}` : ""}
    ${reviewBox(m)}`;
  await inlineSvg($("#fig"), file);
  document.querySelectorAll("[data-view]").forEach((b) => b.onclick = () => location.hash = `#/molecule/${m.id}/${b.dataset.view}`);
}

async function showReaction(r, view = "default") {
  const file = fileOf(r, view);
  const mols = [...new Set([...r.substrates, ...r.products].map((p) => `molecule:${p.id}`))];
  const chip = (p) => { const m = IDX.byKey[`molecule:${p.id}`]; return m ? `<a class="chip" href="#/molecule/${m.id}">${esc(title(m))}</a>` : ""; };
  const eq = `<div class="row eq"><span class="sub">跳到：</span>${r.substrates.map(chip).join('<span class="op">+</span>')}
    <span class="op">${r.reversible ? "⇌" : "→"}</span>${r.products.map(chip).join('<span class="op">+</span>')}</div>`;
  $("#detail").innerHTML = `
    <h1>${esc(r.title)}</h1>
    <div class="sub">${esc(r.enzyme.replace(/_/g, ""))}${r.reversible ? " · 可逆" : ""}</div>
    ${viewButtons(r, view)}
    <div class="fig big" id="fig"></div>
    ${view !== "default" && r.variant_notes[view] ? `<p class="note">${esc(r.variant_notes[view])}</p>` : ""}
    ${eq}
    ${toolbar(r, file)}
    ${cardTable(r.card_zh)}
    ${partsSection(r, view)}
    <h2>相關分子</h2>${thumbs(mols)}
    ${r.pathways.length ? `<h2>所屬路徑</h2>${thumbs(r.pathways.map((p) => `pathway:${p}`))}` : ""}
    ${reviewBox(r)}`;
  await inlineSvg($("#fig"), file);
  document.querySelectorAll("[data-view]").forEach((b) => b.onclick = () => location.hash = `#/reaction/${r.id}/${b.dataset.view}`);
}

async function showPathway(p, view = "default") {
  const file = fileOf(p, view);
  $("#detail").innerHTML = `
    <h1>${esc(p.title)}</h1><div class="sub">版面：${esc(p.layout)} · 點箭頭看該步反應，點分子看分子頁</div>
    ${viewButtons(p, view)}
    <div class="fig big" id="fig"></div>
    ${view !== "default" && p.variant_notes[view] ? `<p class="note">${esc(p.variant_notes[view])}</p>` : ""}
    ${toolbar(p, file)}
    ${p.notes_zh ? `<h2>筆記</h2><p>${esc(p.notes_zh)}</p>` : ""}
    ${partsSection(p, view)}
    <h2>步驟</h2>${thumbs(p.steps.map((s) => `reaction:${s}`))}
    ${reviewBox(p)}`;
  await inlineSvg($("#fig"), file);
  document.querySelectorAll("[data-view]").forEach((b) => b.onclick = () => location.hash = `#/pathway/${p.id}/${b.dataset.view}`);
}

async function showQueue() {
  const items = IDX.all.filter(matches);
  $("#detail").innerHTML = `<h1>審核佇列</h1>
    <div class="sub">共 ${items.length} 項待審。核准 → 狀態變「已核准」；送出意見 → 存到 review/comments.jsonl，之後對 AI 說「處理審核意見」。</div>
    <div id="q-items" style="margin-top:14px"></div>`;
  const box = $("#q-items");
  for (const x of items) {
    const div = document.createElement("div");
    div.className = "qcard";
    div.innerHTML = `<div class="row"><b><a href="#/${x.kind}/${x.id}">${esc(title(x))}</a></b><span class="sub">${KIND_ZH[x.kind]}</span></div>
      <div class="fig"><img src="${svgUrl(fileOf(x))}" style="max-width:100%;max-height:340px"></div>
      ${x.kind === "reaction" ? cardTable(x.card_zh) : x.notes_zh ? `<p class="sub">${esc(x.notes_zh)}</p>` : ""}
      ${reviewBox(x)}`;
    box.appendChild(div);
  }
}

async function route() {
  if (!IDX) await load();
  const [, kind, id, view] = location.hash.split("/");
  if (kind === "review") state.tab = "review";
  renderChrome();
  const x = IDX.byKey[`${kind}:${id}`];
  renderList(x ? `${kind}:${id}` : null);
  if (state.tab === "review" && !x) return showQueue();
  if (!x) {
    $("#detail").innerHTML = `<h1>結構式資料庫</h1><p class="sub">${IDX.molecules.length} 個分子 · ${IDX.reactions.length} 個反應 · ${IDX.pathways.length} 條路徑 · 建置於 ${esc(IDX.generated)}</p>
      <h2>路徑</h2>${thumbs(IDX.pathways.map((p) => `pathway:${p.id}`))}<h2>反應</h2>${thumbs(IDX.reactions.map((r) => `reaction:${r.id}`))}`;
    return;
  }
  if (x.kind === "molecule") return showMolecule(x, view);
  if (x.kind === "reaction") return showReaction(x, view);
  return showPathway(x, view);
}

document.addEventListener("click", async (ev) => {
  const t = ev.target.closest("[data-tab],[data-tag],[data-sub],[data-copy],[data-act]");
  if (!t) return;
  if (t.dataset.tab) {
    state.tab = t.dataset.tab;
    if (state.tab === "review") location.hash = "#/review";
    else if (location.hash === "#/review") location.hash = "#/";
    else route();
  } else if (t.dataset.tag) {
    state.tag = state.tag === t.dataset.tag ? null : t.dataset.tag; state.sub = null; route();
  } else if (t.dataset.sub) {
    state.sub = state.sub === t.dataset.sub ? null : t.dataset.sub; route();
  } else if (t.dataset.copy) {
    await navigator.clipboard.writeText(t.dataset.copy); toast("已複製 Markdown 引用");
  } else if (t.dataset.act) {
    const box = t.closest(".review");
    const [kind, id] = box.dataset.key.split(":");
    const text = $("textarea", box).value;
    const res = await fetch("/api/review", { method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ kind, id, action: t.dataset.act, text }) });
    const j = await res.json();
    if (!res.ok) return toast(j.error || "失敗");
    toast({ approve: "已核准", comment: "意見已儲存", reopen: "已改回待審" }[t.dataset.act]);
    IDX = null; route();
  }
});
$("#q").addEventListener("input", (e) => { state.q = e.target.value.trim(); renderList(); });
$("#base").addEventListener("change", (e) => {
  state.base = e.target.value;
  try { localStorage.setItem("sdb-base", state.base); } catch {}
  route();
});
window.addEventListener("hashchange", route);
route();
