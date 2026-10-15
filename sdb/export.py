"""Build everything under out/: SVG, PNG (3x, transparent), tag galleries (.md), index.json."""
import json
import os
import shutil
from pathlib import Path
from datetime import datetime

import resvg_py
import yaml

from . import review, style
from .compose.molecule import name_block, show_name, structure, views_of, with_name
from .compose.pathway import compose_pathway
from .compose.reaction import compose_reaction
from .graphic import Graphic, document
from .style import OUT, ROOT
from .validate import validate

NOTES_BASE = "../../structure_db/out"     # relative to typora_file/生化/


def svg_name(kind, id_, view=None, default=None):
    stem = id_ if view in (None, default) else f"{id_}.{view}"
    return f"{kind}s/{stem}"


def _write(rel, svg):
    (OUT / "svg" / f"{rel}.svg").write_text(svg, encoding="utf-8")
    png = resvg_py.svg_to_bytes(svg_string=svg, zoom=style.PNG_ZOOM, skip_system_fonts=True)
    (OUT / "png" / f"{rel}.png").write_bytes(bytes(png))


# ---------------------------------------------------------------- stickers split into parts

def molecule_parts(m):
    """Level 1–3 split: structure (text colours, no background) / name / each chain row /
    each background layer (box_R, box_S …) on its own."""
    view = m["default_view"]
    deco = m.get("decorate") or {}
    g = structure(m, view, {k: v for k, v in deco.items() if k not in ("boxes", "box")})
    full = structure(m, view, deco) if deco else g
    parts = {"structure": Graphic(g.body, g.x0, g.y0, g.x1, g.y1)}
    for k, layer in full.parts.items():
        if k.startswith("box_") or k == "highlight":
            parts[k] = layer
    if show_name(m, view):
        tb = name_block(m)
        parts["name"] = Graphic(tb.svg(0, 0), -tb.width / 2, tb.top, tb.width / 2, tb.bottom)
    if "ring_mirror" in views_of(m):                 # left/right twin of a base as its own sticker
        mg = structure(m, "ring_mirror", None)
        parts["structure_mirror"] = Graphic(mg.body, mg.x0, mg.y0, mg.x1, mg.y1)
    if "numbers" in g.parts:
        parts["numbers"] = g.parts["numbers"]
    rows = sorted((k for k in g.parts if k.startswith("row")), key=lambda k: int(k[3:]))
    if len(rows) > 1:
        for k in rows:
            parts[k] = g.parts[k]
    return parts


def _part_order(name):
    for i, prefix in enumerate(("structure", "name", "row", "mol_", "name_", "arrow", "cofactors", "prosthetic",
                                "enzyme", "conditions", "highlight", "box_", "path", "regions", "molecules",
                                "proteins", "labels", "drugs", "electrons", "protons", "inhibitions",
                                "panel_", "axes", "curves", "guides", "arrows")):
        if name.startswith(prefix):
            return i
    return 99


def sheet(parts, max_w=30 * style.S, gap=1.2 * style.S):
    """All parts of one item laid out left→right, wrapping lines: a one-page sticker overview."""
    body, x, y, line_h, w_used = [], 0.0, 0.0, 0.0, 0.0
    for name in sorted(parts, key=lambda n: (_part_order(n), n)):
        g = parts[name]
        if x > 0 and x + g.w > max_w:
            x, y, line_h = 0.0, y + line_h + gap, 0.0
        body.append(g.moved(x - g.x0, y - g.y0).body)
        x += g.w + gap
        line_h = max(line_h, g.h)
        w_used = max(w_used, x - gap)
    return Graphic("".join(body), 0, 0, w_used, y + line_h)


def write_parts(kind, id_, title, parts, variant=None, png=True):
    """Write each part as its own SVG + PNG (SVG only when png=False, i.e. `parts_png: false`; the
    _sheet always gets a PNG); variants go to a sub-folder (…/<id>/isotope/)."""
    if variant:
        id_ = f"{id_}/{variant}"
    folder = OUT / "parts" / f"{kind}s" / id_
    folder.mkdir(parents=True, exist_ok=True)
    out = []
    pad = 0.3 * style.S
    for name, g in sorted(parts.items(), key=lambda kv: (_part_order(kv[0]), kv[0])):
        svg = document(g, f"{title} — {name}", pad=pad)
        (folder / f"{name}.svg").write_text(svg, encoding="utf-8")
        if png:
            (folder / f"{name}.png").write_bytes(bytes(resvg_py.svg_to_bytes(svg_string=svg, zoom=style.PNG_ZOOM, skip_system_fonts=True)))
        out.append({"name": name, "file": f"{kind}s/{id_}/{name}"})
    svg = document(sheet(parts), f"{title} — sticker sheet", pad=pad)
    (folder / "_sheet.svg").write_text(svg, encoding="utf-8")
    (folder / "_sheet.png").write_bytes(bytes(resvg_py.svg_to_bytes(svg_string=svg, zoom=style.PNG_ZOOM, skip_system_fonts=True)))
    return out


def render_all(db, problems):
    shutil.rmtree(OUT / "parts", ignore_errors=True)
    for sub in ("svg", "png"):
        shutil.rmtree(OUT / sub, ignore_errors=True)
        for k in ("molecules", "reactions", "pathways"):
            (OUT / sub / k).mkdir(parents=True, exist_ok=True)

    files = {}
    for m in db.molecules.values():
        files[("molecule", m["id"])] = []
        for view in views_of(m):
            try:
                g = structure(m, view, "default")
                if show_name(m, view):
                    g = with_name(g, m)
                rel = svg_name("molecule", m["id"], view, m["default_view"])
                _write(rel, document(g, m["name"], data={"ref": f"molecule:{m['id']}", "view": view}))
                files[("molecule", m["id"])].append({"view": view, "file": rel})
            except Exception as e:
                problems.append(f"molecule {m['id']} ({view}) 繪圖失敗：{e}")
        try:
            files[("molecule-parts", m["id"])] = {"default": write_parts("molecule", m["id"], m["name"], molecule_parts(m),
                                                                             png=m.get("parts_png", True))}
        except Exception as e:
            problems.append(f"molecule {m['id']} 拆件失敗：{e}")
    for kind, table, fn in (("reaction", db.reactions, compose_reaction), ("pathway", db.pathways, compose_pathway)):
        for x in table.values():
            try:
                title = x.get("title", x["id"])
                files[(kind, x["id"])], files[(f"{kind}-parts", x["id"])] = [], {}
                for variant in [None] + list(x.get("variants") or {}):
                    view = variant or "default"
                    rel = svg_name(kind, x["id"]) + (f".{variant}" if variant else "")
                    g = fn(db, x, variant)
                    _write(rel, document(g, title, data={"ref": f"{kind}:{x['id']}", "view": view}))
                    files[(kind, x["id"])].append({"view": view, "file": rel})
                    if g.parts:
                        files[(f"{kind}-parts", x["id"])][view] = write_parts(kind, x["id"], title, g.parts, variant,
                                                                                          png=x.get("parts_png", True) and not (variant and x["variants"][variant].get("scene")))
            except Exception as e:
                problems.append(f"{kind} {x['id']} 繪圖失敗：{e}")
    return files


def tag_tree():
    """data/tags.yaml → {"tree": {parent: [children]}, "labels": {tag: 中文}} for the web sidebar."""
    f = ROOT / "data" / "tags.yaml"
    d = yaml.safe_load(f.read_text(encoding="utf-8")) if f.exists() else {}
    return {"tree": d.get("tree") or {}, "labels": d.get("labels") or {}}


def build_index(db, files, info):
    comments = review.load()
    open_c = {}
    for c in comments:
        if not c.get("resolved") and c["action"] == "comment":
            open_c.setdefault(f"{c['kind']}:{c['id']}", []).append(c)

    def common(kind, x):
        key = f"{kind}:{x['id']}"
        return {"id": x["id"], "kind": kind, "tags": x["tags"], "status": x["status"],
                "files": files.get((kind, x["id"]), []), "parts": files.get((f"{kind}-parts", x["id"]), {}),
                "variant_notes": {k: v.get("note_zh", "") for k, v in (x.get("variants") or {}).items()},
                "comments": open_c.get(key, [])}

    idx = {"generated": datetime.now().isoformat(timespec="seconds"), "notes_base": NOTES_BASE,
           "abs_out": str(OUT),
           "molecules": [], "reactions": [], "pathways": [], "cofactors": db.cofactors, "tag_tree": tag_tree()}
    for m in db.molecules.values():
        d = common("molecule", m)
        d.update(name=m["name"], aliases=m["aliases"], smiles=m["smiles"], notes_zh=m.get("notes_zh", ""),
                 default_view=m["default_view"], reactions=db.reactions_of(m["id"]),
                 pathways=[p["id"] for p in db.pathways.values() if any(n.get("ref") == m["id"] for n in p.get("nodes", []))],
                 **info.get(("molecule", m["id"]), {}))
        idx["molecules"].append(d)
    for r in db.reactions.values():
        d = common("reaction", r)
        d.update(title=r.get("title", r["id"]), enzyme=r.get("enzyme", ""), reversible=r["reversible"],
                 substrates=r["substrates"], products=r["products"],
                 cofactors_in=r["cofactors_in"], cofactors_out=r["cofactors_out"],
                 prosthetic=r["prosthetic"], card_zh=r.get("card_zh", {}),
                 pathways=[p["id"] for p in db.pathways.values() if any(e.get("step") == r["id"] for e in p.get("edges", []))])
        idx["reactions"].append(d)
    for p in db.pathways.values():
        d = common("pathway", p)
        d.update(title=p.get("title", p["id"]), layout=p["layout"], notes_zh=p.get("notes_zh", ""),
                 steps=[e["step"] for e in p.get("edges", []) if e.get("step")])
        idx["pathways"].append(d)
    (OUT / "index.json").write_text(json.dumps(idx, ensure_ascii=False, indent=1), encoding="utf-8")
    return idx


CARD_LABELS = [("location", "位置"), ("type", "類型"), ("deltaG", "ΔG / 可逆性"), ("cofactors", "輔因子"),
               ("mechanism", "機制"), ("regulation", "調控"), ("deficiency", "缺乏症 / 臨床"),
               ("exam_points", "考點")]


def card_md(card):
    lines = []
    for k, label in CARD_LABELS:
        v = card.get(k)
        if not v:
            continue
        if isinstance(v, list):
            lines.append(f"- **{label}**")
            lines += [f"    - {x}" for x in v]
        else:
            lines.append(f"- **{label}**：{v}")
    return "\n".join(lines)


def build_galleries(db, files):
    md_dir = OUT / "md"
    shutil.rmtree(md_dir, ignore_errors=True)
    md_dir.mkdir(parents=True)
    tags = sorted({t for tbl in (db.molecules, db.reactions, db.pathways) for x in tbl.values() for t in x["tags"]})
    for tag in tags:
        out = [f"# {tag}", "", f"> 自動產生（`python -m sdb build`），請勿手動編輯。", ""]
        for kind, tbl, head in (("pathway", db.pathways, "路徑"), ("reaction", db.reactions, "反應"),
                                ("molecule", db.molecules, "分子")):
            items = [x for x in tbl.values() if tag in x["tags"] and files.get((kind, x["id"]))]
            if not items:
                continue
            out += [f"## {head}", ""]
            for x in items:
                title = x.get("title") or x.get("name") or x["id"]
                out += [f"### {title}", "", f'<img src="../svg/{files[(kind, x["id"])][0]["file"]}.svg" alt="{title}">', ""]
                if kind == "reaction" and x.get("card_zh"):
                    out += [card_md(x["card_zh"]), ""]
                if kind == "molecule" and x.get("notes_zh"):
                    out += [x["notes_zh"], ""]
        (md_dir / f"{tag}.md").write_text("\n".join(out), encoding="utf-8")


TOPIC_DIR = "分類"
KIND_DIR = {"molecule": "分子", "reaction": "反應", "pathway": "路徑"}


def topics_of(db):
    """{(kind, id): {tag, …}} — a molecule belongs to its own tags plus those of every reaction
    and pathway it appears in; a reaction also to the tags of the pathways it is a step of."""
    out = {}
    for p in db.pathways.values():
        out[("pathway", p["id"])] = set(p["tags"])
    for r in db.reactions.values():
        t = set(r["tags"])
        for p in db.pathways.values():
            if any(e.get("step") == r["id"] for e in p.get("edges", [])):
                t |= set(p["tags"])
        out[("reaction", r["id"])] = t
    for m in db.molecules.values():
        t = set(m["tags"])
        for r in db.reactions.values():
            if any(x["id"] == m["id"] for x in r["substrates"] + r["products"]):
                t |= set(r["tags"])
        for p in db.pathways.values():
            if any(n.get("ref") == m["id"] for n in p.get("nodes", [])):
                t |= set(p["tags"])
        out[("molecule", m["id"])] = t
    return out


def topic_folder(tag, tt):
    """胺基酸/鹼性 for a child tag, 糖解 for a top-level one (Chinese names from tags.yaml)."""
    zh = lambda t: tt["labels"].get(t, t)                                   # noqa: E731
    parent = next((p for p, kids in tt["tree"].items() if tag in kids), None)
    return f"{zh(parent)}/{zh(tag)}" if parent else zh(tag)


def _link(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(src, dst)                         # hard link: no extra disk space
    except OSError:
        shutil.copy2(src, dst)


def write_topics(db, files):
    """out/分類/<主題>/{分子,反應,路徑}/<id>[.<view>].png and …/拆件/<id>/…png — the same images
    sorted by topic for browsing on the iPad (out/svg and out/png keep their paths for notes)."""
    root = OUT / TOPIC_DIR
    shutil.rmtree(root, ignore_errors=True)
    tt = tag_tree()
    for (kind, id_), tags in topics_of(db).items():
        # a child tag (酸性) already says where it belongs: skip its parent (胺基酸) as a folder
        kids = {t for t in tags if any(t in v for v in tt["tree"].values())}
        parents = {p for p, v in tt["tree"].items() if kids & set(v)}
        for tag in sorted(tags - parents):
            base = root / topic_folder(tag, tt)
            for f in files.get((kind, id_), []):
                if f["view"] in ("skeletal", "formula") and len(files[(kind, id_)]) > 1:
                    continue                      # automatic extra views: not worth a sticker
                src = OUT / "png" / f"{f['file']}.png"
                if src.exists():
                    _link(src, base / KIND_DIR[kind] / f"{Path(f['file']).name}.png")
            pdir = OUT / "parts" / f"{kind}s" / id_
            if pdir.is_dir():
                for f in pdir.rglob("*"):         # PNG stickers; SVG for `parts_png: false` items
                    if f.suffix == ".png" or (f.suffix == ".svg" and not f.with_suffix(".png").exists()):
                        _link(f, base / "拆件" / id_ / f.relative_to(pdir))


def build(db, force=False):
    errors, warnings, info = validate(db)
    if errors and not force:
        return errors, warnings
    problems = []
    files = render_all(db, problems)
    build_index(db, files, info)
    build_galleries(db, files)
    write_topics(db, files)
    return errors + problems, warnings


def refresh_index(db):
    """Cheap re-index after a review action (no re-rendering)."""
    _, _, info = validate(db)
    files = {}
    for kind, tbl in (("molecule", db.molecules), ("reaction", db.reactions), ("pathway", db.pathways)):
        for x in tbl.values():
            if kind == "molecule":
                files[(kind, x["id"])] = [{"view": v, "file": svg_name(kind, x["id"], v, x["default_view"])}
                                          for v in views_of(x)]
            else:
                files[(kind, x["id"])] = [{"view": "default", "file": svg_name(kind, x["id"])}] + \
                    [{"view": v, "file": f"{svg_name(kind, x['id'])}.{v}"} for v in (x.get("variants") or {})]
            parts = {}
            for view in ["default"] + list(x.get("variants") or {}):
                sub = x["id"] if view == "default" else f"{x['id']}/{view}"
                folder = OUT / "parts" / f"{kind}s" / sub
                if folder.is_dir():
                    names = sorted((f.stem for f in folder.glob("*.svg") if f.stem != "_sheet"),
                                   key=lambda n: (_part_order(n), n))
                    parts[view] = [{"name": n, "file": f"{kind}s/{sub}/{n}"} for n in names]
            files[(f"{kind}-parts", x["id"])] = parts
    build_index(db, files, info)


def build_site(dest=None):
    """Static copy for GitHub Pages: web UI + index.json + SVGs only (no PNG / stickers / review data)."""
    dest = Path(dest or ROOT / "site")
    shutil.rmtree(dest, ignore_errors=True)
    shutil.copytree(ROOT / "web", dest / "web", ignore=shutil.ignore_patterns(".DS_Store"))
    shutil.move(dest / "web" / "index.html", dest / "index.html")
    shutil.copytree(OUT / "svg", dest / "out" / "svg", ignore=shutil.ignore_patterns(".DS_Store", "* 2", "* 2.*"))
    idx = json.loads((OUT / "index.json").read_text(encoding="utf-8"))
    idx["static"] = True
    idx.pop("abs_out", None)
    for kind in ("molecules", "reactions", "pathways"):
        for x in idx[kind]:
            x["comments"], x["parts"] = [], {}
    (dest / "out" / "index.json").write_text(json.dumps(idx, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    (dest / ".nojekyll").write_text("")
    return dest
