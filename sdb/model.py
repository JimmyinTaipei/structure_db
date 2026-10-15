"""Load data/*.yaml into plain dicts with defaults, and check cross references."""
import re
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from .style import ROOT

DATA = ROOT / "data"
KINDS = {"molecule": "molecules", "reaction": "reactions", "pathway": "pathways"}
STATUSES = ("draft", "pending", "changes_requested", "approved")


def participant(p):
    """'pyruvate' or {id: pyruvate, n: 2, view: skeletal} → normalised dict."""
    if isinstance(p, str):
        return {"id": p, "n": 1}
    return {"n": 1, **p}


@dataclass
class DB:
    molecules: dict = field(default_factory=dict)
    reactions: dict = field(default_factory=dict)
    pathways: dict = field(default_factory=dict)
    cofactors: dict = field(default_factory=dict)
    files: dict = field(default_factory=dict)      # (kind, id) → Path
    errors: list = field(default_factory=list)

    def table(self, kind):
        return getattr(self, KINDS[kind])

    def get(self, kind, id_):
        return self.table(kind).get(id_)

    def reactions_of(self, mol_id):
        out = []
        for r in self.reactions.values():
            ids = [p["id"] for p in r["substrates"] + r["products"]]
            if mol_id in ids:
                out.append(r["id"])
        return out


def _load_dir(db, kind):
    for path in sorted((DATA / KINDS[kind]).glob("*.yaml")):
        try:
            doc = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        except yaml.YAMLError as e:
            db.errors.append(f"{path.name}: YAML 語法錯誤 {e}")
            continue
        id_ = doc.get("id")
        if id_ != path.stem:
            db.errors.append(f"{path.name}: id ({id_}) 必須與檔名相同")
            continue
        doc.setdefault("tags", [])
        doc.setdefault("status", "draft")
        if doc["status"] not in STATUSES:
            db.errors.append(f"{path.name}: 未知 status {doc['status']}")
        doc["kind"] = kind
        db.table(kind)[id_] = doc
        db.files[(kind, id_)] = path


def load():
    db = DB()
    cof = DATA / "cofactors.yaml"
    if cof.exists():
        db.cofactors = yaml.safe_load(cof.read_text(encoding="utf-8")) or {}
    for kind in KINDS:
        _load_dir(db, kind)

    for m in db.molecules.values():
        m.setdefault("name", m["id"])
        m.setdefault("aliases", [])
        m.setdefault("default_view", "lipid" if m.get("lipid") else "haworth" if m.get("haworth") else "chair" if m.get("chair") else "ring" if m.get("ring")
                     else "free" if m.get("free") else "chain" if m.get("chain")
                     else ("formula" if m.get("formula") else "skeletal"))
        if not m.get("smiles"):
            db.errors.append(f"molecule {m['id']}: 缺少 smiles")

    for r in db.reactions.values():
        r["substrates"] = [participant(p) for p in r.get("substrates", [])]
        r["products"] = [participant(p) for p in r.get("products", [])]
        for k in ("cofactors_in", "cofactors_out", "prosthetic"):
            r.setdefault(k, [])
        r.setdefault("reversible", False)
        r.setdefault("decorate", {})
        if r.get("text") is not False:                   # names-only versions of every reaction
            v = r.get("variants") or {}
            v.setdefault("text", {"auto": True, "note_zh": R_TEXT_NOTE})
            v.setdefault("text_full", {"auto": True, "note_zh": R_TEXT_FULL_NOTE})
            r["variants"] = v
        ov = DATA / "overrides" / f"{r['id']}.yaml"       # hand tweaks kept out of the main file
        if ov.exists():
            r["layout"] = {**r.get("layout", {}), **(yaml.safe_load(ov.read_text(encoding="utf-8")) or {})}
        for p in r["substrates"] + r["products"]:
            if p["id"] not in db.molecules:
                db.errors.append(f"reaction {r['id']}: 找不到分子 {p['id']}")
            elif p.get("chain") and p["chain"] not in (db.molecules[p["id"]].get("chains") or {}):
                db.errors.append(f"reaction {r['id']}: 分子 {p['id']} 沒有名為 {p['chain']} 的排法（chains）")
        for c in r["cofactors_in"] + r["cofactors_out"]:
            if c not in db.cofactors and c not in db.molecules:
                db.errors.append(f"reaction {r['id']}: 找不到輔因子 {c}（請加到 data/cofactors.yaml）")

    for p in db.pathways.values():
        p.setdefault("layout", "linear")
        ids = set()
        for n in p.get("nodes", []):
            ids.add(n["id"])
            if n.get("ref") and n["ref"] not in db.molecules:
                db.errors.append(f"pathway {p['id']}: 節點 {n['id']} 引用不存在的分子 {n['ref']}")
        for e in p.get("edges", []):
            for end in ("from", "to"):
                v = e.get(end, "")
                if not v.startswith("edge:") and v not in ids:
                    db.errors.append(f"pathway {p['id']}: edge 的 {end}={v} 不是已定義的節點")
            if e.get("step") and e["step"] not in db.reactions:
                db.errors.append(f"pathway {p['id']}: step {e['step']} 不是已定義的反應")
        add_text_variants(p)
    return db


R_TEXT_NOTE = "文字版：只寫物質名稱與箭頭，不放酵素與輔因子。"
R_TEXT_FULL_NOTE = "文字版＋酵素：不畫結構式，輔因子與酵素名稱保留。拆件區可單獨取酵素層。"
TEXT_NOTE = "文字版：只有物質名稱與箭頭，不放酵素與輔因子（指出路徑、接其他路徑用）。"
TEXT_FULL_NOTE = "文字版＋酵素：不畫結構式，酵素名稱與輔因子保留（省空間的完整版）。拆件區可單獨取酵素層。"


def add_text_variants(p):
    """Every pathway with molecule nodes gets text versions for free (opt out with `text: false`).
    text_full only when some node is drawn as a structure (otherwise it equals the default)."""
    refs = [n for n in p.get("nodes", []) if n.get("ref")]
    if not refs or p.get("text") is False:
        return
    v = p.setdefault("variants", {}) or {}
    p["variants"] = v
    v.setdefault("text", {"auto": True, "note_zh": TEXT_NOTE})
    if any(n.get("view") != "name" for n in refs):
        v.setdefault("text_full", {"auto": True, "note_zh": TEXT_FULL_NOTE})


STATUS_RE = re.compile(r"^status:\s*\S+\s*$", re.M)


def set_status(db, kind, id_, status):
    """Rewrite only the `status:` line so comments/formatting in the YAML survive."""
    if status not in STATUSES:
        raise ValueError(status)
    path = db.files[(kind, id_)]
    text = path.read_text(encoding="utf-8")
    if STATUS_RE.search(text):
        text = STATUS_RE.sub(f"status: {status}", text, count=1)
    else:
        text = text.rstrip("\n") + f"\nstatus: {status}\n"
    path.write_text(text, encoding="utf-8")
    db.table(kind)[id_]["status"] = status
