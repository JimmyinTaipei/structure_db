"""Review comments: review/comments.jsonl (one JSON object per line, append-only except `resolved`)."""
import json
from datetime import datetime

from .style import ROOT

FILE = ROOT / "review" / "comments.jsonl"


def load():
    if not FILE.exists():
        return []
    return [json.loads(l) for l in FILE.read_text(encoding="utf-8").splitlines() if l.strip()]


def _save(items):
    FILE.parent.mkdir(exist_ok=True)
    FILE.write_text("".join(json.dumps(c, ensure_ascii=False) + "\n" for c in items), encoding="utf-8")


def add(kind, id_, action, text=""):
    items = load()
    c = {"cid": f"c{len(items) + 1:04d}", "ts": datetime.now().isoformat(timespec="seconds"),
         "kind": kind, "id": id_, "action": action, "text": text, "resolved": action != "comment"}
    items.append(c)
    _save(items)
    return c


def open_comments():
    return [c for c in load() if c["action"] == "comment" and not c.get("resolved")]


def resolve(cid, reply=""):
    items = load()
    for c in items:
        if c["cid"] == cid:
            c["resolved"] = True
            c["reply"] = reply
            c["resolved_ts"] = datetime.now().isoformat(timespec="seconds")
            _save(items)
            return c
    raise KeyError(cid)
