"""python -m sdb {validate|build|serve|site|review|resolve}"""
import argparse
import sys

from . import model, review


def main(argv=None):
    ap = argparse.ArgumentParser(prog="sdb")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("validate", help="檢查 SMILES、反應平衡、Fischer 立體化學、id 引用")
    b = sub.add_parser("build", help="產出 out/ 下所有 SVG、PNG、md 圖庫、index.json")
    b.add_argument("--force", action="store_true", help="即使驗證有錯也照樣出圖")
    s = sub.add_parser("serve", help="啟動本機網頁")
    s.add_argument("--port", type=int, default=8765)
    t = sub.add_parser("site", help="產出 GitHub Pages 用的靜態網頁 site/（只有 SVG，不含審核）")
    t.add_argument("--dest", default=None)
    sub.add_parser("review", help="列出尚未處理的審核意見（給 AI 用）")
    r = sub.add_parser("resolve", help="標記審核意見已處理，並把項目狀態改回 pending")
    r.add_argument("cid")
    r.add_argument("reply", nargs="?", default="")
    a = ap.parse_args(argv)

    if a.cmd == "validate":
        from .validate import validate
        errors, warnings, _ = validate(model.load())
        for w in warnings:
            print("⚠️ ", w)
        for e in errors:
            print("❌", e)
        print(f"{len(errors)} 個錯誤，{len(warnings)} 個警告")
        return 1 if errors else 0

    if a.cmd == "build":
        from .export import build
        errors, warnings = build(model.load(), force=a.force)
        for w in warnings:
            print("⚠️ ", w)
        for e in errors:
            print("❌", e)
        print("建置失敗" if errors and not a.force else "建置完成 → out/", f"（{len(errors)} 個錯誤）")
        return 1 if errors else 0

    if a.cmd == "serve":
        from .server import serve
        serve(a.port)
        return 0

    if a.cmd == "site":
        from .export import build_site
        print("靜態網頁 →", build_site(a.dest))
        return 0

    if a.cmd == "review":
        items = review.open_comments()
        if not items:
            print("沒有待處理的審核意見。")
        for c in items:
            print(f"[{c['cid']}] {c['kind']}:{c['id']}  ({c['ts']})\n    {c['text']}")
        return 0

    if a.cmd == "resolve":
        db = model.load()
        c = review.resolve(a.cid, a.reply)
        if (c["kind"], c["id"]) in db.files:
            model.set_status(db, c["kind"], c["id"], "pending")
        from .export import refresh_index
        refresh_index(db)
        print(f"已處理 {a.cid}，{c['kind']}:{c['id']} 狀態 → pending")
        return 0


if __name__ == "__main__":
    sys.exit(main())
