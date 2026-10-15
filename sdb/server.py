"""Local web UI + review API (standard library only)."""
import json
import mimetypes
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import unquote, urlparse

from . import export, model, review
from .style import ROOT

STATIC = {"/web/": ROOT / "web", "/out/": export.OUT, "/reference/": ROOT / "reference"}
LOCK = threading.Lock()


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):
        pass

    def _send(self, code, body, ctype="application/json; charset=utf-8"):
        if isinstance(body, (dict, list)):
            body = json.dumps(body, ensure_ascii=False)
        if isinstance(body, str):
            body = body.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        p = unquote(urlparse(self.path).path)
        if p == "/":
            p = "/web/index.html"
        for prefix, base in STATIC.items():
            if p.startswith(prefix):
                f = (base / p[len(prefix):]).resolve()
                if base.resolve() in f.parents and f.is_file():
                    ctype = mimetypes.guess_type(f.name)[0] or "application/octet-stream"
                    if ctype.startswith("text/") or ctype in ("application/javascript", "image/svg+xml", "application/json"):
                        ctype += "; charset=utf-8"
                    return self._send(200, f.read_bytes(), ctype)
        self._send(404, {"error": "not found"})

    def do_POST(self):
        p = urlparse(self.path).path
        n = int(self.headers.get("Content-Length") or 0)
        try:
            data = json.loads(self.rfile.read(n) or b"{}")
        except json.JSONDecodeError:
            return self._send(400, {"error": "bad json"})
        with LOCK:
            if p == "/api/review":
                db = model.load()
                kind, id_, action = data.get("kind"), data.get("id"), data.get("action")
                if (kind, id_) not in db.files or action not in ("approve", "comment", "reopen"):
                    return self._send(400, {"error": "bad request"})
                text = (data.get("text") or "").strip()
                if action == "comment" and not text:
                    return self._send(400, {"error": "意見不可空白"})
                c = review.add(kind, id_, action, text)
                status = {"approve": "approved", "comment": "changes_requested", "reopen": "pending"}[action]
                model.set_status(db, kind, id_, status)
                export.refresh_index(db)
                return self._send(200, {"ok": True, "comment": c, "status": status})
            if p == "/api/rebuild":
                errors, warnings = export.build(model.load())
                return self._send(200, {"ok": not errors, "errors": errors, "warnings": warnings})
        self._send(404, {"error": "not found"})


def serve(port=8765):
    if not (export.OUT / "index.json").exists():
        export.build(model.load(), force=True)
    srv = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    print(f"structure_db → http://localhost:{port}   (Ctrl+C 停止)")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
