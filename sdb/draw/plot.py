"""Plots: axes + curves (enzyme kinetics, O₂-binding curves, titration, Ramachandran).

    layout: plot
    plot:                               # or  panels: [{title: "…", <plot keys>}, …]  side by side
      size: [14, 9]                     # plot area, in S
      params: {Vmax: 1, Km: 1}          # names usable in every expression below
      x: {min: 0, max: 10, label: "[S]", ticks: [0, 5, 10], marks: [{at: Km, label: "K_m"}]}
      y: {min: 0, max: 1.2, label: "v_0"}        # axes cross at 0 when 0 is inside the range
      curves:
        - {id: mm, y: "mm(x, Vmax, Km)", from: 0, to: 10, color: black, dash: false,
           label: "No inhibitor", label_at: 8, label_off: [0, -0.7]}
        - {id: tit, t: [0, 14], x: "…", y: "t"}  # parametric curve (t runs over the range)
      points: [{on: mm, x: Km, y: "Vmax/2", drop: true, label: "…", dot: true, leader: false}]
      hlines: [{y: Vmax, label: "V_max"}]        # dashed guides across the plot
      vlines: [{x: 40, label: "Tissues"}]
      bands:  [{x0: 20, x1: 40, label: "Tissues", color: cyan}]       # shaded strip
      areas:  [{points: [[x, y], …], color: cyan, label: "α", label_at: [x, y]}]   # polygons
      arrows: [{from: [x, y], to: [x, y], color: red, label: "…"}]
      texts:  [{at: [x, y], text: "…", color: …, align: left|center|right}]
      legend: {at: [x, y], items: [curve_id | {curve, text, color}], row: 1.25}   # swatch + name rows

Every number may be an expression of the params (e.g. "-1/Km"). A point on a curve needs only
one of x / y (the other is read off the curve); when both are given the validator checks that
the point really lies on the curve — that is what keeps K_m, P50, intercepts … honest.

Built-in functions: mm(s, Vmax, Km), hill(p, P50, n), mwc(p, KR, KT, L, n=4), logit(Y),
bell(x, mu, sigma), exp, log, log10, sqrt.

Pathway-level keys: panel_gap (S), caption: {text, color} — one line under the whole figure.
Parts: axes, curves, guides (dashed lines, bands, areas), labels, arrows (+ panel_<n> per panel).
"""
import math

from .. import style
from ..graphic import BBox, Graphic, line, path
from ..text import TextBlock

S = style.S
CURVE_W = 1.7 * style.LINE_W
N_SAMPLES = 600


def mm(s, Vmax, Km):
    return Vmax * s / (Km + s)


def hill(p, P50, n):
    return 0.0 if p <= 0 else p ** n / (P50 ** n + p ** n)


def mwc(p, KR, KT, L, n=4):
    a, c = p / KR, KR / KT
    return (a * (1 + a) ** (n - 1) + L * c * a * (1 + c * a) ** (n - 1)) / ((1 + a) ** n + L * (1 + c * a) ** n)


def logit(Y):
    return math.log10(Y / (1 - Y))


def bell(x, mu, sigma):
    return math.exp(-0.5 * ((x - mu) / sigma) ** 2)


FUNCS = {"mm": mm, "hill": hill, "mwc": mwc, "logit": logit, "bell": bell, "exp": math.exp,
         "log": math.log, "log10": math.log10, "sqrt": math.sqrt, "pi": math.pi, "abs": abs,
         "min": min, "max": max}


def ev(v, env):
    """Number, or an expression string of the params."""
    if isinstance(v, (int, float)):
        return float(v)
    return float(eval(str(v), {"__builtins__": {}}, {**FUNCS, **env}))   # data files are our own


class Axes:
    def __init__(self, spec):
        self.spec = spec
        self.env = dict(spec.get("params") or {})
        for k, v in list(self.env.items()):               # params may refer to earlier params
            self.env[k] = ev(v, self.env)
        w, h = spec.get("size", [14, 9])
        self.W, self.H = w * S, h * S
        X, Y = spec["x"], spec["y"]
        self.x0, self.x1 = ev(X["min"], self.env), ev(X["max"], self.env)
        self.y0, self.y1 = ev(Y["min"], self.env), ev(Y["max"], self.env)
        # where the axes cross: 0 when inside the range, else the low end
        self.cx = ev(Y.get("at"), self.env) if Y.get("at") is not None else (0.0 if self.x0 <= 0 <= self.x1 else self.x0)
        self.cy = ev(X.get("at"), self.env) if X.get("at") is not None else (0.0 if self.y0 <= 0 <= self.y1 else self.y0)

    def X(self, x):
        return (x - self.x0) / (self.x1 - self.x0) * self.W

    def Y(self, y):
        return self.H - (y - self.y0) / (self.y1 - self.y0) * self.H

    def P(self, x, y):
        return self.X(x), self.Y(y)

    def e(self, v):
        return ev(v, self.env)

    # ------------------------------------------------------------ curves as sampled polylines
    def samples(self, c):
        """Unclipped [(x, y)] in data units."""
        env = self.env
        out = []
        if "t" in c:
            t0, t1 = (self.e(v) for v in c["t"])
            for i in range(N_SAMPLES + 1):
                t = t0 + (t1 - t0) * i / N_SAMPLES
                try:
                    out.append((ev(c["x"], {**env, "t": t}), ev(c["y"], {**env, "t": t})))
                except (ValueError, ZeroDivisionError, OverflowError):
                    out.append(None)
            return out
        a = self.e(c.get("from", self.x0))
        b = self.e(c.get("to", self.x1))
        for i in range(N_SAMPLES + 1):
            x = a + (b - a) * i / N_SAMPLES
            try:
                out.append((x, ev(c["y"], {**env, "x": x})))
            except (ValueError, ZeroDivisionError, OverflowError):
                out.append(None)
        return out

    def crossing(self, c, along, value, which=0):
        """Point of curve c where x (along='x') or y equals value; linear interpolation."""
        if along == "x" and "t" not in c:
            return value, ev(c["y"], {**self.env, "x": value})
        k = 0 if along == "x" else 1
        pts = self.samples(c)
        hits = []
        for p, q in zip(pts, pts[1:]):
            if p is None or q is None:
                continue
            if (p[k] - value) * (q[k] - value) <= 0 and p[k] != q[k]:
                f = (value - p[k]) / (q[k] - p[k])
                hits.append((p[0] + f * (q[0] - p[0]), p[1] + f * (q[1] - p[1])))
        if not hits:
            raise ValueError(f"曲線 {c.get('id')} 沒有通過 {along} = {value:g}")
        return hits[min(which, len(hits) - 1)]

    def distance(self, c, x, y):
        """Smallest distance (in S) from the point to curve c."""
        px, py = self.P(x, y)
        best = float("inf")
        pts = [q for q in self.samples(c) if q]
        for (ax, ay), (bx, by) in zip(pts, pts[1:]):
            (ax, ay), (bx, by) = self.P(ax, ay), self.P(bx, by)
            dx, dy = bx - ax, by - ay
            L2 = dx * dx + dy * dy or 1e-12
            t = max(0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / L2))
            best = min(best, math.hypot(ax + t * dx - px, ay + t * dy - py))
        return best / S


def _clip(p, q, box):
    """Liang–Barsky: the part of segment p→q inside box (x0, y0, x1, y1), or None."""
    x0, y0, x1, y1 = box
    dx, dy = q[0] - p[0], q[1] - p[1]
    t0, t1 = 0.0, 1.0
    for pp, qq in ((-dx, p[0] - x0), (dx, x1 - p[0]), (-dy, p[1] - y0), (dy, y1 - p[1])):
        if pp == 0:
            if qq < 0:
                return None
            continue
        r = qq / pp
        if pp < 0:
            t0 = max(t0, r)
        else:
            t1 = min(t1, r)
        if t0 > t1:
            return None
    return (p[0] + t0 * dx, p[1] + t0 * dy), (p[0] + t1 * dx, p[1] + t1 * dy)


def _fmt(v):
    return f"{v:g}".replace("-", "−")


def _dash(svg, pattern="5 4"):
    return svg.replace("/>", f' stroke-dasharray="{pattern}"/>', 1)


def _head(x, y, dx, dy, color=None, size=0.45 * S):
    n = math.hypot(dx, dy) or 1
    ux, uy = dx / n, dy / n
    a = math.radians(22)
    pts = []
    for s in (1, -1):
        c, si = math.cos(a), s * math.sin(a)
        pts.append((x - size * (ux * c - uy * si), y - size * (uy * c + ux * si)))
    return path(f"M{pts[0][0]:.2f},{pts[0][1]:.2f} L{x:.2f},{y:.2f} L{pts[1][0]:.2f},{pts[1][1]:.2f}", stroke=color)


class Layers:
    NAMES = ("guides", "axes", "curves", "arrows", "labels")

    def __init__(self):
        self.d = {k: ([], BBox()) for k in self.NAMES}

    def add(self, layer, svg, x0, y0, x1, y1):
        self.d[layer][0].append(svg)
        self.d[layer][1].add(x0, y0, x1, y1)

    def text(self, tb, cx, cy, fill=None, align="center", layer="labels"):
        """TextBlock with its visual middle at (cx, cy); align left/right puts (cx, cy) on that edge
        and aligns every line of a multi-line block to it."""
        base = cy - (tb.top + tb.bottom) / 2
        tb.align = align
        x0 = cx if align == "left" else cx - tb.width if align == "right" else cx - tb.width / 2
        self.add(layer, tb.svg(cx, base, fill), x0, base + tb.top, x0 + tb.width, base + tb.bottom)

    def graphic(self):
        bb = BBox()
        body = []
        for k in self.NAMES:
            svg, b = self.d[k]
            if svg:
                body += svg
                bb.add(b.x0, b.y0, b.x1, b.y1)
        parts = {k: Graphic("".join(svg), b.x0, b.y0, b.x1, b.y1) for k, (svg, b) in self.d.items() if svg}
        return Graphic("".join(body), bb.x0, bb.y0, bb.x1, bb.y1, {}, parts)


def label_tb(text, size=None, kind="sans"):
    return TextBlock(str(text), kind, size or style.ENZYME_SIZE)


def draw_axes(ax, L):
    spec = ax.spec
    X, Y = spec["x"], spec["y"]
    ox, oy = ax.X(ax.cx), ax.Y(ax.cy)
    over = 0.6 * S                                       # arrow tips stick out past the range
    L.add("axes", line(0, oy, ax.W + over, oy) + _head(ax.W + over, oy, 1, 0), 0, oy - 0.3 * S, ax.W + over, oy + 0.3 * S)
    L.add("axes", line(ox, ax.H, ox, -over) + _head(ox, -over, 0, -1), ox - 0.3 * S, -over, ox + 0.3 * S, ax.H)
    tick = 0.25 * S
    num = style.ENZYME_SIZE * 0.85
    below = oy + tick + 0.55 * S
    for v in X.get("ticks") or []:
        x = ax.X(ax.e(v))
        L.add("axes", line(x, oy, x, oy + tick), x, oy, x, oy + tick)
        if X.get("tick_labels", True):
            L.text(label_tb(_fmt(ax.e(v)), num), x, below, layer="axes")
    left = ox - tick - 0.3 * S
    for v in Y.get("ticks") or []:
        y = ax.Y(ax.e(v))
        L.add("axes", line(ox - tick, y, ox, y), ox - tick, y, ox, y)
        if Y.get("tick_labels", True):
            L.text(label_tb(_fmt(ax.e(v)), num), left, y, align="right", layer="axes")
    for m in X.get("marks") or []:                      # named positions: K_m, −1/K_m, P50 …
        x = ax.X(ax.e(m["at"]))
        L.add("axes", line(x, oy - tick * 0.6, x, oy + tick), x, oy, x, oy + tick)
        dy = m.get("dy", 0) * S
        side = -1 if m.get("above") else 1
        y = oy + side * (tick + 0.6 * S) + dy
        L.text(label_tb(m["label"]), x + m.get("dx", 0) * S, y, m.get("color"), layer="labels")
    for m in Y.get("marks") or []:
        y = ax.Y(ax.e(m["at"]))
        L.add("axes", line(ox - tick, y, ox + tick * 0.6, y), ox - tick, y, ox, y)
        ly = y + m.get("dy", 0) * S
        if m.get("right"):
            L.text(label_tb(m["label"]), ox + tick + 0.3 * S + m.get("dx", 0) * S, ly, m.get("color"), align="left")
        else:
            L.text(label_tb(m["label"]), ox - tick - 0.3 * S + m.get("dx", 0) * S, ly, m.get("color"), align="right")
    lab = style.NAME_SIZE
    if X.get("label"):
        tb = label_tb(X["label"], lab)
        if X.get("label_pos", "center" if ax.cy == ax.y0 else "end") == "center":
            ty = oy + tick + (1.6 * S if X.get("ticks") else 0.9 * S) + (0.6 * S if X.get("marks") else 0)
            L.text(tb, ax.W / 2, ty, layer="axes")
        else:
            L.text(tb, ax.W + over + 0.35 * S, oy, align="left", layer="axes")
    if Y.get("label"):
        tb = label_tb(Y["label"], lab)
        if Y.get("label_pos", "end") == "center":           # rotated along the axis, left of the numbers
            widest = max([label_tb(_fmt(ax.e(v)), num).width for v in Y.get("ticks") or []] + [0])
            cx = left - widest - 0.5 * S - (tb.bottom - tb.top) / 2
            cy = ax.H / 2
            base = -(tb.top + tb.bottom) / 2
            svg = f'<g transform="translate({cx:.2f},{cy:.2f}) rotate(-90)">{tb.svg(0, base)}</g>'
            h = tb.bottom - tb.top
            L.add("axes", svg, cx - h / 2, cy - tb.width / 2, cx + h / 2, cy + tb.width / 2)
        else:
            L.text(tb, ox, -over - 0.45 * S - (tb.bottom - tb.top) / 2, layer="axes")


def draw_plot(spec):
    ax = Axes(spec)
    L = Layers()
    box = (0, 0, ax.W, ax.H)
    curves = {c["id"]: c for c in spec.get("curves") or []}

    for b in spec.get("bands") or []:
        xa, xb = sorted((ax.X(ax.e(b["x0"])), ax.X(ax.e(b["x1"]))))
        col = style.color(b.get("color", "cyan"))
        L.add("guides", f'<rect x="{xa:.2f}" y="0" width="{xb - xa:.2f}" height="{ax.H:.2f}" fill="{col}" '
                        f'fill-opacity="{b.get("opacity", 0.12)}"/>', xa, 0, xb, ax.H)
        if b.get("label"):
            L.text(label_tb(b["label"]), (xa + xb) / 2, -0.15 * S + b.get("dy", 0) * S - 0.5 * S, b.get("label_color", "gray"))
    for a in spec.get("areas") or []:
        pts = [ax.P(ax.e(x), ax.e(y)) for x, y in a["points"]]
        d = "M" + " L".join(f"{x:.2f},{y:.2f}" for x, y in pts) + " Z"
        col = style.color(a.get("color", "cyan"))
        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        L.add("guides", f'<path d="{d}" fill="{col}" fill-opacity="{a.get("opacity", 0.22)}" stroke="{col}" '
                        f'stroke-width="{style.LINE_W:.2f}"/>', min(xs), min(ys), max(xs), max(ys))
        if a.get("label"):
            lx, ly = ax.P(*(ax.e(v) for v in a["label_at"])) if a.get("label_at") else (sum(xs) / len(xs), sum(ys) / len(ys))
            L.text(label_tb(a["label"]), lx, ly, a.get("label_color"))
    for g in spec.get("hlines") or []:
        y = ax.Y(ax.e(g["y"]))
        x0 = ax.X(ax.e(g.get("from", ax.x0)))
        x1 = ax.X(ax.e(g.get("to", ax.x1)))
        L.add("guides", _dash(line(x0, y, x1, y, stroke=g.get("color", "gray"))), x0, y, x1, y)
        if g.get("label"):
            tb = label_tb(g["label"])
            off = g.get("label_off", [0, 0])
            if g.get("label_pos") == "left":
                L.text(tb, x0 - 0.3 * S + off[0] * S, y + off[1] * S, g.get("label_color"), align="right")
            else:
                L.text(tb, x1 + 0.3 * S + off[0] * S, y + off[1] * S, g.get("label_color"), align="left")
    for g in spec.get("vlines") or []:
        x = ax.X(ax.e(g["x"]))
        y0 = ax.Y(ax.e(g.get("from", ax.y0)))
        y1 = ax.Y(ax.e(g.get("to", ax.y1)))
        L.add("guides", _dash(line(x, y0, x, y1, stroke=g.get("color", "gray"))), x, y1, x, y0)
        if g.get("label"):
            off = g.get("label_off", [0, 0])
            L.text(label_tb(g["label"]), x + off[0] * S, y1 - 0.55 * S + off[1] * S, g.get("label_color"))

    draw_axes(ax, L)

    for c in curves.values():
        segs, cur = [], []
        pts = [ax.P(*q) if q else None for q in ax.samples(c)]
        for p, q in zip(pts, pts[1:]):
            s = _clip(p, q, box) if p and q else None
            if s is None:
                if cur:
                    segs.append(cur)
                cur = []
                continue
            if cur and math.hypot(cur[-1][0] - s[0][0], cur[-1][1] - s[0][1]) < 1e-6:
                cur.append(s[1])
            else:
                if cur:
                    segs.append(cur)
                cur = [s[0], s[1]]
        if cur:
            segs.append(cur)
        col = c.get("color")
        for sg in segs:
            d = "M" + " L".join(f"{x:.2f},{y:.2f}" for x, y in sg)
            svg = path(d, stroke=col, w=c.get("width", 1) * CURVE_W)
            if c.get("dash"):
                svg = _dash(svg, "7 5")
            xs, ys = [p[0] for p in sg], [p[1] for p in sg]
            L.add("curves", svg, min(xs), min(ys), max(xs), max(ys))
        if c.get("label"):
            if "label_at" in c:
                lx, ly = ax.P(*ax.crossing(c, "x", ax.e(c["label_at"]))) if "t" not in c else \
                    ax.P(*ax.crossing(c, "y", ax.e(c["label_at"])))
            else:
                lx, ly = segs[-1][-1]
            off = c.get("label_off", [0.3, -0.6])
            L.text(label_tb(c["label"]), lx + off[0] * S, ly + off[1] * S, c.get("label_color", col),
                   c.get("label_align", "left"))

    for p in spec.get("points") or []:
        c = curves.get(p.get("on"))
        if c is not None and ("x" not in p or "y" not in p):
            along = "x" if "x" in p else "y"
            x, y = ax.crossing(c, along, ax.e(p[along]), p.get("which", 0))
        else:
            x, y = ax.e(p["x"]), ax.e(p["y"])
        px, py = ax.P(x, y)
        col = p.get("color")
        drop = p.get("drop", False)
        if drop in (True, "x", "both"):
            oy = ax.Y(ax.cy)
            L.add("guides", _dash(line(px, py, px, oy, stroke=p.get("drop_color", "gray"))), px, min(py, oy), px, max(py, oy))
        if drop in (True, "y", "both"):
            ox = ax.X(ax.cx)
            L.add("guides", _dash(line(ox, py, px, py, stroke=p.get("drop_color", "gray"))), min(ox, px), py, max(ox, px), py)
        if p.get("dot", True):
            r = 0.17 * S
            L.add("curves", f'<circle cx="{px:.2f}" cy="{py:.2f}" r="{r:.2f}" fill="{style.color(col)}"/>',
                  px - r, py - r, px + r, py + r)
        if p.get("label"):
            off = p.get("label_off", [0.4, -0.5])
            if p.get("leader"):                          # thin line from the dot to a far-away label
                lx, ly = px + off[0] * S, py + off[1] * S
                L.add("labels", line(px, py, lx - 0.15 * S * (1 if off[0] >= 0 else -1), ly, stroke="gray",
                                     w=0.6 * style.LINE_W), min(px, lx), min(py, ly), max(px, lx), max(py, ly))
            L.text(label_tb(p["label"]), px + off[0] * S, py + off[1] * S, p.get("label_color", col),
                   p.get("label_align", "left"))

    for a in spec.get("arrows") or []:
        (fx, fy), (tx, ty) = ax.P(*(ax.e(v) for v in a["from"])), ax.P(*(ax.e(v) for v in a["to"]))
        col = a.get("color")
        svg = line(fx, fy, tx, ty, stroke=col, w=1.3 * style.LINE_W) + _head(tx, ty, tx - fx, ty - fy, col)
        L.add("arrows", svg, min(fx, tx), min(fy, ty), max(fx, tx), max(fy, ty))
        if a.get("label"):
            off = a.get("label_off", [0, -0.6])
            L.text(label_tb(a["label"]), (fx + tx) / 2 + off[0] * S, (fy + ty) / 2 + off[1] * S, a.get("label_color", col),
                   a.get("label_align", "center"))
    lg = spec.get("legend")
    if lg:                                               # swatch + name per curve, one row each
        x, y = ax.P(*(ax.e(v) for v in lg["at"]))
        for item in lg["items"]:
            item = {"curve": item} if isinstance(item, str) else item
            c = curves.get(item.get("curve"), {})
            col = item.get("color", c.get("color"))
            sw = line(x, y, x + 1.6 * S, y, stroke=col, w=c.get("width", 1) * CURVE_W)
            if c.get("dash"):
                sw = _dash(sw, "5 3")
            L.add("labels", sw, x, y - 0.1 * S, x + 1.6 * S, y + 0.1 * S)
            L.text(label_tb(item.get("text", c.get("label", ""))), x + 2.0 * S, y, item.get("text_color", col), "left")
            y += lg.get("row", 1.25) * S
    for t in spec.get("texts") or []:
        x, y = ax.P(*(ax.e(v) for v in t["at"]))
        tb = TextBlock(t["text"], t.get("kind", "sans"), t.get("size", 1) * style.ENZYME_SIZE)
        L.text(tb, x, y, t.get("color"), t.get("align", "center"))
    return L.graphic()


def check_plot(spec, tol=0.12):
    """Problems with one plot spec: expressions that fail, points that miss their curve."""
    ax = Axes(spec)
    curves = {c["id"]: c for c in spec.get("curves") or []}
    out = []
    for c in curves.values():
        pts = [q for q in ax.samples(c) if q]
        if len(pts) < 2:
            out.append(f"曲線 {c['id']} 算不出任何點")
    for p in spec.get("points") or []:
        on = p.get("on")
        if on is None:
            continue
        if on not in curves:
            out.append(f"點 {p.get('label', '')} 指向不存在的曲線 {on}")
            continue
        if "x" in p and "y" in p:
            x, y = ax.e(p["x"]), ax.e(p["y"])
            d = ax.distance(curves[on], x, y)
            if d > p.get("tol", tol):
                out.append(f"點 ({x:g}, {y:g}) {p.get('label', '')} 不在曲線 {on} 上（差 {d:.2f} S）")
    return out


def compose_plot(db, p, variant=None):
    """One plot, or several panels side by side (each with an optional title)."""
    spec_list = p.get("panels") or [p["plot"]]
    over = (p.get("variants") or {}).get(variant, {}) if variant else {}
    if over.get("panels"):
        spec_list = over["panels"]
    elif over.get("plot"):
        spec_list = [{**spec_list[0], **over["plot"]}]
    gs = []
    for spec in spec_list:
        g = draw_plot(spec)
        if spec.get("title"):
            tb = TextBlock(spec["title"], "sans", style.NAME_SIZE * 1.05)
            y = g.y0 - 0.7 * S - tb.bottom
            g = Graphic(g.body + tb.svg((g.x0 + g.x1) / 2, y), min(g.x0, (g.x0 + g.x1) / 2 - tb.width / 2),
                        y + tb.top, max(g.x1, (g.x0 + g.x1) / 2 + tb.width / 2), g.y1, {},
                        {**g.parts, "labels": _merge(g.parts.get("labels"), Graphic(tb.svg((g.x0 + g.x1) / 2, y),
                         (g.x0 + g.x1) / 2 - tb.width / 2, y + tb.top, (g.x0 + g.x1) / 2 + tb.width / 2, y + tb.bottom))})
        gs.append(g)
    if len(gs) == 1:
        return _caption(gs[0], p.get("caption"))
    gap = p.get("panel_gap", 3) * S
    body, parts, bb, x = [], {}, BBox(), 0.0
    top = min(g.y0 for g in gs)
    for i, g in enumerate(gs):
        m = g.moved(x - g.x0, 0)
        body.append(m.body)
        bb.add_g(m)
        parts[f"panel_{i + 1}"] = Graphic(m.body, m.x0, m.y0, m.x1, m.y1)
        for k, sub in m.parts.items():
            parts[k] = _merge(parts.get(k), sub)
        x = m.x1 + gap
    return _caption(Graphic("".join(body), bb.x0, min(bb.y0, top), bb.x1, bb.y1, {}, parts), p.get("caption"))


def _caption(g, cap):
    """caption: {text, color} — one line centred under the whole figure (its own part)."""
    if not cap:
        return g
    tb = TextBlock(cap["text"], "sans", style.NAME_SIZE * 1.05)
    y = g.y1 + 1.0 * S - tb.top
    svg = tb.svg(g.cx, y, cap.get("color"))
    cg = Graphic(svg, g.cx - tb.width / 2, y + tb.top, g.cx + tb.width / 2, y + tb.bottom)
    return Graphic(g.body + svg, min(g.x0, cg.x0), g.y0, max(g.x1, cg.x1), cg.y1, {},
                   {**g.parts, "labels": _merge(g.parts.get("labels"), cg)})


def _merge(a, b):
    if a is None:
        return b
    if b is None:
        return a
    return Graphic(a.body + b.body, min(a.x0, b.x0), min(a.y0, b.y0), max(a.x1, b.x1), max(a.y1, b.y1))
