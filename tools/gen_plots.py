"""Round 13: plots (enzyme kinetics, O2 binding, titration, Ramachandran) + peptide bond.

    python tools/gen_plots.py data

Plots are pathways with `layout: plot` (see sdb/draw/plot.py). Every labelled point that sits on
a curve gives both x and y, so `validate` checks it (K_m at V_{max}/2, intercepts −1/K_m and 1/V_{max},
P50 at half saturation, pK at half-equivalence …).
"""
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent))
from reaction_colors import COLORS  # noqa: E402

M, C = "magenta", "cyan"
RED, GREEN, ORANGE, PURPLE, GRAY = "red", "green", "orange", "purple", "gray"


# ---------------------------------------------------------------- enzyme kinetics

MM_PARAMS = {"Vmax": 1, "Km": 2}


def mm_axes(xmax=20, ymax=1.2, marks=True):
    return {
        "x": {"min": 0, "max": xmax, "label": "[S]",
              "marks": [{"at": "Km", "label": "K_m"}] if marks else []},
        "y": {"min": 0, "max": ymax, "label": "v_0",
              "marks": [{"at": "Vmax", "label": "V_{max}"}, {"at": "Vmax/2", "label": "V_{max}/2"}] if marks else []},
    }


MICHAELIS_MENTEN = {
    "size": [15, 9.5],
    "params": MM_PARAMS,
    **mm_axes(),
    "hlines": [{"y": "Vmax"}],
    "curves": [{"id": "mm", "y": "mm(x, Vmax, Km)", "width": 1.1}],
    "points": [{"on": "mm", "x": "Km", "y": "Vmax/2", "drop": True, "dot": True}],
    "texts": [
        {"at": [3.0, 0.16], "text": "First order (low [S])\nv_0 ≈ (V_{max}/K_m)[S]", "align": "left", "color": C},
        {"at": [15.5, 0.7], "text": "Zero order (saturated)\nv_0 ≈ V_{max}", "color": M},
        {"at": [12.5, 0.38], "text": "v_0 = V_{max}[S] / (K_m + [S])", "size": 1.05},
    ],
}

LB_PARAMS = {"Vmax": 1, "Km": 2}
LINEWEAVER_BURK = {
    "size": [15, 9.5],
    "params": LB_PARAMS,
    "x": {"min": -0.8, "max": 1.6, "label": "1/[S]",
          "marks": [{"at": "-1/Km", "label": "−1/K_m"}]},
    "y": {"min": -0.6, "max": 4.4, "label": "1/v_0",
          "marks": [{"at": "1/Vmax", "label": "1/V_{max}", "dy": -0.6}]},
    "curves": [{"id": "lb", "y": "Km/Vmax*x + 1/Vmax", "from": "-1/Km*1.25", "to": 1.6, "width": 1.1,
                "label": "slope = K_m/V_{max}", "label_at": 1.15, "label_off": [0.3, 0.9]}],
    "points": [{"on": "lb", "x": "-1/Km", "y": 0, "dot": True},
               {"on": "lb", "x": 0, "y": "1/Vmax", "dot": True}],
    "texts": [{"at": [1.0, 0.7], "text": "1/v_0 = (K_m/V_{max})(1/[S]) + 1/V_{max}", "size": 1.05}],
}


def inhibition(kind):
    """MM panel + Lineweaver–Burk panel: no inhibitor (black) vs + inhibitor (red)."""
    # apparent Km and Vmax with the inhibitor present
    km_i, vm_i = {"competitive": ("2.5*Km", "Vmax"), "noncompetitive": ("Km", "Vmax/2.2"),
                  "uncompetitive": ("Km/3", "Vmax/3")}[kind]
    params = {"Vmax": 1, "Km": 2, "KmI": km_i, "VmI": vm_i}
    hl = [{"y": "Vmax", "to": 11}] + ([{"y": "VmI", "color": RED, "to": 11}] if kind != "competitive" else [])
    xm = {"competitive": [{"at": "Km", "label": "K_m"}, {"at": "KmI", "label": "K_m^{app}", "color": RED}],
          "noncompetitive": [{"at": "Km", "label": "K_m (unchanged)", "dx": 1.2}],
          "uncompetitive": [{"at": "KmI", "label": "K_m^{app}", "color": RED, "dx": -0.6},
                            {"at": "Km", "label": "K_m", "dx": 0.5}]}[kind]
    ym = [{"at": "Vmax", "label": "V_{max}"}] + \
        ([{"at": "VmI", "label": "V_{max}^{app}", "color": RED}] if kind != "competitive" else [])
    inh_lab = {"competitive": (16, [-1.0, 0.85]), "noncompetitive": (16, [-1.0, -0.65]),
               "uncompetitive": (16, [-1.0, -0.65])}[kind]
    mmp = {
        "title": "Michaelis–Menten",
        "size": [12, 8.5],
        "params": params,
        "x": {"min": 0, "max": 20, "label": "[S]", "marks": xm},
        "y": {"min": 0, "max": 1.2, "label": "v_0", "marks": ym},
        "hlines": hl,
        "curves": [
            {"id": "no", "y": "mm(x, Vmax, Km)", "label": "No inhibitor", "label_at": 16, "label_off": [-1.0, -0.7]},
            {"id": "inh", "y": "mm(x, VmI, KmI)", "color": RED, "label": "+ Inhibitor", "label_at": inh_lab[0],
             "label_off": inh_lab[1]},
        ],
        "points": [
            {"on": "no", "x": "Km", "y": "Vmax/2", "drop": "x", "dot": True},
            {"on": "inh", "x": "KmI", "y": "VmI/2", "drop": "x", "dot": True, "color": RED},
        ],
    }
    lbx = {"competitive": -0.8, "noncompetitive": -0.8, "uncompetitive": -1.9}[kind]
    ymax = {"competitive": 6, "noncompetitive": 6, "uncompetitive": 4}[kind]
    red_lab = {"competitive": (0.85, [-0.4, 0], "right"), "noncompetitive": (0.5, [0.5, 0.4], "left"),
               "uncompetitive": (-0.45, [-0.4, -0.8], "right")}[kind]
    lbp = {
        "title": "Lineweaver–Burk",
        "size": [12, 8.5],
        "params": params,
        "x": {"min": lbx, "max": 1.4, "label": "1/[S]", "marks": []},
        "y": {"min": -0.6, "max": ymax, "label": "1/v_0", "marks": []},
        "curves": [
            {"id": "no", "y": "Km/Vmax*x + 1/Vmax", "from": "-1/Km*1.3", "to": 1.4, "label": "No inhibitor",
             "label_at": 1.15, "label_off": [0.3, 0.8]},
            {"id": "inh", "y": "KmI/VmI*x + 1/VmI", "from": "-1/KmI*1.3", "to": 1.4, "color": RED,
             "label": "+ Inhibitor", "label_at": red_lab[0], "label_off": red_lab[1], "label_align": red_lab[2]},
        ],
        "points": [
            {"on": "no", "x": "-1/Km", "y": 0, "dot": True},
            {"on": "no", "x": 0, "y": "1/Vmax", "dot": True},
            {"on": "inh", "x": "-1/KmI", "y": 0, "dot": True, "color": RED},
            {"on": "inh", "x": 0, "y": "1/VmI", "dot": True, "color": RED},
        ],
    }
    if kind == "competitive":
        lbp["x"]["marks"] = [{"at": "-1/Km", "label": "−1/K_m", "dy": 1.0},
                             {"at": "-1/KmI", "label": "−1/K_m^{app}", "color": RED, "dx": -0.9}]
        lbp["y"]["marks"] = [{"at": "1/Vmax", "label": "1/V_{max} (same)", "dy": -0.75}]
    elif kind == "noncompetitive":
        lbp["x"]["marks"] = [{"at": "-1/Km", "label": "−1/K_m (same)", "dx": -1.8, "dy": 0.6}]
        lbp["y"]["marks"] = [{"at": "1/Vmax", "label": "1/V_{max}", "right": True, "dx": 0.2, "dy": 0.6},
                             {"at": "1/VmI", "label": "1/V_{max}^{app}", "color": RED, "dy": -0.35}]
    else:
        lbp["x"]["marks"] = [{"at": "-1/Km", "label": "−1/K_m", "dx": 0.9},
                             {"at": "-1/KmI", "label": "−1/K_m^{app}", "color": RED, "dy": 0.9}]
        lbp["y"]["marks"] = [{"at": "1/Vmax", "label": "1/V_{max}", "right": True, "dx": 0.2, "dy": 0.6},
                             {"at": "1/VmI", "label": "1/V_{max}^{app}", "color": RED, "dy": -0.45}]
    return [mmp, lbp]


INHIBITION_NOTES = {
    "competitive": ("Competitive inhibition", "K_m ↑   V_{max} unchanged",
                    "競爭性抑制：抑制劑和受質搶同一個活性位（只結合 E）。K_m^app = αK_m 變大、V_{max} 不變；"
                    "加大量受質可以克服。雙倒數圖：各直線交於 y 軸（1/V_{max} 相同），斜率變大。"
                    "例：methotrexate（DHFR）、statin（HMG-CoA reductase）、malonate（SDH）、ethanol 治 methanol 中毒。"),
    "noncompetitive": ("Noncompetitive inhibition (pure)", "K_m unchanged   V_{max} ↓",
                       "（純）非競爭性抑制：抑制劑結合在別處，對 E 與 ES 親和力相同（α = α'）。"
                       "V_{max}^app = V_{max}/α 變小、K_m 不變；加受質無法克服。雙倒數圖：直線交於 x 軸（−1/K_m 相同），y 截距變大。"
                       "例：重金屬、許多別位抑制劑。考題的 noncompetitive 通常指這種。"),
    "uncompetitive": ("Uncompetitive inhibition", "K_m ↓   V_{max} ↓ (same factor) → parallel lines",
                      "反競爭性抑制：抑制劑只結合 ES 複合物。V_{max} 與 K_m 以同樣倍數下降（V_{max}/K_m 比值不變）"
                      "→ 雙倒數圖為平行線（斜率 K_m/V_{max} 相同）。例：lithium 抑制 inositol monophosphatase（口訣：parallel）。"),
}


def inhibition_pathway(kind):
    title, change, notes = INHIBITION_NOTES[kind]
    panels = inhibition(kind)
    return {"id": f"inhibition_{kind}", "title": title, "layout": "plot", "panels": panels,
            "caption": {"text": f"{title}:  {change}", "color": RED},
            "tags": ["enzyme_kinetics", "inhibition"], "notes": notes}


ALLOSTERIC = {
    "size": [15, 10],
    "params": {"Vmax": 1, "K": 6, "n": 2.8},
    "x": {"min": 0, "max": 24, "label": "[S]", "marks": [{"at": "K", "label": "K_{0.5}"}]},
    "y": {"min": 0, "max": 1.32, "label": "v_0", "marks": [{"at": "Vmax", "label": "V_{max}"},
                                                           {"at": "Vmax/2", "label": "V_{max}/2"}]},
    "hlines": [{"y": "Vmax"}],
    "curves": [
        {"id": "act", "y": "Vmax*hill(x, K/2.2, n)", "color": C},
        {"id": "none", "y": "Vmax*hill(x, K, n)"},
        {"id": "inh", "y": "Vmax*hill(x, K*2, n)", "color": RED},
        {"id": "mm", "y": "mm(x, Vmax, K)", "color": GRAY, "dash": True},
    ],
    "points": [{"on": "none", "x": "K", "y": "Vmax/2", "drop": "x", "dot": True}],
    "legend": {"at": [25.5, 1.05], "items": [
        {"curve": "act", "text": "+ Activator (K_{0.5} ↓)"}, {"curve": "none", "text": "No effector (sigmoid)"},
        {"curve": "inh", "text": "+ Inhibitor (K_{0.5} ↑)"}, {"curve": "mm", "text": "Hyperbolic (M–M), same K"}]},
}

PH_TEMP = [
    {"title": "pH optimum", "size": [12, 7.5],
     "x": {"min": 0, "max": 12, "label": "pH", "ticks": [0, 2, 4, 6, 8, 10, 12]},
     "y": {"min": 0, "max": 1.25, "label": "Activity"},
     "curves": [
         {"id": "pepsin", "y": "bell(x, 1.8, 0.9)", "color": M, "label": "Pepsin", "label_at": 1.8, "label_off": [0, -0.7],
          "label_align": "center"},
         {"id": "amylase", "y": "bell(x, 6.8, 1.0)", "color": GREEN, "label": "Salivary amylase", "label_at": 5.6,
          "label_off": [-0.3, -0.3], "label_align": "right"},
         {"id": "trypsin", "y": "bell(x, 8.0, 1.0)", "color": C, "label": "Trypsin", "label_at": 8, "label_off": [0, -0.7], "label_align": "center"},
         {"id": "alp", "y": "bell(x, 10.0, 0.9)", "color": PURPLE, "label": "Alkaline phosphatase", "label_at": 10.9,
          "label_off": [0.3, -0.2]},
     ],
     "points": [{"on": "pepsin", "x": 1.8, "y": 1, "dot": False}]},
    {"title": "Temperature", "size": [12, 7.5],
     "params": {"T0": 37},
     "x": {"min": 0, "max": 70, "label": "Temperature (°C)", "ticks": [0, 10, 20, 30, 40, 50, 60, 70]},
     "y": {"min": 0, "max": 1.25, "label": "Activity"},
     "curves": [{"id": "t", "y": "2**((x-T0)/10) / (1 + exp((x-44)/2.6)) / 0.935", }],
     "arrows": [{"from": [55, 0.95], "to": [50.5, 0.55], "color": RED, "label": "Denaturation", "label_off": [1.2, -1.0]},
                {"from": [8, 0.62], "to": [24, 0.5], "color": C, "label": "Rate ×2 per 10 °C (Q_{10})",
                 "label_off": [-1.0, -0.9]}],
     "points": [{"on": "t", "x": 37, "drop": "x", "label": "≈ 37 °C", "label_off": [-0.4, -0.75], "label_align": "right"}]},
]


# ---------------------------------------------------------------- O2 binding

O2P = {"P50m": 2.8, "P50": 26, "n": 2.8}
O2_AXES = {
    "x": {"min": 0, "max": 110, "label": "pO_2 (mmHg)", "ticks": [0, 20, 40, 60, 80, 100]},
    "y": {"min": 0, "max": 1.12, "label": "Fractional saturation (Y)", "label_pos": "center",
          "ticks": [0, 0.5, 1]},
}

O2_MB_HB = {
    "size": [16, 10],
    "params": O2P,
    **O2_AXES,
    "bands": [{"x0": 20, "x1": 40, "label": "Tissues", "color": C},
              {"x0": 90, "x1": 105, "label": "Lungs", "color": M}],
    "curves": [
        {"id": "mb", "y": "hill(x, P50m, 1)", "color": C, "label": "Myoglobin (hyperbolic, n_H = 1)",
         "label_at": 60, "label_off": [-1.2, -0.75], "label_align": "left"},
        {"id": "hb", "y": "hill(x, P50, n)", "color": RED, "label": "Hemoglobin (sigmoid, n_H ≈ 2.8)",
         "label_at": 50, "label_off": [0.5, 0.75], "label_align": "left"},
    ],
    "points": [
        {"on": "mb", "x": "P50m", "y": 0.5, "drop": "x", "color": C, "label": "P_{50} ≈ 2.8", "label_off": [0.4, 1.5]},
        {"on": "hb", "x": "P50", "y": 0.5, "drop": "both", "color": RED, "label": "P_{50} ≈ 26", "label_off": [0.4, 0.6]},
    ],
}


def hb_curves(extra):
    return [{"id": "normal", "y": "hill(x, P50, n)", "label": "Normal (P_{50} 26)", "label_at": 33,
             "label_off": [0.4, 0.15]}] + extra


HB_SHIFT = {
    "size": [20, 10],
    "params": {**O2P, "Pl": 17, "Pr": 36},
    "x": {"min": 0, "max": 150, "label": "pO_2 (mmHg)", "ticks": [0, 20, 40, 60, 80, 100],
          "marks": [{"at": "Pl", "label": "17", "color": C, "dy": 0.9}, {"at": "P50", "label": "26 (normal)", "dy": 1.8},
                    {"at": "Pr", "label": "36", "color": RED, "dy": 0.9}]},
    "y": O2_AXES["y"],
    "curves": [{"id": "left", "y": "hill(x, Pl, n)", "color": C}, {"id": "normal", "y": "hill(x, P50, n)"},
               {"id": "right", "y": "hill(x, Pr, n)", "color": RED}],
    "points": [{"on": "left", "x": "Pl", "y": 0.5, "color": C, "drop": "x"},
               {"on": "normal", "x": "P50", "y": 0.5, "drop": "x"},
               {"on": "right", "x": "Pr", "y": 0.5, "color": RED, "drop": "x"}],
    "arrows": [{"from": [21, 0.38], "to": [11, 0.38], "color": C}, {"from": [31, 0.38], "to": [42, 0.38], "color": RED}],
    "texts": [
        {"at": [56, 0.4], "align": "left", "color": RED,
         "text": "Right shift (P_{50} ↑)\naffinity ↓, unloading ↑\n↑ H^+ (↓ pH, Bohr)\n↑ CO_2\n↑ 2,3-BPG\n↑ Temperature"},
        {"at": [118, 0.4], "align": "left", "color": C,
         "text": "Left shift (P_{50} ↓)\nO_2 affinity ↑\nHbF, CO, MetHb\n↓ H^+, ↓ CO_2\n↓ 2,3-BPG\n↓ Temperature"},
    ],
}

BPG = {
    "size": [18, 10],
    "params": {"n": 2.8, "P0": 8, "P5": 26, "P8": 31},
    "x": {"min": 0, "max": 130, "label": "pO_2 (mmHg)", "ticks": [0, 20, 40, 60, 80, 100, 120]},
    "y": O2_AXES["y"],
    "bands": [{"x0": 20, "x1": 40, "label": "Tissues", "color": C}, {"x0": 90, "x1": 105, "label": "Lungs", "color": M}],
    "curves": [
        {"id": "stripped", "y": "hill(x, P0, 2.4)", "color": C},
        {"id": "normal", "y": "hill(x, P5, n)"},
        {"id": "altitude", "y": "hill(x, P8, n)", "color": RED},
    ],
    "points": [{"on": "stripped", "x": "P0", "y": 0.5, "color": C, "drop": "x"},
               {"on": "normal", "x": "P5", "y": 0.5, "drop": "x"},
               {"on": "altitude", "x": "P8", "y": 0.5, "color": RED, "drop": "x"}],
    "legend": {"at": [46, 0.62], "items": [
        {"curve": "stripped", "text": "Stripped Hb, no BPG (P_{50} ≈ 8)"},
        {"curve": "normal", "text": "Normal, 5 mM BPG (P_{50} ≈ 26)"},
        {"curve": "altitude", "text": "High altitude, 8 mM BPG (P_{50} ≈ 31)"}]},
}

CO_CURVES = {
    "size": [18, 10],
    "params": {"n": 2.8, "P50": 26},
    "x": {"min": 0, "max": 130, "label": "pO_2 (mmHg)", "ticks": [0, 20, 40, 60, 80, 100, 120]},
    "y": {"min": 0, "max": 1.12, "label": "O_2 content (× normal capacity)", "label_pos": "center",
          "ticks": [0, 0.5, 1]},
    "curves": [
        {"id": "normal", "y": "hill(x, P50, n)"},
        {"id": "anemia", "y": "0.5*hill(x, P50, n)", "color": C},
        {"id": "co", "y": "0.5*hill(x, 10, 1.6)", "color": RED},
    ],
    "points": [{"on": "anemia", "x": "P50", "y": 0.25, "color": C, "dot": True},
               {"on": "co", "x": 10, "y": 0.25, "color": RED, "dot": True}],
    "legend": {"at": [52, 0.36], "items": [
        {"curve": "normal", "text": "Normal Hb"},
        {"curve": "anemia", "text": "Anemia, 50 % Hb (same P_{50})"},
        {"curve": "co", "text": "50 % COHb (capacity ↓ + left shift)"}]},
}

HILL_PLOT = {
    "size": [15, 10],
    "params": {"P50m": 2.8, "P50": 26, "KR": 1.5, "KT": 250, "L": 8.5e4},   # MWC fit: P50 ≈ 26, n_H ≈ 3
    "x": {"min": -1, "max": 3, "label": "log pO_2", "ticks": [-1, 0, 1, 2, 3], "at": -4},
    "y": {"min": -4, "max": 4, "label": "log [Y/(1−Y)]", "ticks": [-4, -2, 0, 2, 4], "at": -1},
    "hlines": [{"y": 0, "label": "Y = 0.5"}],
    "curves": [
        {"id": "mb", "y": "x - log10(P50m)", "color": C},
        {"id": "hb", "y": "logit(mwc(10**x, KR, KT, L))", "color": RED},
        {"id": "r", "y": "x - log10(KR)", "color": GRAY, "dash": True, "from": 1.6, "to": 3},
        {"id": "t", "y": "x - log10(KT)", "color": GRAY, "dash": True, "from": -1, "to": 1.9},
    ],
    "points": [{"on": "mb", "x": "log10(P50m)", "y": 0, "color": C},
               {"on": "hb", "x": "log10(P50)", "y": 0, "color": RED, "tol": 0.2}],
    "legend": {"at": [1.3, -1.75], "items": [
        {"curve": "mb", "text": "Mb: slope n_H = 1"}, {"curve": "hb", "text": "Hb: slope n_H ≈ 3 at Y = 0.5"},
        {"curve": "r", "text": "slope 1 asymptotes (all-T / all-R)"}]},
    "texts": [{"at": [2.3, 2.75], "text": "R state", "size": 0.9, "color": GRAY},
              {"at": [0.0, -2.9], "text": "T state", "size": 0.9, "color": GRAY}],
}


# ---------------------------------------------------------------- titration / Ramachandran

TITRATION_GLY = {
    "size": [14, 10],
    "params": {"pK1": 2.34, "pK2": 9.60, "pI": "(pK1+pK2)/2"},
    "x": {"min": 0, "max": 2.15, "label": "OH^- added (equivalents)", "ticks": [0, 0.5, 1, 1.5, 2]},
    "y": {"min": 0, "max": 13, "label": "pH", "ticks": [0, 3, 6, 9, 12]},
    "curves": [{"id": "gly", "t": [0.6, 12.6], "x": "1/(1+10**(pK1-t)) + 1/(1+10**(pK2-t))", "y": "t", "width": 1.1}],
    "points": [
        {"on": "gly", "x": 0.5, "y": "pK1", "drop": "y", "label": "pK_1 = 2.34 (α-COOH)", "label_off": [0.4, 0.7]},
        {"on": "gly", "x": 1.0, "y": "pI", "drop": "y", "color": M, "label": "pI = (pK_1 + pK_2)/2 = 5.97",
         "label_off": [0.4, 0.1]},
        {"on": "gly", "x": 1.5, "y": "pK2", "drop": "y", "label": "pK_2 = 9.60 (α-NH_3^+)", "label_off": [-0.4, -0.65],
         "label_align": "right"},
    ],
    "texts": [
        {"at": [0.05, 4.6], "text": "^+H_3N–CH_2–COOH", "align": "left", "kind": "serif", "size": 0.8},
        {"at": [0.9, 7.7], "text": "^+H_3N–CH_2–COO^-\n(zwitterion, net 0)", "align": "right", "kind": "serif", "size": 0.8},
        {"at": [1.97, 12.25], "text": "H_2N–CH_2–COO^-", "align": "right", "kind": "serif", "size": 0.8},
        {"at": [1.6, 4.0], "text": "Buffer regions:\npK ± 1", "color": C, "size": 0.9},
    ],
}

# Ramachandran: allowed (dark) and partially allowed regions, schematic after Lehninger Fig. 4-8
RAMA = {
    "size": [12, 12],
    "x": {"min": -180, "max": 180, "label": "φ (degrees)", "ticks": [-180, -90, 0, 90, 180], "at": -180,
          "label_pos": "center"},
    "y": {"min": -180, "max": 180, "label": "ψ (degrees)", "ticks": [-180, -90, 0, 90, 180], "at": -180,
          "label_pos": "center"},
    "areas": [
        {"points": [[-180, 180], [-45, 180], [-45, 100], [-60, 70], [-100, 60], [-150, 90], [-180, 90]], "color": C,
         "label": "β", "label_at": [-165, 165], "opacity": 0.3},
        {"points": [[-180, -180], [-45, -180], [-45, -170], [-180, -165]], "color": C, "opacity": 0.3},
        {"points": [[-160, -70], [-110, -70], [-45, -70], [-35, -40], [-45, -15], [-90, 0], [-160, -20]],
         "color": M, "label": "α_R", "label_at": [-135, -40], "opacity": 0.3},
        {"points": [[45, 20], [75, 25], [85, 60], [65, 85], [45, 70]], "color": GREEN, "label": "α_L",
         "label_at": [100, 55], "opacity": 0.3},
    ],
    "points": [
        {"x": -51, "y": 153, "color": ORANGE, "label": "Collagen triple helix", "label_off": [0.5, -0.3], "leader": True},
        {"x": -139, "y": 135, "color": C, "label": "Antiparallel β sheet", "label_off": [3.5, 0.15], "leader": True},
        {"x": -119, "y": 113, "color": C, "label": "Parallel β sheet", "label_off": [2.9, 0.95], "leader": True},
        {"x": -57, "y": -47, "color": M, "label": "Right-handed α helix", "label_off": [1.3, 0.2], "leader": True},
        {"x": 57, "y": 47, "color": GREEN, "label": "Left-handed α helix", "label_off": [0.6, 1.5], "leader": True},
    ],
    "texts": [{"at": [95, -125], "text": "Gly (no side chain)\ncan go almost anywhere;\nPro: φ fixed ≈ −60°", "size": 0.85,
               "color": GRAY}],
}


# ---------------------------------------------------------------- peptide molecules

def physio_tag(tags):
    return tags


MOLECULES = [
    {
        "id": "glycylalanine", "name": "Glycylalanine", "label": "Gly-Ala",
        "aliases": ["Gly-Ala", "GA dipeptide", "雙肽"],
        "smiles": "[NH3+]CC(=O)N[C@@H](C)C(=O)[O-]",
        "free": {
            "atoms": {
                "n1": ["^+H_3N", -3.3, 0],
                "ca1": ["CH_2", -2.2, 0.55],
                "c1": ["C", -1.1, 0],
                "o1": ["O", -1.1, -1.05],
                "n2": ["N", 0, 0.55],
                "h2": ["H", 0, 1.55],
                "ca2": ["C", 1.1, 0],
                "ha": ["H", 1.1, -1.0],
                "cb": ["CH_3", 1.1, 1.05],
                "c2": ["COO^-", 2.4, 0.55],
            },
            "bonds": [["n1", "ca1"], ["ca1", "c1"], ["c1", "o1", 2], ["c1", "n2"], ["n2", "h2"], ["n2", "ca2"],
                      ["ca2", "ha"], ["ca2", "cb"], ["ca2", "c2"]],
            "notes": [["Gly residue", -2.4, 1.65, "gray"], ["Ala residue", 1.5, 2.1, "gray"],
                      ["N-terminus", -3.3, -0.85, "gray"], ["C-terminus", 3.1, -0.35, "gray"]],
        },
        "decorate": {"boxes": [{"atoms": ["c1", "o1", "n2", "h2"], "name": "peptide", "color": "salmon"}]},
        "tags": ["peptide"],
        "notes": "Gly-Ala 雙肽：胜肽依慣例由 N 端（左）寫到 C 端（右）。底色框為胜肽鍵（–CO–NH–，醯胺鍵）。"
                 "n 個胺基酸的直鏈胜肽有 n − 1 個胜肽鍵、脫去 n − 1 個 H₂O。",
    },
    {
        "id": "peptide_unit", "name": "Peptide bond", "label": "Peptide bond",
        "aliases": ["peptide bond", "amide plane", "胜肽鍵", "肽鍵"],
        "smiles": "CC(=O)NC",
        "free": {
            "schematic": True,
            "atoms": {
                "ca1": ["C_α", -1.1, 0.55],
                "c": ["C", 0, 0],
                "o": ["O", 0, -1.05],
                "n": ["N", 1.1, 0.55],
                "h": ["H", 1.1, 1.55],
                "ca2": ["C_α", 2.2, 0],
            },
            "bonds": [["ca1", "c"], ["c", "o", 2], ["c", "n"], ["n", "h"], ["n", "ca2"]],
        },
        "decorate": {"colors": {"o": "red", "n": "cyan"}},
        "tags": ["peptide"],
        "notes": "胜肽鍵的一般寫法（C=O 雙鍵、C–N 單鍵）。真實結構是兩個共振式的混成，見 peptide_bond_resonance。",
    },
    {
        "id": "peptide_unit_resonance", "name": "Peptide bond (resonance form)", "label": "Resonance form",
        "aliases": ["peptide bond resonance", "胜肽鍵共振式"],
        "smiles": "C/C([O-])=[NH+]/C",
        "free": {
            "schematic": True,
            "atoms": {
                "ca1": ["C_α", -1.1, 0.55],
                "c": ["C", 0, 0],
                "o": ["O^-", 0, -1.05],
                "n": ["N^+", 1.1, 0.55],
                "h": ["H", 1.1, 1.55],
                "ca2": ["C_α", 2.2, 1.1],
            },
            "bonds": [["ca1", "c"], ["c", "o"], ["c", "n", 2], ["n", "h"], ["n", "ca2"]],
        },
        "decorate": {"colors": {"o": "red", "n": "cyan"}},
        "tags": ["peptide"],
        "notes": "胜肽鍵的共振式：N 的孤對電子推向 C=O，使 C–N 帶部分雙鍵性質（約 40 %）。",
    },
    {
        "id": "peptide_backbone", "name": "Polypeptide backbone (φ, ψ)", "label": "Peptide planes (φ, ψ)",
        "aliases": ["phi psi", "dihedral angle", "二面角", "peptide plane"],
        "smiles": "CC(=O)NC(C)C(=O)NC",
        "free": {
            "schematic": True,
            "atoms": {
                "ca1": ["C_α", -3.3, 0.55],
                "c1": ["C", -2.2, 0],
                "o1": ["O", -2.2, -1.05],
                "n1": ["N", -1.1, 0.55],
                "h1": ["H", -1.1, 1.55],
                "ca2": ["C_α", 0, 0],
                "r": ["R", 0, -1.05],
                "c2": ["C", 1.1, 0.55],
                "o2": ["O", 1.1, 1.6],
                "n2": ["N", 2.2, 0],
                "h2": ["H", 2.2, -1.0],
                "ca3": ["C_α", 3.3, 0.55],
            },
            "bonds": [["ca1", "c1"], ["c1", "o1", 2], ["c1", "n1"], ["n1", "h1"], ["n1", "ca2"], ["ca2", "r"],
                      ["ca2", "c2"], ["c2", "o2", 2], ["c2", "n2"], ["n2", "h2"], ["n2", "ca3"]],
            "notes": [["φ", -0.62, -0.05, "magenta"], ["ψ", 0.62, -0.05, "cyan"],
                      ["Peptide plane", -2.2, 2.45, "gray"], ["Peptide plane", 2.75, -1.95, "gray"],
                      ["C–N 1.32 Å (partial double bond)", 0, 3.3, "gray"]],
        },
        "decorate": {"boxes": [{"atoms": ["ca1", "c1", "o1", "n1", "h1"], "name": "plane", "color": "#E3F1FA"},
                               {"atoms": ["c2", "o2", "n2", "h2", "ca3"], "name": "plane", "color": "#E3F1FA"}]},
        "tags": ["peptide", "protein_structure"],
        "notes": "多肽主鏈：每個胜肽鍵連同兩側 Cα 共平面（剛性、幾乎都是 trans）。只有 Cα 兩側的單鍵可旋轉："
                 "φ（N–Cα）與 ψ（Cα–C）。允許的 φ/ψ 組合見 Ramachandran plot。"
                 "C–N 1.32 Å，介於單鍵 1.49 Å 與雙鍵 1.27 Å 之間。",
    },
    {
        "id": "glutathione", "name": "Glutathione (GSH)", "label": "Glutathione (GSH)",
        "aliases": ["GSH", "γ-Glu-Cys-Gly", "穀胱甘肽", "reduced glutathione"],
        "smiles": "[NH3+][C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(=O)[O-])C(=O)[O-]",
        "free": {
            "atoms": {
                "a1": ["COO^-", -4.4, -0.55],
                "ca1": ["CH", -4.4, 0.55],
                "n1": ["^+H_3N", -5.6, 0.55],
                "cb1": ["CH_2", -3.3, 1.1],
                "cg1": ["CH_2", -2.2, 0.55],
                "cd1": ["C", -1.1, 1.1],
                "od1": ["O", -1.1, 2.15],
                "n2": ["NH", 0, 0.55],
                "ca2": ["CH", 1.1, 1.1],
                "cb2": ["CH_2", 1.1, 2.15],
                "sg": ["SH", 1.1, 3.2],
                "c2": ["C", 2.2, 0.55],
                "o2": ["O", 2.2, -0.5],
                "n3": ["NH", 3.3, 1.1],
                "ca3": ["CH_2", 4.4, 0.55],
                "c3": ["COO^-", 5.6, 1.1],
            },
            "bonds": [["a1", "ca1"], ["ca1", "n1"], ["ca1", "cb1"], ["cb1", "cg1"], ["cg1", "cd1"], ["cd1", "od1", 2],
                      ["cd1", "n2"], ["n2", "ca2"], ["ca2", "cb2"], ["cb2", "sg"], ["ca2", "c2"], ["c2", "o2", 2],
                      ["c2", "n3"], ["n3", "ca3"], ["ca3", "c3"]],
            "notes": [["γ-Glu", -3.9, 2.4, "gray"], ["Cys", 2.3, 2.9, "gray"], ["Gly", 4.6, 2.2, "gray"],
                      ["γ-peptide bond", -1.0, -0.55, "magenta"]],
        },
        "decorate": {"colors": {"sg": "orange"},
                     "boxes": [{"atoms": ["cd1", "od1", "n2"], "name": "gamma", "color": "#F8D3E6"}]},
        "tags": ["peptide"],
        "notes": "GSH：γ-Glu-Cys-Gly 三胜肽。Glu 以側鏈 γ-COOH（不是 α-COOH）接 Cys → γ-胜肽鍵，不被一般 protease 切。"
                 "–SH 是還原劑：2 GSH → GSSG（glutathione peroxidase 清除 H₂O₂，含 Se）；"
                 "glutathione reductase 用 NADPH（來自 PPP，G6PD 缺乏 → 溶血）把 GSSG 還原回 GSH。",
    },
]


REACTIONS = [
    {
        "id": "dipeptide_formation", "title": "Peptide bond formation (condensation)",
        "substrates": ["glycine", "alanine"], "products": ["glycylalanine"], "cofactors_out": ["H2O"],
        "enzyme": "Ribosome (peptidyl transferase)",
        "reversible": False,
        "card": {
            "type": "縮合（脫水）反應：前一個胺基酸的 α-COO⁻ 與下一個的 α-NH₃⁺ 形成醯胺鍵，脫去 1 H₂O",
            "deltaG": "直接縮合在熱力學上不利（ΔG > 0）；細胞內先把胺基酸活化成 aminoacyl-tRNA（消耗 ATP），由核糖體 peptidyl transferase（23S/28S rRNA，ribozyme）催化",
            "mechanism": "胜肽鍵的逆反應（水解）雖然有利但很慢，所以蛋白質在水中穩定；酸、鹼或 protease 可水解",
            "exam_points": ["n 個胺基酸的直鏈胜肽：n − 1 個胜肽鍵、脫 n − 1 個 H₂O；環狀胜肽則為 n 個",
                            "胜肽鍵 = C(=O)–N 鍵，部分雙鍵性質、平面、幾乎都是 trans（Pro 前可為 cis）",
                            "寫法：N 端在左、C 端在右"],
        },
        "tags": ["peptide"],
    },
]


PATHWAYS = [
    {"id": "michaelis_menten", "title": "Michaelis–Menten kinetics", "layout": "plot", "plot": MICHAELIS_MENTEN,
     "tags": ["enzyme_kinetics"],
     "notes": "Michaelis–Menten：v₀ = V_{max}[S]/(K_m + [S])。K_m 是 v₀ = ½V_{max} 時的 [S]（K_m 小 → 親和力高）。"
              "[S] ≪ K_m 時近似一級反應（v₀ ∝ [S]），[S] ≫ K_m 時近似零級（v₀ ≈ V_{max}，酵素飽和）。"
              "V_{max} = k_cat[E]_t；催化效率 = k_cat/K_m。"},
    {"id": "lineweaver_burk", "title": "Lineweaver–Burk (double reciprocal) plot", "layout": "plot",
     "plot": LINEWEAVER_BURK, "tags": ["enzyme_kinetics"],
     "notes": "雙倒數圖：1/v₀ = (K_m/V_{max})(1/[S]) + 1/V_{max}。y 截距 = 1/V_{max}，x 截距 = −1/K_m，斜率 = K_m/V_{max}。"
              "把雙曲線拉成直線，方便讀出 K_m、V_{max} 和判斷抑制型式。"},
    inhibition_pathway("competitive"),
    inhibition_pathway("noncompetitive"),
    inhibition_pathway("uncompetitive"),
    {"id": "allosteric_kinetics", "title": "Allosteric enzyme kinetics (sigmoid)", "layout": "plot",
     "plot": ALLOSTERIC, "tags": ["enzyme_kinetics", "allostery"],
     "notes": "別位酵素（如 ATCase、PFK-1）的 v₀–[S] 是 S 形（協同性），不遵守 M–M；用 K_{0.5} 代替 K_m。"
              "正向調節劑（activator）讓曲線左移、K_{0.5} 下降；負向調節劑（inhibitor）讓曲線右移。"
              "K 系統改變 K_{0.5}（本圖），V 系統改變 V_{max}。灰色虛線為同 K 的雙曲線作對照。"},
    {"id": "enzyme_ph_temperature", "title": "Enzyme activity vs pH and temperature", "layout": "plot",
     "panels": PH_TEMP, "tags": ["enzyme_kinetics"],
     "notes": "pH：每個酵素有最適 pH（pepsin ≈ 1.5–2、唾液 amylase ≈ 6.8、trypsin ≈ 8、alkaline phosphatase ≈ 10），"
              "反映活性位殘基的解離狀態。溫度：先依 Q₁₀ ≈ 2 上升，過最適溫度（人體約 37 °C）後因變性快速下降。曲線為示意。"},
    {"id": "o2_binding_mb_hb", "title": "O2 binding: myoglobin vs hemoglobin", "layout": "plot",
     "plot": O2_MB_HB, "tags": ["hemoglobin"],
     "notes": "Mb：單一 heme，雙曲線，P₅₀ ≈ 2.8 mmHg（親和力高，儲氧）。Hb：四聚體、協同結合，S 形曲線，P₅₀ ≈ 26 mmHg，"
              "Hill 係數 ≈ 2.8。肺（~100 mmHg）幾乎飽和、組織（20–40 mmHg）落在曲線陡峭段 → 可大量卸氧；Mb 在組織 pO₂ 仍幾乎飽和，所以不能當運輸者。"},
    {"id": "hb_oxygen_shift", "title": "Hb O2 curve: right and left shifts", "layout": "plot",
     "plot": HB_SHIFT, "tags": ["hemoglobin"],
     "notes": "右移（P₅₀ ↑、親和力 ↓、組織卸氧 ↑）：H⁺ ↑（pH ↓，Bohr effect）、CO₂ ↑（carbamino + H⁺）、2,3-BPG ↑、溫度 ↑ —— 運動中的肌肉。"
              "左移（親和力 ↑）：HbF（γ 鏈不易結合 BPG）、CO、MetHb、pH ↑、低溫、BPG ↓（如儲存血）。"
              "口訣：右移 = 'CADET face Right'（CO₂、Acid、DPG、Exercise、Temperature）。"},
    {"id": "bpg_effect", "title": "2,3-BPG and Hb O2 affinity", "layout": "plot", "plot": BPG,
     "caption": {"text": "2,3-BPG binds the central cavity of deoxy (T-state) Hb → stabilises T → right shift"},
     "tags": ["hemoglobin"],
     "notes": "2,3-BPG 結合在 deoxy-Hb（T 態）兩條 β 鏈之間的中央空腔（帶正電的 His143、Lys82、His2），穩定 T 態 → 曲線右移。"
              "去除 BPG 的 Hb 親和力極高（P₅₀ 約 8 mmHg），幾乎無法在組織卸氧。高海拔、慢性缺氧、貧血時 BPG 上升。"
              "HbF 的 γ 鏈 His143 → Ser，BPG 結合弱 → 親和力高於 HbA，讓胎兒能從母體取得 O₂。"},
    {"id": "co_poisoning", "title": "Carbon monoxide vs anemia", "layout": "plot", "plot": CO_CURVES,
     "caption": {"text": "CO: less O_2 carried AND the rest is released less easily → worse than anemia", "color": RED},
     "tags": ["hemoglobin"],
     "notes": "CO 與 heme 的親和力約為 O₂ 的 200–250 倍：50 % COHb 時攜氧量減半，且剩下的 heme 對 O₂ 親和力上升（左移、曲線變得較不 S 形），"
              "組織更難取得 O₂ —— 比同樣減半的貧血更危險。曲線為示意（O₂ content 以正常最大攜氧量為 1）。"},
    {"id": "hill_plot", "title": "Hill plot", "layout": "plot", "plot": HILL_PLOT, "tags": ["hemoglobin", "allostery"],
     "notes": "Hill plot：log[Y/(1−Y)] 對 log pO₂。斜率 = Hill 係數 n_H。Mb 為斜率 1 的直線（無協同性）；"
              "Hb 在 Y = 0.5 附近斜率約 3（n_H ≈ 2.8–3，正協同），兩端回到斜率 1（全 T 態／全 R 態）。"
              "n_H > 1 正協同、< 1 負協同；n_H 不等於結合位數目（Hb 有 4 個）。Hb 曲線以 MWC 模型計算。"},
    {"id": "titration_glycine", "title": "Titration curve of glycine", "layout": "plot", "plot": TITRATION_GLY,
     "tags": ["amino_acid"],
     "notes": "Glycine 的滴定曲線：兩個緩衝區，中心為 pK₁ = 2.34（α-COOH）與 pK₂ = 9.60（α-NH₃⁺）。"
              "加入 1 當量 OH⁻ 時全部為兩性離子，pH = pI = (pK₁ + pK₂)/2 = 5.97，淨電荷 0（此處斜率最大、緩衝力最差）。"
              "pH < pI 帶正電、往陰極移動；pH > pI 帶負電、往陽極移動。"},
    {"id": "ramachandran_plot", "title": "Ramachandran plot", "layout": "plot", "plot": RAMA,
     "tags": ["protein_structure"],
     "notes": "Ramachandran plot：φ（N–Cα）與 ψ（Cα–C）的允許組合（色塊為允許區，示意）。左上大片為 β 區，左下為右手 α-helix（φ ≈ −57°、ψ ≈ −47°），"
              "右上小區為左手 α-helix（少見，多為 Gly）。Gly 沒有側鏈，幾乎可落在任何區域；Pro 的 φ 固定在約 −60°。"},
    {"id": "peptide_bond_resonance", "title": "Peptide bond resonance", "layout": "free", "cell_size": 7,
     "nodes": [{"id": "a", "ref": "peptide_unit", "pos": [0, 0]}, {"id": "b", "ref": "peptide_unit_resonance", "pos": [1.5, 0]},
               {"id": "c", "text": "C–N bond: ~40 % double-bond character\n→ planar, rigid, no rotation; usually trans",
                "kind": "concept", "pos": [0.75, 0.85]}],
     "edges": [{"from": "a", "to": "b", "resonance": True}],
     "tags": ["peptide", "protein_structure"],
     "notes": "胜肽鍵的兩個共振式：O 帶部分負電、N 帶部分正電，C–N 鍵約 40 % 雙鍵性質 → 六個原子（Cα、C、O、N、H、Cα）共平面、不能自由旋轉。"
              "幾乎都是 trans（兩個 Cα 在對側，立體障礙小）；X–Pro 鍵約 10 % 為 cis（peptidyl prolyl isomerase 催化轉換）。"},
]


def dump(path, d):
    path.write_text(yaml.safe_dump(d, allow_unicode=True, sort_keys=False, width=200), encoding="utf-8")


def with_parent(tags):
    tree = {"enzyme_kinetics": "enzyme", "inhibition": "enzyme", "allostery": "enzyme",
            "hemoglobin": "protein", "peptide": "protein", "protein_structure": "protein"}
    out = []
    for t in tags:
        if t in tree:
            out.append(tree[t])
        out.append(t)
    return list(dict.fromkeys(out))


def main(data):
    data = Path(data)
    for m in MOLECULES:
        d = {"id": m["id"], "name": m["name"]}
        if m.get("label"):
            d["label"] = m["label"]
        d["aliases"] = m.get("aliases", [])
        d["smiles"] = m["smiles"]
        for k in ("chain", "free", "decorate"):
            if m.get(k):
                d[k] = m[k]
        d["tags"] = with_parent(m["tags"])
        d["notes_zh"] = m["notes"]
        d["status"] = "pending"
        dump(data / "molecules" / f"{m['id']}.yaml", d)
    for r in REACTIONS:
        d = {"id": r["id"], "title": r["title"], "substrates": r["substrates"], "products": r["products"]}
        for k in ("cofactors_in", "cofactors_out", "prosthetic"):
            if r.get(k):
                d[k] = r[k]
        d["enzyme"] = r["enzyme"]
        d["reversible"] = r.get("reversible", False)
        deco = COLORS.get(r["id"], r.get("decorate"))
        if deco:
            d["decorate"] = deco
        d["card_zh"] = r["card"]
        d["tags"] = with_parent(r["tags"])
        d["status"] = "pending"
        dump(data / "reactions" / f"{r['id']}.yaml", d)
    for p in PATHWAYS:
        d = {k: v for k, v in p.items() if k not in ("tags", "notes")}
        d["tags"] = with_parent(p["tags"])
        d["notes_zh"] = p["notes"]
        d["status"] = "pending"
        dump(data / "pathways" / f"{p['id']}.yaml", d)
    print(f"{len(MOLECULES)} molecules, {len(REACTIONS)} reactions, {len(PATHWAYS)} pathways")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data")
