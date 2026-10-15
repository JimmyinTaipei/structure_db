"""Global drawing style. Every size is expressed relative to S (structure font size)."""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# Images are written to the iCloud copy so Typora notes and GoodNotes on the iPad keep working;
# the code itself lives outside iCloud. Override with SDB_OUT.
ICLOUD_OUT = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/後醫/08_考古題與模擬考/82_考題整理/structure_db/out"
OUT = Path(os.environ.get("SDB_OUT") or (ICLOUD_OUT if ICLOUD_OUT.parent.is_dir() else ROOT / "out"))
FONT_DIR = ROOT / "fonts"
SERIF_FONT = FONT_DIR / "texgyreschola-regular.otf"   # Century Schoolbook look (structures, cofactors)
SANS_FONT = FONT_DIR / "SourceSans3-Regular.ttf"               # Myriad look (names, enzymes)

S = 20.0                       # structure font size (SVG user units)
NAME_SIZE = S * 0.95           # molecule name under a structure
ENZYME_SIZE = S * 0.85         # enzyme / conditions text under arrows
LINE_W = S * 0.055             # bond & arrow stroke width
PNG_ZOOM = 3

# colours (semantic names are what YAML refers to)
COLORS = {
    "black": "#111111",
    "magenta": "#D6177F",      # transferred fragment / changed part
    "cyan": "#00A0DC",         # enzyme text, traced fragment
    "salmon": "#F4C3AE",       # highlight box behind changed atoms
    "gray": "#777777",
    "red": "#D62828",          # inhibition
    "orange": "#E8731A",       # isotope-tracing marks (* ○ △ ①)
    "yellow": "#F7E04B",       # Mathews-style block behind S (Cys, Met)
    "green": "#2E9E4F",        # extra source colours (atom-origin figures)
    "purple": "#7B4FC0",
    "membrane": "#F6EAC2",     # lipid bilayer band (ETC, shuttles)
    "membrane_head": "#D9B45A",
    "protein": "#DCE9F7",      # protein complexes in membrane schematics
    "protein_line": "#4F7FB5",
    "carrier": "#FCE6CC",      # mobile carriers (Q, cyt c)
    "carrier_line": "#D9822B",
    "compartment": "#9A9A9A",  # organ / organelle boxes
}
INK = COLORS["black"]
ENZYME_COLOR = COLORS["cyan"]
BOX_COLOR = COLORS["salmon"]


def color(name):
    """Accept a semantic name or a raw '#rrggbb'."""
    if name is None:
        return INK
    return COLORS.get(name, name)
