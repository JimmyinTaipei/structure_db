"""Skeletal (line-angle) structures via RDKit, restyled to match the global style."""
import re

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from .. import style
from ..graphic import Graphic

S = style.S
BOND = 1.55 * S            # target bond length
R_LABELS = {1: "CoA", 2: "ADPR", 3: "FAD", 4: "Ado", 5: "Guo", 6: "ADPR-P", 9: "Glycogen", 16: "G", 17: "dR", 18: "MGMT", 19: "Rib", 20: "tRNA",
            22: "Enz", 23: "Pep_N", 24: "Pep_C", 25: "R_1", 26: "R_2"}   # see data/cofactors.yaml


def draw_skeletal(smiles, deco=None):
    """deco: {"smarts": "C=C", "color": "magenta"} recolours the matched bonds and atom labels."""
    deco = deco or {}
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"bad SMILES: {smiles}")
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)

    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 0:
            a.SetProp("atomLabel", R_LABELS.get(a.GetIsotope(), "R"))

    hl_atoms, hl_bonds = [], []
    if deco.get("smarts"):
        patt = Chem.MolFromSmarts(deco["smarts"])
        for match in mol.GetSubstructMatches(patt):
            hl_atoms += list(match)
            for b in patt.GetBonds():
                hl_bonds.append(mol.GetBondBetweenAtoms(match[b.GetBeginAtomIdx()],
                                                        match[b.GetEndAtomIdx()]).GetIdx())

    def render(font_px):
        d = rdMolDraw2D.MolDraw2DSVG(-1, -1)
        o = d.drawOptions()
        o.useBWAtomPalette()
        o.clearBackground = False
        o.bondLineWidth = 1
        o.scaleBondWidth = False
        o.multipleBondOffset = 0.18
        o.fontFile = str(style.SERIF_FONT)
        o.minFontSize = o.maxFontSize = max(6, round(font_px))
        o.padding = 0.03
        d.DrawMolecule(mol)
        d.FinishDrawing()
        lens = []
        for b in mol.GetBonds():
            p, q = d.GetDrawCoords(b.GetBeginAtomIdx()), d.GetDrawCoords(b.GetEndAtomIdx())
            lens.append(((p.x - q.x) ** 2 + (p.y - q.y) ** 2) ** 0.5)
        return d.GetDrawingText(), (sum(lens) / len(lens) if lens else BOND)

    # RDKit's flexible canvas ignores fixedBondLength, so measure and scale ourselves,
    # shrinking the font beforehand so that it ends up at size S after scaling.
    _, b0 = render(S)
    k = BOND / b0
    svg, _ = render(S * 1.1 / k)

    # colour the matched bonds / atom labels directly (no glow): RDKit tags every element
    # with class='bond-N atom-i atom-j' or class='atom-i'
    if hl_atoms:
        hex_col = style.color(deco.get("color", "magenta"))
        bonds, atoms = set(hl_bonds), set(hl_atoms)
        out = []
        for ln in svg.split("\n"):
            m = re.search(r"class='([^']*)'", ln)
            if m:
                cls = m.group(1).split()
                b = [int(c[5:]) for c in cls if c.startswith("bond-")]
                a = [int(c[5:]) for c in cls if c.startswith("atom-")]
                if (b and b[0] in bonds) or (not b and a and a[0] in atoms):
                    ln = ln.replace("#000000", hex_col)
            out.append(ln)
        svg = "\n".join(out)

    w = float(re.search(r"width='([\d.]+)px'", svg).group(1))
    h = float(re.search(r"height='([\d.]+)px'", svg).group(1))
    inner = svg[svg.index(">", svg.index("<svg")) + 1: svg.rindex("</svg>")]
    inner = re.sub(r"<rect[^>]*/>", "", inner, count=1) if "style='opacity:1.0;fill:#FFFFFF" in inner else inner
    inner = inner.replace("#000000", style.INK)
    inner = re.sub(r"stroke-width:[\d.]+px(?=;stroke-linecap)", f"stroke-width:{style.LINE_W / k:.2f}px", inner)
    # RDKit puts the origin at the top-left; centre horizontally so x=0 is the middle
    w, h = w * k, h * k
    body = f'<g transform="translate({-w / 2:.2f},{-h:.2f}) scale({k:.4f})">{inner}</g>'
    return Graphic(body, -w / 2, -h, w / 2, 0, {"base": (0, 0), "axis": (0, -h / 2)})
