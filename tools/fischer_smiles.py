"""Fischer drawing → isomeric SMILES (open-chain aldoses / ketoses)."""
import sys
sys.path.insert(0, ".")
from rdkit import Chem
from sdb.validate import _cip_from


def fischer_smiles(kind, sides):
    """kind: 'aldose' | 'ketose'; sides: 'R'/'L' for each stereocentre, top to bottom
    (the side the OH is drawn on). Returns isomeric SMILES."""
    n_st = len(sides)
    if kind == "aldose":
        base = "O=C" + "C(O)" * n_st + "CO"
        first = 1                       # atom index of the first stereocentre
    elif kind == "alditol":
        base = "OC" + "C(O)" * n_st + "CO"
    else:
        base = "OCC(=O)" + "C(O)" * n_st + "CO"
        first = 3
    mol = Chem.RWMol(Chem.MolFromSmiles(base))
    mol = Chem.AddHs(mol)
    # backbone carbons in order
    chain = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "C"]
    centres = []
    if kind in ("aldose", "alditol"):
        centres = chain[1:1 + n_st]
    else:
        centres = chain[2:2 + n_st]
    want = {}
    for i, c in enumerate(centres):
        atom = mol.GetAtomWithIdx(c)
        k = chain.index(c)
        up, down = chain[k - 1], chain[k + 1]
        oh = next(n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == "O")
        h = next(n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == "H")
        left, right = (oh, h) if sides[i] == "L" else (h, oh)
        pos = {c: (0, 0, 0), up: (0, 1, -0.7), down: (0, -1, -0.7), left: (-1, 0, 0.7), right: (1, 0, 0.7)}
        want[c] = _cip_from(mol, c, pos)
    m = Chem.Mol(mol)
    for c in centres:
        for tag in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            m.GetAtomWithIdx(c).SetChiralTag(tag)
            Chem.AssignStereochemistry(m, cleanIt=True, force=True)
            if m.GetAtomWithIdx(c).GetPropsAsDict().get("_CIPCode") == want[c]:
                break
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    for c in centres:
        assert m.GetAtomWithIdx(c).GetPropsAsDict().get("_CIPCode") == want[c]
    return Chem.MolToSmiles(Chem.RemoveHs(m))


if __name__ == "__main__":
    g = fischer_smiles("aldose", "RLRR")
    pub = "O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO"   # D-glucose already validated in the database
    print(g, Chem.CanonSmiles(g) == Chem.CanonSmiles(pub))
    print(fischer_smiles("ketose", "LRR"))
