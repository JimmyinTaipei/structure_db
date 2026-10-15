"""SMILES helpers for the sugar generator: glycosyl fragments that keep stereo tags."""
from rdkit import Chem

B_GLC = "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O"
A_GLC = "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O"
B_GAL = "C([C@@H]1[C@@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O"
A_GAL = "C([C@@H]1[C@@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O"
A_MAN = "C([C@@H]1[C@H]([C@@H]([C@@H]([C@H](O1)O)O)O)O)O"
B_FRU = "C([C@@H]1[C@H]([C@@H]([C@](O1)(CO)O)O)O)O"
A_FRU = "C([C@@H]1[C@H]([C@@H]([C@@](O1)(CO)O)O)O)O"


def anomeric_oh(mol):
    """Index of the OH oxygen on the anomeric carbon (the ring carbon bonded to two O)."""
    for a in mol.GetAtoms():
        if a.GetSymbol() == "C" and a.IsInRing():
            os_ = [n for n in a.GetNeighbors() if n.GetSymbol() == "O"]
            if len(os_) == 2:
                return next(n.GetIdx() for n in os_ if not n.IsInRing())
    raise ValueError("no anomeric carbon")


def glycosyl(smiles):
    """'…' rooted at the anomeric OH, with that O dropped: the fragment starts at C1 and its
    first neighbour (in SMILES order) is whatever atom precedes it — the glycosidic O."""
    mol = Chem.MolFromSmiles(smiles)
    s = Chem.MolToSmiles(mol, rootedAtAtom=anomeric_oh(mol), canonical=False)
    assert s.startswith("O"), s
    return s[1:]


def check(smiles):
    assert Chem.MolFromSmiles(smiles) is not None, smiles
    return smiles
