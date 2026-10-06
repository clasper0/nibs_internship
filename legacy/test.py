# coding=utf-8
from pymatgen.io.cif import CifParser   # 更完整地读取 cif 文件信息（包括键长等）
from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.core.periodic_table import Element    # 获取元素半径
from rdkit import Chem
from rdkit.Chem import AllChem
from all_frags import Molecule
from all_frags import PROBE_FRAGMENTS, REF_FRAGMENTS
import copy
import numpy as np


mol = Chem.MolFromSmiles("[H]OS(O[H])(c1c([H])c([H])c(C([H])([H])[H])c([H])c1[H])c1c([H])c([H])c([H])c([H])c1[C@@]([H])([C@]([H])(c1c([H])c([H])c(OC([H])([H])[H])c([H])c1[H])N([H])c1c([H])c([H])c(C([H])([H])N([H])[H])c([H])c1[H])S(O[H])(O[H])C([H])([H])[H].[H]OS(O[H])(c1c([H])c([H])c(C([H])([H])[H])c([H])c1[H])c1c([H])c([H])c([H])c([H])c1[C@@]([H])([C@]([H])(c1c([H])c([H])c(OC([H])([H])[H])c([H])c1[H])N([H])c1c([H])c([H])c(C([H])([H])N([H])[H])c([H])c1[H])S(O[H])(O[H])C([H])([H])[H]")
print(Chem.MolToSmiles(mol))

# mol = Chem.MolFromMol2File("/home/clasper/2507/codes/MBZ.mol2")


mol1 = Chem.MolFromSmarts("c1([*])c([*])c([*])c([*])c([*])c1([*])")
mol2 = Chem.MolFromSmarts("c1ccccc1")

print(mol.GetSubstructMatches(mol1, useQueryQueryMatches=True))
print(mol.GetSubstructMatches(mol2, useQueryQueryMatches=True))
