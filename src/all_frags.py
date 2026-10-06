from dataclasses import dataclass

@dataclass(frozen=True)
class Molecule:
    name: str
    smarts: str
    smiles: str
    # heavy_atom: list
    # # num_heavy : int
    # num_all : int
    # heavy_atom_idx_to_saveH: list

ACEM = Molecule(
    name = "ACEM",
    smarts = "C([*])([*])([*])C(N([*])([*]))=O",
    smiles = "CC([NH2])=O"
)
ACET = Molecule(
    name = "ACET",
    smarts = "C([*])([*])([*])C(=O)[O-]",
    smiles = "CC([O-])=O"
)
ETAM = Molecule(
    name = "ETAM",
    smarts = "C([*])([*])([*])C([*])([*])[N+]([*])([*])([*])",
    smiles = "CC[NH3+]"
)
ETOH = Molecule(
    name = "ETOH",
    smarts = "C([*])([*])([*])C([*])([*])O([*])",
    smiles = "CC[OH]"   # 添加极性氢
)
ETSH = Molecule(
    name = "ETSH",
    smarts = "C([*])([*])([*])C([*])([*])S([*])",
    smiles = "CC[SH]"
)
HOH = Molecule(
    name = "HOH",
    smarts = "O([H])([H])", # 水必须把H显示写
    smiles = "O"
)
MBZ = Molecule(
    name = "MBZ",
    smarts = "C([*])([*])([*])c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "Cc1ccccc1"
)
MGDM = Molecule(
    name = "MGDM",
    smarts = "C([*])([*])([*])N([*])C(N([*])([*]))=N([*])([*])",  # 先不加上"+"
    smiles = "CNC(N)=[NH2+]"
)
MIMD = Molecule(
    name = "MIMD",
    smarts = "C([*])([*])([*])c1c([*])nc([*])n([*])1",
    smiles = "Cc1cnc[nH]1"
)
MIME = Molecule(
    name = "MIME",
    smarts = "C([*])([*])([*])c1c([*])n([*])c([*])n1",
    smiles = "Cc1c[nH]cn1"
)
MIMM = Molecule(
    name = "MIMM",
    smarts = "C([*])([*])([*])c1c([*])[n+]([*])c([*])n([*])1",
    smiles = "Cc1c[nH+]c[nH]1"
)
MIND = Molecule(
    name = "MIND",
    smarts = "C([*])([*])([*])c1c([*])n([*])c2c([*])c([*])c([*])c([*])c12",
    smiles = "Cc1c[nH]c2ccccc12"
)
MSM = Molecule(
    name = "MSM",
    smarts = "C([*])([*])([*])SC([*])([*])([*])",
    smiles = "CSC"
)
MTYR = Molecule(
    name = "MTYR",
    smarts = "C([*])([*])([*])c1c([*])c([*])c(O([*]))c([*])c1([*])",
    smiles = "Cc1ccc([OH])cc1"  # 添加极性氢
)
N1PA = Molecule(
    name = "N1PA",
    smarts = "C([*])([*])([*])C(=O)N1C([*])([*])C([*])([*])C([*])([*])C1([*])([*])",
    smiles = "CC(=O)N1CCCC1"
)
NMA = Molecule(
    name = "NMA",
    smarts = "C([*])([*])([*])N([*])C(=O)C([*])([*])([*])", # 暂时将～改成=键
    smiles = "C[NH]C(=O)C"  # 添加极性氢
)
PRPA = Molecule(
    name = "PRPA",
    smarts = "C([*])([*])([*])C([*])([*])C([*])([*])([*])",
    smiles = "CCC"
)

frag_1 = Molecule(
    name = "1",
    smarts = "C([*])([*])([*])([*])", # 少写一个
    smiles = "C"
)
frag_2 = Molecule(
    name = "2",
    smarts = "C([*])([*])([*])C([*])(C([*])([*])([*]))C([*])([*])([*])",
    smiles = "CC(C)C"
)
frag_3 = Molecule(
    name = "3",
    smarts = "C([*])([*])=C([*])([*])",
    smiles = "C=C"
)
frag_4 = Molecule(
    name = "4",
    smarts = "C([*])#CC([*])([*])([*])",
    smiles = "C#CC"
)
frag_5 = Molecule(
    name = "5",
    smarts = "C1([*])([*])C([*])([*])C1([*])([*])",
    smiles = "C1CC1"
)
frag_6 = Molecule(
    name = "6",
    smarts = "C([*])([*])=C([*])C([*])=C([*])C([*])([*])([*])",
    smiles = "C=CC=CC"
)
frag_7 = Molecule(
    name = "7",
    smarts = "C1([*])([*])C([*])([*])C([*])([*])C1([*])([*])",
    smiles = "C1CCC1"
)
frag_8 = Molecule(
    name = "8",
    smarts = "C1([*])([*])C([*])([*])C([*])([*])C([*])([*])C([*])([*])C1([*])([*])",
    smiles = "C1CCCCC1"
)
frag_9 = Molecule(
    name = "9",
    smarts = "C([*])([*])([*])C1([*])C([*])([*])C2([*])C([*])([*])C([*])([*])C1([*])C2([*])([*])",
    smiles = "CC1CC2CCC1C2"
)
frag_10 = Molecule(
    name = "10",
    smarts = "C([*])([*])([*])N([*])([*])",
    smiles = "C[NH2]"
)
frag_11 = Molecule(
    name = "11",
    smarts = "C([*])([*])([*])C#N",
    smiles = "CC#N"
)
frag_12 = Molecule(
    name = "12",
    smarts = "C([*])([*])(F)F",
    smiles = "C(F)F"
)
frag_13 = Molecule(
    name = "13",
    smarts = "C([*])([*])([*])SC([*])([*])([*])",   # 与MSM相同
    smiles = "CSC"
)
frag_14 = Molecule(
    name = "14",
    smarts = "C([*])([*])([*])O([*])",
    smiles = "C[OH]"    # 添加极性氢
)
frag_15 = Molecule(
    name = "15",
    smarts = "C([*])([*])([*])C([*])([*])O([*])",
    smiles = "CC[OH]"  # 添加极性氢
)
frag_16 = Molecule(
    name = "16",
    smarts = "C([*])([*])([*])C([*])(C([*])([*])([*]))O([*])",
    smiles = "CC(C)[OH]"  # 添加极性氢
)
frag_17 = Molecule(
    name = "17",
    smarts = "C([*])([*])([*])C(C([*])([*])([*]))(C([*])([*])([*]))O([*])",
    smiles = "CC(C)(C)[OH]"  # 添加极性氢
)
frag_18 = Molecule(
    name = "18",
    smarts = "C([*])([*])([*])OC([*])([*])([*])",
    smiles = "COC"
)
frag_19 = Molecule(
    name = "19",
    smarts = "C([*])([*])([*])C([*])(C([*])([*])([*]))OC([*])([*])([*])",
    smiles = "CC(C)OC"
)
frag_20 = Molecule(
    name = "20",
    smarts = "C([*])([*])([*])OC([*])([*])C([*])([*])OC([*])([*])([*])",
    smiles = "COCCOC"
)
frag_21 = Molecule(
    name = "21",
    smarts = "C([*])([*])([*])C(=O)C([*])([*])([*])",
    smiles = "CC(=O)C"
)
frag_22 = Molecule(
    name = "22",
    smarts = "C([*])([*])([*])C(=O)O([*])",
    smiles = "CC(=O)[OH]"  # 添加极性氢
)
frag_23 = Molecule(
    name = "23",
    smarts = "C([*])([*])([*])C(=O)OC([*])([*])([*])",
    smiles = "CC(=O)OC"
)
frag_24 = Molecule(
    name = "24",
    smarts = "C([*])([*])([*])C(=O)OC([*])(C([*])([*])([*]))C([*])([*])([*])", # 少一个[*]
    smiles = "CC(=O)OC(C)C"
)
frag_25 = Molecule(
    name = "25",
    smarts = "C([*])([*])(C([*])(C1([*])C([*])(C([*])(C([*])(O([*]))O1)O([*]))O([*]))O([*]))O([*])",
    smiles = "C(C(C1C(C(C([OH])O1)[OH])[OH])[OH])[OH]"  # 添加极性氢
)
frag_26 = Molecule(
    name = "26",
    smarts = "C([*])(=O)N([*])([*])",
    smiles = "C(=O)[NH2]"   # 添加极性氢
)
frag_27 = Molecule(
    name = "27",
    smarts = "C([*])([*])([*])C(=O)N([*])([*])",
    smiles = "CC(=O)[NH2]"  # 添加极性氢
)
frag_28 = Molecule(
    name = "28",
    smarts = "C1([*])([*])C([*])([*])C1([*])C(=O)N([*])([*])",
    smiles = "C1CC1C(=O)[NH2]"  # 添加极性氢
)
frag_29 = Molecule(
    name = "29",
    smarts = "C([*])([*])([*])C(=O)N([*])C([*])([*])([*])",
    smiles = "CC(=O)[NH]C"  # 添加极性氢
)
frag_30 = Molecule(
    name = "30",
    smarts = "C([*])([*])([*])C(=O)N([*])C1([*])C([*])([*])C1([*])([*])",
    smiles = "CC(=O)[NH]C1CC1"  # 添加极性氢
)
frag_31 = Molecule(
    name = "31",
    smarts = "C([*])([*])([*])C(=O)N(C([*])([*])([*]))C([*])([*])([*])",
    smiles = "CC(=O)N(C)C"
)
frag_32 = Molecule(
    name = "32",
    smarts = "C([*])([*])([*])C(=O)N(C([*])([*])([*]))c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "CC(=O)N(C)c1ccccc1"
)
frag_33 = Molecule(
    name = "33",
    smarts = "C([*])([*])([*])C(=O)N([*])O([*])",
    smiles = "CC(=O)[NH][OH]"  # 添加极性氢
)
frag_34 = Molecule(
    name = "34",
    smarts = "C([*])([*])([*])C(=S)N([*])C([*])([*])([*])",
    smiles = "CC(=S)[NH]C"  # 添加极性氢
)
frag_35 = Molecule(
    name = "35",
    smarts = "C([*])([*])([*])N([*])C(=O)N([*])([*])",
    smiles = "C[NH]C(=O)[NH2]"  # 添加极性氢
)
frag_36 = Molecule(
    name = "36",
    smarts = "C([*])([*])([*])N([*])C(=O)N([*])C([*])([*])([*])",
    smiles = "C[NH]C(=O)[NH]C"  # 添加极性氢
)
frag_37 = Molecule(
    name = "37",
    smarts = "C([*])([*])([*])N(C([*])([*])([*]))C(=O)N([*])([*])",
    smiles = "CN(C)C(=O)[NH2]"  # 添加极性氢
)
frag_38 = Molecule(
    name = "38",
    smarts = "C([*])([*])([*])N([*])C(=O)N(C([*])([*])([*]))C([*])([*])([*])",
    smiles = "C[NH]C(=O)N(C)C"
)
frag_39 = Molecule(
    name = "39",
    smarts = "C([*])([*])([*])OC(=O)N([*])([*])",
    smiles = "COC(=O)[NH2]" # 添加极性氢
)
frag_40 = Molecule(
    name = "40",
    smarts = "C([*])([*])([*])N([*])C(=O)OC([*])([*])([*])",
    smiles = "C[NH]C(=O)OC" # 添加极性氢
)
frag_41 = Molecule(
    name = "41",
    smarts = "C([*])([*])([*])N([*])C(=S)N([*])C([*])([*])([*])",
    smiles = "C[NH]C(=S)[NH]C"  # 添加极性氢
)
frag_42 = Molecule(
    name = "42",
    smarts = "C([*])([*])([*])S(=O)C([*])([*])([*])",
    smiles = "CS(=O)C"
)
frag_43 = Molecule(
    name = "43",
    smarts = "C([*])([*])([*])S(=O)(=O)C([*])([*])([*])",
    smiles = "CS(=O)(=O)C"
)
frag_44 = Molecule(
    name = "44",
    smarts = "C([*])([*])([*])S(=O)(=O)N([*])([*])",
    smiles = "CS(=O)(=O)[NH2]"  # 添加极性氢
)
frag_45 = Molecule(
    name = "45",
    smarts = "C([*])([*])([*])N([*])S(=O)(=O)C([*])([*])([*])",
    smiles = "C[NH]S(=O)(=O)C"  # 添加极性氢
)
frag_46 = Molecule(
    name = "46",
    smarts = "c1([*])c([*])c([*])c([*])c([*])c1([*])",
    smiles = "c1ccccc1"
)
frag_47 = Molecule(
    name = "47",
    smarts = "C([*])([*])([*])c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "Cc1ccccc1"
)
frag_48 = Molecule(
    name = "48",
    smarts = "C([*])([*])([*])C([*])(C([*])([*])([*]))c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "CC(C)c1ccccc1"
)
frag_49 = Molecule(
    name = "49",
    smarts = "c1([*])c([*])c([*])c(C2([*])C([*])([*])C2([*])([*]))c([*])c1([*])",
    smiles = "c1ccc(C2CC2)cc1"
)
frag_50 = Molecule(
    name = "50",
    smarts = "C([*])([*])([*])C(C([*])([*])([*]))(C([*])([*])([*]))c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "CC(C)(C)c1ccccc1"
)
frag_51 = Molecule(
    name = "51",
    smarts = "C([*])([*])([*])c1c(C([*])([*])([*]))c([*])c([*])c([*])c1([*])",
    smiles = "Cc1c(C)cccc1"
)
frag_52 = Molecule(
    name = "52",
    smarts = "C([*])([*])=C([*])c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "C=Cc1ccccc1"
)
frag_53 = Molecule(
    name = "53",
    smarts = "C([*])#Cc1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "C#Cc1ccccc1"
)
frag_54 = Molecule(
    name = "54",
    smarts = "FC(F)(F)c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "FC(F)(F)c1ccccc1"
)
frag_55 = Molecule(
    name = "55",
    smarts = "Fc1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "Fc1ccccc1"
)
frag_56 = Molecule(
    name = "56",
    smarts = "Clc1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "Clc1ccccc1"
)
frag_57 = Molecule(
    name = "57",
    smarts = "Brc1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "Brc1ccccc1"
)
frag_58 = Molecule(
    name = "58",
    smarts = "Fc1c([*])c([*])c([*])c([*])c1F",
    smiles = "Fc1ccccc1F"
)
frag_59 = Molecule(
    name = "59",
    smarts = "Clc1c([*])c([*])c([*])c([*])c1Cl",
    smiles = "Clc1ccccc1Cl"
)
frag_60 = Molecule(
    name = "60",
    smarts = "C([*])([*])([*])N([*])c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "C[NH]c1ccccc1"    # 添加极性氢
)
frag_61 = Molecule(
    name = "61",
    smarts = "C([*])([*])([*])Oc1c([*])c([*])c([*])c([*])c1N([*])([*])",
    smiles = "COc1ccccc1[NH2]"  # 添加极性氢
)
frag_62 = Molecule(
    name = "62",
    smarts = "O([*])c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "[OH]c1ccccc1" # 添加极性氢
)
frag_63 = Molecule(
    name = "63",
    smarts = "C([*])([*])([*])Oc1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "COc1ccccc1"
)
frag_64 = Molecule(
    name = "64",
    smarts = "C([*])([*])([*])Oc1c([*])c([*])c([*])c([*])c1O([*])",
    smiles = "COc1ccccc1[OH]"  # 添加极性氢
)
frag_65 = Molecule(
    name = "65",
    smarts = "C([*])([*])([*])Oc1c([*])c([*])c([*])c([*])c1OC([*])([*])([*])",
    smiles = "COc1ccccc1OC"
)
frag_66 = Molecule(
    name = "66",
    smarts = "O([*])C([*])([*])c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "[OH]Cc1ccccc1"    # 添加极性氢
)
frag_67 = Molecule(
    name = "67",
    smarts = "C([*])([*])([*])C([*])(O([*]))c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "CC([OH])c1ccccc1" # 添加极性氢
)
frag_68 = Molecule(
    name = "68",
    smarts = "FC(F)(F)Oc1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "FC(F)(F)Oc1ccccc1"
)
frag_69 = Molecule(
    name = "69",
    smarts = "O([*])c1c([*])c([*])c(F)c(F)c1F",
    smiles = "[OH]c1ccc(F)c(F)c1F"  # 添加极性氢
)
frag_70 = Molecule(
    name = "70",
    smarts = "Fc1c(Cl)c(Cl)c([*])c([*])c1([*])",
    smiles = "Fc1c(Cl)c(Cl)ccc1"
)
frag_71 = Molecule(
    name = "71",
    smarts = "N#Cc1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "N#Cc1ccccc1"
)
frag_72 = Molecule(
    name = "72",
    smarts = "S([*])c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "[SH]c1ccccc1" # 添加极性氢
)
frag_73 = Molecule(
    name = "73",
    smarts = "C([*])([*])([*])Sc1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "CSc1ccccc1"
)
frag_74 = Molecule(
    name = "74",
    smarts = "N([*])([*])C(=O)c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "[NH2]C(=O)c1ccccc1"   # 添加极性氢
)
frag_75 = Molecule(
    name = "75",
    smarts = "C([*])([*])([*])N([*])C(=O)c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "C[NH]C(=O)c1ccccc1"   # 添加极性氢
)
frag_76 = Molecule(
    name = "76",
    smarts = "C([*])([*])([*])N(C([*])([*])([*]))C(=O)c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "CN(C)C(=O)c1ccccc1"
)
frag_77 = Molecule(
    name = "77",
    smarts = "C([*])([*])([*])OC(=O)c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "COC(=O)c1ccccc1"
)
frag_78 = Molecule(
    name = "78",
    smarts = "C([*])([*])([*])C(=O)c1c([*])c([*])c([*])c([*])c1O([*])",
    smiles = "CC(=O)c1ccccc1[OH]"   # 添加极性氢
)
frag_79 = Molecule(
    name = "79",
    smarts = "C([*])([*])([*])OC(=O)c1c([*])c([*])c([*])c([*])c1O([*])",
    smiles = "COC(=O)c1ccccc1[OH]"  # 添加极性氢
)
frag_80 = Molecule(
    name = "80",
    smarts = "O=[N+]([O-])c1c([*])c([*])c([*])c([*])c1([*])", # 不确定要不要这么严格写电荷
    smiles = "O=[N+]([O-])c1ccccc1"
)
frag_81 = Molecule(
    name = "81",
    smarts = "C([*])([*])([*])C(=O)N([*])c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "CC(=O)[NH]c1ccccc1"   # 添加极性氢
)
frag_82 = Molecule(
    name = "82",
    smarts = "C([*])([*])([*])C(=O)N([*])c1c([*])c([*])c([*])c([*])c1Cl",
    smiles = "CC(=O)[NH]c1ccccc1Cl" # 添加极性氢
)
frag_83 = Molecule(
    name = "83",
    smarts = "N([*])([*])C(=O)c1c([*])c([*])c([*])c([*])c1Cl",
    smiles = "[NH2]C(=O)c1ccccc1Cl" # 添加极性氢
)
frag_84 = Molecule(
    name = "84",
    smarts = "N([*])([*])C(=O)c1c(F)c([*])c([*])c([*])c1F",
    smiles = "[NH2]C(=O)c1c(F)cccc1F"   # 添加极性氢
)
frag_85 = Molecule(
    name = "85",
    smarts = "O=C([*])N([*])c1c([*])c([*])c(F)c([*])c1F",
    smiles = "O=C[NH]c1ccc(F)cc1F"  # 添加极性氢
)
frag_86 = Molecule(
    name = "86",
    smarts = "N([*])([*])C(=S)N([*])c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "[NH2]C(=S)[NH]c1ccccc1"   # 添加极性氢
)
frag_87 = Molecule(
    name = "87",
    smarts = "C([*])([*])([*])S(=O)(=O)c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "CS(=O)(=O)c1ccccc1"
)
frag_88 = Molecule(
    name = "88",
    smarts = "N([*])([*])S(=O)(=O)c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "[NH2]S(=O)(=O)c1ccccc1"
)
frag_89 = Molecule(
    name = "89",
    smarts = "C([*])([*])([*])N([*])S(=O)(=O)c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "C[NH]S(=O)(=O)c1ccccc1"   # 添加极性氢
)
frag_90 = Molecule(
    name = "90",
    smarts = "N([*])([*])c1c([*])c([*])c([*])c([*])c1S(=O)(=O)N([*])([*])",
    smiles = "[NH2]c1ccccc1S(=O)(=O)[NH2]"  # 添加极性氢
)
frag_91 = Molecule(
    name = "91",
    smarts = "C([*])([*])([*])C1([*])C([*])([*])O1",
    smiles = "CC1CO1"
)
frag_92 = Molecule(
    name = "92",
    smarts = "C1([*])([*])C([*])([*])C([*])([*])OC1([*])([*])",
    smiles = "C1CCOC1"
)
frag_93 = Molecule(
    name = "93",
    smarts = "C1([*])([*])C([*])([*])OC([*])([*])OC1([*])([*])",
    smiles = "C1COCOC1"
)
frag_94 = Molecule(
    name = "94",
    smarts = "C1([*])([*])C([*])([*])C([*])([*])N([*])C1([*])([*])",
    smiles = "C1CC[NH]C1"   # 添加极性氢
)
frag_95 = Molecule(
    name = "95",
    smarts = "C1([*])([*])C([*])([*])N([*])C1=O",
    smiles = "O=C1CC[NH]1"  # 添加极性氢
)
frag_96 = Molecule(
    name = "96",
    smarts = "C([*])([*])([*])N1C([*])([*])C([*])([*])C1=O",
    smiles = "CN1CCC1=O"
)
frag_97 = Molecule(
    name = "97",
    smarts = "C1([*])([*])C([*])([*])C(=O)N([*])C1([*])([*])",
    smiles = "O=C1CCC[NH]1" # 添加极性氢
)
frag_98 = Molecule(
    name = "98",
    smarts = "C([*])([*])([*])N1C([*])([*])C([*])([*])C([*])([*])C1=O",
    smiles = "CN1CCCC1=O"
)
frag_99 = Molecule(
    name = "99",
    smarts = "C1([*])([*])C([*])([*])C([*])([*])N([*])C(=O)C1([*])([*])",
    smiles = "O=C1CCCC[NH]1"    # 添加极性氢
)
frag_100 = Molecule(
    name = "100",
    smarts = "C([*])([*])([*])N1C([*])([*])C([*])([*])C([*])([*])C([*])([*])C1=O",
    smiles = "CN1CCCCC1=O"
)
frag_101 = Molecule(
    name = "101",
    smarts = "C1([*])([*])C([*])([*])N([*])C(=O)N1([*])",
    smiles = "O=C1[NH]CC[NH]1"  # 添加极性氢
)
frag_102 = Molecule(
    name = "102",
    smarts = "C([*])([*])([*])N1C([*])([*])C([*])([*])N([*])C1=O",
    smiles = "CN1CC[NH]C1=O"    # 添加极性氢
)
frag_103 = Molecule(
    name = "103",
    smarts = "C([*])([*])([*])N1C([*])([*])C([*])([*])N(C([*])([*])([*]))C1=O",
    smiles = "CN1CCN(C)C1=O"
)
frag_104 = Molecule(
    name = "104",
    smarts = "C1([*])([*])C([*])([*])N([*])C(=O)N([*])C1([*])([*])",
    smiles = "O=C1[NH]CCC[NH]1" # 添加极性氢
)
frag_105 = Molecule(
    name = "105",
    smarts = "C([*])([*])([*])N1C([*])([*])C([*])([*])C([*])([*])N([*])C1=O",
    smiles = "CN1CCC[NH]C1=O"
)
frag_106 = Molecule(
    name = "106",
    smarts = "C([*])([*])([*])N1C([*])([*])C([*])([*])C([*])([*])N(C([*])([*])([*]))C1=O",
    smiles = "CN1CCCN(C)C1=O"
)
frag_107 = Molecule(
    name = "107",
    smarts = "C1([*])([*])C([*])([*])C(=O)N([*])C(=O)C1([*])([*])",
    smiles = "O=C1CCCC(=O)[NH]1"    # 添加极性氢
)
frag_108 = Molecule(
    name = "108",
    smarts = "c1([*])c([*])c([*])n([*])c1([*])",
    smiles = "c1cc[nH]c1"
)
frag_109 = Molecule(
    name = "109",
    smarts = "c1([*])c([*])n([*])nc1([*])", 
    smiles = "c1cn[nH]c1"
)
frag_110 = Molecule(
    name = "110",
    smarts = "c1([*])c([*])n([*])c([*])n1",
    smiles = "c1c[nH]cn1"
)
frag_111 = Molecule(
    name = "111",
    smarts = "c1([*])c([*])n([*])nn1",
    smiles = "c1c[nH]nn1"
)
frag_112 = Molecule(
    name = "112",
    smarts = "c1([*])nc([*])nn1([*])",
    smiles = "c1ncn[nH]1"
)
frag_113 = Molecule(
    name = "113",
    smarts = "c1([*])c([*])c([*])oc1([*])",
    smiles = "c1ccoc1"
)
frag_114 = Molecule(
    name = "114",
    smarts = "c1([*])c([*])onc1([*])",
    smiles = "c1conc1"
)
frag_115 = Molecule(
    name = "115",
    smarts = "c1([*])c([*])oc([*])n1",
    smiles = "c1cocn1"
)
frag_116 = Molecule(
    name = "116",
    smarts = "c1([*])nc([*])no1",
    smiles = "c1ncno1"
)
frag_117 = Molecule(
    name = "117",
    smarts = "c1([*])nnc([*])o1",
    smiles = "c1nnco1"
)
frag_118 = Molecule(
    name = "118",
    smarts = "c1([*])nonc1([*])",
    smiles = "c1nonc1"
)
frag_119 = Molecule(
    name = "119",
    smarts = "c1([*])c([*])c([*])sc1([*])",
    smiles = "c1ccsc1"
)
frag_120 = Molecule(
    name = "120",
    smarts = "c1([*])c([*])snc1([*])",
    smiles = "c1csnc1"
)
frag_121 = Molecule(
    name = "121",
    smarts = "c1([*])c([*])sc([*])n1",
    smiles = "c1cscn1"
)
frag_122 = Molecule(
    name = "122",
    smarts = "c1([*])nnc([*])s1",
    smiles = "c1nncs1"
)
frag_123 = Molecule(
    name = "123",
    smarts = "c1([*])nsnc1([*])",
    smiles = "c1nsnc1"
)
frag_124 = Molecule(
    name = "124",
    smarts = "C([*])([*])([*])n1nc([*])c([*])c1([*])",  
    smiles = "Cn1nccc1"
)
frag_125 = Molecule(
    name = "125",
    smarts = "C([*])([*])([*])c1c([*])nc([*])n1([*])",
    smiles = "Cc1cnc[nH]1"
)
frag_126 = Molecule(
    name = "126",
    smarts = "c2([*])c([*])n([*])c(C1([*])C([*])([*])C1([*])([*]))n2",
    smiles = "c2c[nH]c(C1CC1)n2"
)
frag_127 = Molecule(
    name = "127",
    smarts = "O=C1C([*])([*])C([*])([*])C([*])([*])C([*])([*])C1([*])n1c([*])c([*])nn1",
    smiles = "O=C1CCCCC1n1ccnn1"
)
frag_128 = Molecule(
    name = "128",
    smarts = "C([*])([*])([*])n1c([*])nc([*])n1",
    smiles = "Cn1cncn1"
)
frag_129 = Molecule(
    name = "129",
    smarts = "C([*])([*])([*])c1nnc([*])n1([*])",
    smiles = "Cc1nnc[nH]1"
)
frag_130 = Molecule(
    name = "130",
    smarts = "C([*])([*])([*])c1nnc([*])n1C([*])([*])([*])",
    smiles = "Cc1nncn1C"
)
frag_131 = Molecule(
    name = "131",
    smarts = "C([*])([*])([*])C([*])(C([*])([*])([*]))n1c([*])nc([*])n1",
    smiles = "CC(C)n1cncn1"
)
frag_132 = Molecule(
    name = "132",
    smarts = "C([*])([*])([*])C([*])(C([*])([*])([*]))c1nc([*])nn1([*])",
    smiles = "CC(C)c1ncn[nH]1"
)
frag_133 = Molecule(
    name = "133",
    smarts = "C([*])([*])([*])c1c([*])c([*])c([*])o1",
    smiles = "Cc1ccco1"
)
frag_134 = Molecule(
    name = "134",
    smarts = "C([*])([*])([*])c1c([*])c(C([*])([*])([*]))no1",
    smiles = "Cc1cc(C)no1"
)
frag_135 = Molecule(
    name = "135",
    smarts = "C([*])([*])([*])C([*])(C([*])([*])([*]))c1nc([*])no1",
    smiles = "CC(C)c1ncno1"
)
frag_136 = Molecule(
    name = "136",
    smarts = "C([*])([*])([*])C(C([*])([*])([*]))(C([*])([*])([*]))c1noc([*])n1",
    smiles = "CC(C)(C)c1nocn1"
)
frag_137 = Molecule(
    name = "137",
    smarts = "C([*])([*])([*])c1nnc([*])o1",
    smiles = "Cc1nnco1"
)
frag_138 = Molecule(
    name = "138",
    smarts = "C([*])([*])([*])c1c([*])c([*])c([*])s1",
    smiles = "Cc1cccs1"
)
frag_139 = Molecule(
    name = "139",
    smarts = "C([*])([*])([*])c1c([*])nc([*])s1",
    smiles = "Cc1cncs1"
)
frag_140 = Molecule(
    name = "140",
    smarts = "C([*])([*])([*])C([*])(C([*])([*])([*]))c1nsc([*])c1([*])",
    smiles = "CC(C)c1nscc1"
)
frag_141 = Molecule(
    name = "141",
    smarts = "O=c1n([*])c(=O)c([*])c1([*])",
    smiles = "O=c1[nH]c(=O)cc1"
)
frag_142 = Molecule(
    name = "142",
    smarts = "N([*])([*])c1nc([*])c([*])n([*])1",
    smiles = "[NH2]c1ncc[nH]1"
)
frag_143 = Molecule(
    name = "143",
    smarts = "N([*])([*])c1nc([*])c([*])s1",
    smiles = "[NH2]c1nccs1"
)
frag_144 = Molecule(
    name = "144",
    smarts = "C([*])([*])([*])c1c([*])nc(N([*])([*]))s1",
    smiles = "Cc1cnc([NH2])s1"
)
frag_145 = Molecule(
    name = "145",
    smarts = "C([*])([*])([*])c1nnc(N([*])([*]))s1",
    smiles = "Cc1nnc([NH2])s1"
)
frag_146 = Molecule(
    name = "146",
    smarts = "C([*])([*])([*])C(=O)N([*])c1nc([*])c([*])n([*])1",
    smiles = "CC(=O)[NH]c1ncc[nH]1"
)
frag_147 = Molecule(
    name = "147",
    smarts = "C([*])([*])([*])C(=O)N([*])c1nc([*])c([*])s1",
    smiles = "CC(=O)[NH]c1nccs1"
)
frag_148 = Molecule(
    name = "148",
    smarts = "C([*])([*])([*])OC(=O)c1c(C([*])([*])([*]))c([*])c([*])n([*])1",
    smiles = "COC(=O)c1c(C)cc[nH]1"
)
frag_149 = Molecule(
    name = "149",
    smarts = "C([*])([*])([*])N([*])C(=O)c1c([*])c([*])c([*])o1",
    smiles = "C[NH]C(=O)c1ccco1"
)
frag_150 = Molecule(
    name = "150",
    smarts = "N([*])([*])C(=O)c1c([*])nc([*])s1",
    smiles = "[NH2]C(=O)c1cncs1"
)
frag_151 = Molecule(
    name = "151",
    smarts = "C([*])([*])([*])OC(=O)c1c([*])nc([*])s1",
    smiles = "COC(=O)c1cncs1"
)
frag_152 = Molecule(
    name = "152",
    smarts = "C([*])([*])([*])C(=O)c1c([*])sc(C([*])([*])([*]))n1",
    smiles = "CC(=O)c1csc(C)n1"
)
frag_153 = Molecule(
    name = "153",
    smarts = "C([*])([*])([*])S(=O)(=O)c1c([*])c([*])c([*])s1",
    smiles = "CS(=O)(=O)c1cccs1"
)
frag_154 = Molecule(
    name = "154",
    smarts = "c1([*])c([*])c([*])nc([*])c1([*])",
    smiles = "c1ccncc1"
)
frag_155 = Molecule(
    name = "155",
    smarts = "C([*])([*])([*])c1nc([*])c([*])c([*])c1([*])",
    smiles = "Cc1ncccc1"
)
frag_156 = Molecule(
    name = "156",
    smarts = "C([*])([*])([*])C([*])(C([*])([*])([*]))c1nc([*])c([*])c([*])c1([*])",
    smiles = "CC(C)c1ncccc1"
)
frag_157 = Molecule(
    name = "157",
    smarts = "Clc1nc([*])c([*])c([*])c1([*])",
    smiles = "Clc1ncccc1"
)
frag_158 = Molecule(
    name = "158",
    smarts = "N([*])([*])c1nc([*])c([*])c([*])c1([*])",
    smiles = "[NH2]c1ncccc1"
)
frag_159 = Molecule(
    name = "159",
    smarts = "C([*])([*])([*])Oc1c([*])c([*])c([*])c([*])n1",
    smiles = "COc1ccccn1"
)
frag_160 = Molecule(
    name = "160",
    smarts = "C([*])([*])([*])Sc1nc([*])c([*])c([*])c1([*])",
    smiles = "CSc1ncccc1"
)
frag_161 = Molecule(
    name = "161",
    smarts = "N([*])([*])C(=O)c1nc([*])c([*])c([*])c1([*])",
    smiles = "[NH2]C(=O)c1ncccc1"
)
frag_162 = Molecule(
    name = "162",
    smarts = "C([*])([*])([*])N([*])C(=O)c1nc([*])c([*])c([*])c1([*])",
    smiles = "C[NH]C(=O)c1ncccc1"
)
frag_163 = Molecule(
    name = "163",
    smarts = "C([*])([*])([*])OC(=O)c1nc([*])c([*])c([*])c1([*])",
    smiles = "COC(=O)c1ncccc1"
)
frag_164 = Molecule(
    name = "164",
    smarts = "O=C([*])N([*])c1nc([*])c([*])c([*])c1([*])",
    smiles = "O=C[NH]c1ncccc1"
)
frag_165 = Molecule(
    name = "165",
    smarts = "C([*])([*])([*])C(=O)N([*])c1nc([*])c([*])c([*])c1([*])",
    smiles = "CC(=O)[NH]c1ncccc1"
)
frag_166 = Molecule(
    name = "166",
    smarts = "[O-][n+]1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "[O-][n+]1ccccc1"  
)
frag_167 = Molecule(
    name = "167",
    smarts = "c1([*])c([*])nnc([*])c1([*])",
    smiles = "c1cnncc1"
)
frag_168 = Molecule(
    name = "168",
    smarts = "c1([*])c([*])nc([*])nc1([*])",
    smiles = "c1cncnc1"
)
frag_169 = Molecule(
    name = "169",
    smarts = "c1([*])c([*])nc([*])c([*])n1",
    smiles = "c1cnccn1"
)
frag_170 = Molecule(
    name = "170",
    smarts = "c1([*])c([*])nnnc1([*])",
    smiles = "c1cnnnc1"
)
frag_171 = Molecule(
    name = "171",
    smarts = "c1([*])nc([*])nnc1([*])",
    smiles = "c1ncnnc1"
)
frag_172 = Molecule(
    name = "172",
    smarts = "c1([*])nc([*])nc([*])n1",
    smiles = "c1ncncn1"
)
frag_173 = Molecule(
    name = "173",
    smarts = "N([*])([*])c1nc([*])c([*])c([*])n1",
    smiles = "[NH2]c1ncccn1"
)
frag_174 = Molecule(
    name = "174",
    smarts = "O=c1c([*])c([*])c([*])c([*])n([*])1",
    smiles = "O=c1cccc[nH]1"
)
frag_175 = Molecule(
    name = "175",
    smarts = "C([*])([*])([*])n1c([*])c([*])c([*])c([*])c1=O",
    smiles = "Cn1ccccc1=O"
)
frag_176 = Molecule(
    name = "176",
    smarts = "O=c1c([*])c([*])n([*])c([*])c1([*])",
    smiles = "O=c1cc[nH]cc1"
)
frag_177 = Molecule(
    name = "177",
    smarts = "C([*])([*])([*])n1c([*])c([*])c(=O)c([*])c1([*])",
    smiles = "Cn1ccc(=O)cc1"
)
frag_178 = Molecule(
    name = "178",
    smarts = "O=c1n([*])c([*])c([*])c([*])n1",
    smiles = "O=c1[nH]cccn1"
)
frag_179 = Molecule(
    name = "179",
    smarts = "O=c1c([*])c([*])n([*])c([*])n1",
    smiles = "O=c1cc[nH]cn1"
)
frag_180 = Molecule(
    name = "180",
    smarts = "O=c1c([*])c([*])c([*])c([*])o1",
    smiles = "O=c1cccco1"
)
frag_181 = Molecule(
    name = "181",
    smarts = "O=c1c([*])c([*])oc([*])c1([*])",
    smiles = "O=c1ccocc1"
)
frag_182 = Molecule(
    name = "182",
    smarts = "c1([*])c([*])c([*])c2c(c1([*]))C([*])([*])C([*])([*])N([*])2",
    smiles = "c1ccc2c(c1)CC[NH]2"
)
frag_183 = Molecule(
    name = "183",
    smarts = "C([*])([*])([*])N1c2c([*])c([*])c([*])c([*])c2C([*])([*])C([*])([*])C1([*])([*])",
    smiles = "CN1c2ccccc2CCC1"
)
frag_184 = Molecule(
    name = "184",
    smarts = "O=C1C([*])([*])c2c([*])c([*])c([*])c([*])c2N([*])1",
    smiles = "O=C1Cc2ccccc2[NH]1"
)
frag_185 = Molecule(
    name = "185",
    smarts = "C([*])([*])([*])N1C(=O)C([*])([*])c2c([*])c([*])c([*])c([*])c21",
    smiles = "CN1C(=O)Cc2ccccc21"
)
frag_186 = Molecule(
    name = "186",
    smarts = "O=C1N([*])C([*])([*])c2c([*])c([*])c([*])c([*])c21",
    smiles = "O=C1[NH]Cc2ccccc21"
)
frag_187 = Molecule(
    name = "187",
    smarts = "C([*])([*])([*])N1C([*])([*])c2c([*])c([*])c([*])c([*])c2C1=O",
    smiles = "CN1Cc2ccccc2C1=O"
)
frag_188 = Molecule(
    name = "188",
    smarts = "O=C1N([*])C([*])([*])C([*])([*])c2c([*])c([*])c([*])c([*])c21",
    smiles = "O=C1[NH]CCc2ccccc21"
)
frag_189 = Molecule(
    name = "189",
    smarts = "c1([*])c([*])c([*])c2c(c1([*]))OC([*])([*])O2",
    smiles = "c1ccc2c(c1)OCO2"
)
frag_190 = Molecule(
    name = "190",
    smarts = "C1([*])([*])C([*])([*])C([*])([*])c2c(C1([*])([*]))c([*])c([*])n([*])2",
    smiles = "C1CCc2c(C1)cc[nH]2"
)
frag_191 = Molecule(
    name = "191",
    smarts = "C([*])([*])([*])n1c([*])c([*])c(=O)n([*])c1=O",
    smiles = "Cn1ccc(=O)[nH]c1=O"
)
frag_192 = Molecule(
    name = "192",
    smarts = "C([*])([*])([*])c1c([*])n(C([*])([*])([*]))c(=O)n([*])c1=O",
    smiles = "Cc1cn(C)c(=O)[nH]c1=O"
)
frag_193 = Molecule(
    name = "193",
    smarts = "C([*])([*])([*])n1c([*])c([*])c(N([*])([*]))nc1=O",
    smiles = "Cn1ccc([NH2])nc1=O"
)
frag_194 = Molecule(
    name = "194",
    smarts = "C([*])([*])([*])n1c([*])nc2c(N([*])([*]))nc([*])nc21",
    smiles = "Cn1cnc2c([NH2])ncnc21"
)
frag_195 = Molecule(
    name = "195",
    smarts = "C([*])([*])([*])n1c([*])nc2c1nc(N([*])([*]))n([*])c2=O",
    smiles = "Cn1cnc2c1nc([NH2])[nH]c2=O"
)
frag_196 = Molecule(
    name = "196",
    smarts = "c1([*])c([*])c([*])c(-c2c([*])c([*])c([*])c([*])c2([*]))c([*])c1([*])",
    smiles = "c1ccc(-c2ccccc2)cc1"
)
frag_197 = Molecule(
    name = "197",
    smarts = "c1([*])c([*])c([*])c(-c2c([*])c([*])c([*])n([*])2)c([*])c1([*])",
    smiles = "c1ccc(-c2ccc[nH]2)cc1"
)
frag_198 = Molecule(
    name = "198",
    smarts = "C([*])([*])([*])n1c([*])c([*])c([*])c1-c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "Cn1cccc1-c1ccccc1"
)
frag_199 = Molecule(
    name = "199",
    smarts = "c1([*])c([*])c([*])c(-c2n([*])c([*])c([*])n2)c([*])c1([*])",
    smiles = "c1ccc(-c2[nH]ccn2)cc1"
)
frag_200 = Molecule(
    name = "200",
    smarts = "c1([*])c([*])c([*])c(n2c([*])c([*])c([*])n2)c([*])c1([*])",
    smiles = "c1ccc(n2cccn2)cc1"
)
frag_201 = Molecule(
    name = "201",
    smarts = "c1([*])c([*])c([*])c(n2c([*])c([*])nn2)c([*])c1([*])",
    smiles = "c1ccc(n2ccnn2)cc1"
)
frag_202 = Molecule(
    name = "202",
    smarts = "c1([*])c([*])c([*])c(-c2c([*])c([*])c([*])o2)c([*])c1([*])",
    smiles = "c1ccc(-c2ccco2)cc1"
)
frag_203 = Molecule(
    name = "203",
    smarts = "c1([*])c([*])c([*])c(-c2nc([*])c([*])o2)c([*])c1([*])",
    smiles = "c1ccc(-c2ncco2)cc1"
)
frag_204 = Molecule(
    name = "204",
    smarts = "c1([*])c([*])c([*])c(-c2c([*])c([*])c([*])s2)c([*])c1([*])",
    smiles = "c1ccc(-c2cccs2)cc1"
)
frag_205 = Molecule(
    name = "205",
    smarts = "C([*])([*])([*])c1c([*])c([*])sc1-c2c([*])c([*])c([*])c([*])c2([*])",
    smiles = "Cc1ccsc1-c2ccccc2"
)
frag_206 = Molecule(
    name = "206",
    smarts = "c1([*])c([*])oc(c1([*]))-c1nc([*])nn([*])1",
    smiles = "c1coc(c1)-c1ncn[nH]1"
)
frag_207 = Molecule(
    name = "207",
    smarts = "C([*])([*])([*])n1c([*])nnc1-c1c([*])c([*])c([*])s1",
    smiles = "Cn1cnnc1-c1cccs1"
)
frag_208 = Molecule(
    name = "208",
    smarts = "c1([*])c([*])c([*])c(C([*])([*])c2c([*])c([*])c([*])c([*])c2([*]))c([*])c1([*])",
    smiles = "c1ccc(Cc2ccccc2)cc1"
)
frag_209 = Molecule(
    name = "209",
    smarts = "c1([*])c([*])c([*])c(N([*])c2c([*])c([*])c([*])c([*])c2([*]))c([*])c1([*])",
    smiles = "c1ccc([NH]c2ccccc2)cc1"
)
frag_210 = Molecule(
    name = "210",
    smarts = "c1([*])c([*])c([*])c(Oc2c([*])c([*])c([*])c([*])c2([*]))c([*])c1([*])",
    smiles = "c1ccc(Oc2ccccc2)cc1"
)
frag_211 = Molecule(
    name = "211",
    smarts = "O=C(c1c([*])c([*])c([*])c([*])c1([*]))c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "O=C(c1ccccc1)c1ccccc1"
)
frag_212 = Molecule(
    name = "212",
    smarts = "c1([*])c([*])c([*])c2c(c1([*]))c([*])c([*])c([*])c2([*])",
    smiles = "c1ccc2c(c1)cccc2"
)
frag_213 = Molecule(
    name = "213",
    smarts = "c1([*])c([*])nc2c(c1([*]))c([*])c([*])c([*])c2([*])",
    smiles = "c1cnc2c(c1)cccc2"
)
frag_214 = Molecule(
    name = "214",
    smarts = "c1([*])nc([*])c2c(c1([*]))c([*])c([*])c([*])c2([*])",
    smiles = "c1ncc2c(c1)cccc2"
)
frag_215 = Molecule(
    name = "215",
    smarts = "n1c([*])nc2c(c1([*]))c([*])c([*])c([*])c2([*])",
    smiles = "n1cnc2c(c1)cccc2"
)
frag_216 = Molecule(
    name = "216",
    smarts = "c1([*])c([*])nc2c(n1)c([*])c([*])c([*])c2([*])",
    smiles = "c1cnc2c(n1)cccc2"
)
frag_217 = Molecule(
    name = "217",
    smarts = "c1([*])c([*])nc2c(c1([*]))c([*])nc([*])c2([*])",
    smiles = "c1cnc2c(c1)cncc2"
)
frag_218 = Molecule(
    name = "218",
    smarts = "c1([*])c([*])nc2c(c1([*]))c([*])c([*])c([*])n2",
    smiles = "c1cnc2c(c1)cccn2"
)
frag_219 = Molecule(
    name = "219",
    smarts = "c1([*])c([*])nc2c(c1([*]))nc([*])nc2([*])",
    smiles = "c1cnc2c(c1)ncnc2"
)
frag_220 = Molecule(
    name = "220",
    smarts = "c1([*])c([*])nc2c(n1)c([*])nc([*])n2",
    smiles = "c1cnc2c(n1)cncn2"
)
frag_221 = Molecule(
    name = "221",
    smarts = "c1([*])c([*])c([*])c2c(c1([*]))n([*])c([*])c2([*])",
    smiles = "c1ccc2c(c1)[nH]cc2"
)
frag_222 = Molecule(
    name = "222",
    smarts = "c1([*])c([*])c([*])c2c(c1([*]))n([*])nc2([*])",
    smiles = "c1ccc2c(c1)[nH]nc2"
)
frag_223 = Molecule(
    name = "223",
    smarts = "c1([*])c([*])c([*])c2c(c1([*]))c([*])n([*])n2",
    smiles = "c1ccc2c(c1)c[nH]n2"
)
frag_224 = Molecule(
    name = "224",
    smarts = "c1([*])c([*])c([*])c2c(c1([*]))nc([*])n([*])2",
    smiles = "c1ccc2c(c1)nc[nH]2"
)
frag_225 = Molecule(
    name = "225",
    smarts = "c1([*])c([*])c([*])c2n(c1([*]))c([*])c([*])n2",
    smiles = "c1ccc2n(c1)ccn2"
)
frag_226 = Molecule(
    name = "226",
    smarts = "c1([*])c([*])nc2c(c1([*]))n([*])c([*])c2([*])",
    smiles = "c1cnc2c(c1)[nH]cc2"
)
frag_227 = Molecule(
    name = "227",
    smarts = "n1c([*])c([*])c2c(c1([*]))c([*])c([*])n([*])2",
    smiles = "n1ccc2c(c1)cc[nH]2"
)
frag_228 = Molecule(
    name = "228",
    smarts = "c1([*])c([*])nc2c(c1([*]))c([*])c([*])n([*])2",
    smiles = "c1cnc2c(c1)cc[nH]2"
)
frag_229 = Molecule(
    name = "229",
    smarts = "c1([*])nn2c([*])c([*])c([*])c([*])c2c1([*])",
    smiles = "c1nn2ccccc2c1"
)
frag_230 = Molecule(
    name = "230",
    smarts = "c1([*])c([*])c([*])c2n(c1([*]))c([*])nc2([*])",
    smiles = "c1ccc2n(c1)cnc2"
)
frag_231 = Molecule(
    name = "231",
    smarts = "n1c([*])c([*])c2c(c1([*]))nc([*])n([*])2",
    smiles = "n1ccc2c(c1)nc[nH]2"
)
frag_232 = Molecule(
    name = "232",
    smarts = "n1c([*])nc2c(c1([*]))c([*])c([*])n([*])2",
    smiles = "n1cnc2c(c1)cc[nH]2"
)
frag_233 = Molecule(
    name = "233",
    smarts = "c1([*])c([*])nc2n(c1([*]))c([*])c([*])n2",
    smiles = "c1cnc2n(c1)ccn2"
)
frag_234 = Molecule(
    name = "234",
    smarts = "n1c([*])nc2c(c1([*]))n([*])c([*])n2",
    smiles = "n1cnc2c(c1)[nH]cn2"
)
frag_235 = Molecule(
    name = "235",
    smarts = "c1([*])c([*])c([*])c2c(c1([*]))oc([*])c2([*])",
    smiles = "c1ccc2c(c1)occ2"
)
frag_236 = Molecule(
    name = "236",
    smarts = "c1([*])c([*])c([*])c2c(c1([*]))oc([*])n2",
    smiles = "c1ccc2c(c1)ocn2"
)
frag_237 = Molecule(
    name = "237",
    smarts = "n1c([*])nc2c(c1([*]))oc([*])c2([*])",
    smiles = "n1cnc2c(c1)occ2"
)
frag_238 = Molecule(
    name = "238",
    smarts = "c1([*])c([*])c([*])c2c(c1([*]))sc([*])c2([*])",
    smiles = "c1ccc2c(c1)scc2"
)
frag_239 = Molecule(
    name = "239",
    smarts = "c1([*])c([*])c([*])c2c(c1([*]))sc([*])n2",
    smiles = "c1ccc2c(c1)scn2"
)
frag_240 = Molecule(
    name = "240",
    smarts = "C([*])([*])([*])n1c([*])nc2c1c([*])c([*])c([*])c2([*])",
    smiles = "Cn1cnc2c1cccc2"
)
frag_241 = Molecule(
    name = "241",
    smarts = "C([*])([*])([*])n1c([*])nc2c1c([*])c([*])nc2([*])",
    smiles = "Cn1cnc2c1ccnc2"
)
frag_242 = Molecule(
    name = "242",
    smarts = "C([*])([*])([*])n1c([*])nc2c1ncn([*])c2=O",
    smiles = "Cn1cnc2c1nc[nH]c2=O"
)
frag_243 = Molecule(
    name = "243",
    smarts = "C([*])([*])([*])n1c([*])nc2c(c1(=O))c([*])c([*])n([*])2",
    smiles = "Cn1cnc2c(c1(=O))cc[nH]2"
)
frag_244 = Molecule(
    name = "244",
    smarts = "O=c1n([*])c([*])c([*])c2c1c([*])c([*])c([*])c2([*])",
    smiles = "O=c1[nH]ccc2c1cccc2"
)
frag_245 = Molecule(
    name = "245",
    smarts = "C([*])([*])([*])n1c([*])c([*])c2c(c1(=O))c([*])c([*])c([*])c2([*])",
    smiles = "Cn1ccc2c(c1(=O))cccc2"
)
frag_246 = Molecule(
    name = "246",
    smarts = "O=c1nc([*])n([*])c2c1nc([*])c([*])c2([*])",
    smiles = "O=c1nc[nH]c2c1nccc2"
)
frag_247 = Molecule(
    name = "247",
    smarts = "O=c1c([*])c([*])c2c(n([*])1)c([*])c([*])c([*])c2([*])",
    smiles = "O=c1ccc2c([nH]1)cccc2"
)
frag_248 = Molecule(
    name = "248",
    smarts = "O=c1c([*])c([*])c2c(o1)c([*])c([*])c([*])c2([*])",
    smiles = "O=c1ccc2c(o1)cccc2"
)
frag_249 = Molecule(
    name = "249",
    smarts = "C([*])([*])([*])C(=O)[O-]",
    smiles = "CC(=O)[O-]"
)
frag_250 = Molecule(
    name = "250",
    smarts = "C([*])([*])([*])C([*])(C([*])([*])([*]))C(=O)[O-]",
    smiles = "CC(C)C(=O)[O-]"
)
frag_251 = Molecule(
    name = "251",
    smarts = "C([*])([*])([*])C(C([*])([*])([*]))(C([*])([*])([*]))C(=O)[O-]",
    smiles = "CC(C)(C)C(=O)[O-]"
)
frag_252 = Molecule(
    name = "252",
    smarts = "C1([*])([*])C([*])([*])C1([*])C(=O)[O-]",
    smiles = "O=C([O-])C1CC1"
)
frag_253 = Molecule(
    name = "253",
    smarts = "O=C([O-])c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "O=C([O-])c1ccccc1"
)
frag_254 = Molecule(
    name = "254",
    smarts = "O=C([O-])c1c([*])c([*])c([*])c([*])c1F",
    smiles = "O=C([O-])c1ccccc1F"
)
frag_255 = Molecule(
    name = "255",
    smarts = "O=C([O-])c1c([*])c([*])c([*])c([*])c1Cl",
    smiles = "O=C([O-])c1ccccc1Cl"
)
frag_256 = Molecule(
    name = "256",
    smarts = "O=C([O-])c1c([*])c([*])c([*])c([*])n1",
    smiles = "O=C([O-])c1ccccn1"
)
frag_257 = Molecule(
    name = "257",
    smarts = "O=C([O-])c1nc([*])c([*])s1",
    smiles = "O=C([O-])c1nccs1"
)
frag_258 = Molecule(
    name = "258",
    smarts = "C([*])([*])([*])c1n[n-]nn1",
    smiles = "Cc1n[n-]nn1"
)
frag_259 = Molecule(
    name = "259",
    smarts = "C([*])([*])([*])c1c([*])snc1[O-]",
    smiles = "Cc1csnc1[O-]"
)
frag_260 = Molecule(
    name = "260",
    smarts = "C([*])([*])([*])c1c([*])onc1[O-]",
    smiles = "Cc1conc1[O-]"
)
frag_261 = Molecule(
    name = "261",
    smarts = "C([*])([*])([*])c1noc([O-])n1",
    smiles = "Cc1noc([O-])n1"
)
frag_262 = Molecule(
    name = "262",
    smarts = "C([*])([*])([*])c1nsc([O-])n1",
    smiles = "Cc1nsc([O-])n1"
)
frag_263 = Molecule(
    name = "263",
    smarts = "O=c1c([*])c([O-])c([*])([*])o1",
    # NOTE: the aromatic SMILES "O=c1cc([O-])co1" used here originally is not
    # kekulizable by RDKit (an [O-] substituent on an aromatic C cannot be
    # assigned a valid Kekule structure), so construct_template() silently
    # returned None for this fragment. Written in the kekulized form instead.
    smiles = "O=C1C=C([O-])CO1"
)
frag_264 = Molecule(
    name = "264",
    smarts = "O=C([O-])c1c([*])c([*])c([*])c([*])c1-c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "O=C([O-])c1ccccc1-c1ccccc1"
)
frag_265 = Molecule(
    name = "265",
    smarts = "C([*])([*])([*])C(=O)[N-]S(=O)(=O)C([*])([*])([*])",
    smiles = "CC(=O)[N-]S(=O)(=O)C"
)
frag_266 = Molecule(
    name = "266",
    smarts = "[O-]c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "[O-]c1ccccc1"
)
frag_267 = Molecule(
    name = "267",
    smarts = "[O-]c1c([*])c([*])c([*])c([*])c1F",
    smiles = "[O-]c1ccccc1F"
)
frag_268 = Molecule(
    name = "268",
    smarts = "[O-]c1c([*])c([*])c([*])c([*])c1Cl",
    smiles = "[O-]c1ccccc1Cl"
)
frag_269 = Molecule(
    name = "269",
    smarts = "[O-]c1c([*])c([*])c([*])c([*])n1",
    smiles = "[O-]c1ccccn1"
)
frag_270 = Molecule(
    name = "270",
    smarts = "[O-]c1c([*])c([*])c([*])c2c1c([*])c([*])c([*])c2([*])",
    smiles = "[O-]c1cccc2c1cccc2"
)
frag_271 = Molecule(
    name = "271",
    smarts = "[S-]c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "[S-]c1ccccc1"
)
frag_272 = Molecule(
    name = "272",
    smarts = "C([*])([*])([*])[S-]",
    smiles = "C[S-]"
)
frag_273 = Molecule(
    name = "273",
    smarts = "C([*])([*])([*])C(=O)N([*])[O-]",
    smiles = "CC(=O)[NH][O-]"   # 添加极性氢
)
frag_274 = Molecule(
    name = "274",
    smarts = "C([*])(C(F)(F)F)(C(F)(F)F)[O-]",
    smiles = "C(C(F)(F)F)(C(F)(F)F)[O-]"
)
frag_275 = Molecule(
    name = "275",
    smarts = "C([*])([*])([*])S(=O)(=O)[O-]",
    smiles = "CS(=O)(=O)[O-]"
)
frag_276 = Molecule(
    name = "276",
    smarts = "[O-]S(=O)(=O)[O-]",
    smiles = "O=S(=O)([O-])[O-]"
)
frag_277 = Molecule(
    name = "277",
    smarts = "C([*])([*])([*])P(=O)([O-])OC([*])([*])([*])",
    smiles = "COP(=O)(C)[O-]"
)
frag_278 = Molecule(
    name = "278",
    smarts = "C([*])([*])([*])P(=O)(O([*]))[O-]",
    smiles = "CP(=O)([OH])[O-]" # 添加极性氢
)
frag_279 = Molecule(
    name = "279",
    smarts = "C([*])([*])([*])OP(=O)([O-])OC([*])([*])([*])",
    smiles = "COP(=O)([O-])OC"
)
frag_280 = Molecule(
    name = "280",
    smarts = "[Cl-]",
    smiles = "[Cl-]"
)
frag_281 = Molecule(
    name = "281",
    smarts = "C([*])([*])([*])[N+]([*])([*])([*])",
    smiles = "C[NH3+]"
)
frag_282 = Molecule(
    name = "282",
    smarts = "C([*])([*])([*])[N+]([*])([*])C([*])([*])([*])",
    smiles = "C[NH2+]C"
)
frag_283 = Molecule(
    name = "283",
    smarts = "C([*])([*])([*])[N+](C([*])([*])([*]))(C([*])([*])([*]))C([*])([*])([*])",
    smiles = "C[N+](C)(C)C"
)
frag_284 = Molecule(
    name = "284",
    smarts = "C1([*])([*])C([*])([*])C([*])([*])[N+]([*])([*])C([*])([*])C1([*])([*])",
    smiles = "C1CC[NH2+]CC1"
)
frag_285 = Molecule(
    name = "285",
    smarts = "C1([*])([*])C([*])([*])C([*])([*])[N+]([*])([*])C1([*])([*])",
    smiles = "C1CC[NH2+]C1"
)
frag_286 = Molecule(
    name = "286",
    smarts = "C1([*])([*])C([*])([*])[N+]([*])([*])C1([*])([*])",
    smiles = "C1C[NH2+]C1"
)
frag_287 = Molecule(
    name = "287",
    smarts = "C([*])([*])([*])[N+]([*])1C([*])([*])C([*])([*])C([*])([*])C([*])([*])C1([*])([*])",
    smiles = "C[NH+]1CCCCC1"
)
frag_288 = Molecule(
    name = "288",
    smarts = "C([*])([*])([*])[N+]([*])1C([*])([*])C([*])([*])C([*])([*])C1([*])([*])",
    smiles = "C[NH+]1CCCC1"
)
frag_289 = Molecule(
    name = "289",
    smarts = "C([*])([*])([*])[N+]([*])1C([*])([*])C([*])([*])C1([*])([*])",
    smiles = "C[NH+]1CCC1"
)
frag_290 = Molecule(
    name = "290",
    smarts = "c1([*])c([*])[n+]([*])c([*])n([*])1",
    smiles = "c1c[nH+]c[nH]1"
)
frag_291 = Molecule(
    name = "291",
    smarts = "C([*])([*])([*])N([*])c1n([*])c([*])c([*])[n+]([*])1",
    smiles = "CNc1[nH]cc[nH+]1"
)
frag_292 = Molecule(
    name = "292",
    smarts = "c2([*])c([*])c([*])c(c1n([*])c([*])c([*])[n+]([*])1)c([*])c2([*])",
    smiles = "c1ccc(c2[nH]cc[nH+]2)cc1"
)
frag_293 = Molecule(
    name = "293",
    smarts = "c2([*])c([*])c([*])c1[n+]([*])c([*])n([*])c1c2([*])",
    smiles = "c1ccc2[nH+]c[nH]c2c1"
)
frag_294 = Molecule(
    name = "294",
    smarts = "C([*])([*])([*])c2n([*])c1c([*])c([*])c([*])c([*])c1[n+]([*])2",
    smiles = "Cc1[nH]c2ccccc2[nH+]1"
)
frag_295 = Molecule(
    name = "295",
    smarts = "c2([*])c([*])c([*])n1c([*])c([*])[n+]([*])c1c2([*])",
    smiles = "c1ccn2cc[nH+]c2c1"
)
frag_296 = Molecule(
    name = "296",
    smarts = "N([*])([*])C(=[N+]([*])([*]))c1c([*])c([*])c([*])c([*])c1([*])",
    smiles = "[NH2]C(=[NH2+])c1ccccc1"
)
frag_297 = Molecule(
    name = "297",
    smarts = "c1([*])c([*])c([*])[n+]([*])c([*])c1([*])",
    smiles = "c1cc[nH+]cc1"
)
frag_298 = Molecule(
    name = "298",
    smarts = "C([*])([*])([*])c1c([*])c([*])c([*])c([*])[n+]([*])1",
    smiles = "Cc1cccc[nH+]1"
)
frag_299 = Molecule(
    name = "299",
    smarts = "N([*])([*])c1c([*])c([*])c([*])c([*])[n+]([*])1",
    smiles = "[NH2]c1cccc[nH+]1"
)
frag_300 = Molecule(
    name = "300",
    smarts = "N([*])([*])c1c([*])c([*])[n+]([*])c([*])c1([*])",
    smiles = "[NH2]c1cc[nH+]cc1"
)
frag_301 = Molecule(
    name = "301",
    smarts = "C([*])([*])([*])N([*])c1c([*])c([*])[n+]([*])c([*])c1([*])",
    smiles = "C[NH]c1cc[nH+]cc1"
)
frag_302 = Molecule(
    name = "302",
    smarts = "c2([*])c([*])c([*])c1[n+]([*])c([*])c([*])c([*])c1c2([*])",
    smiles = "c1ccc2[nH+]cccc2c1"
)
frag_303 = Molecule(
    name = "303",
    smarts = "[Na+]",
    smiles = "[Na+]"
)
frag_304 = Molecule(
    name = "304",
    smarts = "c1([*])nnnn([*])1",
    smiles = "c1nnn[nH]1"
)
frag_305 = Molecule(
    name = "305",
    smarts = "c2([*])nc1c([*])n([*])c([*])nc1n2",
    smiles = "c1nc2c[nH]cnc2n1"
)
frag_306 = Molecule(
    name = "306",
    smarts = "O=c1c([*])c([*])n([*])c2c([*])c([*])c([*])c([*])c12",
    smiles = "O=c1cc[nH]c2c1cccc2"
)
frag_307 = Molecule(
    name = "307",
    smarts = "C1([*])([*])C([*])([*])C([*])([*])N([*])C([*])([*])C1([*])([*])",
    smiles = "C1CC[NH]CC1"
)
frag_308 = Molecule(
    name = "308",
    smarts = "C1([*])([*])C([*])([*])N([*])C([*])([*])C([*])([*])N([*])1",
    smiles = "C1C[NH]CC[NH]1"
)
frag_309 = Molecule(
    name = "309",
    smarts = "C1([*])([*])C([*])([*])C([*])([*])OC([*])([*])C1([*])([*])",
    smiles = "C1CCOCC1"
)
frag_310 = Molecule(
    name = "310",
    smarts = "C1([*])([*])C([*])([*])C([*])([*])C([*])([*])C1([*])([*])",
    smiles = "C1CCCC1"
)
frag_311 = Molecule(
    name = "311",
    smarts = "C1([*])([*])C([*])([*])OC([*])([*])C([*])([*])N([*])1",
    smiles = "C1COCC[NH]1"
)
frag_312 = Molecule(
    name = "312",
    smarts = "C1([*])=NC([*])([*])C([*])([*])N([*])1",
    smiles = "C1CN=C[NH]1"
)
frag_313 = Molecule(
    name = "313",
    smarts = "C1([*])=C([*])N([*])C([*])=C([*])C1([*])([*])",
    smiles = "C1C=C[NH]C=C1"
)
frag_314 = Molecule(
    name = "314",
    smarts = "C1([*])=C([*])C([*])([*])C([*])([*])C([*])([*])C1([*])([*])",
    smiles = "C1CCC=CC1"
)

# Create the PROBE_FRAGMENTS dictionary with all defined molecules
PROBE_FRAGMENTS = {
    "ACEM": ACEM,
    "ACET": ACET,
    "ETAM": ETAM,
    "ETOH": ETOH,
    "ETSH": ETSH,
    "HOH": HOH,
    "MBZ": MBZ,
    "MGDM": MGDM,
    "MIMD": MIMD,
    "MIME": MIME,
    "MIMM": MIMM,
    "MIND": MIND,
    "MSM": MSM,
    "MTYR": MTYR,
    "N1PA": N1PA,
    "NMA": NMA,
    "PRPA": PRPA,
}

REF_FRAGMENTS = {}
for i in range(1, 315):
    _var = globals()[f"frag_{i}"]
    REF_FRAGMENTS[f"frag_{i}"] = _var