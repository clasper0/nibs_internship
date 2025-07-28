# coding=utf-8
from pymatgen.io.cif import CifParser   # 更完整地读取 cif 文件信息（包括键长等）
from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.core.periodic_table import Element    # 获取元素半径
from rdkit import Chem
from rdkit.Chem import AllChem
from all_frags import Molecule
from all_frags import PROBE_FRAGMENTS, REF_FRAGMENTS


def get_struct_info(cif_path, cod_id):
    """
    从 CIF 文件中提取必要的原子坐标信息、成键信息，返回 atom_dict, bond_dict
    """
    file_path = cif_path + f"{cod_id[0]}/{cod_id[1:3]}/{cod_id[3:5]}/{cod_id}.cif"
    parser = CifParser(file_path)
    structure = parser.parse_structures(primitive=True)[0]
    cif_data = parser.as_dict()[cod_id]

    # atom_dict -> {index: [position: label, symbol]}
    atom_dict = {}
    for i, site in enumerate(structure):
        dict = site.as_dict()
        atom_dict[i] = [site.coords, dict['label'], dict['species'][0]['element']]
    # print(atom_dict[440])

    # 获取成键原子信息
    label1 = cif_data['_geom_bond_atom_site_label_1']
    label2 = cif_data['_geom_bond_atom_site_label_2']
    labels = list(zip(label1, label2))
    # 获取所有原子的相邻成键原子（但在 CIF 文件中通常不包括 H ）
    bond_label_dict = {}
    for label1, label2 in labels:
        if label1 not in bond_label_dict:
            bond_label_dict[label1] = []
        bond_label_dict[label1].append(label2)
        if label2 not in bond_label_dict:
            bond_label_dict[label2] = []
        bond_label_dict[label2].append(label1)
    # print(bond_label_dict)
    # 由于成键两端标签和两端索引是一对多的关系，所以必须把标签字典转换成索引字典 (bond_dict)
    bond_dict = {}
    for label1, label2s in bond_label_dict.items():
        bond_indexes_1 = [index for index, info in atom_dict.items() if label1 in info]
        for bond_index_1 in bond_indexes_1:
            bond_dict[bond_index_1] = set()
            
            # 建立嵌套子对象，存储候选 bond_indexes_2 的各个索引及坐标
            site_1 = structure[bond_index_1] # 成键原子 1 的 site, 用于计算距离，所用方法为: site_1.distance(site_2)
            for label2 in label2s:  # 获得候选成键原子的索引集合
                candidates_index_2 = [index for index, info in atom_dict.items() if label2 in info]
                distances = []
                for candidate_index_2 in candidates_index_2:
                    site_2 = structure[candidate_index_2]
                    distances.append(site_1.distance(site_2))
                bond_index_2 = candidates_index_2[distances.index(min(distances))]  # 取出距离最近的索引作为成键原子
                bond_dict[bond_index_1].add(bond_index_2)
    # print(bond_dict)

    return atom_dict, bond_dict

def construct_mol(atom_dict, bond_dict):
    """
    根据 atom_dict, bond_dict, 从头搭建 Mol 分子
    """
    # 从 atom_dict 中提取元素符号及坐标信息
    positions, symbols = [], []
    for key, item in atom_dict.items():
        positions.append(item[0])
        symbols.append(item[2])
    # print(positions[0][0])
    # print(symbols)

    # 创建空对象
    mol = Chem.Mol()
    mol = Chem.EditableMol(mol)

    # 添加原子
    for symbol in symbols:
        mol.AddAtom(Chem.Atom(symbol))

    # 添加键连信息
    for i, js in bond_dict.items():
        for j in js:
            if i < j:   # 防止双向映射导致成键失败
                mol.AddBond(i, j, order=Chem.BondType.SINGLE)   # 由于缺少键类型等，默认按照单键处理
    mol = mol.GetMol()

    # 添加原子坐标
    conf = Chem.Conformer(len(symbols))
    for i, pos in enumerate(positions):
        x, y, z = pos
        conf.SetAtomPosition(i, (x, y, z))
    mol.AddConformer(conf)

    return mol


def search(mol, substruct1, substruct2):
    """
    进行子结构查询匹配
    """
    match1 = mol.HasSubstructMatch(substruct1)
    match2 = mol.HasSubstructMatch(substruct2)
    if match1 and match2:
        return True
    else:
        return False


if __name__ == '__main__':
    cif_path = '/home/clasper/2507/cif/'
    cod_id = '1004005'
    print(f"Loading file from {cod_id}.cif")
    atom_dict, bond_dict = get_struct_info(cif_path, cod_id)

    print(f"Constructing mol for {cod_id}.cif")
    mol = construct_mol(atom_dict, bond_dict)

    patt1 = "C1CCCCC1"
    patt2 = "CCC"

    substruct1 = Chem.MolFromSmarts("C1CCCCC1")
    substruct2 = Chem.MolFromSmarts("CCC")

    # mol, substruct1, substruct2 的结构都有待优化
    print(f"Searching cod_id {cod_id} for {patt1} and {patt2}...")
    answer = search(mol, substruct1, substruct2)

    print(answer)
