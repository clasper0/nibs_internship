# coding=utf-8
"""
1. 读取 cif 文件
2. 将 cif 文件转换为 RDKit 分子对象
3. 使用 RDKit 对分子进行分析，匹配子结构
4. 将提取到的匹配分子对进行标准化处理（参考 csd_workflow.py）
5. 将标准化后的对象保存为 sdf 文件

csd_workflow.py
1. 读取预定义 --frag1 & --frag2，并定义各自的smart, smiles, heavy_atom及其数目
2. 构建模板片段对
3. 搜索 CSD 数据库，提取匹配的分子片段对
4. 将匹配到的 hits 转换为 mol2 文件
5. 矫正氢原子
6. 对齐
7. 去除重复分子
8. 保存对齐后/未能矫正的分子为 sdf 文件
"""
from pymatgen.io.cif import CifParser   # 更完整地读取 cif 文件信息（包括键长等）
from rdkit import Chem
from rdkit.Chem import AllChem
from all_frags import Molecule
from all_frags import PROBE_FRAGMENTS, REF_FRAGMENTS
import numpy as np
from scipy.spatial.distance import cdist


def SimplySub(mol):
    """
    Change all bonds to single for the substructures. 
    Return mol -> mol. 
    """
    for b in mol.GetBonds():
        if b.GetIsAromatic():
            continue
        b.SetBondType(Chem.BondType.SINGLE)
    
    return mol


def SimplifyMol(mol):
    """Change all bonds to single and discharge/dearomatize all atoms.
    The molecule is modified in-place (no copy is made).
    from https://github.com/oddt/oddt/blob/master/oddt/toolkits/extras/rdkit/fixer.py
    """
    for b in mol.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in mol.GetAtoms():
        a.SetFormalCharge(0)
        a.SetIsAromatic(False)

    return mol


def get_heavy_indices_from_smarts(smarts):
    """
    Parses a SMARTS string and returns a list of indices for atoms that are
    explicitly defined as heavy atoms (i.e., not '*' or '[H]').
    """
    mol = Chem.MolFromSmarts(smarts)
    if not mol:
        print(f"Warning: Could not parse SMARTS string: {smarts}", flush=True)
        return []
    
    heavy_atom_indices = []
    for atom in mol.GetAtoms():
        # Get the SMARTS representation of the single atom
        atom_smarts = atom.GetSmarts()
        # We consider it a core heavy atom if it's not a generic wildcard or explicit hydrogen
        if atom_smarts not in ['*', '[H]', '[#1]']:
            heavy_atom_indices.append(atom.GetIdx())

    return heavy_atom_indices


def construct_template(smiles):
    template = Chem.MolFromSmiles(smiles)
    template = Chem.AddHs(template)
    AllChem.EmbedMolecule(template)
    order = []
    heavy_count = 0
    for a in template.GetAtoms():
        if a.GetAtomicNum() == 1:
            continue
        heavy_count += 1
        order.append(a.GetIdx())
        for n in a.GetNeighbors():
            if n.GetAtomicNum() == 1:
                order.append(n.GetIdx())
    template = Chem.RenumberAtoms(template, order)

    return template


def GetStructInfo(cif_path, cod_id):
    """
    Extracts the necessary information from the .cif file, including the atom and bond information. 
    Returns atom_dict -> {index: coords, label, symbol} and bond_dict -> {index: {index}}. 
    """
    file_path = cif_path + f"{cod_id[0]}/{cod_id[1:3]}/{cod_id[3:5]}/{cod_id}.cif"
    parser = CifParser(file_path)
    structure = parser.parse_structures(primitive=False)[0]
    structure *= (2, 2, 2)
    cif_data = parser.as_dict()[cod_id]

    # atom_dict -> {index: [position: label, symbol]}
    atom_dict = {}
    for i, site in enumerate(structure):
        dict = site.as_dict()
        atom_dict[i] = [site.coords, dict['label'], dict['species'][0]['element']]
    # print(atom_dict[440])

    # 获取成键原子信息
    bond_dict = {}
    try:
        label1 = cif_data['_geom_bond_atom_site_label_1']
        label2 = cif_data['_geom_bond_atom_site_label_2']
        labels = list(zip(label1, label2))
        print(labels)
        # 获取所有原子的相邻成键原子（但在 CIF 文件中通常不包括 H ）
        bond_label_dict = {}
        for label1, label2 in labels:
            if label1 not in bond_label_dict:
                bond_label_dict[label1] = []
            bond_label_dict[label1].append(label2)
            # if label2 not in bond_label_dict:
            #     bond_label_dict[label2] = []
            # bond_label_dict[label2].append(label1)
        print(bond_label_dict)
        # 由于成键两端标签和两端索引是一对多的关系，所以必须把标签字典转换成索引字典 (bond_dict)
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
                        distance = site_1.distance(site_2)
                        distances.append(distance)
                    # if min(distances) < 2:  # 添加距离判断，如果过大则表示两个原子不在同一分子上（位于晶胞的对侧），也就不可能成键
                    #     print(min(distances))
                    #     bond_index_2 = candidates_index_2[distances.index(min(distances))]  # 取出距离最近的索引作为成键原子
                    #     bond_dict[bond_index_1].add(bond_index_2)
                    # print(min(distances))
                    bond_index_2 = candidates_index_2[distances.index(min(distances))]  # 取出距离最近的索引作为成键原子
                    bond_dict[bond_index_1].add(bond_index_2)
        print(bond_dict)
    except KeyError:
        print(f"{cod_id}.cif has no bond information, searching for jmol version...")

    return atom_dict, bond_dict


def ConstructMol(atom_dict, bond_dict):
    """
    Extracts the molecule according to the atom_dict and bond_dict. 
    Return mol -> mol. 
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

    # 添加原子
    mol = Chem.EditableMol(mol)
    for symbol in symbols:
        mol.AddAtom(Chem.Atom(symbol))
    mol = mol.GetMol()

    # 添加原子坐标
    conf = Chem.Conformer(len(symbols))
    for i, pos in enumerate(positions):
        x, y, z = pos
        conf.SetAtomPosition(i, (x, y, z))
    mol.AddConformer(conf)

    # 添加键连信息
    mol = Chem.EditableMol(mol)
    for i, js in bond_dict.items():
        for j in js:
            # if i < j:   # 防止双向映射导致成键失败
            mol.AddBond(i, j, order=Chem.BondType.SINGLE)   # 由于缺少键类型等，默认按照单键处理
    mol = mol.GetMol()

    # 根据键长去除不合理的成键
    bonds_to_remove = []
    conf = mol.GetConformer()
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        pos_i = np.array(conf.GetAtomPosition(i))
        pos_j = np.array(conf.GetAtomPosition(j))
        bond_length = np.linalg.norm(pos_i - pos_j)
        # print(f"Bond {bond.GetIdx()}: {i}({a1.GetSymbol()}) - {j}({a2.GetSymbol()})")
        # print(f"Bond length: {pos_i} - {pos_j} = {bond_length}")
        if bond_length > 2:
            bonds_to_remove.append((i, j))
    if bonds_to_remove:
        rw_mol = Chem.RWMol(mol)
        for i, j in bonds_to_remove:
            rw_mol.RemoveBond(i, j)
        mol = rw_mol.GetMol()

    return mol


def OptimizeMol(mol):
    """
    Optimizes the mol structure, incluing removing the Hs, 
    allocating the aromaticity and Adding the Hs. 
    Return mol -> mol. 
    """
    rw_mol = Chem.RWMol(mol)
    bonds_to_remove = []
    for bond in rw_mol.GetBonds():
        atom1 = rw_mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        atom2 = rw_mol.GetAtomWithIdx(bond.GetEndAtomIdx())
        if atom1.GetAtomicNum() == 1 and atom2.GetAtomicNum() == 1:
            bonds_to_remove.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
    if bonds_to_remove:
        for begin_idx, end_idx in bonds_to_remove:
            rw_mol.RemoveBond(begin_idx, end_idx)
    h_atoms = [atom.GetIdx() for atom in rw_mol.GetAtoms() if atom.GetAtomicNum() == 1]
    for idx in sorted(h_atoms, reverse=True):
        rw_mol.RemoveAtom(idx)
    mol = rw_mol.GetMol()

    tolerance = 0.05
    sssr = Chem.GetSymmSSSR(mol)
    if sssr:
        n = len(sssr)
        i, j = 0, 0
        for ring in sssr:
            coords = []
            conf = mol.GetConformer()
            for atom_idx in ring:
                position = conf.GetAtomPosition(atom_idx)
                # print(f"position: ({position.x},{position.y},{position.z})")
                coords.append(position)
            coords = np.array(coords)
            centroid = np.mean(coords, axis=0)
            coords_centered = coords - centroid
            # SVD 拟合平面
            _, _, vh = np.linalg.svd(coords_centered)
            normal = vh[2, :]
            # 计算距离
            distances = np.abs(np.dot(coords_centered, normal))
            max_deviation = np.max(distances)
            rms_deviation = np.sqrt(np.mean(distances**2))
            if max_deviation < tolerance:
                # print(f"Ring {[atom_idx for atom_idx in ring]} is aromatic. ")
                for atom_idx in ring:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    atom.SetIsAromatic(True)
                for bond in mol.GetBonds():
                    begin_idx = bond.GetBeginAtomIdx()
                    end_idx = bond.GetEndAtomIdx()
                    if begin_idx in ring and end_idx in ring:
                        bond.SetIsAromatic(True)
                        bond.SetBondType(Chem.BondType.AROMATIC)
                i += 1
            else:
                # print(f"Ring {[atom_idx for atom_idx in ring]} is not aromatic. ")
                j += 1
        print(f"Found {n} rings, with {i} aromatic while {j} unaromatic. ")

    # for id in [1, 32, 65, 97, 128, 161, 192, 225, 256]:
    #     a = mol.GetAtomWithIdx(id)
    #     print(a.GetSymbol())

    sanitize_ops = (
        Chem.SanitizeFlags.SANITIZE_FINDRADICALS |
        Chem.SanitizeFlags.SANITIZE_SETAROMATICITY |
        Chem.SanitizeFlags.SANITIZE_SETCONJUGATION |
        Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION |
        Chem.SanitizeFlags.SANITIZE_CLEANUP 
        # Chem.SanitizeFlags.SANITIZE_CHECKVALENCE 
        )
    Chem.SanitizeMol(mol, sanitizeOps=sanitize_ops)
    mol = Chem.AddHs(mol, addCoords=True)

    # Chem.Kekulize(mol, clearAromaticFlags=False)
    # Chem.AssignAtomChiralTagsFromStructure(mol)
    # Chem.AssignStereochemistry(mol)
    # final_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    # print(final_smiles)
    # mol2 = Chem.MolFromSmiles(final_smiles)
    # mol2 = Chem.AddHs(mol2)

    # Chem.rdmolfiles.MolToMolFile(mol, "mol.mol")
    # Chem.rdmolfiles.MolToMolFile(mol2, "mol2.mol")

    return mol


def Search(cod_id, mol, substruct1, substruct2):
    """
    Tries to get the matches substructures, 
    Return the hits -> ((index), (index)), which exists in a molecule. 
    """
    hits = []
    match1s = mol.GetSubstructMatches(substruct1, useQueryQueryMatches=True)
    match2s = mol.GetSubstructMatches(substruct2, useQueryQueryMatches=True)

    if match1s and match2s:
        # 获取第一个子结构的所有重原子索引
        heavy_indices_1s = {}
        for match1 in match1s:
            heavy_indices_1 = set()
            for idx in match1:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 1:
                    heavy_indices_1.add(idx)
            if heavy_indices_1:
                heavy_indices_1s[frozenset(heavy_indices_1)] = match1
        # print(len(heavy_indices_1s))
        # 获取第二个子结构的所有重原子索引
        heavy_indices_2s = {}
        for match2 in match2s:
            heavy_indices_2 = set()
            for idx in match2:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 1:
                    heavy_indices_2.add(idx)
            if heavy_indices_2:
                heavy_indices_2s[frozenset(heavy_indices_2)] = match2
        # print(len(heavy_indices_2s))
        conf = mol.GetConformer()
        for heavy_indices_1, match1 in heavy_indices_1s.items():
            coords_1 = np.array([conf.GetAtomPosition(i) for i in heavy_indices_1])
            for heavy_indices_2, match2 in heavy_indices_2s.items():
                coords_2 = np.array([conf.GetAtomPosition(j) for j in heavy_indices_2])
                distances = cdist(coords_1, coords_2, metric='euclidean')
                if np.min(distances) < 5.0:
                    hits.append((match1, match2))
        # print(len(hits))
        if hits:
            print("Get hits! continue...")
        else:
            print("No matched distances, exiting...")
    elif match1s and not match2s:
        smiles = Chem.MolToSmiles(substruct2)
        print(f"{cod_id} has no substructure: {smiles}, exiting...")
    elif match2s and not match1s:
        smiles = Chem.MolToSmiles(substruct1)
        print(f"{cod_id} has no substructure: {smiles}, exiting...")
    else:
        print("No matched substructures, exiting...")


if __name__ == '__main__':
    # 待完善：暂时跳过 all_frags.py 的读取
    molecule1 = REF_FRAGMENTS['frag_43']
    molecule2 = REF_FRAGMENTS['frag_46']

    smarts_1 = molecule1.smarts.replace("=", "")    # 尝试？？？
    smarts_2 = molecule2.smarts.replace("=", "")
    smiles_1 = molecule1.smiles
    smiles_2 = molecule2.smiles
    heavy_atom_1 = [atom.GetIdx() for atom in Chem.MolFromSmiles(smiles_1).GetAtoms() if atom.GetSymbol() != "H"]
    heavy_atom_2 = [atom.GetIdx() for atom in Chem.MolFromSmiles(smiles_2).GetAtoms() if atom.GetSymbol() != "H"]
    num_heavy_1 = len(heavy_atom_1)
    num_heavy_2 = len(heavy_atom_2)
    heavy_indices_1 = get_heavy_indices_from_smarts(smarts_1)
    heavy_indices_2 = get_heavy_indices_from_smarts(smarts_2)

    print(f"Fragment 1: {smarts_1} -> {smiles_1}", flush=True)
    print(f"Fragment 2: {smarts_2} -> {smiles_2}", flush=True)
    print(f"Heavy atoms: {num_heavy_1} and {num_heavy_2}", flush=True)
    print("\nConstructing templates...", flush=True)
    template_1 = construct_template(smiles_1)
    template_2 = construct_template(smiles_2)
    num_1 = template_1.GetNumAtoms()
    num_2 = template_2.GetNumAtoms()
    print(f"Template atoms: {num_1} and {num_2}", flush=True)
    substruct1 = Chem.MolFromSmarts(smarts_1)
    substruct2 = Chem.MolFromSmarts(smarts_2)
    print("fragments: ", Chem.MolToSmarts(substruct1), Chem.MolToSmarts(substruct2))
    print("=" * 50)
    # 待完善：暂不考虑输入输出文件目录问题

    cif_path = '/home/clasper/2507/cif/'
    cod_id = '1000001'
    print(f"Loading file from {cod_id}.cif...")
    atom_dict, bond_dict = GetStructInfo(cif_path, cod_id)

    if atom_dict and bond_dict:
        print(f"Constructing mol for {cod_id}.cif")
        mol = ConstructMol(atom_dict, bond_dict)
        # print(Chem.MolToSmiles(mol))
    if not bond_dict:
        mol = Chem.MolFromMolFile(f"{cod_id}.mol")
        if not mol:
            print(f"Parse for {cod_id}.mol failed, please check if cif2mol.jmol has runned. ")

    if mol:
        print("Optimizing mol structure...")
        mol = OptimizeMol(mol)
        print("=" * 50)
        print(f"Searching cod_id {cod_id} for {smarts_1} and {smarts_2}...")
        hits = Search(cod_id, mol, substruct1, substruct2)
        if hits:
            print(f"Got {len(hits)} hits! continue...")
        else:
            print("No hits, exit.")
    else:
        print(f"Failed to get mol for {cod_id}.cif, exit. ")
