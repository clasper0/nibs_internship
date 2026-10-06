# coding=utf-8
import os
from pathlib import Path
from pymatgen.io.cif import CifParser   # 更完整地读取 cif 文件信息（包括键长等）
from rdkit import Chem
from rdkit.Chem import rdMolAlign, AllChem
from all_frags import Molecule
from all_frags import PROBE_FRAGMENTS, REF_FRAGMENTS
import numpy as np
from scipy.spatial.distance import cdist
import copy
import argparse


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

    # atom_dict -> {index: {position: np.array, label: str, symbol:str}}
    atom_dict = {}
    for i, site in enumerate(structure):
        dict = site.as_dict()
        species = dict['species'][0]
        atom_dict[i] = {'position': site.coords, 'label': dict['label'], 'symbol': species['element']}
        if 'oxidation_state in species' in species:
            atom_dict[i]['charge'] = species['oxidation_state']
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
            bond_indexes_1 = [index for index, info in atom_dict.items() if label1 == info['label']]
            for bond_index_1 in bond_indexes_1:
                bond_dict[bond_index_1] = set()
                
                # 建立嵌套子对象，存储候选 bond_indexes_2 的各个索引及坐标
                site_1 = structure[bond_index_1] # 成键原子 1 的 site, 用于计算距离，所用方法为: site_1.distance(site_2)
                for label2 in label2s:  # 获得候选成键原子的索引集合
                    candidates_index_2 = [index for index, info in atom_dict.items() if label2 == info['label']]
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
    positions, symbols, charges = [], [], []
    for value in atom_dict.values():
        positions.append(value['position'])
        symbols.append(value['symbol'])
        if 'charge' in value:
            charges.append(value['charge'])
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

    # 添加电荷信息
    for i, atom in enumerate(mol.GetAtoms()):
        charge = charges[i]
        # 先自定义属性，不知是否能利用上
        atom.SetProperty("PartialCharge", str(charge))

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
    i, j = 0, 0
    if bonds_to_remove:
        for begin_idx, end_idx in bonds_to_remove:
            rw_mol.RemoveBond(begin_idx, end_idx)
            i += 1
    h_atoms = [atom.GetIdx() for atom in rw_mol.GetAtoms() if atom.GetAtomicNum() == 1]
    for idx in sorted(h_atoms, reverse=True):
        rw_mol.RemoveAtom(idx)
        j += 1
    mol = rw_mol.GetMol()
    print(f"Removed {i} H-H and {j} Hs. ")

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
    else:
        print("No ring found. ")

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
    # 确定分子对象索引
    molecules = Chem.GetMolFrags(mol, asMols=False)
    mol_num_dict = {}
    for molecule in molecules:
        num = len(molecule)
        if num not in mol_num_dict:
            mol_num_dict[num] = 0
        mol_num_dict[num] += 1
    
    # 定义分子大小：晶胞为相同分子的重复，则重复最多次的应是完整分子（其余为晶胞边缘处的片段）
    valid_num = max(mol_num_dict.values())
    valid_sizes = [size for size, num in mol_num_dict.items() if num == valid_num]

    # 建立有效分子集合和原子 -> 分子映射
    valid_molecules = []
    atom_mol_dict = {}
    for mol_id, molecule in enumerate(molecules):
        if len(molecule) in valid_sizes:
            valid_molecules.append(mol_id)
        for atom_id in molecule:
            atom_mol_dict[atom_id] = mol_id
    
    hits = []
    match1s = mol.GetSubstructMatches(substruct1, useQueryQueryMatches=True)
    match2s = mol.GetSubstructMatches(substruct2, useQueryQueryMatches=True)

    if match1s and match2s:
        # 获取第一个子结构的所有重原子索引
        heavy_indices_1s = {}
        for match1 in match1s:
            if not all([atom_mol_dict[i] == atom_mol_dict[i + 1] for i in range(len(match1) - 1)]):
                print(f"match1 {match1} on different molecule! continue...")
                continue
            if not atom_mol_dict[match1[0]] in valid_molecules:
                print("match1 not on molecules, continue...")
                continue
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
            if not all([atom_mol_dict[i] == atom_mol_dict[i + 1] for i in range(len(match2) - 1)]):
                print(f"match2 {match2} on different molecule! ")
            if not atom_mol_dict[match2[0]] in valid_molecules:
                print("match2 not on molecules, continue...")
                continue
            heavy_indices_2 = set()
            for idx in match2:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 1:
                    heavy_indices_2.add(idx)
            if heavy_indices_2:
                heavy_indices_2s[frozenset(heavy_indices_2)] = match2
        # print(len(heavy_indices_2s))
        conf = mol.GetConformer()
        i = 0
        for heavy_indices_1, match1 in heavy_indices_1s.items():
            coords_1 = np.array([conf.GetAtomPosition(i) for i in heavy_indices_1])
            for heavy_indices_2, match2 in heavy_indices_2s.items():
                if atom_mol_dict[match1[0]] != atom_mol_dict[match2[0]]:
                    coords_2 = np.array([conf.GetAtomPosition(j) for j in heavy_indices_2])
                    distances = cdist(coords_1, coords_2, metric='euclidean')
                    if np.min(distances) < 5.0:
                        hits.append((match1, match2))
                else:
                    i += 1
                    continue
        # print(f"{cod_id} got {len(hits)} hits! Skipped {i} groups as they are of the same molecule. ")
        # print(len(hits))
    elif match1s and not match2s:
        smiles = Chem.MolToSmiles(substruct2)
        print(f"{cod_id} has no substructure: {smiles}, exiting...")
    elif match2s and not match1s:
        smiles = Chem.MolToSmiles(substruct1)
        print(f"{cod_id} has no substructure: {smiles}, exiting...")
    else:
        print("No matched substructures, exiting...")
    
    # 为提高效率，暂定一个晶胞只返回一个匹配结果
    return hits[0]


def construct_frag(mol, indices):
    """
    Constrcts the fragment mol from the total mol according to the indices. 
    Returns the fragment mol -> Mol. 
    """
    idx_map = {}

    # 添加原子
    frag = Chem.EditableMol(Chem.Mol())
    for new_idx, old_idx in enumerate(indices):
        atom = mol.GetAtomWithIdx(old_idx)
        frag.AddAtom(atom)
        idx_map[old_idx] = new_idx
    frag = frag.GetMol()
    
    # 添加坐标
    old_conf = mol.GetConformer()
    new_conf = Chem.Conformer(len(indices))
    for new_idx, old_idx in enumerate(indices):
        position = old_conf.GetAtomPosition(old_idx)
        new_conf.SetAtomPosition(new_idx, position)
    frag.AddConformer(new_conf)

    # 添加键
    frag = Chem.EditableMol(frag)
    for old_idx in indices:
        atom = mol.GetAtomWithIdx(old_idx)
        for neighbor in atom.GetNeighbors():
            neigh_id = neighbor.GetIdx()
            if neigh_id in indices and old_idx < neigh_id:
                bond = mol.GetBondBetweenAtoms(old_idx, neigh_id)
                frag.AddBond(idx_map[old_idx], idx_map[neigh_id], bond.GetBondType())
    
    return frag.GetMol()


def construct_frag_pair(mol, indices1, indices2):
    """
    Constructs the fragment pair within a cell. 
    Return frag_pair -> (Mol, Mol)
    """
    molFrag1 = construct_frag(mol, indices1)
    molFrag2 = construct_frag(mol, indices2)
    return (molFrag1, molFrag2)


def correct_H(mol, template, heavy_indices):
    """
    modify the hydrogen atom coordinates in the molecule based on a template
    """
    copy_mol = copy.deepcopy(mol)
    for i in range(copy_mol.GetNumAtoms()):
        if i not in heavy_indices:
            copy_mol.GetAtomWithIdx(i).SetAtomicNum(1)

    # 检查并移除H-H键
    rw_mol = Chem.RWMol(copy_mol)
    bonds_to_remove = []
    for bond in rw_mol.GetBonds():
        atom1 = rw_mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        atom2 = rw_mol.GetAtomWithIdx(bond.GetEndAtomIdx())
        if atom1.GetAtomicNum() == 1 and atom2.GetAtomicNum() == 1:
            bonds_to_remove.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
    
    if bonds_to_remove:
        for begin_idx, end_idx in bonds_to_remove:
            rw_mol.RemoveBond(begin_idx, end_idx)
    
    copy_mol = rw_mol.GetMol()

    simplified_mol = SimplifyMol(copy.deepcopy(copy_mol))
    simplified_template = SimplifyMol(copy.deepcopy(template))

    match_forward = simplified_template.GetSubstructMatch(simplified_mol)
    match_reverse = simplified_mol.GetSubstructMatch(simplified_template)

    if match_forward and match_reverse:
        template_to_mol_map = {template_idx: mol_idx for mol_idx, template_idx in enumerate(match_reverse)}

        rw_copy_mol = Chem.RWMol(copy_mol)
        for template_idx, mol_idx in template_to_mol_map.items():
            template_atom = template.GetAtomWithIdx(template_idx)
            mol_atom = rw_copy_mol.GetAtomWithIdx(mol_idx)
            
            if template_atom.GetAtomicNum() != 1:  
                # 设置电荷
                mol_atom.SetFormalCharge(template_atom.GetFormalCharge())
                
                # 计算该重原子在template中连接的氢原子数
                h_count = 0
                for neighbor in template_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 1:  
                        h_count += 1
                
                # 设置显式氢原子数
                mol_atom.SetNumExplicitHs(h_count)
                # print(f"set atom {mol_idx} of mol(template {template_idx}): "
                #       f"charge={template_atom.GetFormalCharge()}, H count={h_count}")

        # 设置键类型
        for template_bond in template.GetBonds():
            template_begin_idx = template_bond.GetBeginAtomIdx()
            template_end_idx = template_bond.GetEndAtomIdx()
            
            mol_begin_idx = template_to_mol_map.get(template_begin_idx)
            mol_end_idx = template_to_mol_map.get(template_end_idx)
            
            if mol_begin_idx is not None and mol_end_idx is not None:
                mol_bond = rw_copy_mol.GetBondBetweenAtoms(mol_begin_idx, mol_end_idx)
                if mol_bond is not None:
                    template_bond_type = template_bond.GetBondType()
                    mol_bond.SetBondType(template_bond_type)
                    # print(f"set bond {mol_begin_idx}-{mol_end_idx}: {template_bond_type}")

        copy_mol = rw_copy_mol.GetMol()

        H_cor_dict = {}
        for t_bond in template.GetBonds():
            t_atom1 = t_bond.GetBeginAtom()
            t_atom2 = t_bond.GetEndAtom()

            is_polar_h = False
            if t_atom1.GetSymbol() == 'H' and t_atom2.GetSymbol() != 'H':
                t_h_idx, t_heavy_idx = t_atom1.GetIdx(), t_atom2.GetIdx()
                is_polar_h = True
            elif t_atom2.GetSymbol() == 'H' and t_atom1.GetSymbol() != 'H':
                t_h_idx, t_heavy_idx = t_atom2.GetIdx(), t_atom1.GetIdx()
                is_polar_h = True

            if is_polar_h:
                target_length = get_bond_length(template, t_heavy_idx, t_h_idx)
                
                mol_heavy_idx = template_to_mol_map.get(t_heavy_idx)
                mol_h_idx = template_to_mol_map.get(t_h_idx)

                if mol_heavy_idx is None or mol_h_idx is None:
                    continue

                heavy_pos = np.array(copy_mol.GetConformer().GetAtomPosition(mol_heavy_idx))
                heavy_coord_key = tuple(heavy_pos.tolist())

                h_pos = np.array(copy_mol.GetConformer().GetAtomPosition(mol_h_idx))
                bond_vector = h_pos - heavy_pos
                current_length = np.linalg.norm(bond_vector)
                
                if current_length > 1e-6:  
                    unit_vector = bond_vector / current_length
                    new_h_pos = heavy_pos + unit_vector * target_length
                    
                    if heavy_coord_key not in H_cor_dict:
                        H_cor_dict[heavy_coord_key] = []
                    H_cor_dict[heavy_coord_key].append(new_h_pos)

        mol_rmHs = Chem.RemoveHs(copy.deepcopy(copy_mol))
        mol_addHs = Chem.AddHs(mol_rmHs, addCoords=True)
        # for atom in mol_addHs.GetAtoms():
        #     print(f"Atom: {atom.GetIdx()}, Symbol: {atom.GetSymbol()}, Position: {mol_addHs.GetConformer().GetAtomPosition(atom.GetIdx())}", flush=True)

        for heavy_coord, correct_H_coords in H_cor_dict.items():
            heavy_coord = np.array(heavy_coord)
            
            # 根据坐标找到重原子的索引
            heavy_idx = None
            for i in range(mol_addHs.GetNumAtoms()):
                if np.allclose(mol_addHs.GetConformer().GetAtomPosition(i), heavy_coord, atol=1e-3):
                    heavy_idx = i
                    break
            
            if heavy_idx is None:
                continue
                
            # 找到该重原子的所有氢邻居
            H_indices = [atom.GetIdx() for atom in mol_addHs.GetAtomWithIdx(heavy_idx).GetNeighbors() 
                        if atom.GetSymbol() == "H"]
            
            # 将校正后的坐标应用到氢原子上
            for i, h_idx in enumerate(H_indices):
                mol_addHs.GetConformer().SetAtomPosition(h_idx, correct_H_coords[i])

        return (True, mol_addHs)

    else:
        print("Debug: mol and template do not match, skipping hydrogen correction", flush=True)
        return (False, copy_mol)


def get_bond_length(mol, atom1_idx, atom2_idx):
    """
    calculate the bond length between two atoms in a molecule
    """
    conf = mol.GetConformer()
    pos1 = np.array(conf.GetAtomPosition(atom1_idx))
    pos2 = np.array(conf.GetAtomPosition(atom2_idx))
    return np.linalg.norm(pos2 - pos1)


def align_sidemol(ref_central_mol, qry_mol, mol_idx):
    qry_mols = Chem.GetMolFrags(qry_mol, asMols=True)
    # align central molecule of qry to ref
    atomMap = [(i,i) for i in range(num_heavy_1)]
    rmsd, transform = Chem.rdMolAlign.GetAlignmentTransform(qry_mols[0], ref_central_mol, atomMap=atomMap)
    rotation = transform[:3, :3]
    translation = transform[:3, 3]
    if rmsd > 0.5:
        print(f"{mol_idx+1} align error: {rmsd}", flush=True)
    # for qry_mol in qry_mols:
    #     # apply the transformation to the molecule coordinates
    conf = qry_mol.GetConformer()
    for atom_idx in range(qry_mol.GetNumAtoms()):
        pos = np.array(conf.GetAtomPosition(atom_idx))
        new_pos = np.dot(rotation, pos) + translation
        conf.SetAtomPosition(atom_idx, new_pos)
    return qry_mol


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--frag1", type=str, required=False, help="fragment 1", default="UREA"
    )
    parser.add_argument(
        "--frag2", type=str, required=False, help="fragment 2", default="MSM"
    )
    args = parser.parse_args()
    frag1 = args.frag1
    frag2 = args.frag2

    print(f"Starting workflow for fragment pair: {frag1} - {frag2}", flush=True)
    print("=" * 50, flush=True)
    
    molecule1 = REF_FRAGMENTS[frag1]
    molecule2 = PROBE_FRAGMENTS[frag2]
    # molecule1 = REF_FRAGMENTS['frag_14']
    # molecule2 = REF_FRAGMENTS['frag_46']

    smarts_1 = molecule1.smarts.replace("=", "")    # 目前需要去掉 all_frags.py 中的键类型和电荷信息
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
    
    
    # 设置 cif 路径，sdf 路径和输出文件路径
    home_dir = Path().home()
    work_path = home_dir / "2507"
    cif_dir = work_path / "cif"
    sdf_dir = work_path / "sdf"
    output_path = work_path / "output"
    ref_parent_dir = Path(f"{output_path}/{frag1}")
    pair_name = f"{frag1}_{frag2}"
    pair_name_dir = ref_parent_dir / pair_name
    pair_name_dir.mkdir(parents=True, exist_ok=True)
    corrected_file = f"{pair_name_dir}/{pair_name}.sdf"
    not_corrected_file = f"{pair_name_dir}/{pair_name}_uncorrect.sdf"


    # 待完善：获得 cod_ids 的思路：在外部写脚本遍历数据库下的文件名提取
    cod_ids = ['1100743']
    mol_pairs = []

    for cod_id in cod_ids:
        print(f"Loading file from {cod_id}.cif...")
        atom_dict, bond_dict = GetStructInfo(cif_dir, cod_id)
        bond_dict = []  # 使得程序自动搜索 jmol 转换的 sdf

        if atom_dict and bond_dict:
            print(f"Constructing mol for {cod_id}.cif")
            mol = ConstructMol(atom_dict, bond_dict)
            # print(Chem.MolToSmiles(mol))
        elif not bond_dict:
            try:
                sdf_file = sdf_dir / f"{cod_id}.sdf"
                mol = Chem.SDMolSupplier(sdf_file)[0]
            except OSError:
                print(f"Failed to parse {cod_id}.sdf, please check if cif2mol.jmol has runned. ")
                mol = None
        else:
            print(f"Failed to get mol for {cod_id}, please check. ")
            continue

        if mol:
            print("Optimizing mol structure...")
            mol = OptimizeMol(mol)
            print("=" * 50)
            print(f"Searching cod_id {cod_id} for {smarts_1} and {smarts_2}...")
            hit = Search(cod_id, mol, substruct1, substruct2)
            if hit:
                print(f"{cod_id} got hits! Getting mol pair...")
                match1, match2 = hit
                mol_pair = construct_frag_pair(mol, match1, match2)
                mol_pairs.append(mol_pair)
            else:
                print("No hits, continue...")
                continue
        else:
            print(f"Failed to get mol for {cod_id}.cif, continue...")
            continue


    print("\nCorrecting Hs...", flush=True)
    correctH_mol_pairs = []
    correctH_mol_names = []
    not_corrected_pairs = []
    not_corrected_names = []
    not_corrected_pairs_Hreplace = []
    for i, (frag1_mol, frag2_mol) in enumerate(mol_pairs):
        if i % 50 == 0:
            print(f"Processing molecule pair {i+1}/{len(mol_pairs)}", flush=True)
        if_correct_1, corrected_frag1 = correct_H(frag1_mol, template_1, heavy_indices_1)
        if_correct_2, corrected_frag2 = correct_H(frag2_mol, template_2, heavy_indices_2)
        if not if_correct_1 or not if_correct_2:
            print("Failed to correct. ")
            combined_mol_not_corrected = Chem.CombineMols(frag1_mol, frag2_mol)
            not_corrected_pairs.append(combined_mol_not_corrected)
            not_corrected_names.append(cod_id)  
            combined_mol_not_corrected_Hreplace = Chem.CombineMols(corrected_frag1, corrected_frag2)
            not_corrected_pairs_Hreplace.append(combined_mol_not_corrected_Hreplace)
            continue
        print("Corrected. ")
        combined_mol = Chem.CombineMols(corrected_frag1, corrected_frag2)
        correctH_mol_pairs.append(combined_mol)
        correctH_mol_names.append(cod_id)


    print("\nAligning molecules...", flush=True)
    ref_mol = correctH_mol_pairs[0]
    ref_mols = Chem.GetMolFrags(ref_mol, asMols=True)
    ref_central_mol, ref_side_mol = ref_mols[0], ref_mols[1]
    aligned_mols = []
    aligned_names = []

    for i, qry_mol in enumerate(correctH_mol_pairs):
        if i % 50 == 0:
            print(f"Aligning molecule {i+1}/{len(correctH_mol_pairs)}", flush=True)
        if i == 0:
            aligned_mols.append(qry_mol)
            aligned_names.append(correctH_mol_names[i])
            continue
        qry_mol = align_sidemol(ref_central_mol, qry_mol, i)
        aligned_mols.append(qry_mol)
        aligned_names.append(correctH_mol_names[i])

    print("Alignment completed.", flush=True)

    # delete similar molecules
    atomMap = [[(i,i) for i in range(num_heavy_2)]]
    print("\nRemoving redundant molecules...", flush=True)
    nonredun_mols = []
    nonredun_names = []

    for i, qry_mol in enumerate(aligned_mols):
        if i % 50 == 0:
            print(f"Checking redundancy for molecule {i+1}/{len(aligned_mols)}", flush=True)
        redun_flag = False
        for j, ref_mol in enumerate(aligned_mols):
            if j <= i:
                continue
            qry_side_mol = Chem.GetMolFrags(qry_mol, asMols=True)[1]
            ref_side_mol = Chem.GetMolFrags(ref_mol, asMols=True)[1]
            rmsd = Chem.rdMolAlign.CalcRMS(ref_side_mol, qry_side_mol, map=atomMap)
            if rmsd < 0.00001:
                print(f"Too similar, delete {i+1} and {j+1}, RMSD={rmsd:.6f}", flush=True)
                redun_flag = True
                break
        if not redun_flag:
            nonredun_mols.append(qry_mol)
            nonredun_names.append(aligned_names[i])

    print(f"Redundancy removal completed. Saving {len(nonredun_mols)} molecules", flush=True)

    # save the aligned molecule with names
    with Chem.SDWriter(str(corrected_file)) as writer:
        for mol, name in zip(nonredun_mols, nonredun_names):
            mol.SetProp("_Name", name)  
            writer.write(mol)
    print(f"Total molecules saved: {len(nonredun_mols)}", flush=True)

    # save the not corrected pairs with names
    with Chem.SDWriter(str(not_corrected_file)) as writer:
        for mol, name in zip(not_corrected_pairs, not_corrected_names):
            try:
                Chem.SanitizeMol(mol)
                mol.SetProp("_Name", name)  
                writer.write(mol)
            except Exception as e:
                print(f"Error sanitizing molecule {name}: {e}", flush=True)
                continue
    print(f"Not corrected pairs saved to: {not_corrected_file}", flush=True)

    print(f"\nWorkflow completed! Final output: {corrected_file}, Not corrected pairs: {not_corrected_file}", flush=True)
    print("=" * 50, flush=True)

