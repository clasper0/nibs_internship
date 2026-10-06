from rdkit import Chem
import numpy as np
from typing import Tuple, Optional, Union
import copy
import os
from pathlib import Path
import sys


# Combine all the unicon-generated sdfs
def CombineSdfs():
    sdf_dir = "/home/clasper/2507/tmp"
    sdf_files = [file for file in os.listdir(sdf_dir) if file.endswith(".sdf")]
    final_mol = ""
    i = 0
    for sdf_file in sdf_files:
        try:
            mol = Chem.SDMolSupplier(f"{sdf_dir}/{sdf_file}", sanitize=False, removeHs=False, strictParsing=False)[0]
            if not final_mol:
                final_mol = mol
            final_mol = Chem.CombineMols(final_mol, mol)
            i += 1
        except OSError:
            pass
    if i:
        return final_mol
    else:
        print("No valid sdf files. continue. ")
        return None
    

# CIF文件晶胞参数解析函数
def parse_cif_cell_params(cif_file_path):
    """
    从CIF文件中解析晶胞长度和角度参数
    
    参数:
        cif_file_path: CIF文件路径
    
    返回:
        tuple: (cell_params, success)
            - cell_params: [a, b, c, alpha, beta, gamma] 或 None
            - success: 是否成功解析
    """
    
    # 定义需要查找的CIF标签
    cell_tags = {
        '_cell_length_a': 0,
        '_cell_length_b': 1, 
        '_cell_length_c': 2,
        '_cell_angle_alpha': 3,
        '_cell_angle_beta': 4,
        '_cell_angle_gamma': 5
    }
    
    try:
        with open(cif_file_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except Exception as e:
        print(f"错误：无法读取CIF文件 {cif_file_path}: {e}")
        return None, False
    
    # 初始化参数数组
    cell_params = [None] * 6
    found_params = 0
    
    # 按行处理
    lines = content.split('\n')
    
    for line in lines:
        line = line.strip()
        
        # 跳过空行和注释行
        if not line or line.startswith('#'):
            continue
            
        # 查找晶胞参数
        for tag, index in cell_tags.items():
            if line.lower().startswith(tag.lower()):
                # 分割行获取值
                parts = line.split()
                if len(parts) >= 2:
                    value_str = parts[1]
                    
                    # 解析数值（处理误差标记，如 5.123(4)）
                    try:
                        # 移除括号中的误差值
                        if '(' in value_str:
                            value_str = value_str.split('(')[0]
                        
                        value = float(value_str)
                        cell_params[index] = value
                        found_params += 1
                        
                    except ValueError as e:
                        print(f"警告：无法解析参数 {tag}: {value_str}")
                        continue
    
    # 检查是否找到了所有6个参数
    if found_params == 6:
        print(f"成功解析CIF文件晶胞参数:")
        print(f"  a = {cell_params[0]:.6f} Å")
        print(f"  b = {cell_params[1]:.6f} Å") 
        print(f"  c = {cell_params[2]:.6f} Å")
        print(f"  α = {cell_params[3]:.6f}°")
        print(f"  β = {cell_params[4]:.6f}°")
        print(f"  γ = {cell_params[5]:.6f}°")
        return cell_params, True
    else:
        print(f"错误：只找到 {found_params}/6 个晶胞参数")
        print("缺失的参数:")
        for tag, index in cell_tags.items():
            if cell_params[index] is None:
                print(f"  {tag}")
        return None, False


def parse_cif_cell_params_advanced(cif_file_path):
    """
    高级CIF文件晶胞参数解析（支持多个数据块）
    
    参数:
        cif_file_path: CIF文件路径
    
    返回:
        dict: {block_name: [a, b, c, alpha, beta, gamma]} 或空字典
    """
    
    cell_tags = [
        '_cell_length_a', '_cell_length_b', '_cell_length_c',
        '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma'
    ]
    
    try:
        with open(cif_file_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except Exception as e:
        print(f"错误：无法读取CIF文件 {cif_file_path}: {e}")
        return {}
    
    # 分析数据块
    blocks = {}
    current_block = None
    current_data = {}
    
    lines = content.split('\n')
    
    for line in lines:
        line = line.strip()
        
        # 跳过空行和注释
        if not line or line.startswith('#'):
            continue
        
        # 检查数据块开始
        if line.lower().startswith('data_'):
            # 保存之前的块
            if current_block is not None:
                if len([v for v in current_data.values() if v is not None]) == 6:
                    blocks[current_block] = list(current_data.values())
            
            # 开始新块
            current_block = line[5:].strip()  # 移除 'data_'
            current_data = {tag: None for tag in cell_tags}
            continue
        
        # 在当前块中查找晶胞参数
        if current_block is not None:
            for tag in cell_tags:
                if line.lower().startswith(tag.lower()):
                    parts = line.split()
                    if len(parts) >= 2:
                        value_str = parts[1]
                        try:
                            # 处理误差标记
                            if '(' in value_str:
                                value_str = value_str.split('(')[0]
                            current_data[tag] = float(value_str)
                        except ValueError:
                            continue
    
    # 处理最后一个块
    if current_block is not None:
        if len([v for v in current_data.values() if v is not None]) == 6:
            blocks[current_block] = list(current_data.values())
    
    # 输出结果
    if blocks:
        print(f"从CIF文件中找到 {len(blocks)} 个数据块的晶胞参数:")
        for block_name, params in blocks.items():
            print(f"\n数据块: {block_name}")
            print(f"  a = {params[0]:.6f} Å")
            print(f"  b = {params[1]:.6f} Å")
            print(f"  c = {params[2]:.6f} Å") 
            print(f"  α = {params[3]:.6f}°")
            print(f"  β = {params[4]:.6f}°")
            print(f"  γ = {params[5]:.6f}°")
    else:
        print("未在CIF文件中找到完整的晶胞参数")
    
    return blocks


# Generating supercell
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from typing import Tuple, Optional, Union
import copy


def make_mol_supercell(mol, P=None, cell_params=None, wrap=True, order="cell-major", tol=1e-5):
    """
    为RDKit Mol对象生成超胞结构
    
    参数:
        mol: RDKit Mol对象，必须包含3D坐标
        P: 3x3 numpy数组，超胞变换矩阵。如果为None，默认为单位矩阵
        cell_params: 晶胞参数 [a, b, c, alpha, beta, gamma]，如果为None则自动计算
        wrap: 是否将原子坐标包装到晶胞内
        order: 原子排序方式，"cell-major" 或 "atom-major"
        tol: 容差值
    
    返回:
        新的RDKit Mol对象，包含超胞结构
    """
    
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("输入的Mol对象为空或不包含原子")
    
    # 检查是否有3D坐标
    conf = mol.GetConformer()
    if conf is None:
        raise ValueError("Mol对象不包含构象信息")
    
    # 默认超胞变换矩阵
    if P is None:
        P = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=int)
    
    # 获取原子坐标
    positions = []
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        positions.append([pos.x, pos.y, pos.z])
    positions = np.array(positions)
    
    # 计算或使用给定的晶胞参数
    print("cell_params: ", cell_params)
    if cell_params is None:
        cell_params = _estimate_cell_params(positions)
    
    # 从晶胞参数构建晶胞矩阵
    cell_matrix = _cellpar_to_cell(cell_params)
    
    # 计算超胞矩阵
    supercell_matrix = P
    supercell_cell = supercell_matrix @ cell_matrix
    
    # 获取晶格点
    lattice_points_frac = _lattice_points_in_supercell(supercell_matrix)
    lattice_points = np.dot(lattice_points_frac, supercell_cell)
    N = len(lattice_points)
    
    print(f"生成超胞: {P[0,0]}x{P[1,1]}x{P[2,2]}")
    print(f"原分子原子数: {mol.GetNumAtoms()}")
    print(f"晶格点数: {N}")
    print(f"超胞总原子数: {mol.GetNumAtoms() * N}")
    
    # 创建新的分子对象
    supercell_mol = Chem.RWMol()
    
    # 存储原始原子到新原子的映射
    atom_mapping = {}  # {(copy_idx, orig_atom_idx): new_atom_idx}
    
    # 为每个晶格点添加分子副本
    for copy_idx, lattice_point in enumerate(lattice_points):
        for orig_atom_idx in range(mol.GetNumAtoms()):
            orig_atom = mol.GetAtomWithIdx(orig_atom_idx)
            
            # 创建新原子，复制所有属性
            new_atom = Chem.Atom(orig_atom.GetAtomicNum())
            
            # 复制原子属性
            if orig_atom.HasProp('_Name'):
                new_atom.SetProp('_Name', orig_atom.GetProp('_Name'))
            if orig_atom.GetFormalCharge() != 0:
                new_atom.SetFormalCharge(orig_atom.GetFormalCharge())
            if orig_atom.GetNumRadicalElectrons() != 0:
                new_atom.SetNumRadicalElectrons(orig_atom.GetNumRadicalElectrons())
            if orig_atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                new_atom.SetChiralTag(orig_atom.GetChiralTag())
            if orig_atom.GetHybridization() != Chem.HybridizationType.UNSPECIFIED:
                new_atom.SetHybridization(orig_atom.GetHybridization())
            
            # 复制其他属性
            for prop_name in orig_atom.GetPropNames():
                if not prop_name.startswith('_'):  # 跳过内部属性
                    try:
                        prop_val = orig_atom.GetProp(prop_name)
                        new_atom.SetProp(prop_name, prop_val)
                    except:
                        pass
            
            # 添加原子到分子
            new_atom_idx = supercell_mol.AddAtom(new_atom)
            atom_mapping[(copy_idx, orig_atom_idx)] = new_atom_idx
    
    # 复制键连接
    for copy_idx in range(N):
        for bond in mol.GetBonds():
            begin_atom = bond.GetBeginAtomIdx()
            end_atom = bond.GetEndAtomIdx()
            
            # 获取新分子中对应的原子索引
            new_begin_idx = atom_mapping[(copy_idx, begin_atom)]
            new_end_idx = atom_mapping[(copy_idx, end_atom)]
            
            # 复制键属性
            bond_type = bond.GetBondType()
            bond_dir = bond.GetBondDir()
            is_aromatic = bond.GetIsAromatic()
            
            # 添加键
            new_bond_idx = supercell_mol.AddBond(new_begin_idx, new_end_idx, bond_type) - 1
            new_bond = supercell_mol.GetBondWithIdx(new_bond_idx)

            
            if bond_dir != Chem.BondDir.NONE:
                new_bond.SetBondDir(bond_dir)
            if is_aromatic:
                new_bond.SetIsAromatic(True)
            
            # 复制其他键属性
            for prop_name in bond.GetPropNames():
                if not prop_name.startswith('_'):
                    try:
                        prop_val = bond.GetProp(prop_name)
                        new_bond.SetProp(prop_name, prop_val)
                    except:
                        pass
    
    # 设置3D坐标
    mol_with_coords = supercell_mol.GetMol()
    new_conf = Chem.Conformer(mol_with_coords.GetNumAtoms())
    
    # 计算新的原子坐标
    for copy_idx, lattice_point in enumerate(lattice_points):
        for orig_atom_idx in range(mol.GetNumAtoms()):
            new_atom_idx = atom_mapping[(copy_idx, orig_atom_idx)]
            
            # 获取原始坐标
            orig_pos = conf.GetAtomPosition(orig_atom_idx)
            orig_coords = np.array([orig_pos.x, orig_pos.y, orig_pos.z])
            
            # 计算新坐标
            if order == "cell-major":
                new_coords = orig_coords + lattice_point
            elif order == "atom-major":
                new_coords = orig_coords + lattice_point
            else:
                raise ValueError(f"无效的排序方式: {order}")
            
            # 设置新坐标
            new_conf.SetAtomPosition(new_atom_idx, 
                                   (float(new_coords[0]), 
                                    float(new_coords[1]), 
                                    float(new_coords[2])))
    
    mol_with_coords.AddConformer(new_conf)
    

    
    # 坐标包装（如果需要）
    if wrap:
        _wrap_coordinates(mol_with_coords, supercell_cell, tol)
    
    return mol_with_coords


def _estimate_cell_params(positions):
    """从原子坐标估算晶胞参数"""
    min_coords = np.min(positions, axis=0)
    max_coords = np.max(positions, axis=0)
    
    # 添加一些padding
    padding = 0  # Angstrom
    lengths = max_coords - min_coords + 2 * padding
    
    # 默认为正交晶系
    a, b, c = lengths[0], lengths[1], lengths[2]
    alpha, beta, gamma = 90.0, 90.0, 90.0
    
    return [a, b, c, alpha, beta, gamma]


def _cellpar_to_cell(cellpar):
    """将晶胞参数转换为晶胞矩阵"""
    if len(cellpar) == 6:
        a, b, c, alpha, beta, gamma = cellpar
        alpha_r = np.radians(alpha)
        beta_r = np.radians(beta)
        gamma_r = np.radians(gamma)
        
        cos_alpha = np.cos(alpha_r)
        cos_beta = np.cos(beta_r)
        cos_gamma = np.cos(gamma_r)
        sin_gamma = np.sin(gamma_r)
        
        cell = np.zeros((3, 3))
        cell[0, 0] = a
        cell[1, 0] = b * cos_gamma
        cell[1, 1] = b * sin_gamma
        cell[2, 0] = c * cos_beta
        cell[2, 1] = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
        cell[2, 2] = c * np.sqrt(1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 
                               + 2 * cos_alpha * cos_beta * cos_gamma) / sin_gamma
        return cell
    else:
        return np.array(cellpar)


def _lattice_points_in_supercell(supercell_matrix):
    """找到超胞中包含的所有晶格点"""
    diagonals = np.array([
        [0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1],
        [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1],
    ])
    d_points = np.dot(diagonals, supercell_matrix)

    mins = np.min(d_points, axis=0)
    maxes = np.max(d_points, axis=0) + 1

    ar = np.arange(mins[0], maxes[0])[:, None] * np.array([1, 0, 0])[None, :]
    br = np.arange(mins[1], maxes[1])[:, None] * np.array([0, 1, 0])[None, :]
    cr = np.arange(mins[2], maxes[2])[:, None] * np.array([0, 0, 1])[None, :]

    all_points = ar[:, None, None] + br[None, :, None] + cr[None, None, :]
    all_points = all_points.reshape((-1, 3))

    frac_points = np.dot(all_points, np.linalg.inv(supercell_matrix))
    tvects = frac_points[np.all(frac_points < 1 - 1e-10, axis=1)
                         & np.all(frac_points >= -1e-10, axis=1)]
    return tvects


def _wrap_coordinates(mol, cell_matrix, tol):
    """将原子坐标包装到晶胞内"""
    conf = mol.GetConformer()
    inv_cell = np.linalg.inv(cell_matrix)
    
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        coords = np.array([pos.x, pos.y, pos.z])
        
        # 转换到分数坐标
        frac_coords = np.dot(coords, inv_cell.T)
        
        # 包装到[0, 1)范围内
        frac_coords = frac_coords % 1.0
        
        # 转换回笛卡尔坐标
        new_coords = np.dot(frac_coords, cell_matrix)
        
        # 更新坐标
        conf.SetAtomPosition(i, (float(new_coords[0]), 
                                float(new_coords[1]), 
                                float(new_coords[2])))


# 便捷函数：生成常见大小的超胞
def make_2x2x2_supercell(mol, **kwargs):
    """生成2x2x2超胞"""
    P = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]], dtype=int)
    return make_mol_supercell(mol, cell_params=cell_params, P=P, **kwargs)


def make_3x3x3_supercell(mol, **kwargs):
    """生成3x3x3超胞"""
    P = np.array([[3, 0, 0], [0, 3, 0], [0, 0, 3]], dtype=int)
    return make_mol_supercell(mol, cell_params=cell_params, P=P, **kwargs)


def make_custom_supercell(mol, nx=1, ny=1, nz=1, **kwargs):
    """生成自定义大小的超胞"""
    P = np.array([[nx, 0, 0], [0, ny, 0], [0, 0, nz]], dtype=int)
    return make_mol_supercell(mol, cell_params=cell_params, P=P, **kwargs)


if __name__ == '__main__':
    cif_dir = Path("/home/clasper/2507/cif")
    sdf_dir = Path("/home/clasper/2507/sdf")
    cif_name = sys.argv[1].split("/")[-1]
    cod_id = cif_name[0:7]
    cif_path = cif_dir / cif_name[0] / cif_name[1:3] / cif_name[3:5] / cif_name
    sdf_path = sdf_dir / cif_name[0] / cif_name[1:3] / cif_name[3:5] / f"{cod_id}.sdf"
    cell_params, success = parse_cif_cell_params(cif_path)

    if not success:
        cell_params = None

    mol = CombineSdfs()
    
    print("原分子信息:")
    print(f"原子数: {mol.GetNumAtoms()}")
    print(f"键数: {mol.GetNumBonds()}")
    
    # 根据原子数生成不过大的超胞
    if mol.GetNumAtoms() < 370:
        supercell = make_3x3x3_supercell(mol)
    elif mol.GetNumAtoms() < 1250:
        supercell = make_2x2x2_supercell(mol)
    else:
        supercell = make_custom_supercell(mol)

    # 去除不合理成键
    bonds_to_remove = []
    conf = supercell.GetConformer()
    for bond in supercell.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        pos_i = np.array(conf.GetAtomPosition(i))
        pos_j = np.array(conf.GetAtomPosition(j))
        bond_length = np.linalg.norm(pos_i - pos_j)
        if bond_length > 5:
            bonds_to_remove.append((i, j))
    if bonds_to_remove:
        rw_mol = Chem.RWMol(supercell)
        for i, j in bonds_to_remove:
            rw_mol.RemoveBond(i, j)
        supercell = rw_mol.GetMol()
    
    print(f"\n超胞信息:")
    print(f"原子数: {supercell.GetNumAtoms()}")
    print(f"键数: {supercell.GetNumBonds()}")
    print(f"预期原子数: {mol.GetNumAtoms() * 8}")  # 2^3 = 8
    print(f"预期键数: {mol.GetNumBonds() * 8}")
    
    writer = Chem.SDWriter(sdf_path)
    writer.write(supercell)
    writer.close()