#!/usr/bin/env python3
"""
CIF到超胞SDF转换器

该脚本从CIF文件开始，执行以下步骤：
1. 解析CIF文件并生成完整的单晶胞结构
2. 将晶胞内的分子分离为独立的PDB文件
3. 使用Unicon为每个分子分配键连信息
4. 保留原始原子坐标
5. 将处理后的分子合并
6. 生成3×3×3超胞
7. 输出为RDKit兼容的SDF格式

使用方法:
    python cif_to_supercell_sdf.py input.cif [options]

依赖:
    - numpy
    - unicon (外部程序)
    - rdkit (可选，用于验证)
"""

import sys
import os
import subprocess
import argparse
import tempfile
import shutil
import time
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Set, Union
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
import re
import shlex
import warnings
import collections.abc

# 导入第一个脚本的CIF解析功能
CIFDataValue = Union[str, int, float]
CIFData = Union[CIFDataValue, List[CIFDataValue]]

class Atoms:
    """简单的原子结构类"""
    
    def __init__(self, symbols=None, positions=None, cell=None, pbc=None):
        self.symbols = symbols or []
        self.positions = np.array(positions) if positions is not None else np.array([])
        self.cell = np.array(cell) if cell is not None else np.zeros((3, 3))
        self.pbc = pbc or [True, True, True]
        self.info = {}
        self.arrays = {}
        
    def __len__(self):
        return len(self.symbols)
    
    def get_chemical_formula(self):
        """获取化学式"""
        if not self.symbols:
            return ""
        
        from collections import Counter
        counts = Counter(self.symbols)
        formula = ""
        for symbol in sorted(counts.keys()):
            count = counts[symbol]
            if count == 1:
                formula += symbol
            else:
                formula += f"{symbol}{count}"
        return formula
    
    def get_scaled_positions(self, wrap=False):
        """获取分数坐标"""
        if self.cell.shape == (3, 3) and np.linalg.det(self.cell) != 0:
            inv_cell = np.linalg.inv(self.cell)
            scaled_pos = np.dot(self.positions, inv_cell.T)
            if wrap:
                scaled_pos = scaled_pos % 1.0
            return scaled_pos
        return self.positions.copy()
    
    def set_scaled_positions(self, scaled_positions):
        """从分数坐标设置位置"""
        if self.cell.shape == (3, 3) and np.linalg.det(self.cell) != 0:
            self.positions = np.dot(scaled_positions, self.cell)
        else:
            self.positions = np.array(scaled_positions)
    
    def set_array(self, name, array):
        """为原子设置数组"""
        self.arrays[name] = np.array(array)
    
    def wrap(self, eps=1e-7):
        """将原子包装到晶胞内"""
        scaled = self.get_scaled_positions()
        scaled = scaled % 1.0
        self.set_scaled_positions(scaled)


class Cell:
    """简单的晶胞类"""
    
    def __init__(self, cell_array):
        self.array = np.array(cell_array)
        self.rank = 3 if self.array.shape == (3, 3) and np.linalg.det(self.array) != 0 else 0
    
    def cellpar(self):
        """获取晶胞参数 [a, b, c, alpha, beta, gamma]"""
        if self.rank != 3:
            return [0, 0, 0, 0, 0, 0]
        
        cell = self.array
        a = np.linalg.norm(cell[0])
        b = np.linalg.norm(cell[1])
        c = np.linalg.norm(cell[2])
        
        alpha = np.arccos(np.dot(cell[1], cell[2]) / (b * c)) * 180 / np.pi
        beta = np.arccos(np.dot(cell[0], cell[2]) / (a * c)) * 180 / np.pi
        gamma = np.arccos(np.dot(cell[0], cell[1]) / (a * b)) * 180 / np.pi
        
        return [a, b, c, alpha, beta, gamma]
    
    @classmethod
    def new(cls, cellpar):
        """从参数创建晶胞"""
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
            return cls(cell)
        else:
            return cls(np.array(cellpar))


def handle_subscripts(text):
    """处理下标"""
    return text


def convert_value(value: str) -> CIFDataValue:
    """将CIF值字符串转换为对应的python类型"""
    value = value.strip()
    if re.match('(".*")|(\'.*\')', value):
        return handle_subscripts(value[1:-1])
    elif re.match(r'[+-]?\d+$', value):
        return int(value)
    elif re.match(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$', value):
        return float(value)
    elif re.match(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\(\d+\)$', value):
        return float(value[:value.index('(')])
    elif re.match(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\(\d+', value):
        warnings.warn(f'Badly formed number: "{value}"')
        return float(value[:value.index('(')])
    else:
        return handle_subscripts(value)

def parse_multiline_string(lines: List[str], line: str) -> str:
    """解析分号包围的多行字符串"""
    assert line[0] == ';'
    strings = [line[1:].lstrip()]
    while True:
        line = lines.pop().strip()
        if line[:1] == ';':
            break
        strings.append(line)
    return '\n'.join(strings).strip()


def parse_singletag(lines: List[str], line: str) -> Tuple[str, CIFDataValue]:
    """解析CIF标签，返回键值对"""
    kv = line.split(None, 1)
    if len(kv) == 1:
        key = line
        line = lines.pop().strip()
        while not line or line[0] == '#':
            line = lines.pop().strip()
        if line[0] == ';':
            value = parse_multiline_string(lines, line)
        else:
            value = line
    else:
        key, value = kv
    return key, convert_value(value)


def parse_cif_loop_headers(lines: List[str]):
    while lines:
        line = lines.pop()
        tokens = line.split()

        if len(tokens) == 1 and tokens[0].startswith('_'):
            header = tokens[0].lower()
            yield header
        else:
            lines.append(line)
            return


def parse_cif_loop_data(lines: List[str], ncolumns: int) -> List[List[CIFDataValue]]:
    columns: List[List[CIFDataValue]] = [[] for _ in range(ncolumns)]

    tokens = []
    while lines:
        line = lines.pop().strip()
        lowerline = line.lower()
        if (not line or
            line.startswith('_') or
            lowerline.startswith('data_') or
                lowerline.startswith('loop_')):
            lines.append(line)
            break

        if line.startswith('#'):
            continue

        line = line.split(' #')[0]

        if line.startswith(';'):
            moretokens = [parse_multiline_string(lines, line)]
        else:
            if ncolumns == 1:
                moretokens = [line]
            else:
                moretokens = shlex.split(line, posix=False)

        tokens += moretokens
        if len(tokens) < ncolumns:
            continue
        if len(tokens) == ncolumns:
            for i, token in enumerate(tokens):
                columns[i].append(convert_value(token))
        else:
            warnings.warn(f'Wrong number {len(tokens)} of tokens, '
                          f'expected {ncolumns}: {tokens}')

        tokens = []

    return columns


def parse_loop(lines: List[str]) -> Dict[str, List[CIFDataValue]]:
    """解析CIF循环，返回以列标签名为键的字典"""
    headers = list(parse_cif_loop_headers(lines))
    columns = parse_cif_loop_data(lines, len(headers))

    columns_dict = {}
    for i, header in enumerate(headers):
        if header in columns_dict:
            warnings.warn(f'Duplicated loop tags: {header}')
        else:
            columns_dict[header] = columns[i]
    return columns_dict


def parse_items(lines: List[str], line: str) -> Dict[str, CIFData]:
    """解析CIF数据项并返回包含所有标签的字典"""
    tags: Dict[str, CIFData] = {}

    while True:
        if not lines:
            break
        line = lines.pop().strip()
        if not line:
            continue
        lowerline = line.lower()
        if not line or line.startswith('#'):
            continue
        elif line.startswith('_'):
            key, value = parse_singletag(lines, line)
            tags[key.lower()] = value
        elif lowerline.startswith('loop_'):
            tags.update(parse_loop(lines))
        elif lowerline.startswith('data_'):
            if line:
                lines.append(line)
            break
        elif line.startswith(';'):
            parse_multiline_string(lines, line)
        else:
            raise ValueError(f'Unexpected CIF file entry: "{line}"')
    return tags


def parse_symmetry_operation(symop_str):
    """解析对称操作字符串，如'x,y,z'或'x+1/2,-y,z'"""
    parts = [p.strip() for p in symop_str.replace(' ', '').split(',')]
    
    if len(parts) != 3:
        raise ValueError(f"Invalid symmetry operation: {symop_str}")
    
    rotation = np.zeros((3, 3))
    translation = np.zeros(3)
    
    for i, part in enumerate(parts):
        for j, var in enumerate(['x', 'y', 'z']):
            if var in part:
                var_idx = part.find(var)
                if var_idx == 0 or part[var_idx-1] == '+':
                    coeff = 1.0
                elif part[var_idx-1] == '-':
                    coeff = -1.0
                else:
                    start = var_idx - 1
                    while start > 0 and part[start-1] not in '+-':
                        start -= 1
                    coeff_str = part[start:var_idx]
                    if coeff_str == '' or coeff_str == '+':
                        coeff = 1.0
                    elif coeff_str == '-':
                        coeff = -1.0
                    else:
                        coeff = float(coeff_str)
                
                rotation[i, j] = coeff
        
        trans_part = part
        for var in ['x', 'y', 'z']:
            trans_part = re.sub(r'[+-]?[\d]*' + var, '', trans_part)
        
        if trans_part:
            trans_value = 0.0
            terms = re.split(r'(?=[+-])', trans_part)
            for term in terms:
                if term.strip():
                    if '/' in term:
                        if '+' in term or '-' in term:
                            sign = 1 if term.startswith('+') else -1 if term.startswith('-') else 1
                            frac_part = term.lstrip('+-')
                        else:
                            sign = 1
                            frac_part = term
                        
                        if '/' in frac_part:
                            num, den = frac_part.split('/')
                            trans_value += sign * float(num) / float(den)
                        else:
                            trans_value += sign * float(frac_part)
                    else:
                        if term.strip():
                            trans_value += float(term)
            
            translation[i] = trans_value
    
    return rotation, translation


class CIFBlock:
    """晶体学信息文件中的一个块"""

    cell_tags = ['_cell_length_a', '_cell_length_b', '_cell_length_c',
                 '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma']

    def __init__(self, name: str, tags: Dict[str, CIFData]):
        self.name = name
        self._tags = tags

    def __getitem__(self, key: str) -> CIFData:
        return self._tags[key]

    def __contains__(self, key):
        return key in self._tags

    def get(self, key, default=None):
        return self._tags.get(key, default)

    def get_cellpar(self) -> Optional[List]:
        try:
            return [self[tag] for tag in self.cell_tags]
        except KeyError:
            return None

    def get_cell(self) -> Cell:
        cellpar = self.get_cellpar()
        if cellpar is None:
            return Cell(np.zeros((3, 3)))
        return Cell.new(cellpar)

    def _raw_scaled_positions(self) -> Optional[np.ndarray]:
        coords = [self.get(name) for name in ['_atom_site_fract_x',
                                              '_atom_site_fract_y',
                                              '_atom_site_fract_z']]
        if None in coords:
            return None
        return np.array(coords).T

    def get_symbols(self) -> List[str]:
        labels = self._get_any(['_atom_site_type_symbol', '_atom_site_label'])
        if labels is None:
            raise ValueError('No symbols found')

        symbols = []
        for label in labels:
            if label == '.' or label == '?':
                raise ValueError('Symbols are undetermined')
            match = re.search(r'([A-Z][a-z]?)', str(label))
            if match:
                symbol = match.group(0)
                if symbol == 'D':
                    symbol = 'H'
                symbols.append(symbol)
            else:
                symbols.append('X')
        return symbols

    def _get_any(self, names):
        for name in names:
            if name in self:
                return self[name]
        return None
    
    def get_symmetry_operations(self):
        """从CIF获取对称操作"""
        symop_tags = [
            '_symmetry_equiv_pos_as_xyz',
            '_space_group_symop_operation_xyz',
            '_symmetry_equiv_pos_xyz'
        ]
        
        symops = None
        for tag in symop_tags:
            symops = self.get(tag)
            if symops is not None:
                break
        
        if symops is None:
            print("Warning: No symmetry operations found, using identity operation")
            return [np.eye(3), np.zeros(3)]
        
        if isinstance(symops, str):
            symops = [symops]
        
        operations = []
        for symop in symops:
            try:
                rotation, translation = parse_symmetry_operation(symop)
                operations.append((rotation, translation))
            except Exception as e:
                print(f"Warning: Could not parse symmetry operation '{symop}': {e}")
                continue
        
        if not operations:
            operations = [(np.eye(3), np.zeros(3))]
        
        return operations

    def get_occupancies(self):
        """获取位点占有率"""
        occs = self.get('_atom_site_occupancy')
        if occs is None:
            symbols = self.get_symbols()
            return [1.0] * len(symbols)
        return occs

    def get_atoms(self) -> Atoms:
        """从CIF标签返回包含所有对称生成原子的Atoms对象"""
        cell = self.get_cell()
        
        base_symbols = self.get_symbols()
        base_positions = self._raw_scaled_positions()
        if base_positions is None:
            raise ValueError('No atomic positions found')
        
        occupancies = self.get_occupancies()
        symops = self.get_symmetry_operations()
        
        print(f"Found {len(base_symbols)} unique atomic sites")
        print(f"Found {len(symops)} symmetry operations")
        
        all_symbols = []
        all_positions = []
        tolerance = 1e-5
        
        for rotation, translation in symops:
            for i, (symbol, pos, occ) in enumerate(zip(base_symbols, base_positions, occupancies)):
                if abs(occ) < tolerance:
                    continue
                
                new_pos = np.dot(rotation, pos) + translation
                new_pos = new_pos % 1.0
                
                is_duplicate = False
                for existing_pos in all_positions:
                    diff = new_pos - existing_pos
                    diff = diff - np.round(diff)
                    if np.linalg.norm(diff) < tolerance:
                        is_duplicate = True
                        break
                
                if not is_duplicate:
                    all_symbols.append(symbol)
                    all_positions.append(new_pos)
        
        print(f"Generated {len(all_symbols)} atoms in the complete unit cell")
        
        atoms = Atoms(symbols=all_symbols, cell=cell.array)
        atoms.set_scaled_positions(np.array(all_positions))
        
        return atoms


def parse_block(lines: List[str], line: str) -> CIFBlock:
    assert line.lower().startswith('data_')
    blockname = line.split('_', 1)[1].rstrip()
    tags = parse_items(lines, line)
    return CIFBlock(blockname, tags)


def parse_cif(fileobj):
    """解析CIF文件"""
    if isinstance(fileobj, str):
        with open(fileobj, 'rb') as f:
            data = f.read()
    else:
        data = fileobj.read()

    if isinstance(data, bytes):
        data = data.decode('latin1')
    
    lines = [e for e in data.split('\n') if len(e) > 0]
    lines = [''] + lines[::-1]

    while lines:
        line = lines.pop().strip()
        if not line or line.startswith('#'):
            continue
        yield parse_block(lines, line)


def read_cif(filename):
    """读取CIF文件并返回Atoms对象"""
    blocks = list(parse_cif(filename))
    if not blocks:
        raise ValueError("No data blocks found in CIF file")
    
    for block in blocks:
        try:
            atoms = block.get_atoms()
            return atoms
        except (ValueError, KeyError) as e:
            print(f"Warning: Skipping block {block.name}: {e}")
            continue
    
    raise ValueError("No valid structure data found in CIF file")


class Atom:
    """原子类，存储原子信息"""
    
    def __init__(self, symbol, x, y, z, atom_id=None):
        self.symbol = symbol
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.atom_id = atom_id or 1
        
    def get_coordinates(self) -> np.ndarray:
        """获取原子坐标"""
        return np.array([self.x, self.y, self.z])
    
    def distance_to(self, other_atom: 'Atom') -> float:
        """计算到另一个原子的距离"""
        coord1 = self.get_coordinates()
        coord2 = other_atom.get_coordinates()
        return np.linalg.norm(coord1 - coord2)
    
    def to_pdb_line(self, atom_number: int = None) -> str:
        """转换为PDB格式行"""
        if atom_number is None:
            atom_number = self.atom_id
        
        atom_name = f"{self.symbol:<4}"
        res_name = "UNL"
        chain_id = "A"
        res_seq = 1
        
        return (f"HETATM{atom_number:5d} {atom_name} {res_name} {chain_id}{res_seq:4d}    "
               f"{self.x:8.3f}{self.y:8.3f}{self.z:8.3f}"
               f"{1.00:6.2f}{0.00:6.2f}          {self.symbol:>2}")


class Molecule:
    """分子类，包含一个分子的所有原子"""
    
    def __init__(self, atoms: List[Atom], molecule_id: int):
        self.atoms = atoms
        self.molecule_id = molecule_id
    
    def __len__(self) -> int:
        return len(self.atoms)
    
    def get_center(self) -> np.ndarray:
        """获取分子的几何中心"""
        coords = np.array([atom.get_coordinates() for atom in self.atoms])
        return np.mean(coords, axis=0)
    
    def get_formula(self) -> str:
        """获取分子式"""
        from collections import Counter
        elements = [atom.symbol for atom in self.atoms]
        element_count = Counter(elements)
        
        formula = ""
        order = ['C', 'H', 'N', 'O'] + sorted([e for e in element_count.keys() if e not in ['C', 'H', 'N', 'O']])
        
        for element in order:
            if element in element_count:
                count = element_count[element]
                if count == 1:
                    formula += element
                else:
                    formula += f"{element}{count}"
        
        return formula
    
    def to_pdb_string(self, title: str = None) -> str:
        """转换为完整的PDB字符串"""
        lines = []
        
        if title:
            lines.append(f"TITLE     {title}")
        else:
            lines.append(f"TITLE     Molecule {self.molecule_id} ({self.get_formula()})")
        
        lines.append("REMARK    Extracted from unit cell structure")
        lines.append(f"REMARK    Molecule ID: {self.molecule_id}")
        lines.append(f"REMARK    Number of atoms: {len(self.atoms)}")
        lines.append(f"REMARK    Molecular formula: {self.get_formula()}")
        
        for i, atom in enumerate(self.atoms, 1):
            lines.append(atom.to_pdb_line(i))
        
        lines.append("END")
        
        return "\n".join(lines)


class MoleculeSplitter:
    """分子分离器"""
    
    def __init__(self, bond_distance_cutoff: float = 2.0, hydrogen_cutoff: float = 1.3):
        self.bond_distance_cutoff = bond_distance_cutoff
        self.hydrogen_cutoff = hydrogen_cutoff
        self.atoms = []
        self.molecules = []
    
    def load_from_atoms_object(self, atoms: Atoms):
        """从Atoms对象加载原子"""
        print(f"Loading {len(atoms)} atoms from unit cell")
        
        self.atoms = []
        for i, (symbol, pos) in enumerate(zip(atoms.symbols, atoms.positions)):
            atom = Atom(symbol, pos[0], pos[1], pos[2], i+1)
            self.atoms.append(atom)
    
    def _are_bonded(self, atom1: Atom, atom2: Atom) -> bool:
        """判断两个原子是否成键"""
        distance = atom1.distance_to(atom2)
        
        if atom1.symbol == 'H' or atom2.symbol == 'H':
            return distance <= self.hydrogen_cutoff
        else:
            return distance <= self.bond_distance_cutoff
    
    def _find_connected_atoms(self, start_atom_idx: int, visited: Set[int]) -> List[int]:
        """使用深度优先搜索找到连接的原子"""
        if start_atom_idx in visited:
            return []
        
        connected = [start_atom_idx]
        visited.add(start_atom_idx)
        
        start_atom = self.atoms[start_atom_idx]
        
        for i, other_atom in enumerate(self.atoms):
            if i not in visited and self._are_bonded(start_atom, other_atom):
                connected.extend(self._find_connected_atoms(i, visited))
        
        return connected
    
    def split_molecules(self) -> List[Molecule]:
        """分离分子"""
        print("Splitting molecules using connectivity analysis...")
        
        visited = set()
        molecules = []
        molecule_id = 1
        
        for i, atom in enumerate(self.atoms):
            if i not in visited:
                connected_indices = self._find_connected_atoms(i, visited)
                
                if connected_indices:
                    molecule_atoms = [self.atoms[idx] for idx in connected_indices]
                    molecule = Molecule(molecule_atoms, molecule_id)
                    molecules.append(molecule)
                    molecule_id += 1
        
        self.molecules = molecules
        print(f"Found {len(molecules)} molecules")
        
        for mol in molecules:
            center = mol.get_center()
            print(f"  Molecule {mol.molecule_id}: {len(mol)} atoms, "
                  f"formula: {mol.get_formula()}, "
                  f"center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
        
        return molecules
    
    def save_molecules(self, output_dir: str, filename_pattern: str = "molecule_{:03d}.pdb", 
                      filter_hydrogen_only: bool = True) -> List[str]:
        """保存分子到单独的PDB文件"""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        saved_files = []
        skipped_count = 0
        
        print(f"\nSaving molecules to {output_dir}:")
        for molecule in self.molecules:
            # 检查是否只包含氢原子
            if filter_hydrogen_only:
                heavy_atoms = [atom for atom in molecule.atoms if atom.symbol != 'H']
                if len(heavy_atoms) == 0:
                    print(f"  Skipped: Molecule {molecule.molecule_id} (only hydrogen atoms)")
                    skipped_count += 1
                    continue
                elif len(heavy_atoms) == 1:
                    print(f"  Warning: Molecule {molecule.molecule_id} has only 1 heavy atom")
            
            filename = filename_pattern.format(molecule.molecule_id)
            filepath = output_path / filename
            
            with open(filepath, 'w') as f:
                f.write(molecule.to_pdb_string())
            
            saved_files.append(str(filepath))
            print(f"  Saved: {filename} ({len(molecule)} atoms, {molecule.get_formula()})")
        
        if skipped_count > 0:
            print(f"  Skipped {skipped_count} molecules with only hydrogen atoms")
        
        return saved_files


class UniconProcessor:
    """Unicon处理器"""
    
    def __init__(self, unicon_path: str = "unicon", timeout: int = 300):
        self.unicon_path = unicon_path
        self.timeout = timeout
    
    def check_unicon_availability(self) -> bool:
        """检查unicon是否可用"""
        try:
            # 尝试运行unicon，检查是否存在
            result = subprocess.run(
                [self.unicon_path], 
                capture_output=True, 
                text=True, 
                timeout=5
            )
            # 如果命令存在，即使返回错误码也说明程序可用
            return True
            
        except FileNotFoundError:
            return False
        except subprocess.TimeoutExpired:
            # 超时也说明程序存在
            return True
        except Exception:
            return False
    
    def process_molecule(self, input_file: str, output_file: str = None) -> Tuple[bool, str, str]:
        """使用unicon处理单个分子"""
        if output_file is None:
            input_path = Path(input_file)
            output_file = str(input_path.with_suffix('.sdf'))
        
        # 构建unicon命令 - 使用标准参数格式
        cmd = [
            self.unicon_path,
            "-i", input_file,
            "-o", output_file,
            "-t", "sdf"  # 明确指定输出格式为SDF
        ]
        
        try:
            print(f"Processing {Path(input_file).name} with unicon...")
            print(f"Command: {' '.join(cmd)}")
            
            # 确保输出目录存在
            os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=self.timeout,
                cwd=os.path.dirname(os.path.abspath(input_file))  # 在输入文件目录运行
            )
            
            print(f"Unicon return code: {result.returncode}")
            if result.stdout:
                print(f"Unicon stdout: {result.stdout[:500]}")
            if result.stderr:
                print(f"Unicon stderr: {result.stderr[:500]}")
            
            # 检查输出文件是否存在且不为空
            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                # 验证SDF文件内容
                if self._validate_sdf_file(output_file):
                    if result.returncode == 0:
                        return True, output_file, ""
                    elif result.returncode == 9:
                        warning_msg = "Unicon completed with warnings, but SDF file created"
                        return True, output_file, warning_msg
                    else:
                        # 其他返回码但文件有效
                        warning_msg = f"Unicon returned code {result.returncode}, but valid SDF created"
                        return True, output_file, warning_msg
                else:
                    # SDF文件无效
                    error_msg = f"Unicon generated invalid SDF file\n"
                    if result.stderr:
                        error_msg += f"STDERR: {result.stderr[:200]}\n"
                    return False, "", error_msg
            else:
                # 输出文件不存在或为空
                error_msg = f"Unicon failed - no output file generated (return code: {result.returncode})\n"
                
                # 分析错误原因
                if result.stderr:
                    stderr_lower = result.stderr.lower()
                    if "no heavy atoms" in stderr_lower:
                        error_msg += "Reason: Molecule contains only hydrogen atoms\n"
                    elif "invalid" in stderr_lower or "error" in stderr_lower:
                        error_msg += "Reason: Invalid input structure\n"
                    elif "not found" in stderr_lower:
                        error_msg += "Reason: Input file or dependency not found\n"
                    error_msg += f"STDERR: {result.stderr[:200]}\n"
                
                if result.stdout:
                    error_msg += f"STDOUT: {result.stdout[:200]}\n"
                    
                return False, "", error_msg
                
        except subprocess.TimeoutExpired:
            return False, "", f"Unicon process timed out after {self.timeout} seconds"
        except Exception as e:
            return False, "", f"Error running unicon: {str(e)}"
    
    def _validate_sdf_file(self, sdf_file: str) -> bool:
        """验证SDF文件是否有效"""
        try:
            with open(sdf_file, 'r') as f:
                content = f.read()
                
            # 检查基本SDF结构
            if not content.strip():
                print(f"Warning: Empty SDF file {sdf_file}")
                return False
                
            lines = content.strip().split('\n')
            if len(lines) < 4:
                print(f"Warning: SDF file {sdf_file} too short")
                return False
                
            # 查找计数行和M  END标记
            has_count_line = False
            has_end_marker = False
            n_atoms = 0
            n_bonds = 0
            
            for i, line in enumerate(lines):
                # 寻找计数行格式: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
                if len(line) >= 6 and not has_count_line:
                    try:
                        # 尝试解析前6个字符作为原子数和键数
                        atoms_str = line[:3].strip()
                        bonds_str = line[3:6].strip()
                        if atoms_str.isdigit() and bonds_str.isdigit():
                            n_atoms = int(atoms_str)
                            n_bonds = int(bonds_str)
                            if n_atoms > 0:
                                has_count_line = True
                                print(f"Found {n_atoms} atoms, {n_bonds} bonds in {sdf_file}")
                    except (ValueError, IndexError):
                        continue
                
                # 检查M  END标记
                if line.strip() == "M  END":
                    has_end_marker = True
            
            if not has_count_line:
                print(f"Warning: No valid count line found in {sdf_file}")
                return False
                
            if not has_end_marker:
                print(f"Warning: No M  END marker found in {sdf_file}")
                return False
            
            # 如果没有键，给出警告但仍然接受文件
            if n_bonds == 0:
                print(f"Warning: No bonds found in {sdf_file} - unicon may not have generated bond information")
                
            return True
            
        except Exception as e:
            print(f"Error validating SDF file {sdf_file}: {e}")
            return False
    

    
    def process_molecules_batch(self, input_files: List[str], output_dir: str = None, 
                               max_workers: int = 4) -> Dict[str, Tuple[bool, str, str]]:
        """批量处理分子"""
        if output_dir:
            Path(output_dir).mkdir(exist_ok=True)
        
        results = {}
        
        print(f"\nProcessing {len(input_files)} molecules with unicon (max_workers={max_workers}):")
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_file = {}
            for input_file in input_files:
                if output_dir:
                    input_path = Path(input_file)
                    output_file = Path(output_dir) / f"{input_path.stem}.sdf"
                else:
                    output_file = None
                
                future = executor.submit(self.process_molecule, input_file, str(output_file) if output_file else None)
                future_to_file[future] = input_file
            
            completed = 0
            for future in as_completed(future_to_file):
                input_file = future_to_file[future]
                try:
                    success, output_file, error_msg = future.result()
                    results[input_file] = (success, output_file, error_msg)
                    
                    completed += 1
                    status = "✓" if success else "✗"
                    filename = Path(input_file).name
                    print(f"  [{completed}/{len(input_files)}] {status} {filename}")
                    
                    if not success and error_msg:
                        print(f"    Error: {error_msg}")
                
                except Exception as e:
                    results[input_file] = (False, "", f"Exception: {str(e)}")
                    completed += 1
                    print(f"  [{completed}/{len(input_files)}] ✗ {Path(input_file).name} (Exception)")
        
        return results



def parse_sdf_file(sdf_file: str) -> Tuple[List[str], List[np.ndarray]]:
    """解析SDF文件，返回原子符号和坐标"""
    symbols = []
    positions = []
    
    try:
        with open(sdf_file, 'r') as f:
            lines = f.readlines()
        
        # 找到计数行
        count_line_idx = -1
        for i, line in enumerate(lines):
            if len(line.strip()) >= 6:
                try:
                    atoms_str = line[:3].strip()
                    bonds_str = line[3:6].strip()
                    if atoms_str.isdigit() and bonds_str.isdigit():
                        n_atoms = int(atoms_str)
                        if n_atoms > 0:
                            count_line_idx = i
                            break
                except (ValueError, IndexError):
                    continue
        
        if count_line_idx == -1:
            print(f"Warning: Could not find count line in {sdf_file}")
            return symbols, positions
        
        n_atoms = int(lines[count_line_idx][:3].strip())
        
        # 解析原子块
        atom_start = count_line_idx + 1
        for i in range(atom_start, min(atom_start + n_atoms, len(lines))):
            line = lines[i].strip()
            if len(line) < 31:  # 最小SDF原子行长度
                continue
                
            try:
                x = float(line[:10].strip())
                y = float(line[10:20].strip())
                z = float(line[20:30].strip())
                symbol = line[31:34].strip()
                
                symbols.append(symbol)
                positions.append(np.array([x, y, z]))
                
            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse atom line in {sdf_file}: {line}")
                continue
        
        print(f"Parsed {len(symbols)} atoms from {Path(sdf_file).name}")
        return symbols, positions
        
    except Exception as e:
        print(f"Error parsing SDF file {sdf_file}: {e}")
        return symbols, positions


def reconstruct_unit_cell_from_sdf_files(sdf_files: List[str], original_cell: np.ndarray, 
                                        original_molecules: List[Molecule]) -> Atoms:
    """从SDF文件重建单晶胞结构，保留原始坐标"""
    all_symbols = []
    all_positions = []
    
    print(f"Reconstructing unit cell from {len(sdf_files)} SDF files...")
    
    # 为每个成功处理的SDF文件解析原子信息
    for i, sdf_file in enumerate(sdf_files):
        symbols, positions = parse_sdf_file(sdf_file)
        
        if symbols and positions:
            # 找到对应的原始分子（通过文件名中的分子ID）
            try:
                # 从文件名提取分子ID，例如 molecule_001.sdf -> 1
                filename = Path(sdf_file).stem
                if 'molecule_' in filename:
                    mol_id_str = filename.split('molecule_')[1]
                    mol_id = int(mol_id_str)
                    
                    # 找到对应的原始分子
                    original_mol = None
                    for mol in original_molecules:
                        if mol.molecule_id == mol_id:
                            original_mol = mol
                            break
                    
                    if original_mol:
                        # 使用原始分子的坐标，但保留SDF中的原子顺序和符号
                        print(f"  Using original coordinates for molecule {mol_id} ({len(symbols)} atoms)")
                        
                        # 确保原子数量匹配
                        if len(symbols) == len(original_mol.atoms):
                            for j, atom in enumerate(original_mol.atoms):
                                all_symbols.append(atom.symbol)
                                all_positions.append(atom.get_coordinates())
                        else:
                            print(f"  Warning: Atom count mismatch for molecule {mol_id}")
                            print(f"    SDF: {len(symbols)} atoms, Original: {len(original_mol.atoms)} atoms")
                            # 使用SDF坐标作为备选
                            all_symbols.extend(symbols)
                            all_positions.extend(positions)
                    else:
                        print(f"  Warning: Could not find original molecule {mol_id}")
                        # 使用SDF坐标
                        all_symbols.extend(symbols)
                        all_positions.extend(positions)
                        
                else:
                    print(f"  Warning: Could not extract molecule ID from {filename}")
                    # 使用SDF坐标
                    all_symbols.extend(symbols)
                    all_positions.extend(positions)
                    
            except (ValueError, IndexError) as e:
                print(f"  Warning: Error processing {sdf_file}: {e}")
                # 使用SDF坐标作为备选
                all_symbols.extend(symbols)
                all_positions.extend(positions)
    
    if not all_symbols:
        raise ValueError("No atoms found in SDF files")
    
    # 创建重建的Atoms对象
    reconstructed_atoms = Atoms(
        symbols=all_symbols,
        positions=np.array(all_positions),
        cell=original_cell,
        pbc=[True, True, True]
    )
    
    return reconstructed_atoms


def merge_sdf_files(sdf_files: List[str], output_file: str) -> bool:
    """合并多个SDF文件为一个"""
    print(f"\nMerging {len(sdf_files)} SDF files into {output_file}")
    
    try:
        with open(output_file, 'w') as outf:
            for i, sdf_file in enumerate(sdf_files):
                if os.path.exists(sdf_file):
                    with open(sdf_file, 'r') as inf:
                        content = inf.read().strip()
                        if content:
                            outf.write(content)
                            if not content.endswith('\n'):
                                outf.write('\n')
                            # 确保正确的SDF分隔符
                            if not content.endswith('$$$$'):
                                outf.write('$$$$\n')
                            else:
                                outf.write('\n')
                    print(f"  Merged: {Path(sdf_file).name}")
                else:
                    print(f"  Warning: File not found: {sdf_file}")
        
        print(f"✓ Successfully merged into {output_file}")
        return True
        
    except Exception as e:
        print(f"Error merging SDF files: {e}")
        return False


def make_supercell(prim, P, wrap=True, order="cell-major", tol=1e-5):
    """通过应用变换矩阵P生成超胞"""
    
    def lattice_points_in_supercell(supercell_matrix):
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

    supercell_matrix = P
    supercell_cell = supercell_matrix @ prim.cell

    lattice_points_frac = lattice_points_in_supercell(supercell_matrix)
    lattice_points = np.dot(lattice_points_frac, supercell_cell)
    N = len(lattice_points)

    if order == "cell-major":
        shifted = prim.positions[None, :, :] + lattice_points[:, None, :]
    elif order == "atom-major":
        shifted = prim.positions[:, None, :] + lattice_points[None, :, :]
    else:
        raise ValueError(f"invalid order: {order}")
    
    shifted_reshaped = shifted.reshape(-1, 3)

    superatoms = Atoms(
        symbols=prim.symbols * N,
        positions=shifted_reshaped,
        cell=supercell_cell,
        pbc=prim.pbc
    )

    for name, arr in prim.arrays.items():
        if name == "positions":
            continue
        shape = (N * arr.shape[0], *arr.shape[1:])
        if order == "cell-major":
            new_arr = np.repeat(arr[None, :], N, axis=0).reshape(shape)
        elif order == "atom-major":
            new_arr = np.repeat(arr[:, None], N, axis=1).reshape(shape)
        superatoms.set_array(name, new_arr)

    if wrap:
        superatoms.wrap(eps=tol)

    return superatoms
def atoms_to_sdf_string(atoms: Atoms, title: str = "Supercell Structure") -> str:
    """将Atoms对象转换为SDF格式字符串"""
    lines = []
    
    # SDF头部
    lines.append(title)
    lines.append("  Generated by CIF to Supercell SDF converter")
    lines.append("")
    
    # 计数行：原子数 键数 0 0 0 0 0 0 0 0999 V2000
    num_atoms = len(atoms)
    lines.append(f"{num_atoms:3d}  0  0  0  0  0  0  0  0  0999 V2000")
    
    # 原子块
    for symbol, pos in zip(atoms.symbols, atoms.positions):
        lines.append(f"{pos[0]:10.4f}{pos[1]:10.4f}{pos[2]:10.4f} {symbol:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
    
    # 结束标记
    lines.append("M  END")
    
    return "\n".join(lines)


def write_pdb_with_bonds(filename: str, atoms: Atoms, sdf_files: List[str] = None):
    """写入PDB格式文件，包含键连信息注释"""
    
    with open(filename, 'w') as f:
        f.write("TITLE     Generated supercell structure with bond information\n")
        f.write("REMARK    Created by CIF to Supercell converter with Unicon processing\n")
        f.write("REMARK    Original coordinates preserved from unit cell\n")
        
        if sdf_files:
            f.write(f"REMARK    Bond information processed for {len(sdf_files)} molecules\n")
            f.write("REMARK    Individual molecule SDF files contain detailed bond information\n")
        else:
            f.write("REMARK    No bond information available (Unicon processing skipped)\n")
        
        f.write(f"REMARK    Total atoms: {len(atoms)}\n")
        f.write(f"REMARK    Formula: {atoms.get_chemical_formula()}\n")
        
        # 写入CRYST1记录
        if atoms.cell is not None and atoms.cell.shape == (3, 3):
            cell_obj = Cell(atoms.cell)
            cellpar = cell_obj.cellpar()
            f.write(f"CRYST1{cellpar[0]:9.3f}{cellpar[1]:9.3f}{cellpar[2]:9.3f}"
                   f"{cellpar[3]:7.2f}{cellpar[4]:7.2f}{cellpar[5]:7.2f} P 1           1\n")
        
        # 写入原子
        for i, (symbol, pos) in enumerate(zip(atoms.symbols, atoms.positions)):
            atom_name = f"{symbol:<4}"
            res_name = "UNL"
            chain_id = "A"
            res_seq = 1
            
            f.write(f"HETATM{i+1:5d} {atom_name} {res_name} {chain_id}{res_seq:4d}    "
                   f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}"
                   f"{1.00:6.2f}{0.00:6.2f}          {symbol:>2}\n")
        
        f.write("END\n")


def write_sdf(filename: str, atoms: Atoms, title: str = None):
    """将atoms写入SDF格式文件"""
    if title is None:
        title = f"Supercell Structure ({atoms.get_chemical_formula()})"
    
    sdf_content = atoms_to_sdf_string(atoms, title)
    
    with open(filename, 'w') as f:
        f.write(sdf_content)
        f.write('\n$$$$\n')  # SDF分隔符


def process_cif_to_supercell_sdf(cif_path: str, output_path: str = None, 
                                supercell_size: Tuple[int, int, int] = (3, 3, 3),
                                bond_cutoff: float = 2.0, hydrogen_cutoff: float = 1.3,
                                unicon_path: str = "unicon", max_workers: int = 4,
                                keep_temp: bool = False, timeout: int = 300,
                                skip_unicon: bool = False) -> bool:
    """主处理函数"""
    
    print(f"Loading CIF file: {cif_path}")
    
    # 1. 读取CIF文件并生成完整的单晶胞
    try:
        unit_cell_atoms = read_cif(cif_path)
        print(f"Successfully loaded unit cell structure:")
        print(f"  - Formula: {unit_cell_atoms.get_chemical_formula()}")
        print(f"  - Number of atoms: {len(unit_cell_atoms)}")
        
        if unit_cell_atoms.cell is not None and unit_cell_atoms.cell.shape == (3, 3):
            cell_obj = Cell(unit_cell_atoms.cell)
            cellpar = cell_obj.cellpar()
            print(f"  - Cell parameters: {cellpar}")
        
    except Exception as e:
        print(f"Error loading CIF file: {e}")
        return False
    
    # 2. 分离分子
    print(f"\n{'='*50}")
    print("Step 1: Splitting molecules from unit cell")
    print(f"{'='*50}")
    
    splitter = MoleculeSplitter(
        bond_distance_cutoff=bond_cutoff,
        hydrogen_cutoff=hydrogen_cutoff
    )
    
    try:
        splitter.load_from_atoms_object(unit_cell_atoms)
        molecules = splitter.split_molecules()
        
        if not molecules:
            print("No molecules found!")
            return False
        
        # 创建临时目录保存分子PDB文件
        temp_dir = tempfile.mkdtemp(prefix="cif_to_sdf_")
        print(f"Using temporary directory: {temp_dir}")
        
        pdb_files = splitter.save_molecules(temp_dir)
        
    except Exception as e:
        print(f"Error during molecule splitting: {e}")
        return False
    
    # 3. 使用Unicon处理分子
    print(f"\n{'='*50}")
    print("Step 2: Processing molecules with Unicon")
    print(f"{'='*50}")
    
    if skip_unicon:
        print("Skipping Unicon processing as requested")
        successful_sdf_files = []
    else:
        processor = UniconProcessor(unicon_path=unicon_path, timeout=timeout)
        
        if not processor.check_unicon_availability():
            print(f"Error: unicon not found at '{unicon_path}'")
            print("Please make sure unicon is installed and accessible")
            print("Common solutions:")
            print("  1. Install unicon and add it to your PATH")
            print("  2. Use --unicon-path to specify the full path to unicon")
            print("  3. Use --skip-unicon to generate supercell without bond information")
            sys.exit(1)
    
    if not skip_unicon:
        try:
            unicon_output_dir = os.path.join(temp_dir, "unicon_output")
            
            results = processor.process_molecules_batch(
                pdb_files, 
                unicon_output_dir,
                max_workers
            )
            
            successful = sum(1 for success, _, _ in results.values() if success)
            failed = len(results) - successful
            
            print(f"\nUnicon processing completed:")
            print(f"  Successful: {successful}")
            print(f"  Failed: {failed}")
            
            if failed > 0:
                print("\nFailed files:")
                for input_file, (success, _, error_msg) in results.items():
                    if not success:
                        print(f"  {Path(input_file).name}: {error_msg}")
            
            # 收集成功处理的SDF文件
            successful_sdf_files = []
            for input_file, (success, output_file, _) in results.items():
                if success and os.path.exists(output_file):
                    successful_sdf_files.append(output_file)
            
            if successful == 0:
                print("No molecules were successfully processed by Unicon!")
                print("This might be due to:")
                print("  - Molecules containing only hydrogen atoms")
                print("  - Invalid molecular structures")
                print("  - Unicon configuration issues")
                print("Continuing without bond information...")
            
        except Exception as e:
            print(f"Error during unicon processing: {e}")
            print("Continuing without bond information...")
            successful_sdf_files = []
    
    # 4. 合并处理后的分子
    print(f"\n{'='*50}")
    print("Step 3: Merging processed molecules")
    print(f"{'='*50}")
    
    merged_sdf_path = None
    if successful_sdf_files:
        try:
            merged_sdf_path = os.path.join(temp_dir, "merged_molecules.sdf")
            if not merge_sdf_files(successful_sdf_files, merged_sdf_path):
                print("Failed to merge SDF files!")
                merged_sdf_path = None
            
        except Exception as e:
            print(f"Error merging molecules: {e}")
            merged_sdf_path = None
    else:
        print("No processed molecules to merge (Unicon was skipped or failed)")
    
    # 5. 读取合并的分子并重建Atoms对象
    print(f"\n{'='*50}")
    print("Step 4: Reconstructing unit cell from processed molecules")
    print(f"{'='*50}")
    
    try:
        if successful_sdf_files:
            # 从成功处理的SDF文件重建单晶胞结构
            reconstructed_atoms = reconstruct_unit_cell_from_sdf_files(
                successful_sdf_files, unit_cell_atoms.cell, molecules
            )
            print(f"Successfully reconstructed unit cell with {len(reconstructed_atoms)} atoms")
            print(f"  - Original unit cell: {len(unit_cell_atoms)} atoms")
            print(f"  - Reconstructed: {len(reconstructed_atoms)} atoms")
            print(f"  - Formula: {reconstructed_atoms.get_chemical_formula()}")
            
            # 使用重建的结构生成超胞
            unit_cell_for_supercell = reconstructed_atoms
        else:
            print("No processed molecules available, using original unit cell")
            print("(No bond information will be included)")
            unit_cell_for_supercell = unit_cell_atoms
        
    except Exception as e:
        print(f"Error reconstructing unit cell: {e}")
        print("Falling back to original unit cell")
        unit_cell_for_supercell = unit_cell_atoms
    
    # 6. 生成超胞
    print(f"\n{'='*50}")
    print(f"Step 5: Generating {supercell_size[0]}×{supercell_size[1]}×{supercell_size[2]} supercell")
    print(f"{'='*50}")
    
    try:
        P = np.array([
            [supercell_size[0], 0, 0], 
            [0, supercell_size[1], 0], 
            [0, 0, supercell_size[2]]
        ], dtype=int)
        
        supercell = make_supercell(unit_cell_for_supercell, P)
        
        print(f"Supercell created successfully:")
        print(f"  - Unit cell atoms: {len(unit_cell_for_supercell)}")
        print(f"  - Supercell atoms: {len(supercell)}")
        print(f"  - Expected atoms: {len(unit_cell_for_supercell) * np.prod(supercell_size)}")
        print(f"  - Formula: {supercell.get_chemical_formula()}")
        
    except Exception as e:
        print(f"Error creating supercell: {e}")
        return False
    
    # 7. 输出为SDF格式
    print(f"\n{'='*50}")
    print("Step 6: Writing final SDF file")
    print(f"{'='*50}")
    
    if output_path is None:
        cif_path_obj = Path(cif_path)
        output_path = cif_path_obj.with_suffix('').name + f'_supercell_{supercell_size[0]}x{supercell_size[1]}x{supercell_size[2]}.sdf'
    
    try:
        write_sdf(output_path, supercell, 
                 f"Supercell {supercell_size[0]}×{supercell_size[1]}×{supercell_size[2]} from {Path(cif_path).name}")
        
        print(f"✓ Successfully exported supercell to: {output_path}")
        
        # 生成PDB文件，包含所有成功处理的分子
        pdb_output_path = Path(output_path).with_suffix('.pdb')
        write_pdb_with_bonds(pdb_output_path, supercell, successful_sdf_files if not skip_unicon else [])
        print(f"✓ Supercell PDB file with bond information saved to: {pdb_output_path}")
        
        # 同时保存合并的分子SDF文件到最终输出目录（如果存在）
        if merged_sdf_path and os.path.exists(merged_sdf_path):
            merged_output_path = Path(output_path).with_suffix('').name + '_molecules_with_bonds.sdf'
            shutil.copy2(merged_sdf_path, merged_output_path)
            print(f"✓ Molecules with bond information saved to: {merged_output_path}")
        else:
            print("Note: No bond information available (Unicon processing was skipped or failed)")
        
    except Exception as e:
        print(f"Error writing output files: {e}")
        return False
    
    # 8. 清理临时文件
    if not keep_temp:
        print(f"\nCleaning up temporary files...")
        try:
            shutil.rmtree(temp_dir)
            print("✓ Temporary files cleaned up")
        except Exception as e:
            print(f"Warning: Could not clean up temporary files: {e}")
    else:
        print(f"\nTemporary files kept in: {temp_dir}")
    
    return True


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="Convert CIF file to supercell SDF with bond information via Unicon",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python cif_to_supercell_sdf.py input.cif
  python cif_to_supercell_sdf.py input.cif -o output.sdf
  python cif_to_supercell_sdf.py input.cif --supercell 2 2 2
  python cif_to_supercell_sdf.py input.cif --bond-cutoff 2.5 --keep-temp

Process:
  1. Parse CIF file and generate complete unit cell
  2. Split unit cell into individual molecules
  3. Process each molecule with Unicon to assign bonds
  4. Merge processed molecules
  5. Generate supercell
  6. Output as RDKit-compatible SDF format
        """
    )
    
    parser.add_argument("input_cif", help="Input CIF file")
    parser.add_argument("-o", "--output", help="Output SDF file (auto-generated if not specified)")
    
    parser.add_argument("--supercell", nargs=3, type=int, default=[3, 3, 3],
                       metavar=('X', 'Y', 'Z'),
                       help="Supercell dimensions (default: 3 3 3)")
    
    parser.add_argument("--bond-cutoff", type=float, default=2.0,
                       help="Bond distance cutoff for non-hydrogen atoms in Å (default: 2.0)")
    parser.add_argument("--hydrogen-cutoff", type=float, default=1.3,
                       help="Bond distance cutoff for hydrogen atoms in Å (default: 1.3)")
    
    parser.add_argument("--unicon-path", default="unicon",
                       help="Path to unicon executable (default: unicon)")
    parser.add_argument("--max-workers", type=int, default=4,
                       help="Maximum number of parallel workers for unicon (default: 4)")
    parser.add_argument("--timeout", type=int, default=300,
                       help="Timeout for each unicon process in seconds (default: 300)")
    
    parser.add_argument("--keep-temp", action="store_true",
                       help="Keep temporary files for debugging")
    parser.add_argument("--skip-unicon", action="store_true",
                       help="Skip Unicon processing and generate supercell without bond information")

    
    args = parser.parse_args()
    
    # 检查输入文件
    if not os.path.exists(args.input_cif):
        print(f"Error: Input file '{args.input_cif}' not found!")
        sys.exit(1)
    
    print("CIF to Supercell SDF Converter with Unicon Bond Assignment")
    print("=" * 60)
    print(f"Input CIF: {args.input_cif}")
    print(f"Output SDF: {args.output if args.output else 'Auto-generated'}")
    print(f"Supercell size: {args.supercell[0]}×{args.supercell[1]}×{args.supercell[2]}")
    print(f"Bond cutoffs: {args.bond_cutoff} Å (non-H), {args.hydrogen_cutoff} Å (H)")
    print(f"Unicon path: {args.unicon_path}")
    print(f"Max workers: {args.max_workers}")
    print()
    
    start_time = time.time()
    
    try:
        success = process_cif_to_supercell_sdf(
            cif_path=args.input_cif,
            output_path=args.output,
            supercell_size=tuple(args.supercell),
            bond_cutoff=args.bond_cutoff,
            hydrogen_cutoff=args.hydrogen_cutoff,
            unicon_path=args.unicon_path,
            max_workers=args.max_workers,
            keep_temp=args.keep_temp,
            timeout=args.timeout,
            skip_unicon=getattr(args, 'skip_unicon', False)
        )
        
        elapsed_time = time.time() - start_time
        
        if success:
            print(f"\n{'='*60}")
            print(f"✓ Conversion completed successfully in {elapsed_time:.2f} seconds!")
            print(f"✓ Generated files:")
            print(f"  - Supercell SDF: Ready for RDKit analysis")
            print(f"  - Supercell PDB: Contains all successfully processed molecules with original coordinates")
            if not getattr(args, 'skip_unicon', False):
                print(f"  - Individual molecule SDF files: Contain detailed bond information")
        else:
            print(f"\n{'='*60}")
            print(f"✗ Conversion failed after {elapsed_time:.2f} seconds!")
            sys.exit(1)
            
    except KeyboardInterrupt:
        print(f"\n\nProcess interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nUnexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()