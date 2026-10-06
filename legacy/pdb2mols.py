#!/usr/bin/env python3
"""
PDB分子分离器和Unicon处理器

该脚本从含有多个分子的PDB文件中分离出单独的分子，
然后使用unicon为每个分子添加键连信息。

功能特点：
1. 基于距离算法自动检测分子边界
2. 为每个分子生成单独的PDB文件
3. 自动调用unicon处理每个分子
4. 支持批量处理和并行执行
5. 提供详细的进度信息和统计

使用方法:
    python pdb_splitter_unicon.py input.pdb [options]

依赖:
    - numpy
    - unicon (外部程序)
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Set
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
import tempfile
import shutil
import time


class Atom:
    """原子类，存储原子信息"""
    
    def __init__(self, line: str):
        """从PDB行解析原子信息"""
        self.record_type = line[0:6].strip()
        self.atom_number = int(line[6:11].strip())
        self.atom_name = line[12:16].strip()
        self.residue_name = line[17:20].strip()
        self.chain_id = line[21:22].strip()
        self.residue_number = int(line[22:26].strip()) if line[22:26].strip() else 1
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip())
        self.occupancy = float(line[54:60].strip()) if line[54:60].strip() else 1.0
        self.temperature_factor = float(line[60:66].strip()) if line[60:66].strip() else 0.0
        self.element = line[76:78].strip() if len(line) > 76 else self.atom_name[0]
        
        # 原始行，用于输出
        self.original_line = line.rstrip()
    
    def get_coordinates(self) -> np.ndarray:
        """获取原子坐标"""
        return np.array([self.x, self.y, self.z])
    
    def distance_to(self, other_atom: 'Atom') -> float:
        """计算到另一个原子的距离"""
        coord1 = self.get_coordinates()
        coord2 = other_atom.get_coordinates()
        return np.linalg.norm(coord1 - coord2)
    
    def to_pdb_line(self, atom_number: Optional[int] = None) -> str:
        """转换为PDB格式行"""
        if atom_number is not None:
            # 更新原子编号
            line = self.original_line
            new_line = f"{line[:6]}{atom_number:5d}{line[11:]}"
            return new_line
        return self.original_line


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
        elements = [atom.element for atom in self.atoms]
        element_count = Counter(elements)
        
        formula = ""
        # 通常的元素顺序：C, H, N, O, 然后按字母顺序
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
        
        # 添加头部信息
        if title:
            lines.append(f"TITLE     {title}")
        else:
            lines.append(f"TITLE     Molecule {self.molecule_id} ({self.get_formula()})")
        
        lines.append("REMARK    Extracted from supercell structure")
        lines.append(f"REMARK    Molecule ID: {self.molecule_id}")
        lines.append(f"REMARK    Number of atoms: {len(self.atoms)}")
        lines.append(f"REMARK    Molecular formula: {self.get_formula()}")
        
        # 添加原子信息，重新编号
        for i, atom in enumerate(self.atoms, 1):
            lines.append(atom.to_pdb_line(i))
        
        lines.append("END")
        
        return "\n".join(lines)


class PDBMoleculeSplitter:
    """PDB分子分离器"""
    
    def __init__(self, bond_distance_cutoff: float = 2.0, hydrogen_cutoff: float = 1.3):
        """
        初始化分离器
        
        Args:
            bond_distance_cutoff: 非氢原子间的键长阈值
            hydrogen_cutoff: 氢原子的键长阈值
        """
        self.bond_distance_cutoff = bond_distance_cutoff
        self.hydrogen_cutoff = hydrogen_cutoff
        self.atoms = []
        self.molecules = []
    
    def read_pdb(self, filename: str) -> None:
        """读取PDB文件"""
        print(f"Reading PDB file: {filename}")
        
        self.atoms = []
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith(('HETATM', 'ATOM')):
                    try:
                        atom = Atom(line)
                        self.atoms.append(atom)
                    except (ValueError, IndexError) as e:
                        print(f"Warning: Could not parse line: {line.strip()} - {e}")
        
        print(f"Successfully read {len(self.atoms)} atoms")
    
    def _are_bonded(self, atom1: Atom, atom2: Atom) -> bool:
        """判断两个原子是否成键"""
        distance = atom1.distance_to(atom2)
        
        # 根据原子类型选择不同的阈值
        if atom1.element == 'H' or atom2.element == 'H':
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
        
        # 检查与其他未访问原子的连接
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
                # 找到这个连通分量的所有原子
                connected_indices = self._find_connected_atoms(i, visited)
                
                if connected_indices:
                    # 创建分子
                    molecule_atoms = [self.atoms[idx] for idx in connected_indices]
                    molecule = Molecule(molecule_atoms, molecule_id)
                    molecules.append(molecule)
                    molecule_id += 1
        
        self.molecules = molecules
        print(f"Found {len(molecules)} molecules")
        
        # 显示分子统计信息
        for mol in molecules:
            center = mol.get_center()
            print(f"  Molecule {mol.molecule_id}: {len(mol)} atoms, "
                  f"formula: {mol.get_formula()}, "
                  f"center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
        
        return molecules
    
    def save_molecules(self, output_dir: str, filename_pattern: str = "molecule_{:03d}.pdb") -> List[str]:
        """保存分子到单独的PDB文件"""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        saved_files = []
        
        print(f"\nSaving molecules to {output_dir}:")
        for molecule in self.molecules:
            filename = filename_pattern.format(molecule.molecule_id)
            filepath = output_path / filename
            
            with open(filepath, 'w') as f:
                f.write(molecule.to_pdb_string())
            
            saved_files.append(str(filepath))
            print(f"  Saved: {filename} ({len(molecule)} atoms, {molecule.get_formula()})")
        
        return saved_files


class UniconProcessor:
    """Unicon处理器"""
    
    def __init__(self, unicon_path: str = "unicon", 
                 output_format: str = "pdb",
                 generate_3d: bool = True,
                 max_conformers: int = 1,
                 verbosity: int = 1,
                 timeout: int = 300):
        """
        初始化Unicon处理器
        
        Args:
            unicon_path: unicon可执行文件路径
            output_format: 输出格式 (sdf, mol2等)
            generate_3d: 是否生成3D坐标
            max_conformers: 最大构象数量
            verbosity: 详细程度 (0-3)
            timeout: 超时时间（秒）
        """
        self.unicon_path = unicon_path
        self.output_format = output_format
        self.generate_3d = generate_3d
        self.max_conformers = max_conformers
        self.verbosity = verbosity
        self.timeout = timeout
    
    def check_unicon_availability(self) -> bool:
        """检查unicon是否可用"""
        try:
            result = subprocess.run(
                [self.unicon_path, "-h"], 
                capture_output=True, 
                text=True, 
                timeout=10
            )
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
            return False
    
    def process_molecule(self, input_file: str, output_file: str = None) -> Tuple[bool, str, str]:
        """
        使用unicon处理单个分子
        
        Args:
            input_file: 输入PDB文件路径
            output_file: 输出文件路径（可选）
        
        Returns:
            (成功标志, 输出文件路径, 错误信息)
        """
        if output_file is None:
            input_path = Path(input_file)
            output_file = str(input_path.with_suffix(f'.{self.output_format}'))
        
        # 构建unicon命令
        cmd = [
            self.unicon_path,
            "-i", input_file,
            "-o", output_file,
            # "--outFormat", self.output_format.upper(),
            "-t", "single",
            "-p", "single", 
            "-g", "1"
        ]
        
        # if self.generate_3d:
        #     cmd.extend(["-g", "3"])  # 生成3D坐标
        
        try:
            print(f"Processing {Path(input_file).name} with unicon...")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=self.timeout
            )
            
            if result.returncode == 0:
                if os.path.exists(output_file):
                    return True, output_file, ""
                else:
                    return False, "", f"Output file not created: {output_file}"
            else:
                error_msg = f"Unicon failed with return code {result.returncode}\n"
                error_msg += f"STDERR: {result.stderr}\n"
                error_msg += f"STDOUT: {result.stdout}"
                return False, "", error_msg
                
        except subprocess.TimeoutExpired:
            return False, "", f"Unicon process timed out after {self.timeout} seconds"
        except Exception as e:
            return False, "", f"Error running unicon: {str(e)}"
    
    def process_molecules_batch(self, input_files: List[str], output_dir: str = None, 
                               max_workers: int = 4) -> Dict[str, Tuple[bool, str, str]]:
        """
        批量处理分子
        
        Args:
            input_files: 输入文件列表
            output_dir: 输出目录
            max_workers: 最大并行工作线程数
        
        Returns:
            处理结果字典 {输入文件: (成功标志, 输出文件, 错误信息)}
        """
        if output_dir:
            Path(output_dir).mkdir(exist_ok=True)
        
        results = {}
        
        print(f"\nProcessing {len(input_files)} molecules with unicon (max_workers={max_workers}):")
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # 提交所有任务
            future_to_file = {}
            for input_file in input_files:
                if output_dir:
                    input_path = Path(input_file)
                    output_file = Path(output_dir) / f"{input_path.stem}.{self.output_format}"
                else:
                    output_file = None
                
                future = executor.submit(self.process_molecule, input_file, str(output_file) if output_file else None)
                future_to_file[future] = input_file
            
            # 收集结果
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


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="Split PDB file into individual molecules and process with unicon"
    )
    
    parser.add_argument("input_pdb", help="Input PDB file")
    parser.add_argument("-o", "--output-dir", default="molecules", 
                       help="Output directory for molecules (default: molecules)")
    parser.add_argument("--unicon-dir", default="unicon_output",
                       help="Output directory for unicon results (default: unicon_output)")
    parser.add_argument("--bond-cutoff", type=float, default=2.0,
                       help="Bond distance cutoff for non-hydrogen atoms (default: 2.0)")
    parser.add_argument("--hydrogen-cutoff", type=float, default=1.3,
                       help="Bond distance cutoff for hydrogen atoms (default: 1.3)")
    parser.add_argument("--unicon-path", default="unicon",
                       help="Path to unicon executable (default: unicon)")
    parser.add_argument("--output-format", default="sdf", choices=["sdf", "mol2"],
                       help="Output format for unicon (default: sdf)")
    parser.add_argument("--max-workers", type=int, default=4,
                       help="Maximum number of parallel workers for unicon (default: 4)")
    parser.add_argument("--timeout", type=int, default=300,
                       help="Timeout for each unicon process in seconds (default: 300)")
    parser.add_argument("--skip-unicon", action="store_true",
                       help="Skip unicon processing, only split molecules")
    parser.add_argument("--keep-temp", action="store_true",
                       help="Keep temporary PDB files after processing")
    
    args = parser.parse_args()
    
    # 检查输入文件
    if not os.path.exists(args.input_pdb):
        print(f"Error: Input file '{args.input_pdb}' not found!")
        sys.exit(1)
    
    print("PDB Molecule Splitter and Unicon Processor")
    print("=" * 50)
    print(f"Input PDB: {args.input_pdb}")
    print(f"Output directory: {args.output_dir}")
    print(f"Unicon output directory: {args.unicon_dir}")
    print(f"Bond cutoffs: {args.bond_cutoff} Å (non-H), {args.hydrogen_cutoff} Å (H)")
    print()
    
    start_time = time.time()
    
    # 1. 分离分子
    splitter = PDBMoleculeSplitter(
        bond_distance_cutoff=args.bond_cutoff,
        hydrogen_cutoff=args.hydrogen_cutoff
    )
    
    try:
        splitter.read_pdb(args.input_pdb)
        molecules = splitter.split_molecules()
        
        if not molecules:
            print("No molecules found!")
            sys.exit(1)
        
        # 保存分子
        pdb_files = splitter.save_molecules(args.output_dir)
        
    except Exception as e:
        print(f"Error during molecule splitting: {e}")
        sys.exit(1)
    
    # 2. 使用unicon处理
    if not args.skip_unicon:
        processor = UniconProcessor(
            unicon_path=args.unicon_path,
            output_format=args.output_format,
            generate_3d=True,
            max_conformers=1,
            verbosity=1,
            timeout=args.timeout
        )
        
        # 检查unicon可用性
        if not processor.check_unicon_availability():
            print(f"Warning: unicon not found at '{args.unicon_path}'")
            print("Please make sure unicon is installed and in your PATH")
            print("Skipping unicon processing...")
        else:
            try:
                # 批量处理
                results = processor.process_molecules_batch(
                    pdb_files, 
                    args.unicon_dir,
                    args.max_workers
                )
                
                # 统计结果
                successful = sum(1 for success, _, _ in results.values() if success)
                failed = len(results) - successful
                
                print(f"\nUnicon processing completed:")
                print(f"  Successful: {successful}")
                print(f"  Failed: {failed}")
                
                # 显示失败的文件
                if failed > 0:
                    print("\nFailed files:")
                    for input_file, (success, _, error_msg) in results.items():
                        if not success:
                            print(f"  {Path(input_file).name}: {error_msg}")
            
            except Exception as e:
                print(f"Error during unicon processing: {e}")
    
    # 3. 清理临时文件
    if not args.keep_temp and not args.skip_unicon:
        print(f"\nCleaning up temporary PDB files in {args.output_dir}...")
        try:
            shutil.rmtree(args.output_dir)
            print("✓ Temporary files cleaned up")
        except Exception as e:
            print(f"Warning: Could not clean up temporary files: {e}")
    
    # 总结
    elapsed_time = time.time() - start_time
    print(f"\n" + "=" * 50)
    print(f"Processing completed in {elapsed_time:.2f} seconds")
    print(f"Found and processed {len(molecules)} molecules")
    
    if not args.skip_unicon:
        print(f"Unicon output directory: {args.unicon_dir}")
    else:
        print(f"Molecule PDB files saved in: {args.output_dir}")


if __name__ == "__main__":
    main()