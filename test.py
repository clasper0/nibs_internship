# coding=utf-8
from openbabel import pybel
import os

def cif_to_mol(cif_file_path, mol_file_path):
    """
    将CIF文件转换为MOL文件
    cif_file_path: CIF文件的路径
    mol_file_path: MOL文件的保存路径
    """
    mol = list(pybel.readfile("cif", cif_file_path))  # 读取CIF文件
    if mol:
        if len(mol) > 1:
            print(f"Warning: More than one molecule found in {cif_file_path}. Only the first will be saved.")
        # 只处理第一个分子
        mol = mol[0]
        mol.write("mol", mol_file_path, overwrite=True)  # 写入MOL文件


# 示例调用
cif_directory = "/home/clasper/2507/cif/1/00/00"  # CIF文件所在的目录
mol_directory = "/home/clasper/2507/mol/00"  # MOL文件保存的目录

# 确保MOL目录存在
if not os.path.exists(mol_directory):
    os.makedirs(mol_directory)

# 遍历CIF目录中的每个文件
mol_file_num = 0
for cif_file in os.listdir(cif_directory):
    if cif_file.endswith(".cif"):
        cif_file_path = os.path.join(cif_directory, cif_file)   # 补全 cif 文件路径
        mol_file_name = os.path.splitext(cif_file)[0] + ".mol"
        mol_file_path = os.path.join(mol_directory, mol_file_name)  # 补全 mol 文件路径
        cif_to_mol(cif_file_path, mol_file_path)    # 进行转换
    if os.path.exists(mol_file_path):
        mol_file_num += 1
print(f"converted {mol_file_num} .cif -> .mol. ")
