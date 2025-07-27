from ase.io import read

cif_file = '/home/clasper/2507/cif/1/00/40/1004002.cif'
atoms = read(cif_file)

symbols = atoms.get_chemical_symbols()  # 原子类型列表
positions = atoms.get_positions()

print(positions)