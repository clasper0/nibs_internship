from rdkit import Chem
from pathlib import Path
import sys

def split_sdf(sdf_file, cod_id):
    tmp_dir = Path("/home/clasper/2507/tmp")
    structure = Chem.SDMolSupplier(sdf_file, removeHs=False, sanitize=False, strictParsing=False)[0]

    mols = Chem.GetMolFrags(structure, asMols=True, sanitizeFrags=Chem.SanitizeFlags.SANITIZE_NONE)
    i = 1
    for mol in mols:
        Chem.MolToPDBFile(mol, f"{str(tmp_dir)}/{cod_id}_{i}.pdb")
        i += 1
    print(f"Generated {i} pdb files. ")


if __name__ == '__main__':
    sdf_dir = Path("/home/clasper/2507/sdf")
    sdf_name = sys.argv[1].split("/")[-1]
    cod_id = sdf_name[0:7]
    sdf_file = sdf_dir / sdf_name[0] / sdf_name[1:3] / sdf_name[3:5] / sdf_name
    split_sdf(sdf_file, cod_id)