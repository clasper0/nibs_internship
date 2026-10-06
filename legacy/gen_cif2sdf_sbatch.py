from pathlib import Path
import os

home_dir = Path.home()
cif_dir = home_dir / "2507/cif"
scripts_dir = home_dir / "2507/codes"
tmp_dir = home_dir / "2507/tmp"
cif2pdb_file = scripts_dir / "cif2pdb.py"
sdf2pdbs_file = scripts_dir / "sdf2pdbs.py"
unicon_file = scripts_dir / "unicon-sdf.sh"
getfinalsdf_file = scripts_dir / "getfinalsdf.py"
cod_id_file = scripts_dir / "cod_ids.txt"

cod_ids = []
with open(cod_id_file, "r", encoding="utf-8") as f:
    for line in f.readlines():
        cod_id = line.strip()
        cod_ids.append(cod_id)

for cod_id in cod_ids:
    cif_file = cif_dir / cod_id[0] / cod_id[1:3] / cod_id[3:5] / f"{cod_id}.cif"
    sdf_dir = home_dir / f"2507/sdf/{cod_id[0]}/{cod_id[1:3]}/{cod_id[3:5]}"
    pdb_file = sdf_dir / f"{cod_id}.pdb"
    sdf_file = sdf_dir / f"{cod_id}.sdf"
    os.makedirs(sdf_dir, exist_ok=True)

    csd_sbatch_script = ''
    csd_sbatch_script += f'#!/bin/bash\n'
    csd_sbatch_script += f'#SBATCH --job-name={cod_id}\n'
    csd_sbatch_script += f'#SBATCH --output={sdf_dir}/{cod_id}.log\n'
    csd_sbatch_script += f'#SBATCH --partition=honda\n'
    csd_sbatch_script += f'#SBATCH --ntasks=1\n'
    csd_sbatch_script += f'#SBATCH --qos=short_many\n'
    csd_sbatch_script += f'#SBATCH --time=200:00\n'
    csd_sbatch_script += f'source ~/.bashrc\n'
    csd_sbatch_script += f'source /mnt/e/Users/clasper/tmp/schrodinger2021-2/env.sh\n'
    csd_sbatch_script += f'conda activate csd\n'
    csd_sbatch_script += f'mkdir -p /home/clasper/2507/tmp\n'
    csd_sbatch_script += f'date\n'
    csd_sbatch_script += f'hostname\n'
    csd_sbatch_script += f'python {cif2pdb_file} {cif_file} {pdb_file}\n'
    csd_sbatch_script += f'structconvert {pdb_file} {sdf_file}\n'
    csd_sbatch_script += f'python {sdf2pdbs_file} {sdf_file}\n'
    csd_sbatch_script += f'{unicon_file}\n'
    csd_sbatch_script += f'python {getfinalsdf_file} {cif_file}\n'
    csd_sbatch_script += f'rm {pdb_file}\n'
    csd_sbatch_script += f'rm -r {tmp_dir}\n'
    csd_sbatch_script += f'date\n'

    sbatch_path = Path("/mnt/e/Users/clasper/tmp/sdf") / f"{cod_id[0]}/{cod_id[1:3]}/{cod_id[3:5]}"
    sbatch_path.mkdir(parents=True, exist_ok=True)
    cif2sdf_sbatch_script_path = Path("/mnt/e/Users/clasper/tmp/sdf") / f"{cod_id[0]}/{cod_id[1:3]}/{cod_id[3:5]}/{cod_id}.sbatch"
    with open(cif2sdf_sbatch_script_path, 'w') as f:
        f.write(csd_sbatch_script)

print("COD workflow scripts generated successfully!")
