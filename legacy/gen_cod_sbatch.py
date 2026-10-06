from pathlib import Path
import os
from all_frags import REF_FRAGMENTS, PROBE_FRAGMENTS

home_dir = Path.home()
splinter_frag_pairs_dir = home_dir / "2507/output"
scripts_dir = home_dir / "2507/codes"

for ref_frag in REF_FRAGMENTS:
    ref_frag_dir = splinter_frag_pairs_dir / ref_frag
    os.makedirs(ref_frag_dir, exist_ok=True)

    for probe_frag in PROBE_FRAGMENTS:
        pair_name = f"{ref_frag}_{probe_frag}"
        pair_dir = ref_frag_dir / pair_name
        os.makedirs(pair_dir, exist_ok=True)

        csd_sbatch_script = ''
        csd_sbatch_script += f'#!/bin/bash\n'
        csd_sbatch_script += f'#SBATCH --job-name={pair_name}\n'
        csd_sbatch_script += f'#SBATCH --output={pair_dir}/cod_workflow.log\n'
        csd_sbatch_script += f'#SBATCH --partition=honda\n'
        csd_sbatch_script += f'#SBATCH --nodelist=n152\n'
        csd_sbatch_script += f'#SBATCH --ntasks=1\n'
        csd_sbatch_script += f'#SBATCH --qos=short_many\n'
        csd_sbatch_script += f'#SBATCH --time=200:00\n'
        csd_sbatch_script += f'source ~/.bashrc\n'
        csd_sbatch_script += f'conda activate csd\n'
        csd_sbatch_script += f'date\n'
        csd_sbatch_script += f'hostname\n'
        csd_sbatch_script += f'python {scripts_dir}/cod_workflow.py --frag1 {ref_frag} --frag2 {probe_frag}\n'
        csd_sbatch_script += f'date\n'  

        csd_sbatch_script_path = pair_dir / 'cod_workflow.sbatch'
        with open(csd_sbatch_script_path, 'w') as f:
            f.write(csd_sbatch_script)

print("COD workflow scripts generated successfully!")
