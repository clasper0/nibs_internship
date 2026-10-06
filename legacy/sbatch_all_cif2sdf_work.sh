#!/bin/bash
cwd=`pwd`
# find "/mnt/e/Users/clasper/tmp/sdf" -maxdepth 4 -type f -name '*.sbatch' > sbatch_cif2sdf.list
cat sbatch_cif2sdf.list | while read -r line; do
    cif2sdf_sbatch_path=`readlink -f $line`
    cif2sdf_path=`dirname $cif2sdf_sbatch_path`
    cif2sdf_sbatch=`basename $cif2sdf_sbatch_path`
    cd $cif2sdf_path

    # if [[ -n "$(find . -maxdepth 1 -type f -name '*.json.gz')" ]]; then
    #     cd $cwd
    #     continue
    # fi

    # while [[ $(squeue -u "ylzhang" | grep -E 'honda' | wc -l) -ge 50 ]]; do
    #     sleep 30s
    # done

    sbatch $cif2sdf_sbatch

    cd $cwd

done
