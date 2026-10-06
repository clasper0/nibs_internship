#!/bin/bash
cwd=`pwd`
# find "/home/yqwu/2507/output" -maxdepth 3 -type f -name 'cod_workflow.sbatch' | LC_COLLATE=C sort > sbatch_cod.list
cat sbatch_cod.list | while read -r line; do
    cal_cod_workflow_path=`readlink -f $line`
    cal_cod_workflow_dir=`dirname $cal_cod_workflow_path`
    cd $cal_cod_workflow_dir

    # if [[ -n "$(find . -maxdepth 1 -type f -name '*.json.gz')" ]]; then
    #     cd $cwd
    #     continue
    # fi

    # while [[ $(squeue -u "ylzhang" | grep -E 'honda' | wc -l) -ge 50 ]]; do
    #     sleep 30s
    # done

    sbatch csd_workflow.sbatch

    cd $cwd

done
