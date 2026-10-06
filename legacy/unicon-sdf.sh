#!/bin/bash

find /home/clasper/2507/tmp -type f -name '*.pdb' | LC_COLLATE=C sort | while read -r line; do
    line=`readlink -f $line`
    pdb_path=`dirname $line`
    pdb_file=`basename $line`
    sdf_file=`echo ${pdb_file%.pdb}.sdf`
    cd $pdb_path
    echo $pdb_file
    echo $sdf_file

    unicon -i $pdb_file -o $sdf_file -t single -p single -g 1
done

echo "Generated unicon-sdf files. "