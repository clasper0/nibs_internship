#!/bin/bash

# 配置
CIF2PDB_SCRIPT="/home/clasper/2507/codes/cif2pdb.py"
CIF_ROOT_DIR="/mnt/e/Users/clasper/tmp/cif"
OUTPUT_DIR="/mnt/e/Users/clasper/tmp/sdf"

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 日志
LOG_SUCCESS="$OUTPUT_DIR/success.log"
LOG_FAILED="$OUTPUT_DIR/failed.log"
LOG_SKIPPED="$OUTPUT_DIR/skipped.log"

# 清空日志
> "$LOG_SUCCESS"
> "$LOG_FAILED"
> "$LOG_SKIPPED"

# 准备环境
eval "$(micromamba shell hook --shell bash)"
micromamba activate csd
source /mnt/e/Users/clasper/tmp/schrodinger2021-2/env.sh

# 准备转换
echo "converting CIF -> SDF"
echo "CIF dir: $CIF_ROOT_DIR"
echo "SDF dir: $OUTPUT_DIR"
printf '%.0s=' {1..50}; echo

# 开始转换
cat ~/2507/codes/cod_ids.txt | while read cod_id; do
    cif_file="$CIF_ROOT_DIR/${cod_id:0:1}/${cod_id:1:2}/${cod_id:3:2}/$cod_id.cif"
    pdb_file="$OUTPUT_DIR/$cod_id.pdb"
    sdf_dir="$OUTPUT_DIR/${cod_id:0:1}/${cod_id:1:2}/${cod_id:3:2}"
    sdf_file="$sdf_dir/$cod_id.sdf"

    mkdir -p "$sdf_dir"

    # 重复运行脚本用
    if [[ -f "$sdf_file" ]]; then
        echo "skip $sdf_file"
        echo "$cod_id" >> "$LOG_SKIPPED"
        continue
    fi

    # cif 转 pdb
    python "$CIF2PDB_SCRIPT" "$cif_file" "$pdb_file"

    # pdb 转 sdf
    if structconvert "$pdb_file" "$sdf_file" > /dev/null 2>&1; then
        echo "success: $cod_id"
        echo "$cod_id" >> "$LOG_SUCCESS"
        # 删除临时产生的 pdb 文件
        rm "$pdb_file"
        printf '%.0s=' {1..50}; echo
    else
        echo "failed: $cod_id"
        echo "$cod_id" >> "$LOG_FAILED"
        printf '%.0s=' {1..50}; echo
    fi

    python 

    # 防止内存累积，小间隔
    sleep 0.1
done

printf '%.0s=' {1..50}; echo
success_count=$(wc -l < "$LOG_SUCCESS")
failed_count=$(wc -l < "$LOG_FAILED")
skipped_count=$(wc -l < "$LOG_SKIPPED")
total_count=$(($success_count + $failed_count + $skipped_count))
echo "Convertion complete! "
echo "Success: $success_count, skipped: $skipped_count, failed: $failed_count, total: $total_count"
echo "Logs saved to $OUTPUT_DIR. "