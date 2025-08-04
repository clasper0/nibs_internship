#!/bin/bash

# 配置
JMOL_JAR="/home/clasper/softwares/jmol/jmol-16.3.31/Jmol.jar"
JAVA_BIN="/home/clasper/micromamba/envs/jmol/bin/java"
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

eval "$(micromamba shell hook --shell bash)"
micromamba activate jmol

# 准备转换
echo "converting CIF -> SDF"
echo "CIF dir: $CIF_ROOT_DIR"
echo "SDF dir: $OUTPUT_DIR"
printf '%.0s=' {1..50}; echo

# 开始转换
cat /mnt/e/Users/clasper/tmp/cif/cod_ids.txt | while read cod_id; do
    cif_file="$CIF_ROOT_DIR/${cod_id:0:1}/${cod_id:1:2}/${cod_id:3:2}/$cod_id.cif"
    sdf_dir="$OUTPUT_DIR/${cod_id:0:1}/${cod_id:1:2}/${cod_id:3:2}"
    old_sdf_file="$OUTPUT_DIR/$cod_id.sdf"
    new_sdf_file="$sdf_dir/$cod_id.sdf"
    jmol_file="$OUTPUT_DIR/$cod_id.jmol"

    mkdir -p "$sdf_dir"

    if [[ -f "$old_sdf_file" ]]; then
        echo "move $old_sdf_file"
        mv "$old_sdf_file" "$new_sdf_file"
        echo "$cod_id" >> "$LOG_SKIPPED"
        continue
    fi

    # 生成 jmol 脚本
    echo 'set defaultLoadScript "unitcell"' > $jmol_file
    echo "load '"$cif_file"' {3 3 3}" >> $jmol_file
    echo "write V3000 '"$new_sdf_file"'" >> $jmol_file
    echo "quit" >> $jmol_file

    # 执行 jmol 脚本
    echo "$cif_file -> $new_sdf_file..."
    if timeout 30 "$JAVA_BIN" -jar "$JMOL_JAR" -n -s $jmol_file > /dev/null 2>&1; then
        echo "Success: $cod_id"
        echo "$cod_id" >> "$LOG_SUCCESS"
    else
        # 如果超时自动跳过当前
        if [[ $? -eq 124 ]]; then 
            echo "timeout: $cod_id"
            echo "$cod_id" >> "$LOG_FAILED"
        else
            echo "Failure: $cod_id"
            echo "$cod_id" >> "$LOG_FAILED"
        fi

        # 删除可能的不完整文件
        rm -f "$new_sdf_file"
    fi

    # 删除 jmol 脚本
    rm "$jmol_file"

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