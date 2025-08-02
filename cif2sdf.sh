#!/bin/bash

# 配置
JMOL_JAR="/home/clasper/softwares/jmol/jmol-16.3.31/Jmol.jar"
JAVA_BIN="/home/clasper/micromamba/envs/jmol/bin/java"
CIF_ROOT_DIR="/mnt/e/Users/clasper/tmp/cif"
OUTPUT_DIR="/home/clasper/2507/sdf"

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 日志
LOG_SUCCESS="$OUTPUT_DIR/success.log"
LOG_FAILED="$OUTPUT_DIR/failed.log"

eval "$(micromamba shell hook --shell bash)"
micromamba activate jmol

# 清空日志
> "$LOG_SUCCESS"
> "$LOG_FAILED"

# 准备转换
echo "converting CIF -> SDF"
echo "CIF dir: $CIF_ROOT_DIR"
echo "SDF dir: $OUTPUT_DIR"
printf '%.0s=' {1..50}; echo

# 开始转换
cat /mnt/e/Users/clasper/tmp/cif/cod_ids.txt | while read cod_id; do
    cif_file="$CIF_ROOT_DIR/${cod_id:0:1}/${cod_id:1:2}/${cod_id:3:2}/$cod_id.cif"
    sdf_file="$OUTPUT_DIR/$cod_id.sdf"
    jmol_file="$OUTPUT_DIR/$cod_id.jmol"

    # 生成 jmol 脚本
    echo 'set defaultLoadScript "unitcell"' > $jmol_file
    echo "load '"$cif_file"' {3 3 3}" >> $jmol_file
    echo "write V3000 '"$sdf_file"'" >> $jmol_file
    echo "quit" >> $jmol_file

    # 执行 jmol 脚本
    echo "$cif_file -> $sdf_file..."
    if "$JAVA_BIN" -jar "$JMOL_JAR" -n -s $jmol_file > /dev/null 2>&1; then
        echo "Success: $sdf_file"
        echo "$cif_file -> $sdf_file" >> "$LOG_SUCCESS"
    else
        echo "Failure: $sdf_file"
        echo "$cif_file" >> "$LOG_FAILED"
    fi

    # 删除 jmol 脚本
    rm "$jmol_file"

    # 防止内存累积，小间隔
    sleep 0.1
done

printf '%.0s=' {1..50}; echo
success_count=$(wc -l < "$LOG_SUCCESS")
failed_count=$(wc -l < "$LOG_FAILED")
total_count=$(($success_count + $failed_count))
echo "Convertion complete! "
echo "Converted $success_count/$total_count cif files. "
echo "Logs saved to $OUTPUT_DIR. "