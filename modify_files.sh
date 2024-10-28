#!/bin/bash

# 检查是否提供了目录地址参数
if [ -z "$1" ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# 检查提供的目录是否存在
if [ ! -d "$1" ]; then
    echo "Directory does not exist."
    exit 1
fi

# 切换到指定的目录
cd "$1" || exit

# 遍历当前目录下所有以.all.txt为后缀的文件
for file in *.all.txt; do
    # 提取文件名（不包括扩展名）
    sample_name=$(basename "$file" .all.txt)
    
    # 临时文件名
    temp_file=".temp.txt"

    # 使用sed命令替换第一行
    sed "1s/.*/gene_id\t$sample_name/" "$file" > "$temp_file"

    # 将临时文件替换原文件
    mv "$temp_file" "$file"
done
