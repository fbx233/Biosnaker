#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
source /public3/soft/modules/module.sh
module load anaconda/3-Python3.7.7-wjl
source activate py3.9.10

# --input_path 输入数据所在地址。
# --output_path 输出结果文件所在地址。
# --suffix 设置所需整合文件的后缀。

python -u Count_matrix.py \
--input_path "/public3/home/scg8403/Data/zhy/result/read_counts/result/" \
--output_path "/public3/home/scg8403/Data/zhy/result/counts_matrix" \
--suffix ".all.txt"
