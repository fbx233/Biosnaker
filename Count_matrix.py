import pandas as pd
import argparse
import os
from functools import reduce

def merge_files(input_files, output_file):
    data_list = []
    for file in input_files:
        try:
            # 先读取文件的前几行，打印出来以进行调试
            with open(file, 'r') as f:
                lines = [next(f) for _ in range(10)]
                print(f"Preview of {file}:\n{''.join(lines)}")

            # 跳过注释行并读取文件
            temp = pd.read_csv(file, sep='\t', comment='#', index_col=0)
            temp = temp.iloc[:, [-1]]  # 只保留最后一列
            temp.columns = [os.path.splitext(os.path.basename(file))[0]]  # 使用文件名作为列名
            data_list.append(temp)
        except Exception as e:
            print(f"Error processing file {file}: {e}")
            return

    if not data_list:
        print("No valid data to merge.")
        return

    # 使用reduce函数和join来合并所有DataFrame
    matrix_C = reduce(lambda left, right: left.join(right, how='outer'), data_list)

    # 确保目标目录存在
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # 输出到文件
    matrix_C.to_csv(output_file, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge featureCounts files.')
    parser.add_argument('--input_files', nargs='+', required=True, help='List of input files')
    parser.add_argument('--output_file', required=True, help='Output file path')
    args = parser.parse_args()

    merge_files(args.input_files, args.output_file)
