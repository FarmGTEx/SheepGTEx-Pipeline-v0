#!/bin/bash

# 定义文件名
input_file1="${1}_without_age.txt"        # 假设这是原始文件名
input_file2="${1}_without_predict_age.txt" # 预测结果文件
output_file="${1}_final.txt"               # 最终结果文件

# 提取每个文件的第一列
# 提取第一个文件的第一列，删除第一行
tail -n +2 "$input_file1" | cut -f1 > temp1.txt
cut -f1 "$input_file2" > temp2.txt

# 合并这两个临时文件的内容到最终文件
paste temp1.txt temp2.txt > "$output_file"

# 清理临时文件
rm temp1.txt temp2.txt

