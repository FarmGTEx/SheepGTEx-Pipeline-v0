sed 's/ /,/g' $1_transposed.txt > $1_temp_file.txt

# 第二步：使用sed修改第一行的第一个字符串
sed '1s/^[^,]*/sample/' $1_temp_file.txt > $1_transposed_comma.txt

# 清理临时文件
rm $1_temp_file.txt

awk -F"\t" 'NR==FNR{a[$1]=$3; next} {split($1, b, ","); if(b[1] in a) print $0","a[b[1]]; else print $0",NA"}' sample_age.list FS="," OFS="," $1_transposed_comma.txt > $2.with_age.txt

