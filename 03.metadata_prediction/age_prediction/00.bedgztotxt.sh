#使用方式：sh .sh 组织名

gunzip -k $1.expression.bed.gz 

cp $1.expression.bed $1.expression.txt

cut -f4- $1.expression.txt > $2.txt
#转置文件
awk '
{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' $2.txt > $2_transposed.txt


