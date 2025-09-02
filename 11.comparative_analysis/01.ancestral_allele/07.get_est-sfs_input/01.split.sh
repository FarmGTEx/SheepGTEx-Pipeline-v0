chrom=${1}
total_lines=$(wc -l < chr${chrom}_sfs_input.txt)  
lines_per_file=$(( (total_lines +29) /30 ))  
split -l "$lines_per_file" -d -a 2 chr${chrom}_sfs_input.txt split/chr${chrom}_sfs_input_ --additional-suffix=.txt
