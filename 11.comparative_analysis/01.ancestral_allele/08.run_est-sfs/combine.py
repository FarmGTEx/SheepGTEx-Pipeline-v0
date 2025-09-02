for c in [str(i) for i in range(1,27)] + ['X', 'Y']:
    outfile_path = f'chr{c}_sfs_pvalue.txt'
    with open(outfile_path, 'w') as f_out:
        # 写入空头部 
        for _ in range(8):
            f_out.write('0\n')
        # 合并子文件 
        for sub in range(0,30): 
            sub_str = '{:02d}'.format(sub)
            subfile_path = f'chr{c}/chr{c}_sfs_pvalue_{sub_str}.txt'
            try:
                with open(subfile_path, 'r') as f_sub:
                    # 跳过头
                    for _ in range(8):
                        f_sub.readline()
                        # 将文件内容写入输出文件
                    for line in f_sub:
                        f_out.write(line)
            except FileNotFoundError:
                print(f'Warning: {subfile_path} not found, skipping...')
