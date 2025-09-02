for c in ['1', '2', '3','4','5','6','7','8','9','10','11', '12', '13','14','15','16','17','18','19','20','21', '22', '23','24','25','26','X','Y','MT','S']:
    for species in ['Cattle','Pig', 'Human']:
        with open('04.param/02.pairwise_maf2vcf/chr{:}_Sheep_{:}.vcf.param'.format(c, species), 'w') as f_out:
            f_out.write('input.file={:}/split/chr{:}.maf'.format(species,c) + '\n')
            f_out.write('input.file.compression=none' + '\n')
            f_out.write('output.log=log/04.chr{:}_Sheep_{:}.maffilter.vcf.log'.format(c, species) + '\n')
            f_out.write('maf.filter=VcfOutput(file=04.pairwise_align_vcf/chr{:}_Sheep_{:}.vcf.gz,\ '.format(c, species) + '\n')
            f_out.write('compression=gzip,\ ' + '\n')
            f_out.write('reference=Sheep,\ ' + '\n')
            f_out.write('genotypes=(Sheep,{:}),\ '.format(species) + '\n')
            f_out.write('all=yes)' + '\n')
