### Software required (on linux):

## Miniconda3 (add mirrors to speed up: https://mirror.tuna.tsinghua.edu.cn/help/anaconda/)

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

## Snakemake

conda install -n base -c conda-forge mamba
mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal==7.30.1 pandas scipy click

## sratoolkit (prefetch)/ascp for data downloading

# For example, ascp:

wget https://download.asperasoft.com/download/sw/connect/3.10.0/ibm-aspera-connect-3.10.0.180973-linux-g2.12-64.tar.gz --no-check-certificate
tar zxvf ibm-aspera-connect-3.10.0.180973-linux-g2.12-64.tar.gz
bash ibm-aspera-connect-3.10.0.180973-linux-g2.12-64.sh

## bedtools v2.31.0

wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static
mv bedtools.static.binary bedtools ; chmod a+x bedtools ; mv bedtools ~/bin

## STAR v2.7.10b

wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz
tar -xzf 2.7.10b.tar.gz
cd STAR-2.7.10b/bin/Linux_x86_64_static
ln -s $PWD/STAR ~/bin

## gatk v4.3.0.0

# Download gatk from the github releases page https://github.com/broadinstitute/gatk/releases

unzip gatk-4.3.0.0.zip
cd gatk-4.3.0.0
ln -s $PWD/gatk ~/bin

## fastp v0.23.2

wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
mv fastp ~/bin

## download GLIMPSE2 v2.0.0 (GLIMPSE2_*_static) from github release page https://github.com/odelaneau/glimpse/releases/

chmod a+x GLIMPSE2_*_static ; mv GLIMPSE2_*_static ~/bin

## stringtie v2.2.1

wget https://github.com/gpertea/stringtie/releases/download/v2.2.1/stringtie-2.2.1.Linux_x86_64.tar.gz
tar -xzf stringtie-2.2.1.Linux_x86_64.tar.gz
cd stringtie-2.2.1.Linux_x86_64
ln -s $PWD/stringtie ~/bin

## subread v2.0.5

wget https://udomain.dl.sourceforge.net/project/subread/subread-2.0.5/subread-2.0.5-Linux-x86_64.tar.gz --no-check-certificate
tar -xzf subread-2.0.5-Linux-x86_64.tar.gz
cd subread-2.0.5-Linux-x86_64/bin

# Add the software to the PATH in ~/.bashrc

export PATH=$PWD:$PATH

## bcftools, samtools and htslib v1.7+.

# According to http://www.htslib.org/download/

## Salmon v1.10.0

# Download salmon from github released page https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz

tar -xzf salmon-1.10.0_linux_x86_64.tar.gz
cd salmon-latest_linux_x86_64/bin
ln -s $PWD/salmon ~/bin

## gffread v0.12.7

# Download gffread from github released page https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz

tar -xzf gffread-0.12.7.Linux_x86_64.tar.gz
cd gffread-0.12.7.Linux_x86_64
ln -s $PWD/gffread ~/bin

## Download liftOver, gtfToGenePred, genePredToBed in https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64

wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
chmod 755 liftOver gtfToGenePred genePredToBed
mv liftOver gtfToGenePred genePredToBed ~/bin

## DaPars2 for 3'UTR polya detection

# Downdload DaPars2 codes from github https://github.com/3UTR/DaPars2

unzip DaPars2-main.zip

## WASP for unbiased allele-specific read mapping

# NOTE: It's better to change the output mode into 'wb' in line 139-142 in mapping/find_intersecting_snps.py to output BAM file.

git clone https://github.com/bmvdgeijn/WASP.git
conda create -n WASP hdf5 pysam scipy
conda activate WASP
pip install tables -i https://pypi.tuna.tsinghua.edu.cn/simple
cd WASP/snp2h5 ; make
conda deactivate
ln -s $PWD/snp2h5 ~/bin
ln -s $PWD/fasta2h5 ~/bin

## LeafCutter v0.2.7 for splicing detection

# installation from http://davidaknowles.github.io/leafcutter/articles/Installation.html for splicing detection.

git clone https://github.com/davidaknowles/leafcutter

# install RegTools (https://regtools.readthedocs.io/en/latest/#installation)

git clone https://github.com/griffithlab/regtools
cd regtools/
mkdir build
cd build/
cmake ..
make
ln -s $PWD/regtools ~/bin

# install 

## phASER v1.1.1 for ASE (https://github.com/secastel/phaser)

git clone https://github.com/secastel/phaser.git
conda create -n phaser python=2.7 pysam scipy cython
conda activate phaser
cd phaser/phaser/
python setup.py build_ext --inplace
conda deactivate
