# Analysis pipelines for the pilot phase (V0) of SheepGTEx
## 1. Introduction

As an essential component of the ongoing Farm animal Genotype-Tissue Expression (**FarmGTEx**, https://www.farmgtex.org/), the **SheepGTEx** (http://sheepgtex.farmgtex.org/) aims to establish a comprehensive public catalogue of genetic regulatory variants across diverse biological contexts (e.g., tissues, developmental stages, sex, environmental exposures, and genetic backgrounds) in sheep. This resource will serve as a foundational reference for understanding the genetic regulation of basic biological processes and complex traits. It will also support research in comparative transcriptomics, livestock evolution and domestication, precision breeding, sustainable agriculture, and human biomedicine.

![SheepGTEx-V0](https://raw.githubusercontent.com/FarmGTEx/SheepGTEx-Pipeline-v0/master/SheepGTEx-V0.jpg)

## 2. Analysis pipelines

This repository contains codes used in data analysis of the SheepGTEx pilot phase, including:

- 00.preparation: software and reference genome preparation.
- 01.WGS: whole genome sequencing raw data alignment, SNP calling and quality control.
- 02.RNA-seq: RNA-seq raw data processing, quantification and normalization of gene/exon/isoform/enhancer expression, splicing, RNA stability and 3'UTR APA, as well as ASE analysis.
- 03.metadata_prediction: sex, age and breed prediction for RNA-seq samples
- 04.gene_expression_analyses: gene expression analyses including variance decomposition, sample clustering, imputation, tissue specificity, co-expression, and sex- and developmental stage-biased expression.
- 05.molQTL_mapping: sample deduplication and molecular quantitative trait locus (molQTL) discovery, annotation and validation.
- 06.molQTL_colocalization: comparison between eQTL and other types of molQTL with the same gene.
- 07.molQTL_tissue_specificity: analyses of molQTL tissue sharing and effect size specificity. 
- 08.molQTL_context_specificity: context-specific (ancestry/sex/developmental stage) molQTL mapping.
- 09.GWAS_and_molQTL: GWAS meta-analysis, expression-mediated heritability, TWAS, and co-localization.
- 10.population_genetics: analyses of population structure (ADMIXTURE) and signatures of selection (FST, lnπ ratio, Tajima's D).
- 11.comparative_analysis: comparative analyses between sheep and cattle, pig and human, including ancestral allele inference, conservation of gene expression and eGene effect size.
- 12.ancient_process: ancient genome raw data alignment, imputation, quality control and gene expresison prediction.
- script: scripts used in the analyses above.

# 3. Citation

#### **A multi-tissue atlas of genetic regulatory effects in sheep**
 (Comming soon)
