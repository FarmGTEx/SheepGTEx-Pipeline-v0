library(data.table)
library(dplyr)
library(purrr)
library(coloc)

ARGS <- commandArgs(trailingOnly = TRUE)
file1 = ARGS[1] # eQTL nominal file
file2 = ARGS[2] # other QTL nominal file
tissue = ARGS[3] # tissue name
chrom = ARGS[4]  # chromosome
outfile = ARGS[5] # output file path

# 数据读取
if (!file.exists(file1)) stop("找不到文件 file1")
if (!file.exists(file2)) stop("找不到文件 file2")

gwas <- fread(file1)
qtl  <- fread(file2)

# 检查列名（调试用）
cat("GWAS列名：\n"); print(colnames(gwas))
cat("QTL列名：\n"); print(colnames(qtl))

# 衍生 maf 列（使用 Freq 列）
gwas$maf <- ifelse(gwas$Freq <= 0.5, gwas$Freq, 1 - gwas$Freq)

# 估算 GWAS 样本量
gwas$N_gwas_est <- 1 / (gwas$se^2 * 2 * gwas$maf * (1 - gwas$maf))
sample_size_gwas <- round(median(gwas$N_gwas_est, na.rm = TRUE))
cat("估算的 GWAS 样本量:", sample_size_gwas, "\n")

# 合并 GWAS 和 QTL 数据
df <- inner_join(gwas, qtl, by = c("SNP" = "variant_id"))
head(df)
# 检查并清洗 QTL 中 ma_count 和 af 字段
cat("检查 QTL 中 ma_count 和 af 的缺失值...\n")
cat("NA in ma_count:", sum(is.na(df$ma_count)), "\n")
cat("NA in af:", sum(is.na(df$af)), "\n")

df_qtl_clean <- df[!is.na(df$ma_count) & !is.na(df$af) & df$af != 0, ]
head(df_qtl_clean)
if (nrow(df_qtl_clean) == 0) {
  cat("⚠️ 警告：清洗后 QTL 数据为空，无法估算样本量。\n")
  sample_size_qtl <- NA
} else {
  sample_size_qtl <- round(max(df_qtl_clean$ma_count / df_qtl_clean$af, na.rm = TRUE) / 2)
  cat("估算的 QTL 样本量:", sample_size_qtl, "\n")
}

# 运行 coloc 分析
unique_phenotypes <- unique(df$phenotype_id)

for (phenotype in unique_phenotypes) {
  df_target <- subset(df, phenotype_id == phenotype)
  
  # 检查缺失值
  if (any(is.na(df_target[, .(b, p, slope, pval_nominal)]))) next
  SNP = unique(df_target$SNP)
  df1_coloc <- list(
    beta = df_target$b,
    varbeta = df_target$se^2,
    snp = df_target$SNP,
    type = "quant",
    N = rep(sample_size_gwas, nrow(df_target)),
    MAF = ifelse(df_target$maf <= 0.5, df_target$maf, 1 - df_target$maf)
  )
  
  df2_coloc <- list(
    beta = df_target$slope,
    # pvalues = df_target$pval_nominal,
    varbeta = df_target$slope_se^2,
    snp = df_target$SNP,
    type = "quant",
    N = rep(sample_size_qtl, nrow(df_target)),
    MAF = ifelse(df_target$af <= 0.5, df_target$af, 1 - df_target$af)
  )
  
  cat("Running coloc for phenotype:", phenotype, "\n")
  my.res <- coloc.abf(dataset1 = df1_coloc, dataset2 = df2_coloc)
  
  lead_snp <- df_target$SNP[which.min(df_target$p)]  # GWAS 最显著 SNP
  
  results <- as.data.frame(t(unlist(my.res$summary)))
  # 添加 SNP 和 phenotype 信息到结果中
  results$SNP <- lead_snp  
  results$phenotype <- phenotype
  results$chrom <- chrom
  results$tissue <- tissue
  
  # 输出结果
  fwrite(results, outfile, append = TRUE, col.names = !file.exists(outfile), sep = "\t")
}
