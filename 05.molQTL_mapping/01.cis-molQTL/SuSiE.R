library(data.table)
library(dplyr)
library(tidyr)
library(snpStats)
library(susieR)
library(MESuSiE)
#library(Rfast)


ARGS <- commandArgs(trailingOnly = TRUE)
sumfile = ARGS[1] # summary file 
#NOTE: summary data required columns: CHR_POS, CHR, POS, REF, ALT, MAF, BETA, SE, Z and N.
plinkprefix = ARGS[2] # plink file prefix
outprefix = ARGS[3] # output file prefix
if (!file.exists(sumfile)) { stop("Can not find the summary file") }
if (!file.exists(paste0(plinkprefix, ".bed"))) { stop("Can not find the plink file") }

# read inputs
UKBB <- fread(sumfile)
WB_bim <- fread(paste0(plinkprefix, ".bim"))
WB_bim <- WB_bim %>%
  rename(CHR = V1, POS = V4, REF = V6, ALT = V5, CHR_POS=V2)
## MAF filtering
MAF_threshold=0.001
UKBB <- UKBB %>% filter(MAF > MAF_threshold, MAF < 1 - MAF_threshold)
common_SNP <- find_common_snps(UKBB, WB_bim)
UKBB_used <- UKBB %>%
    filter(CHR_POS %in% common_SNP) %>%
    arrange(match(CHR_POS, common_SNP))

# Check for LD Mismatches
create_correlation_matrix <- function(geno_name, SNPs_in_region) {
  ## Read in the genotype file selecting only the SNPs in the specified region
  geno_data <- read.plink(paste0(geno_name, ".bed"), select.snps = SNPs_in_region)
  
  ## Convert genotypes to numeric and match them with SNPs_in_region
  plink_geno <- as(geno_data$genotypes, "numeric")
  plink_geno <- plink_geno[, match(SNPs_in_region, geno_data$map$snp.name)]
  
  ## Replace NA values with the mean of the respective column (ignoring NA values)
  plink_geno <- apply(plink_geno, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
  })
  
  ## Calculate the correlation matrix
  cov <- cov2cor(crossprod(scale(plink_geno)))
  return(cov)
}

# Generate correlation matrices
WB_cov <- create_correlation_matrix(plinkprefix, common_SNP)

# Run susieR
susie_WB<-susie_rss(UKBB_used$Z, WB_cov, n=UKBB$N[1])
cat("susie_rss done!\n")

# Save results
df = as.data.frame(cbind(row.names(WB_cov), NA, susie_WB$pip)) %>%
  rename(SNP = V1, cs = V2, PIP = V3)
for (name in names(susie_WB$sets$cs)) {
  snps = as.vector(row.names(WB_cov)[susie_WB$sets$cs[[name]]])
  df$cs[df$SNP %in% snps] = name
}
fwrite(df, paste0(outprefix, ".txt"), col.names=T, sep="\t")
cat("susie results saved!\n")
