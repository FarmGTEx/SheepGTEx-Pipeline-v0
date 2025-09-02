library(MESuSiE)
library(data.table)
library(dplyr)
library(tidyr)
library(snpStats)
library(ggplot2)
library(cowplot)

ARGS <- commandArgs(trailingOnly = TRUE)
sumfile1 = ARGS[1] # summary file of the first population
sumfile2 = ARGS[2] # summary file of the second population
#NOTE: summary data required columns: CHR_POS, CHR, POS, REF, ALT, MAF, BETA, SE, Z and N.
plinkprefix1 = ARGS[3] # plink file prefix of the first population
plinkprefix2 = ARGS[4] # plink file prefix of the second population
outprefix = ARGS[5] # output file prefix
if (!file.exists(sumfile1)) { stop("Can not find the summary file1") }
if (!file.exists(sumfile2)) { stop("Can not find the summary file2") }
if (!file.exists(paste0(plinkprefix1, ".bed"))) { stop("Can not find the plink file1") }
if (!file.exists(paste0(plinkprefix2, ".bed"))) { stop("Can not find the plink file2") }

# https://borangao.github.io/meSuSie_Analysis/GWAS_QC.html
# 1. GWAS Preparations 
## MAF filtering
MAF_threshold=0.001
UKBB <- fread(sumfile1)
GLGC <- fread(sumfile2)
UKBB <- UKBB %>% filter(MAF > MAF_threshold, MAF < 1 - MAF_threshold)
GLGC <- GLGC %>% filter(MAF > MAF_threshold, MAF < 1 - MAF_threshold)

## Subset of SNPs
common_SNP <- find_common_snps(UKBB, GLGC)
UKBB_subset <- UKBB %>%
  filter(CHR_POS %in% common_SNP) %>%
  arrange(match(CHR_POS, common_SNP))
GLGC_subset <- GLGC %>%
  filter(CHR_POS %in% common_SNP) %>%
  arrange(match(CHR_POS, common_SNP))

## Identify Variation & Adjust BETA & Z Score
GLGC_subset_flip <- allele_flip(UKBB_subset, GLGC_subset)

# 2. Reference Panel Preparation
## Read in the bim file of both ancestries
WB_bim <- fread(paste0(plinkprefix1, ".bim"))
BB_bim <- fread(paste0(plinkprefix2, ".bim"))
WB_bim <- WB_bim %>%
  rename(CHR = V1, POS = V4, REF = V6, ALT = V5, CHR_POS=V2)
BB_bim <- BB_bim %>%
  rename(CHR = V1, POS = V4, REF = V6, ALT = V5, CHR_POS=V2)

## Find a common set of the SNPs across ancestries
common_SNP_bim <- find_common_snps(WB_bim, BB_bim)
WB_bim_subset <- WB_bim %>%
  filter(CHR_POS %in% common_SNP_bim)

## Find a common set of the SNPs with GWAS summary statistics
common_SNP_GWAS <- find_common_snps(UKBB_subset, WB_bim_subset)
UKBB_used <- UKBB_subset %>%
  filter(CHR_POS %in% common_SNP_GWAS) %>%
  arrange(match(CHR_POS, common_SNP_GWAS))
GLGC_used <- GLGC_subset_flip %>%
  filter(CHR_POS %in% common_SNP_GWAS) %>%
  arrange(match(CHR_POS, common_SNP_GWAS))

# 3. Check for LD Mismatches
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

## Generate correlation matrices for WB and BB genotype datasets
WB_cov <- create_correlation_matrix(plinkprefix1, common_SNP_GWAS)
BB_cov <- create_correlation_matrix(plinkprefix2, common_SNP_GWAS)
#WB_diagnostic <- kriging_rss(UKBB_used$Z, WB_cov)
#BB_diagnostic <- kriging_rss(GLGC_used$Z, BB_cov)
#WB_diagnostic$plot

# 4. Data Formatting for MESuSiE
summ_stat_list<-organize_gwas(UKBB_used%>%rename(SNP = CHR_POS),
                              GLGC_used%>%rename(SNP = CHR_POS),
                              c("EUR","CEA"))
#colnames(WB_cov)<-UKBB_used$CHR_POS
#colnames(BB_cov)<-UKBB_used$CHR_POS
LD_list<-organize_ld(WB_cov,BB_cov,summ_stat_list)

# 5. Run MESuSiE
MESuSiE_res<-meSuSie_core(LD_list,summ_stat_list,L=10)

# 6. save results
df = as.data.frame(cbind(row.names(WB_cov), NA, NA, MESuSiE_res$pip, MESuSiE_res$pip_config)) %>%
     rename(SNP = V1, cs = V2, category = V3, PIP = V4)
for (name in names(MESuSiE_res$cs$cs)) {
  snps = as.vector(row.names(WB_cov)[MESuSiE_res$cs$cs[[name]]])
  category = MESuSiE_res$cs$cs_category[[name]]
  df$cs[df$SNP %in% snps] = name
  df$category[df$SNP %in% snps] = category
}
fwrite(df, paste0(outprefix, ".txt"), col.names=T, sep="\t")

## plot results
# 提取 PIP 值
pip_1 <- get("pip", envir = MESuSiE_res)
pip_plot <- data.frame(Index = seq_along(pip_1), Value = pip_1)
pip_config_plot <- get("pip_config", envir = MESuSiE_res)

# 筛选 PIP > 0.5 的 SNP
filtered_indices <- which(pip_1 > 0.5)
filtered_values <- pip_1[filtered_indices]
result_pip <- data.frame(Index = filtered_indices, Value = filtered_values)

# 匹配 pip_config，判断变异类型
pip_config_filtered <- pip_config_plot[filtered_indices, , drop = FALSE]
max_values <- apply(pip_config_filtered, 1, max)
column_names <- apply(pip_config_filtered, 1, function(x) names(x)[which.max(x)])
classification <- ifelse(max_values > 0.5, column_names, "none")


classification_pip_result <- data.frame(
  Index = filtered_indices,
  MaxValue = max_values,
  ColumnName = column_names,
  Classification = classification
)

# 匹配SNP
common_SNP_index <- data.frame(
  Index = 1:nrow(UKBB_subset),
  CHR_POS = UKBB_subset$CHR_POS
)
merged_result <- merge(classification_pip_result, common_SNP_index, by = "Index", all.x = TRUE)
classification_pip_result <- merged_result
classification_pip_result

#作图——————————————————————————————————————————————————————————————————————————————————————————————————————————————————

#图表风格
custom_theme <- function() {
  theme(
    text = element_text( size = 12),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )
}


pip_plot <- function(res, R_mat, summary_data, classification_pip_result) {
  name_vec <- names(summary_data)
  
  if (!("POS" %in% colnames(summary_data[[name_vec[1]]]))) {
    summary_data <- lapply(summary_data, function(x) {
      x$POS = seq_len(nrow(x))
      return(x)
    })
  }
  
  Baseline_data <- summary_data[[names(summary_data)[1]]] %>% 
    select(SNP, POS, pos)
  
  Z_data <- Reduce(cbind, lapply(name_vec, function(x) {
    Z <- matrix(summary_data[[x]] %>% pull(Z), ncol = 1)
    colnames(Z) <- paste0("Z_", x)
    return(Z)
  }))
  
  PIP_data <- data.frame(PIP = res$pip)
  
  all_data <- bind_cols(Baseline_data, Z_data, PIP_data)
  
  lead_SNP_Ancestry <- all_data %>% select(SNP, starts_with("Z_")) %>% 
    pivot_longer(cols = -SNP, names_to = "ancestry", values_to = "z_val") %>% 
    arrange(desc(abs(z_val))) %>% slice(1) %>% select(SNP, ancestry)
  
  lead_SNP <- lead_SNP_Ancestry %>% pull(SNP)
  lead_Ancestry <- gsub("Z_", "", lead_SNP_Ancestry %>% pull(ancestry))
  
  r_data <- Reduce(cbind, lapply(name_vec, function(x) {
    r <- matrix(R_mat[[x]][, lead_SNP], ncol = 1)
    colnames(r) <- paste0("r_", x)
    return(r)
  }))
  
  all_data <- bind_cols(all_data, r_data)
  
  plotlist <- list()
  for (ancestry_name in name_vec) {
    gwas_plot_data <- all_data %>% select(SNP, pos, paste0("Z_", ancestry_name), paste0("r_", ancestry_name)) %>% 
      mutate(P = -log10(2 * pnorm(-abs(!!sym(paste0("Z_", ancestry_name))))))
    
    locus_zoom <- ggplot(gwas_plot_data, aes(x = pos, y = P, color = !!sym(paste0("r_", ancestry_name)))) + 
      geom_point(size = 1.5) +
      scale_color_stepsn(colors = c("navy", "lightskyblue", "green", "orange", "red"),
                         breaks = seq(0.2, 0.8, by = 0.2),
                         limits = c(0, 1),
                         show.limits = TRUE,
                         na.value = "grey50",
                         name = expression(R^2)) +
      custom_theme()+
      xlab(ancestry_name) +
      ylab("-log10 P-value") +
      geom_hline(yintercept = -log10(5e-08), linetype = "dashed") +
      custom_theme() +
      scale_x_continuous(labels = function(x) paste0(x / 1000000, "MB")) + # 修改 x 轴标签
      theme(legend.position = "none") +
      geom_text(aes(label = ifelse(SNP == lead_SNP, lead_SNP, "")),color="black", hjust = -0.1,size = 3)
    
    
    plotlist[[ancestry_name]] <- locus_zoom
  }
  
  PIP_plot_data <- all_data %>% 
    select(SNP, pos, PIP, paste0("r_", lead_Ancestry))
  
  # 将 PIP > 0.5 的 SNP 与 classification_pip_result 中的 CHR_POS 匹配，获取 Classification
  PIP_plot_data <- all_data %>% 
    select(SNP, pos, PIP, paste0("r_", lead_Ancestry)) %>% 
    left_join(classification_pip_result, by = c("SNP" = "CHR_POS")) %>% 
    mutate(Ancestry = ifelse(is.na(Classification), "none", Classification)) 
  
  # 绘制 PIP 图
  PIP_Plot <- ggplot(PIP_plot_data, aes(x = pos, y = PIP, color = !!sym(paste0("r_", lead_Ancestry)), shape = Ancestry)) +
    geom_point(size = 1.5,stroke=1) +
    scale_color_stepsn(colors = c("navy", "lightskyblue", "green", "orange", "red"),
                       breaks = seq(0.2, 0.8, by = 0.2),
                       limits = c(0, 1),
                       show.limits = TRUE,
                       na.value = "grey50",
                       name = expression(R^2)) + 
    scale_shape_manual(values = c(2, 6, 5, 16)) + 
    custom_theme()+
    theme(
      legend.key.size = unit(1, "lines"), 
      legend.justification = c(0, 0.5)
    )+
    xlab("MESuSiE") +
    ylab("PIP") +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    custom_theme() +
    theme(legend.position = "bottom") +
    scale_x_continuous(labels = function(x) paste0(x / 1000000, "MB")) + 
    geom_text(aes(label = ifelse(SNP == lead_SNP, lead_SNP, "")),color="black", hjust = -0.1,size = 3)
  plotlist[["PIP_Plot"]] <- PIP_Plot
  # 获取 PIP 图例
  PIP_plot_for_legend <- PIP_Plot + theme(legend.position = "bottom") 
  PIP_legend <- get_legend(PIP_plot_for_legend)
  
  combined_plot <- cowplot::plot_grid(plotlist = plotlist,
                                      ncol = 1,rel_heights = c(rep(1, length(plotlist) - 1), 1.5),align = "v")
  combined_plot_out <- plot_grid(combined_plot,
                                 PIP_legend,
                                 ncol = 1, rel_heights = c(1, 0.1), axis = "lr",align = "v")
  
  return(combined_plot_out)
}

cb_plot <- pip_plot(MESuSiE_res,LD_list,summ_stat_list,classification_pip_result)
ggsave(paste0(outprefix, ".png"), plot = cb_plot, width = 6, height = 12, units = "in", dpi = 300)

