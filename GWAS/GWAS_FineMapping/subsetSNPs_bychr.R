## NOTE : Date - 11/11/2022
# Project: Pick's Disease + AD
# @Ze

# https://github.com/boxiangliu/locuscomparer

###====================================
### NOTE : Date - 11/11/2022
###====================================
library(dplyr)

# NOTE PATH
FTD_PATH <- "/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Data/GWAS/GWAS_summary_statistics/International_Frontotemporal_Dementia_Genetics_Consortium/FTD_GWAS/SummaryData_Downloaded/formatted/"
AD_PATH <- "/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Data/GWAS/GWAS_summary_statistics/AD_GWAS_summary_statistics/GRCh37/GCST90012877/formatted/"

dir_out <- '/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/PAINTOR/SNPsbyChr_n_bedfiles/'

###====================================
# NOTE AD
###====================================
AD_GWAS <- read.table(paste0(AD_PATH, 'AD_GWAS_GCST90012877_munged.tsv'), header = TRUE)
head(AD_GWAS)


chr_number <- c(1:22)

# cur_chr = 7

for(i in 1:length(chr_number)){
  cur_chr=chr_number[i]
  print(paste("In Progress: chr -> ", cur_chr))

  AD_GWAS_chr <- subset(AD_GWAS, CHR == cur_chr)
  head(AD_GWAS_chr)
  table(AD_GWAS_chr$CHR)

  AD_GWAS_chr_SNPs <- AD_GWAS_chr %>% select(SNP)
  head(AD_GWAS_chr_SNPs)
  nrow(AD_GWAS_chr_SNPs)

  print(paste("Saving In Progress: chr -> ", cur_chr))
  write.table(AD_GWAS_chr_SNPs, file= paste0(dir_out, 'AD_SNPs_bychr/',"AD_GWAS_",cur_chr,"_SNPs.txt"), col.names=FALSE, row.names=FALSE, quote = FALSE)

}


###====================================
# NOTE FTD
###====================================
FTD_GWAS <- read.table(paste0(FTD_PATH, 'FTD_GWAS_META_Munged_11052022.txt'), header = TRUE)
head(FTD_GWAS)

chr_number <- c(1:22)

for(i in 1:length(chr_number)){
  cur_chr=chr_number[i]
  print(paste("In Progress: chr -> ", cur_chr))

  FTD_GWAS_chr <- subset(FTD_GWAS, CHR == cur_chr)
  head(FTD_GWAS_chr)
  table(FTD_GWAS_chr$CHR)

  FTD_GWAS_chr_SNPs <- FTD_GWAS_chr %>% select(SNP)
  head(FTD_GWAS_chr_SNPs)
  nrow(FTD_GWAS_chr_SNPs)

  print(paste("Saving In Progress: chr -> ", cur_chr))
  write.table(FTD_GWAS_chr_SNPs, file= paste0(dir_out, 'FTD_SNPs_bychr/',"FTD_GWAS_",cur_chr,"_SNPs.txt"), col.names=FALSE, row.names=FALSE, quote = FALSE)

}
