###====================================
### NOTE : Date - 11/02/2022
###====================================
# NOTE Standardise the format of GWAS summary statistics
# Project: Pick's Disease + AD
# @Ze
# https://github.com/neurogenomics/MungeSumstats

# # NOTE library
# if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("MungeSumstats")
BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")


library(MungeSumstats)

# NOTE FTD_GWAS_META.txt
FTD_PATH <- "/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Data/GWAS/GWAS_summary_statistics/International_Frontotemporal_Dementia_Genetics_Consortium/FTD_GWAS/SummaryData_Downloaded/"
AD_PATH <- "/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Data/GWAS/GWAS_summary_statistics/AD_GWAS_summary_statistics/GRCh37/"

dir_out <- "/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Data/GWAS/GWAS_summary_statistics/TEST/"

###====================================
# NOTE FTD
###====================================
#
FTD_GWAS <- read.table(paste0(FTD_PATH, 'FTD_GWAS_META.txt'), header = TRUE)

head(FTD_GWAS)

FTD_GWAS_Munged <- MungeSumstats::format_sumstats(path = FTD_GWAS, ref_genome="GRCh37", save_path = paste0(FTD_PATH, 'formatted/', 'FTD_GWAS_META_Munged_11052022.txt'), bi_allelic_filter = FALSE, allele_flip_drop = FALSE)



###====================================
# NOTE AD
###====================================

AD_NG2022_PATH <- '/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Data/GWAS/GWAS_summary_statistics/AD_GWAS_Bellenguez_2022_NatureGenetics/GRCh38/GCST90027158/Data_Download/'

AD_NG2022_GWAS <- read.table(paste0(AD_NG2022_PATH, 'GCST90027158_buildGRCh38.tsv'), header = TRUE)
head(AD_NG2022_GWAS)

dir.create(paste0(AD_NG2022_PATH, 'formatted/'))

# library(MungeSumstats)
BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")

#
AD_NG2022_GWAS_Munged <- MungeSumstats::format_sumstats(path = AD_NG2022_GWAS, ref_genome="GRCh38", save_path = paste0(AD_NG2022_PATH, 'formatted/', 'AD_NG2022_GWAS_Munged_12072022.txt'), bi_allelic_filter = TRUE, allele_flip_drop = FALSE)
