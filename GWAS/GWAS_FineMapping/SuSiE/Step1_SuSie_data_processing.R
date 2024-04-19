## NOTE : Date - 11/22/2022
# Project: Pick's Disease + AD
# @Ze

# Coloc: a package for colocalisation analyses
# https://chr1swallace.github.io/coloc/articles/a01_intro.html
conda activate echolocator

###====================================
### NOTE : Date - 12/16/2022
###====================================
# NOTE the following loop is used for generating the shared SNPs between GWAS and cis-eQTL
# NOTE and also generating SNP-list for LD matrix

library(dplyr)
# library(coloc)
# library(locuscomparer)

####################################################################################
AD_GWAS_PATH <- "/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/SNPs_within_1MB_lead_SNP/AD_GWAS_Bellenguez_2022_NatureGenetics_LeadSNPs_wPositions/"
# QTL_path <- '/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Data/cis_eQTL/Summary_statistics_celltype_specific_cis_eQTLs_eight_brain_cell_types/'

# NOTE
# dir_out <- '/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/COLOC/AD/'
data_out <- '/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/COLOC_SuSiE/COLOC_SuSiE_Results/AD_NG2022/Data_susie/'
# analysis_out <- '/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/COLOC/locuscompare/FTD/'

####################################################################################
# NOTE I changed the Folder name to match the celltype cis-eQTL inside the folder
SNPid <- c('ABCA1_rs1800978',
'ABCA7_rs12151021',
'ABI3_rs616338',
'ACE_rs4277405',
'ADAM10_rs602602',
'ADAM17_rs72777026',
'ADAMTS1_rs2830489',
'ANK3_rs7068231',
'ANKH_rs112403360',
'APH1B_rs117618017',
'APP_rs2154481',
'BIN1_rs6733839',
'BLNK_rs6584063',
'CASS4_rs6014724',
'CD2AP_rs7767350',
'CLNK_rs6846529',
'CLU_rs11787077',
'COX7C_rs62374257',
'CR1_rs679515',
'CTSB_rs1065712',
'CTSH_rs12592898',
'DOC2A_rs1140239',
'EPDR1_rs6966331',
'EPHA1_rs11771145',
'FERMT2_rs17125924',
'FOXF1_rs16941239',
'GRN_rs5848',
'HLA_DQA1_rs6605556',
'HS3ST5_rs785129',
'ICA1_rs10952097',
'IDUA_rs3822030',
'IGH_rs7157106',
'IL34_rs4985556',
'INPP5D_rs10933431',
'JAZF1_rs1160871',
'KAT8_rs889555',
'KLF16_rs149080927',
'LILRB2_rs587709',
'MAF_rs450674',
'MAPT_rs199515',
'MME_rs61762319',
'MS4A4A_rs1582763',
'MYO15A_rs2242595',
'NCK2_rs143080277',
'PICALM_rs3851179',
'PLCG2_rs12446759',
'PLEKHA1_rs7908662',
'PRDM7_rs56407236',
'PRKD3_rs17020490',
'PTK2B_rs73223431',
'RASGEF1C_rs113706587',
'RBCK1_rs1358782',
'RHOH_rs2245466',
'SCIMP_rs7225151',
'SEC61G_rs76928645',
'SHARPIN_rs34173062',
'SIGLEC11_rs9304690',
'SLC24A4_rs7401792',
'SLC2A4RG_rs6742',
'SNX1_rs3848143',
'SORL1_rs74685827',
'SORT1_rs141749679',
'SPDYE3_rs7384878',
'SPI1_rs10437655',
'SPPL2A_rs8025980',
'TMEM106B_rs13237518',
'TNIP1_rs871269',
'TPCN1_rs6489896',
'TREM2_rs143332484',
'TSPAN14_rs6586028',
'TSPOAP1_rs2526377',
'UMAD1_rs6943429',
'UNC5CL_rs10947943',
'USP6NL_rs7912495',
'WDR12_rs139643391',
'WDR81_rs35048651')
# NOTE Testing:
# cur_SNPid <- 'APP_rs2154481'

##Run a loop for each cell-type
####################################################################################

for(i in 1:length(SNPid)){
  cur_SNPid=SNPid[i]
  print(paste("In Progress:", cur_SNPid))

  # NOTE loading the GWAS SNPs
  GWAS_SNPs <- read.table(paste0(AD_GWAS_PATH, cur_SNPid,'/', 'AD_NG2022_GWAS_', cur_SNPid,'_LDSNPs'), header = TRUE)
  # nrow(GWAS_SNPs)
  # head(GWAS_SNPs)
  GWAS_SNPs$chr_pos <- paste0(GWAS_SNPs$chr, ':', GWAS_SNPs$position)
  # GWAS_SNPs$chr_pos <- gsub('chr','',GWAS_SNPs$chr_pos)

  # NOTE extract the Lead SNPs
  # GWAS_SNPs_list <- GWAS_SNPs_filtered$chr_pos
  GWAS_SNPs_list <- GWAS_SNPs$chr_pos

  dir.create(paste0(data_out, cur_SNPid))

  # write.table(GWAS_SNPs_filtered, file= paste0(data_out, cur_SNPid,'/',"AD_NG2022_GWAS_",cur_SNPid,"_LDSNPs.txt"), col.names=TRUE, row.names=FALSE, quote = FALSE)
  write.table(GWAS_SNPs, file= paste0(data_out, cur_SNPid,'/',"AD_NG2022_GWAS_",cur_SNPid,"_LDSNPs.txt"), col.names=TRUE, row.names=FALSE, quote = FALSE)
  # write.table(cis_eQTL_min_subset, file= paste0(data_out, SNPid,'/',"cis_eQTL_",cur_celltype,'_',SNPid,".txt"), col.names=TRUE, row.names=FALSE, quote = FALSE)
  write.table(GWAS_SNPs_list, file= paste0(data_out, cur_SNPid,'/',"SNPslist_for_LD_",cur_SNPid,".txt"), col.names=FALSE, row.names=FALSE, quote = FALSE)

}

####################################################################################
