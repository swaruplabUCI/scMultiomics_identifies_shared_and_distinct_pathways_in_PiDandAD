## NOTE : Date - 11/24/2022
# Project: Pick's Disease + AD
# @Ze

# Coloc: a package for colocalisation analyses
# https://chr1swallace.github.io/coloc/articles/a01_intro.html
# conda activate echolocator

###====================================
### NOTE : Date - 12/18/2022
###====================================
library(dplyr)
library(coloc)
# library(locuscomparer)

####################################################################################
# NOTE Loading PATH
# FTD_GWAS_PATH <- "/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/COLOC/COLOC_SuSiE_Results/FTD/Data_susie/"
AD_NG2022_GWAS_PATH <- "/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/COLOC_SuSiE/COLOC_SuSiE_Results/AD_NG2022/Data_susie/"

# NOTE Saving Dir
dir_out <- '/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/COLOC_SuSiE/COLOC_SuSiE_Results/AD_NG2022/Data_susie/'
analysis_out <- '/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/COLOC_SuSiE/COLOC_SuSiE_Results/AD_NG2022/output_susie/'

###====================================
# NOTE
###====================================
# NOTE GWAS SNPs

# celltypes <- c('Astrocytes', 'Endothelial.cells', 'Excitatory.neurons', 'Inhibitory.neurons', 'Microglia', 'Oligodendrocytes', 'OPCs...COPs', 'Pericytes')
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
# cur_SNPid <- c('APP_rs2154481')

for(i in 1:length(SNPid)){
  cur_SNPid=SNPid[i]
  print(paste("In Progress:", cur_SNPid))

  # NOTE constructing the PATH
  SNP_LD_PATH <- paste0(AD_NG2022_GWAS_PATH, cur_SNPid, '/LD_Matrix/')

  # if (cur_celltype == 'Excitatory.neurons'){
  #   SNPslist_celltype = 'Excitatory_neurons'} else if (cur_celltype == 'Inhibitory.neurons'){
  #   SNPslist_celltype = 'Inhibitory_neurons'} else if (cur_celltype == 'OPCs...COPs'){
  #   SNPslist_celltype = 'OPC'} else {SNPslist_celltype = cur_celltype}
  #
  # print(paste("Loading SNPslist In Progress:", SNPslist_celltype))
  # NOTE loading the GWAS SNPs list

  SNPslist <- read.table(paste0(SNP_LD_PATH, cur_SNPid, '_matrix.snplist'), header = FALSE)
  # head(SNPslist)
  nrow(SNPslist)

  print(paste("Loading GWAS_SNPs In Progress"))
  # NOTE loading the GWAS SNPs
  GWAS_SNPs <- read.table(paste0(AD_NG2022_GWAS_PATH, cur_SNPid,'/', 'AD_NG2022_GWAS_', cur_SNPid,'_LDSNPs.txt'), header = TRUE)
  #nrow(GWAS_SNPs)
  # GWAS_SNPs_subset <- GWAS_SNPs[GWAS_SNPs$snp %in% SNPslist$V1,]
  GWAS_SNPs_subset <- GWAS_SNPs[GWAS_SNPs$chr_pos %in% SNPslist$V1,]
  # nrow(GWAS_SNPs_subset)
  all.equal(GWAS_SNPs_subset$chr_pos, SNPslist$V1)
  identical(GWAS_SNPs_subset$chr_pos, SNPslist$V1)
  # NOTE add varbeta from se
  GWAS_SNPs_subset$varbeta=GWAS_SNPs_subset$se^2

  # NOTE loading LD matrix
  LD <- as.matrix(read.table(paste0(SNP_LD_PATH, cur_SNPid, '_matrix.ld')))

  colnames(LD) <- GWAS_SNPs_subset$snp
  rownames(LD) <- GWAS_SNPs_subset$snp
  LD[is.nan(LD)] <- 0

  ###====================================
  # NOTE Constructing the data matrix
  # NOTE
  GWAStype <- "cc"
  # eQTLtype <- "quant"
  GWAS_SNPs_data_list <- c(as.list(GWAS_SNPs_subset), LD = list(LD), type = GWAStype, N = 788989)
  # eQTL_SNPs_data_list <- c(as.list(eQTL_subset), LD = list(LD), sdY=1, type = eQTLtype, N = 192) # , N = 192 from the paper

  print(check_dataset(GWAS_SNPs_data_list))
  # print(check_dataset(eQTL_SNPs_data_list))
  # NOTE

  ###====================================
  # NOTE saving the data list of GWAS and eQTL for analysis
  # NOTE
  # save(GWAS_SNPs_data_list, eQTL_SNPs_data_list, file = paste0(dir_out, SNPid, '/', 'COLOC_data_', SNPid, '_', SNPslist_celltype, '.rda'))
  save(GWAS_SNPs_data_list, file = paste0(dir_out, cur_SNPid, '/', 'susie_data_AD_GWAS_2022_input_', cur_SNPid, '.rda'))

  ##====================================
  # NOTE Analysis
  print(paste0("SuSie Analysis between SNPs in GWAS: ", cur_SNPid))
  ###=============
  # NOTE SuSie
  print("coloc.susie: In Progress")
  S_GWAS=runsusie(GWAS_SNPs_data_list)
  print(summary(S_GWAS))

  susie_AD_summary <- summary(S_GWAS)
  # str(S_GWAS)

  susie_AD_snp <- S_GWAS$set$cs

  # NOTE saving
  dir.create(paste0(analysis_out, cur_SNPid))
  # # save(S_GWAS, S_QTL, susie.res, susie.res_result, my.res, my.res_result, file = paste0(analysis_out, SNPid, '/', 'COLOC_susie_analysis_', SNPid, '_', SNPslist_celltype, '.rda'))
  save(S_GWAS, susie_AD_snp, file = paste0(analysis_out, cur_SNPid, '/', 'susie_analysis_using_colocwarpper_', cur_SNPid,'.rda'))
  # NOTE capture.output save data list in csv or txt
  capture.output(susie_AD_summary, file= paste0(analysis_out, cur_SNPid,'/',"susie_results_summary_",cur_SNPid,".txt"))
  capture.output(susie_AD_snp, file= paste0(analysis_out, cur_SNPid,'/',"susie_results_snps_",cur_SNPid,".txt"))

  # NOTE remove to clean the data stored for the loop
  rm(S_GWAS, susie_AD_summary, susie_AD_snp)

}
