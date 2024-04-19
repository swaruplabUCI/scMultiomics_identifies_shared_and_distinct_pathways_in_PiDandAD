###====================================
### NOTE : Date - 12/19/2022
###====================================
# Project: Pick's Disease + AD
# @Ze
# https://www.cog-genomics.org/plink/2.0/

# https://zenodo.org/record/3359882#.Y2FVdezMITv
# The 1000 Genomes Project Consortium. (2019). 1000 Genomes Project (1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.3359882

###====================================
# Working Directory
# PATH:
# /dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/PAINTOR


cd /dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/COLOC_SuSiE/COLOC_SuSiE_Results/AD_NG2022/Data_susie/


###===================
# NOTE SORT1
# Variant	Chr	Pos	Gene	Known_locus
# rs141749679	1	109345810	SORT1	New
# AD_NG2022_GWAS[AD_NG2022_GWAS$SNP == "rs141749679",]

# GeneName <- 'SORT1'
# SNPid <- 'rs141749679'
# SNPposition <- 109345810
# chrN <- 1

rsid=rs141749679
Genename_rsid=SORT1_rs141749679
chrN=1

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

###===================
# NOTE CR1
# Variant	Chr	Pos	Gene	Known_locus
# rs679515	1	207577223	CR1	CR1

# GeneName <- 'CR1'
# SNPid <- 'rs679515'
# SNPposition <- 207577223
# chrN <- 1

rsid=rs679515
Genename_rsid=CR1_rs679515
chrN=1

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

###===================
# NOTE ADAM17
# Variant	Chr	Pos	Gene	Known_locus
# rs72777026	2	9558882	ADAM17	New

# GeneName <- 'ADAM17'
# SNPid <- 'rs72777026'
# SNPposition <- 9558882
# chrN <- 2

# NOTE NEXT

rsid=rs72777026
Genename_rsid=ADAM17_rs72777026
chrN=2

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

###===================
# NOTE PRKD3
# Variant	Chr	Pos	Gene	Known_locus
# rs17020490	2	37304796	PRKD3	New

# GeneName <- 'PRKD3'
# SNPid <- 'rs17020490'
# SNPposition <- 37304796
# chrN <- 2

rsid=rs17020490
Genename_rsid=PRKD3_rs17020490
chrN=2

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

###===================
# NOTE NCK2
# Variant	Chr	Pos	Gene	Known_locus
# rs143080277	2	105749599	NCK2	New

# GeneName <- 'NCK2'
# SNPid <- 'rs143080277'
# SNPposition <- 105749599
# chrN <- 2

rsid=rs143080277
Genename_rsid=NCK2_rs143080277
chrN=2

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

###===================
# NOTE BIN1
# Variant	Chr	Pos	Gene	Known_locus
# rs6733839	2	127135234	BIN1	BIN1

# GeneName <- 'BIN1'
# SNPid <- 'rs6733839'
# SNPposition <- 127135234
# chrN <- 2

rsid=rs6733839
Genename_rsid=BIN1_rs6733839
chrN=2

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

###===================
# NOTE WDR12
# Variant	Chr	Pos	Gene	Known_locus
# rs139643391	2	202878716	WDR12	New

# GeneName <- 'WDR12'
# SNPid <- 'rs139643391'
# SNPposition <- 202878716
# chrN <- 2

rsid=rs139643391
Genename_rsid=WDR12_rs139643391
chrN=2

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

###===================
# NOTE INPP5D
# Variant	Chr	Pos	Gene	Known_locus
# rs10933431	2	233117202	INPP5D	INPP5D

# GeneName <- 'INPP5D'
# SNPid <- 'rs10933431'
# SNPposition <- 233117202
# chrN <- 2

rsid=rs10933431
Genename_rsid=INPP5D_rs10933431
chrN=2

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..


# ###===================
# # NOTE MME
# # Variant	Chr	Pos	Gene	Known_locus
# # rs16824536	3	155069722	MME	New
#
# GeneName <- 'MME'
# SNPid <- 'rs16824536'
# SNPposition <- 155069722
# chrN <- 3
#
# rsid=rs16824536
# Genename_rsid=MME_rs16824536
# chrN=3
#
# cd ${Genename_rsid}
# mkdir LD_Matrix/
#
# PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
# bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles
#
# ${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix
#
# cd ..

###===================
# NOTE MME
# Variant	Chr	Pos	Gene	Known_locus
# rs61762319	3	155084189	MME	New

# GeneName <- 'MME'
# SNPid <- 'rs61762319'
# SNPposition <- 155084189
# chrN <- 3

rsid=rs61762319
Genename_rsid=MME_rs61762319
chrN=3

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..


###===================
# NOTE IDUA
# Variant	Chr	Pos	Gene	Known_locus
# rs3822030	4	993555	IDUA	New

# GeneName <- 'IDUA'
# SNPid <- 'rs3822030'
# SNPposition <- 993555
# chrN <- 4

rsid=rs3822030
Genename_rsid=IDUA_rs3822030
chrN=4

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..


###===================
# NOTE CLNK
# Variant	Chr	Pos	Gene	Known_locus
# rs6846529	4	11023507	CLNK	CLNK/HS3ST1

# GeneName <- 'CLNK'
# SNPid <- 'rs6846529'
# SNPposition <- 11023507
# chrN <- 4

rsid=rs6846529
Genename_rsid=CLNK_rs6846529
chrN=4

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

###===================
# NOTE RHOH
# Variant	Chr	Pos	Gene	Known_locus
# rs2245466	4	40197226	RHOH	New

# GeneName <- 'RHOH'
# SNPid <- 'rs2245466'
# SNPposition <- 40197226
# chrN <- 4

rsid=rs2245466
Genename_rsid=RHOH_rs2245466
chrN=4

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..


###===================
# NOTE ANKH
# Variant	Chr	Pos	Gene	Known_locus
# rs112403360	5	14724304	ANKH	New

# GeneName <- 'ANKH'
# SNPid <- 'rs112403360'
# SNPposition <- 14724304
# chrN <- 5

rsid=rs112403360
Genename_rsid=ANKH_rs112403360
chrN=5

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE COX7C
# Variant	Chr	Pos	Gene	Known_locus
# rs62374257	5	86927378	COX7C	New

# GeneName <- 'COX7C'
# SNPid <- 'rs62374257'
# SNPposition <- 86927378
# chrN <- 5

rsid=rs62374257
Genename_rsid=COX7C_rs62374257
chrN=5

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE TNIP1
# Variant	Chr	Pos	Gene	Known_locus
# rs871269	5	151052827	TNIP1	New

# GeneName <- 'TNIP1'
# SNPid <- 'rs871269'
# SNPposition <- 151052827
# chrN <- 5

rsid=rs871269
Genename_rsid=TNIP1_rs871269
chrN=5

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE RASGEF1C
# Variant	Chr	Pos	Gene	Known_locus
# rs113706587	5	180201150	RASGEF1C	New

# GeneName <- 'RASGEF1C'
# SNPid <- 'rs113706587'
# SNPposition <- 180201150
# chrN <- 5

rsid=rs113706587
Genename_rsid=RASGEF1C_rs113706587
chrN=5

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE HLA_DQA1
# Variant	Chr	Pos	Gene	Known_locus
# rs6605556	6	32615322	HLA-DQA1	HLA

# GeneName <- 'HLA_DQA1'
# SNPid <- 'rs6605556'
# SNPposition <- 32615322
# chrN <- 6

rsid=rs6605556
Genename_rsid=HLA_DQA1_rs6605556
chrN=6

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE UNC5CL
# Variant	Chr	Pos	Gene	Known_locus
# rs10947943	6	41036354	UNC5CL	TREM2

# GeneName <- 'UNC5CL'
# SNPid <- 'rs10947943'
# SNPposition <- 41036354
# chrN <- 6

rsid=rs10947943
Genename_rsid=UNC5CL_rs10947943
chrN=6

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE TREM2
# Variant	Chr	Pos	Gene	Known_locus
# rs143332484	6	41161469	TREM2	TREM2
##===================
# NOTE TREM2 or TREML2
# Variant	Chr	Pos	Gene	Known_locus
# rs75932628	6	41161514	TREM2	TREM2
# rs60755019	6	41181270	TREML2	TREM2

# TREM2 are within in the 200 bp of the above lead SNPs under the same Genes, so let's use the same
##===================

# GeneName <- 'TREM2'
# SNPid <- 'rs143332484'
# SNPposition <- 41161469
# chrN <- 6

rsid=rs143332484
Genename_rsid=TREM2_rs143332484
chrN=6

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE CD2AP
# Variant	Chr	Pos	Gene	Known_locus
# rs7767350	6	47517390	CD2AP	CD2AP

# GeneName <- 'CD2AP'
# SNPid <- 'rs7767350'
# SNPposition <- 47517390
# chrN <- 6

rsid=rs7767350
Genename_rsid=CD2AP_rs7767350
chrN=6

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE HS3ST5
# Variant	Chr	Pos	Gene	Known_locus
# rs785129	6	114291731	HS3ST5	New

# GeneName <- 'HS3ST5'
# SNPid <- 'rs785129'
# SNPposition <- 114291731
# chrN <- 6

rsid=rs785129
Genename_rsid=HS3ST5_rs785129
chrN=6

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE UMAD1
# Variant	Chr	Pos	Gene	Known_locus
# rs6943429	7	7817263	UMAD1	New

# GeneName <- 'UMAD1'
# SNPid <- 'rs6943429'
# SNPposition <- 7817263
# chrN <- 7

rsid=rs6943429
Genename_rsid=UMAD1_rs6943429
chrN=7

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE ICA1
# Variant	Chr	Pos	Gene	Known_locus
# rs10952097	7	8204382	ICA1	New

# GeneName <- 'ICA1'
# SNPid <- 'rs10952097'
# SNPposition <- 8204382
# chrN <- 7

rsid=rs10952097
Genename_rsid=ICA1_rs10952097
chrN=7

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE TMEM106B
# Variant	Chr	Pos	Gene	Known_locus
# rs13237518	7	12229967	TMEM106B	New

# GeneName <- 'TMEM106B'
# SNPid <- 'rs13237518'
# SNPposition <- 12229967
# chrN <- 7

rsid=rs13237518
Genename_rsid=TMEM106B_rs13237518
chrN=7

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..


##===================
# NOTE JAZF1
# Variant	Chr	Pos	Gene	Known_locus
# rs1160871	7	28129126	JAZF1	New

# GeneName <- 'JAZF1'
# SNPid <- 'rs1160871'
# SNPposition <- 28129126
# chrN <- 7

rsid=rs1160871
Genename_rsid=JAZF1_rs1160871
chrN=7

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE EPDR1
# Variant	Chr	Pos	Gene	Known_locus
# rs6966331	7	37844191	EPDR1	NME8

# GeneName <- 'EPDR1'
# SNPid <- 'rs6966331'
# SNPposition <- 37844191
# chrN <- 7

rsid=rs6966331
Genename_rsid=EPDR1_rs6966331
chrN=7

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..


##===================
# NOTE SEC61G
# Variant	Chr	Pos	Gene	Known_locus
# rs76928645	7	54873635	SEC61G	New

# GeneName <- 'SEC61G'
# SNPid <- 'rs76928645'
# SNPposition <- 54873635
# chrN <- 7

rsid=rs76928645
Genename_rsid=SEC61G_rs76928645
chrN=7

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..


##===================
# NOTE SPDYE3
# Variant	Chr	Pos	Gene	Known_locus
# rs7384878	7	100334426	SPDYE3	ZCWPW1/NYAP1

# GeneName <- 'SPDYE3'
# SNPid <- 'rs7384878'
# SNPposition <- 100334426
# chrN <- 7

rsid=rs7384878
Genename_rsid=SPDYE3_rs7384878
chrN=7

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..


##===================
# NOTE EPHA1
# Variant	Chr	Pos	Gene	Known_locus
# rs11771145	7	143413669	EPHA1	EPHA1

# GeneName <- 'EPHA1'
# SNPid <- 'rs11771145'
# SNPposition <- 143413669
# chrN <- 7

rsid=rs11771145
Genename_rsid=EPHA1_rs11771145
chrN=7

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE CTSB
# Variant	Chr	Pos	Gene	Known_locus
# rs1065712	8	11844613	CTSB	New

# GeneName <- 'CTSB'
# SNPid <- 'rs1065712'
# SNPposition <- 11844613
# chrN <- 8

rsid=rs1065712
Genename_rsid=CTSB_rs1065712
chrN=8

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE PTK2B
# Variant	Chr	Pos	Gene	Known_locus
# rs73223431	8	27362470	PTK2B	PTK2B

# GeneName <- 'PTK2B'
# SNPid <- 'rs73223431'
# SNPposition <- 27362470
# chrN <- 8

rsid=rs73223431
Genename_rsid=PTK2B_rs73223431
chrN=8

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE CLU
# Variant	Chr	Pos	Gene	Known_locus
# rs11787077	8	27607795	CLU	CLU

# GeneName <- 'CLU'
# SNPid <- 'rs11787077'
# SNPposition <- 27607795
# chrN <- 8

rsid=rs11787077
Genename_rsid=CLU_rs11787077
chrN=8

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE SHARPIN
# Variant	Chr	Pos	Gene	Known_locus
# rs34173062	8	144103704	SHARPIN	New

# GeneName <- 'SHARPIN'
# SNPid <- 'rs34173062'
# SNPposition <- 144103704
# chrN <- 8

rsid=rs34173062
Genename_rsid=SHARPIN_rs34173062
chrN=8

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE ABCA1
# Variant	Chr	Pos	Gene	Known_locus
# rs1800978	9	104903697	ABCA1	New

# GeneName <- 'ABCA1'
# SNPid <- 'rs1800978'
# SNPposition <- 104903697
# chrN <- 9

rsid=rs1800978
Genename_rsid=ABCA1_rs1800978
chrN=9

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE USP6NL
# Variant	Chr	Pos	Gene	Known_locus
# rs7912495	10	11676714	USP6NL	ECHDC3

# GeneName <- 'USP6NL'
# SNPid <- 'rs7912495'
# SNPposition <- 11676714
# chrN <- 10

rsid=rs7912495
Genename_rsid=USP6NL_rs7912495
chrN=10

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE ANK3
# Variant	Chr	Pos	Gene	Known_locus
# rs7068231	10	60025170	ANK3	New

# GeneName <- 'ANK3'
# SNPid <- 'rs7068231'
# SNPposition <- 60025170
# chrN <- 10

rsid=rs7068231
Genename_rsid=ANK3_rs7068231
chrN=10

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE TSPAN14
# Variant	Chr	Pos	Gene	Known_locus
# rs6586028	10	80494228	TSPAN14	New

# GeneName <- 'TSPAN14'
# SNPid <- 'rs6586028'
# SNPposition <- 80494228
# chrN <- 10

rsid=rs6586028
Genename_rsid=TSPAN14_rs6586028
chrN=10

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE BLNK
# Variant	Chr	Pos	Gene	Known_locus
# rs6584063	10	96266650	BLNK	New

# GeneName <- 'BLNK'
# SNPid <- 'rs6584063'
# SNPposition <- 96266650
# chrN <- 10

rsid=rs6584063
Genename_rsid=BLNK_rs6584063
chrN=10

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE PLEKHA1
# Variant	Chr	Pos	Gene	Known_locus
# rs7908662	10	122413396	PLEKHA1	New

# GeneName <- 'PLEKHA1'
# SNPid <- 'rs7908662'
# SNPposition <- 122413396
# chrN <- 10

rsid=rs7908662
Genename_rsid=PLEKHA1_rs7908662
chrN=10

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE SPI1
# Variant	Chr	Pos	Gene	Known_locus
# rs10437655	11	47370397	SPI1	CELF1/SPI1

# GeneName <- 'SPI1'
# SNPid <- 'rs10437655'
# SNPposition <- 47370397
# chrN <- 11

rsid=rs10437655
Genename_rsid=SPI1_rs10437655
chrN=11

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..


##===================
# NOTE MS4A4A
# Variant	Chr	Pos	Gene	Known_locus
# rs1582763	11	60254475	MS4A4A	MS4A

# GeneName <- 'MS4A4A'
# SNPid <- 'rs1582763'
# SNPposition <- 60254475
# chrN <- 11

rsid=rs1582763
Genename_rsid=MS4A4A_rs1582763
chrN=11

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE PICALM
# Variant	Chr	Pos	Gene	Known_locus
# rs3851179	11	86157598	EED	PICALM

# GeneName <- 'PICALM'
# SNPid <- 'rs3851179'
# SNPposition <- 86157598
# chrN <- 11

rsid=rs3851179
Genename_rsid=PICALM_rs3851179
chrN=11

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE SORL1
# Variant	Chr	Pos	Gene	Known_locus
# rs74685827	11	121482368	SORL1	SORL1
#
# rs11218343	11	121564878	SORL1	SORL1

# GeneName <- 'SORL1'
# SNPid <- 'rs74685827'
# SNPposition <- 121482368
# chrN <- 11

rsid=rs74685827
Genename_rsid=SORL1_rs74685827
chrN=11

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE TPCN1
# Variant	Chr	Pos	Gene	Known_locus
# rs6489896	12	113281983	TPCN1	New

# GeneName <- 'TPCN1'
# SNPid <- 'rs6489896'
# SNPposition <- 113281983
# chrN <- 12

rsid=rs6489896
Genename_rsid=TPCN1_rs6489896
chrN=12

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..


##===================
# NOTE FERMT2
# Variant	Chr	Pos	Gene	Known_locus
# rs17125924	14	52924962	FERMT2	FERMT2

# GeneName <- 'FERMT2'
# SNPid <- 'rs17125924'
# SNPposition <- 52924962
# chrN <- 14

rsid=rs17125924
Genename_rsid=FERMT2_rs17125924
chrN=14

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE SLC24A4
# Variant	Chr	Pos	Gene	Known_locus
# rs7401792	14	92464917	SLC24A4	SLC24A4/RIN3
#
# rs12590654	14	92472511	SLC24A4	SLC24A4/RIN3

# GeneName <- 'SLC24A4'
# SNPid <- 'rs7401792'
# SNPposition <- 92464917
# chrN <- 14

rsid=rs7401792
Genename_rsid=SLC24A4_rs7401792
chrN=14

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE IGH
# Variant	Chr	Pos	Gene	Known_locus
# rs7157106	14	105761758	IGH gene cluster	New
#
# rs10131280	14	106665591	IGH gene cluster	New

# GeneName <- 'IGH'
# SNPid <- 'rs7157106'
# SNPposition <- 105761758
# chrN <- 14

rsid=rs7157106
Genename_rsid=IGH_rs7157106
chrN=14

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE SPPL2A
# Variant	Chr	Pos	Gene	Known_locus
# rs8025980	15	50701814	SPPL2A	SPPL2A

# GeneName <- 'SPPL2A'
# SNPid <- 'rs8025980'
# SNPposition <- 50701814
# chrN <- 15

rsid=rs8025980
Genename_rsid=SPPL2A_rs8025980
chrN=15

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..


##===================
# NOTE ADAM10
# Variant	Chr	Pos	Gene	Known_locus
# rs602602	15	58764824	MINDY2	ADAM10

# GeneName <- 'ADAM10'
# SNPid <- 'rs602602'
# SNPposition <- 58764824
# chrN <- 15

rsid=rs602602
Genename_rsid=ADAM10_rs602602
chrN=15

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE APH1B
# Variant	Chr	Pos	Gene	Known_locus
# rs117618017	15	63277703	APH1B	APH1B

# GeneName <- 'APH1B'
# SNPid <- 'rs117618017'
# SNPposition <- 63277703
# chrN <- 15

rsid=rs117618017
Genename_rsid=APH1B_rs117618017
chrN=15

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE SNX1
# Variant	Chr	Pos	Gene	Known_locus
# rs3848143	15	64131307	SNX1	New

# GeneName <- 'SNX1'
# SNPid <- 'rs3848143'
# SNPposition <- 64131307
# chrN <- 15

rsid=rs3848143
Genename_rsid=SNX1_rs3848143
chrN=15

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE CTSH
# Variant	Chr	Pos	Gene	Known_locus
# rs12592898	15	78936857	CTSH	New

# GeneName <- 'CTSH'
# SNPid <- 'rs12592898'
# SNPposition <- 78936857
# chrN <- 15

rsid=rs12592898
Genename_rsid=CTSH_rs12592898
chrN=15

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE DOC2A
# Variant	Chr	Pos	Gene	Known_locus
# rs1140239	16	30010081	DOC2A	New

# GeneName <- 'DOC2A'
# SNPid <- 'rs1140239'
# SNPposition <- 30010081
# chrN <- 16

rsid=rs1140239
Genename_rsid=DOC2A_rs1140239
chrN=16

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE KAT8
# Variant	Chr	Pos	Gene	Known_locus
# rs889555	16	31111250	BCKDK	KAT8

# GeneName <- 'KAT8'
# SNPid <- 'rs889555'
# SNPposition <- 31111250
# chrN <- 16

rsid=rs889555
Genename_rsid=KAT8_rs889555
chrN=16

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE IL34
# Variant	Chr	Pos	Gene	Known_locus
# rs4985556	16	70660097	IL34	IL34

# GeneName <- 'IL34'
# SNPid <- 'rs4985556'
# SNPposition <- 70660097
# chrN <- 16

rsid=rs4985556
Genename_rsid=IL34_rs4985556
chrN=16

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE MAF
# Variant	Chr	Pos	Gene	Known_locus
# rs450674	16	79574511	MAF	New

# GeneName <- 'MAF'
# SNPid <- 'rs450674'
# SNPposition <- 79574511
# chrN <- 16

rsid=rs450674
Genename_rsid=MAF_rs450674
chrN=16

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..


##===================
# NOTE PLCG2
# Variant	Chr	Pos	Gene	Known_locus
# rs12446759	16	81739398	PLCG2	PLCG2
#
# rs72824905	16	81908423	PLCG2	PLCG2

# GeneName <- 'PLCG2'
# SNPid <- 'rs12446759'
# SNPposition <- 81739398
# chrN <- 16

rsid=rs12446759
Genename_rsid=PLCG2_rs12446759
chrN=16

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE FOXF1
# Variant	Chr	Pos	Gene	Known_locus
# rs16941239	16	86420604	FOXF1	New

# GeneName <- 'FOXF1'
# SNPid <- 'rs16941239'
# SNPposition <- 86420604
# chrN <- 16

rsid=rs16941239
Genename_rsid=FOXF1_rs16941239
chrN=16

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE PRDM7
# Variant	Chr	Pos	Gene	Known_locus
# rs56407236	16	90103687	PRDM7	New

# GeneName <- 'PRDM7'
# SNPid <- 'rs56407236'
# SNPposition <- 90103687
# chrN <- 16

rsid=rs56407236
Genename_rsid=PRDM7_rs56407236
chrN=16

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE WDR81
# Variant	Chr	Pos	Gene	Known_locus
# rs35048651	17	1728046	WDR81	New

# GeneName <- 'WDR81'
# SNPid <- 'rs35048651'
# SNPposition <- 1728046
# chrN <- 17

rsid=rs35048651
Genename_rsid=WDR81_rs35048651
chrN=17

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE SCIMP
# Variant	Chr	Pos	Gene	Known_locus
# rs7225151	17	5233752	SCIMP	SCIMP/RABEP1

# GeneName <- 'SCIMP'
# SNPid <- 'rs7225151'
# SNPposition <- 5233752
# chrN <- 17

rsid=rs7225151
Genename_rsid=SCIMP_rs7225151
chrN=17

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE MYO15A
# Variant	Chr	Pos	Gene	Known_locus
# rs2242595	17	18156140	MYO15A	New

# GeneName <- 'MYO15A'
# SNPid <- 'rs2242595'
# SNPposition <- 18156140
# chrN <- 17

rsid=rs2242595
Genename_rsid=MYO15A_rs2242595
chrN=17

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE GRN
# Variant	Chr	Pos	Gene	Known_locus
# rs5848	17	44352876	GRN	New

# GeneName <- 'GRN'
# SNPid <- 'rs5848'
# SNPposition <- 44352876
# chrN <- 17

rsid=rs5848
Genename_rsid=GRN_rs5848
chrN=17

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE MAPT
# Variant	Chr	Pos	Gene	Known_locus
# rs199515	17	46779275	WNT3	MAPT

# GeneName <- 'MAPT'
# SNPid <- 'rs199515'
# SNPposition <- 46779275
# chrN <- 17

rsid=rs199515
Genename_rsid=MAPT_rs199515
chrN=17

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE ABI3
# Variant	Chr	Pos	Gene	Known_locus
# rs616338	17	49219935	ABI3	ABI3

# GeneName <- 'ABI3'
# SNPid <- 'rs616338'
# SNPposition <- 49219935
# chrN <- 17

rsid=rs616338
Genename_rsid=ABI3_rs616338
chrN=17

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE TSPOAP1
# Variant	Chr	Pos	Gene	Known_locus
# rs2526377	17	58332680	TSPOAP1	TSPOAP1

# GeneName <- 'TSPOAP1'
# SNPid <- 'rs2526377'
# SNPposition <- 58332680
# chrN <- 17

rsid=rs2526377
Genename_rsid=TSPOAP1_rs2526377
chrN=17

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE ACE
# Variant	Chr	Pos	Gene	Known_locus
# rs4277405	17	63471557	ACE	ACE

# GeneName <- 'ACE'
# SNPid <- 'rs4277405'
# SNPposition <- 63471557
# chrN <- 17

rsid=rs4277405
Genename_rsid=ACE_rs4277405
chrN=17

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE ABCA7
# Variant	Chr	Pos	Gene	Known_locus
# rs12151021	19	1050875	ABCA7	ABCA7

# GeneName <- 'ABCA7'
# SNPid <- 'rs12151021'
# SNPposition <- 1050875
# chrN <- 19

rsid=rs12151021
Genename_rsid=ABCA7_rs12151021
chrN=19

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE KLF16
# Variant	Chr	Pos	Gene	Known_locus
# rs149080927	19	1854254	KLF16	New

# GeneName <- 'KLF16'
# SNPid <- 'rs149080927'
# SNPposition <- 1854254
# chrN <- 19

rsid=rs149080927
Genename_rsid=KLF16_rs149080927
chrN=19

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE SIGLEC11
# Variant	Chr	Pos	Gene	Known_locus
# rs9304690	19	49950060	SIGLEC11	New

# GeneName <- 'SIGLEC11'
# SNPid <- 'rs9304690'
# SNPposition <- 49950060
# chrN <- 19

rsid=rs9304690
Genename_rsid=SIGLEC11_rs9304690
chrN=19

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE LILRB2
# Variant	Chr	Pos	Gene	Known_locus
# rs587709	19	54267597	LILRB2	New

# GeneName <- 'LILRB2'
# SNPid <- 'rs587709'
# SNPposition <- 54267597
# chrN <- 19

rsid=rs587709
Genename_rsid=LILRB2_rs587709
chrN=19

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE RBCK1
# Variant	Chr	Pos	Gene	Known_locus
# rs1358782	20	413334	RBCK1	New

# GeneName <- 'RBCK1'
# SNPid <- 'rs1358782'
# SNPposition <- 413334
# chrN <- 20

rsid=rs1358782
Genename_rsid=RBCK1_rs1358782
chrN=20

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE CASS4
# Variant	Chr	Pos	Gene	Known_locus
# rs6014724	20	56423488	CASS4	CASS4

# GeneName <- 'CASS4'
# SNPid <- 'rs6014724'
# SNPposition <- 56423488
# chrN <- 20

rsid=rs6014724
Genename_rsid=CASS4_rs6014724
chrN=20

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE SLC2A4RG
# Variant	Chr	Pos	Gene	Known_locus
# rs6742	20	63743088	SLC2A4RG	New

# GeneName <- 'SLC2A4RG'
# SNPid <- 'rs6742'
# SNPposition <- 63743088
# chrN <- 20

rsid=rs6742
Genename_rsid=SLC2A4RG_rs6742
chrN=20

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE APP
# Variant	Chr	Pos	Gene	Known_locus
# rs2154481	21	26101558	APP	New

# GeneName <- 'APP'
# SNPid <- 'rs2154481'
# SNPposition <- 26101558
# chrN <- 21

rsid=rs2154481
Genename_rsid=APP_rs2154481
chrN=21

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..

##===================
# NOTE ADAMTS1
# Variant	Chr	Pos	Gene	Known_locus
# rs2830489	21	26775872	ADAMTS1	ADAMTS1

# GeneName <- 'ADAMTS1'
# SNPid <- 'rs2830489'
# SNPposition <- 26775872
# chrN <- 21

rsid=rs2830489
Genename_rsid=ADAMTS1_rs2830489
chrN=21

cd ${Genename_rsid}
mkdir LD_Matrix/

PLINK=/dfs7/swaruplab/zechuas/tools_software/PLINK
bedfile=/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/SNPsbyChr_n_bedfiles/AD_NG2022_SNPs_bedfiles

${PLINK}/plink --bfile ${bedfile}/AD_NG2022_GWAS.chr${chrN}_20181129.GRCh38 --write-snplist --extract SNPslist_for_LD_${Genename_rsid}.txt --keep-allele-order --r square --out LD_Matrix/${Genename_rsid}_matrix

cd ..
