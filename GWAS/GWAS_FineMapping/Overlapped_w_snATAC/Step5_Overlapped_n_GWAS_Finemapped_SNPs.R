################################
### NOTE : Date - 08/29/2023
################################
# Project: Pick's Disease + AD
# Ze Tristan SHI

# conda activate scRNAnATAC_R
#############################################
setwd('/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/Overlapped_wsnATACseq/')
dir_out <- '/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/Overlapped_wsnATACseq/'

# library(CMplot)
library(stringr)
library(tidyverse)
library(ggplot2)

##============ ##============
# NOTE Read in GWAS data
# gwas_data <- read.table("path/to/gwas_data.txt", header = TRUE)

AD_NG_GWAS_2022 <- read.table("/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Data/GWAS/GWAS_summary_statistics/AD_GWAS_Bellenguez_2022_NatureGenetics/GRCh38/GCST90027158/Data_Download/formatted/AD_NG2022_GWAS_Munged_12132022.txt", header = T)

head(AD_NG_GWAS_2022)
gwas_data_subset <- AD_NG_GWAS_2022 %>% select(SNP, CHR, BP, P)
head(gwas_data_subset)
colnames(gwas_data_subset) <- c('SNP', 'Chromosome', 'Position', 'P')

##============ ##============
# NOTE Read in FineMapped SNPs data
FineMapping_withSNPsetNumbers <- read.csv("/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/Figures/Figure3/FineMapping.csv")

# NOTE only kept AD_2022 finemapped SNPs
GWAS_finmapped <- FineMapping_withSNPsetNumbers %>% subset(AD_2022 > 0 | FTD_2014 > 0)
head(GWAS_finmapped)

# class(GWAS_finmapped$Credible_Set)

# NOTE reformat the columns to keep only a list of SNPs
# GWAS_finmapped$SNP_list <- gsub(" $", " ", GWAS_finmapped$Credible_Set)
GWAS_finmapped$SNP_list <- GWAS_finmapped$Credible_Set
GWAS_finmapped$SNP_list <- gsub(" ", ", ", GWAS_finmapped$SNP_list)
GWAS_finmapped$SNP_list <- gsub("\\s*,+\\s*", ", ", gsub("(^|\\s)[^rs]+", "", GWAS_finmapped$SNP_list))

# GWAS_finmapped$Credible_Set <- NULL

head(GWAS_finmapped)

GWAS_finmapped[GWAS_finmapped$GeneName=='TREM2', ]

# NOTE remove the last ", $" in rows in the column
# NOTE ", $" is a regular expression pattern that matches a comma followed by a space at the end of a string.
GWAS_finmapped$SNP_list <- gsub(", $", "", GWAS_finmapped$SNP_list)

class(GWAS_finmapped$SNP_list)
# [1] "character"

write.csv(GWAS_finmapped, file = paste0(dir_out, 'GWAS_finmapped_FTD2014_AD2022_wDup.csv'), quote=TRUE, row.names=FALSE)


##============ ##============
## NOTE Find duplicated rows based on the SNP_list column
duplicates <- GWAS_finmapped[duplicated(GWAS_finmapped$SNP_list) | duplicated(GWAS_finmapped$SNP_list, fromLast = TRUE), ]

# Print the duplicated rows
print(duplicates)

write.csv(duplicates, file = paste0(dir_out, 'GWAS_finmapped_duplicates_FTD2014_AD2022.csv'), quote=TRUE, row.names=FALSE)




# NOTE
#############################################
setwd('/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/Overlapped_wsnATACseq/')
dir_out <- '/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/Overlapped_wsnATACseq/'

# library(CMplot)
library(stringr)
library(tidyverse)
library(ggplot2)
library(dplyr)

##============ ##============
# NOTE Read in GWAS data
# gwas_data <- read.table("path/to/gwas_data.txt", header = TRUE)

GWAS_finmapped <- read.csv(paste0(dir_out, 'GWAS_finmapped_FTD2014_AD2022_wDup.csv'), header = T)

head(GWAS_finmapped, 1)

##============ ##============
## NOTE
list <- str_split(GWAS_finmapped$SNP_list, ", ")
length(list)
# To remove empty strings "" from a list that contains multiple lists inside,
# snp_lists <- lapply(list, function(x) x[nzchar(x)])
snp_vec <- unlist(list, recursive = TRUE)
snp_vec <- snp_vec[nzchar(snp_vec)]

## NOTE add Lead_SNPs  and get the unique
snp_vec_wleadSNP <- unique(c(unique(snp_vec), unique(unlist(str_split(GWAS_finmapped$Lead_SNPs, ", ")))))
snp_vec_wleadSNP
length(snp_vec_wleadSNP)
length(unique(snp_vec_wleadSNP))


# ##============ ##============
# ## NOTE
SNPs <- snp_vec_wleadSNP
length(SNPs)


# ##============ ##============
# ## NOTE  get location for all SNPs
library("biomaRt")
# snp_mart = useMart(biomart = "ENSEMBL_MART_SNP", dataset="hsapiens_snp", host='may2017.archive.ensembl.org')
snp_mart = useMart(biomart = "ENSEMBL_MART_SNP", dataset="hsapiens_snp")

# snp_ids = c("rs16828074", "rs17232800")
snp_ids <- SNPs
class(snp_ids)

snp_attributes = c("refsnp_id", "chr_name", "chrom_start")

snp_locations = getBM(attributes=snp_attributes, filters="snp_filter",
                       values=snp_ids, mart=snp_mart)


# List of values to exclude
values_to_exclude = c('HG2232_PATCH', 'HG708_PATCH', 'HSCHR14_3_CTG1', 'HSCHR17_1_CTG5', 'HSCHR17_2_CTG5', 'HSCHR19_4_CTG3_1', 'HSCHR19LRC_COX1_CTG3_1','HSCHR19LRC_PGF1_CTG3_1', 'HSCHR19LRC_PGF2_CTG3_1', 'HSCHR6_MHC_APD_CTG1','HSCHR6_MHC_COX_CTG1', 'HSCHR6_MHC_DBB_CTG1', 'HSCHR6_MHC_MANN_CTG1','HSCHR6_MHC_MCF_CTG1', 'HSCHR6_MHC_QBL_CTG1', 'HSCHR6_MHC_SSTO_CTG1')

# Filter the DataFrame to exclude values in the 'chr_name' column
SNP_dflocations <- snp_locations[!(snp_locations$chr_name %in% values_to_exclude), ]

print(SNP_dflocations)

length(unique(SNP_dflocations$refsnp_id))

identical(sort(SNP_dflocations$refsnp_id), sort(snp_ids))
# [1] TRUE

##============ ##============
# NOTE Read in GWAS data
# gwas_data <- read.table("path/to/gwas_data.txt", header = TRUE)

AD_NG_GWAS_2022 <- read.table("/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Data/GWAS/GWAS_summary_statistics/AD_GWAS_Bellenguez_2022_NatureGenetics/GRCh38/GCST90027158/Data_Download/formatted/AD_NG2022_GWAS_Munged_12132022.txt", header = T)

head(AD_NG_GWAS_2022)
gwas_data_subset <- AD_NG_GWAS_2022 %>% dplyr::select(SNP, CHR, BP, P)
# Error in (function (classes, fdef, mtable)  :
#   unable to find an inherited method for function ‘select’ for signature ‘"data.frame"’
head(gwas_data_subset)
colnames(gwas_data_subset) <- c('SNP', 'Chromosome', 'Position', 'P')

##============
# NOTE Check if SNP Position match
gwas_data_subset[gwas_data_subset$SNP=='rs17047661', ]
SNP_dflocations[SNP_dflocations$refsnp_id=='rs17047661', ]

##============
colnames(SNP_dflocations) <- c('SNP', 'Chromosome', 'Position')
head(SNP_dflocations)


# NOTE match the Credible_Set and GeneName
##============ ##============ ##============ ##============ ##============ ##============ ##============

# test <- SNP_dflocations

# Initialize SNP_dflocations$Credible_Set as a list
SNP_dflocations$Credible_Set <- vector("list", length(SNP_dflocations$SNP))
SNP_dflocations$GeneName <- vector("list", length(SNP_dflocations$SNP))

# NOTE library(stringr)
# Loop through each row in SNP_dflocations
for (i in 1:nrow(SNP_dflocations)) {
  # Split the SNP_list in GWAS_finmapped by comma and space
  split_snps <- str_split(GWAS_finmapped$SNP_list, ", ")

  # Find the rows in GWAS_finmapped$SNP_list where SNP_dflocations$SNP[i] is present
  match_rows <- sapply(split_snps, function(x) {
    any(SNP_dflocations$SNP[i] %in% unlist(x))
  })

  # Extract corresponding information from GWAS_finmapped$X for the matched rows
  matching_info <- GWAS_finmapped$X[match_rows]

  # Store the matching information in SNP_dflocations$Credible_Set and SNP_dflocations$GeneName as lists
  if (length(matching_info) > 0) {
    SNP_dflocations$Credible_Set[[i]] <- matching_info
    SNP_dflocations$GeneName[[i]] <- unique(GWAS_finmapped$GeneName[match_rows])
  } else {
    # When Credible_Set is NA, find unique GeneNames based on SNP matching Lead_SNPs
    lead_snp_match <- sapply(GWAS_finmapped$Lead_SNPs, function(x) {
      any(SNP_dflocations$SNP[i] %in% unlist(x))
    })
    if (any(lead_snp_match)) {
      SNP_dflocations$GeneName[[i]] <- unique(GWAS_finmapped$GeneName[lead_snp_match])
    } else {
      SNP_dflocations$Credible_Set[[i]] <- NA
      SNP_dflocations$GeneName[[i]] <- NA
    }
  }
}


head(SNP_dflocations)

# NOTE
# Flatten the list columns in SNP_dflocations
SNP_dflocations$Credible_Set <- sapply(SNP_dflocations$Credible_Set, function(x) if (is.null(x)) NA_character_ else paste(x, collapse = ", "))
SNP_dflocations$GeneName <- sapply(SNP_dflocations$GeneName, function(x) if (is.null(x)) NA_character_ else paste(x, collapse = ", "))


# Print the result
head(SNP_dflocations, 20)

SNP_dflocations[SNP_dflocations$SNP == 'rs9394764', ]

# Write the flattened data frame to a CSV file
write.csv(SNP_dflocations, file = paste0(dir_out, 'GWAS_finmapped_SNPs_withlocations_FTD2014_AD2022.csv'), quote = TRUE, row.names = FALSE)

# print(SNP_dflocations)




# NOTE
##============ ##============ ##============ ##============ ##============ ##============ ##============
#############################################

setwd('/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/Overlapped_wsnATACseq/')
dir_out <- '/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/Overlapped_wsnATACseq/'

library(dplyr)
library(data.table)
options(stringsAsFactors=F)
# library(CMplot)
library(stringr)
# library(tidyverse)
library(ggplot2)


SNP_dflocations <- read.csv(paste0(dir_out, 'GWAS_finmapped_SNPs_withlocations_FTD2014_AD2022.csv'))
# write.csv(SNP_dflocations, file = paste0(dir_out, 'GWAS_finmapped_SNPs_withlocations_FTD2014_AD2022.csv'), quote = TRUE, row.names = FALSE)
head(SNP_dflocations)
nrow(SNP_dflocations)

# Add characters to a numeric column in dataframe
SNP_dflocations$Chromosome <- sub("^", "chr", SNP_dflocations$Chromosome)

# SNP_dflocations$start <- SNP_dflocations$Position - 100
# SNP_dflocations$end <- SNP_dflocations$Position + 100
SNP_dflocations$start <- SNP_dflocations$Position - 150
SNP_dflocations$end <- SNP_dflocations$Position + 150
# SNP_dflocations$start <- SNP_dflocations$Position - 200
# SNP_dflocations$end <- SNP_dflocations$Position + 200
SNP_dflocations$range <- SNP_dflocations$end - SNP_dflocations$start


SNP_GRange_locations <- GenomicRanges::makeGRangesFromDataFrame(SNP_dflocations, seqnames.field = "Chromosome", start.field="start", end.field="end", keep.extra.columns=TRUE)
head(SNP_GRange_locations)

##============ ##============
# NOTE Use the FindMarkers Calculated Differential Open Chromatin regions to overlap with SNP regions

load('/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/All_Celltype_Markers_DARs/FindAllMarkers/DAR_allcelltypes_nogenename.rda')

head(DAR_allcelltypes_nogenename)

# Split the "X" column into "chr," "start," and "end" columns
DAR_allcelltypes_nogenename <- DAR_allcelltypes_nogenename %>%
  tidyr::separate(X, into = c("chr", "start", "end"), sep = "-")
#
head(DAR_allcelltypes_nogenename)

class(DAR_allcelltypes_nogenename$start)
# [1] "character"

# Convert "start" and "end" columns to numeric (if needed)
DAR_allcelltypes_nogenename$start <- as.numeric(DAR_allcelltypes_nogenename$start)
DAR_allcelltypes_nogenename$end <- as.numeric(DAR_allcelltypes_nogenename$end)


DAR_allcelltypes_nogenename_GRange <- GenomicRanges::makeGRangesFromDataFrame(DAR_allcelltypes_nogenename, seqnames.field = "chr", start.field="start", end.field="end", keep.extra.columns=TRUE)
head(DAR_allcelltypes_nogenename_GRange)


# PiD_DAR_GRange_hg19 <- XGR::xLiftOver(PiD_DAR_GRange_hg38, format.file = 'GRanges', build.conversion = "hg38.to.hg19")
# test <- which(overlapsAny(PiD_GRange_hg19, FTD_GRange_hg19))
snps.in.peak <- SNP_GRange_locations[which(IRanges::overlapsAny(SNP_GRange_locations, DAR_allcelltypes_nogenename_GRange))]
# length(snps.in.peak)
head(snps.in.peak)

peaks_wSNPs <- DAR_allcelltypes_nogenename_GRange[which(IRanges::overlapsAny(DAR_allcelltypes_nogenename_GRange, SNP_GRange_locations))]
head(peaks_wSNPs)


##============ ##============
# NOTE

load('/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/All_Celltype_Markers_DARs/Seurat_CallPeaks/peaks_incelltypes.rda')

head(peaks)

peaks$peak_range <- paste0(seqnames(peaks), "-", start(peaks), "-", end(peaks))

peaks_wSNPs <- peaks[which(IRanges::overlapsAny(peaks, SNP_GRange_locations))]
head(peaks_wSNPs)


snps.in.peak <- SNP_GRange_locations[which(IRanges::overlapsAny(SNP_GRange_locations, peaks))]
# length(snps.in.peak)
head(snps.in.peak)

snps.in.peak


# Find SNPs that overlap with peaks
snps.in.peak <- IRanges::subsetByOverlaps(SNP_GRange_locations, peaks)
# Use findOverlaps to identify overlaps
overlap_indices <- IRanges::findOverlaps(SNP_GRange_locations, peaks)

# Add the peak information to snps.in.peak
library(GenomicRanges)
snps.in.peak$peak_called_in <- unlist(peaks$peak_called_in[subjectHits(overlap_indices)])
snps.in.peak$peak_range <- unlist(peaks$peak_range[subjectHits(overlap_indices)])


# Convert the GRanges object to a DataFrame
snps_in_peak_df <- as.data.frame(snps.in.peak)


# Save the DataFrame as a CSV file
write.csv(snps_in_peak_df, file = paste0(dir_out, "snps_range300_in_peak.csv"), row.names = FALSE)
save(snps.in.peak, file= paste0(dir_out, "snps_range300_in_peak.rda"))
