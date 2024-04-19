################################
### NOTE : Date - 08/29/2023
################################
# Project: Pick's Disease + AD
# Ze Tristan SHI

# conda activate scRNAnATAC_R

# conda activate scRNAnATAC_R
#############################################

library(dplyr)
library(data.table)
options(stringsAsFactors=F)
# library(CMplot)
library(stringr)
# library(tidyverse)
library(ggplot2)
library(ggrepel)

# install.packages("MetBrewer")
library("MetBrewer")
# met.brewer("Hokusai1", n=100)

setwd('/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/Overlapped_wsnATACseq/')
dir_out <- '/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/Overlapped_wsnATACseq/'

#======== NOTE loading File gnomadv2
# /dfs7/swaruplab/zechuas/resources/gnomAD/gnomad.v2.1.1/gnomad.v2.1.1.lof_metrics.by_gene.txt
gnomadv2 <- read.delim("/dfs7/swaruplab/zechuas/resources/gnomAD/gnomad.v2.1.1/gnomad.v2.1.1.lof_metrics.by_gene.txt", header=TRUE)
head(gnomadv2)
names(gnomadv2)

length(gnomadv2$gene)
# [1] 19704
length(unique(gnomadv2$gene))
# [1] 19658

gene_pli_oelof <- gnomadv2 %>% dplyr::select(gene, pLI, oe_lof, oe_lof_upper)

#======== NOTE Load GWAS
# GWAS_finmapped_FTD2014_AD2022_wDup <- read.csv(paste0(dir_out, 'GWAS_finmapped_FTD2014_AD2022_wDup.csv'))
# head(GWAS_finmapped_FTD2014_AD2022_wDup, 1)

GWAS_finmapped <- read.csv('/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/Figures/Figure3/FineMapping.csv')

#======== NOTE Some key genes of interest to be labeled
# key.gns <- unique(GWAS_finmapped_FTD2014_AD2022_wDup$GeneName)
# key.gns <- key.gns[-41]

key.gns <- unique(GWAS_finmapped$GeneName)
key.gns <- key.gns[-56]
length(key.gns)
key.gns

# subet missing genename
subset_for_AD_GWASGenes<- subset(gene_pli_oelof, gene %in% key.gns)
nrow(subset_for_AD_GWASGenes)
gene_missing <- setdiff(key.gns, subset_for_AD_GWASGenes$gene)
gene_missing
# > gene_missing
# [1] "UMAD1"   "TSPOAP1"

key.gns <- key.gns[!key.gns %in% gene_missing]


#======== NOTE Plotting Color
# met.brewer("Hokusai1", n=10)
# pal2 = c("#007a8b", "#224b5e", "#94b594", "#edc775",  "#e09351", "#df7e66", "#b75347", "#6d2f20")
pal2 = c("#94b594", "#007a8b", "#95c36e", "#edc775",  "#e09351", "#df7e66", "#b75347", "#6d2f20")


#======== NOTE Plotting AD GWAS genes in all genes' pLI score
# Create the histogram with labels for key.gns using geom_text_repel
gg <- ggplot(gene_pli_oelof, aes(x = pLI, fill = cut(pLI, breaks = 8))) +
  geom_histogram(show.legend = FALSE) +
  scale_fill_manual(values = pal2) +  # Use custom colors
  theme_minimal() +
  labs(x = "pLI Score", y = "Gene Count") +
  ggtitle("Histogram of pLI Score with All GWAS Gene Count - labled with AD GWAS Genes") +
  geom_text_repel(data = subset(gene_pli_oelof, gene %in% key.gns),
                  aes(x = pLI, y = 0, label = gene),
                  box.padding = 0.5,  # Adjust as needed
                  size = 2,
                  max.overlaps = 100)



#
pdf(paste0(dir_out, 'Figures/', 'pLI_plot_AD_GWAS_genes_inAllGene_pLI_Sep05_2023.pdf'), width=20, height=10)
print(gg)
dev.off()



#======== NOTE Plotting AD GWAS genes in all genes' pLI score

subset_gene_pli_oelof <- subset(gene_pli_oelof, gene %in% key.gns)
length(key.gns)
nrow(subset_gene_pli_oelof)
# [1] 84

# subset_gene_pli_oelof$gene
# pal1 = c("#94b594", "#007a8b", "#95c36e", "#edc775",  "#e09351", "#df7e66", "#b75347", "#6d2f20")

#======== NOTE Plotting AD GWAS genes in all genes' pLI score
# Create the histogram with labels for key.gns using geom_text_repel
gg2 <- ggplot(subset_gene_pli_oelof, aes(x = pLI, fill = cut(pLI, breaks = 8))) +
  geom_histogram(show.legend = FALSE) +
  scale_fill_manual(values = pal2) +  # Use custom colors
  theme_minimal() +
  labs(x = "pLI Score", y = "Gene Count") +
  ggtitle("Histogram of pLI Score with AD GWAS Gene Count") +
  geom_text_repel(data = subset(subset_gene_pli_oelof, gene %in% key.gns),
                  aes(x = pLI, y = 0, label = gene),
                  box.padding = 0.5,  # Adjust as needed
                  size = 2,
                  max.overlaps = 100)



#
pdf(paste0(dir_out, 'Figures/', 'pLI_plot_AD_GWAS_genes_w_pLI_Sep05_2023.pdf'), width=20, height=10)
print(gg2)
dev.off()

#
pdf(paste0(dir_out, 'Figures/', 'pLI_plot_for_GeneLabels.pdf'), width=40, height=10)
print(gg2)
dev.off()



#############################################
# NOTE Prep for Heatmaps
#############################################


# Save the DataFrame as a CSV file
setwd('/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/Overlapped_wsnATACseq/')
dir_out <- '/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/Overlapped_wsnATACseq/'

snps_in_peak_df <- read.csv(paste0(dir_out, "snps_range300_in_peak.csv"))

head(snps_in_peak_df)


Credible_Set_list <- unique(unlist(snps_in_peak_df$Credible_Set))
Credible_Set_list

# Remove "NA" elements
Credible_Set_list <- Credible_Set_list[!is.na(Credible_Set_list)]

# Split elements containing commas and unlist the result
Credible_Set_list <- unlist(strsplit(Credible_Set_list, ", "))
Credible_Set_list

table(snps_in_peak_df$Credible_Set)
table(snps_in_peak_df$GeneName)


#======== NOTE Load GWAS
GWAS_finmapped_FTD2014_AD2022_wDup <- read.csv(paste0(dir_out, 'GWAS_finmapped_FTD2014_AD2022_wDup.csv'))
head(GWAS_finmapped_FTD2014_AD2022_wDup, 1)


overlapped_finemapped_CS <- subset(GWAS_finmapped_FTD2014_AD2022_wDup, X %in% Credible_Set_list)
nrow(overlapped_finemapped_CS)

head(overlapped_finemapped_CS[, c(1, 2, 3, 5, 6, 7)])
# > head(overlaped_finemapped_CS[, c(1, 2, 3, 5, 6, 7)])
#       X GeneName  Lead_SNPs  FTD_2014 AD_2021   AD_2022
# 3   CS3  SLC30A8  rs6984249 0.8100579       0 0.0000000
# 4   CS4     GLDN  rs2446406 0.7519578       0 0.0000000
# 5   CS5   SCUBE2 rs10840161 0.9192060       0 0.0000000
# 6   CS6    ZFHX3 rs12934137 0.6755328       0 0.0000000
# 10 CS12      CR1   rs679515 0.0000000       0 0.8237380
# 11 CS20   ADAM17 rs72777026 0.0000000       0 0.9209054

head(snps_in_peak_df)

# library(dplyr)
#======== NOTE
# Group and summarize snps_in_peak_df to get unique characters per Credible_Set
snps_summary <- snps_in_peak_df %>%
  group_by(Credible_Set) %>%
  summarize(cs_celltypes = paste(unique(peak_called_in), collapse = ","))

#
unique(snps_summary$Credible_Set)
# > unique(snps_summary$Credible_Set)
#  [1] "CS103"      "CS104"      "CS106"      "CS107"      "CS113"
#  [6] "CS114"      "CS12"       "CS121"      "CS122"      "CS127"
# [11] "CS130"      "CS131"      "CS132"      "CS133"      "CS137"
# [16] "CS138"      "CS142"      "CS143"      "CS144"      "CS145"
# [21] "CS146"      "CS151"      "CS155"      "CS157"      "CS158"
# [26] "CS165"      "CS167"      "CS168"      "CS172"      "CS190"
# [31] "CS196"      "CS20"       "CS200"      "CS204"      "CS205"
# [36] "CS21"       "CS3"        "CS33"       "CS34"       "CS35"
# [41] "CS36"       "CS4"        "CS43"       "CS49"       "CS5"
# [46] "CS51"       "CS53"       "CS58, CS63" "CS59, CS65" "CS6"
# [51] "CS60, CS64" "CS65"       "CS72"       "CS74, CS77" "CS75"
# [56] "CS76, CS79" "CS80"       "CS86"       "CS88"       "CS89"
# [61] NA

# Split the multiple Credible_Set values into separate rows
snps_summary <- snps_summary %>%
  tidyr::separate_rows(Credible_Set, sep = ", ") %>%
  filter(!is.na(Credible_Set)) %>%
  distinct()

#
unique(snps_summary$Credible_Set)

# NOTE check if they are the same
snps_summary[snps_summary$Credible_Set =="CS59", ]
snps_summary[snps_summary$Credible_Set =="CS65", ]

identical(snps_summary$cs_celltypes[snps_summary$Credible_Set =="CS58"], snps_summary$cs_celltypes[snps_summary$Credible_Set =="CS63"])
identical(snps_summary$cs_celltypes[snps_summary$Credible_Set =="CS59"], snps_summary$cs_celltypes[snps_summary$Credible_Set =="CS65"])

head(snps_summary)
nrow(snps_summary)


# NOTE Merge the summary back into overlaped_finemapped_CS based on Credible_Set
colnames(snps_summary)[colnames(snps_summary) == "Credible_Set"] <- "X"

overlapped_finemapped_CS <- overlapped_finemapped_CS %>%
  left_join(snps_summary, by = "X")

head(overlapped_finemapped_CS, 1)

class(overlapped_finemapped_CS$cs_celltypes)

# Split the cs_celltypes column, remove duplicates, and collapse unique values
overlapped_finemapped_CS$cs_celltypes <- sapply(strsplit(overlapped_finemapped_CS$cs_celltypes, ","), function(x) paste(unique(x), collapse = ","))

head(overlapped_finemapped_CS, 1)


write.csv(overlapped_finemapped_CS, file = paste0(dir_out, 'GWAS_finmapped_overlapped_finemapped_CS_FTD2014_AD2022.csv'), quote = TRUE, row.names = FALSE)


###====================================###
# NOTE # Project: Pick's Disease + AD
###====================================###

library(dplyr)
library(data.table)
options(stringsAsFactors=F)
# library(CMplot)
library(stringr)
# library(tidyverse)
library(ggplot2)
library(ggrepel)

setwd('/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/Overlapped_wsnATACseq/')
dir_out <- '/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/Overlapped_wsnATACseq/'

#======== NOTE loading File gnomadv2
# /dfs7/swaruplab/zechuas/resources/gnomAD/gnomad.v2.1.1/gnomad.v2.1.1.lof_metrics.by_gene.txt
gnomadv2 <- read.delim("/dfs7/swaruplab/zechuas/resources/gnomAD/gnomad.v2.1.1/gnomad.v2.1.1.lof_metrics.by_gene.txt", header=TRUE)
head(gnomadv2)
names(gnomadv2)

length(gnomadv2$gene)
# [1] 19704
length(unique(gnomadv2$gene))
# [1] 19658

gene_pli_oelof <- gnomadv2 %>% dplyr::select(gene, pLI, oe_lof, oe_lof_upper)

# NOTE loading more
# write.csv(overlapped_finemapped_CS, file = paste0(dir_out, 'GWAS_finmapped_overlapped_finemapped_CS_FTD2014_AD2022.csv'), quote = TRUE, row.names = FALSE)
overlapped_finemapped_CS <- read.csv(paste0(dir_out, 'GWAS_finmapped_overlapped_finemapped_CS_FTD2014_AD2022.csv'))
head(overlapped_finemapped_CS, 1)

names(overlapped_finemapped_CS)

keygenes <- unique(overlapped_finemapped_CS$GeneName)
keygenes
keygenes <- keygenes[-31]
length(keygenes)
# 47
keygenes

# subet missing genename
subset_for_AD_GWASGenes<- subset(gene_pli_oelof, gene %in% keygenes)
head(subset_for_AD_GWASGenes)
nrow(subset_for_AD_GWASGenes)
# 46 -- missing "UMAD1"

for(meta in names(subset_for_AD_GWASGenes)){
  overlapped_finemapped_CS[[meta]] <- subset_for_AD_GWASGenes[match(as.character(overlapped_finemapped_CS$GeneName), subset_for_AD_GWASGenes$gene), meta]
}


# NOTE reorder the cs_celltypes column
library(stringr)

# It seems that the sapply function is not appropriate for this task because it expects atomic vectors. Instead, you can use the apply function with MARGIN = 1 to apply your custom function to each row of the dataframe.
overlapped_finemapped_CS$cs_celltypes <- apply(overlapped_finemapped_CS, 1, function(row) {
  sorted_values <- str_sort(unlist(str_split(row["cs_celltypes"], ",")))
  return(paste(sorted_values, collapse = ","))
})

head(overlapped_finemapped_CS, 10)

overlapped_finemapped_CS[overlapped_finemapped_CS$GeneName == "UMAD1", ]

###====== NOTE saving
write.csv(overlapped_finemapped_CS, file = paste0(dir_out, 'GWAS_reordered_overlapped_finemapped_CS_FTD2014_AD2022_w_pLI_oe_lof.csv'), quote = TRUE, row.names = FALSE)




###====================================###
# NOTE # Project: Pick's Disease + AD
# Ze Tristan SHI
# conda activate scRNAnATAC_R
###====================================###
#==================
# NOTE HEATMAP
#==================
##============ ##============
# NOTE
# install.packages("NatParksPalettes")
# install.packages("MetBrewer")

# NOTE loading library
set.seed(5)

library(ComplexHeatmap)
library(viridis)
library(dplyr)
library(MetBrewer)
library(NatParksPalettes)

# NOTE PATH /dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/Overlapped_wsnATACseq/
data_in <- "/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/Overlapped_wsnATACseq//"
GWAS_fm_w_pLI_oe_lof <- read.csv(paste0(data_in, "GWAS_Update_reordered_overlapped_finemapped_CS_FTD2014_AD2022_w_pLI_oe_lof.csv"))

##== NOTE checking the dataset
names(GWAS_fm_w_pLI_oe_lof)

head(GWAS_fm_w_pLI_oe_lof)

length(unique(GWAS_fm_w_pLI_oe_lof$GeneName))
# [1] 48
GWAS_fm_w_pLI_oe_lof[GWAS_fm_w_pLI_oe_lof$GeneName =='IGH gene cluster', ]
# X         GeneName             Lead_SNPs
# 44 CS133 IGH gene cluster rs7157106, rs10131280
GWAS_fm_w_pLI_oe_lof[GWAS_fm_w_pLI_oe_lof$GeneName =='UMAD1', ]

GWAS_fm_w_pLI_oe_lof <- GWAS_fm_w_pLI_oe_lof[-44,]

##== NOTE checking the dataset
subset_finemapped <- GWAS_fm_w_pLI_oe_lof %>% select(X, GeneName, FTD_2014, AD_2022, pLI, oe_lof_upper)
#subset_finemapped[subset_finemapped$GeneName =='IGH gene cluster', ]
#subset_finemapped <- subset_finemapped[-44,]
#length(unique(subset_finemapped$GeneName))

head(subset_finemapped)
# names(subset_finemapped)

matrix_finemapped <- subset_finemapped
rownames(matrix_finemapped) <- matrix_finemapped$X
matrix_finemapped$GeneName <- NULL
matrix_finemapped$X <- NULL

# mat <- matrix_finemapped %>% select(FTD_2014, AD_2022)
# head(mat)

# NOTE Finamapped SuSiE
mat1 <- matrix_finemapped %>% select(FTD_2014, AD_2022) %>% data.matrix()
#mat1

# NOTE pLI
mat2 <- matrix_finemapped %>% select(pLI) %>% data.matrix()
#mat2

# NOTE LOEUF
mat3 <- matrix_finemapped %>% select(oe_lof_upper) %>% data.matrix()
#mat2


###==============
# NOTE:  ht1  --- Finamapped SuSiE
###==============
# row annotation
row_ha <- rowAnnotation(
  Gene_Name = anno_simple(as.character(subset_finemapped$GeneName), pch=as.character(subset_finemapped$GeneName),
                          simple_anno_size = unit(10, "mm"))
  #module = anno_simple(as.character(rownames(subset_finemapped)))
)

#ht1 = Heatmap(mat1, col=viridis(100), name = "Finamapped")
# rev(met.brewer("Hiroshige", n=100))
# met.brewer("Morgenstern", n=100)
ht1 <- Heatmap(mat1,col=viridis(100), name = "Finamapped",
               width = ncol(mat1)*unit(40, "mm"),
               # right_annotation=ha,column_names_rot = 45,
               left_annotation = row_ha,
               #row_names_side = 'left',
               #row_dend_side = 'right',
               #row_split = coloc_all_format_per_gene_direction[, 1],
               # row_title_rot = 0,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               layer_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.2f", pindex(mat1, i, j)),
                           x, y, gp = gpar(fontsize = 10,col="black"))
               },
               #row_title_gp = gpar(col = 'red', font = 2),
               #row_title = rep('',length(unique(coloc_all_format_per_gene_direction[, 1]))),
               heatmap_legend_param = list(title="r2",at = c(0,0.5,1),
                                           lables = c(0,0.5,1)))


#draw(ht1,heatmap_legend_side = "right")

###==============
# NOTE:  OPTIONAL --- make is.na disappear
###==============
# #Replace na values with 0 using is.na()
# mat2[is.na(mat2)] = 0
# mat3[is.na(mat3)] = 0



###==============
# NOTE:  ht2  --- pLI
###==============
ht2 = Heatmap(mat2, name = "pLI", col=met.brewer("VanGogh3", n=100), width = ncol(mat2)*unit(10, "mm"),
              heatmap_legend_param = list(title="pLI",at = c(0,0.5,1),
                                          lables = c(0,0.5,1)))
#draw(ht2,heatmap_legend_side = "right")

###==============
# NOTE:  Colors
###==============
# rev(met.brewer("Hiroshige", n=100))
# natparks.pals("Arches2", n=10)
# met.brewer("Homer1", n=100)


###==============
# NOTE:  ht3  --- LOEUF
###==============
ht3 = Heatmap(mat3, name = "LOEUF", col=rev(natparks.pals("Arches2", n=10)), width = ncol(mat2)*unit(10, "mm"),
              heatmap_legend_param = list(title="LOEUF",at = c(0,0.5, 1, 1.5, 2),
                                          lables = c(0,0.5,1, 1.5, 2)))
#draw(ht3,heatmap_legend_side = "right")

###==============
# NOTE:  ht4 -- GeneNames -- this is used to be deleted later but for split the rows
###==============

le <- as.character(subset_finemapped$GeneName)
ht4 = Heatmap(le, name = "Gene Name", width = 1*unit(10, "mm"))
# draw(ht4)


###==============
# NOTE:  ht5 -- overlapped with celltypes
###==============
# Load the required library for plotting
library(gplots)

# Define the color mapping
color_mapping <- data.frame(
  Category = c("ASC", "EX", "INH", "MG", "ODC", "OPC", "PER.END" , "grey"),
  Color = c("#F3746C", "#C1992D", "#51B348", "#29B891", "#11B4E9", "#355C7D", "#D56CAA" , "grey")
)

# # Define the color mapping
# color_mapping <- data.frame(
#   Category = c("ASC", "EX", "INH", "MG", "ODC", "OPC", "PER.END" ),
#   Color = c("#F3746C", "#C1992D", "#51B348", "#29B891", "#11B4E9", "#355C7D", "#D56CAA")
# )

# # NOTE: used a different category than color mapping since Gery should not be another column
categories <- c("ASC", "EX", "INH", "MG", "ODC", "OPC", "PER.END")

# Initialize the matrix with grey
num_samples <- nrow(GWAS_fm_w_pLI_oe_lof)
num_categories <- length(categories)
heatmap_matrix <- matrix("grey", nrow = num_samples, ncol = num_categories)

# Loop through each row and update the matrix
for (i in 1:num_samples) {
  cell_types <- unlist(strsplit(GWAS_fm_w_pLI_oe_lof$cs_celltypes[i], ","))
  for (category in cell_types) {
    if (category %in% categories) {
      column_index <- match(category, categories)
      heatmap_matrix[i, column_index] <- category  # Store category names, not colors
    }
  }
}

rownames(heatmap_matrix) <- GWAS_fm_w_pLI_oe_lof$X
colnames(heatmap_matrix) <- c("ASC", "EX", "INH", "MG", "ODC", "OPC", "PER.END")
head(heatmap_matrix)

# Create a numeric matching matrix
numeric_matrix <- matrix(0, nrow = nrow(heatmap_matrix), ncol = ncol(heatmap_matrix))
rownames(numeric_matrix) <- GWAS_fm_w_pLI_oe_lof$X
colnames(numeric_matrix) <- c("ASC", "EX", "INH", "MG", "ODC", "OPC", "PER.END")
head(numeric_matrix)


for (i in 1:nrow(heatmap_matrix)) {
  for (j in 1:ncol(heatmap_matrix)) {
    category <- heatmap_matrix[i, j]
    color <- color_mapping$Color[match(category, color_mapping$Category)]
    numeric_matrix[i, j] <- match(color, color_mapping$Color)
  }
}

# Create the heatmap using the numeric matrix
ht5 = Heatmap(numeric_matrix,col = color_mapping$Color,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              name = "Overlap", width = ncol(numeric_matrix)*unit(5, "mm"))

#ht5 <- heatmap(
#  numeric_matrix,
#  col = color_mapping$Color,
# Rowv = NA,
# Colv = NA,
# scale = "none"#, width = ncol(numeric_matrix)*unit(2, "mm")
##margins = c(5, 10),
##main = "Heatmap with Color Mapping"
#)

# ht5


###==============
# NOTE:  Plotting
###==============
ht_list = ht1 + ht2 + ht3 + ht5 + ht4
# draw(ht_list, cluster_rows = FALSE)
#draw(ht_list,row_split = le, cluster_rows = FALSE, ht_gap = unit(c(0.01, 0), "mm"))
# draw(ht_list,row_split = le, cluster_rows = FALSE)


pdf(paste0(data_in, 'Figures/', 'Heatmap_FTD2014_AD2022_SuSie_Finemapped_pLI_LOEUF_Sep12_2023.pdf'), width = 20, height=20)
# draw(ht2,heatmap_legend_side = "right")
draw(ht_list,row_split = le, cluster_rows = FALSE)
dev.off()




##============ ##============
# NOTE

#================ #================ #================ #================ #================ #================
# NOTE loading package
library(scCustomize)
library(Seurat)
library(Signac)
#================ NOTE loading dataset for processing
seurat_AD <- readRDS(file="/dfs7/swaruplab/smorabit/analysis/ADDS_2021/splitseq/integration/data/AD_integrated.rds" )
table(seurat_AD$cell_type)

# head(seurat_AD@meta.data)
# table(seurat_AD$Tissue)
# table(seurat_AD$Diagnosis)

table(seurat_AD$cell_type[seurat_AD$Diagnosis == "Control"])
 #  ASC   END    EX   FBR   INH    MG   ODC   OPC   PER   SMC
 # 5483   214 15794    77  6950  3269 35406  2730   185    91

seuratRNA <- subset(seurat_AD, Diagnosis == 'Control')
table(seuratRNA$Diagnosis)
table(seuratRNA$cell_type)

# NOTE run this one first
seuratRNA$cell_type <- gsub('PER', 'place_holder', seuratRNA$cell_type) # NOTE this stop the replace PER in all PER.END BUT I have to run this first


seuratRNA$cell_type <- gsub('END', 'PER.END', seuratRNA$cell_type)
seuratRNA$cell_type <- gsub('FBR', 'PER.END', seuratRNA$cell_type)
seuratRNA$cell_type <- gsub('SMC', 'PER.END', seuratRNA$cell_type)

seuratRNA$cell_type <- gsub('place_holder', 'PER.END', seuratRNA$cell_type)

table(seuratRNA$cell_type)

 seuratRNA$cell_type <- factor(
   as.character(seuratRNA$cell_type),
   levels = c("ASC", "EX", "INH", "MG", "ODC", "OPC", "PER.END") # An example
 )

table(Idents(seuratRNA))
# > table(Idents(seuratRNA))
#
# Mathys_2019 Swarup_2021   Zhou_2020
#       12153       21010       37036

Idents(seuratRNA) <- seuratRNA$cell_type


####################################
# NOTE rescale
# Normalize and FindVariableFeatures
seuratRNA <- NormalizeData(seuratRNA, normalization.method = "LogNormalize", scale.factor = 10000)

# NucSeq <- FindVariableFeatures(NucSeq, selection.method = "vst", nfeatures = 2000)
seuratRNA <- FindVariableFeatures(
  seuratRNA,
  selection.method = "vst",
  nfeatures = 4000
)
# scale data:
seuratRNA <- ScaleData(seuratRNA, features = rownames(seuratRNA))
# seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
#?ScaleData
####################################

AD_genenames <-  unique(GWAS_fm_w_pLI_oe_lof$GeneName)

# Sort the gene list in alphabetical order
sorted_gene_list <- sort(AD_genenames)

pdf(paste0(data_in, 'Figures/', 'DotPlot_scCustom_FTD2014_AD2022_SuSie_Finemapped_pLI_LOEUF_Sep12_2023.pdf'), width = 5, height=20)
DotPlot_scCustom(seurat_object = seuratRNA, features = rev(sorted_gene_list), flip_axes = T, x_lab_rotate = TRUE)
dev.off()
