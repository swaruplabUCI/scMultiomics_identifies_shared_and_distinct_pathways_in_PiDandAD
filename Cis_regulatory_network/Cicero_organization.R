##################################
# NOTE: Date: 04/19/2023 # Cicero Running
##################################
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(cicero)

# theme_set(theme_cowplot())
# Monocle3
# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                        'limma', 'S4Vectors', 'SingleCellExperiment',
#                        'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

# for Cicero
# BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer"))
# install.packages("devtools")
# devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")

# convertion of Seurat_Object
data_dir <- '/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Data/'

# seurat_obj <- readRDS(paste0(data_dir, "PiD_AD_integrated_input_12_16_2021.rds"))
seurat_obj_AD_only <- readRDS(paste0(data_dir, "Seurat_obj_AD_only_Dataset_wUpdatedFragPath_06_24_2022.rds"))
seurat_obj_PiD_only <- readRDS(paste0(data_dir, "Seurat_obj_PiD_only_Dataset_wUpdatedFragPath_06_24_2022.rds"))

##=========================================
### NOTE
hg38.genome.chrom.sizes <- seqlengths(seurat_obj_AD_only)
hg38.genome.chrom.sizes <- data.frame("chr" = names(hg38.genome.chrom.sizes), "length" = hg38.genome.chrom.sizes)
write.table(hg38.genome.chrom.sizes, file=paste0(data_dir, 'hg38.genome.chrom.sizes'), row.names = FALSE, sep="\t", quote = FALSE)

AD_cds <- as.cell_data_set(seurat_obj_AD_only)
PiD_cds <- as.cell_data_set(seurat_obj_PiD_only)

saveRDS(AD_cds, file = paste0(data_dir, "Monocle3_cds_ADnControl_04_19_2023.rds"))
saveRDS(PiD_cds, file = paste0(data_dir, "Monocle3_cds_PiDnControl_04_19_2023.rds"))
