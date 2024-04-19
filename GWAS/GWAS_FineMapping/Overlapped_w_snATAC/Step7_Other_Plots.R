# NOTE GWAS finemapped SNPs' count
# Ze

# conda activate scRNAnATAC_R

library(ggplot2)
library(stringr)

# NOTE
# GWAS_df <- read.csv("~/Library/CloudStorage/GoogleDrive-zechuas@uci.edu/My Drive/Research/Project/AD_PiD_2021/Code/GWAS/output/FineMapped_SNPs/Overlapped_wsnATACseq/GWAS_Update_reordered_overlapped_finemapped_CS_FTD2014_AD2022_w_pLI_oe_lof.csv")

# loading PATH
data_dir <- '/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/GWAS_finemapping/Analysis/Overlapped_wsnATACseq/'
GWAS_df <- read.csv(paste0(data_dir, 'GWAS_Update_reordered_overlapped_finemapped_CS_FTD2014_AD2022_w_pLI_oe_lof.csv'))
head(GWAS_df)
class(GWAS_df$cs_celltypes)

#GWAS_df$cs_celltypes[1]
#length(GWAS_df$cs_celltypes[1])


# NOTE Function to count the items in a row
count_items <- function(row) {
  return(str_count(row, ",") + 1)  # Count commas and add 1 to get the number of items
}

# Apply the function to each row of the column
GWAS_df$ItemCount <- sapply(GWAS_df$cs_celltypes, count_items)


GWAS_df$GWAS_SNPs <- ifelse(GWAS_df$FTD_2014 > 0, "FTD GWAS Risk Loci", "AD GWAS Risk Loci")

# Display the result
# print(GWAS_df)

#
count_plot <- ggplot(GWAS_df,aes(x = ItemCount)) +
  geom_histogram(bins = 7,aes(fill = GWAS_SNPs)) +
  geom_histogram(bins = 7, fill = NA, color = 'black') +
  theme_minimal() +
  # scale_color_brewer(palette="Dark2")
  # scale_fill_manual(values=c("#E69F00", "#56B4E9"))
  scale_fill_manual(values=c("#a0522d", "#a572e8"))



pdf(paste0(data_dir, 'Figures/', 'Count_SNPs_forCelltype.pdf'), width=15, height=10)
count_plot
dev.off()
