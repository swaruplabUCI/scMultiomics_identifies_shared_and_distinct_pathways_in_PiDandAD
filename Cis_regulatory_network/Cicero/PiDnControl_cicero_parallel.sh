#!/bin/bash
#SBATCH --job-name=PiD_Cicero
#SBATCH -p highmem
#SBATCH -A vswarup_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --error=error_%A_%a.txt ## error log file name: %A is job id, %a is array task id
#SBATCH --output=out_%A_%a.txt  ## output filename
#SBATCH --mem 256G
#SBATCH --array=1-7
#SBATCH --time=72:00:00

source /data/homezvol0/zechuas/.bashrc
conda activate Monocle3


SEURAT="/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Data/Monocle3_cds_PiDnControl_04_19_2023.rds"
OUTDIR="/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/Cis_regulatory_network/Processing_window3e5_sampleN200/PiDnControl/"
CLUST="celltype"
GENOME="/dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Data/hg38.genome.chrom.sizes"

Rscript --vanilla /dfs7/swaruplab/zechuas/Projects/PiD_2021/scATAC_data/Analysis/Cis_regulatory_network/Processing_window3e5_sampleN200/scripts/PiDnControl_cicero_parallel.R \
    -f $SEURAT \
    -o $OUTDIR \
    -c $CLUST \
    -g $GENOME \
    -n $SLURM_ARRAY_TASK_ID
