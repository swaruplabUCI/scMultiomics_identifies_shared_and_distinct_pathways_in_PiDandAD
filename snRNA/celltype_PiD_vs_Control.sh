#!/bin/bash
#SBATCH --job-name=ct_PiD
#SBATCH -p standard
#SBATCH -A vswarup_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --error=slurm-%J.err
#SBATCH --mem 96G
#SBATCH --array=1-8
#SBATCH --time=12:00:00

source ~/.bashrc
conda activate hdWGCNA

# seurat object path:
seurat="/dfs7/swaruplab/smorabit/analysis/PiD_2021/snRNA/data/PiD_snRNA_seurat.rds"
outdir="/dfs7/swaruplab/smorabit/analysis/PiD_2021/snRNA/DEGs/data/celltype_DX/"
type="conditions"
name="celltype_PiD_vs_Control"
condition="DX"
group1="PiD"
group2="Control"
cluster="cell_type"
latent="Sex,total_counts,SampleID"

mkdir $outdir

# launch R script:
Rscript --vanilla /dfs7/swaruplab/smorabit/analysis/PiD_2021/snRNA//bin/parallel_DEGs.R \
    --seurat $seurat \
    --outdir $outdir \
    --type $type \
    --cluster $cluster \
    --condition $condition \
    --name $name \
    --index $SLURM_ARRAY_TASK_ID \
    --group1 $group1 \
    --group2 $group2 \
    --test "MAST" \
    --pos "FALSE" \
    --pct 0 \
    --logfc 0 \
    --verbose "TRUE" \
    --latent $latent \
    --cores 16 \
    --overwrite "TRUE"
