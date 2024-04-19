#!/bin/bash
#SBATCH --job-name=FindMarkers
#SBATCH -p standard
#SBATCH -A vswarup_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --error=slurm-%J.err
#SBATCH --mem 96G
#SBATCH --array=1-8 # number of clusters
#SBATCH --time=4:00:00

source ~/.bashrc
conda activate hdWGCNA # activate your conda environment

# seurat object path:
seurat="/dfs7/swaruplab/smorabit/analysis/PiD_2021/snRNA/data/PiD_snRNA_seurat.rds"
outdir="/dfs7/swaruplab/smorabit/analysis/PiD_2021/snRNA/DEGs/data/celltype_markers/"
type="markers" # type of DEGs to run. "markers" or "condition"
name="cell_type" # name to append to output files
cluster="cell_type" # column in the seurat object with cluster / cell type information
latent="total_counts,SampleID"

mkdir $outdir

# launch R script:
Rscript --vanilla /dfs7/swaruplab/smorabit/analysis/PiD_2021/snRNA/bin/parallel_DEGs.R \
    --seurat $seurat \
    --outdir $outdir \
    --type $type \
    --cluster $cluster \
    --name $name \
    --index $SLURM_ARRAY_TASK_ID \
    --test "MAST" \
    --pos "FALSE" \
    --pct 0.1 \
    --logfc 0.1 \
    --verbose "TRUE" \
    --latent $latent \
    --cores 16
