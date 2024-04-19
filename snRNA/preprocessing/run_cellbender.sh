#!/bin/bash
#SBATCH --job-name=cb_pid
#SBATCH -A vswarup_lab_gpu
#SBATCH -p gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:A30:1
#SBATCH --cpus-per-task=4
#SBATCH --mem 128G
#SBATCH --array=0-7
#SBATCH --error=slurm-%J.err ## error log file
#SBATCH --output=out_%A_%a.txt  ## output filename
#SBATCH --time=8:00:00

# activate appropriate env
source ~/.bashrc
conda activate cellbender

# get a list of all samples:
counts_dir="/dfs7/swaruplab/zechuas/Projects/AD_GM_WM_2022/snRNA_seq/rawdata/ParseBio/WT2_27samples/WT2_combined/"
out_dir="/dfs7/swaruplab/smorabit/analysis/PiD_2021/snRNA/cellbender/"

let index="$SLURM_ARRAY_TASK_ID"
#echo $index

# get array of all samples
cd $counts_dir

samples=($(ls -d */ | grep "PiD" | uniq))
samples_ctrl=($(ls -d */ | grep "Ctrl" | uniq))

# combine PiD & Ctrl into one list
samples+=(${samples_ctrl[*]})
sample_id=${samples[$index]}
echo $sample_id

echo "${samples[*]}"

# make the output dir:
mkdir $out_dir/$sample_id/

infile="$counts_dir/$sample_id/DGE_unfiltered/"
outfile="$infile/cellbender.h5"
echo $outfile
echo $infile

touch $infile/aaaaa.eeeee

#srun -A vswarup_lab --cpus-per-task=8 --mem 96G --pty bash -i --time=24:00:00

#cd $out_dir/$sample_name/
cd $infile 


if ! [ -f $infile/cellbender.h5 ]; then
    cellbender remove-background \
    --input $infile/adata.h5ad \
    --output $outfile \
    --expected-cells 10000 \
    --total-droplets-included 25000 \
    --epochs 150 \
    --cuda
fi



