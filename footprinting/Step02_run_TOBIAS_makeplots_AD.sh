#!/bin/bash
#SBATCH --job-name=tobias      ## Name of the job.
#SBATCH -A vswarup_lab     ## account to charge (all labs have a separate account, if not request HPC)
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=6    ## number of cores the job needs
#SBATCH --mem=32G    ## number of cores the job needs
#SBATCH --error=slurm-%J.err ## error log file

source /data/homezvol2/vswarup/project/miniconda3/etc/profile.d/conda.sh
conda activate tobias

BAMFILE=$1
NAME=$2
INDIR=/dfs7/swaruplab/smorabit/analysis/AD_NucSeq_2019/atac_analysis/all_data/trackhubs/data/bams/merged/celltypes
OUTDIR=/dfs7/swaruplab/smorabit/analysis/PiD_2021/tobias/AD
BAMFILE=${INDIR}/${BAMFILE}

GENOME=/dfs7/swaruplab/shared_lab/ExternalDatasets/tobias_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa ## must be indexed, samtools faidx
MERGED_PEAKS=/dfs7/swaruplab/smorabit/analysis/PiD_2021/trackhubs/data/bams/merged/pid_peaks_merged.bed #### From Seurat
MERGED_PEAKS_HEADER=/dfs7/swaruplab/smorabit/analysis/PiD_2021/trackhubs/data/bams/merged/merged_peaks_header.txt ## header
BLACKLIST=/dfs7/swaruplab/shared_lab/ExternalDatasets/tobias_ref/ENCODE_hg38_blacklist.bed #from https://sites.google.com/site/anshulkundaje/projects/blacklists
JASPER_MOTIFS=/dfs7/swaruplab/shared_lab/ExternalDatasets/tobias_ref/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.motifs
TFBS_List=${INDIR}/TF_list2022_JASPAR.txt ##generated in R using the output of directories from TOBIAS BINDetect

TFBS_List=${OUTDIR}/TF_list2022_JASPAR.txt ##generated in R using the output of directories from TOBIAS BINDetect

cd ${OUTDIR}/${NAME}_ATACorrect

for TF in $(cat ${TFBS_List}); do
  echo $TF;
  TF_BED=${OUTDIR}/${NAME}_ATACorrect/BINDetect_single_output/${TF}/beds/${TF}_all.bed;
  BOUND_BED=${OUTDIR}/${NAME}_ATACorrect/BINDetect_single_output/${TF}/beds/${TF}_*_bound.bed;
  UNBOUND_BED=${OUTDIR}/${NAME}_ATACorrect/BINDetect_single_output/${TF}/beds/${TF}_*_unbound.bed;
  TOBIAS PlotAggregate --TFBS ${TF_BED}  --signals ${NAME}_corrected.bw --output ${TF}_footprint_${NAME}.pdf --share_y both --plot_boundaries --signal-on-x;
  TOBIAS PlotAggregate --TFBS ${TF_BED} ${BOUND_BED} ${UNBOUND_BED} --signals ${NAME}_uncorrected.bw ${NAME}_expected.bw ${NAME}_corrected.bw --output ${TF}_footprint_split_${NAME}.pdf --share_y sites --plot_boundaries;
done

mkdir ${NAME}_footprint_pdfs
mv *.pdf ${NAME}_footprint_pdfs

##submit

cd /dfs7/swaruplab/smorabit/analysis/AD_NucSeq_2019/atac_analysis/all_data/trackhubs/data/bams/merged/celltypes
for folder in *.bam;do
echo $folder;
name=`basename $folder .bam`;
echo $name;
cd /dfs7/swaruplab/smorabit/analysis/PiD_2021/tobias/AD;
sbatch Step02_run_tobias_makeplots.sh $folder $name;
done
