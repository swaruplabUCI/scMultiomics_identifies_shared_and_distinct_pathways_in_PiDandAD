#!/bin/bash
#SBATCH --job-name=tobias      ## Name of the job.
#SBATCH -A class-nb227     ## account to charge (all labs have a separate account, if not request HPC)
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=12    ## number of cores the job needs
#SBATCH --mem=32G    ## number of cores the job needs
#SBATCH --error=slurm-%J.err ## error log file

source /data/homezvol2/vswarup/miniconda3/etc/profile.d/conda.sh
conda activate tobias

NAME1=WT.STR.8W.Pos
NAME2=HD.STR.8W.Pos
COMPARISON=HDvsWT_STR.8W.Pos

INDIR=/data/homezvol2/vswarup/swaruplab/shared_lab/Collaborations/Leslie/R6_ATACseq_Footprint
OUTDIR=${INDIR}/footprint_compare
SIGNAL1=${INDIR}/footprints/${NAME1}_ATACorrect/${NAME1}_footprints.bw
SIGNAL2=${INDIR}/footprints/${NAME2}_ATACorrect/${NAME2}_footprints.bw
CORRECTED1=${INDIR}/footprints/${NAME1}_ATACorrect/${NAME1}_corrected.bw
CORRECTED2=${INDIR}/footprints/${NAME2}_ATACorrect/${NAME2}_corrected.bw

GENOME=/data/homezvol2/vswarup/project/bin/TOBIAS/tobias/mm10.fa ## must be indexed, samtools faidx
MERGED_PEAKS=${INDIR}/Merged_bed/merged_peaks.bed ## merging all bedfiles from macs2 run_macs2.sh file
MERGED_PEAKS_HEADER=${INDIR}/Merged_bed/merged_peaks_header.txt ## header
JASPER_MOTIFS=${INDIR}/motifs_2022.jaspar
TFBS_List=${INDIR}/TF_list2022.txt ##generated in R using the output of directories from TOBIAS BINDetect


cd ${OUTDIR}

TOBIAS BINDetect --motifs ${JASPER_MOTIFS} --signals ${SIGNAL1} ${SIGNAL2} --genome ${GENOME} --peaks ${MERGED_PEAKS} --peak_header ${MERGED_PEAKS_HEADER} --outdir BINDetect_output_${COMPARISON} --cores 12

cd ${OUTDIR}/BINDetect_output_${COMPARISON}

for TF in $(cat ${TFBS_List}); do
  echo $TF;
  TF_BED=${OUTDIR}/BINDetect_output_${COMPARISON}/${TF}/beds/${TF}_all.bed;
  BOUND_BED1=${OUTDIR}/BINDetect_output_${COMPARISON}/${TF}/beds/${TF}_${NAME1}_footprints_bound.bed;
  BOUND_BED2=${OUTDIR}/BINDetect_output_${COMPARISON}/${TF}/beds/${TF}_${NAME2}_footprints_bound.bed;

  TOBIAS PlotAggregate --TFBS ${TF_BED} --signals ${CORRECTED1} ${CORRECTED2} --output ${TF}_footprint_${COMPARISON}.pdf --share_y both --plot_boundaries --signal-on-x;
  TOBIAS PlotAggregate --TFBS ${BOUND_BED1} ${BOUND_BED2} --signals ${CORRECTED1} ${CORRECTED2} --output ${TF}_footprint_${COMPARISON}_boundOnly.pdf --share_y both --plot_boundaries;

done

mkdir ${COMPARISON}_footprint_pdfs
mv *.pdf ${COMPARISON}_footprint_pdfs
