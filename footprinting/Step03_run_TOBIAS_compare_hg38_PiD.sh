#!/bin/bash
#SBATCH --job-name=tobias      ## Name of the job.
#SBATCH -A vswarup_lab     ## account to charge (all labs have a separate account, if not request HPC)
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=12    ## number of cores the job needs
#SBATCH --mem=32G    ## number of cores the job needs
#SBATCH --error=slurm-%J.err ## error log file

source /data/homezvol2/vswarup/project/miniconda3/etc/profile.d/conda.sh
conda activate tobias


CellType=$1
NAME1=${CellType}_control
NAME2=${CellType}_PiD
COMPARISON=PiDvscontrol_${CellType}

INDIR=/dfs7/swaruplab/smorabit/analysis/PiD_2021/tobias
OUTDIR=${INDIR}/footprint_compare
SIGNAL1=${INDIR}/${NAME1}_ATACorrect/${NAME1}_footprints.bw
SIGNAL2=${INDIR}/${NAME2}_ATACorrect/${NAME2}_footprints.bw
CORRECTED1=${INDIR}/${NAME1}_ATACorrect/${NAME1}_corrected.bw
CORRECTED2=${INDIR}/${NAME2}_ATACorrect/${NAME2}_corrected.bw

GENOME=/dfs7/swaruplab/shared_lab/ExternalDatasets/tobias_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa ## must be indexed, samtools faidx
MERGED_PEAKS=/dfs7/swaruplab/smorabit/analysis/PiD_2021/trackhubs/data/bams/merged/pid_peaks_merged.bed ## merging all bedfiles from macs2 run_macs2.sh file
MERGED_PEAKS_HEADER=/dfs7/swaruplab/smorabit/analysis/PiD_2021/trackhubs/data/bams/merged/merged_peaks_header.txt ## header
JASPER_MOTIFS=/dfs7/swaruplab/shared_lab/ExternalDatasets/tobias_ref/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.motifs

TFBS_List=${INDIR}/TF_list2022_JASPAR.txt ##generated in R using the output of directories from TOBIAS BINDetect


cd ${OUTDIR}

TOBIAS BINDetect --motifs ${JASPER_MOTIFS} --signals ${SIGNAL1} ${SIGNAL2} --genome ${GENOME} --peaks ${MERGED_PEAKS} --peak_header ${MERGED_PEAKS_HEADER} --outdir BINDetect_output_${COMPARISON} --cores 12

mkdir ${OUTDIR}/BINDetect_output_${COMPARISON}
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



####Submit code

for folder in *_control_ATACorrect; do
echo $folder;
name=`basename $folder _control_ATACorrect`;
echo $name;
sbatch ./Step03_run_tobias_compare.sh $name
done




ASC_control_ATACorrect
ASC
Submitted batch job 15479772
EX_control_ATACorrect
EX
Submitted batch job 15479773
INH_control_ATACorrect
INH
Submitted batch job 15479774
MG_control_ATACorrect
MG
Submitted batch job 15479775
ODC_control_ATACorrect
ODC
Submitted batch job 15479776
OPC_control_ATACorrect
OPC
Submitted batch job 15479777
PER.END_control_ATACorrect
PER.END
Submitted batch job 15479778
