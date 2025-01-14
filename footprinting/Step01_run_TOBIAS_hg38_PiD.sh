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

BAMFILE=$1
NAME=$2
INPUT_DIR=/dfs7/swaruplab/smorabit/analysis/PiD_2021/trackhubs/data/bams/merged
OUTPUT_DIR=/dfs7/swaruplab/smorabit/analysis/PiD_2021/tobias
BAMFILE=${INPUT_DIR}/${BAMFILE}

GENOME=/dfs7/swaruplab/shared_lab/ExternalDatasets/tobias_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa ## must be indexed, samtools faidx
MERGED_PEAKS=${INPUT_DIR}/pid_peaks_merged.bed #### From Seurat

MERGED_PEAKS_HEADER=${INPUT_DIR}/merged_peaks_header.txt ## header
BLACKLIST=/dfs7/swaruplab/shared_lab/ExternalDatasets/tobias_ref/ENCODE_hg38_blacklist.bed #from https://sites.google.com/site/anshulkundaje/projects/blacklists
JASPER_MOTIFS=/dfs7/swaruplab/shared_lab/ExternalDatasets/tobias_ref/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.motifs


cd ${OUTPUT_DIR}


TOBIAS ATACorrect --bam ${BAMFILE} --genome ${GENOME} --peaks ${MERGED_PEAKS} --blacklist ${BLACKLIST} --outdir ${NAME}_ATACorrect --cores 12
cd ${NAME}_ATACorrect
TOBIAS FootprintScores --signal ${NAME}_corrected.bw --regions ${MERGED_PEAKS} --output ${NAME}_footprints.bw --cores 12
TOBIAS BINDetect --motifs ${JASPER_MOTIFS} --signals ${NAME}_footprints.bw --genome ${GENOME} --peaks ${MERGED_PEAKS} --peak_header ${MERGED_PEAKS_HEADER} --outdir BINDetect_single_output --cores 12



##submit

INPUT_DIR=/dfs7/swaruplab/smorabit/analysis/PiD_2021/trackhubs/data/bams/merged
OUTPUT_DIR=/dfs7/swaruplab/smorabit/analysis/PiD_2021/tobias
cd ${INPUT_DIR}
for folder in *.bam;do
echo $folder;
name=`basename $folder .bam`;
echo $name;
cd ${OUTPUT_DIR};
sbatch run_tobias.sh $folder $name;
done

##make bamfile list in R
options(stringsAsFactors=F)
targets=read.csv("metadata_footprinting_r6.csv")
filename="/share/crsp/lab/alsppg/share/MIT_R62_ATAC/bams/"
endname="_1.merged.nodup.bam"

group=unique(targets$Group)
for (i in 1:length(group)){
  gr1=group[i]
  samples=targets[targets$Group==gr1,"title"]
  sample.names=paste0(filename,samples,endname)
  write.table(sample.names,file=paste0(gr1,"_bamlist.txt"),row.names=F,col.names =F, quote=F)
}



###Merge bamfiles
#!/bin/bash
#SBATCH --job-name=bamtools      ## Name of the job.
#SBATCH -A class-nb227     ## account to charge (all labs have a separate account, if not request HPC)
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=12    ## number of cores the job needs
#SBATCH --mem=32G    ## number of cores the job needs
#SBATCH --error=slurm-%J.err ## error log file

module load bamtools/2.5.1

OUTPATH=/data/homezvol2/vswarup/swaruplab/shared_lab/Collaborations/Leslie/R6_ATACseq_Footprint/merged_bams
LIST1=$1
BASENAME=$2
MERGELIST=/data/homezvol2/vswarup/swaruplab/shared_lab/Collaborations/Leslie/R6_ATACseq_Footprint/${LIST1}

cd ${OUTPATH}

bamtools merge -list ${MERGELIST} -out ${BASENAME}.bam
bamtools index -in ${BASENAME}.bam


for file in *_bamlist.txt; do
  name=`basename $file _bamlist.txt`;
  echo $name;
  sbatch mergeBamfiles.sh $file $name
done




##cat WT_peaks.bed treatment_peaks.bed | bedtools sort | bedtools merge > merged_peaks.bed
#######merging bedfiles
# cd /dfs7/swaruplab/smorabit/analysis/PiD_2021/trackhubs/data/bams/merged
# for files in *.bam; do
#   echo $files
#   name=`basename $files .bam`
#   echo $name
#   bedtools bamtobed -i $files |cut -f1,2,3 > ${name}.bed
# done
#
# cat ASC_control.bed ASC_PiD.bed EX_control.bed EX_PiD.bed INH_control.bed INH_PiD.bed MG_control.bed MG_PiD.bed ODC_control.bed ODC_PiD.bed OPC_control.bed OPC_PiD.bed PER.END_control.bed PER.END_PiD.bed | bedtools sort | bedtools merge > merged_peaks.bed
# #bedtools sort -i merged.bed > merged_sort.bed
# cd /dfs7/swaruplab/smorabit/analysis/PiD_2021/trackhubs/data/bams/merged
# sort -S 60% --parallel=10 -k 1,1 -k2,2n merged.bed > merged_sort.bed
# bedtools merge -i merged_sort.bed > merged_final.bed

#cat pid_peaks.bed | bedtools sort | bedtools merge > pid_peaks_merged.bed
