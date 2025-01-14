  #!/bin/bash
  #SBATCH --job-name=bam2fastq      ## Name of the job.
  #SBATCH -A class-nb227     ## account to charge (all labs have a separate account, if not request HPC)
  #SBATCH -p standard          ## partition/queue name
  #SBATCH --nodes=1            ## (-N) number of nodes to use
  #SBATCH --ntasks=1           ## (-n) number of tasks to launch
  #SBATCH --cpus-per-task=12    ## number of cores the job needs
  #SBATCH --mem=32G    ## number of cores the job needs
  #SBATCH --error=slurm-%J.err ## error log file

  module load bedtools2/2.29.2
  module load samtools/1.10

  BASENAME=$1
  INPUT_DIR=/share/crsp/lab/alsppg/share/MIT_R62_ATAC/bams
  OUTPUT_DIR=/data/homezvol2/vswarup/swaruplab/shared_lab/Collaborations/Leslie/R6_ATACseq_Footprint
  BAMFILE=/share/crsp/lab/alsppg/share/MIT_R62_ATAC/bams/${BASENAME}_1.merged.nodup.bam

  cd $OUTPUT_DIR
  mkdir ${OUTPUT_DIR}/${BASENAME}

  cd ${OUTPUT_DIR}/${BASENAME}

  samtools sort -n -o ${BASENAME}_sorted.bam ${BAMFILE}

  bedtools bamtofastq -i ${BASENAME}_sorted.bam -fq ${BASENAME}.R1.fq -fq2 ${BASENAME}.R2.fq
  gzip ${BASENAME}.R1.fq
  gzip ${BASENAME}.R2.fq

##submit

cd /share/crsp/lab/alsppg/share/MIT_R62_ATAC/bams
for folder in *.bam;do
echo $folder;
name=`basename $folder _1.merged.nodup.bam`;
echo $name;
cd /data/homezvol2/vswarup/swaruplab/shared_lab/Collaborations/Leslie/R6_ATACseq_Footprint;
sbatch bam2fastq.sh $name;
done
