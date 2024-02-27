#!/bin/bash

# Slurm Parameters
#SBATCH -p medium
#SBATCH --time 2-00:00:00
#SBATCH --job-name=STAR_index
#SBATCH -c 10
#SBATCH --mem 40G
#SBATCH -e STAR_index_%A.err
#SBATCH -o STAR_index_%A.log
#SBATCH -C scratch2



module load star/2.6.0c
module load rsem/1.3.0
rnadir='/scratch/users/txie/RNA-seq'

ifs=(GSC457B-Arima-allReps-filtered.mcool   \
     GSC390B-Arima-allReps-filtered.mcool  \
      GSC28-Arima-allReps-filtered.mcool  \
      GSC275B-Arima-allReps-filtered.mcool  \
      GSC275-Arima-allReps-filtered.mcool  \
      GSC23p-Arima-allReps-filtered.mcool  \
      GSC213R1-Arima-allReps-filtered.mcool  \
      GSC208R1-Arima-allReps-filtered.mcool  \
      GSC1-Arima-allReps-filtered.mcool  \
      GSC163-Arima-allReps-filtered.mcool  \
      GSC120-Arima-allReps-filtered.mcool  \
      GSC148-Arima-allReps-filtered.mcool  \
      GSC171-Arima-allReps-filtered.mcool  \
      GSC181-Arima-allReps-filtered.mcool  \
      GSC208M-Arima-allReps-filtered.mcool  \
      GSC213M-Arima-allReps-filtered.mcool  \
      GSC318-Arima-allReps-filtered.mcool  \
      GSC323-Arima-allReps-filtered.mcool  \
      GSC394B-Arima-allReps-filtered.mcool \
      GSC412-Arima-allReps-filtered.mcool \
      GSC450-Arima-allReps-filtered.mcool  \
      GSC452C-Arima-allReps-filtered.mcool \
      GSC452P-Arima-allReps-filtered.mcool  \
      GSC61-Arima-allReps-filtered.mcool  \
      GSC62-Arima-allReps-filtered.mcool  \
      GSC83-Arima-allReps-filtered.mcool  \
      GSC208R2-Arima-allReps-filtered.mcool  \
      GSC213R2-Arima-allReps-filtered.mcool  \
      GSC351-Arima-allReps-filtered.mcool  \
      GSC486-Arima-allReps-filtered.mcool  \
      GSC428-Arima-allReps-filtered.mcool  \
      GSC402-Arima-allReps-filtered.mcool  )




if=${ifs[${SLURM_ARRAY_TASK_ID}]}
id=${if%%-*}


###############################################bulid STAR index

STAR --runThreadN 10  \
--runMode genomeGenerate  \
--genomeDir  $rnadir/V22/   \
--genomeFastaFiles $rnadir/0ref/hg38.fa  \
--sjdbGTFfile $rnadir/0ref/merged_annotation_V29.gtf    \
--sjdbOverhang 150  

###############################################align reads with STAR


STAR --runMode alignReads  --runThreadN 12    \
--readFilesCommand zcat  \
--genomeDir   ${path}/v29    \
--readFilesIn  $rnadir/1rawdata/${id}.fastq.gz  \
--outFileNamePrefix  ${path}/2align_out/${id##*/}.   \
--outSAMtype BAM SortedByCoordinate   \
--outBAMsortingThreadN 12    \
--quantMode TranscriptomeSAM GeneCounts  


###############################################Estimate gene and isoform expression from RNA-Seq data.

rsem-calculate-expression  --no-bam-output --alignments  \
--strandedness   reverse  \
-p 8  --append-names        \
#--paired-end                 \
${rnadir}/2align_out/${id}.Aligned.toTranscriptome.out.bam   \
${rnadir}/0ref/RSEM_hg38/RSEM_hg38  \
${rnadir}/3rsem_out/${id}  

