#!/bin/bash

# Slurm Parameters
#SBATCH -p medium
#SBATCH --time 8:00:00
#SBATCH --job-name=neoloop
#SBATCH -o  neoloop_%A.log
#SBATCH -e  neoloop_%A.err
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem 20G
#SBATCH -C scratch
#SBATCH --array=0-14

#THIS SCRIPT CAN PRODUCE ALL CNV/CNV-bw/SV/filtedSV/NEO-LOOP/NEOTADS/deleted or duplicated or fusion genes.###################
#################################things to put########################
ifs=(GSC457B-Arima-allReps-filtered.mcool   \
     GSC390B-Arima-allReps-filtered.mcool  \
      GSC28-Arima-allReps-filtered.mcool  \
      GSC275B-Arima-allReps-filtered.mcool  \
      GSC275-Arima-allReps-filtered.mcool  \
      GSC23p-Arima-allReps-filtered.mcool  \
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
      GSC351-Arima-allReps-filtered.mcool  \
      GSC486-Arima-allReps-filtered.mcool  \
      GSC428-Arima-allReps-filtered.mcool  \
      GSC402-Arima-allReps-filtered.mcool  )



if=${ifs[${SLURM_ARRAY_TASK_ID}]}

coolerpath='/scratch1/users/txie/workplace_Pan/coolers-hg38'
SVoutpath='/home/uni08/txie/analysis/Pan/analysis/SV'

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate EagleC


predictSV    \
--hic-5k $coolerpath/$if::/resolutions/5000   \
--hic-10k  $coolerpath/$if::/resolutions/10000   \
--hic-50k  $coolerpath/$if::/resolutions/50000    \
-O $SVoutpath/${if%%-*}_NOkb -g hg38     \
--balance-type CNV  --output-format full    \
--prob-cutoff-5k 0.8 --prob-cutoff-10k 0.8 --prob-cutoff-50k 0.99999    

