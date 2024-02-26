#!/bin/bash

# Slurm Parameters
#SBATCH -p medium
#SBATCH --time 8:00:00
#SBATCH --job-name=neoloop
#SBATCH -o  Neoloop_%A.log
#SBATCH -e  Neoloop_%A.err
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem 50G
#SBATCH -C scratch
#SBATCH --array=2

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

R=$1
res1=25000
res2=50000
res3=10000
res4=5000
coolerpath='/scratch1/users/txie/workplace_Pan/coolers-hg38'
SVoutpath='/home/uni08/txie/analysis/Pan/analysis/SV'
plotoutpath='/home/uni08/txie/analysis/Pan/analysis/SV'

######################################################call SV####################################
source ~/miniconda3/etc/profile.d/conda.sh 

conda activate neoloop043

for R in $res1 $res2 $res3 $res4;   \
do    \
calculate-cnv -H $coolerpath/$if::/resolutions/$R  \
-g hg38 -e Arima  \
--output $SVoutpath/${if%%-*}_$[R/1000]kb.1 --cachefolder $SVoutpath/  \
--logFile $SVoutpath/cnv-calculation1.log


segment-cnv --cnv-file $SVoutpath/${if%%-*}_$[R/1000]kb.1   \
--bins $R --nproc 16 \
--output $SVoutpath/${if%%-*}_$[R/1000]kb.2 \
--logFile $SVoutpath/lastk.log    \
  -f    \




done
command

correct-cnv -H $coolerpath/$if::/resolutions/$R  --cnv-file $SVoutpath/${if%%-*}_$[R/1000]kb.2  \
--nproc 16    \
--logFile $SVoutpath/cnv-norm1.log   


source ~/miniconda3/etc/profile.d/conda.sh 
conda activate base
sort -k1,1 -k2,2n  $SVoutpath/${if%%-*}_$[res/1000]kb.1  >$SVoutpath/${if%%-*}_$[res/1000]kb.1.sorted

bedGraphToBigWig $SVoutpath/${if%%-*}_$[res/1000]kb.1.sorted      \
/scratch/users/txie/data/hg38/hg38.chrom.sizes $SVoutpath/${if%%-*}_$[res/1000]kb.bw    

rm $SVoutpath/${if%%-*}_$[res/1000]kb.1.sorted


#############################################################this part is to map the CNV profile to gene #######to find deleted or duplicated gene##########3

source ~/miniconda3/etc/profile.d/conda.sh 
conda deactivate 

python 10cal_CNV.py   ${if%%-*} ${SVoutpath}

command
