#!/bin/bash

# Slurm Parameters
#SBATCH -p medium
#SBATCH --time 8:00:00
#SBATCH --job-name=neoloop
#SBATCH -o  Neoloop_%A.log
#SBATCH -e  Neoloop_%A.err
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem 80G
#SBATCH -C scratch
#SBATCH --array=13

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

res=$1
coolerpath='/scratch1/users/txie/workplace_Pan/coolers-hg38'
SVoutpath='/home/uni08/txie/analysis/Pan/analysis/SV'
plotoutpath='/home/uni08/txie/analysis/Pan/analysis/SV'


######################################call SVs #####################

source ~/miniconda3/etc/profile.d/conda.sh 
conda activate EagleC



merge-redundant-SVs --full-sv-files   \
$SVoutpath/${if%%-*}.CNN_SVs.5K_combined.txt     \
-O $SVoutpath/${if%%-*}.CNN_SVs.5K_combined.merge    \
--output-format NeoLoopFinder


annotate-gene-fusion --sv-file $SVoutpath/${if%%-*}.CNN_SVs.5K_combined.txt    \
 --output-file $SVoutpath/${if%%-*}.gene-fusions_5k.txt    \
--buff-size 5000 --skip-rows 1    \
--ensembl-release 93 --species human




python 11plot_SV.py  $coolerpath $SVoutpath $plotoutpath


################################END##call SV #####################


####################################call neoloop and neotad################################
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate neoloop043


assemble-complexSVs     \
-H $coolerpath/$if::/resolutions/5000   \
$coolerpath/$if::/resolutions/10000   \
$coolerpath/$if::/resolutions/25000   \
-B $SVoutpath/${if%%-*}.CNN_SVs.5K_combined.merge \
--balance-type CNV   \
--nproc 10   \
--minimum-size  10  \
-O $SVoutpath/${if%%-*}_CNV  \
--logFile $SVoutpath/assembleSVs_cnv.log 


neotad-caller -H $coolerpath/$if::/resolutions/25000   \
--assembly $SVoutpath/${if%%-*}_CNV.assemblies.txt   \
--balance-type CNV   \   ####important
--nproc 15   \
-O $SVoutpath/${if%%-*}_25K.neotad  \   #####
--logFile $SVoutpath/neotad.log 


neoloop-caller \
-H    \
$coolerpath/$if::/resolutions/5000   \
$coolerpath/$if::/resolutions/10000   \
$coolerpath/$if::/resolutions/25000   \
--assembly $SVoutpath/${if%%-*}_CNV.assemblies.txt   \
--balance-type CNV  --nproc 15 --prob 0.95  \
-O $SVoutpath/${if%%-*}.'95'.Mneoloop  \
--logFile $SVoutpath/neoloop.log 


################################END##call neoloop and neotad################################






