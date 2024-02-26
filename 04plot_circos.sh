#!/bin/bash

# Slurm Parameters
#SBATCH -p gpu
#SBATCH --time 00:30:00
#SBATCH --job-name=plot_circos
#SBATCH -o plot_circos_%A.log
#SBATCH -e plot_circos_%A.err
#SBATCH -c 1
#SBATCH --mem 20G
#SBATCH -C scratch
#SBATCH --array=10-14
#SBATCH --gpus 1

#######################THIS SCRIPT CAN PRODUCE CIRCOS PLOT AS WELL AS INTER OR INTRA SVS DETECED BY EAGLEC###############

coolerpath='/scratch1/users/txie/workplace_Pan/coolers-hg38'
plotoutpath='/home/uni08/txie/analysis/Pan/analysis/SV'
SVdir='/home/uni08/txie/analysis/Pan/analysis/SV'

ifs=(Pan3T-Arima-allReps-filtered.mcool   \
     Pan3T_ICM-Arima-allReps-filtered.mcool  \
    Pan51T-Arima-allReps-filtered.mcool  \
    Pan74T-Arima-allReps-filtered.mcool  \
    Pan84T-Arima-allReps-filtered.mcool  \
    Pan35T-Arima-allReps-filtered.mcool  \
    Pan51T_ICM-Arima-allReps-filtered.mcool  \
    Pan54T-Arima-allReps-filtered.mcool    \
    Pan61T-Arima-allReps-filtered.mcool   \
    Pan2T-Arima-allReps-filtered.mcool   \
    Pan51T_T-Arima-allReps-filtered.mcool  \
    Pan61T_T-Arima-allReps-filtered.mcool  \
    Pan54T_T-Arima-allReps-filtered.mcool  \
    Pan74T_T-Arima-allReps-filtered.mcool  \
    Pan84T_T-Arima-allReps-filtered.mcool  )


if=${ifs[${SLURM_ARRAY_TASK_ID}]}




sample=${if%%-*}


###circos plot##########

python 04generate_cnv.py $SVdir  ####### to convert neoloop CNV segment to circos format CNV #####################
python 04ass2circos_type.py $SVdir   ########to convert EagleC SV to circos format SV##############


source ~/miniconda3/etc/profile.d/conda.sh 
conda activate base
conda deactivate 
module load circos/0.69
##########the cicos folder with other .conf files and circos.conf should be in the same folder

circos  -conf 04circos.conf     \
--param links/link/file=${SVdir}/${sample}.forcircos.type   \
--param plots/plot/file=${SVdir}/${sample}_25kb.3   \
--param ideogram/spacing/show_labels=yes   \
--param show_ticks=yes  \
--param show_tick_labels=no  \
-outputfile  ${sample}


module unload circos/0.69

