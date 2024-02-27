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
      GSC402-Arima-allReps-filtered.mcool  \
      astr-Arima-allReps-filtered.mcool)


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

