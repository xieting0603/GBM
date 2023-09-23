#!/bin/bash

# Slurm Parameters
#SBATCH -p medium
#SBATCH --time 2-00:00:00
#SBATCH --job-name=predictSV
#SBATCH -o predictSV_%A.log
#SBATCH -e predictSV_%A.err
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem 80G
#SBATCH -C scratch
#SBATCH --array=0-31

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

transdir=/home/uni08/txie/analysis/CP/8trans
loc=/scratch1/users/txie/workplace_CP/coolers-hg38
dir=/home/uni08/txie/analysis/CP/8trans/cnv

################calculate cnv and do cnv-based normalization
calculate-cnv -H $loc/$if::/resolutions/$res  \
-g hg38 -e Arima  \
--output $dir/${if%%-*}_$[res/1000]kb.1 --cachefolder $dir/  \
--logFile $dir/cnv-calculation.log

segment-cnv --cnv-file $dir/${if%%-*}_$[res/1000]kb.1   \
--bins $res --nproc 16 \
--output $dir/${if%%-*}_$[res/1000]kb.2 \
--logFile $dir/cnv-seg.log

correct-cnv -H $loc/$if::/resolutions/$res  --cnv-file $dir/${if%%-*}_$[res/1000]kb.2  \
--nproc 16  -f  \
--logFile $dir/cnv-norm.log

sort -k1,1 -k2,2n  $dir/${if%%-*}_$[res/1000]kb.1  >$dir/${if%%-*}_$[res/1000]kb.1.sorted
bedGraphToBigWig $dir/${if%%-*}_$[res/1000]kb.1.sorted  /scratch2/txie/data/hg38/hg38.chrom.sizes $dir/${if%%-*}_$[res/1000]kb.bw
rm $dir/${if%%-*}_$[res/1000]kb.1.sorted


###########################detecte SVs and fusion events

predictSV --hic-5k $loc/$if::/resolutions/5000 \
            --hic-10k  $loc/$if::/resolutions/10000 \
            --hic-50k  $loc/$if::/resolutions/50000 -O ${if%%-*} -g hg38 --balance-type CNV --output-format full \
            --prob-cutoff-5k 0.8 --prob-cutoff-10k 0.8 --prob-cutoff-50k 0.99999

merge-redundant-SVs --full-sv-files ${if%%-*}.CNN_SVs.5K_combined.txt -O ${if%%-*}.CNN_SVs.5K_combined.merge  --output-format NeoLoopFinder

annotate-gene-fusion --sv-file $transdir/translocation/raw/${if%%-*}.CNN_SVs.5K_combined.txt \
                       --output-file $transdir/translocation/${if%%-*}.gene-fusions_5k.txt \
                       --buff-size 5000 --skip-rows 1 --ensembl-release 93 --species human

###########################assemble and filteri SVs 


assemble-complexSVs     \
-H $loc/$if::/resolutions/5000   \
$loc/$if::/resolutions/10000   \
$loc/$if::/resolutions/25000   \
-B $transdir/translocation/${if%%-*}.CNN_SVs.5K_combined.merge \
--balance-type CNV --nproc 15  \
--minimum-size  10  \
-O $transdir/ass_neoloop/${if%%-*}_CNV  \
--logFile $transdir/ass_neoloop/assembleSVs_cnv.log 

###########################detecte neo-loops and neo-tads


neoloop-caller \
-H    \
$loc/$if::/resolutions/5000   \
$loc/$if::/resolutions/10000   \
$loc/$if::/resolutions/25000   \
--assembly $transdir/ass_neoloop_beforepeakachu/${if%%-*}_CNV.assemblies.txt   \
--balance-type CNV  --nproc 15 --prob 0.95  \
-O $transdir/ass_neoloop_beforepeakachu/${if%%-*}.'95'.Mneoloop  \
--logFile $transdir/ass_neoloop_beforepeakachu/neoloop.log 
command

neotad-caller -H $loc/$if::/resolutions/25000   \
--assembly $transdir/ass_neoloop_beforepeakachu/${if%%-*}_CNV.assemblies.txt   \
--balance-type CNV  --nproc 15 \
-O $transdir/ass_neoloop_beforepeakachu/${if%%-*}_25K.neotad  \
--logFile $transdir/ass_neoloop_beforepeakachu/neotad.log 



