#!/bin/bash

# Slurm Parameters
#SBATCH -p medium
#SBATCH --time 8:00:00
#SBATCH --job-name=neoloop
#SBATCH -o HIC_pipeline%A.log
#SBATCH -e HIC_pipeline_%A.err
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem 50G
#SBATCH -C scratch
#SBATCH --array=0-32

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

########################################### generate cooler files with runhic
conda activate runHiC86
cd /scratch1/users/txie/workplace_CP
runHiC pileup -p ../data/ -g hg38 -f tumor_seq -F FASTQ -A bwa-mem -t 10 -O BAM --logFile runhic22.log --tmpdir /scratch1/users/txie  --chunkSize 60000000 --memory 50G
runHiC quality -L filtered-hg38
##################################################get CNV and CNV-normalized cooler
conda activate neoloop
if=${ifs[${SLURM_ARRAY_TASK_ID}]}
res=5000
Coolerloc=/scratch1/users/txie/workplace_CP/coolers-hg38
cnvdir=/home/uni08/txie/analysis/CP/8trans/cnv
Compartmentdir=/home/uni08/txie/analysis/CP/04Compartment/cooltools/tmp
CNVCoolerdir=/scratch1/users/txie/workplace_CP/coolers-hg38/sweight
taddir=/home/uni08/txie/analysis/CP/05TAD
loopsdir=/home/uni08/txie/analysis/CP/06dot/peakchu


calculate-cnv -H $Coolerloc/$if::/resolutions/$res  \
-g hg38 -e Arima  \
--output $cnvdir/${if%%-*}_$[res/1000]kb.1 --cachefolder $cnvdir/  \
--logFile $cnvdir/cnv-calculation.log

segment-cnv --cnv-file $cnvdir/${if%%-*}_$[res/1000]kb.1   \
--bins $res --nproc 16 \
--output $cnvdir/${if%%-*}_$[res/1000]kb.2 \
--logFile $cnvdir/cnv-seg.log

correct-cnv -H $Coolerloc/$if::/resolutions/$res  --cnv-file $cnvdir/${if%%-*}_$[res/1000]kb.2  \
--nproc 16  -f  \
--logFile $cnvdir/cnv-norm.log

sort -k1,1 -k2,2n  $cnvdir/${if%%-*}_$[res/1000]kb.1  >$cnvdir/${if%%-*}_$[res/1000]kb.1.sorted
bedGraphToBigWig $cnvdir/${if%%-*}_$[res/1000]kb.1.sorted  /scratch2/txie/data/hg38/hg38.chrom.sizes $cnvdir/${if%%-*}_$[res/1000]kb.bw
rm $cnvdir/${if%%-*}_$[res/1000]kb.1.sorted
###################################################################

###################################################identify compartment
binsize=50000
binsizek=$[binsize/1000]

cooltools call-compartments   \
$CNVCoolerdir/$if::/resolutions/$binsize  \
--reference-track  ../0base_files/hg38_${binsizek}kb.gc  \
--bigwig  \
-o $Compartmentdir/${if%%-*}_${binsizek}kb

###################################################identify TADs and calculate insulation score
binsize=25000
binsizek=$[binsize/1000]
cooltools insulation    \
$CNVCoolerdir/$if::/resolutions/$binsize 125000 \
 --threshold  Li  \
-o $TADdir/${if%%-*}_${binsizek}kb.tsv
#############################################################################call loops
binsize=5000
binsizek=$[binsize/1000]
thre1=0.95



A=$(peakachu depth -p ${CNVCoolerdir}/${if}::/resolutions/1000000 | gawk '{print $0}' )
a=${A##*:}
b=${a:1:3}
peakachu score_genome -r $binsize --balance -p ${CNVCoolerdir}/${if}::/resolutions/$binsize -O ${loopsdir}/${if%%-*}_${binsizek}kb -m ./peakachu_models/high-confidence.${b}million.${binsizek}kb.w6.pkl

peakachu pool -r $binsize -i  ${loopsdir}/${if%%-*}_${binsizek}kb  -o ${loopsdir}/${if%%-*}_${binsizek}kb.${thre1}.loops -t $thre1



