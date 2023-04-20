#!/bin/bash

# Slurm Parameters
#SBATCH -p medium
#SBATCH --time 20:00:00
#SBATCH --job-name=cuttag_gz
#SBATCH -o cuttag_gz_%A.log
#SBATCH -e cuttag_gz_%A.err
#SBATCH -c 10
#SBATCH --mem 100G
#SBATCH -C scratch
#SBATCH --array=0

#### original pipeline and tutorial available at https://yezhengstat.github.io/CUTTag_tutorial

# throw an error if -name argument is not provided
if [ $# -eq 0 ]; then
	echo "No argument is provided"
	echo "call the script with sbatch CutandTag.sh sampleName ReadsName"
	exit 1 # general erro; if sucessful type exit 0
fi

#### set processor number
export cpu_per_task=20

bt2_ref_genome="/scratch/users/txie/data/hg38.bt/hg38" #depending on data should be changed.
spikeInRef="/scratch/users/txie/data/ecoil/ecoil"
chromSize="/scratch/users/txie/data/hg38/hg38.chrom.sizes" #depending on data should be changed. 

# set the base name of files without the extension
#readsName=$(find . -type f -name "*_merged_R1_001.fastq.gz" -printf '%f\n' | sed -e 's/merged_R1_001.fastq.gz//') ##merging, may not need this
sampleName=$1
readsName=$2

#### load modules
#module load cutadapt/2.3
module load fastqc/0.11.4
#module load bowtie/2.3.4.1
#module load samtools/1.9
module load bedtools
module load r/4.0.3
module load deeptools
######################### setting the directories
mkdir -p $(pwd)/$sampleName"_results" && resultPath=$_
mkdir -p ${resultPath}/Merged_fastq_files && fastqPath=$_
mkdir -p ${fastqPath}/FastQC && fastQCPath=$_
mkdir -p ${resultPath}/Aligned_files && alignedPath=$_
mkdir -p ${alignedPath}/Bam_files
mkdir -p ${alignedPath}/Bed_files
mkdir -p ${alignedPath}/BedGraph_files
mkdir -p ${alignedPath}/Bigwig_scaled
mkdir -p ${alignedPath}/picard_summary && picardresult=$_
mkdir -p ${resultPath}/peakCalling && PeakPath=$_

############################# Pre-processing steps

##### If the samples are sequenced in-house they will be split into sequencing lanes
##### Pulling the lanes or techical (!) replicates together with:
echo 'Step 1 : Merging the reads from separate lanes'
echo 'Step 1 : File name of processed reads is' $readsName

# set the path for merged files
#read1_uncut=${fastqPath}/${sampleName}"_merged_R1_001.fastq.gz"
#read2_uncut=${fastqPath}/${sampleName}"_merged_R2_001.fastq.gz"
read1_uncut=/home/uni08/txie/analysis/CP/07Chips/cuttag/raw/${sampleName}".R1.fastq.gz"
read2_uncut=/home/uni08/txie/analysis/CP/07Chips/cuttag/raw/${sampleName}".R2.fastq.gz"
# check if merged reads exist from previous run to skip this step
if [[ -f "$read1_uncut" && -f "$read2_uncut" ]]; then
    echo 'Step 1 : Merged reads already exist, moving to Step 2.'
else
    # concatenate separate lanes and output to merged files folder
    cat /home/uni08/txie/rawdata/cuttag/*${readsName}*_1.fq.gz > $read1_uncut
    cat /home/uni08/txie/rawdata/cuttag/*${readsName}*_2.fq.gz > $read2_uncut
fi

# check if merging is done
if [[ -f "$read1_uncut" && -f "$read2_uncut" ]]; then
    echo 'Files from separate lanes are merged into' 
    echo $read1_uncut 'and' 
    echo $read2_uncut
    echo 'Step 1 : Done!'
else
    echo 'Step 1 : Failed.'
    exit 1
fi

############################# Quality check with FASTQC

fastqc1=${fastQCPath}/$sampleName"_trimmed_R1_fastqc.html"
fastqc2=${fastQCPath}/$sampleName"_trimmed_R2_fastqc.html"

read1=$read1_uncut
read2=$read2_uncut


echo 'Step 3 : fastQC on trimmed reads'

if [[ -f "$fastqc1" && -f "$fastqc2" ]]; then
    echo 'Step 3 : FastQC reports already exist, moving to Step 4.'
else
    fastqc -o ${fastQCPath}/ -t $cpu_per_task -f fastq $read1
    fastqc -o ${fastQCPath}/ -t $cpu_per_task -f fastq $read2
fi

# check if fastqc was done 

if [[ -f "$fastqc1" && -f "$fastqc2" ]]; then
    echo 'Step 3 : Done!'
else
    echo 'Step 3 : Failed, but moving on.'
fi

############################## Alignment of reads
echo 'Step 4 : Read alignment to Reference genome as sample and E.Coli as spike-in'

bamFile=${alignedPath}/Bam_files/$sampleName"_bowtie2.bam"
spikeFile=${alignedPath}/Bam_files/$sampleName"_bowtie2_spikeIn.seqDepth"

if [[ -f "$bamFile" &&  -f "$spikeFile" ]]; then
    echo 'Step 4 : Reads are already mapped, moving on to Step 5.'
else
    ## Map to Reference genome: trimmed reads
#    bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p $cpu_per_task -x $bt2_ref_genome -1 $read1 -2 $read2 -S ${alignedPath}/$sampleName"_bowtie2.sam"  &> ${alignedPath}/Bam_files/$sampleName"_bowtie2.txt"

	## Map to Reference genome: untrimmed reads
	bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p $cpu_per_task -x $bt2_ref_genome -1 $read1 -2 $read2 -S ${alignedPath}/$sampleName"_bowtie2.sam"  &> ${alignedPath}/Bam_files/$sampleName"_bowtie2.txt"

    samtools view -bS -F 0x04 ${alignedPath}/$sampleName"_bowtie2.sam" > $bamFile
    ## Extract the 9th column from the alignment sam file which is the fragment length
    samtools view -F 0x04 ${alignedPath}/$sampleName"_bowtie2.sam" | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $sampleName/2}' > ${alignedPath}/Bam_files/$sampleName"_fragmentLen.txt"

    ## Align the reads to e coli genome : untrimmed reads
    bowtie2 --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10 -X 700 -p $cpu_per_task -x $spikeInRef -1 $read1 -2 $read2 -S ${alignedPath}/$sampleName"_bowtie2_spikeIn.sam" &> ${alignedPath}/Bam_files/$sampleName"_bowtie2_spikeIn.txt"

	## Align the reads to e coli genome : trimmed reads
#	bowtie2 --local --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10 -X 700 -p $cpu_per_task -x $spikeInRef -1 $read1 -2 $read2 -S ${alignedPath}/$sampleName"_bowtie2_spikeIn.sam" &> ${alignedPath}/Bam_files/$sampleName"_bowtie2_spikeIn.txt"

    ## extract sequencing info from e.coli
    seqDepthDouble=`samtools view -F 0x04 ${alignedPath}/$sampleName"_bowtie2_spikeIn.sam" | wc -l`
    seqDepth=$((seqDepthDouble/2))
    echo $seqDepth > ${alignedPath}/Bam_files/$sampleName"_bowtie2_spikeIn.seqDepth"
fi
############################## MQ filter and remove duplication (PCR and not PCR duplicates) 

echo 'Step 5 : removing PCR duplication and low mapping quality'

bamFileFilter=${alignedPath}/Bam_files/$sampleName"_bowtie2.filter.bam"
bamFileFilterSorted=${alignedPath}/Bam_files/$sampleName"_bowtie2.filter.sorted.bam"
minQualityScore=30
picardCMD="java  -jar -Djava.io.tmpdir=${resultPath}/tmp /usr/product/bioinfo/PICARD/2.20.2/picard.jar   "

if [[ -f "$bamFileFilterSorted" ]]; then
	echo "Step 5 : removing Duplication and filtering mapping quality have done."
else
	## removing duplication, pcr and not PCR.
	### sorting
	$picardCMD SortSam I=${bamFile} O=${alignedPath}/$sampleName"_bowtie2.sorted.bam" SORT_ORDER=coordinate # sam format seems doesn't work.
	### rmdup
	$picardCMD MarkDuplicates I=${alignedPath}/$sampleName"_bowtie2.sorted.bam" O=${alignedPath}/$sampleName"_bowtie2.rmdup.bam" REMOVE_DUPLICATES=true METRICS_FILE=${picardresult}/$sampleName"_picard.txt" TMP_DIR=${resultPath}/tmp
	samtools view -h -q $minQualityScore ${alignedPath}/$sampleName"_bowtie2.rmdup.bam"  > $bamFileFilter
	samtools sort -n -o $bamFileFilterSorted -O BAM --threads $cpu_per_task $bamFileFilter
fi

#### check if done and remove .sam files
if [[ -f "$bamFile" && -f "$spikeFile" && -f "$bamFileFilter" && -f $bamFileFilterSorted ]]; then
	rm -rf ${alignedPath}/$sampleName"_bowtie2.sam"
	rm -rf ${alignedPath}/$sampleName"_bowtie2_spikeIn.sam"
	rm -rf ${alignedPath}/$sampleName"_bowtie2.rmdup.bam"
	rm -rf ${alignedPath}/$sampleName"_bowtie2.sorted.bam"
    echo 'Step 5 : Done!'
else
    echo 'Step 5 : Failed for Reference genome or E.Coli'
    exit 1
fi

################################## File conversion
echo 'Step 6 : Filtering and conversion to bed and bedgraph format'

bedFile=${alignedPath}/Bed_files/$sampleName".fragments.bed"
bedgraphFile=${alignedPath}/BedGraph_files/$sampleName".fragments.scaled.bedgraph"
#bedgraphFile=${alignedPath}/BedGraph_files/$sampleName".fragments.raw.bedgraph"
if [[ -f "$bedFile" && -f "$bedgraphFile" ]]; then
    echo 'Step 6 : bed and bedgraph files are found, moving to Step 6.'
else
	##keep MQ30 and rumdup
    bedtools bamtobed -i $bamFileFilterSorted  -bedpe > ${alignedPath}/Bed_files/$sampleName".bed"
    #### Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
    awk '$1==$4 && $6-$2 < 1000 {print $0}' ${alignedPath}/Bed_files/$sampleName".bed" > ${alignedPath}/Bed_files/$sampleName".clean.bed"
    #### Only extract the fragment related columns
    cut -f 1,2,6 ${alignedPath}/Bed_files/$sampleName".clean.bed" | sort -k1,1 -k2,2n -k3,3n  > ${alignedPath}/Bed_files/$sampleName".fragments.bed"
    #### Create bedgraph file for peak calling
    bedtools genomecov -bg -i ${alignedPath}/Bed_files/$sampleName".fragments.bed" -g $chromSize  > ${alignedPath}/BedGraph_files/$sampleName".fragments.raw.bedgraph"

    scale_factor=`echo "10000 / $seqDepth" | bc -l`
    echo "Scaling factor for $sampleName is: $scale_factor!"
    bedtools genomecov -bg -scale $scale_factor -i ${alignedPath}/Bed_files/$sampleName".fragments.bed" -g $chromSize > ${alignedPath}/BedGraph_files/$sampleName".fragments.scaled.bedgraph"
fi

if [[ -f "$bedFile" && -f "$bedgraphFile" ]]; then
#if [[ -f "$bedFile" ]]; then
    echo 'Step 6 : Done!'
else
    echo 'Step 6 : Failed.'
    exit 1
fi


################################## SEACR peak calling
echo 'Step 7 : Peak calling with SEACR'

bash /home/uni08/txie/scripts/07Chips/cuttag/SEACR_1.3.sh $bedgraphFile 1 non stringent $PeakPath/$sampleName"_seacr_all.peaks"
bash /home/uni08/txie/scripts/07Chips/cuttag/SEACR_1.3.sh $bedgraphFile 0.01 norm stringent $PeakPath/$sampleName"_seacr_top0.01.peaks"
bash /home/uni08/txie/scripts/07Chips/cuttag/SEACR_1.3.sh $bedgraphFile 0.02 norm stringent $PeakPath/$sampleName"_seacr_top0.02.peaks"
bash /home/uni08/txie/scripts/07Chips/cuttag/SEACR_1.3.sh $bedgraphFile 0.03 norm stringent $PeakPath/$sampleName"_seacr_top0.03.peaks"


peakFile=$PeakPath/$sampleName"_seacr_top0.01.peaks.stringent.bed"

if [ -f "$peakFile" ]; then
    echo 'Step 7 : Done!'
else
    echo 'Step 7 : Failed.'
    exit 1
fi

################################## generate scaled bigwig file for IGV
echo 'Step 8 : BigWig file generate'
bigwigFile=${alignedPath}/Bigwig_scaled/$sampleName".bw"

if [[ -f "$bigwigFile" ]]; then
	echo 'Step 8 : scaled bigwig files are found.'
else
	samtools sort -o ${alignedPath}/Bigwig_scaled/$sampleName"_bowtie2.filter.sorted.sorted.bam" $bamFileFilterSorted
	samtools index ${alignedPath}/Bigwig_scaled/$sampleName"_bowtie2.filter.sorted.sorted.bam"
	bamCoverage -b ${alignedPath}/Bigwig_scaled/$sampleName"_bowtie2.filter.sorted.sorted.bam"  -o $bigwigFile
	bamCoverage --normalizeUsing RPGC  --effectiveGenomeSize 2913022398 -b ${alignedPath}/Bigwig_scaled/$sampleName"_bowtie2.filter.sorted.sorted.bam"  -o ${bigwigFile}_RPGC.bw

fi
echo 'All done!'
echo "$sampleName is done"


