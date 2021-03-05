#!/bin/env bash
#SBATCH -J SCV
#SBATCH -n 10
#SBATCH -p batch
#SBATCH --output scv2.o
#SBATCH --error scv2.e
#SBATCH --mail-user s.henson@cgiar.org
#SBATCH --mail-type ALL

module load cutadapt/v1.8.1 ## required for TrimGalore
module load fastqc/0.11.7 ## required for TrimGalore
module load ivar/1.3 ## load samtools and htslib as well
module load bwa/0.7.17
module load bedtools/2.29.0
module load bcftools/1.11

module list

## Required data files
samplelist=$1 ## File with a list of sample IDs each on a new line
BASEDIR=~/SARS-CoV/ILRI_data
ref=${BASEDIR}/../db/nCoV-2019.reference.fasta
ref_gff=${BASEDIR}/../db/MN908947.3.gff3
primers_fa=${BASEDIR}/../db/nCoV-2019.primers.fa 
primers_bed=${BASEDIR}/../db/nCoV-2019.primer.bed ## if already exists
ampl_bed=${BASEDIR}/../db/nCoV-2019.amplicons.bed

while read sample;
do
	sample_start="$(date +%s)"
	echo "##### Working on sample $sample #####"
	
	## sample files
	prefix=$sample
	R1=${BASEDIR}/fastq/${sample}_*_R1_001.fastq.gz
	R2=${BASEDIR}/fastq/${sample}_*_R2_001.fastq.gz
	outdir=${BASEDIR}/${sample}_output
	mkdir -p $outdir
	
	## Trim reads
	## Defaults:
	##		-q 20 (Trims from ends)
	##		--phre33
	##		--adapter (autodetects)
	##		--length 20 (reads shorter than this will be discarded)
	## --paired performs trimming in pair-aware manner
	## --fastqc run fastqc after trimming
	TRIM_START="$(date +%s)"
	~/soft/TrimGalore-0.6.6/trim_galore --paired -o $outdir --basename ${prefix} $R1 $R2
	TRIM_END="$(date +%s)"
	echo "##### Trim_galore ran in $((TRIM_END-TRIM_START))s."

	cd $outdir
	
	if [ ! -f "${ref}.bwt" ]
	then
		bwa index $ref
	fi
	
	echo "##### Mapping $sample reads to $ref #####"
	## exclude unmapped reads
	## exclude supplementary alignment
	## Map trimmed reads to reference
	## Write mapped (-F 4) and primary alignment (-F 2048) reads only to bam
	BWA_START="$(date +%s)"
	bwa mem -t 10 $ref ${prefix}_val_1.fq.gz ${prefix}_val_2.fq.gz | samtools view -b -F 4 -F 2048 | samtools sort -o ${prefix}.sorted.bam   
	samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
	## Get amplicon coverage
	bedtools coverage -mean -a $ampl_bed -b ${prefix}.sorted.bam > ${prefix}.sorted.amplicon.cvrg
	BWA_END="$(date +%s)"
	echo "##### BWA ran in $((BWA_END-BWA_START))s."

	## Map primers to reference
	## Uncomment below 2 lines if primer_bed doesn't exist
	#bwa mem -k 5 -T 16 $ref $primers_fa | samtools view -b -F 4 > primers.bam
	#bedtools bamtobed -i primers.bam > primers.bed

	## Get coverage of primers
	bedtools coverage -mean -a $primers_bed -b ${prefix}.sorted.bam > ${prefix}.sorted.primers.cvrg

	echo "##### Trim primers from reads using ivar #####"
	## Trim primers from reads
	## Following that, reads < Q20 (-q) are also trimmed using a sliding window of 4 (default)
	## Reads shorter than 30 bp (-m) are removed (default)
	## -e includes reads without primers. Used for Nextera. If ligation prep remove this flag
	PRIMERS_START="$(date +%s)"
	time ivar trim -e -b $primers_bed -p ${prefix}.trimmed -i ${prefix}.sorted.bam 2>&1 | tee -a primer_trim.log
	echo "##### Primer trimmed output: ${prefix}.trimmed.bam"
	PRIMERS_END="$(date +%s)"
	echo "##### Primer trimming ran in $((PRIMERS_END-PRIMERS_START))s."

	echo "##### Sort and calculate depth of trimmed reads #####"
	# Sort and index trimmed BAM file.
	samtools sort -o ${prefix}.trimmed.sorted.bam ${prefix}.trimmed.bam
	samtools index ${prefix}.trimmed.sorted.bam

	## Get the depth of mapped reads
	samtools depth -a ${prefix}.sorted.bam > ${prefix}.sorted.depth
	samtools depth -a ${prefix}.trimmed.sorted.bam > ${prefix}.trimmed.sorted.depth

	## Get genome coverage in BedGraph format for visualizations
	bedtools genomecov -bg -ibam ${prefix}.sorted.bam > ${prefix}.sorted.BedGraph
	bedtools genomecov -bg -ibam ${prefix}.trimmed.sorted.bam > ${prefix}.trimmed.sorted.BedGraph

	echo "##### Make a consensus #####"
	CONSENSUS_START="$(date +%s)"
	samtools mpileup -aa -A -d 0 -Q 0 ${prefix}.trimmed.sorted.bam | ivar consensus -t 0.75 -m 10 -n N -p ${prefix}.consensus 2>&1 | tee -a consensus_call.log
	bwa index -p ${prefix}.consensus ${prefix}.consensus.fa
	CONSENSUS_END="$(date +%s)"
	echo "##### ivar consensus calling ran in $((CONSENSUS_END-CONSENSUS_START))s."

	echo "##### Calling variants against reference #####"
	## Call iSNVs from the BAMS without reads from the masked amplicons.
	## -t Min. frequency threshold = 25%
	## -m Min. read depth = 10 
	## -q Min. quality = 20 (default)
	VAR_START="$(date +%s)"
	samtools mpileup -A -d 0 --reference $ref -Q 0 ${prefix}.trimmed.sorted.bam | ivar variants -m 10 -p ${prefix} -t 0.25 -r $ref -g $ref_gff 
	VAR_END="$(date +%s)"
	echo "##### ivar variants ran in $((VAR_END-VAR_START))s."

	echo "##### Genotyping #####" 
	python3 type_vcf.py -i $prefix -y ${BASEDIR}/../db/SARS-CoV-2.types.yaml -ov ${prefix}_typed.vcf -ot ${prefix}_typed.csv -os ${prefix}_typed.summary.csv -dp 10 -t ${prefix}.tsv $ref_gff $ref
	
	cd $BASEDIR
	sample_end="$(date +%s)"
	echo "##### $sample ran in $((sample_end-sample_start))s."
	echo "##################################################"
done < $samplelist

echo "All done!!"
