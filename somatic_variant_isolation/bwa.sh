#!/bin/bash

module load bwa-mem
module load samtools

THREADS=16 #Threads for bwa & samtools
LIBRARY_PATH=/path/to/clonal_hematopoiesis #Path to all fastqs
OUTPUT_PATH=/path/to/Mapped_libraries #Path to output directory
REFERENCE_PATH=/path/to/GRCh38.d1.vd1.fa #Path to reference genome (with associated indices for BWA & GATK)

echo "LIBRARY DIRECTORY: ${LIBRARY_PATH}"
echo "OUTPUT DIRECTORY: ${OUTPUT_PATH}"
echo "REFERENCE PATH: ${REFERENCE_PATH}"

while read COLONY_ID #Loop through names in uniq_libraries.txt
do
	if [[ `ls ${LIBRARY_PATH}/${COLONY_ID}*R1_001.fastq.gz | wc -l` == 4 ]]; then #If sample was run on four lanes
		FWD=`echo -n "'<zcat ${LIBRARY_PATH}/${COLONY_ID}_L001_R1_001.fastq.gz ${LIBRARY_PATH}/${COLONY_ID}_L002_R1_001.fastq.gz ${LIBRARY_PATH}/${COLONY_ID}_L003_R1_001.fastq.gz ${LIBRARY_PATH}/${COLONY_ID}_L004_R1_001.fastq.gz'"`
		REV=`echo -n "'<zcat ${LIBRARY_PATH}/${COLONY_ID}_L001_R2_001.fastq.gz ${LIBRARY_PATH}/${COLONY_ID}_L002_R2_001.fastq.gz ${LIBRARY_PATH}/${COLONY_ID}_L003_R2_001.fastq.gz ${LIBRARY_PATH}/${COLONY_ID}_L004_R2_001.fastq.gz'"`
	fi
	if [[ `ls ${LIBRARY_PATH}/${COLONY_ID}*R1_001.fastq.gz | wc -l` == 2 ]]; then #If sample was run on two lanes
		FWD=`echo -n "'<zcat ${LIBRARY_PATH}/${COLONY_ID}_L001_R1_001.fastq.gz ${LIBRARY_PATH}/${COLONY_ID}_L002_R1_001.fastq.gz'"`
		REV=`echo -n "'<zcat ${LIBRARY_PATH}/${COLONY_ID}_L001_R2_001.fastq.gz ${LIBRARY_PATH}/${COLONY_ID}_L002_R2_001.fastq.gz'"`
	fi

	echo "${FWD}"
	echo "${REV}"

	CMD=`echo "bwa mem -t $THREADS ${REFERENCE_PATH} ${FWD} ${REV}"`
	
	#Run BWA & convert SAM-->BAM
	eval ${CMD} | samtools view -Shb -@ $THREADS -o ${OUTPUT_PATH}/${COLONY_ID}.merged.bam
	
	#Add read-mate coordinates to BAM file
	samtools fixmate -m \
		-@ $THREADS ${OUTPUT_PATH}/${COLONY_ID}.merged.bam \
		${OUTPUT_PATH}/${COLONY_ID}.merged.fm.bam
	
	#Coordinate sort BAM file
	samtools sort -@ $THREADS -m 4G \
		-o ${OUTPUT_PATH}/${COLONY_ID}.merged.sort.bam \
		${OUTPUT_PATH}/${COLONY_ID}.merged.fm.bam
	
	#Flag PCR duplicates
	samtools markdup -@ $THREADS \
		-f ${OUTPUT_PATH}/${COLONY_ID}.merged.sort.markdup_metrics.txt \
		-d 2500 \
		${OUTPUT_PATH}/${COLONY_ID}.merged.sort.bam \
		${OUTPUT_PATH}/${COLONY_ID}.merged.sort.markdup.bam
	
	#Index BAM file
	samtools index -@ $THREADS ${OUTPUT_PATH}/${COLONY_ID}.merged.sort.markdup.bam 
	
	#Generate read-mapping stats
	samtools stats -@ $THREADS \
		-r ${REFERENCE_PATH} \
		--reference ${REFERENCE_PATH} \
		${OUTPUT_PATH}/${COLONY_ID}.merged.sort.markdup.bam > ${OUTPUT_PATH}/${COLONY_ID}.merged.sort.markdup.bam.stats

	#Clean-up intermediate files
	rm ${OUTPUT_PATH}/${COLONY_ID}.merged.bam
	rm ${OUTPUT_PATH}/${COLONY_ID}.merged.fm.bam 
	rm ${OUTPUT_PATH}/${COLONY_ID}.merged.sort.bam

done < uniq_libraries.txt #Text file containing a list unique sample identifiers


################ ERRATA ##################
# Example structure of $LIBRARY_PATH:
#
# >>> ls ${LIBRARY_PATH}
# SAMPLEA1_S10_L001_R1_001.fastq.gz
# SAMPLEA1_S10_L001_R2_001.fastq.gz
# SAMBLEB1_S3_L001_R1_001.fastq.gz
# SAMPLEB1_S3_L001_R2_001.fastq.gz

# Example structure of uniq_libraries.txt:
#
# >>> cat uniq_libraries.txt
# SAMPLEA1_S10
# SAMPLEB1_S3
##########################################