#!/bin/bash

LIB_DIR=/path/to/Mapped_libraries_with_readgroups
REF=/path/to/GRCh38.d1.vd1.fa

while read CLONE
do
	for i in {1..24} #Set i to chromosome number; if i is 23/24, set to X or Y, respectively.
	do
		if [[ $i == 23 ]]; then
			i="X"
		fi
		if [[ $i == 24 ]]; then
			i="Y"
		fi

		#Run haplotype caller in GVCF mode per chromosome (originally run as an array job with for every chromosome in parallel)
		gatk HaplotypeCaller --emit-ref-confidence GVCF \
			--output ${CLONE}.merged.sort.markdup.RG.g.chr$i.vcf.gz \
			--input ${LIB_DIR}/${CLONE}.merged.sort.markdup.RG.bam \
			--reference ${REF} \
			--intervals chr${i}
	done
done < lib_list.txt #List of library ID's

################ ERRATA ##################
# Example structure of $LIB_DIR:
#
# >>> ls ${LIB_DIR}
# SAMPLEA1_S10.merged.sort.markdup.RG.bam
# SAMPLEA1_S10.merged.sort.markdup.RG.bai
# SAMPLEB1_S3.merged.sort.markdup.RG.bam
# SAMPLEB1_S3.merged.sort.markdup.RG.bai
#
# Example structure of lib_list.txt:
#
# >>> cat lib_list.txt
# SAMPLEA1_S10
# SAMPLEB1_S3
##########################################