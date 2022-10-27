#!/bin/bash
ml gatk

while read BAM #Loop through each line (path to output from bwa.sh) of bam_list.txt
do
	FIRST_HEADER=`samtools view ${BAM} | head -1 | awk '{print $1}'` #Extract data from first read alignment header
	OUTPUT=`echo "${BAM}" | awk -F"/" '{print $(NF)}' | sed 's/\.bam/\.RG.bam/'` #Set output name variable

	RG_SAMPLE=`echo "${BAM}" | awk -F"/" '{print $(NF)}' | awk -F"." '{print $1}' | awk -F"_" '{print $1}' | sed 's/JH189//'` #Extract sample ID from ${BAM}
	RG_LIBRARY=`echo "${BAM}" | awk -F"/" '{print $(NF)}' | awk -F"." '{print $1}'` #Extract library ID from ${BAM}
	RG_PLATFORM=illumina #Set sequencing platform variable
	RG_FLOWCELL=`echo "${FIRST_HEADER}" | awk -F":" '{print $3}'` #Extract flowcell ID from first header 

	#stdout sanity-check	
	echo "File: ${BAM}"
	echo "Output: Mapped_libraries_with_readgroups/${OUTPUT}"
	echo "Sample: ${RG_SAMPLE}"
	echo "Library: ${RG_LIBRARY}"
	echo "Platform: ${RG_PLATFORM}"
	echo "Flowcell: ${RG_FLOWCELL}"
	echo ""

	#Add readgroup information to each BAM
	gatk AddOrReplaceReadGroups -I ${BAM} \
		-O  Mapped_libraries_with_readgroups/${OUTPUT} \
		-SO coordinate \
		-LB ${RG_LIBRARY} \
		-PL ${RG_PLATFORM} \
		-PU ${RG_FLOWCELL} \
		-SM ${RG_SAMPLE} \
		--CREATE_INDEX
	
done < bam_list.txt #List of bam output from bwa.sh with path (e.g., "/path/to/library.merged.sort.markdup.bam")
