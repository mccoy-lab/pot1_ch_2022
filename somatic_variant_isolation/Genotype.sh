#!/bin/bash

ml bcftools
ml gatk

SAMPLE_LIST=/path/to/Individual_Samples.txt
GRCh=/path/to/GRCh38.d1.vd1.fa

for i in {1..24}
do
	if (( $i == 23 ));then
		i=X
	fi
	if (( $i == 24 ));then
		i=Y
	fi

	#Subset full gVCF file to focal samples
	bcftools view --threads 8 \
		--with-header \
		-O z \
		-o A.chr${i}.g.vcf.gz \
		--samples-file ${SAMPLE_LIST} \
		All_samples.chr${i}.g.vcf.gz

	#Index subsetted gVCF file
	bcftools index --threads 8 \
		-t \
		A.chr${i}.g.vcf.gz

	#Call genotypes for every sample in gVCF
	gatk GenotypeGVCFs \
		--output Genotypes/A.chr${i}.genotype.vcf.gz \
		--variant A.chr${i}.g.vcf.gz \
		--reference ${GRCh}
done

################ ERRATA ##################
# Example structure of $SAMPLE_LIST:
#
# >>> ls ${SAMPLE_LIST}_A
# SAMPLEA1
# SAMPLEA2
# SAMPLEA3
# SAMBLEA4
#
# Example of a different ${SAMPLE_LIST}
# >>> ls ${SAMPLE_LIST}_B
# SAMPLEB1
# SAMPLEB2
# SAMPLEB3
# SAMBLEB4
#
# These samples should be present in the 
# header each gVCF; assigned from read-
# group information during variant calling
##########################################