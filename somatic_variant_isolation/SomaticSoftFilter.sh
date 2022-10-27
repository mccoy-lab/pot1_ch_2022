#!/bin/bash
ml bcftools
ml gatk

#Output from VQSR.sh was merged bcftools merge

SAMPLE=`echo "${1}" | awk -F"." '{print $1}'` #Create sample ID from input VCF file, e.g., "A"
OUTBASE=`echo "${1}" | sed 's/\.vcf\.gz//'` #Create an output base ID from input VCF file

TGP_UNRELIABLE=/path/to/1kgp_masked_grch38.unreliablesites_mask.bed
GNOMAD_AF_FILTER=/path/to/af-only-gnomad.hg38.AF001.hg38.vcf.gz
GERMLINE=/path/to/SampleA_Bulk.matched_normal.hetsites.vcf.gz

bcftools view ${1} |\
bcftools filter \
	--threads 4 \
	--exclude "COUNT(INFO/AC) > 1" \
	-s MultiAllelic -m + |\
bcftools filter \
	--threads 4 \
	--mask-file ${TGP_UNRELIABLE} \
	-s UnreliableSite -m + |\
bcftools filter \
	--threads 4 \
	--mask-file ${GNOMAD_AF_FILTER} \
	-s CommonPopVariant -m + |\
bcftools filter \
	--threads 4 \
	--mask-file ${GERMLINE} \
	-s GermlineHet -m + |\
bcftools filter \
	--threads 4 \
	--exclude "F_MISSING > 0.2" \
	-s MissingData -m + |\
bcftools filter \
	--threads 4 \
	--exclude 'COUNT(FORMAT/GT == "het" & FORMAT/DP >= 11) < COUNT(FORMAT/GT == "het")*0.2' \
	-s DepthFilter -m + |\
bcftools filter \
	--threads 4 \
	--exclude "INFO/AF == 1" \
	-s HomALT -m + |\
bcftools filter \
	--threads 4 \
	--exclude 'COUNT(FORMAT/GT == "het") == 0' \
	-s NoHET -m + |\
bcftools filter \
	--threads 4 \
	--exclude 'COUNT(FORMAT/GT == "hom") < COUNT(FORMAT/GT != "mis")*0.2' \
	-s HighHet -m + |\
bcftools filter \
	--threads 4 \
	--exclude 'COUNT(FORMAT/GT == "het" & FORMAT/DP >= 11 & sMIN(FORMAT/AD)/FORMAT/DP >= 0.33) < COUNT(FORMAT/GT == "het" & FORMAT/DP >= 11)*0.5' \
	-s AlleleBalance -m + \
	--output-type z \
	--output ${OUTBASE}.tmp.vcf.gz
ClusterFilter.R ${OUTBASE}.tmp.vcf.gz ${OUTBASE}.softfiltered.vcf.gz
gzip -d ${OUTBASE}.softfiltered.vcf.gz
bgzip ${OUTBASE}.softfiltered.vcf
bcftools index --threads 4 --tbi ${OUTBASE}.softfiltered.vcf.gz

rm ${OUTBASE}.tmp.vcf.gz

