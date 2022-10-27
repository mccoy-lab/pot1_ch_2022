#!/bin/bash
ml gatk

#This command is written specifically to process all genotype *vcf.gz files for a single individual - in this case "A"

#All resources obtained from gs://gcp-public-data--broad-references/hg38/v0 as presented in https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels 
MILLS=/path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
AXIOM=/path/to/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
DBSNP=/path/to/Homo_sapiens_assembly38.dbsnp138.vcf
TGP=/hpath/to/1000G_phase1.snps.high_confidence.hg38.vcf.gz
HAPMAP=/path/to/hapmap_3.3.hg38.vcf.gz
OMNI=/path/to/1000G_omni2.5.hg38.vcf.gz

for i in {1..24}
do

	#Obtain base name
	if (( $i == 23 ));then
		i=X
	fi
	if (( $i == 24 ));then
		i=Y
	fi

	BASENAME=`echo "A.chr${i}.genotype.vcf.gz" | sed 's/\.vcf\.gz//'`

	#Remove sample information
	gatk MakeSitesOnlyVcf \
		-I A.chr${i}.genotype.vcf.gz \
		-O ${BASENAME}.sites_only.vcf.gz

	if [[ ! -f ${BASENAME}.sites_only.vcf.gz ]]; then
		echo "ERROR: MakeSitesOnlyVcf output not found."
		exit
	fi

	#Indel filter
	gatk --java-options "-Xmx64g" VariantRecalibrator \
		-V ${BASENAME}.sites_only.vcf.gz \
		--trust-all-polymorphic \
		-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
		-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP -an BaseQRankSum \
		-mode INDEL \
		--max-gaussians 4 \
		--resource:mills,known=false,training=true,truth=true,prior=12 ${MILLS} \
	 	--resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${AXIOM} \
	 	--resource:dbsnp,known=true,training=false,truth=false,prior=2 ${DBSNP} \
	 	-O ${BASENAME}.indels.recal \
	 	--tranches-file ${BASENAME}.indels.tranches

	 #SNP filter
	 gatk --java-options "-Xmx64g" VariantRecalibrator \
		-V ${BASENAME}.sites_only.vcf.gz \
		--trust-all-polymorphic \
		-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
		-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP -an BaseQRankSum \
		-mode SNP \
		--max-gaussians 6 \
	    --resource:hapmap,known=false,training=true,truth=true,prior=15 ${HAPMAP} \
	    --resource:omni,known=false,training=true,truth=true,prior=12 ${OMNI} \
	    --resource:1000G,known=false,training=true,truth=false,prior=10 ${TGP} \
	    --resource:dbsnp,known=true,training=false,truth=false,prior=7 ${DBSNP} \
	 	-O ${BASENAME}.snps.recal \
	 	--tranches-file ${BASENAME}.snps.tranches

	if [[ ! -f ${BASENAME}.indels.recal ]]; then
		echo "ERROR: VariantRecalibrator Indel output not found."
		exit
	fi

	if [[ ! -f ${BASENAME}.snps.recal ]]; then
		echo "ERROR: VariantRecalibrator Snps output not found."
		exit
	fi

	#Apply Indel VQSR
	gatk --java-options "-Xmx64g -Xms64g" ApplyVQSR \
		-V A.chr${i}.genotype.vcf.gz \
		--recal-file ${BASENAME}.indels.recal \
		--tranches-file ${BASENAME}.indels.tranches \
		--truth-sensitivity-filter-level 99.7 \
		--create-output-variant-index true \
		-mode INDEL \
		-O ${BASENAME}.indel.recalibrated.vcf.gz

	#Apply SNP VQSR
	gatk --java-options "-Xmx64g -Xms64g" ApplyVQSR \
		-V ${BASENAME}.indel.recalibrated.vcf.gz \
		--recal-file ${BASENAME}.snps.recal \
		--tranches-file ${BASENAME}.snps.tranches \
		--truth-sensitivity-filter-level 99.7 \
		--create-output-variant-index true \
		-mode SNP \
		-O ${BASENAME}.VQSR.vcf.gz

	rm ${BASENAME}.sites_only.vcf.gz
	rm ${BASENAME}.indel.recalibrated.vcf.gz 
	rm ${BASENAME}.indel.recalibrated.vcf.gz.tbi
done