#!/bin/bash
ml bcftools
ml samtools
ml emboss

#PROGRAM TAKES A SINGLE COMMAND LINE ARGUMENT: LIST OF SAMPLES CONTAINED WITHIN A VCF

REF_FASTA=/path/to/GRCh38.d1.vd1.fa #Reference genome
INDIVIDUAL=`echo "${1}" | awk -F"." '{print $1}'` #For format of {INDIVIDUAL}.clone_list.txt, extract INDIVIDUAL
NOW=`date '+%F_%H-%M-%S'` #Create unique date for tmp directory

REFBASE_REGIONS=IndividualA.full.REFbase.regions #BED format of SNP sites where germline = homREF; determined by clone genotyping sites where AF <= 0.5
ALTBASE_REGIONS=IndividualA.full.ALTbase.regions #BED format of SNP sites where germline = homALT; determined by clone genotyping sites where AF > 0.5
VCF=IndividualA.full.genotype.VQSR.PASS.vcf.gz #VCF with PASS variant sitesper clone; subset output from SomaticSoftFilter.sh

mkdir fasta_tmp_${NOW} #Create tmp directory

while read SAMPLE #Loop through clone list
do
	#Index focal sites in reference fasta, obtain consensus haplotype where het is set to ALT, merge all single site sequences under a single fasta header, rename sequence, and output to a temporary file
	samtools faidx ${REF_FASTA} -r ${REFBASE_REGIONS} |\
		bcftools consensus --haplotype A --sample ${SAMPLE} ${VCF} |\
		union -filter |\
		sed "/^>/s/>.*/>${SAMPLE}/" > fasta_tmp_${NOW}/${SAMPLE}.REFbase.fa
	#Do the same as above, but set heterozygous sites to REF (homALT -> het mutations)
	samtools faidx ${REF_FASTA} -r ${ALTBASE_REGIONS} |\
		bcftools consensus --haplotype R --sample ${SAMPLE} ${VCF} |\
		union -filter |\
		sed "/^>/s/>.*/>${SAMPLE}/" > fasta_tmp_${NOW}/${SAMPLE}.ALTbase.fa
	
	#Consolidate all fasta sequences from ${SAMPLE} (each site output from above is has a separate header) into a concatenated sequence under a single header
	cat fasta_tmp_${NOW}/${SAMPLE}.REFbase.fa fasta_tmp_${NOW}/${SAMPLE}.ALTbase.fa |\
		union -filter > fasta_tmp_${NOW}/${SAMPLE}.${INDIVIDUAL}.variants.fasta
	rm fasta_tmp_${NOW}/${SAMPLE}.REFbase.fa fasta_tmp_${NOW}/${SAMPLE}.ALTbase.fa

done < $1 #List of samples to extract from vcf

#CREATE A GERMLINE REFERENCE SEQUENCE
REFBASE_LOCATIONS=IndividualA.full.REFbase.regions.tab #Same as $REFBASE_REGIONS, but tab delimited for bcftools
ALTBASE_LOCATIONS=IndividualA.full.ALTbase.regions.tab #Same as $ALTBASE_REGIONS, but tab delimited for bcftools
VCF=IndividualA.full.genotype.VQSR.PASS.vcf.gz 

#Construct synthetic germline reference haplotype where heterozygous sites with AF <0.5 are set to the REF allele
bcftools query -f '>%CHROM\_%POS\n%REF\n' -R ${REFBASE_LOCATIONS} ${VCF} |\
	union -filter |\
	sed "/^>/s/>.*/>${INDIVIDUAL}Bulk/" > fasta_tmp_${NOW}/${INDIVIDUAL}Bulk.REFbase.fa
#Construct synthetic germline reference haplotype where heterozygous sites with AF >=0.5 are set to the ALT allele
bcftools query -f '>%CHROM\_%POS\n%ALT\n' -R ${ALTBASE_LOCATIONS} ${VCF} |\
	union -filter |\
	sed "/^>/s/>.*/>${INDIVIDUAL}Bulk/" > fasta_tmp_${NOW}/${INDIVIDUAL}Bulk.ALTbase.fa

#Consolidate all fasta sequences generated for ${INDIVIDUAL}Bulk
cat fasta_tmp_${NOW}/${INDIVIDUAL}Bulk.REFbase.fa fasta_tmp_${NOW}/${INDIVIDUAL}Bulk.ALTbase.fa |\
	union -filter > fasta_tmp_${NOW}/${INDIVIDUAL}Bulk.${INDIVIDUAL}.variants.fasta
rm fasta_tmp_${NOW}/${INDIVIDUAL}Bulk.REFbase.fa fasta_tmp_${NOW}/${INDIVIDUAL}Bulk.ALTbase.fa

#Consolidate all single sample fasta files into a single joint fasta
cat fasta_tmp_${NOW}/*.${INDIVIDUAL}.variants.fasta >> ${INDIVIDUAL}.all_samples.variant_sites.fasta
rm fasta_tmp_${NOW}/*.${INDIVIDUAL}.variants.fasta

rmdir fasta_tmp_${NOW}