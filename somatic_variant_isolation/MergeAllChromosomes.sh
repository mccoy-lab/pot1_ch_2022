#!/bin/bash
ml gatk

for i in {1..24}
do
	if (( $i == 23 ));then
		i=X
	fi
	if (( $i == 24 ));then
		i=Y
	fi

	#Create a variable which contains a new-line delimited list of all chromosome vcf files output from HaplotypeCaller_wrapper
	gVCFs=`ls *.chr${i}.vcf.gz | sed 's/^/--variant /' | tr '\n' ' '`

	gatk CombineGVCFs -R ~/path/to/GRCh38.d1.vd1.fa \
		`echo "${gVCFs}"` \
		--output Merged_Chromosomes_Full/All_samples.chr${i}.g.vcf.gz
done