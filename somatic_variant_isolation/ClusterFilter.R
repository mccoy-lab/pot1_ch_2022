#!/usr/bin/env Rscript

#Read input VCF from arg1; read output name from arg2
args=commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("\n>>>Syntax: ClusterFilter.R [in.vcf.gz] [out.vcf.gz]", call.=FALSE)
} 

#Load required libraries
library(GenomicRanges,quietly=TRUE)
library(stringr,quietly=TRUE)
library(readr,quietly=TRUE)
library(dplyr,quietly=TRUE)
library(tidyr,quietly=TRUE)
library(data.table,quietly=TRUE)
library(vcfR,quietly=TRUE)

in_vcf<-read.vcfR(args[1])
out_vcf<-args[2]

chrom <- unlist(in_vcf@fix[, 1]) #create chromosome name vector from input vcf
pos <- as.numeric(unlist(in_vcf@fix[, 2])) #create variant position vector from input vcf
grange <- makeGRangesFromDataFrame(data.table(chr = chrom, start = pos, end = pos)) #Create a grange object
dists <- distanceToNearest(grange)@elementMetadata$distance + 1 #Compute distance to nearest variant 

filter_field<-in_vcf@fix[,"FILTER"] #Obtain filter field for each variant
filter_field[dists <= 10]<-paste(filter_field[dists <= 10],";ClusteredEvents",sep="") %>% str_replace("PASS;","") #If within 10bp of variant, add "ClusteredEvents" to filter; remove "PASS" if present
in_vcf@fix[,"FILTER"]<-filter_field #Assign new filter field to FILTER slot

in_vcf@meta<-append(in_vcf@meta,"##FILTER=<ID=ClusteredEvents,Description=\"Variant within 10bp of another variant\">") #Add info for filter information to meta slot

write.vcf(in_vcf,file=out_vcf) #Outputs in gzip; must be converted to bgzip to be compatible with bcftools indexing