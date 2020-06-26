#!/usr/bin/env bash

while read CHROM LINK
do
    bcftools view -R targets.bed -S samples.txt --force-samples $LINK -Oz > $CHROM.vcf.gz
done < remotes.tsv

rm *.tbi

for i in *.gz
do 
    tabix $i
done

#bcftools concat *.vcf.gz -Ob > test_samples.bcf
#bcftools index test_samples.bcf

#rm *.vcf.g*
