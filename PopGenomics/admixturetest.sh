#!/bin/bash

# Set the path and filename to the vcf file filtered for minor allele frequency (maf)

VCF=/data/project_data/PopGenomics/poplar_hybrids.maf05.vcf.gz

cd ~/myresults


plink2 --vcf $VCF \
--threads 3 \
--double-id \ 
--allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \ 
--out poplar_hybrids_noLD


mkdir Admixture

FILE=poplar_hybrids

plink2 --vcf $VCF \
--threads 3 \
--allow-extra-chr \
--make-bed \
--out Admixture/$FILE 

plink2 --bfile Admixture/$FILE \
--threads 3 \
--allow-extra-chr \
--set-missing-var-ids @:# \
--extract poplar_hybrids_noLD.prune.in \
--make-bed \
--out Admixture/$FILE.LDpruned 


cd Admixture

FILE2=poplar_hybrids.LDpruned

awk '{$1=0;print $0}' $FILE2.bim > $FILE2.bim.tmp
mv $FILE2.bim.tmp $FILE2.bim


K=3

admixture -j3 --cv $FILE2.bed $K >log${K}.out