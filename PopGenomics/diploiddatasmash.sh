#!/bin/bash

## From within your LAI/ directory:

CHR="Chr04"

echo $CHR # Is it right?

Nsites=`tail -n +2 ${CHR}.kept.sites | wc -l | sed 's/\s/\t/' | cut -f1` # calculates and stores the number of SNP sites for your chromosome

echo $Nsites # For Chr02, Nsites=281887 SNPs

Ninds=`wc -l Admixed.Inds | sed 's/\s/\t/' | cut -f1` # calculates and stores the number of admixed individuals you previously identified

echo $Ninds  # Should be 442 individuals

touch ${CHR}_matrix.diploid

for file in Loter_out/*.txt
do
datamash --field-separator=" " sum 1-${Nsites} <$file >>${CHR}_matrix.diploid
done


#change into plink format
seq -f "snp%02g" 1 $Nsites >sites

printf 'A\n%.0s' $(seq $Nsites) >allele1  # Create a dummy column of 'A' the length of your Nsites file

printf "T\n%.0s" $(seq $Nsites) >allele2 # Create a dummy column of 'T' the length of your Nsites file

mkdir Plink

paste sites allele1 allele2 ${CHR}_matrix.diploid.tr >Plink/${CHR}_matrix.diploid.tr.forPlink

#fam file for our samples - using transect as family
cat /data/project_data/PopGenomics/Combined_Transect_Sampling_Data_2020.txt | \
cut -f1-2 | \
grep -w -f Admixed.Inds - | \
cut -f2 | \
paste - Admixed.Inds >FID_IID

printf '0\t0\t0\t-9\n%.0s' $(seq $Ninds) >dummy

paste FID_IID dummy  >Plink/${CHR}_fam.forPlink


cd Plink/ 

plink2 --import-dosage ${CHR}_matrix.diploid.tr.forPlink noheader \
--fam ${CHR}_fam.forPlink \
--make-bed \
--out ${CHR}_Admixed_FAI

plink2 --bfile ${CHR}_Admixed_FAI --freq --out ${CHR}_LAI_freq
