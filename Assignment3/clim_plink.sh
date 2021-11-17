#!/bin/bash

CHR="Chr04"

echo $CHR  # Does it look right?

cd /data/project_data/PopGenomics/shared/${CHR}_sweeps/LAI/

# Get phenotype data from meta file for the admixed samples:

tail -n +2 /data/project_data/PopGenomics/climDat.txt | \
grep -w -f Admixed.Inds - >clim.pheno

printf "#FID\tIID\tmean_finalFreeze\tmean_cGDDfreeze\tmed_DD0\n%.0s" >Plink/assignment3/pheno.forPlink

cat clim.pheno >>Plink/assignment3/pheno.forPlink

#larger positive number - more infected relative to the mean
# Get K=2 ADMIX to use as covariate; Need to use tail -n +2 to skip the header in the metadata file before pasting to the Q file

tail -n +2 /data/project_data/PopGenomics/climDat.txt | \
cut -f1 | \
paste - /data/project_data/PopGenomics/shared/poplar_hybrids.LDpruned.5.Q | \
grep -w -f Admixed.Inds - | \
sed 's/\s/\t/g' | \
cut -f 1,5 >Plink/assignment3/Admixed_KBals

cat /data/project_data/PopGenomics/climDat.txt | \
cut -f1-2 | \
grep -w -f Admixed.Inds - | \
sed 's/\s/\t/g' | \
cut -f2 >Plink/assignment3/Transect

# Create the cov file with KBals

printf "#FID\tIID\tK2Q\n%.0s" >Plink/assignment3/cov.forPlink
paste Plink/assignment3/Transect Plink/assignment3/Admixed_KBals >>Plink/assignment3/cov.forPlink

#run plink
cd Plink/assignment3/

plink2 --bfile ${CHR}_Admixed_FAI \
--pheno pheno.forPlink \
--covar cov.forPlink \
--glm omit-ref
