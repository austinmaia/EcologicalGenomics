#!/bin/bash

myChr=Chr04  

# myChr=Chr04  # I used Chr02 for my test run...


cd /data/project_data/PopGenomics

# Run VCFtools to subset the big vcf file for just your chromosome

vcftools --gzvcf poplar_hybrids.maf05.vcf.gz \
--chr $myChr \
--out shared/$myChr \
--recode

# Extract the centromere coordinates for your chromosome so you can exclude those regions from your sweep analysis

grep $myChr poplar_centromeres.txt > shared/${myChr}_centromere.txt # grab the centromere location for your chromosome

cd shared/

mkdir ${myChr}_sweeps  # make a new directory for your chromosome analyses

mv *${myChr}* ${myChr}_sweeps # clean up the space by moving all files into your the directory you just made

cd ${myChr}_sweeps

# Test for selective sweeps

RAiSD -n $myChr \
-I ${myChr}.recode.vcf \
-f -t -R -P -D -A 0.99 \
-X ${myChr}_centromere.txt


# First, need to subset the metadata file for just those individuals with balsamifera ancestry
# We can do this using an interactive R session at the commandline. 
# An alternative is to put these R commands in a script, save it with the ".r" extension, 
# and at the commandline type "Rscript myscript.r"

R # Opens an interactive R session within Unix...
Qscores <- read.table("../poplar_hybrids.LDpruned.5.Q", sep=" ",header=F)
names(Qscores) = c("K1","K2","K3","K4","K5")

meta <- read.table("../../Combined_Transect_Sampling_Data_2020.txt",sep="\t",header=T)

merged <- cbind(meta,Qscores)
str(merged)

Bals_Inds <- merged[which(merged$K4>0.5),1]  #Bals - all with greater than 50% K4 ancestry
length(Bals_Inds) # Should net you 188 individuals

Tricho_Inds <- merged[which(merged$K4<=0.5),1] #Tricho - everything with less than 50% K4
length(Tricho_Inds) # Should net you 388 individuals

# Write out your Bals and Tricho lists as tab-delimited text files
write.table(Bals_Inds, "Bals_Inds.txt", quote=F, row.names=F, col.names=F)

write.table(Tricho_Inds, "Tricho_Inds.txt", quote=F, row.names=F, col.names=F)

quit()

# Calculate Fst between Balsam and Tricho using sliding windows of 50kb

vcftools --vcf ${myChr}.recode.vcf \
--weir-fst-pop Bals_Inds.txt \
--weir-fst-pop Tricho_Inds.txt \
--fst-window-size 50000 \
--out Bals_Tricho_All
# When prompted with: "Save workspace image? [y/n/c]"  choose: n