library(ggplot2)

getwd()
setwd("~/")
Qscores_new <- read.table("poplar_hybrids.LDpruned.5.Q", sep=" ",header=F) # Rename to the new file!
names(Qscores_new) <- c("K1","K2", "K3","K4","K5")  # Customize to the level of K!
meta <- read.table("Combined_Transect_Sampling_Data_2020.txt", sep="\t",header=T) 

head(Qscores_new)
# K1       K2       K3       K4       K5    ShDiv
# 1 0.099116 0.640200 0.001462 0.087846 0.171376 1.040106
# 2 0.037871 0.412992 0.004180 0.404622 0.140335 1.153771
# 3 0.080780 0.640749 0.002973 0.069082 0.206416 1.016065
# 4 0.080562 0.390881 0.000137 0.408208 0.120213 1.191723
# 5 0.157144 0.591063 0.004248 0.086425 0.161120 1.130561
# 6 0.026937 0.372067 0.000010 0.444475 0.156512 1.116006

# Calculate Shannon Diversity across the K different Q-scores per individual
  #shannon diversity = summation of (p(k) / ln(p(k)))
  #higher shannon diversity with higher k

K=5  # Change X to the level of K we're investigating

tmp=numeric()

for(i in 1:nrow(Qscores_new)){
  for(j in 1:K){
    tmp[j] = Qscores_new[i,j]*log(Qscores_new[i,j])
  }
  Qscores_new$ShDiv[i] = -1*sum(tmp)
}

het <- read.table("het_maf05.het", sep="\t",header=T)

str(het) # What's in this dataframe?  
# $ INDV   : chr  "201" "202" "203" "204" ...
# $ O.HOM. : int  3298804 2904842 3374260 3045988 3368363 3168390 3490846 3356708 3355916 3297833 ...
# $ E.HOM. : num  3066656 3066656 3066656 3066656 3066656 ...
# $ N_SITES: int  4803704 4803704 4803704 4803704 4803704 4803704 4803704 4803704 4803704 4803704 ...
# $ F      : num  0.1336 -0.0931 0.1771 -0.0119 0.1737 ...

# same order as ancestry .q file and metadata file
#observed homozygosity in SNPs
#expected homozygosity in SNPs assumed HWE (2pq)
# F - inbreeding coefficient; inverse of heterozygosity
    # F = 1 - (obs heterozygosity / exp heterozygosity)

# Combine the meta data, heterozygosity, and admixture data

het2 <- cbind(meta,het,Qscores_new) #bind the het results with the meta data

# How does F vary within each transect?

ggplot(het2, aes(x=Transect, y=F, color=Transect)) +
  geom_dotplot(binaxis='y', binwidth=0.01, stackdir='center')

# Plot admixture diversity spatially
ggplot(het2, aes(x=Longitude, y=Latitude, color=ShDiv)) +
  geom_point(size=4, shape=20)

# And finally, are more admixed individuals more heterozygous in their genomes?

ggplot(het2, aes(x=F, y=ShDiv, color=Transect)) +
  geom_point(size=4, shape=20)

cor.test(het2$F,het2$ShDiv)
