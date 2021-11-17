
library(ggplot2)

setwd("~/Documents/GitHub/EcologicalGenomics")

# Read in list of positions
snps <- read.table("Chr04.kept.sites",sep="\t", header=T)

# Read in the local ancestry frequencies from Plink
AF <- read.table("Chr04_LAI_freq.afreq", skip=1,sep="\t",header=F)  

# Note the skip=1 here.  
# This skips the first line, since Plink's header line doesn't play well with R.  
# We'll define our own header line below.

names(AF) = c("CHROM",  "ID",   "REF",  "ALT",  "ALT_FREQS",    "OBS_CT")

AF2 <- cbind(snps,AF)

str(AF2) # How does it look?

# A simple plot:

p1 <- ggplot(AF2[,-3],aes(x=POS,y=ALT_FREQS)) +
  geom_line(size=0.25, color="blue") + 
  xlab("Position (bp) along chromosome") +
  ylab("Frequency P. trichocarpa ancestry")

p1

