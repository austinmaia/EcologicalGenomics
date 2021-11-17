# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("GenomicRanges")
# BiocManager::install("GenomicFeatures")

library(GenomicRanges)
library(GenomicFeatures)

####### Randomization test for Fst and sweep outliers #########

setwd("~/Documents/GitHub/EcologicalGenomics/PopGenomics/Chr04") # Customize to where your Fst, Pi, and RAiSD files are located on your laptop

fst <- read.table("Bals_Tricho_All.windowed.weir.fst", sep="\t",header=T) # Import the Fst results

cent <- read.table("Chr04_centromere.txt", sep="\t",header=F) # Import the centromere coordinates

fst <- fst[-which(fst$BIN_START>cent$V2 & fst$BIN_END<cent$V3),] # Mask Fst windows in the centromere region


# Calculate Genomic Ranges for outlier bins

CHR="Chr04"  # Customize to your chromosome.  Make sure syntax is exact!

# Define the genomic ranges of the Fst bins
fstGR <- GRanges(CHR,IRanges(fst[,"BIN_START"],fst[,"BIN_END"]))

# GRanges object with 481 ranges and 1 metadata column:
#   seqnames            ranges strand |
#   <Rle>         <IRanges>  <Rle> |
#   [1]    Chr04           1-50000      * |
#   [2]    Chr04      50001-100000      * |
#   [3]    Chr04     100001-150000      * |
#   [4]    Chr04     150001-200000      * |
#   [5]    Chr04     200001-250000      * |
#   ...      ...               ...    ... .
# [477]    Chr04 23900001-23950000      * |


# Define what should be a "low" value of Fst.  Here, we'll try the lowest 10% of windows that don't incude the centromere.  Can play with this if you want (choosing the last value in the quantile function call). 
fstThreshold = quantile(fst$MEAN_FST,0.1)
# 10% 
# 0.0948187

# Call outliers (=1) or non-outliers (=2) based on fst Threshold
fstGR$outlier <- ifelse(fst$MEAN_FST<fstThreshold,1,2)

# Grab just the outlier Fst windows
fstCand <- subset(fstGR, outlier==1)
#GRanges object with 48 ranges and 1 metadata column:
# seqnames            ranges strand |
#   <Rle>         <IRanges>  <Rle> |
#   [1]    Chr04     500001-550000      * |
#   [2]    Chr04     600001-650000      * |
#   [3]    Chr04     650001-700000      * |
#   [4]    Chr04     700001-750000      * |
#   [5]    Chr04     900001-950000      * |
#   ...      ...               ...    ... .
# [44]    Chr04 20700001-20750000      * |
#   [45]    Chr04 20750001-20800000      * |
#   [46]    Chr04 21150001-21200000      * |
#   [47]    Chr04 23550001-23600000      * |
#   [48]    Chr04 23600001-23650000      * |


raisd <- read.table(paste0("RAiSD_Report.",CHR,".",CHR), sep="\t",header=F) # Import the RAiSD results

# Define RAiSD genomic ranges based on first and last SNP within the RAiSD windows of 50 SNPs each.   How does this relate to the Fst window size?
raisdGR <- GRanges(CHR,IRanges(raisd[,2],raisd[,3]))

# GRanges object with 278028 ranges and 1 metadata column:
#   seqnames            ranges strand
# <Rle>         <IRanges>  <Rle>
#   [1]    Chr04       10911-15403      *
#   [2]    Chr04       11054-15454      *
#   [3]    Chr04       11173-15463      *
#   [4]    Chr04       11177-15475      *
#   [5]    Chr04       11180-15514      *
#   ...      ...               ...    ...
# [278024]    Chr04 24123691-24133757      *

# Define RAiSD outliers, based on the upper 1% of SNPs
raisdThreshold = quantile(raisd$V7,0.99)
# 99% 
# 1.183e-10

raisdGR$outlier <- ifelse(raisd$V7>raisdThreshold,"outlier","non")

raisdGR_out <- unlist(reduce(split(raisdGR, ~outlier)))
raisdGR_out$outlier <- names(raisdGR_out)
raisdCand <- subset(raisdGR_out, outlier=="outlier")

# GRanges object with 65 ranges and 1 metadata column:
#   seqnames            ranges strand |
#   <Rle>         <IRanges>  <Rle> |
#   outlier    Chr04      95359-111075      * |
#   outlier    Chr04     633615-734233      * |
#   outlier    Chr04     757109-791640      * |
#   outlier    Chr04     906670-936465      * |
#   outlier    Chr04     947558-976241      * |
#   ...      ...               ...    ... .
# outlier    Chr04 23373528-23398662      * |

# Use GenomicRanges to pull out the RAiSD outlier windows that overlap the lowest Fst windows
overlap <- subsetByOverlaps(raisdCand, fstCand)

length(overlap) # Number of the RAiSD sweep candidate loci that overlap with the lowest n% of Fst windows
#[1] 18


# A little code I wrote to permute Fst values randomly among the 50kb windows, and re-calculate the overlap with the RAiSD outliers. 
# Can set number of permutation replicates (NumPerm)
# 1000 permutations seems to give a decent randomized distribution

NumPerm=1000

Perm = NA
for(i in 1:NumPerm){
  FstSamp = cbind(fst[,c(2,3)], fst[sample(nrow(fst)),6])
  names(FstSamp) = c("Start", "Stop", "FST")
  FstRand = FstSamp[which(FstSamp[3]<fstThreshold),]  
  FstRand_ranges <- GRanges(seqnames=CHR, ranges=IRanges(FstRand[,1],FstRand[,2]))
  FstRand_ranges_red <- reduce(FstRand_ranges)
  Perm[i] = sum(countOverlaps(raisdCand, FstRand_ranges_red))
}

# Plot the random distribution with a blue line marking the observed overlap
hist(Perm, col="gray", main="Chr04 Randomization test", xlab="Overlapping windows of Fst and RAiSD outliers")
abline(v=length(overlap), col="blue", lwd=3)

# Calculate the p-pvalue based on the rank of the observed value in the randomized distirbution.
# Note the 1-tailed test here (alpha = 0.05)

p_value = 1-ecdf(Perm)(length(overlap))
p_value
#[1] 0.002

# Import the GFF annotation file and make a transcript database
txdb <- makeTxDbFromGFF("Ptrichocarpa_533_v4.1.gene.gff3.gz", format="gff3")

txdb

#TxDb object:
  # Db type: TxDb
  # Supporting package: GenomicFeatures
  # Data source: Ptrichocarpa_533_v4.1.gene.gff3.gz
  # Organism: NA
  # Taxonomy ID: NA
  # miRBase build ID: NA
  # Genome: NA
  # Nb of transcripts: 52400
  # Db created by: GenomicFeatures package from Bioconductor
  # Creation time: 2021-11-01 10:04:15 -0400 (Mon, 01 Nov 2021)
  # GenomicFeatures version at creation time: 1.42.3
# RSQLite version at creation time: 2.2.8
# DBSCHEMAVERSION: 1.2

# How many chromosomes are present?
head(seqlevels(txdb))

# Subset the database for just your chromosome of interest
seqlevels(txdb) <- CHR # subset for just your chromosome

# Reduce the transcript database to just the non-redundant gene names, instead of multiple entries for all the variant transcript types per gene
genes <- unlist(reduce(transcriptsBy(txdb, by="gene"))) 
#GRanges object with 2078 ranges and 1 metadata column
genes$geneID <- names(genes)
# [1] "Potri.004G000150" "Potri.004G000300"
# [3] "Potri.004G000400" "Potri.004G000600"
# [5] "Potri.004G000750" "Potri.004G000900"
# [7] "Potri.004G001000" "Potri.004G001200"
# [9] "Potri.004G001300" "Potri.004G001400"
# [11] "Potri.004G001500" "Potri.004G001600"
# [13] "Potri.004G001701" "Potri.004G001800"
# [15] "Potri.004G001900" "Potri.004G002000"

candGenes <- subsetByOverlaps(genes, overlap)
# GRanges object with 84 ranges and 1 metadata column:
#   seqnames            ranges strand |
#   <Rle>         <IRanges>  <Rle> |
#   Potri.004G011101    Chr04     653019-654971      + |
#   Potri.004G011201    Chr04     658506-658842      - |
#   Potri.004G011300    Chr04     657197-659949      + |
#   Potri.004G011375    Chr04     683629-684275      + |
#   Potri.004G011450    Chr04     697784-699295      + |
write.table(candGenes$geneID, paste0("candGenes",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")
# GO- Biological Process
# GO:0019538	3.051e-05	14/17 | 3727/14903	protein metabolic process
# GO:0016310	1.118e-04	10/17 | 1964/14903	phosphorylation
# GO:0006468	1.308e-04	10/17 | 1911/14903	protein phosphorylation
# GO:0043170	1.494e-04	14/17 | 4689/14903	macromolecule metabolic process
# GO:0006793	1.657e-04	10/17 | 2319/14903	phosphorus metabolic process
# GO:0006464	1.733e-04	10/17 | 2227/14903	cellular protein modification process
# GO:0036211	1.733e-04	10/17 | 2227/14903	protein modification process
# GO:0006796	1.836e-04	10/17 | 2315/14903	phosphate-containing compound metabolic process
# GO:0043412	1.890e-04	10/17 | 2288/14903	macromolecule modification
# GO:0044267	1.435e-03	10/17 | 3010/14903	cellular protein metabolic process
# GO:0044238	5.272e-03	14/17 | 6540/14903	primary metabolic process
# GO:0071704	6.252e-03	14/17 | 6836/14903	organic substance metabolic process

# GO - Molecular function
# GO:0045735	1.707e-06	5/22 | 92/19622	nutrient reservoir activity
# GO:0008234	5.275e-05	4/22 | 92/19622	cysteine-type peptidase activity
# GO:0004672	2.149e-04	10/22 | 1946/19622	protein kinase activity
# GO:0016301	3.094e-04	10/22 | 2150/19622	kinase activity
# GO:0016773	3.242e-04	10/22 | 2107/19622	phosphotransferase activity, alcohol group as acceptor
# GO:0030554	6.268e-04	12/22 | 3643/19622	adenyl nucleotide binding
# GO:0016772	6.541e-04	10/22 | 2395/19622	transferase activity, transferring phosphorus-containing groups
# GO:0032559	6.956e-04	12/22 | 3638/19622	adenyl ribonucleotide binding
# GO:0005524	6.956e-04	12/22 | 3638/19622	ATP binding
# GO:0032553	8.122e-04	12/22 | 3989/19622	ribonucleotide binding
# GO:0017076	8.217e-04	12/22 | 3969/19622	purine nucleotide binding
# GO:0001882	8.700e-04	12/22 | 3966/19622	nucleoside binding
# GO:0032549	9.299e-04	12/22 | 3965/19622	ribonucleoside binding
# GO:0001883	9.794e-04	12/22 | 3956/19622	purine nucleoside binding
# GO:0032550	9.794e-04	12/22 | 3956/19622	purine ribonucleoside binding
# GO:0032555	9.794e-04	12/22 | 3956/19622	purine ribonucleotide binding
# GO:0035639	9.794e-04	12/22 | 3956/19622	purine ribonucleoside triphosphate binding

# Pfam
# PF00190	5.072e-08	5/22 | 92/29014	Cupin
# PF00112	1.952e-07	4/22 | 50/29014	Papain family cysteine protease
# PF01657	7.546e-07	4/22 | 77/29014	Domain of unknown function DUF26
# PF00954	7.823e-03	2/22 | 132/29014	S-locus glycoprotein family