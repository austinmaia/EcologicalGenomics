# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("DESeq2")

# practicing deseq2 with transplant (RT) tonsa data
# 10/6/21

library(DESeq2) 
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)
library(hexbin)

# import counts matrix

countsTable <- read.table("DE_counts_F1.txt", header=T, row.names = 1)
head(countsTable)
dim(countsTable) # 24362  x  16

countsTableRound <- round(countsTable) #removing decimals b/c dseq2 doesn't like them
head(countsTableRound)

# import the sample description table
conds <- read.delim("RT_tonsa_F1_samples.txt", header=T, stringsAsFactors = T, row.names = 1)
head(conds)
