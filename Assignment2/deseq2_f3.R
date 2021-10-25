# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("DESeq2")

# practicing deseq2 with transplant (RT) tonsa data
# 10/6/21

setwd("~/Documents/GitHub/EcologicalGenomics/Assignment2/tonsaF3")

library(DESeq2) 
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)
library(hexbin)
library(eulerr)
library(pheatmap)

# import counts matrix

countsTable <- read.table("DE_counts_F3.txt", header=T, row.names = 1)
head(countsTable)
dim(countsTable) # 25279  x  16

countsTableRound <- round(countsTable) #removing decimals b/c dseq2 doesn't like them
head(countsTableRound)

# import the sample description table
conds <- read.delim("RT_tonsa_F3_samples.txt", header=T, stringsAsFactors = T, row.names = 1)
head(conds)

### Visualizing Data

#check how many reads we have from each sample
colSums(countsTableRound)
mean(colSums(countsTableRound))
# 17.7 million reads per sample
barplot(colSums(countsTableRound), names.arg=colnames(countsTableRound), cex.names = 0.5, las=3,ylim=c(0,20000000)) #cex.names changes size, las changes orientation
abline(h=mean(colSums(countsTableRound)), col="blue",lwd=2)

#The average number counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) #11220.54
median(rowSums(countsTableRound)) #2144
#median much lower than mean; expected with distribution with long tail

apply(countsTableRound, 1, mean) #across rows
apply(countsTableRound, 2, mean) #across columns
hist(apply(countsTableRound, 1, mean), xlim=c(0,1000), breaks=10000)
hist(apply(countsTableRound, 1, mean), xlim=c(50000,150000),ylim=c(0,10), breaks=100)

##########################################################################

### create a DESeq object and define experimental design here with the tilda

dds <- DESeqDataSetFromMatrix(countData=countsTableRound, 
                              colData=conds, 
                              design=~line+environment+line:environment) #[cond]:[cond] = interaction

dim(dds) #25279  16

### Filter out genes with too few reads - keep reads with average > 10 reads per sample

dds <- dds[rowSums(counts(dds)) > 160] #10 * 16 samples
dim(dds)
     
### Run DESeq Model to test for differential gene expression
dds <- DESeq(dds)

### List the results you've generated
resultsNames(dds)
# [1] "Intercept"                 
# [2] "line_combined_vs_ambient"  
# [3] "environment_HH_vs_AA"      
# [4] "linecombined.environmentHH"

### Generate PCA to look at global gene expression patterns

vsd <- vst(dds, blind=F)

data <- plotPCA(vsd, intgroup=c("line","environment"),returnData=T)
percentVar <- round(100 *attr(data,"percentVar"))

ggplot(data, aes(PC1,PC2,color=environment,shape=line)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()

### Order and summarize results from specific contrasts

#### INTERACTION EFFECT ####

resInteraction <- results(dds, alpha=0.05)
resInteraction <- resInteraction[order(resInteraction$padj),]
head(resInteraction)
summary(resInteraction)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 271, 1.1%
# LFC < 0 (down)     : 60, 0.24%
# outliers [1]       : 5, 0.02%
# low counts [2]     : 2451, 9.7%
# (mean count < 26)

# significant interaction: cross of lines b/w conditions
# significant line: parallel but diff lines b/w conditions
# significant environment:parallel/equal lines (points clustered)


### ENVIRONMENT EFFECT ###
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, 
                              colData = conds, 
                              design = ~ line + environment)
dds <- DESeq(dds, test="LRT", reduced=~line)
resultsNames(dds)
# [1] "Intercept"               
# [2] "line_combined_vs_ambient"
# [3] "environment_HH_vs_AA"   

resEnv <- results(dds, alpha = 0.05)
resEnv <- resEnv[order(resEnv$padj),] #ordering rows by adjusted p-value
head(resEnv)

summary(resEnv)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 513, 2%
# LFC < 0 (down)     : 315, 1.2%
# outliers [1]       : 16, 0.063%
# low counts [2]     : 491, 1.9%
# (mean count < 19)

resEnv <- resEnv[!is.na(resEnv$padj),] #filters data table to exclude NA results; lost about 500 genes

degsEnv <- row.names(resEnv[resEnv$padj < 0.05,]) #differentially expressed genes from the environment
length(degsEnv) 
#828 genes differentially expressed


### LINE EFFECT ###
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, 
                              colData = conds, 
                              design = ~ environment + line)

dds <- DESeq(dds, test="LRT", reduced=~environment)
resultsNames(dds)
# [1] "Intercept"               
# [2] "environment_HH_vs_AA"    
# [3] "line_combined_vs_ambient"

resLine <- results(dds, alpha = 0.05)
resLine <- resLine[order(resLine$padj),]
head(resLine)


summary(resLine)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 808, 3.2%
# LFC < 0 (down)     : 837, 3.3%
# outliers [1]       : 16, 0.063%
# low counts [2]     : 981, 3.9%
# (mean count < 21)

resLine <- resLine[!is.na(resLine$padj),]

degsline <- row.names(resLine[resLine$padj < 0.05,])
length(degsline)
#1645 differentially expressed gene


### INTERACTION EFFECT PT2 - using LRT ###
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, 
                              colData = conds, 
                              design = ~ environment + line + environment:line)

dds <- DESeq(dds, test="LRT", reduced=~environment + line)
resultsNames(dds)
# [1] "Intercept"                 
# [2] "environment_HH_vs_AA"      
# [3] "line_combined_vs_ambient"  
# [4] "environmentHH.linecombined"
resInt <- results(dds, alpha = 0.05)
resInt <- resInt[order(resInt$padj),]
head(resInt)

summary(resInt)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 235, 0.93%
# LFC < 0 (down)     : 48, 0.19%
# outliers [1]       : 5, 0.02%
# low counts [2]     : 2941, 12%
# (mean count < 28)

resInt <- resInt[!is.na(resInt$padj),]

degsInt <- row.names(resInt[resInt$padj < 0.05,])
length(degsInt)
#283 differentially expressed genes

#### Data Visualization

# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds, gene="TRINITY_DN115950_c0_g1", intgroup = (c("line","environment")), returnData=TRUE)
d

p <-ggplot(d, aes(x=environment, y=count, color=line, shape=line, group=line)) + 
  theme_minimal() + theme(text = element_text(size=15), panel.grid.major=element_line(colour="grey"))+
  geom_point(position=position_jitter(w=0.2,h=0), size=3)+ 
  stat_summary(fun = mean, geom = "line") + 
  stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
print(p)

#counts for top environment gene
d1 <-plotCounts(dds, gene="TRINITY_DN138549_c1_g2", intgroup = (c("line","environment")), returnData=TRUE)
d1

p1 <-ggplot(d1, aes(x=environment, y=count, color=line, shape=line, group=line)) + 
  theme_minimal() + theme(text = element_text(size=15), panel.grid.major=element_line(colour="grey"))+
  geom_point(position=position_jitter(w=0.2,h=0), size=3)+ 
  stat_summary(fun = mean, geom = "line") + 
  stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
print(p1)

## VENN DIAGRAM ##
# Total
length(degsEnv)  # 828
length(degsline)  # 1645
length(degsInt)  # 283

# Intersections
length(intersect(degsEnv,degsline))  # 141
length(intersect(degsEnv,degsInt))  # 14
length(intersect(degsInt,degsline))  # 32

intEL <- intersect(degsEnv,degsline)
length(intersect(degsInt,intEL)) # 7

# Number unique
828-14-141-7 # 666
1645-141-32-7 # 1465
283-14-32-7 # 230

# unique, present in 2, present in 3
fit1 <- euler(c("Env" = 666, "Line" = 1465, "Interaction" = 230, "Env&Line" = 141, "Env&Interaction" = 14, "Line&Interaction" = 32, "Env&Line&Interaction" = 7))

plot(fit1,  lty = 1:3, quantities = TRUE)

plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))

## HEATMAP ##
# Heatmap of top 20 genes sorted by pvalue
# Interaction

topgenes <- head(rownames(resInt),50)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df, fontsize_row = 6)
#dendogram - cluster by similarity across the samples
#heatmap colors - up and downregulation

# By line

topgenes <- head(rownames(resLine),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat) # scaling - subtract row mean from matrix
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)
#clustered by line, and then by environment
#see genes that are differentially regulated in different lines

# By environment

topgenes <- head(rownames(resEnv),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat) # scaling - subtract row mean from matrix
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)
#clustered by line, and then by environment
#see genes that are differentially regulated in different lines

