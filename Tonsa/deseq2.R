# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("DESeq2")

# practicing deseq2 with transplant (RT) tonsa data
# 10/6/21

setwd("~/Documents/GitHub/EcologicalGenomics/Tonsa")

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

countsTable <- read.table("DE_counts_F1.txt", header=T, row.names = 1)
head(countsTable)
dim(countsTable) # 24362  x  16

countsTableRound <- round(countsTable) #removing decimals b/c dseq2 doesn't like them
head(countsTableRound)

# import the sample description table
conds <- read.delim("RT_tonsa_F1_samples.txt", header=T, stringsAsFactors = T, row.names = 1)
head(conds)

### Visualizing Data

#check how many reads we have from each sample
colSums(countsTableRound)
mean(colSums(countsTableRound))
# just over 18 million reads per sample
barplot(colSums(countsTableRound), names.arg=colnames(countsTableRound), cex.names = 0.5, las=3,ylim=c(0,20000000)) #cex.names changes size, las changes orientation
abline(h=mean(colSums(countsTableRound)), col="blue",lwd=2)

#The average number counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) #11930.81
median(rowSums(countsTableRound)) #2226
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

dim(dds)

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
# out of 24362 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2839, 12% (NOTE: Increased Log Fold Change)
# LFC < 0 (down)     : 1053, 4.3% (NOTE: Decreased Log Fold Change)
# outliers [1]       : 9, 0.037%
# low counts [2]     : 473, 1.9%
# (mean count < 18)

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
# out of 24362 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 213, 0.87%
# LFC < 0 (down)     : 235, 0.96%
# outliers [1]       : 41, 0.17%
# low counts [2]     : 473, 1.9%
# (mean count < 18)

resEnv <- resEnv[!is.na(resEnv$padj),] #filters data table to exclude NA results; lost about 500 genes

degsEnv <- row.names(resEnv[resEnv$padj < 0.05,]) #differentially expressed genes from the environment
length(degsEnv) 
#448 genes differentially expressed


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
# out of 24362 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 72, 0.3%
# LFC < 0 (down)     : 154, 0.63%
# outliers [1]       : 41, 0.17%
# low counts [2]     : 2362, 9.7%
# (mean count < 25)

resLine <- resLine[!is.na(resLine$padj),]

degsline <- row.names(resLine[resLine$padj < 0.05,])
length(degsline)
#226 differentially expressed gene


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
# out of 24362 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2802, 12%
# LFC < 0 (down)     : 1052, 4.3%
# outliers [1]       : 9, 0.037%
# low counts [2]     : 945, 3.9%
# (mean count < 20)

resInt <- resInt[!is.na(resInt$padj),]

degsInt <- row.names(resInt[resInt$padj < 0.05,])
length(degsInt)
#3854 differentially expressed genes

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
length(degsEnv)  # 448
length(degsline)  # 226
length(degsInt)  # 3854

# Intersections
length(intersect(degsEnv,degsline))  # 37
length(intersect(degsEnv,degsInt))  # 44
length(intersect(degsInt,degsline))  # 34

intEL <- intersect(degsEnv,degsline)
length(intersect(degsInt,intEL)) # 7

# Number unique
448-44-37-7 # 360 
226-37-34-7 # 148
3854-44-34-7 # 3769

# unique, present in 2, present in 3
fit1 <- euler(c("Env" = 360, "Line" = 148, "Interaction" = 3769, "Env&Line" = 37, "Env&Interaction" = 44, "Line&Interaction" = 34, "Env&Line&Interaction" = 7))

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
pheatmap(mat, annotation_col=df)
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