library(ggplot2)
library(gridExtra)

setwd("~/Documents/GitHub/EcologicalGenomics/Assignment3")

# Get the list of admixed individuals:
Admixed <- read.table("Admixed.Inds",header=F)

# Get the meta data:
meta <- read.table("Combined_Transect_Sampling_Data_2020.txt", sep="\t",header=T)

# merge them together:
meta_admx <- merge(meta, Admixed, by.x="ID", by.y="V1")
str(meta_admx)  

# Read in the Admixture coefficients for KBals that we made from the K=5 file:
KBals <- read.table("Admixed_KBals", sep="\t", header=F)
names(KBals) = c("ID","KBals")

# Second merge:
meta_admx_KBals <- merge(meta_admx,KBals,by="ID")


# Bring in phenotype data:
pheno <- read.table("climDat.txt",sep="\t",header=T)


# Merge pheno data with meta and KBals:
meta_admx_KBals_pheno <- merge(meta_admx_KBals,pheno,by="ID")

# Read in list of positions
snps <- read.table("Chr04.kept.sites",sep="\t", header=T)

#  cGDDfreeze
plotGDD <- ggplot(meta_admx_KBals_pheno,aes(x=KBals,y=mean_cGDDfreeze, color=Transect)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Growing Degree Days") 

plotGDD

# Final freeze
plotff <- ggplot(meta_admx_KBals_pheno,aes(x=KBals,y=mean_finalFreeze, color=Transect)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Final Freeze") 

plotff


# DD0
plotDD0 <- ggplot(meta_admx_KBals_pheno,aes(x=KBals,y=med_DD0, color=Transect)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Chilling Degree Days") 

plotDD0
#DD0 - strongest relationship
#Tricho ancestry take longer to break bud in spring than bals (in same environment)
#populations that have evolved in warmer climates tend to require more heat in spring to bud

grid.arrange(plotGDD, plotff, plotDD0, nrow = 3)

# linear models testing trait ~ genome-wide admixture association
summary(lm(mean_cGDDfreeze~KBals + Transect, data=meta_admx_KBals_pheno))
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -186.28  -60.24    0.52   57.18  417.07 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           389.22      15.66  24.860  < 2e-16 ***
#   KBals                  61.67      20.43   3.019  0.00271 ** 
#   TransectCassiar       122.73      19.67   6.238 1.16e-09 ***
#   TransectChilcotin     221.77      16.38  13.541  < 2e-16 ***
#   TransectCrowsnest     180.61      18.10   9.978  < 2e-16 ***
#   TransectIndian_Head    66.42      22.36   2.971  0.00316 ** 
#   TransectJasper        222.87      16.77  13.290  < 2e-16 ***
#   TransectWyoming       366.00      18.05  20.276  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 91.3 on 388 degrees of freedom
# Multiple R-squared:  0.5671,	Adjusted R-squared:  0.5593 
# F-statistic:  72.6 on 7 and 388 DF,  p-value: < 2.2e-16

summary(lm(mean_finalFreeze~KBals + Transect, data=meta_admx_KBals_pheno))

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -23.272  -5.652  -0.322   4.862  33.427 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          130.313      1.519  85.767  < 2e-16 ***
#   KBals                 14.532      1.982   7.330 1.35e-12 ***
#   TransectCassiar       14.228      1.909   7.452 6.03e-13 ***
#   TransectChilcotin      6.254      1.589   3.935 9.87e-05 ***
#   TransectCrowsnest      2.796      1.757   1.592   0.1123    
# TransectIndian_Head   -4.466      2.170  -2.058   0.0402 *  
#   TransectJasper         4.136      1.627   2.541   0.0114 *  
#   TransectWyoming       10.857      1.752   6.198 1.46e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 8.861 on 388 degrees of freedom
# Multiple R-squared:  0.2849,	Adjusted R-squared:  0.272 
# F-statistic: 22.08 on 7 and 388 DF,  p-value: < 2.2e-16

summary(lm(med_DD0~KBals + Transect, data=meta_admx_KBals_pheno))

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -49.979 -11.523   0.884  10.361  43.325 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          123.319      2.917  42.280  < 2e-16 ***
#   KBals                 53.953      3.806  14.177  < 2e-16 ***
#   TransectCassiar        9.423      3.665   2.571   0.0105 *  
#   TransectChilcotin    -20.535      3.051  -6.730 6.09e-11 ***
#   TransectCrowsnest    -31.361      3.372  -9.301  < 2e-16 ***
#   TransectIndian_Head  -27.314      4.165  -6.557 1.75e-10 ***
#   TransectJasper       -28.623      3.124  -9.162  < 2e-16 ***
#   TransectWyoming      -28.040      3.363  -8.339 1.32e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 17.01 on 388 degrees of freedom
# Multiple R-squared:  0.6679,	Adjusted R-squared:  0.6619 
# F-statistic: 111.5 on 7 and 388 DF,  p-value: < 2.2e-16


######  Bring in Association results from Plink   ######

mean_cGDDfreeze <- read.table("plink2.mean_cGDDfreeze.glm.linear",skip=1,sep="\t",header=F)
names(mean_cGDDfreeze) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")


#Row 1 - Beta 1 (additive test)
#Row 2 - Beta 2 (genome-wide test)

mean_cGDDfreeze2 <- mean_cGDDfreeze[which(mean_cGDDfreeze$TEST=="ADD"),]
#parsing out just beta 1/additive

# Define association outliers as the upper 0.1% of p-values

#########  mean_cGDDfreeze  #########
mean_cGDDfreeze2 <- cbind(snps, mean_cGDDfreeze2[,-c(1:2)])
mean_cGDDfreeze2$outlier = ifelse(mean_cGDDfreeze2$P<quantile(mean_cGDDfreeze2$P,0.001),2,1)

p1 <- ggplot(mean_cGDDfreeze2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=mean_cGDDfreeze2$outlier, color=mean_cGDDfreeze2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("mean_cGDDfreeze")

p1

####### mean_finalFreeze  #########
mean_finalFreeze <- read.table("plink2.mean_finalFreeze.glm.linear",skip=1,sep="\t",header=F)
names(mean_finalFreeze) = c("CHROM",  "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
mean_finalFreeze2 <- mean_finalFreeze[which(mean_finalFreeze$TEST=="ADD"),]
mean_finalFreeze2 <- cbind(snps, mean_finalFreeze2[,-c(1,2)])
mean_finalFreeze2$outlier = ifelse(mean_finalFreeze2$P<quantile(mean_finalFreeze2$P,0.001),2,1)

p2 <- ggplot(mean_finalFreeze2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=mean_finalFreeze2$outlier, color=mean_finalFreeze2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("mean_finalFreeze")

p2

#########  Bud med_DD0  #########
med_DD0 <- read.table("plink2.med_DD0.glm.linear",skip=1,sep="\t",header=F)
names(med_DD0) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
med_DD0 <- med_DD0[which(med_DD0$TEST=="ADD"),]
med_DD02 <- cbind(snps, med_DD0[,-c(1,2)])
med_DD02$outlier = ifelse(med_DD02$P<quantile(med_DD02$P,0.001),2,1)

p3 <- ggplot(med_DD02,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=med_DD02$outlier, color=med_DD02$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Bud med_DD0")

p3

grid.arrange(p1, p2, p3, nrow = 3)

# Get outliers for a given trait association:

mean_cGDDfreeze_outliers <- mean_cGDDfreeze2[which(mean_cGDDfreeze2$outlier==2),c(2,3,9)]
mean_finalFreeze_outliers <- mean_finalFreeze2[which(mean_finalFreeze2$outlier==2),c(2,3,9)]
med_DD0_outliers <- med_DD02[which(med_DD02$outlier==2),c(2,3,9)]


# Plot freq of LAI along chr
AF <- read.table("Chr04_LAI_freq.afreq", skip=1,sep="\t",header=F)
names(AF) = c("CHROM",  "ID",   "REF",  "ALT",  "ALT_FREQS",    "OBS_CT")
str(AF)

AF2 <- cbind(snps,AF)

windows <- seq(1,max(AF2$POS),5e4)
AF_windows <- numeric()

for(i in 1:length(windows)){
  tmp=AF2[which(AF2$POS>windows[i] & AF2$POS<windows[i+1]),"ALT_FREQS"]
  ancfreq=mean(tmp)
  AF_windows[i] = ancfreq
}

AF3 <- as.data.frame(cbind(windows,AF_windows))
names(AF3) = c("window","AvgAncFreq")

upper = mean(AF3$AvgAncFreq,na.rm=T) + 2*sd(AF3$AvgAncFreq,na.rm=T)
lower = mean(AF3$AvgAncFreq,na.rm=T) - 2*sd(AF3$AvgAncFreq,na.rm=T)

outliers_upper = AF3[which(AF3$AvgAncFreq>upper),]
outliers_lower = AF3[which(AF3$AvgAncFreq<lower),]

# Print the outlier regions out
outliers_upper
outliers_lower

# And finally, make the 4-panel plot with the trait associations
p4 <- ggplot(AF3[,-3],aes(x=window,y=AvgAncFreq)) +
  geom_line(size=0.8, color="blue") + 
  xlab("Position (bp) along chromosome") +
  ylab("Frequency P. trichocarpa ancestry") +
  geom_hline(yintercept=mean(AF2$ALT_FREQS), color = "red") + 
  geom_hline(yintercept=upper, linetype="dashed", color = "red") + 
  geom_hline(yintercept=lower, linetype="dashed", color = "red") +
  ggtitle("Chr04: Local ancestry")

p4


grid.arrange(p1, p2, p3, p4, nrow = 4)

# Get the betas from each trait and look at pleiotropy between traits
betas <- cbind(mean_cGDDfreeze[,c(1:3,9)],mean_finalFreeze2[,9],med_DD02[,9])
names(betas) = c("CHROM","POS","ID","beta_mean_cGDDfreeze","beta_mean_finalFreeze","beta_med_DD0")
str(betas)

cor(betas[,4:6],betas[4:6])
# beta_mean_cGDDfreeze    beta_mean_finalFreeze   beta_med_DD0
# beta_mean_cGDDfreeze   1.000000000 -0.05717757 -0.008154436
# beta_mean_finalFreeze   -0.057177569  1.00000000  0.782707741
# beta_med_DD0 -0.008154436  0.78270774  1.000000000

plot(beta$beta_mean_finalFreeze,betas$beta_med_DD0)

p5 <- ggplot(betas,aes(x=beta_med_DD0,y=beta_mean_finalFreeze)) +
  geom_point(color="darkgray") + 
  xlab("Beta med_DD0") +
  ylab("Beta mean_finalFreeze") +
  ggtitle("Correlation of mean_finalFreeze and med_DD0 effect sizes")

p5

##Enchrichment analysis

# Import the GFF annotation file and make a transcript database
txdb <- makeTxDbFromGFF("Ptrichocarpa_533_v4.1.gene.gff3.gz", format="gff3")
CHR="Chr04"

head(seqlevels(txdb))

# Subset the database for just your chromosome of interest
seqlevels(txdb) <- CHR # subset for just your chromosome

# Reduce the transcript database to just the non-redundant gene names, instead of multiple entries for all the variant transcript types per gene
genes <- unlist(reduce(transcriptsBy(txdb, by="gene"))) 
#GRanges object with 2078 ranges and 1 metadata column
genes$geneID <- names(genes)

mean_cGDDfreeze_outliers2 <- mean_cGDDfreeze2[which(mean_cGDDfreeze2$outlier==2),c(1,2,3,9)]
mean_finalFreeze_outliers2 <- mean_finalFreeze2[which(mean_finalFreeze2$outlier==2),c(1,2,3,9)]
med_DD0_outliers2 <- med_DD02[which(med_DD02$outlier==2),c(1,2,3,9)]

outlier_range_cGDD <- makeGRangesFromDataFrame(mean_cGDDfreeze_outliers2, seqnames = 'CHROM', start.field = 'POS', end.field = 'POS')
candGenesCGGD <- subsetByOverlaps(genes, outlier_range_cGDD)
write.table(candGenesCGGD$geneID, paste0("candGenesCGDD",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")

outlier_range_finalFreeze <- makeGRangesFromDataFrame(mean_finalFreeze_outliers2, seqnames = 'CHROM', start.field = 'POS', end.field = 'POS')
candGenesFF <- subsetByOverlaps(genes, outlier_range_finalFreeze)
write.table(candGenesFF$geneID, paste0("candGenesFF",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")

outlier_range_dd0 <- makeGRangesFromDataFrame(med_DD0_outliers2, seqnames = 'CHROM', start.field = 'POS', end.field = 'POS')
candGenesDD <- subsetByOverlaps(genes, outlier_range_dd0)
write.table(candGenesDD$geneID, paste0("candGenesDD",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")
