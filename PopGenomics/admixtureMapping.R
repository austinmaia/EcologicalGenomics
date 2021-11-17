library(ggplot2)
library(gridExtra)

setwd("~/Documents/GitHub/EcologicalGenomics/PopGenomics/Chr04")

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
pheno <- read.table("VT_Garden_Phenotypes_2021.txt",sep="\t",header=T)

# Merge pheno data with meta and KBals:
meta_admx_KBals_pheno <- merge(meta_admx_KBals,pheno,by="ID")
#  Rust
plotRust <- ggplot(meta_admx_KBals_pheno,aes(x=KBals,y=RUST, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Rust susceptibility") 

# Read in list of positions
snps <- read.table("Chr04.kept.sites",sep="\t", header=T)


plotRust

# Bud set
plotBudset <- ggplot(meta_admx_KBals_pheno,aes(x=KBals,y=SET, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Bud set") 

plotBudset
#balsamifera tend to winterize later (?) in same environment

# Bud flush
plotBudflush <- ggplot(meta_admx_KBals_pheno,aes(x=KBals,y=FLUSH, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Bud flush") 

plotBudflush
#Bud flush - strongest relationship
#Tricho ancestry take longer to break bud in spring than bals (in same environment)
#populations that have evolved in warmer climates tend to require more heat in spring to bud

grid.arrange(plotRust, plotBudset, plotBudflush, nrow = 3)

# linear models testing trait ~ genome-wide admixture association
summary(lm(RUST~KBals + Transect.x, data=meta_admx_KBals_pheno))
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.62811 -0.09974  0.04646  0.12916  0.38583 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           -0.09572    0.03414  -2.804  0.00530 ** 
#   KBals                 -0.08407    0.04454  -1.887  0.05984 .  
# Transect.xCassiar      0.11772    0.04290   2.744  0.00634 ** 
#   Transect.xChilcotin    0.10476    0.03571   2.934  0.00355 ** 
#   Transect.xCrowsnest    0.15826    0.03947   4.010 7.28e-05 ***
#   Transect.xIndian_Head  0.13497    0.04875   2.769  0.00590 ** 
#   Transect.xJasper       0.10661    0.03656   2.916  0.00376 ** 
#   Transect.xWyoming      0.26945    0.03936   6.846 2.97e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1991 on 388 degrees of freedom
# Multiple R-squared:  0.1725,	Adjusted R-squared:  0.1576 
# F-statistic: 11.56 on 7 and 388 DF,  p-value: 2.357e-13


summary(lm(SET~KBals + Transect.x, data=meta_admx_KBals_pheno))

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -33.685  -7.920  -0.665   8.342  36.650 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             -9.746      2.025  -4.812 2.14e-06 ***
#   KBals                  -26.539      2.642 -10.044  < 2e-16 ***
#   Transect.xCassiar        3.188      2.545   1.253    0.211    
# Transect.xChilcotin     14.403      2.118   6.799 4.00e-11 ***
#   Transect.xCrowsnest     30.544      2.341  13.046  < 2e-16 ***
#   Transect.xIndian_Head   24.354      2.892   8.421 7.33e-16 ***
#   Transect.xJasper        22.768      2.169  10.496  < 2e-16 ***
#   Transect.xWyoming       30.276      2.335  12.968  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 11.81 on 388 degrees of freedom
# Multiple R-squared:  0.6482,	Adjusted R-squared:  0.6419 
# F-statistic: 102.1 on 7 and 388 DF,  p-value: < 2.2e-16

summary(lm(FLUSH~KBals + Transect.x, data=meta_admx_KBals_pheno))

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -9.8893 -3.1791 -0.5628  2.1129 17.5910 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           -1.29280    0.77350  -1.671  0.09546 .  
# KBals                 -4.79023    1.00922  -4.746 2.92e-06 ***
#   Transect.xCassiar      0.07921    0.97196   0.081  0.93509    
# Transect.xChilcotin    1.70354    0.80912   2.105  0.03590 *  
#   Transect.xCrowsnest    5.46074    0.89423   6.107 2.47e-09 ***
#   Transect.xIndian_Head  3.15875    1.10464   2.860  0.00447 ** 
#   Transect.xJasper       1.62067    0.82851   1.956  0.05117 .  
# Transect.xWyoming      6.81336    0.89176   7.640 1.71e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 4.511 on 388 degrees of freedom
# Multiple R-squared:  0.3675,	Adjusted R-squared:  0.3561 
# F-statistic: 32.21 on 7 and 388 DF,  p-value: < 2.2e-16


######  Bring in Association results from Plink   ######

rust <- read.table("plink2.RUST.glm.linear",skip=1,sep="\t",header=F)
names(rust) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")

# CHROM POS    ID  REF ALT A1 TEST OBS_CT       BETA        SE   T_STAT
# 1     0   0 snp01 TRUE   A  A  ADD    396  0.0353984 0.0170634  2.07452
# 2     0   0 snp01 TRUE   A  A  K2Q    396 -0.1080120 0.0469478 -2.30068
# 3     0   0 snp02 TRUE   A  A  ADD    396  0.0353984 0.0170634  2.07452
# 4     0   0 snp02 TRUE   A  A  K2Q    396 -0.1080120 0.0469478 -2.30068
#Row 1 - Beta 1 (additive test)
#Row 2 - Beta 2 (genome-wide test)

rust2 <- rust[which(rust$TEST=="ADD"),]
#parsing out just beta 1/additive

# Define association outliers as the upper 0.1% of p-values

#########  rust  #########
rust2 <- cbind(snps, rust2[,-c(1:2)])
rust2$outlier = ifelse(rust2$P<quantile(rust2$P,0.001),2,1)

p1 <- ggplot(rust2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=rust2$outlier, color=rust2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Rust infection")

p1

####### Bud set  #########
budset <- read.table("plink2.SET.glm.linear",skip=1,sep="\t",header=F)
names(budset) = c("CHROM",  "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
budset2 <- budset[which(budset$TEST=="ADD"),]
budset2 <- cbind(snps, budset2[,-c(1,2)])
budset2$outlier = ifelse(budset2$P<quantile(budset2$P,0.001),2,1)

p2 <- ggplot(budset2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=budset2$outlier, color=budset2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Bud set")

p2

#########  Bud flush  #########
budflush <- read.table("plink2.FLUSH.glm.linear",skip=1,sep="\t",header=F)
names(budflush) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
budflush <- budflush[which(budflush$TEST=="ADD"),]
budflush2 <- cbind(snps, budflush[,-c(1,2)])
budflush2$outlier = ifelse(budflush2$P<quantile(budflush2$P,0.001),2,1)

p3 <- ggplot(budflush2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=budflush2$outlier, color=budflush2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Bud flush")

p3

grid.arrange(p1, p2, p3, nrow = 3)

# Get outliers for a given trait association:

rust_outliers <- rust2[which(rust2$outlier==2),c(2,3,9)]
set_outliers <- budset2[which(budset2$outlier==2),c(2,3,9)]
flush_outliers <- budflush2[which(budflush2$outlier==2),c(2,3,9)]


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
betas <- cbind(rust[,c(1:3,9)],budset2[,9],budflush2[,9])
names(betas) = c("CHROM","POS","ID","beta_rust","beta_set","beta_flush")
str(betas)

cor(betas[,4:6],betas[4:6])
# beta_rust    beta_set   beta_flush
# beta_rust   1.000000000 -0.05717757 -0.008154436
# beta_set   -0.057177569  1.00000000  0.782707741
# beta_flush -0.008154436  0.78270774  1.000000000

plot(beta$beta_set,betas$beta_flush)

p5 <- ggplot(betas,aes(x=beta_flush,y=beta_set)) +
  geom_point(color="darkgray") + 
  xlab("Beta bud flush") +
  ylab("Beta bud set") +
  ggtitle("Correlation of bud set and flush effect sizes")

p5

