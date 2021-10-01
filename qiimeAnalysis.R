library(qiime2R)
library(tidyverse)
library(gplots)
library("ggsci")
library(lemon)

shannon <- read_qza("diversityMetrics/shannon_vector.qza")
metadata <- read.table("diversityMetrics/pyc_manifest", header=T)
uwunifrac<-read_qza("diversityMetrics/unweighted_unifrac_pcoa_results.qza")
bcmat <- read_qza("diversityMetrics/bray_curtis_distance_matrix.qza")

shannon <- shannon$data %>% rownames_to_column("SampleID")

metadata <- merge(x=metadata, y=shannon,by='SampleID') #combine shannon and metadata

p1 <- ggplot(metadata,aes(x=siteanimalhealth,y=shannon_entropy,fill=siteanimalhealth)) + #plotting site/animal health vs shannon
  stat_summary(geom="bar", fun.data=mean_se, color="black") + #bar outline 
  geom_jitter(shape=21, width=0.1, height=0)+ #arrange points
  coord_cartesian(ylim=c(2,8.5)) + # adjust y-axis
  facet_grid(~site) + #make seperate plots for each site
  theme_q2r() +
  theme(legend.position="none")+
  scale_fill_npg()+
  xlab("Site Animal Health") +
  ylab("Shannon Diversity") +
  ggsave("Shannon_barplot.pdf", height=3, width=5, device="pdf")
print(p1)

vec <- uwunifrac$data$Vectors
pc <- select(vec, SampleID, PC1, PC2)
metadata <- merge(x=metadata, y=pc, by='SampleID')

p2 <- ggplot(metadata, aes(x=PC1, y=PC2, color=siteanimalhealth, size=shannon_entropy)) +
  geom_point(alpha=0.7) + #alpha controls transparency and helps when points are overlapping
  theme_q2r()+
  scale_color_npg(name="Site Animal Health")+
  scale_size_continuous(name="Shannon Diversity") +
  ggsave("Shannon_unifrac.pdf", height=3, width=5, device="pdf")
print(p2)

p3 <- ggplot(metadata,aes(x=PC1,y=PC2, color=siteanimalhealth)) + #plotting site/animal health vs shannon
  geom_point(alpha=0.7) +
  facet_rep_wrap(~site) + #make seperate plots for each site
  coord_fixed(1) +
  theme_q2r() +
  theme(legend.position="none")+
  scale_color_npg(name="Site Animal Health:")+
  xlab("PC1") +
  ylab("PC2") +
  theme(legend.position="right") +
  ggsave("unifracSite.pdf", height=3, width=5, device="pdf")
print(p3)
