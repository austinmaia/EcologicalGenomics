library(qiime2R)
library(tidyverse)
library(gplots)
library("ggsci")
library(lemon)
library(wesanderson)
library(patchwork)

shannon <- read_qza("diversityMetrics/shannon_vector.qza")
metadata <- read.table("diversityMetrics/pyc_manifest", header=T)
uwunifrac<-read_qza("diversityMetrics/unweighted_unifrac_pcoa_results.qza")

shannon <- shannon$data %>% rownames_to_column("SampleID")

metadata <- merge(x=metadata, y=shannon,by='SampleID') #combine shannon and metadata


p1 <- ggplot(metadata,aes(x=siteanimalhealth,y=shannon_entropy,fill=siteanimalhealth)) + #plotting site/animal health vs shannon
  stat_summary(geom="bar", fun.data=mean_se) + #bar outline 
  geom_jitter(shape=21, alpha=0.7, width=0.1, height=0, color="transparent")+ #arrange points
  coord_cartesian(ylim=c(2,8.5)) + # adjust y-axis
  theme_q2r()+
  facet_grid(~site) + #make seperate plots for each site
  scale_fill_manual(values=wes_palette(n=3, "FantasticFox1"),name="Site and Animal Health")  + #change to color-blind friendly
  xlab("Site Animal Health") +
  ylab("Shannon Diversity") +
  ggsave("Shannon_barplot.jpg", height=3, width=5, device="jpeg")
print(p1)

vec <- uwunifrac$data$Vectors 
pc <- select(vec, SampleID, PC1, PC2) #get PCA axes
metadata <- merge(x=metadata, y=pc, by='SampleID') #add PCA axes to metadata
proportionEx <- uwunifrac$data$ProportionExplained

p2 <- ggplot(metadata, aes(x=PC1, y=PC2, color=siteanimalhealth, size=shannon_entropy)) + #PCA plot with unifrac distance, putting shannon diversity as variable for size
  theme_q2r()+
  geom_point(alpha=0.7, aes(color = siteanimalhealth)) +
  scale_color_manual(values=wes_palette(n=3, "FantasticFox1"), name="Site and Animal Health")+ 
  xlab("PC1 (12.1%)") +
  ylab("PC2 (7.3%)") +
  stat_ellipse(type='t',size =0.25) +
  scale_size_continuous(name="Shannon Diversity") +
  ggsave("Shannon_unifrac.jpg", height=3, width=5, device="jpeg")
print(p2)

p3 <- ggplot(metadata,aes(x=PC1,y=PC2, color=siteanimalhealth)) + #plotting shannon index for each site
  geom_point(alpha=0.7) +
  facet_rep_wrap(~site) + #make seperate plots for each site
  coord_fixed(1) +
  theme_q2r() +
  theme(legend.position="none")+
  scale_color_manual(values=wes_palette(n=3, "FantasticFox1"), name="Site and Animal Health")+
  xlab("PC1 (12.1%)") +
  ylab("PC2 (7.3%)") +
  theme(legend.position="right") +
  ggsave("unifracSite.jpg", height=3, width=5, device="jpeg")
print(p3)


p4 <- ggplot(metadata,aes(x=siteanimalhealth,y=shannon_entropy,fill=siteanimalhealth)) + #plotting site/animal health vs shannon
  stat_summary(geom="bar", fun.data=mean_se) + #bar outline 
  geom_jitter(shape=21, alpha=0.7, width=0.1, height=0, color="transparent")+ #arrange points
  coord_cartesian(ylim=c(2,8.5)) + # adjust y-axis
  theme_q2r()+
  scale_fill_manual(values=wes_palette(n=3, "FantasticFox1"))  + #change to color-blind friendly
  theme(legend.position="none")+ #add legend
  xlab("Site Animal Health") +
  ylab("Shannon Diversity") +
  ggsave("Shannon_barplot_full.jpg", height=3, width=5, device="jpeg")
print(p4)

shannonplots <- p4 + p1 + plot_layout(ncol=2,widths=c(1,3)) + plot_annotation(tag_levels = 'A')
ggsave("shannon_combo.jpg", shannonplots, device="jpeg")

# T-Tests
sh_hh = subset(metadata, siteanimalhealth=="HH", select=shannon_entropy)
sh_sh = subset(metadata, siteanimalhealth=="SH", select=shannon_entropy)
sh_ss = subset(metadata, siteanimalhealth=="SS", select=shannon_entropy)
sh_h = subset(metadata, animalhealth = 'healthy', select = shannon_entropy)

t.test(sh_sh, sh_hh) # t-test between healthy at naive vs impact
t.test(sh_ss, sh_hh) # t-test between naive healthy and sick
t.test(sh_ss, sh_sh) # t-test between sick and healthy at impacted
t.test(sh_ss, sh_h) # t-test between all healthy and all sick

sh_h12 = subset(metadata, animalhealth == 'healthy' & site == 'site_12', select=shannon_entropy)
sh_s12 = subset(metadata, animalhealth == 'sick' & site == 'site_12', select=shannon_entropy)

t.test(sh_h12, sh_s12) #t-test between healthy and sick at site 12

sh_h13 = subset(metadata, animalhealth == 'healthy' & site == 'site_13', select=shannon_entropy)
sh_s13 = subset(metadata, animalhealth == 'sick' & site == 'site_13', select=shannon_entropy)

t.test(sh_h13, sh_s13) #t-test between healthy and sick at site 13

sh_h15 = subset(metadata, animalhealth == 'healthy' & site == 'site_15', select=shannon_entropy)
sh_s15 = subset(metadata, animalhealth == 'sick' & site == 'site_15', select=shannon_entropy)

t.test(sh_h15, sh_s15) #t-test between healthy and sick at site 15
