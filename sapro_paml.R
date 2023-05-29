library(tidyverse)
require(dplyr)
require(reshape2)
require(ggpubr)

##############################

setwd("~/Desktop/2023.01.20_Sapro/dNdS_dIdS/")

dat <- read.delim("2023.01.24_allPrank85_cga_paml_v2.txt")
di_dat <- read.delim("pairwiseDifferences.txt")
colnames(di_dat)

# filter to exclude when LWL85m_dN == None
dat <- filter(dat,LWL85m_dN != "None")

# select only relevant data from dat
dat_mini <- select(dat,"isolate","comparedTo_iso","YN00_dN", "YN00_dS")
colnames(dat_mini)

# merge pairwise differences of IGR data
dat_di_mini <- left_join(dat_mini, di_dat, by = c("isolate"="iso1","comparedTo_iso"="iso2"))

# build lists containing iso names of each group 
clade1 <- read.csv("../metadata/clade1_isolates.csv",header=F)$V1
clade2 <- read.csv("../metadata/clade2_isolates.csv",header=F)$V1

# calc dn/ds 
dat_di_mini$YN00_dNdS <- dat_di_mini$YN00_dN / dat_di_mini$YN00_dS

# calc di/ds by dividing the number of pairwise differences by the length of 
# the alignment, and then divide by the synonymous variation estimates 
# from paml (dS) 
# IGR aln length for sapro is 261952
dat_di_mini$dids <- (dat_di_mini$numSnps / 261952) / dat_di_mini$YN00_dS
colnames(dat_di_mini)

within_clade1 <- filter(dat_di_mini, isolate %in% clade1 & comparedTo_iso %in% clade1)
within_clade2 <- filter(dat_di_mini, isolate %in% clade2 & comparedTo_iso %in% clade2)
between_clades <- filter(dat_di_mini, isolate %in% clade1 & comparedTo_iso %in% clade2 | isolate %in% clade2 & comparedTo_iso %in% clade1)

# plot with comparison to SNAP piNpiS values
sapro <- rbind(data.frame(Group = "Clade 1", piNpiS = within_clade1$YN00_dNdS, method = "dnds"),
             data.frame(Group = "Clade 1", piNpiS = within_clade1$dids, method = "dids"),
             data.frame(Group = "Clade 2", piNpiS = within_clade2$YN00_dNdS, method = "dnds"),
             data.frame(Group = "Clade 2", piNpiS = within_clade2$dids, method = "dids"),
             data.frame(Group = "between clades", piNpiS = between_clades$YN00_dNdS, method = "dnds"),
             data.frame(Group = "between clades", piNpiS = between_clades$dids, method = "dids")
)

p <- ggplot(sapro, aes(Group, piNpiS,color = method)) + geom_boxplot() + theme_bw() + theme(text = element_text(size=15)) +labs( y = "value",x = "")+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_log10()
p

ggplot(sapro[sapro$method == "dnds",],aes(x=Group,y=piNpiS)) + geom_boxplot() + theme_bw() +
  ylim(0,0.5) + ylab("dN/dS")

ggplot(sapro[sapro$method == "dids",],aes(x=Group,y=piNpiS)) + geom_boxplot() + theme_bw() +
  ylab("dI/dS") + scale_y_log10()

##############################################################################

source.colors <- c("#F89441FF","#F0F921FF","dodgerblue","#0D0887FF","magenta")

iso <- read.csv("../metadata/2023.05.16_metadata.csv",header=T)

animal.iso <- iso$Isolate[iso$Source.general == "Animal"]
human.iso <- iso$Isolate[iso$Source.general == "Human"]
food.iso <- iso$Isolate[iso$Source.general == "Food"]
nenvironment.iso <- iso$Isolate[iso$Source.general == "Built environment"]
benvironment.iso <- iso$Isolate[iso$Source.general == "Natural environment"]

#dat_di_mini <- dat_di_mini %>% mutate(source = case_when((dat_di_mini$isolate %in% animal.iso & dat_di_mini$comparedTo_iso %in% animal.iso) ~ "animal",
                                                         #isolate %in% human.iso & comparedTo_iso %in% human.iso ~ "human",
                                                         #isolate %in% food.iso & comparedTo_iso %in% food.iso ~ "food",
                                                         #isolate %in% environment.iso & comparedTo_iso %in% environment.iso ~ "environment",
#                                                         .default = "between"))

within_animal <- filter(dat_di_mini, isolate %in% animal.iso & comparedTo_iso %in% animal.iso)
within_human <- filter(dat_di_mini, isolate %in% human.iso & comparedTo_iso %in% human.iso)
within_food <- filter(dat_di_mini, isolate %in% food.iso & comparedTo_iso %in% food.iso)
within_benvironment <- filter(dat_di_mini, isolate %in% benvironment.iso & comparedTo_iso %in% benvironment.iso)
within_nenvironment <- filter(dat_di_mini, isolate %in% nenvironment.iso & comparedTo_iso %in% nenvironment.iso)

source_dnds <- rbind(data.frame(Group = "human", dNdS = within_human$YN00_dNdS, dIdS = within_human$dids, dS = within_human$YN00_dS),
                     data.frame(Group = "animal", dNdS = within_animal$YN00_dNdS, dIdS = within_animal$dids, dS = within_animal$YN00_dS),
                     data.frame(Group = "food", dNdS = within_food$YN00_dNdS, dIdS = within_food$dids, dS = within_food$YN00_dS),
                     data.frame(Group = "built-environment", dNdS = within_benvironment$YN00_dNdS, dIdS = within_benvironment$dids, dS = within_benvironment$YN00_dS),
                     data.frame(Group = "natural-environment", dNdS = within_nenvironment$YN00_dNdS, dIdS = within_nenvironment$dids, dS = within_nenvironment$YN00_dS))



dnds <- ggplot(source_dnds) + geom_boxplot(aes(x=Group,y=dNdS,fill=Group)) + theme_bw() + scale_y_log10() +
  scale_fill_manual(values=source.colors)
dnds

dnds_stats <- compare_means(data=source_dnds,dNdS~Group,p.adjust.method = "bonferroni")

dids <- ggplot(source_dnds) + geom_boxplot(aes(x=Group,y=dIdS,fill=Group)) + theme_bw() + scale_y_log10() +
  scale_fill_manual(values=source.colors)
dids

dids_stats <- compare_means(data=source_dnds,dIdS~Group,p.adjust.method = "bonferroni")

ExportPlot(dnds,"Figures/dnds_source",width=6,height=4)
ExportPlot(dids,"Figures/dids_source",width=6,height=4)

source_dnds.melt <- melt(source_dnds,measure.vars=c("dNdS","dIdS"),id.vars = "Group",variable.name = "Type",value.name = "Value")
supp.labs <- c("dN/dS","dI/dS")
names(supp.labs) <- c("dNdS","dIdS")

both <- ggplot(source_dnds.melt) + geom_boxplot(aes(x=Group,y=Value,fill=Group)) + theme_bw() + scale_y_log10() +
  scale_fill_manual(name = "Isolation source",values=source.colors,
                    breaks=c("animal","built-environment","food","human","natural-environment"),
                    labels=c("Animal","Built environment","Food","Human","Natural environment")) + 
  facet_wrap(~Type, scales="free_y",labeller=labeller(Type=supp.labs))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())
both

ExportPlot(both,"Figures/dNdS_dIdS_bySource",width=10,height=4)

ExportPlot <- function(gplot, filename, width=2, height=1.5) {
  # Export plot in PDF and EPS.
  # Notice that A4: width=11.69, height=8.27
  ggsave(paste(filename, '.pdf', sep=""), gplot, width = width, height = height)
  postscript(file = paste(filename, '.eps', sep=""), width = width, height = height, family = "sans")
  print(gplot)
  dev.off()
  png(file = paste(filename, '_.png', sep=""), width = width * 100, height = height * 100)
  print(gplot)
  dev.off()
}


