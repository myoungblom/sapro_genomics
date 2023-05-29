require(ggplot2)
require(viridis)
require(dplyr)
require(ggpubr)
require(scales)
library(gridExtra)

setwd("~/Desktop/2023.01.20_Sapro/clonalframeml/")

#### stats ####
# fragment lengths
# avg frags per isolate

clade.colors <- viridis(n=2,option="A",begin=0.3,end=0.8)

clade1 <- read.delim("DataFiles/2023.02.13_clade1Prank85_cga_CFML.importation_status.txt")
clade2 <- read.delim("DataFiles/2023.02.13_clade2Prank85_cga_CFML.importation_status.txt")

all.clades <- rbind(data.frame(beg=clade1$Beg,end=clade1$End,node=clade1$Node,clade="Clade 1"),
                    data.frame(beg=clade2$Beg,end=clade2$End,node=clade2$Node,clade="Clade 2"))
all.clades$length <- all.clades$end - all.clades$beg

length <- ggplot(all.clades,aes(x=clade,y=length)) + geom_boxplot(aes(fill=clade)) + theme_minimal() +
  scale_y_log10() + xlab(NULL) + ylab("Fragment length (bp)")+
  scale_fill_manual(name=NULL,breaks = c("Clade 1","Clade 2"), values = clade.colors)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=10))+
  stat_compare_means(label="p.signif")
length

counts1 <- read.delim("DataFiles/clade1_frags_per_iso.txt")
counts2 <- read.delim("DataFiles/clade2_frags_per_iso.txt")

counts <- rbind(data.frame(iso=counts1$Strain,count=counts1$NumFrags,clade="Clade 1"),
                data.frame(iso=counts2$Strain,count=counts2$NumFrags,clade="Clade 2"))

count <- ggplot(counts,aes(x=clade,y=count)) + geom_boxplot(aes(fill=clade)) + theme_minimal() +
  xlab(NULL) + ylab("Fragments per strain")+
  scale_fill_manual(name=NULL,breaks = c("Clade 1","Clade 2"), values = clade.colors)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=10))+
  stat_compare_means(label="p.signif")
count

full <- ggarrange(count,length)
full

ExportPlot(full,"Figures/CFML_Clades_v2",width=8,height=6)

### Clade 1 -- plot with tree ###
clade1 <- read.delim("DataFiles/2023.02.13_clade1Prank85_cga_CFML.importation_status.txt")

# read tree
tre <- read.tree("../trees/CORE/2023.02.13_clade1Prank85_cga_CFML.labelled_tree.newick")
r <- ggtree(tre,layout = "rectangular", size = 0.5)

# order to plot strains (on y-axis, matches tree)
isoOrder <- get_taxa_name(r)

currY <- 0
id <- vector()
X <- vector() # x-axis position of recombination block
Y <- vector() # y-axis position of recombination bloc

#building data for recent
scratchdf <- data.frame() #place we're gonna dump stuff

currY = 645 # number of isolates in sample-1
for(iso in 1:length(isoOrder)){  #loop through isolates in lineage
  strain.name <- toString(isoOrder[iso])
  # recombination events where recipient strain is "iso"
  scratchdf <- clade1[clade1$Node == strain.name,]
  #loop through fragments to plot for isolate assuming there are any
  if(nrow(scratchdf)>0){
    for(i in (1:nrow(scratchdf))){
      id <- c(id,rep(strain.name,each=5))
      X <- c(X,scratchdf[i,2],scratchdf[i,3],scratchdf[i,3],scratchdf[i,2],scratchdf[i,2]) # start/stop of frag
      Y <- c(Y,currY,currY,currY+1,currY+1,currY) # y-position of frag
    }}
  currY <- currY-1 #decrement row
}

# make data frames for plotting
datapoly <- data.frame(X=X, Y=Y, strain=id)

#background goes in before recomb fragments
options(scipen = 999)
recentplot <- ggplot() + 
  geom_polygon(data=datapoly, aes(x=X, y=Y, group=strain),fill="black") + 
  xlab("Core Genome Position")+
  theme_classic() +
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = NA, colour = NA))

recentplot

full <- ggarrange(r,recentplot,widths=c(2,5))
full

ExportPlot(recentplot,"Figures/Clade1_CFMLfrags_v3",width=12,height=8)
ExportPlot(full,"Figures/Clade1_CFMLfrags_wTree_v3",width=12,height=16)

### Clade 2 -- plot with tree ###
clade2 <- read.delim("DataFiles/2023.02.13_clade2Prank85_cga_CFML.importation_status.txt")

# read tree
tre <- read.tree("../trees/CORE/2023.02.13_clade2Prank85_cga_CFML.labelled_tree.newick")
r <- ggtree(tre,layout = "rectangular", size = 0.5)

# order to plot strains (on y-axis, matches tree)
isoOrder <- get_taxa_name(r)

currY <- 0
id <- vector()
X <- vector() # x-axis position of recombination block
Y <- vector() # y-axis position of recombination bloc

#building data for recent
scratchdf <- data.frame() #place we're gonna dump stuff

currY = 133 # number of isolates in sample-1
for(iso in 1:length(isoOrder)){  #loop through isolates in lineage
  strain.name <- toString(isoOrder[iso])
  # recombination events where recipient strain is "iso"
  scratchdf <- clade2[clade2$Node == strain.name,]
  #loop through fragments to plot for isolate assuming there are any
  if(nrow(scratchdf)>0){
    for(i in (1:nrow(scratchdf))){
      id <- c(id,rep(strain.name,each=5))
      X <- c(X,scratchdf[i,2],scratchdf[i,3],scratchdf[i,3],scratchdf[i,2],scratchdf[i,2]) # start/stop of frag
      Y <- c(Y,currY,currY,currY+1,currY+1,currY) # y-position of frag
    }}
  currY <- currY-1 #decrement row
}

# make data frames for plotting
datapoly <- data.frame(X=X, Y=Y, strain=id)

options(scipen = 999)
recentplot <- ggplot() + 
  geom_polygon(data=datapoly, aes(x=X, y=Y, group=strain),fill="black") + 
  xlab("Core Genome Position")+
  xlim(0,3553895) + 
  theme_classic() +
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = NA, colour = NA))

recentplot

full <- ggarrange(r,recentplot,widths=c(1,3))
full

ExportPlot(recentplot,"Figures/Clade2_CFMLfrags_v2",width=12,height=8)
ExportPlot(full,"Figures/Clade2_CFMLfrags_wTree",width=12,height=8)

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
