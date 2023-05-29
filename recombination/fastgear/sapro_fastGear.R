library(ggplot2)
library(gtable)
library(gridExtra)
library(viridis)
library(wesanderson)
library(tidyverse)
library(stringr)

# length of core genome alignment
alnLength <- 942474
  
# recent recombination events
recent <- read.table("recombinations_recent.txt",header = TRUE)
# lineage & cluster designations for each strain
lineage <- read.table("lineage_information.txt", header = TRUE)
# order to plot strains (on y-axis, matches tree)
isoOrder <- read.table("../allPrank85_RAxML_isolateOrder.txt", header = FALSE) # this might need to be reversed to plot in the correct order
isoOrder$Name <- isoOrder$V1
# re-order lineage info by strain order
isoReOrder <- left_join(isoOrder,lineage,by="Name")

# clade 1, clade 2
clade.colors <- c(viridis(n=2,option="A",begin=0.3,end=0.8),"white"); clade.colors

### recent recombination events ###

llist <- as.data.frame(table(lineage$Lineage)) #get number of isolates per lineage
isoframe <- list()
for(idx in llist[[1]]){
  isoframe[[idx]] <- lineage[lineage$Lineage == idx,c(1)]
}#build a list of the isolates in each lineage (by isolate index)

currY <- 0
X <- vector() # x-axis position of recombination block
Y <- vector() # y-axis position of recombination block
rcolor <- vector() # color to fill recombination block
id <- vector() # strain ID for recombination block (donating strain?)
bgX <- vector() # x-axis position of background block (recipient strain)
bgY <- vector() # y-axis position of background block
bgid <- vector() # strain ID of background (recipient) strain
bgcolor <- vector() # color to fill background block

# individual recombination events
recipients <- c()
dlin <- c()
rlin <- c()
starts <- c()
stops <- c()

#building data for recent
scratchdf <- data.frame() #place we're gonna dump stuff

currY = 779 # number of isolates in sample
for(iso in (1:length(isoOrder$V1))){   #loop through isolates in lineage
  strain.lineage <- as.integer(isoReOrder$Lineage[iso])
  strain.index <- isoReOrder$StrainIndex[iso]
  strain.name <- toString(isoReOrder$Name[iso])
  
  # each isolate gets background block 1xalnLength
  bgid <- c(bgid,rep(strain.name,each=5))
  bgY <- c(bgY,currY,currY,currY+1,currY+1,currY)
  bgX <- c(bgX,0,alnLength,alnLength,0,0)
  bgcolor <- c(bgcolor,rep(clade.colors[strain.lineage],each=5)) # number of colors = number of lineages + 1
  
  # recombination events where recipient strain is "iso"
  scratchdf <- recent[recent$RecipientStrain == strain.index,c(1:6)]
  #loop through fragments to plot for isolate assuming there are any
  if(nrow(scratchdf)>0){
  for(i in (1:nrow(scratchdf))){
    #if(scratchdf[i,5]>1){ # only include fragments with log(BF) > 1
      recp.strain <- toString(scratchdf[i,6])
      rlin <- c(rlin,strain.lineage)
      recipients <- c(recipients,recp.strain)
      donor.lin <- as.integer(scratchdf[i,3])
      dlin <- c(dlin,donor.lin)
      starts <- c(starts,scratchdf[i,1])
      stops <- c(stops,scratchdf[i,2])
      id <- c(id,rep(recp.strain,each=5)) # id of recipient strain
      X <- c(X,scratchdf[i,1],scratchdf[i,2],scratchdf[i,2],scratchdf[i,1],scratchdf[i,1]) # start/stop of frag
      Y <- c(Y,currY,currY,currY+1,currY+1,currY) # y-position of frag
      rcolor <- c(rcolor,rep(clade.colors[donor.lin],each=5))
    }}
  currY <- currY-1 #decrement row
}

# data frame of individual recombination events
recombs <- data.frame(r.strain=recipients,r.lin=rlin,d.lin=dlin,start=starts,stop=stops)

# make data frames for plotting
datapoly <- data.frame(ids=id, value=id, X=X, Y=Y, color=rcolor)
nrow(datapoly)
length(rcolor)

rbgpoly <- data.frame(ids=bgid, value=bgid, X=bgX, Y=bgY, color=bgcolor)
nrow(rbgpoly)
length(bgcolor)

#c2bg <- data.frame(x=c(0,alnLength,alnLength,0,0),y=c(0,0,134,134,0))
#c1bg <- data.frame(x=c(0,alnLength,alnLength,0,0),y=c(134,134,780,780,134))

#background goes in before recomb fragments
recentplot <- ggplot() + 
  geom_polygon(data=rbgpoly, aes(x=bgX, y=bgY, group=ids),fill=bgcolor) +
  #geom_polygon(data=c2bg,aes(x=x,y=y),fill=clade.colors[1])+
  #geom_polygon(data=c1bg,aes(x=x,y=y),fill=clade.colors[2])+
  geom_polygon(data=datapoly, aes(x=X, y=Y, group=ids),fill=rcolor) + 
  scale_fill_manual(values=clade.colors, name = "Clades")+
  xlab("Core Genome Position")+
  xlim(0,alnLength) + 
  theme_classic(base_size = 10) +
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = NA, colour = NA),
        plot.margin = margin(0, 0, 0, 0))

recentplot

#ExportPlot(recentplot,"../Figures/recent_recombinations",width=12,height=8)

tre <- read.raxml("../../trees/CORE/2023.01.24_allPrank85_cga_RAxML.newick")
r <- ggtree(tre,layout = "rectangular", size = 0.5)

full <- ggarrange(r,recentplot,widths=c(1,3))
full

ExportPlot(full,"../Figures/recent_recombinations_wTree_v3",width=12,height=8)

### getting color legend ###
color_dat <- data.frame(clade=c("Clade 1","Clade 2"),height=c(10,10),
                        color=clade.colors[-3])

testplot <- ggplot(data=color_dat, aes(x=clade,y=height))+
  geom_point(aes(color=clade)) +
  scale_color_manual( values = clade.colors, name = NULL, labels = c("Clade 1","Clade 2"))+
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size=5)))

testplot

library(grid)
leg <- cowplot::get_legend(testplot)
grid.newpage()
grid.draw(leg)


#function to export plot
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
