library(viridis)
library(ggplot2)

####
# Takes in 2 gene_presence_absence.Rtab files from Roary and calculates rarefaction
# and accumulation curves while sub-setting to the size of the smaller group.
####


#read the FIRST input file - the one to be subsampled. 
#rows corespond to genes, columns corespond to isolates
data1 <- read.delim("sapro_clade1_gpa.Rtab", sep="\t",row.names = 1)
ncol(data1)

#read the SECOND input file - the one that IS NOT subsampled
data2 <- read.delim("sapro_clade2_gpa.Rtab", sep="\t",row.names = 1)
ncol(data2)

#size of smaller dataset, size to subsample to
sub.n <- ncol(data2)

#number of iterations
n.iter <- 100

#data 1 accumulation curve
data1.t <- t(data1)
subset <- data1.t[sample(nrow(data1.t),sub.n),]
n.rows <- nrow(subset) #isolates
n.cols <- ncol(subset) #genes
core_results <- matrix(data=NA,nrow=n.rows,ncol=n.iter) #rows = isolates, columns = iterations
presence_line <- matrix(data=0,nrow=1,ncol=n.cols) #columns = genes, single row
for(times in 1: n.iter){ #for each iteration
  for (i in 1: n.rows){ #for each isolate
    sum <- 0; #total genes for each isolate
    if(i==1){ #if its the first isolate
      for(j in 1: n.cols){ #for each gene
        presence_line[1,j]<-subset[i,j] #add 0/1 for each gene
        sum <- sum + subset[i,j]; #add 0/1 to total for each isolate
      }
      core_results[i,times] <- sum #add isolate gene total for each iteration
      next;
    }
    for(j in 1: n.cols){ #(if its not the first isolate) for each gene
      if(presence_line[1,j] || subset[i,j]){ #if either the first isolate or the current
  #                                          isolate has the gene
        presence_line[1,j] <- 1; #add 1
      }
      else{presence_line[1,j] <- 0;} #else add 0
    }
    for(j in 1: n.cols){ #(if its not the first isolate) for each gene
      sum <- sum + presence_line[1,j]; #add up all genes
    }
    core_results[i,times] <- sum; #for this number of isolates, in this iteration, add gene total
  }
  subset <- data1.t[sample(nrow(data1.t),sub.n), ] #shuffle data
}
data1.acc <- as.data.frame(t(core_results))

#data 1 rarefaction curve
data1.t <- t(data1)
subset <- data1.t[sample(nrow(data1.t),sub.n),]
n.rows <- nrow(subset) #isolates
n.cols <- ncol(subset) #genes
core_results <- matrix(data=NA,nrow=n.rows,ncol=n.iter) #rows = isolates, columns = iterations
presence_line <- matrix(data=0,nrow=1,ncol=n.cols) #columns = genes, single row
for(times in 1:n.iter){ #for each iteration
  for (i in 1: n.rows){ #for each isolate
    sum <- 0; #total genes for each isolate?
    if(i==1){ #if its the first isolate
      for(j in 1: n.cols){
        presence_line[1,j]<-subset[i,j]
        sum <- sum + subset[i,j]; 
      }
      core_results[i,times] <- sum
      next;
    }
    for(j in 1: n.cols){
      if(presence_line[1,j] && subset[i,j]){
        presence_line[1,j] <- 1;
      }
      else{presence_line[1,j] <- 0;}
    }
    for(j in 1: n.cols){
      sum <- sum + presence_line[1,j]; 
    }
    core_results[i,times] <- sum;
  }
  subset <- data1.t[sample(nrow(data1.t),sub.n), ]
}
data1.rare <- as.data.frame(t(core_results))

#data 2 accumulation curve
data2.t <- t(data2)
subset <- data2.t[sample(nrow(data2.t),sub.n),]
n.rows <- nrow(subset)
n.cols <- ncol(subset)
core_results <- matrix(data=NA,nrow=n.rows,ncol=n.iter)
presence_line <- matrix(data=0,nrow=1,ncol=n.cols)
for(times in 1: n.iter){
  for (i in 1: n.rows){
    sum <- 0;
    if(i==1){
      for(j in 1: n.cols){
        presence_line[1,j]<-subset[i,j]
        sum <- sum + subset[i,j]; 
      }
      core_results[i,times] <- sum
      next;
    }
    for(j in 1: n.cols){
      if(presence_line[1,j] || subset[i,j]){
        presence_line[1,j] <- 1;
      }
      else{presence_line[1,j] <- 0;}
    }
    for(j in 1: n.cols){
      sum <- sum + presence_line[1,j]; 
    }
    core_results[i,times] <- sum;
  }
  subset <- data2.t[sample(nrow(data2.t),sub.n), ]
}
data2.acc <- as.data.frame(t(core_results))

#data 2 rarefaction curve
data2.t <- t(data2)
subset <- data2.t[sample(nrow(data2.t),sub.n),]
n.rows <- nrow(subset)
n.cols <- ncol(subset)
core_results <- matrix(data=NA,nrow=n.rows,ncol=n.iter)
presence_line <- matrix(data=0,nrow=1,ncol=n.cols)
for(times in 1:n.iter){
  for (i in 1: n.rows){
    sum <- 0;
    if(i==1){
      for(j in 1: n.cols){
        presence_line[1,j]<-subset[i,j]
        sum <- sum + subset[i,j]; 
      }
      core_results[i,times] <- sum
      next;
    }
    for(j in 1: n.cols){
      if(presence_line[1,j] && subset[i,j]){
        presence_line[1,j] <- 1;
      }
      else{presence_line[1,j] <- 0;}
    }
    for(j in 1: n.cols){
      sum <- sum + presence_line[1,j]; 
    }
    core_results[i,times] <- sum;
  }
  subset <- data2.t[sample(nrow(data2.t),sub.n), ]
}
data2.rare <- as.data.frame(t(core_results))

data1.c <- apply(data1.rare, 2, median)
data2.c <- apply(data2.rare, 2, median)
data1.t <- apply(data1.acc, 2, median)
data2.t <- apply(data2.acc, 2, median)


data1.genes <- data.frame( genes_to_genomes1 = c(data1.c,data1.t),
                            genomes1 = c(c(1:length(data1.c)),c(1:length(data1.c))),
                            Key = c(rep("Core genes",length(data1.c)), rep("Total genes",length(data1.t))) )
data2.genes <- data.frame( genes_to_genomes1 = c(data2.c,data2.t),
                            genomes1 = c(c(1:length(data2.c)),c(1:length(data2.c))),
                            Key = c(rep("Core genes",length(data2.c)), rep("Total genes",length(data2.t))) )

colors <- viridis(n=2,option="A",begin=0.3,end=0.8)

rarefaction <- ggplot() + 
  geom_line(data = data1.genes, aes(x = genomes1, y = genes_to_genomes1, group = Key, linetype=Key, color = "Clade 1"), linewidth = 1.5, alpha = 0.8) +
  geom_line(data = data2.genes, aes(x = genomes1, y = genes_to_genomes1, group = Key, linetype=Key, color = "Clade 2"), linewidth = 1.5, alpha = 0.8)+
  scale_color_manual(name=NULL,breaks = c("Clade 1","Clade 2"), values = colors)+
  theme_classic() +
  labs(linetype = "Genes")+
  xlab("\nNumber of genomes") +
  ylab("Number of genes\n")+ theme_bw(base_size = 10) +
  theme(legend.justification=c(0,1),legend.position=c(0.01,.99), legend.box = "vertical",
        legend.text = element_text(color = "grey27", size = 12), 
        axis.text = element_text(color = "grey27", face = "bold",size=12),
        axis.title = element_text(color = "grey27", face = "bold", size = 14))+
  guides( color = guide_legend(order = 1), linetype = guide_legend(order = 2))

rarefaction

ExportPlot(rarefaction, "Figures/sapro_rarefaction_median", width=6, height=8)

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
