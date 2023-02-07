# Makes a heatmap of mean DMR methylation values either per sample or by group

# load in pacakges
library(gplots)
library(RColorBrewer)
library(viridis)
methMean <- c()

# load in data
load("dmrseq.hc.rdata")
sigRegions <- regions[which(regions$pval<0.05),]

# get methylation values from available data
methSm <- bsseq::getMeth(bs.filtered,type="raw") # probably should used smoothed values, but this works either way
chrs <- as.data.frame(seqnames(bs.filtered))
pos <- as.data.frame(start(bs.filtered))
methWork <- cbind(chrs,pos,as.data.frame(methSm))
colnames(methWork) <- c("chr","pos",colnames(bs.filtered))
meth.gr <- with(methWork,GRanges(chr,IRanges(pos,pos+1)))
pal <- colorpanel(100,"dodgerblue3","goldenrod1","firebrick3")

# filter methylation matrix for only CpGs in DMRs
dmrs.gr <- sigRegions
for (i in 1:length(dmrs.gr)) {
  x <- as.data.frame(findOverlaps(meth.gr,dmrs.gr[i]))
  methWork.subset <- methWork[as.numeric(x$queryHits),3:ncol(methWork)]
  dmrMean <- colMeans(methWork.subset)
  methMean <- rbind(methMean,dmrMean)
  prog <- paste0("Done with DMR ",i,"\n")
  prog <- noquote(prog)
  cat(prog)
  }
#colnames(methMean) <- c("KO1","KO2","KO3","WT1","WT2","WT3")

# plot results for "by sample" instead of group means
heatmap.2(as.matrix(methMean),col=pal,trace="none",dendrogram="column",labRow=F,key.title=NULL,key.xlab="Methylation Level",density.info="none")

# Move on to group means
wBlood <- c()
hBlood <- c()
wHC <- c()
hHC <- c()
for (i in 1:nrow(methMean)) {
x <- mean(methMean[i,c(1,3,5)])
wBlood <- rbind(wBlood,x)
x <- mean(methMean[i,c(7,9,11)])
hBlood <- rbind(hBlood,x)
x <- mean(methMean[i,c(2,4,6)])
wHC <- rbind(wHC,x)
x <- mean(methMean[i,c(8,10,12)])
hHC <- rbind(hHC,x)
}
groupMeans <- cbind(wBlood,wHC,hBlood,hHC)
colnames(groupMeans) <- c("wtBlood","wtHC","koBlood","koHC")

# plot results for "by group" mean values
heatmap.2(as.matrix(groupMeans),col=pal,trace="none",dendrogram="column",labRow=F,key.title=NULL,key.xlab="Methylation Level",density.info="none")
