# Overlap significance testing

# load some packages
library(regioneR)
library(BSgenome.Mmusculus.UCSC.mm10.masked)

# parse the genome for chromosomes of interest
genome <- BSgenome.Mmusculus.UCSC.mm10.masked
toKeep <- paste0("chr",c(1:19,"X"))
genome@user_seqnames <- setNames(toKeep,toKeep)
genome@seqinfo <- genome@seqinfo[toKeep]

# read in hippocampal results
load("dmrseq.hc.rdata")
sigRegions.hc <- regions[which(regions$pval<0.05),]
regions.hc <- regions


# read in blodo results
load("dmrseq.blood.rdata")
sigRegions.blood <- regions[which(regions$pval<0.05),]
regions.blood <- regions

# run the permutation test
permRes <- overlapPermTest(A=sigRegions.hc,B=sigRegions.blood,genome=genome,count.once=T,ntimes=10000)

# plot the results
pdf("permResults.pdf")
plot(permRes)
dev.off()
