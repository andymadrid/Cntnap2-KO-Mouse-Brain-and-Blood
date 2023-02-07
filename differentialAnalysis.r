# Run differential methylation analysis with dmrseq

# load packages
library(DSS)
library(dmrseq)
library(ChIPseeker)
library(BiocParallel)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(DOSE)
library(GenomicRanges)
library(fdrtool)
library(ggplot2)
library(ggfortify)
library(methylCC)
library(RColorBrewer)
library(liftOver)
library(rtracklayer)
library(DOSE)
library(DMRichR)
library(viridis)
library(enrichplot)
library(wesanderson)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
mmAnno <- getAnnot("mm10")

# read in data
infile <- list.files(pattern="cov$")

# generate matrices
bs <- read.bismark(files = infile,rmZeroCov=TRUE,strandCollapse=T,verbose=T)
x <- colnames(bs)
x <- gsub(".CpG_report.merged_CpG_evidence.cov","",x)
colnames(bs) <- x
pData(bs)$Genotype <- rep(c("HOMO","WT"),c(6,6))
pData(bs)$Tissue <- rep(c("Blood","Hippo"),times=6)
pData(bs)$Condition <- paste0(pData(bs)$Genotype,"_",pData(bs)$Tissue)
pData(bs)$Mouse <- rep(c("KO1","KO2","KO3","WT1","WT2","WT3"),c(2,2,2,2,2,2))

# remove chrY from dataset
bs <- chrSelectBSseq(bs,seqnames=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX"))

# separate the groups
bloodSamples <- c(1,3,5,7,9,11)
hippoSamples <- c(2,4,6,8,10,12)
homoSamples <- c(1,3,5,2,4,6)
wtSamples <- c(7,9,11,8,10,12)

# differential analysis

### Blood samples
bs2 <- bs[,bloodSamples]
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs2, type="Cov")==0) == 0)
bs.filtered <- bs2[loci.idx,]
cores <- 40
regions <- dmrseq(bs.filtered,testCovariate="Condition",maxPerms=10,minNumRegion=5,cutoff=0.05,BPPARAM=BiocParallel::MulticoreParam(workers=cores))
save(regions,bs.filtered,file="dmrseq.blood.rdata")

### Hippocampus samples
bs2 <- bs[,hippoSamples]
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs2, type="Cov")==0) == 0)
bs.filtered <- bs2[loci.idx,]
cores <- 40
regions <- dmrseq(bs.filtered,testCovariate="Condition",maxPerms=10,minNumRegion=5,cutoff=0.05,BPPARAM=BiocParallel::MulticoreParam(workers=cores))
save(regions,bs.filtered,file="dmrseq.hc.rdata")
