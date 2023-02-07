# DMR annotator?
# I hardly know her!

# load in results
load("dmrseq.hc.rdata")

# filter for significant regions (pval < 0.05 in this case)
sigRegions <- regions[which(regions$pval<0.05),]

# annotate to mouse genes/features
dmrs2 <- sigRegions
peaks <- annotatePeak(dmrs2,tssRegion=c(-2000,2000),TxDb=txdb,annoDb="org.Mm.eg.db")
peaks
peaks <- as.data.frame(peaks)
dmrs2 <- as.data.frame(dmrs2)
dmrs <- cbind(dmrs2,peaks$distanceToTSS,peaks$ENSEMBL,peaks$annotation,peaks$SYMBOL,peaks$GENENAME)
names(dmrs) <- c(names(dmrs)[1:14],"distanceToTSS","ENSEMBL","annotation","SYMBOL","GENENAME")

# remove gene symbols from intergenic regions
dmrs[which(dmrs$annotation=="Distal Intergenic"),c("ENSEMBL","SYMBOL","GENENAME")] <- NA

# get the annotations out
sigRegions <- dmrs
