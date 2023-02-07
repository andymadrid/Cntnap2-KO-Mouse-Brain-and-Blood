# Make a volcano plot from dmrseq results
# Things is heating up in here

# load in data
load("dmrseq.hc.rdata")
x <- as.data.frame(regions)

# parse for significant regions
sigUp <- x[which(x$beta<0 & x$pval<0.05),]
sigDown <- x[which(x$beta>0 & x$pval<0.05),]
notSig <- x[which(x$pval>0.05),]

# add in some pretty colours
sigUp$color <- "indianred"
sigDown$color <- "deepskyblue3"
notSig$color <- "ivory3"

# bring it all back together
x <- rbind(sigUp,sigDown,notSig)

# plot, baby, plot!
plot((-1*x$beta),-log10(x$pval),col=x$color,pch=16,xlab="Beta",ylab="-log10(p)")
