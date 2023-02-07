# Make barplots of GO results obtained upstairs (aka another script)

# load pacakges
library(forcats)
library(gplots)

# set up the GO results
egoAll <- as.data.frame(egoAll)
egoHyper <- as.data.frame(egoHyper)
egoHypo <- as.data.frame(egoHypo)

# this keeps ggplot from changing the order to alphanumeric (ugh)
egoAll$Description <- factor(egoAll$Description,levels=egoAll$Description)
egoHyper$Description <- factor(egoHyper$Description,levels=egoHyper$Description)
egoHypo$Description <- factor(egoHypo$Description,levels=egoHypo$Description)

# plot it out
pdf("GO.all.pdf",height=5,width=12)
ggplot(egoAll[1:10,],aes(x=-log10(p.adjust),y=fct_rev(Description),fill=Description)) + geom_bar(stat="identity") + scale_fill_manual(values=colorpanel(10,"chocolate4","chocolate1")) + theme_classic() + theme(legend.position="none") + ggtitle("Differentially Methylated Genes") + theme(axis.title.y=element_blank(), text=element_text(size=20))
dev.off()

pdf("GO.hyper.pdf",height=5,width=12)
ggplot(egoHyper[1:10,],aes(x=-log10(p.adjust),y=fct_rev(Description),fill=Description)) + geom_bar(stat="identity") + scale_fill_manual(values=colorpanel(10,"deeppink4","pink")) + theme_classic() + theme(legend.position="none") + ggtitle("Hypermethylated Genes") + theme(axis.title.y=element_blank(), text=element_text(size=20))
dev.off()

pdf("GO.hypo.pdf",height=5,width=12)
ggplot(egoHypo[1:10,],aes(x=-log10(p.adjust),y=fct_rev(Description),fill=Description)) + geom_bar(stat="identity") + scale_fill_manual(values=colorpanel(10,"cadetblue4","cadetblue1")) + theme_classic() + theme(legend.position="none") + ggtitle("Hypomethylated Genes") + theme(axis.title.y=element_blank(), text=element_text(size=20))
dev.off()
