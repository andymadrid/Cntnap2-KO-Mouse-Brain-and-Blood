# GO analysis from DMR results

# note: the sigRegions object was annotated upstairs (aka another script)

# filter for hyper/hypo DMRs
hyper <- sigRegions[which(sigRegions$areaStat>0 & sigRegions$annotation != "Distal Intergenic"),]
hypo <- sigRegions[which(sigRegions$areaStat<0 & sigRegions$annotation != "Distal Intergenic"),]
all <- sigRegions[which(sigRegions$annotation != "Distal Intergenic"),]

# convert symbols to Entrez IDs
e <- bitr(hyper$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db")
genesHyper <- e[,2]
e <- bitr(hypo$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db")
genesHypo <- e[,2]
e <- bitr(all$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db")
genesAll <- e[,2]

# run the GO
egoHyper <- enrichGO(gene=genesHyper,ont="BP",readable=T,OrgDb="org.Mm.eg.db")
egoHypo <- enrichGO(gene=genesHypo,ont="BP",readable=T,OrgDb="org.Mm.eg.db")
egoAll <- enrichGO(gene=genesAll,ont="BP",readable=T,OrgDb="org.Mm.eg.db")

# remove redundant terms
egoAll <- simplify(egoAll, cutoff=0.7, by="p.adjust", select_fun=min)
egoHyper <- simplify(egoHyper, cutoff=0.7, by="p.adjust", select_fun=min)
egoHypo <- simplify(egoHypo, cutoff=0.7, by="p.adjust", select_fun=min)
