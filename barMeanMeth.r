# Bar plot of mean methylation values

# load pacakges
library(ggplot2)

# read in mean methylation values (by sample)
x <- read.table("meansMethyl.txt",header=T)

# make data frame with mean methylation values per group
# I guess I hard coded this...don't remember doing that...
methMeans <- as.data.frame(c(77,77.23,81.4,81.13))
methMeans$SEM <- c(0.1,0.0333,0.0577,0.35)
colnames(methMeans) <- c("Mean","SEM")

# switch up the order of the tissues for plotting
x$Group <- factor(x$Group,levels=unique(x$Group))
x <- x[c(7:12,1:6),]
x$Group <- rep(c("WT HC","HOMO HC","WT Bld","HOMO Bld"),c(3,3,3,3))
methMeans$Group <- c("WT HC","HOMO HC","WT Bld","HOMO Bld")
methMeans$Color <- c("goldenrod3","dodgerblue2","forestgreen","firebrick2")
methMeans$Group <- factor(methMeans$Group,levels=methMeans$Group)

# plot it out
ggplot(methMeans, aes(x=Group, y=Mean,fill=Group)) + geom_bar(stat="identity",width=0.5,fill=methMeans$Color) + geom_point(data=x,aes(x=Group,y=Mean),size=3,colour="grey") + theme(text=element_text(size=20)) + ylab(label="Mean Methylation (%)") + xlab(label="") + theme(legend.position="none") + geom_errorbar( aes(x=Group, ymin=Mean-SEM, ymax=Mean+SEM),width=0.4,colour="black",alpha=0.9,size=1.3) + theme_classic() + theme(text=element_text(size=20)) + coord_cartesian(ylim=c(75,82)) + theme(axis.text.x=element_text(angle=45,vjust = 0.5)) + theme(legend.position="none")
