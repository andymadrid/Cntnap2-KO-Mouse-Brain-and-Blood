# Make a basic Circos plot

# load in packages
library(BioCircos)
library(RColorBrewer)
library(wesanderson)

# load in results
load("dmrseq.hc.rdata")
sigRegions <- as.data.frame(regions[which(regions$pval<0.05),])

# set up the colour palette
pal <- as.vector(wes_palette("GrandBudapest1",21,type="continuous"))

# set up chromosome sizes (mm10)
genome <- list('1'=195471971,'2'=182113224,'3'=160039680,'4'=156508116,'5'=151834684,'6'=149736546,'7'=145441459,'8'=129401213,'9'=124595110,'10'=130694993,'11'=122082543,'12'=120129022,'13'=120421639,'14'=124902244,'15'=104043685,'16'=98207768,'17'=94987271,'18'=90702639,'19'=61431566,'X'=171031299,'Y'=91744698)

# get which DMRs are hyper/hypo
hyper <- sigRegions[which(sigRegions$areaStat>0),]
hypo <- sigRegions[which(sigRegions$areaStat<0),]

# keep only relevant columns for plot
hyper <- hyper[,c(1,2,3)]
hypo <- hypo[,1:3]

# remove superfluous "chr"
hyper$chr<-gsub("chr","",hyper$chr)
hypo$chr<-gsub("chr","",hypo$chr)

# set up chromsomes, starts, and ends (ends are artifically large to actually show up in the plot - can edit as needed to look pretty)
hyper_chr <- as.character(hyper[,1])
hyper_start <- hyper[,2]
hyper_end <- hyper[,3]+1000000
hypo_chr <- as.character(hypo[,1])
hypo_start <- hypo[,2]
hypo_end <- hypo[,3]+1000000

# set up tracks for plotting along chromosomes
tracklist <- BioCircosArcTrack('myArcTracks',hypo_chr,hypo_start,hypo_end,opacities=1,colors='blue',minRadius=0.55,maxRadius=0.7)
tracklist <- tracklist + BioCircosArcTrack('myArcTracks',hyper_chr,hyper_start,hyper_end,opacities=1,colors='red',minRadius=0.75,maxRadius=0.9)

# plot the plot
BioCircos(genome=genome,tracklist,genomeFillColor=pal,  genomeTicksLen = 2, genomeTicksTextSize = 0, genomeTicksScale = 1e+8)

