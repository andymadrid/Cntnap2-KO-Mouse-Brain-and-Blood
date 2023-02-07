# Let's make a Sankey plot

# set up a function to do this automatically
makeSankey <- function(x,y) {
  # load some pacakges
  library(networkD3)
  library(dplyr)

  # get the total number of DMRs found
  toUse <- x
  totalNum <- nrow(toUse)
  
  # get the percent of DMRs falling into each feature
  intronPer <- 100*(length(grep("Intron",x$annotation))/totalNum)
  exonPer <- 100*(length(grep("Exon",x$annotation))/totalNum)
  promoterPer <- 100*(length(grep("Promoter",x$annotation))/totalNum)
  utr5Per <- 100*(length(grep("5' UTR",x$annotation))/totalNum)
  utr3Per <- 100*(length(grep("3' UTR",x$annotation))/totalNum)
  downstreamPer <- 100*(length(grep("Downstream",x$annotation))/totalNum)
  interPer <- 100*(length(grep("Intergenic",x$annotation))/totalNum)
  
  # round up percentages if they're greater than 1%
  if(intronPer > 1) { intronPer <- round(intronPer)}
  if(exonPer > 1) { exonPer <- round(exonPer)}
  if(promoterPer > 1) { promoterPer <- round(promoterPer)}
  if(utr5Per > 1) { utr5Per <- round(utr5Per)}
  if(utr3Per > 1) { utr3Per <- round(utr3Per)}
  if(downstreamPer > 1) { downstreamPer <- round(downstreamPer)}
  if(interPer > 1) { interPer <- round(interPer)}
  
  # set percentages for features with <1%
  if(intronPer < 1) { intronPer <- signif(intronPer,digits=1)}
  if(exonPer < 1) { exonPer <- signif(exonPer,digits=1)}
  if(promoterPer < 1) { promoterPer <- signif(promoterPer,digits=1)}
  if(utr5Per < 1) { utr5Per <- signif(utr5Per,digits=1)}
  if(utr3Per < 1) { utr3Per <- signif(utr3Per,digits=1)}
  if(downstreamPer < 1) { downstreamPer <- signif(downstreamPer,digits=1)}
  if(interPer < 1) { interPer <- signif(interPer,digits=1)}
  
  # name the nodes
  name1 <- paste0("Differentially Methylated Regions (DMRs) (N = ",totalNum,")")
  name2 <- paste0("Promoter (",promoterPer,"%)")
  name3 <- paste0("5' UTR (",utr5Per,"%)")
  name4 <- paste0("Exonic (",exonPer,"%)")
  name5 <- paste0("Intronic (",intronPer,"%)")
  name6 <- paste0("3' UTR (",utr3Per,"%)")
  name7 <- paste0("Downstream (",downstreamPer,"%)")
  name8 <- paste0("Intergenic (",interPer,"%)")
  
  # set the links
  links <- data.frame(source=c(name1,name1,name1,name1,name1,name1,name1),
    target=c(name2,name3,name4,name5,name6,name7,name8),
    value=c(promoterPer,utr5Per,exonPer,intronPer,utr3Per,downstreamPer,interPer))
    nodes <- data.frame(name=c(as.character(links$source),as.character(links$target)) %>% unique())

  # set up groups for connections
  links$group <- as.factor(c("type_a","type_b","type_c","type_d","type_e","type_f","type_g"))

  # set up groups for nodes
  nodes$group <- as.factor(c("group1"))

  # set up color scheme
  my_color <- 'd3.scaleOrdinal() .domain(["type_a","type_b","type_c","type_d","type_e","type_f","type_g","group1"]) .range(["#FF0000","#556A5B","#50A45C","#F2AD00","#F69100","#C49647","#5BBCD6","grey"])'
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  
  # plot the Sankey
  sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget",Value = "value", NodeID = "name",colourScale=my_color, LinkGroup="group", NodeGroup="group",fontSize=20,sinksRight=FALSE)
}

# note: dmrs object should be already annotated to regions
makeSankey(sigRegions)
