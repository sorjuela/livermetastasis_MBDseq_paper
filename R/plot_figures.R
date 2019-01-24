#!/usr/bin/env Rscript
# chmod +x 
# run as [R < scriptName.R --no-save]

#########################################################################################
# R script to plot anotations from regions used/detected in csaw run
#
# MBD-seq with 7 metastatis untreated, 5 metastasis treated, 6 CRC, 3 normal mucosa
# Collaboration with Mirco Menigatti and Giancarlo Marra
#
# Stephany Orjuela, December 2018
#########################################################################################

library("ggplot2")
library("tibble")
library("UpSetR")
library("cowplot")

load("MBD_csaw_verify_mod_3grps2_nomapqfilt.RData")

#### Annotation plots ####
#From a total of 289522 DMRs

pdf("myFigs/annotation_plots_mapq.pdf")

#intergenic Vs intragenic
inter <- length(res$gene[grep("^intergenic\\.NA$", res$gene, perl = T)])
intra <- length(res$gene) - inter
d <- data.frame(gene.location = c("intergenic", "intragenic"), number.regions = c(inter, intra))

ggplot(data=d, aes(x=gene.location, y=number.regions)) + 
  geom_bar(stat="identity", width=0.5, fill="steelblue") + 
  geom_text(aes(label = number.regions, x=gene.location, y=number.regions + 2500)) +
  theme_classic() + 
  labs(x = "Location", y = "Number of Regions")

#within intragenic
UTR3 <- grep("3UTR", res$gene)
res_no3UTR <- res[-UTR3,]
proms <- grep("(promoter|1to5kb|5UTR)", res_no3UTR$gene, perl =T)
res_no3UTR_noproms <- res_no3UTR[-proms,]
gene_body <- grep("(exon|intron|CDS)", res_no3UTR_noproms$gene, perl =T)
#res_no3UTR_nogeneb_noproms <- res_no3UTR_noproms[-gene_body,]

d <- tibble(gene.location = factor(c("Regulatory region", "Gene body", "3UTR")),
            number.regions = c(length(proms), length(gene_body), length(UTR3)))
d$gene.location <- factor(d$gene.location, c("Regulatory region", "Gene body", "3UTR"))

ggplot(data=d, aes(x=gene.location, y=number.regions)) + 
  geom_bar(stat="identity", width=0.5, fill="steelblue") + 
  geom_text(aes(label = number.regions, x=gene.location, y=number.regions + 2500)) +
  theme_classic() + 
  labs(x = "Location", y = "Number of Regions")


#open sea Vs close to a CpG Island
inter <- length(res$CpG[grep("^inter:[0-9]+$", res$CpG, perl = T)])
intra <- length(res$CpG) - inter
d <- data.frame(CpG.location = c("Outside SSI", "SSI"), number.regions = c(inter, intra))

ggplot(data=d, aes(x=CpG.location, y=number.regions)) + 
  geom_bar(stat="identity", width=0.5, fill="#F8766D") + 
  geom_text(aes(label = number.regions, x=CpG.location, y=number.regions + 3500)) +
  theme_classic() +
  labs(x = "Location", y = "Number of Regions")

#shelf Vs shore, Vs island
islands <- grep("island", res$CpG)
res_noisl <- res[-islands,]
shore <- grep("shore", res_noisl$CpG)
res_noisl_noshore <- res_noisl[-shore,]
shelf <- grep("shelf", res_noisl_noshore$CpG)
res_noisl_noshore_noshelf <- res_noisl_noshore[-shelf,]
#inter <- length(res_noisl_noshore_noshelf$CpG)

d <- tibble(CpG.location = factor(c("CpG island", "CpG shore", "CpG shelf")),
            number.regions = c(length(islands), length(shore), length(shelf)))
d$CpG.location <- factor(d$CpG.location, c("CpG island", "CpG shore", "CpG shelf"))

ggplot(data=d, aes(x=CpG.location, y=number.regions)) + 
  geom_bar(stat="identity", width=0.5, fill="#F8766D") + 
  geom_text(aes(label = number.regions, x=CpG.location, y=number.regions + 300)) +
  theme_classic() + 
  labs(x = "Location", y = "Number of Regions")

dev.off()

#### Upset plots ####

#choose differential regions for each comparison:

cn.regs <- res[res$cn.de == 1,]
cn.regsGR <- GRanges(cn.regs$seqnames, IRanges(cn.regs$start, cn.regs$end))
mcols(cn.regsGR)$state <- ifelse(cn.regs$cn.direction == "up", "hyper", "hypo")

mn.regs <- res[res$mn.de == 1,]
mn.regsGR <- GRanges(mn.regs$seqnames, IRanges(mn.regs$start, mn.regs$end))
mcols(mn.regsGR)$state <- ifelse(mn.regs$mn.direction == "up", "hyper", "hypo")

#mc.regs <- res[res$mc.de == 1,]
#mc.regsGR <- GRanges(mc.regs$seqnames, IRanges(mc.regs$start, mc.regs$end))


makeUpsetTable <- function(orientationA, orientationS, SSADMRs, AdenDMRs){
  x <- AdenDMRs[mcols(AdenDMRs)$state == orientationA]
  y <- SSADMRs[mcols(SSADMRs)$state == orientationS]
  
  hits <- findOverlaps(x, y)
  comm <- x[queryHits(hits)]
  
  uncomS <- y[-subjectHits(hits)]
  uncomA <- x[-queryHits(hits)]
  
  hyperm <- matrix(0, nrow = sum(length(comm), length(uncomA), length(uncomS)), ncol = 2)
  hyperm[0:length(comm),] <- 1
  hyperm[(length(comm)+1):(length(comm)+length(uncomA)),1] <- 1
  hyperm[(length(comm)+length(uncomA)+1):(length(comm)+length(uncomA)+length(uncomS)),2] <- 1
  
  
  colnames(hyperm)=c("metastasis-normal", "cancer-normal")
  return(hyperm)
  #l <- list(uncomA, uncomS)
  #return(uncomS)
}

hypoM <- makeUpsetTable("hypo","hypo", cn.regsGR, mn.regsGR)
hyperM <- makeUpsetTable("hyper","hyper", cn.regsGR, mn.regsGR)


#Make upset
pdf("myFigs/upset.plots_mapq.pdf")
upset(as.data.frame(hypoM), 
      sets = c("metastasis-normal","cancer-normal"),
      keep.order = T,
      mb.ratio = c(0.8, 0.2), 
      mainbar.y.label = "No. of DMRs (hypomethylated)",
      main.bar.color = "#3b56d8",
      point.size = 3,
      line.size = 1,
      matrix.color = "black",
      mainbar.y.max = 2000, 
      sets.bar.color = "gray23",
      text.scale = c(2,1.5,1.5,1,2,2.5),
      sets.x.label = "Total number of DMRs")

upset(as.data.frame(hyperM), 
      sets = c("metastasis-normal","cancer-normal"),
      keep.order = T,
      mb.ratio = c(0.8, 0.2), 
      mainbar.y.label = "No. of DMRs (hypermethylated)",
      main.bar.color = "#e1cf22",
      point.size = 3,
      line.size = 1,
      matrix.color = "black",
      mainbar.y.max = 2000, 
      sets.bar.color = "gray23",
      text.scale = c(2,1.5,1.5,1,2,2.5),
      sets.x.label = "Total number of DMRs")
dev.off()

#### More annotation ####

cn.regs$state <- ifelse(cn.regs$cn.direction == "up", "hyper", "hypo")
mn.regs$state <- ifelse(mn.regs$mn.direction == "up", "hyper", "hypo")

pdf("myFigs/DMR_annotation_mapq.pdf")

# draw_annot.perc <- function(res){
#   inter <- length(res$gene[grep("^intergenic\\.NA$", res$gene, perl=T)]) 
#   proc.inter <- inter /  length(res$gene) * 100
#   intra <- ((length(res$gene) - inter )/ length(res$gene)) * 100
#   d <- data.frame(location = c("intergenic", "intragenic"), DMRs = c(proc.inter, intra))
#   return(d)
# }

#converge info from hypo and hyper
build_plot <- function(tab1, tab2, orientation, orientation2, drawfunc){
  d1 <- cbind(drawfunc(tab1[tab1$state == orientation,]), state = orientation, comparison = "c-n")
  d2 <- cbind(drawfunc(tab2[tab2$state == orientation,]), state = orientation, comparison = "m-n")
  d3 <- cbind(drawfunc(tab1[tab1$state == orientation2,]), state = orientation2, comparison = "c-n")
  d4 <- cbind(drawfunc(tab2[tab2$state == orientation2,]), state = orientation2, comparison = "m-n")
  
  d <- rbind(d1,d2,d3,d4)
  
  p1 <- ggplot(data=d, aes(x=location, y=number.DMRs, fill = comparison)) + #change number for percentage
    geom_bar(stat="identity", position=position_dodge()) + 
    theme(legend.position="bottom") +
    labs(x = "Location", y = "Number of DMRs") +
    facet_grid(~state)
  return(p1)
}

#genic annotation
draw_annot_genic <- function(res){
  inter <- length(res$gene[grep("^intergenic\\.NA$", res$gene, perl=T)]) 
  #proc.inter <- inter /  length(res$gene) * 100
  intra <- length(res$gene) - inter 
  d <- data.frame(location = c("intergenic", "intragenic"), number.DMRs = c(inter, intra))
  return(d)
}

#intragenic annotation
draw_annot_intragenic <- function(res){
  UTR3 <- grep("3UTR", res$gene)
  res_no3UTR <- res[-UTR3,]
  proms <- grep("(promoter|1to5kb|5UTR)", res_no3UTR$gene, perl =T)
  res_no3UTR_noproms <- res_no3UTR[-proms,]
  gene_body <- grep("(exon|intron|CDS)", res_no3UTR_noproms$gene, perl =T)
  #res_no3UTR_nogeneb_noproms <- res_no3UTR_noproms[-gene_body,]
  d <- tibble(location = factor(c("Regulatory region", "Gene body", "3UTR")),
              number.DMRs = c(length(proms), length(gene_body), length(UTR3)))
  d$location <- factor(d$location, c("Regulatory region", "Gene body", "3UTR"))
  return(d)
}

#CpG related annotation
draw_annot_cpg_general <- function(res){
  inter <- length(res$CpG[grep("^inter:[0-9]+$", res$CpG, perl = T)])
  intra <- length(res$CpG) - inter
  d <- data.frame(location = c("Outside SSI", "SSI"), number.DMRs = c(inter, intra))
  return(d)
}

#CpG island related annotation
draw_annot_cpg <- function(res){
  islands <- grep("island", res$CpG)
  res_noisl <- res[-islands,]
  shore <- grep("shore", res_noisl$CpG)
  res_noisl_noshore <- res_noisl[-shore,]
  shelf <- grep("shelf", res_noisl_noshore$CpG)
  res_noisl_noshore_noshelf <- res_noisl_noshore[-shelf,]
  d <- tibble(location = factor(c("CpG island", "CpG shore", "CpG shelf")),
              number.DMRs = c(length(islands), length(shore), length(shelf)))
  d$location <- factor(d$location, c("CpG island", "CpG shore", "CpG shelf"))
  return(d)
  
}

build_plot(cn.regs, mn.regs, "hyper", "hypo", draw_annot_genic)
build_plot(cn.regs, mn.regs, "hyper", "hypo", draw_annot_intragenic)
build_plot(cn.regs, mn.regs, "hyper", "hypo", draw_annot_cpg_general)
build_plot(cn.regs, mn.regs, "hyper", "hypo", draw_annot_cpg)

dev.off()
