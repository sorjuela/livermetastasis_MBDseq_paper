#!/usr/bin/env Rscript
# chmod +x 
# run as [R < scriptName.R --no-save]

#########################################################################################
# R script to calculate CG rates based on CGI determination (Gardiner-Garden & Frommer, 1987)
#
# MBD-seq with 7 metastatis untreated, 5 metastasis treated, 6 CRC, 3 normal mucosa
# Collaboration with Mirco Menigatti and Giancarlo Marra
#
# Stephany Orjuela, October 2018
#########################################################################################

library("csaw")
library("BSgenome.Hsapiens.UCSC.hg19")  

#Get metadata
samples <- read.table("../metadata.txt", header =T)
samples <- samples[order(samples$sample),] 
samples$file <- paste0("dedups/", list.files("dedups/", "_s.bam$"))

# Counting Number of Reads in Each Bin
# not include low quality reads and eliminate empty bins
# Not use overlapping windows

param <- readParam(restrict = paste0("chr",seq(1:22)), minq=50, dedup=TRUE)

data <- windowCounts(as.character(samples$file), ext = 180,
                     width = 200, spacing = 200, shift = 0,
                     param=param, filter=30, bin = T)

grdata <- rowRanges(data)

#Get all CpG sites from hg19. Idea from https://support.bioconductor.org/p/95239/

CGrates <- function(grdata){
  chrs <- names(Hsapiens)[1:22]
  cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
  cpgr <- do.call(c, lapply(1:22, 
                          function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 2))))

  #Count number of CpG sites in window
  over <- findOverlaps(grdata, cpgr)
  sHits <- queryHits(over)
  grdata$numCpGsites <- 0
  grdata$numCpGsites[rle(sHits)$values] <- rle(sHits)$lengths

  #Get all Cs fom hg19
  cs <- lapply(chrs, function(x) start(matchPattern("C", Hsapiens[[x]])))
  cgr <- do.call(c, lapply(1:22, 
                          function(x) GRanges(names(Hsapiens)[x], IRanges(cs[[x]], width = 2))))

  over <- findOverlaps(grdata, cgr)
  sHits <- queryHits(over)
  grdata$numCsites <- 0
  grdata$numCsites[rle(sHits)$values] <- rle(sHits)$lengths

  #Get all Gs fom hg19
  gs <- lapply(chrs, function(x) start(matchPattern("G", Hsapiens[[x]])))
  ggr <- do.call(c, lapply(1:22, 
                         function(x) GRanges(names(Hsapiens)[x], IRanges(gs[[x]], width = 2))))

  over <- findOverlaps(grdata, ggr)
  sHits <- queryHits(over)
  grdata$numGsites <- 0
  grdata$numGsites[rle(sHits)$values] <- rle(sHits)$lengths

  GCcont <- (grdata$numCpGsites / 200) / ((grdata$numCsites / 200) * (grdata$numGsites / 200))
  return(GCcont)
}

#save(grdata, file = "GCcontent_per200bpbins.RData")

### Do same as above for regions not covered by MBD-seq

grnotcov <- gaps(grdata)
w <- strand(grnotcov) == "*"
grnotcov <- grnotcov[w]
grnotcov.tiles <- tile(grnotcov, width = 200)
grnotcov.untiles <- unlist(grnotcov.tiles)

grnotcov.rates <- CGrates(grnotcov.untiles)

png("GContent_noncovered.png")
hist(grnotcov.rates, xlim = c(0,1), breaks = 1000)
dev.off()

#### filter main gr by abundance ####

abundances <- aveLogCPM(asDGEList(data))
keep.simple <- abundances > -1 
grdata_filt <- grdata[keep.simple]
#save(grdata_filt, file = "GCcontent_per200bpbins_filt.RData")

mmbdfilt.rates <- CGrates(grdata_filt)

grnotcov <- gaps(grdata_filt)
w <- strand(grnotcov) == "*"
grnotcov <- grnotcov[w]
grnotcov.tiles <- tile(grnotcov, width = 200)
grnotcov.untiles <- unlist(grnotcov.tiles)

grnotcov.rates <- CGrates(grnotcov.untiles)

png("GContent_noncovered_filt.png")
hist(grnotcov.rates, xlim = c(0,2), breaks = 1000)
dev.off()

#### Probe density ####

probes <- read.table(file="../../probes/130912_HG19_CpGiant_4M_EPI.bed") #Available online from Roche
probes <- GRanges(probes$V1, IRanges(start = probes$V2, end = probes$V3)) 
probes <- probes[!seqnames(probes) %in% c("chrX", "chrY", "chrM")]

png("GContent_probes.png")
hist(probes.rates, xlim = c(0,2), breaks = 50)
dev.off()

#make tiles from probes
probes.tiles <- tile(probes, width = 200)
probes.untiles <- unlist(probes.tiles)

probes.tiles.rates <- CGrates(probes.untiles)

png("GContent_probestiles.png")
hist(probes.tiles.rates, xlim = c(0,2), breaks = 50)
dev.off()

grnotcov <- gaps(probes)
w <- strand(grnotcov) == "*"
grnotcov <- grnotcov[w]
grnotcov.tiles <- tile(grnotcov, width = 200)
grnotcov.untiles <- unlist(grnotcov.tiles)

grnotcov.rates <- CGrates(grnotcov.untiles)

png("GContent_noncovered_probes.png")
hist(grnotcov.rates, xlim = c(0,2), breaks = 1000)
dev.off()

#### Draw overlapping densities ####
d <- data.frame(cpg_rates = c(mmbdfilt.rates, grnotcov.rates), 
                cov=rep(c("covered", "not-covered"), 
                         c(length(mmbdfilt.rates), length(grnotcov.rates))))

pdf("MBD_GCcontent_density.pdf")
ggplot(d) + 
  geom_density(aes(x=cpg_rates, colour=cov, fill=cov), alpha=0.5) +
  theme_classic() + 
  xlab("CpG rates") + 
  coord_cartesian(xlim = c(0,2))
dev.off()

d <- data.frame(cpg_rates = c(probes.tiles.rates, grnotcovprobes.rates), 
                cov=rep(c("covered", "not-covered"), 
                        c(length(probes.tiles.rates), length(grnotcovprobes.rates))))

pdf("probes_GCcontent_density.pdf")
ggplot(d) + 
  geom_density(aes(x=cpg_rates, colour=cov, fill=cov), alpha=0.5) +
  theme_classic() + 
  xlab("CpG rates") + 
  coord_cartesian(xlim = c(0,2))

dev.off()


d <- data.frame(cpg_rates = c(mmbdfilt.rates, probes.tiles.rates), 
                cov=rep(c("MBD", "TE"), 
                        c(length(mmbdfilt.rates), length(probes.tiles.rates))))

pdf("MBDandTE_GCcontent_density.pdf")
ggplot(d) + 
  geom_density(aes(x=cpg_rates, colour=cov, fill=cov), alpha=0.5) +
  theme_classic() + 
  xlab("CpG rates") + 
  coord_cartesian(xlim = c(0,2))
dev.off()

d <- data.frame(cpg_rates = c(grnotcov.rates, grnotcovprobes.rates), 
                cov=rep(c("MBD", "TE"), 
                        c(length(grnotcov.rates), length(grnotcovprobes.rates))))

pdf("nonMBDandnonTE_GCcontent_density.pdf")
ggplot(d) + 
  geom_density(aes(x=cpg_rates, colour=cov, fill=cov), alpha=0.5) +
  theme_classic() + 
  xlab("CpG rates") + 
  coord_cartesian(xlim = c(0,2))
dev.off()

#### CpG density from Repitools ####

#this is probably wrong...
window <- width(csawGR)

CGfinds <- sequenceCalc(csawGR, Hsapiens, DNAString("CG"), positions = TRUE)

CGfinds <- Map(function(winList, CGList){
  if(!is.null(CGList)) abs(CGList-winList/2) else CGList
}, window, CGfinds)

cpgDensity <- mapply(function(winList, CGList){ 
  sum(1-(CGList/(winList/2)))
}, window, CGfinds)

pdf("MBDregions_cpgdensity_abundancefilt.pdf")
hist(cpgDensity, breaks = 200, xlim = c(0,200))
dev.off()

#i think this plot only gives you a
# summary of the distribution of CpG density among reads/fragments
#it changes if you change the seq.len
csawlist <- GRangesList(csawGR)
cpgDensityPlot(csawlist, organism = Hsapiens, w.function = "none", seq.len = 180)
