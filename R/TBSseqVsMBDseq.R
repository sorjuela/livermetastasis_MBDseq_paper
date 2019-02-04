#!/usr/bin/env Rscript
# chmod +x 
# run as [R < scriptName.R --no-save]

#########################################################################################
# R script to determine "overlap" between tBS-seq and MBD-seq
#
# MBD-seq with 7 metastatis untreated, 5 metastasis treated, 6 CRC, 3 normal mucosa
# Collaboration with Mirco Menigatti and Giancarlo Marra
#
# Stephany Orjuela, January 2019
#########################################################################################

library("GenomicRanges")
library("BSgenome.Hsapiens.UCSC.hg19")
library("csaw")

#### Overlap probes with MBD regions ####
#Load Roche probes
probes <- read.table(file="../../probes/130912_HG19_CpGiant_4M_EPI.bed") #Available online from Roche
probes <- GRanges(probes$V1, IRanges(start = probes$V2, end = probes$V3)) #240.131
probes <- probes[!seqnames(probes) %in% c("chrX", "chrY", "chrM")] #234943

#number of somatic bp covered by probes 79,175,849 
sum(width(probes))/1000 #79,175.85 kb

#hg19 cpg sites, based on https://support.bioconductor.org/p/95239/)
chrs <- names(Hsapiens)[1:24]
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr <- do.call(c, lapply(1:24, 
                          function(x) GRanges(names(Hsapiens)[x], 
                                              IRanges(cgs[[x]], width = 2)))) #28,217,009

#number of somatic CpG sites covered
overprobes <- subsetByOverlaps(cpgr, probes) #2,768,494

#number of CpG islands covered by Roche
annotations <- build_annotations(genome = 'hg19', annotations = "hg19_cpg_islands")
over <- subsetByOverlaps(annotations, probes)

#load covered MBDseq regions
load("csaw_data_nomapqfilt.RData")
merged <- mergeWindows(rowRanges(data), tol=50L)
csawGR <- merged$region #overall covered MBD regions, including low mapq 810,374
sum(width(csawGR)) #923,681,320, not really informative though

#look at filtered regions
abundances <- aveLogCPM(asDGEList(data))
keep.simple <- abundances > -1 
data2 <- data[keep.simple,]
merged <- mergeWindows(rowRanges(data2), tol=50L)
csawGR <- merged$region #filtered by abundance 322,551
sum(width(csawGR)) #249,547,640

#number of somatic CpG sites covered
overcsaw <- subsetByOverlaps(cpgr, csawGR) #16,693,457, 7,623,905 for filtered data

#CpG sites covered by both technologies
over <- findOverlaps(overprobes, overcsaw) #1,828,983

#"regions" covered by both, kb
over <- GenomicRanges::intersect(probes, csawGR)
sum(width(over)) #26,395,873

#probes not covered by csaw
probesNocov <- subsetByOverlaps(probes,csawGR, invert = T) #152516/234943 65%
#266,057/322,551, 82% of MBD not covered by probes

#### DMR stuff ####

#load MLH1 proficient samples results (DMCs)
load("nonCIMP.allclusters.trimmed.RData")
library(BiSeq)
epiprobesDMCs <- GRanges(allClusters.trimmed$chr, IRanges(allClusters.trimmed$pos, width = 1)) #227,824
#actually hit by probes
epiprobesDMCsinProbes <- subsetByOverlaps(epiprobesDMCs,probes) #124,213
#hypermethylated
sum(epiprobesDMCsinProbes$meth.diff > 0) #30,337

over <- findOverlaps(epiprobesDMCs,probes)
#allClustersinProbes <- allClusters.trimmed[unique(sort(queryHits(over))),]
#DMRs <- findDMRs(allClustersinProbes, max.dist = 200, diff.dir = T) #11,555
DMRs <- findDMRs(allClusters.trimmed, max.dist = 200, diff.dir = F) #19,478

DMRs <- subsetByOverlaps(DMRs,probes)
DMRsfilt <- DMRs[width(DMRs) >= 3] #10,600
sum(DMRsfilt$median.meth.diff > 0) #1,495
sum(DMRsfilt$median.meth.diff < 0) #9,105
DMRsfilthyper <- DMRsfilt[DMRsfilt$median.meth.diff > 0]
DMRsfilthypo <- DMRsfilt[DMRsfilt$median.meth.diff < 0]


#load MBD regions
load("MBD_csaw_verify_mod_3grps2_nomapqfilt.RData")
resGR <- GRanges(res$seqnames, IRanges(res$start, res$end), direction = res$cn.direction)
resGR <- resGR[res$cn.de == 1] #2155
sum(resGR$direction == "up") #2020
sum(resGR$direction == "down") #135
resGRup <- resGR[resGR$direction == "up"]
resGRdown <- resGR[resGR$direction == "down"]

#unique hyper MBD
over <- subsetByOverlaps(resGRup, DMRsfilthyper, invert =T) #1584

#unique hyper probes
over <- subsetByOverlaps(DMRsfilthyper,resGRup,  invert =T) #827

#overlap
over <- subsetByOverlaps(resGRup, DMRsfilthyper) #436

#unique hypo MBD 
over <- subsetByOverlaps(resGRdown, DMRsfilthypo, invert =T) #127

#unique hypo probes
over <- subsetByOverlaps(DMRsfilthypo,resGRdown,  invert =T) #8918

#overlap
over <- subsetByOverlaps(resGRdown, DMRsfilthypo) #8

#### Repitools test ####

library(GenomicRanges)
library(Repitools)
data(samplesList)

library(BSgenome.Hsapiens.UCSC.hg18)
cpgDensityPlot(samples.list.subset, 
               organism = Hsapiens, 
               w.function = "none", 
               seq.len = 300,
               cols = c("black", "green", "orange", "red"), 
               xlim = c(0, 30), lwd = 2)