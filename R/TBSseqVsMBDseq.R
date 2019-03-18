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
library("BiSeq")
library("annotatr")

#setwd("/run/user/1000/gvfs/sftp:host=imlssherborne.uzh.ch/home/Shared_penticton/data/seq/mirco_mets_mbdseq/R")

#### Overlap probes with MBD regions ####
#Load Roche probes
probes <- read.table(file="../../../../../Shared_taupo/steph/reference/130912_HG19_CpGiant_4M_EPI.bed") #Available online from Roche
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
load("MBD_csaw_verify_mod_3grps2_nomapqfilt2.RData")
resGR <- GRanges(res$seqnames, IRanges(res$start, res$end), 
                 pval = res$cn.PValue,
                 zscore = qnorm(1-(res$cn.PValue/2)),
                 meanlogFC = res$cn.meanlogFC,
                 FDR = res$cn.FDR)
#csawGR <-  
#filtered by abundance 322,551

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

#### Overlap TE DMRs with MBD DMRs ####

#load MLH1 proficient samples results (DMCs)
load("../../../../../sorjuela/serrated_pathway_paper/BiSulf/Data/nonCIMP.allclusters.trimmed.RData")

#over <- findOverlaps(epiprobesDMCs,probes)
#allClustersinProbes <- allClusters.trimmed[unique(sort(queryHits(over))),]
#DMRs <- findDMRs(allClustersinProbes, max.dist = 200, diff.dir = T) #11,555
DMRs <- findDMRs(allClusters.trimmed, max.dist = 200, diff.dir = F) #19,478

DMRs <- subsetByOverlaps(DMRs,probes) #10,608
DMRsfilt <- DMRs[width(DMRs) >= 3] #10,232

DMRsfilt <- get_table_with_annots(DMRsfilt)

#export tables
DMRstab <- (as(DMRsfilt, "data.frame"))
write.table(DMRstab, "noncimpCRCsVsNORM_DMRs_Te.csv", row.names=F, quote=FALSE, sep="\t")
DMRsbed <- DMRstab[,1:3]
write.table(DMRsbed, "noncimpCRCsVsNORM_DMRs_Te.bed",col.names= F, row.names=F, quote=FALSE, sep="\t")

#count
sum(DMRsfilt$median.meth.diff > 0) #1,306
sum(DMRsfilt$median.meth.diff < 0) #8,926
DMRsfilthyper <- DMRsfilt[DMRsfilt$median.meth.diff > 0]
DMRsfilthypo <- DMRsfilt[DMRsfilt$median.meth.diff < 0]

#get hyper MBD regions
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
over <- findOverlaps(resGRup, DMRsfilthyper)
length(unique(sort(queryHits(over)))) # 436
length(unique(sort(subjectHits(over)))) #479

#unique hypo MBD 
over <- subsetByOverlaps(resGRdown, DMRsfilthypo, invert =T) #127

#unique hypo probes
over <- subsetByOverlaps(DMRsfilthypo,resGRdown,  invert =T) #8918

#overlap
over <- findOverlaps(resGRdown, DMRsfilthypo) #8
length(unique(sort(queryHits(over)))) # 8
length(unique(sort(subjectHits(over)))) #8

#### load adenoma samples results (DMCs) ####
load("../../../../../sorjuela/serrated_pathway_paper/BiSulf/Data/AdenVsNorm.DMRs.RData")

DMRsaden_crc <- subsetByOverlaps(biseqDMRs, DMRs) #87,529 / 19,478 / over: 9745

#get recall
recall <- function(grreg){
  tp <- length(subsetByOverlaps(grreg, resGR)) #all hits on the truth 1435
  fn <- length(subsetByOverlaps(grreg, resGR, invert = T)) #all non hits on truth 7622
  tpr <- tp/length(grreg) #0.15
  tpr
}
recall(DMRsaden_crc)

#get the meth.diff val for each overlapping DMR, and build the recall curve??

over <- findOverlaps(resGR, DMRsaden_crc) #7779 TPs
#subDMRsaden_crc <- DMRsaden_crc[subjectHits(over)]
#subresGR <- resGR[queryHits(over)]

#nonover <- subsetByOverlaps(DMRsaden_crc, resGR, invert = T) #1968 FN?
nonover <- subsetByOverlaps(resGR, DMRsaden_crc, invert = T) #319,424 FP

#Try to make a table
mrecal <- data.frame(matrix(NA, nrow = 7779+1968+319424, ncol = 4 ))
colnames(mrecal) <- c("score", "fdr", "zscore", "labels")

#all D=1 7779+1968 = 9747
mrecal$labels[1:9747] <- 1 #all D=1
mrecal$score[1:7779] <- resGR$pval[queryHits(over)] #the ones that overlap with the truth, which should be TPs
mrecal$fdr[1:7779] <- resGR$FDR[queryHits(over)]
mrecal$zscore[1:7779] <- resGR$zscore[queryHits(over)]

#the ones that overlap with the truth, which should be TPs

mrecal$labels[is.na(mrecal$labels)] <- 0

mrecal$score[mrecal$labels == 0] <- nonover$pval
mrecal$fdr[mrecal$labels == 0] <- nonover$FDR
mrecal$zscore[mrecal$labels == 0] <- nonover$zscore

mrecal$score[is.na(mrecal$score)] <- 1
mrecal$fdr[is.na(mrecal$fdr)] <- 1
mrecal$zcore[is.na(mrecal$zscore)] <- 1

#Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
#library(benchmarkR)
library(ROCR)

# methods <- list(
#   score = as.matrix(data.frame( 
#     score = mrecal$score
#   )), 
#   padj = as.matrix(data.frame(score = mrecal$fdr)),
#   labels = mrecal$labels)
# 
# re <- SimResults(pval = methods$score, padj= methods$padj, labels = methods$labels)

png("precREC_TEastruth.png")
pred <- prediction(mrecal$zscore,mrecal$labels)
perf <- performance(pred,"prec","rec")
plot(perf, colorize=T)
dev.off()

png("roc_TEastruth.png")
pred <- prediction(mrecal$zscore,mrecal$labels)
perf <- performance(pred,"tpr","fpr")
plot(perf, colorize=T)
dev.off()



#### where are the non-overlapping dmrs? ####
resnonover <- subsetByOverlaps(resGRup, DMRsfilthyper, invert = T)

pdf("nonoverlapping_MBDhyperDMRs_annot_plots.pdf")
annot_plots(resnonover)
dev.off()

