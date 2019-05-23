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
library("ROCR")

#setwd("/run/user/1000/gvfs/sftp:host=imlssherborne.uzh.ch/home/Shared_penticton/data/seq/mirco_mets_mbdseq/R")

#### Overlap probes with MBD regions ####
#Load Roche probes
probes <- read.table(file="../../../../../Shared_taupo/steph/reference/130912_HG19_CpGiant_4M_EPI.bed") #Available online from Roche
probes <- GRanges(probes$V1, IRanges(start = probes$V2, end = probes$V3)) #240.131
probes <- probes[!seqnames(probes) %in% c("chrX", "chrY", "chrM")] #234943

#number of somatic bp covered by probes 79,175,849 (almost 3% of somatic genome)
sum(width(probes))/1000 #79,175.85 kb

#hg19 cpg sites, based on https://support.bioconductor.org/p/95239/)
chrs <- names(Hsapiens)[1:22]
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr <- do.call(c, lapply(1:22, 
                          function(x) GRanges(names(Hsapiens)[x], 
                                              IRanges(cgs[[x]], width = 2)))) #26,752,702

#number of somatic CpG sites covered
overprobes <- subsetByOverlaps(cpgr, probes) #2,768,494 (10% of somatic CpGs in the entire genome)

#number of CpG islands covered by Roche
annotations <- build_annotations(genome = 'hg19', annotations = "hg19_cpg_islands")
over <- subsetByOverlaps(annotations, probes)

#load covered MBDseq regions
load("MBD_csaw_verify_mod_3grps2_nomapqfilt3.RData") #<-- full res 
load("MBD_csaw_verify_mod_3grps2_nomapqfilt_final.RData")
#load("unsorted.res.3Vs3crc.RData") #<-- crc vs norm
#load("unsorted.res.3Vs3crc.from0.RData") #<-- crc vs norm from 0

resGR <- GRanges(res$seqnames, IRanges(res$start, res$end),
                 pval = res$cn.PValue,
                 zscore = qnorm(1-(res$cn.PValue/2)),
                 meanlogFC = res$cn.meanlogFC,
                 FDR = res$cn.FDR,
                 direction = res$cn.direction)


#### get numbers ####

sum(width(resGR)) #249,547,640

#number of somatic CpG sites covered by MBD
overcsaw <- subsetByOverlaps(cpgr, resGR) #7,623,905 for filtered data (29%)

#CpG sites covered by both technologies
over <- findOverlaps(overprobes, overcsaw) #1,144,600 (4%)

#"regions" covered by both, kb
over <- GenomicRanges::intersect(probes, resGR)
sum(width(over)) #26,395,873

#probes not covered by csaw
probesNocov <- subsetByOverlaps(probes,resGR, invert = T) #152516/234943 65%
#266,057/322,551, 82% of MBD not covered by probes

#### Overlap TE DMRs with MBD DMRs ####

#load MLH1 proficient samples results (DMCs)
load("../../../../../sorjuela/serrated_pathway_paper/BiSulf/Data/nonCIMP.allclusters.trimmed.RData")

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
DMRstab <- (as(over, "data.frame"))
write.table(DMRstab, "unique_hyper_DMRs_MBDe.csv", row.names=F, quote=FALSE, sep="\t")

#unique hyper probes
over <- subsetByOverlaps(DMRsfilthyper,resGRup,  invert =T) #827
DMRstab <- (as(over, "data.frame"))
write.table(DMRstab, "unique_hyper_DMRs_Te.csv", row.names=F, quote=FALSE, sep="\t")


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

#### Get recall and ROC curves ####

DMRsaden_crc <- DMRs

# make.recall <- function(resGR, DMRsaden_crc){
#   over <- findOverlaps(resGR, DMRsaden_crc) #7779 TPs
#   tps <- length(over)
#   nonover <- subsetByOverlaps(DMRsaden_crc, resGR, invert = T) #1968 FN?
#   fns <- length(nonover)
#   nonover <- subsetByOverlaps(resGR, DMRsaden_crc, invert = T) #FPs?
#   fps <- length(nonover)
#   
#   mrecal <- data.frame(matrix(NA, nrow = tps+fns+fps, ncol = 4 ))
#   colnames(mrecal) <- c("score", "fdr", "zscore", "labels")
#   
#   #all D=1 7779+1968 = 9747
#   mrecal$labels[1:(tps+fns)] <- 1 #all D=1
#   mrecal$score[1:tps] <- resGR$pval[queryHits(over)] #the ones that overlap with the truth, which should be TPs
#   mrecal$fdr[1:tps] <- resGR$FDR[queryHits(over)]
#   mrecal$zscore[1:tps] <- resGR$zscore[queryHits(over)]
#   
#   # D=0
#   mrecal$labels[is.na(mrecal$labels)] <- 0
#   
#   mrecal$score[mrecal$labels == 0] <- nonover$pval
#   mrecal$fdr[mrecal$labels == 0] <- nonover$FDR
#   mrecal$zscore[mrecal$labels == 0] <- nonover$zscore
#   
#   #add lowest score to NA vals
#   mrecal$score[is.na(mrecal$score)] <- 0
#   mrecal$fdr[is.na(mrecal$fdr)] <- 0
#   mrecal$zcore[is.na(mrecal$zscore)] <- 0
#   
#   png("myFigs/precREC_TEastruth_cor.png")
#   pred <- prediction(mrecal$zscore,mrecal$labels)
#   perf <- performance(pred,"prec","rec")
#   plot(perf, colorize=T)
#   dev.off()
#   
#   png("myFigs/roc_TEastruth_cor.png")
#   pred <- prediction(mrecal$zscore,mrecal$labels)
#   perf <- performance(pred,"tpr","fpr")
#   plot(perf, colorize=T)
#   dev.off()
# }
# 
# make.recall(resGR, DMRsaden_crc)

#### get recall for CpG sites ####
#all sites captured by probes assay are considered the tested units

load("betaRegs_Te_nocimpCRCsVsNorm/betaregAll_table.RData") #3 crc Vs 3 norm
dmcs <- GRanges(betaregAll$chr,  #4,686,717
                IRanges(betaregAll$pos, width =1), 
                meth.diff=betaregAll$meth.diff,
                pval = betaregAll$p.val,
                p.li = p.adjust(betaregAll$p.val,method="holm"))

dmcs <- subsetByOverlaps(dmcs, probes) #2,153,119

mrecal <- data.frame(matrix(0, nrow = length(dmcs), ncol = 2))
colnames(mrecal) <- c("zscore", "labels")
mrecal$labels <- ifelse(dmcs$p.li <= 0.05 & abs(dmcs$meth.diff >= 0.2 ), 1, 0)

over <- findOverlaps(dmcs, resGR)
resGRcs <- dmcs[queryHits(over)] #2,296,797
resGRcs$zscore <- resGR$zscore[subjectHits(over)]
#resGRcs$logfc <- resGR$meanlogFC[subjectHits(over)]

mrecal$zscore[queryHits(over)] <- resGRcs$zscore

png("myFigs/precREC_ROC_ss.png", width = 800, height = 480)
par(mfrow=c(1,2))

myColor <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
pred <- prediction(mrecal$zscore,mrecal$labels)
perf <- performance(pred,"prec","rec")
plot(perf, colorize = T, colorize.palette = myColor, lwd = 5)

perf <- performance(pred,"tpr","fpr")
plot(perf, colorize = T, colorize.palette = myColor, lwd = 5)
dev.off()


#### where are the non-overlapping dmrs? ####
resnonover <- subsetByOverlaps(resGRup, DMRsfilthyper, invert = T)

pdf("nonoverlapping_MBDhyperDMRs_annot_plots.pdf")
annot_plots(resnonover)
dev.off()

