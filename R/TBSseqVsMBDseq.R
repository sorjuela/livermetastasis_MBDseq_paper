#!/usr/bin/env Rscript

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
source("get_table_with_annots.R")

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
load("MBD_csaw_verify_mod_3grps2_nomapqfilt_final.RData")

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

#load MLH1 proficient samples results (DMRs) from dmrseq_DMRanalysis script
load("noncimpCRCsVsNORM_DMRs_dmrseq.RData") #105,580

#DMRsfilt <- subsetByOverlaps(DMRs,probes)
DMRsfilt <- DMRs[DMRs$qval < 0.05,]  #13,220

DMRsfilt <- get_table_with_annots(DMRsfilt)

DMRstab <- (as(DMRsfilt, "data.frame"))
DMRstab$beta <- DMRstab$beta * -1
write.table(DMRstab[,-c(4:7,9,12:14)], "noncimpCRCsVsNORM_DMRs_Te_dmrseq.csv", 
            row.names=F, quote=FALSE, sep="\t")
DMRsbed <- DMRstab[,1:3]
write.table(DMRsbed, "noncimpCRCsVsNORM_DMRs_Te_dmrseq.bed",col.names= F, 
            row.names=F, quote=FALSE, sep="\t")

#count
# sum(DMRsfilt$median.meth.diff > 0) #1306 , 2406
# sum(DMRsfilt$median.meth.diff < 0) #8926 , 13144
# DMRsfilthyper <- DMRsfilt[DMRsfilt$median.meth.diff > 0]
# DMRsfilthypo <- DMRsfilt[DMRsfilt$median.meth.diff < 0]

#positive beta is hypo
sum(DMRsfilt$beta < 0 & DMRsfilt$qval < 0.05) #1306, dmrseq= 3187
sum(DMRsfilt$beta > 0 & DMRsfilt$qval < 0.05) #8926 , dmrseq = 10,033
DMRsfilthyper <- DMRsfilt[DMRsfilt$beta < 0 & DMRsfilt$qval < 0.05]
DMRsfilthypo <- DMRsfilt[DMRsfilt$beta > 0 & DMRsfilt$qval < 0.05]


#get hyper MBD regions
resGR <- resGR[res$cn.de == 1] #2155
sum(resGR$direction == "up") #2020
sum(resGR$direction == "down") #135
resGRup <- resGR[resGR$direction == "up"]
resGRdown <- resGR[resGR$direction == "down"]

#unique hyper MBD
over <- subsetByOverlaps(resGRup, DMRsfilthyper, invert =TRUE) #1584 , dmrseq = 617
DMRstab <- (as(over, "data.frame"))
write.table(DMRstab, "unique_hyper_DMRs_MBDe.csv", row.names=F, quote=FALSE, sep="\t")

#unique hyper probes
over <- subsetByOverlaps(DMRsfilthyper,resGRup,  invert =T) #827 , dmrse = 1843
DMRstab <- (as(over, "data.frame"))
write.table(DMRstab, "unique_hyper_DMRs_Te.csv", row.names=F, quote=FALSE, sep="\t")


#overlap
over <- findOverlaps(resGRup, DMRsfilthyper) #dmrseq = 1415
length(unique(sort(queryHits(over)))) # 436 , dmrseq = 1403
length(unique(sort(subjectHits(over)))) #479, dmrseq = 1344

#unique hypo MBD 
over <- subsetByOverlaps(resGRdown, DMRsfilthypo, invert =T) #127

#unique hypo probes
over <- subsetByOverlaps(DMRsfilthypo,resGRdown,  invert =T) #8918

#overlap
over <- findOverlaps(resGRdown, DMRsfilthypo) #8, dmrseq = 25
length(unique(sort(queryHits(over)))) # 8, 25
length(unique(sort(subjectHits(over)))) #8, 24

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

