#!/usr/bin/env Rscript

#########################################################################################
# R script to get DMRs using dmrseq (tBS-seq) 
#
# MBD-seq with 6 CRC, 3 normal mucosa
# Collaboration with Mirco Menigatti and Giancarlo Marra
#
# Stephany Orjuela, June 2019
#########################################################################################

library(dmrseq)

load("../../../../../sorjuela/serrated_pathway_paper/BiSulf/Data/CRCdata/BSrawCRC.RData")
bs <- BSseq(chr = seqnames(BSr), 
            pos = start(BSr),
            M = methReads(BSr)[,c(3,4,9,10,11,12)], 
            Cov = totalReads(BSr)[,c(3,4,9,10,11,12)], 
            sampleNames = colnames(BSr)[c(3,4,9,10,11,12)])
pData(bs)$group <- colData(BSr)$group[c(3,4,9,10,11,12)]
pData(bs)$group <- factor(pData(bs)$group, levels = c("non","NORMAL.non"))

pData(bs)$patient <- rep(1:3, each = 2)

bs <- sort(bs)
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov")==0) == 0)
bs <- bs[loci.idx,]
save(bs, file = "BSobj.RData")

#for some reason the old bsseq obj is more heavy than it should be
#so I re-run the BSseq object from above

DMRs <- dmrseq(bs=bs_new, testCovariate="group", cutoff = 0.05, 
               BPPARAM = MulticoreParam(1),
               adjustCovariate = "patient")

save(DMRs, file = "noncimpCRCsVsNORM_DMRs_dmrseq.RData")