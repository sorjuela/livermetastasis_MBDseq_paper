#!/usr/bin/env Rscript
# chmod +x 
# run as [R < scriptName.R --no-save]

#########################################################################################
# R script to draw delta-meth in both methods (tBS-seq and MBD-seq) relative to CpG density
#
# MBD-seq with 7 metastatis untreated, 5 metastasis treated, 6 CRC, 3 normal mucosa
# Collaboration with Mirco Menigatti and Giancarlo Marra
#
# Stephany Orjuela, January 2019
#########################################################################################

library("BSgenome.Hsapiens.UCSC.hg19")
library("ggplot2")
library("BiSeq")

#### Get hg19 CpG sites and calculate CpG rates ####
chrs <- names(Hsapiens)[1:22]

#Get all CpGs fom hg19
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr <- do.call(c, lapply(1:22, 
                          function(x) GRanges(names(Hsapiens)[x], 
                                              IRanges(cgs[[x]], width = 2))))

slopped.cpgr <- GenomicRanges::promoters(cpgr, upstream=100, downstream=100)


#Get all Cs fom hg19
cs <- lapply(chrs, function(x) start(matchPattern("C", Hsapiens[[x]])))
cr <- do.call(c, lapply(1:22, 
                         function(x) GRanges(names(Hsapiens)[x], IRanges(cs[[x]], width = 2))))


#Get all Gs fom hg19
gs <- lapply(chrs, function(x) start(matchPattern("G", Hsapiens[[x]])))
gr <- do.call(c, lapply(1:22, 
                         function(x) GRanges(names(Hsapiens)[x], IRanges(gs[[x]], width = 2))))

#calculate O/E cpg rate per window, then asign that value to each cpg (to further classify in plot)

calc_rates <- function(granges, cpgr, cr, gr){
  numCpGs <- countOverlaps(granges, cpgr)
  numCs <- countOverlaps(granges, cr)
  numGs <- countOverlaps(granges, gr)
  rates <- (numCpGs / 200) / ((numCs / 200) * (numGs / 200))
  rates
}


#add CpG rate
cpgr$rates <- calc_rates(slopped.cpgr, cpgr, cr, gr)
save(cpgr, file = "AllCpGs_GRanges.RData")

#### overlay CpG sites covered by both techs ####
#load("../../../../../sorjuela/serrated_pathway_paper/BiSulf/Data/nonCIMP.allclusters.trimmed.RData")
#hist(qnorm(1-(allClusters.trimmed$p.val/2)))

load("betaRegs_Te_nocimpCRCsVsNorm/betaregAll_table.RData")
dmcs <- GRanges(betaregAll$chr, 
                IRanges(betaregAll$pos, width =1), 
                meth.diff=betaregAll$meth.diff,
                pval = betaregAll$p.val,
                #p.li = allClusters.trimmed$p.li,
                #zscore= allClusters.trimmed$z.score,
                myzscore=qnorm(1-(betaregAll$p.val/2))) #4,686,717
#dmcshyper <- dmcs[dmcs$meth.diff > 0]
#hist(dmcshyper$myzscore)
#hist(dmcshyper$zscore)

load("MBD_csaw_verify_mod_3grps2_nomapqfilt2.RData")
resGR <- GRanges(res$seqnames, IRanges(res$start, res$end), 
                 pval = res$cn.PValue,
                 zscore = qnorm(1-(res$cn.PValue/2)),
                 meanlogFC = res$cn.meanlogFC)


hist(resGR$zscore, xlim = c(0,4))

#resGR <- resGR[res$cn.de == 1] #2155
#resGRup <- resGR[resGR$direction == "up"]
#resGRup <- resGR[resGR$meanlogFC > 0]
#hist(resGRup$zscore)

#get CpGsite from MBD
over <- findOverlaps(cpgr, resGR)
resGRcs <- cpgr[queryHits(over)]
resGRcs$zscore <- resGR$zscore[subjectHits(over)]
#resGRcs$logFC <- resGRup$logFCup[subjectHits(over)]
resGRcs$logFC <- resGR$meanlogFC[subjectHits(over)] #<-- this is the meanlogFC from script csaw_analysis

#get overlaps between techs
over <- findOverlaps(resGRcs, dmcs)
shared <- resGRcs[queryHits(over)] #30,164
shared$te.zscore <- dmcs$myzscore[subjectHits(over)]
shared$meth.diff <- dmcs$meth.diff[subjectHits(over)]

#add rate category
shared.df <- as.data.frame(shared)
shared.df$rate.category <- ifelse(shared.df$rates < 0.6, "Medium", "High")
shared.df$rate.category <- ifelse(shared.df$rates < 0.3, "Low", shared.df$rate.category)
save(shared.df, file="AllCs.sharedtechs.withrates.RData")

#shared.df$meth.diff <- shared.df$meth.diff / 100
setwd("/run/user/1000/gvfs/sftp:host=imlssherborne.uzh.ch/home/Shared_penticton/data/seq/mirco_mets_mbdseq/R")
load("AllCs.sharedtechs.withrates.RData") #2,297,602
head(shared.df)

#filter out infinite vals for zscores
#shared.df <- shared.df[!is.infinite(shared.df$te.zscore),]
#shared.df <- shared.df[!is.infinite(shared.df$zscore),]

#refactor rate levels
shared.df$rate.category <- factor(shared.df$rate.category, levels = c("Low", "Medium", "High"))
shared.df$orient <- ifelse(shared.df$meth.diff >= 0.2 & shared.df$logFC >= 1, "hyper", "non")
shared.df$orient <- ifelse(shared.df$meth.diff < -0.2 & shared.df$logFC < -1, "hypo", shared.df$orient)

#zscore scatter...fix to make negative pvals
#ggplot(shared.df, aes(x = te.zscore, y = zscore)) + 
  #geom_point(color = "#313695", alpha = 1/15) +
  #geom_smooth(method=lm) +
  #geom_density_2d() +
  #geom_bin2d(bins = 60) +
  #tat_bin_2d(geom = "polygon") +
  #facet_grid(~rate.category) +
  #geom_abline(slope=slp, intercept=int) +
  #theme_classic() 

myColor <- rev(RColorBrewer::brewer.pal(11, "Spectral"))

#png("MBDlogFCVsTEmethdiff_full_ratefacet.png")
#meth diff scatter
ggplot(shared.df, aes(x = meth.diff, y = logFC)) + #, color = orient)) + 
  geom_point(alpha = 1/10) +
  #geom_bin2d(bins = 60) +
  #geom_hex(bins = 50, color = "black") +
  facet_grid(~rate.category) +
  #scale_color_manual(values = c("#ABDDA4", "#D53E4F", "gray")) +
  #scale_fill_gradientn(colours = myColor) +
  geom_hline(yintercept = 0, color = "white") +
  geom_vline(xintercept = 0, color = "white") +
  theme_classic() 
#dev.off()

#

#### rerun BiSeq for CRC samples ####
#up tp getting b vals since we dont need the zscores?
for(i in 2:22){ 
  
  load(file=paste0("../../../../../sorjuela/serrated_pathway_paper/BiSulf/Data/CRCdata/", 
                   paste0("chr",i,"BSraw.CRC.RData")))
  #w <- which(seqnames(BSr) == paste0("chr",i))
  #BSrchr <- BSr[w, ]
  #rm(BSr)
  #save(BSrchr, file = paste0("chr",i,"BSraw.RData"))
  
  #add patient effect
  colData(BSrchr)$patient <- rep(1:6, each = 2)
  
  ### Define clusters ###----------------------------------------------------------------------------- 
  
  BSr.clust.unlim <- clusterSites(BSrchr,
                                  groups = factor(colData(BSrchr)$group),
                                  perc.samples = 0.8,
                                  min.sites = 3,
                                  max.dist =  500,
                                  mc.cores = 10)
  
  ### smooth the methylation values of CpG sites within the clusters ###----------------------------- 
  predictedMeth <- predictMeth(object = BSr.clust.unlim, h=1000, mc.cores = 10)
  #save(predictedMeth, file=paste0("predictedMethchr",i,".RData"))
  
  ### Model methylation within a beta regression for each comparison ###-----------------------------
  
  #subset for groups, change depending on group of interest
  s <- grep("^non$", colData(predictedMeth)$group)
  n <- grep("NORMAL.non", colData(predictedMeth)$group)
  
  SNpMeth <- predictedMeth[,c(s,n)]
  colData(SNpMeth)$group <- factor(colData(SNpMeth)$group, levels = c("non", "NORMAL.non"))
  
  betaResultsSN <- betaRegression(formula = ~group+patient, link = "probit", object = SNpMeth, type = "BR", mc.cores = 10) 
  save(betaResultsSN, file=paste0("betaRegs_Te_nocimpCRCsVsNorm/betaRes.",i,".RData"))
}

load("betaRegs_Te_nocimpCRCsVsNorm/betaRes.1.RData")
betaregAll <- betaResultsSN
for(i in 2:22){
  load(paste0("betaRegs_Te_nocimpCRCsVsNorm/betaRes.",i,".RData"))
  betaregAll <- rbind(betaregAll, betaResultsSN)
}
betaregAll <- betaregAll[!is.na(betaregAll$meth.diff),]
save(betaregAll, file = "betaRegs_Te_nocimpCRCsVsNorm/betaregAll_table.RData")

#after getting this betaregAll table, run the script from the start using this table instead
#of the hyper subset
    
#### Make this plot for regions (cause this will never work?) ####
DMRs <- findDMRs(allClusters.trimmed, max.dist = 200, diff.dir = F) #19,478

DMRs <- subsetByOverlaps(DMRs,probes) #10,608
DMRsfilt <- DMRs[width(DMRs) >= 3]
DMRsfilthyper <- DMRsfilt[DMRsfilt$median.meth.diff > 0]

over <- findOverlaps(resGRup, DMRsfilthyper)

DMRsfilthyper$rates <- calc_rates(DMRsfilthyper, cpgr, cr, gr)


sharedreg <- data.frame(MBDzcore = resGRup$zscore[queryHits(over)],
                        TEzscore = DMRsfilthyper$zscore.median[subjectHits(over)],
                        MBDrates = resGRup$rates[queryHits(over)],
                        TErates = DMRsfilthyper$rates[subjectHits(over)])

# sharedreg$rate.categoryMBD <- ifelse(sharedreg$MBDrates < 0.6, "Medium", "High")
# sharedreg$rate.categoryMBD <- ifelse(sharedreg$MBDrates < 0.3, "low", sharedreg$rate.categoryMBD)
# 
# sharedreg$rate.categoryTE <- ifelse(sharedreg$TErates < 0.6, "Medium", "High")
# sharedreg$rate.categoryTE <- ifelse(sharedreg$TErates < 0.3, "low", sharedreg$rate.categoryTE)

save(sharedreg, file = "hypermethDMRs.sharedtechs.onlyzscores.RData")

load("hypermethDMRs.sharedtechs.onlyzscores.RData")
head(sharedreg)
sharedreg <- sharedreg[!is.infinite(sharedreg$TEzscore),]

ggplot(sharedreg, aes(x = MBDzcore, y = TEzscore, color = MBDrates)) + 
  geom_point(alpha = 1/15) +
  #myColor_scale_fill +
  #geom_smooth(method=lm) +
  #geom_density_2d() +
  #geom_bin2d() +
  #tat_bin_2d(geom = "polygon") +
  facet_grid(~rate.categoryMBD) +
  #geom_abline(slope=slp, intercept=int) +
  #coord_cartesian(ylim = c(0,10))
  theme_classic() 

