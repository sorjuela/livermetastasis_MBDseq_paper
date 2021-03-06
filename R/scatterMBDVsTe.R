#!/usr/bin/env Rscript

#########################################################################################
# R script to draw delta-meth in both methods (tBS-seq and MBD-seq) relative to CpG density
#
# MBD-seq with 6 CRC, 3 normal mucosa
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
#save(cpgr, file = "AllCpGs_GRanges.RData")
#load("AllCpGs_GRanges.RData") #26,752,702
#total CpGs in low density areas  = 9,881,944
#total CpGs in medium density areas  = 11,572,301
#total CpGs in high density areas  = 5,298,457

#### overlay CpG sites covered by both techs ####

load("betaRegs_Te_nocimpCRCsVsNorm/betaregAll_table.RData") #too heavy
dmcs <- GRanges(betaregAll$chr, 
                IRanges(betaregAll$pos, width =1), 
                meth.diff=betaregAll$meth.diff,
                pval = betaregAll$p.val)

dmcs <- subsetByOverlaps(dmcs, probes) #2,153,119

load("MBD_csaw_allregions.RData")

resGR <- GRanges(res$seqnames, IRanges(res$start, res$end), 
                 pval = res$cn.PValue,
                 meanlogFC = res$cn.meanlogFC)


hist(resGR$zscore, xlim = c(0,4))

#get CpGsite from MBD
over <- findOverlaps(cpgr, resGR)
resGRcs <- cpgr[queryHits(over)]
resGRcs$logFC <- resGR$meanlogFC[subjectHits(over)]
resGRcs$pval <- resGR$pval[subjectHits(over)]

#get overlaps between techs
over <- findOverlaps(resGRcs, dmcs)
shared <- resGRcs[queryHits(over)] #30,164
shared$meth.diff <- dmcs$meth.diff[subjectHits(over)]
shared$te.pval <- dmcs$pval[subjectHits(over)]

#add rate category
shared.df <- as.data.frame(shared)
shared.df$rate.category <- ifelse(shared.df$rates < 0.6, "Medium", "High")
shared.df$rate.category <- ifelse(shared.df$rates < 0.3, "Low", shared.df$rate.category)

#relevel rate
shared.df$rate.category <- factor(shared.df$rate.category, levels = c("Low", "Medium", "High"))

#add coloring vector (ugly way of doing this...)
shared.df$orient <- ifelse(shared.df$meth.diff > 0 & 
                             shared.df$logFC > 0 & 
                             shared.df$pval <= 0.05 &
                             shared.df$te.pval <= 0.05,
                           "hyper", "non")
shared.df$orient <- ifelse(shared.df$meth.diff < 0 & 
                             shared.df$logFC < 0 &
                             shared.df$pval <= 0.05 &
                             shared.df$te.pval <= 0.05, 
                           "hypo", shared.df$orient)

shared.df$orient <- ifelse(shared.df$meth.diff > 0 & 
                             shared.df$pval > 0.05 &
                             shared.df$te.pval <= 0.05, 
                           "hyperTe", shared.df$orient)

shared.df$orient <- ifelse(shared.df$logFC > 0 &
                             shared.df$pval <= 0.05 &
                             shared.df$te.pval > 0.05, 
                           "hyperMBD", shared.df$orient)

shared.df$orient <- ifelse(shared.df$meth.diff < 0 & 
                             shared.df$pval > 0.05 &
                             shared.df$te.pval <= 0.05, 
                           "hypoTe", shared.df$orient)

shared.df$orient <- ifelse(shared.df$logFC < 0 &
                             shared.df$pval <= 0.05 &
                             shared.df$te.pval > 0.05, 
                           "hypoMBD", shared.df$orient)


# Calculate correlation for each facet
library("plyr")
cors <- ddply(shared.df, .(rate.category), summarise, cor = round(cor(meth.diff, logFC), 2))

#### plot ####

myColor <- RColorBrewer::brewer.pal(9, "Set1")[c(1:6,9)]

#barplots
d <- table(shared.df$rate.category, shared.df$orient)
ex <- expand.grid(rownames(d),colnames(d))
shared.counts <- data.frame(orient = ex$Var2, rate = ex$Var1, number.sites = as.vector(d))

p1 <- ggplot(data=shared.counts, aes(x=orient, y=number.sites, fill = orient)) + 
  geom_bar(stat="identity") + 
  geom_text(aes(label = number.sites, x=orient, y=number.sites + 10000), size = 3) +
  scale_fill_manual(values = myColor) +
  labs(y = "Number of sites") +
  facet_grid(~rate) +
  theme_classic() +
  theme(legend.position="bottom")
  
#main scatters
tempd <- table(shared.df$rate.category)
d <- data.frame(size = paste0("n = ",as.vector(tempd)), rate.category = names(tempd))

p2 <- ggplot() + 
  geom_point(data = shared.df, aes(x = meth.diff, y = logFC, color = orient)) +
  geom_smooth(data = shared.df, aes(x = meth.diff, y = logFC), method = lm, se = FALSE, 
              color = "white") +
  geom_text(data = d, aes(label = size, x=-0.4, y=7), size = 3) +
  geom_text(data = cors, aes(label = paste0("r = ", cor), x=-0.4, y=6.5), size = 3) +
  facet_grid(~rate.category) +
  scale_color_manual(values = myColor) +
  geom_hline(yintercept = 0, color = "white") +
  geom_vline(xintercept = 0, color = "white") +
  theme_classic() +
  theme(legend.position="none")

p3 <- cowplot::plot_grid(p2, p1, ncol=1, nrow = 2, align="v", rel_widths = 1, rel_heights = c(2,1))
ggplot2::ggsave("myFigs/MBDscatter_barplot_predFC_lm.png", p3, width = 12, height = 12)

