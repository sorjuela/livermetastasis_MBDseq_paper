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

library(csaw)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

# Counting Number of Reads in Each Bin
# not include low quality reads and eliminate empty bins
# Not use overlapping windows
# 
# param <- readParam(restrict = paste0("chr",seq(1:22)), minq=50, dedup=TRUE)
# 
# data <- windowCounts(as.character(samples$file), ext = 180,
#                      width = 200, spacing = 200, shift = 0,
#                      param=param, filter=30, bin = T)
# 
# grdata <- rowRanges(data)

#save(grdata, file = "GCcontent_per200bpbins.RData")

### Do same as above for regions not covered by MBD-seq

# grnotcov <- gaps(grdata)
# w <- strand(grnotcov) == "*"
# grnotcov <- grnotcov[w]
# grnotcov.tiles <- tile(grnotcov, width = 200)
# grnotcov.untiles <- unlist(grnotcov.tiles)
# 
# grnotcov.rates <- CGrates(grnotcov.untiles)
# 
# png("GContent_noncovered.png")
# hist(grnotcov.rates, xlim = c(0,1), breaks = 1000)
# dev.off()

#### Probe density ####

probes <- read.table(file="../../../../../Shared_taupo/steph/reference/130912_HG19_CpGiant_4M_EPI.bed") #Available online from Roche
probes <- GRanges(probes$V1, IRanges(start = probes$V2, end = probes$V3)) #240.131
probes <- probes[!seqnames(probes) %in% c("chrX", "chrY", "chrM")]

# png("GContent_probes.png")
# hist(probes.rates, xlim = c(0,2), breaks = 50)
# dev.off()
# 
# #make tiles from probes
# probes.tiles <- tile(probes, width = 200)
# probes.untiles <- unlist(probes.tiles)
# 
# probes.tiles.rates <- CGrates(probes.untiles)
# 
# png("GContent_probestiles.png")
# hist(probes.tiles.rates, xlim = c(0,2), breaks = 50)
# dev.off()
# 
# grnotcov <- gaps(probes)
# w <- strand(grnotcov) == "*"
# grnotcov <- grnotcov[w]
# grnotcov.tiles <- tile(grnotcov, width = 200)
# grnotcov.untiles <- unlist(grnotcov.tiles)
# 
# grnotcov.rates <- CGrates(grnotcov.untiles)
# 
# png("GContent_noncovered_probes.png")
# hist(grnotcov.rates, xlim = c(0,2), breaks = 1000)
# dev.off()


#### Do overlaping densities for each CpG site ####

load("MBD_csaw_verify_mod_3grps2_nomapqfilt3.RData")

#filt test
#filt <- (res$cn.logFC.down == 0 & res$cn.logFC.up == 0)
#res.red <- res[!filt,]
#res <- res.red

resGR <- GRanges(res$seqnames, IRanges(res$start, res$end), 
                 pval = res$cn.PValue,
                 zscore = qnorm(1-(res$cn.PValue/2)),
                 meanlogFC = res$cn.meanlogFC)

#get CpGsite from MBD
load("AllCpGs_GRanges.RData") #with calculated rates
resGRcs <- subsetByOverlaps(cpgr, resGR)
notcovres <- subsetByOverlaps(cpgr, resGRcs, invert = TRUE)

#get CpG rates for probes
overprobes <- subsetByOverlaps(cpgr, probes) #2,768,494
notoverprobes <- subsetByOverlaps(cpgr, probes, invert = TRUE) #23,984,208

overtechs <- subsetByOverlaps(resGRcs) #1,298,194
uniquembd <- subsetByOverlaps(resGRcs, overprobes, invert = TRUE)
uniquete <- subsetByOverlaps(overprobes, resGRcs, invert = TRUE)

#plot

d <- data.frame(cpg_rates = c(resGRcs$rates, notcovres$rates), 
                State=rep(c("covered", "not covered"), 
                        c(length(resGRcs), length(notcovres))))


a1 <- ggplot(d) + 
  geom_density(aes(x=cpg_rates, colour=State, fill=State), alpha=0.5) +
  scale_fill_manual(values = brewer.pal(n = 3, name = "Set1")) +
  scale_color_manual(values = brewer.pal(n = 3, name = "Set1")) +
  theme_classic(base_size = 6) + 
  labs(x = "CpG rates", title = "MBDe") + 
  coord_cartesian(xlim = c(0,2))

d <- data.frame(cpg_rates = c(overprobes$rates, notoverprobes$rates), 
                State=rep(c("covered", "not covered"), 
                        c(length(overprobes), length(notoverprobes))))

a2 <- ggplot(d) + 
  geom_density(aes(x=cpg_rates, colour=State, fill=State), alpha=0.5) +
  scale_fill_manual(values = brewer.pal(n = 3, name = "Set1")) +
  scale_color_manual(values = brewer.pal(n = 3, name = "Set1")) +
  theme_classic(base_size = 6) + 
  labs(x = "CpG rates", title = "Te") + 
  coord_cartesian(xlim = c(0,2))


d <- data.frame(cpg_rates = c(resGRcs$rates, overprobes$rates), 
                State=rep(c("MBDe", "Te"), 
                        c(length(resGRcs), length(overprobes))))

a3 <- ggplot(d) + 
  geom_density(aes(x=cpg_rates, colour=State, fill=State), alpha=0.5) +
  scale_fill_manual(values = brewer.pal(n = 3, name = "Set2")) +
  scale_color_manual(values = brewer.pal(n = 3, name = "Set2")) +
  theme_classic(base_size = 6) + 
  labs(x = "CpG rates", title = "MBDe vs Te") + 
  coord_cartesian(xlim = c(0,2))

d <- data.frame(cpg_rates = c(uniquembd$rates, uniquete$rates, overtechs$rates), 
                State=factor(rep(c("MBDe unique", "Te unique", "shared"), 
                          c(length(uniquembd), length(uniquete), length(overtechs))), 
                          levels = c("MBDe unique", "Te unique", "shared")))


a4 <- ggplot(d) + 
  geom_density(aes(x=cpg_rates, colour=State, fill=State), alpha=0.5) +
  scale_fill_manual(values = brewer.pal(n = 3, name = "Set2")) +
  scale_color_manual(values = brewer.pal(n = 3, name = "Set2")) +
  theme_classic(base_size = 6) + 
  labs(x="CpG rates", title = "MBDe vs Te - intersections") + 
  coord_cartesian(xlim = c(0,2))

a5 <- gridExtra::grid.arrange(a1,a2,a3,a4, ncol = 4)

ggplot2::ggsave("myFigs/GCrates_perTech.pdf", a5, width = 15, height = 3)
