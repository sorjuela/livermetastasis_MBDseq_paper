#!/usr/bin/env Rscript

#########################################################################################
# R script to calculate CG rates based on CGI determination (Gardiner-Garden & Frommer, 1987)
#
# MBD-seq with 7 metastatis untreated, 5 metastasis treated, 6 CRC, 3 normal mucosa
# Collaboration with Mirco Menigatti and Giancarlo Marra
#
# Stephany Orjuela, October 2018
#########################################################################################

library(GenomicRanges)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

#### Probe density ####
#Available online from Roche
probes <- read.table(file="../../../../../Shared_taupo/steph/reference/130912_HG19_CpGiant_4M_EPI.bed")
probes <- GRanges(probes$V1, IRanges(start = probes$V2, end = probes$V3)) #240.131
probes <- probes[!seqnames(probes) %in% c("chrX", "chrY", "chrM")]

#### Do overlaping densities for each CpG site ####

load("MBD_csaw_verify_mod_3grps2_nomapqfilt3.RData")

resGR <- GRanges(res$seqnames, IRanges(res$start, res$end), 
                 pval = res$cn.PValue,
                 zscore = qnorm(1-(res$cn.PValue/2)),
                 meanlogFC = res$cn.meanlogFC)

#get CpGsite from MBD
load("AllCpGs_GRanges.RData") #with calculated rates
resGRcs <- subsetByOverlaps(cpgr, resGR)
notcovres <- subsetByOverlaps(cpgr, resGRcs, invert = TRUE)
ks.test(resGRcs$rates, notcovres$rates)


#get CpG rates for probes
overprobes <- subsetByOverlaps(cpgr, probes) #2,768,494
notoverprobes <- subsetByOverlaps(cpgr, probes, invert = TRUE) #23,984,208

ks.test(resGRcs$rates,overprobes$rates)
# Two-sample Kolmogorov-Smirnov test
# 
# data:  resGRcs$rates and overprobes$rates
# D = 0.18229, p-value < 2.2e-16
# alternative hypothesis: two-sided
wilcox.test(resGRcs$rates,overprobes$rates)
# Wilcoxon rank sum test with continuity correction
# 
# data:  resGRcs$rates and overprobes$rates
# W = 8.5941e+12, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

overtechs <- subsetByOverlaps(resGRcs, overprobes) #1,298,194
uniquembd <- subsetByOverlaps(resGRcs, overprobes, invert = TRUE)
uniquete <- subsetByOverlaps(overprobes, resGRcs, invert = TRUE)
wilcox.test(overtechs$rates,uniquete$rates)
# Wilcoxon rank sum test with continuity correction
# 
# data:  overtechs$rates and uniquete$rates
# W = 9.5287e+11, p-value = 0.02402
# alternative hypothesis: true location shift is not equal to 0
ks.test(overtechs$rates,uniquete$rates)
# Two-sample Kolmogorov-Smirnov test
# 
# data:  overtechs$rates and uniquete$rates
# D = 0.1592, p-value < 2.2e-16
# alternative hypothesis: two-sided

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
