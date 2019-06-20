#Pipeline

Scripts in order of execution:

1. `csaw_DMRanalysis.R` Generate DMRs with csaw from MBD-seq bam files.
2. `dmrseq_DMRanalysis.R` Generate DMRs with dmrseq from pre-generated BSraw object (too large for repo).
3. `BiSeq_DMCanalysis.R` Generate DMCs with BiSeq from pre-generated BSraw objects (too large for repo).
4. `plot_annotation.R` Generate Figures 1-2.
4. `scatterMBDVsTe.R` Generate Figure 3.C.
5. `CGdensity.R` Generate Figure 3.B.