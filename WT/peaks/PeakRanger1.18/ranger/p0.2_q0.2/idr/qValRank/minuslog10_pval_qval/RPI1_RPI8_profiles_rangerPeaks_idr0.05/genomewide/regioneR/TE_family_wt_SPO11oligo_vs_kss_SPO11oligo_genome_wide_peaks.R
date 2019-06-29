#!/applications/R/R-3.3.2/bin/Rscript

# Plot bar chart of log2(observed:expected) peaks overlapping TEs within each family

# Usage:
# Rscript TE_family_wt_SPO11oligo_vs_kss_SPO11oligo_genome_wide_peaks.R "SPO11-1-oligo hotspots" "wt SPO11-1" "kss SPO11-1" wt_SPO11oligo kss_SPO11oligo 10000

library(ggplot2)
library(ggthemes)

dataName <- "SPO11-1-oligo hotspots" 
pt1Name <- "wt SPO11-1" 
pt2Name <-  "kss SPO11-1" 
pt1LibName <- "wt_SPO11oligo"
pt2LibName <- "kss_SPO11oligo"
# Number of permutations (randomisations) performed
perms <- 10000

args <- commandArgs(trailingOnly = T)
dataName <- as.character(args[1])
pt1Name <- as.character(args[2])
pt2Name <- as.character(args[3])
pt1LibName <- as.character(args[4])
pt2LibName <- as.character(args[5])
# Number of permutations (randomisations) performed
perms <- as.numeric(args[6])

REC8Dir <-  "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1/genome_wide/regioneR/noMinWidth_mergedOverlaps/"
inDir1 <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/regioneR/"
inDir2 <- "/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI34_RPI35_profiles_rangerPeaks_idr0.05/genome_wide/regioneR/"
plotDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/regioneR/plots/"

famNames <- c("dna", "heli", "ptmari", "mudr", "enspm", "hat", "harbinger",
              "rna", "gypsy", "copia", "linel1", "sine")
famNamesPlot <- c("DNA", "Helitrons", "Pogo/Tc1/Mariner", "MuDR", "EnSpm", "hAT", "Harbinger",
                  "RNA", "Gypsy LTR", "Copia LTR", "LINE-1", "SINE")

load(paste0(REC8Dir, "permTest_REC8_HA_Rep1_rangerPeaks_vs_TEsDNA.RData"))
ptREC8_DNA <- ptPeaksTEsDNAPerChrom
ptPeaksTEsDNAPerChrom <- NULL
load(paste0(REC8Dir, "permTest_REC8_HA_Rep1_rangerPeaks_vs_TEsRNA.RData"))
ptREC8_RNA <- ptPeaksTEsRNAPerChrom
ptPeaksTEsRNAPerChrom <- NULL
ptREC8 <- c(ptREC8_DNA, ptREC8_RNA)

load(paste0(inDir1, "permTest_SPO11_RPI1_RPI8_rangerPeaks_vs_TEsDNA.RData"))
pt1_DNA <- ptPeaksTEsDNAPerChrom
ptPeaksTEsDNAPerChrom <- NULL
load(paste0(inDir1, "permTest_SPO11_RPI1_RPI8_rangerPeaks_vs_TEsRNA.RData"))
pt1_RNA <- ptPeaksTEsRNAPerChrom
ptPeaksTEsRNAPerChrom <- NULL
pt1 <- c(pt1_DNA, pt1_RNA)

load(paste0(inDir2, "permTest_kss_SPO11_RPI34_RPI35_rangerPeaks_vs_TEsDNA.RData"))
pt2_DNA <- ptPeaksTEsDNAPerChrom
ptPeaksTEsDNAPerChrom <- NULL
load(paste0(inDir2, "permTest_kss_SPO11_RPI34_RPI35_rangerPeaks_vs_TEsRNA.RData"))
pt2_RNA <- ptPeaksTEsRNAPerChrom
ptPeaksTEsRNAPerChrom <- NULL
pt2 <- c(pt2_DNA, pt2_RNA)

# ptREC8
ptREC8_Pval <- lapply(seq_along(ptREC8), function(x) {
  ptREC8[[x]]$numOverlaps$pval
})
ptREC8_Obs <- lapply(seq_along(ptREC8), function(x) {
  ptREC8[[x]]$numOverlaps$observed
})
ptREC8_Perm <- lapply(seq_along(ptREC8), function(x) {
  ptREC8[[x]]$numOverlaps$permuted
})
ptREC8_Exp <- lapply(seq_along(ptREC8), function(x) {
  mean(ptREC8[[x]]$numOverlaps$permuted)
})
ptREC8_log2ObsExp <- lapply(seq_along(ptREC8_Obs), function(x) {
  log2(ptREC8_Obs[[x]]/ptREC8_Exp[[x]])
})
ptREC8_Zscore <- lapply(seq_along(ptREC8), function(x) {
  ptREC8[[x]]$numOverlaps$zscore
})
ptREC8_AltHyp <- lapply(seq_along(ptREC8), function(x) {
  ptREC8[[x]]$numOverlaps$alternative
})
ptREC8_alpha0.05 <- lapply(seq_along(ptREC8_Perm), function(x) {
  if(ptREC8_AltHyp[[x]] == "greater") {
    quantile(ptREC8_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(ptREC8_Perm[[x]], 0.05)[[1]]
  }
})
ptREC8_log2alpha0.05 <- lapply(seq_along(ptREC8_alpha0.05), function(x) {
  log2(ptREC8_alpha0.05[[x]]/ptREC8_Exp[[x]])
})

ptREC8_log2ObsExp_sorted <- sort.int(unlist(ptREC8_log2ObsExp), decreasing = T)
ptREC8_log2alpha0.05_sorted <- unlist(ptREC8_log2alpha0.05[sort.int(unlist(ptREC8_log2ObsExp), decreasing = T, index.return = T)$ix])
ptREC8_famNames_sorted <- famNames[sort.int(unlist(ptREC8_log2ObsExp), decreasing = T, index.return = T)$ix]
ptREC8_famNamesPlot_sorted <- famNamesPlot[sort.int(unlist(ptREC8_log2ObsExp), decreasing = T, index.return = T)$ix]

# pt1
pt1_Pval <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$pval
})
pt1_Obs <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$observed
})
pt1_Perm <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$permuted
})
pt1_Exp <- lapply(seq_along(pt1), function(x) {
  mean(pt1[[x]]$numOverlaps$permuted)
})
pt1_log2ObsExp <- lapply(seq_along(pt1_Obs), function(x) {
  log2(pt1_Obs[[x]]/pt1_Exp[[x]])
})
pt1_Zscore <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$zscore
})
pt1_AltHyp <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$alternative
})
pt1_alpha0.05 <- lapply(seq_along(pt1_Perm), function(x) {
  if(pt1_AltHyp[[x]] == "greater") {
    quantile(pt1_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt1_Perm[[x]], 0.05)[[1]]
  }
})
pt1_log2alpha0.05 <- lapply(seq_along(pt1_alpha0.05), function(x) {
  log2(pt1_alpha0.05[[x]]/pt1_Exp[[x]])
})

pt1_log2ObsExp_sorted <- unlist(pt1_log2ObsExp[sort.int(unlist(ptREC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt1_log2alpha0.05_sorted <- unlist(pt1_log2alpha0.05[sort.int(unlist(ptREC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt1_famNames_sorted <- famNames[sort.int(unlist(ptREC8_log2ObsExp), decreasing = T, index.return = T)$ix]
pt1_famNamesPlot_sorted <- famNamesPlot[sort.int(unlist(ptREC8_log2ObsExp), decreasing = T, index.return = T)$ix]

# pt2
pt2_Pval <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$pval
})
pt2_Obs <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$observed
})
pt2_Perm <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$permuted
})
pt2_Exp <- lapply(seq_along(pt2), function(x) {
  mean(pt2[[x]]$numOverlaps$permuted)
})
pt2_log2ObsExp <- lapply(seq_along(pt2_Obs), function(x) {
  log2(pt2_Obs[[x]]/pt2_Exp[[x]])
})
pt2_Zscore <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$zscore
})
pt2_AltHyp <- lapply(seq_along(pt2), function(x) {
  pt2[[x]]$numOverlaps$alternative
})
pt2_alpha0.05 <- lapply(seq_along(pt2_Perm), function(x) {
  if(pt2_AltHyp[[x]] == "greater") {
    quantile(pt2_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt2_Perm[[x]], 0.05)[[1]]
  }
})
pt2_log2alpha0.05 <- lapply(seq_along(pt2_alpha0.05), function(x) {
  log2(pt2_alpha0.05[[x]]/pt2_Exp[[x]])
})

pt2_log2ObsExp_sorted <- unlist(pt2_log2ObsExp[sort.int(unlist(ptREC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt2_log2alpha0.05_sorted <- unlist(pt2_log2alpha0.05[sort.int(unlist(ptREC8_log2ObsExp), decreasing = T, index.return = T)$ix])
pt2_famNames_sorted <- famNames[sort.int(unlist(ptREC8_log2ObsExp), decreasing = T, index.return = T)$ix]
pt2_famNamesPlot_sorted <- famNamesPlot[sort.int(unlist(ptREC8_log2ObsExp), decreasing = T, index.return = T)$ix]

df <- data.frame(Sample = rep(c(pt1Name, pt2Name),
                              each = length(pt1_log2ObsExp_sorted)),
                 Transposon_family = rep(pt1_famNamesPlot_sorted, 2),
                 log2ObsExp = c(pt1_log2ObsExp_sorted,
                                pt2_log2ObsExp_sorted),
                 log2alpha0.05 = c(pt1_log2alpha0.05_sorted, pt2_log2alpha0.05_sorted))

df$Transposon_family <- factor(df$Transposon_family,
                               levels = pt1_famNamesPlot_sorted)
df$Sample <- factor(df$Sample,
                    levels = c(pt1Name, pt2Name))

bp <- ggplot(data = df,
             mapping = aes(x = Transposon_family,
                           y = log2ObsExp,
                           fill = Sample)) +
      geom_bar(stat = "identity",
               position = position_dodge()) +
      scale_fill_manual(name = "Genotype",
                        values = c("black",
                                   "dodgerblue3"),
                        labels = c(pt1Name,
                                   pt2Name)) +
      geom_point(mapping = aes(Transposon_family, log2alpha0.05),
                 position = position_dodge(0.9),
                 shape = "-", colour  = "red", size = 4.25) +
      labs(x = "Transposon family",
           y = expression("Log"[2]*"(observed:expected) peak overlap")) +
      theme_bw() +
      theme(axis.line.y = element_line(size = 0.5, colour = "black"),
            axis.ticks.y = element_line(size = 0.25, colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 7),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0(dataName, " (", as.character(perms), " permutations)"))
ggsave(paste0(plotDir, "barplot_TE_families_permTestResults_",
              as.character(perms), "perms_",
              "log2_Observed_Expected_", pt1LibName, "_", pt2LibName, "_peaks.pdf"),
       plot = bp,
       height = 4, width = 5)
save(bp,
     file = paste0(plotDir, "barplot_TE_families_permTestResults_",
                   as.character(perms), "perms_",
                  "log2_Observed_Expected_", pt1LibName, "_", pt2LibName, "_peaks.RData"))

 
