#!/applications/R/R-3.3.2/bin/Rscript

# Use permutation test function in regioneR to determine if peaks overlap features
# of interest (e.g., other peaks, genes, TEs) more or less than expected by chance

# Usage on hydrogen node7:
# csmit -m 50G -c 48 "/applications/R/R-3.3.2/bin/Rscript permTest_SPO11_RPI1_RPI8_rangerPeaks_idr0.05_vs_MSH4_Rep1_peaks.R"

library(regioneR)

# Genomic definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))

inDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/"
outDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/regioneR/"
plotDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/regioneR/plots/"

# MSH4_Rep1
load(paste0("/home/ajt200/analysis/160902_Sasha_ChIP_MSH4_Rep1/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
            "MSH4_Rep1_armrangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
load(paste0("/home/ajt200/analysis/160902_Sasha_ChIP_MSH4_Rep1/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
            "MSH4_Rep1_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
MSH4_Rep1GR <- sort(c(armrangerPeaksGRmergedOverlaps, perirangerPeaksGRmergedOverlaps))
armrangerPeaksGRmergedOverlaps <- NULL
perirangerPeaksGRmergedOverlaps <- NULL
strand(MSH4_Rep1GR) <- "*"
print("***********peaks***********")
print(MSH4_Rep1GR)
print(length(MSH4_Rep1GR))
#[1] 40556

# SPO11-1-oligo hotspots
load(paste0(inDir,
            "armrangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData"))
load(paste0(inDir,
            "perirangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData"))
peaksGR <- sort(c(armrangerPeaksGRmerge, perirangerPeaksGRmerge))
armrangerPeaksGRmerge <- NULL
perirangerPeaksGRmerge <- NULL
strand(peaksGR) <- "*"
print(length(peaksGR))
#[1] 5914

otherNames <- c("MSH4_Rep1GR")

grl <- GRangesList("MSH4_Rep1GR" = MSH4_Rep1GR)

# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions in B as in A
set.seed(123)
ptPeaksOtherPerChrom <- lapply(seq_along(grl), function(x) {
  permTest(A = peaksGR, B = grl[[x]], genome = genome,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE, per.chromosome = TRUE,
           evaluate.function = numOverlaps, count.once = TRUE,
           ntimes = 10000, mc.set.seed = FALSE, mc.cores = detectCores())
})

for(i in 1:length(ptPeaksOtherPerChrom)) {
  assign(paste0(otherNames[i]), ptPeaksOtherPerChrom[[i]])
}
save(ptPeaksOtherPerChrom,
     file = paste0(outDir,
                   "permTest_SPO11_RPI1_RPI8_rangerPeaks_vs_MSH4_Rep1.RData"))

# Summarise results in a table
noOfFeatures <- NULL
expected <- NULL
observed <- NULL
pval <- NULL
zscore <- NULL
for(i in 1:length(ptPeaksOtherPerChrom)) {
  noOfFeaturesi <- print(length(grl[[i]]))
  noOfFeatures <- c(noOfFeatures, noOfFeaturesi)
  expectedi <- print(round(mean(ptPeaksOtherPerChrom[[i]]$numOverlaps$permuted)))
  expected <- c(expected, expectedi)
  observedi <- print(ptPeaksOtherPerChrom[[i]]$numOverlaps$observed)
  observed <- c(observed, observedi)
  pvali <- print(round(ptPeaksOtherPerChrom[[i]]$numOverlaps$pval, 4))
  pval <- c(pval, pvali)
  zscorei <- print(round(ptPeaksOtherPerChrom[[i]]$numOverlaps$zscore, 4))
  zscore <- c(zscore, zscorei)
}
ptPeaksOtherPerChromDataFrame <- cbind(noOfFeatures, expected, observed, pval, zscore)
write.table(ptPeaksOtherPerChromDataFrame,
            file = paste0(outDir,
                          "permTest_SPO11_RPI1_RPI8_rangerPeaks_vs_MSH4_Rep1_DataFrame.txt"),
            sep = "\t", row.names = F)

# plot graphical summaries of results
for(i in 1:length(ptPeaksOtherPerChrom)) {
  pdf(file = paste0(plotDir, otherNames[i], "_permTest_nperm10000_SPO11_RPI1_RPI8_rangerPeaks_perChrom.pdf"), width = 10, height = 7)
  plot(ptPeaksOtherPerChrom[[i]], main = paste0("SPO11_RPI1_RPI8_rangerPeaks vs ", otherNames[i]), xlab = "Number of overlaps", ylab = "Density")
  dev.off()

  # Using the localZScore() function, evaluate whether the association between peaks and other is highly dependent on their exact position
  lz_1kb <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = grl[[i]],
                        window = 1000, step = 50, count.once = TRUE)
  lz_10kb <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = grl[[i]],
                         window = 10000, step = 500, count.once = TRUE)
  lz_custom <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = grl[[i]],
                           window = 10*mean(width(peaksGR)), step = mean(width(peaksGR))/2, count.once = TRUE)
  win <- as.character(round((10*mean(width(peaksGR)))/1000))
  step <- as.character(round(mean(width(peaksGR))/2))
  pdf(file = paste0(plotDir, otherNames[i], "_localZscore_permTest_nperm10000_SPO11_RPI1_RPI8_rangerPeaks_w1kb_s50bp_w10kb_s500bp_w", win ,"kb_s", step, "bp_perChrom.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(lz_1kb, main = paste0("SPO11_RPI1_RPI8_rangerPeaks vs ", otherNames[i], " (1-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("SPO11_RPI1_RPI8_rangerPeaks vs ", otherNames[i], " (1-kb shift)"))
  plot(lz_10kb, main = paste0("SPO11_RPI1_RPI8_rangerPeaks vs ", otherNames[i], " (10-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("SPO11_RPI1_RPI8_rangerPeaks vs ", otherNames[i], " (10-kb shift)"))
  plot(lz_custom, main = paste0("SPO11_RPI1_RPI8_rangerPeaks vs ", otherNames[i], " (~", win, "-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("SPO11_RPI1_RPI8_rangerPeaks vs ", otherNames[i], " (~", win, "-kb shift)"))
  dev.off()
}


