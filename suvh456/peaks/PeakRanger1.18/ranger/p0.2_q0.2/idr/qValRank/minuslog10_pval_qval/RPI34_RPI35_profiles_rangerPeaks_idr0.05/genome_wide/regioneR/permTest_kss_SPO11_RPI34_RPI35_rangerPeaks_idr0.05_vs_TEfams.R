#!/applications/R/R-3.3.2/bin/Rscript

# Use permutation test function in regioneR to determine if peaks overlap features
# of interest (e.g., other peaks, genes, TEs) more or less than expected by chance

# Usage on hydrogen node7:
# csmit -m 50G -c 47 "Rscript permTest_kss_SPO11_RPI34_RPI35_rangerPeaks_idr0.05_vs_TEfams.R"

library(regioneR)

# Genomic definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))

inDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/"

outDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI34_RPI35_profiles_rangerPeaks_idr0.05/genome_wide/regioneR/"
DNAplotDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI34_RPI35_profiles_rangerPeaks_idr0.05/genome_wide/regioneR/plots/TEsDNA/"
RNAplotDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI34_RPI35_profiles_rangerPeaks_idr0.05/genome_wide/regioneR/plots/TEsRNA/"

# kss SPO11-1-oligo peaks
load(paste0(inDir,
            "armrangerPeaksGRmerge_RPI34_RPI35_idr0.05_noMinWidth.RData"))
load(paste0(inDir,
            "perirangerPeaksGRmerge_RPI34_RPI35_idr0.05_noMinWidth.RData"))
peaksGR <- sort(c(armrangerPeaksGRmerge, perirangerPeaksGRmerge))
strand(peaksGR) <- "*"
print("***********peaks***********")
print(peaksGR)

DNAfamNames <- c("dna", "heli", "ptmari", "mudr", "enspm", "hat", "harbinger")
RNAfamNames <- c("rna", "gypsy", "copia", "linel1", "sine")
DNAdir <- "/projects/ajt200/TAIR10/TE_classes/DNA/"
RNAdir <- "/projects/ajt200/TAIR10/TE_classes/RNA/"

### DNA TEs

TEsDNAGR <- lapply(seq_along(DNAfamNames), function(x) {
  TEsDNA <- read.table(file = paste0(DNAdir, "TAIR10_Buisine_TEs_strand_tab_ann_", DNAfamNames[x], ".txt"), header = T)
  GRanges(seqnames = TEsDNA$chr, ranges = IRanges(start = TEsDNA$start, end = TEsDNA$end), strand = "*")
})

# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions in B as in A
set.seed(123)
ptPeaksTEsDNAPerChrom <- lapply(seq_along(TEsDNAGR), function(x) {
  permTest(A = peaksGR, B = TEsDNAGR[[x]], genome = genome,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE, per.chromosome = TRUE,
           evaluate.function = numOverlaps, count.once = TRUE,
           ntimes = 10000, mc.set.seed = FALSE, mc.cores = 47)
})

for(i in 1:length(ptPeaksTEsDNAPerChrom)) {
  assign(paste0(DNAfamNames[i]), ptPeaksTEsDNAPerChrom[[i]])
}
save(ptPeaksTEsDNAPerChrom,
     file = paste0(outDir,
                   "permTest_kss_SPO11_RPI34_RPI35_rangerPeaks_vs_TEsDNA.RData"))

# Summarise results in a table
noOfFeatures <- NULL
expected <- NULL
observed <- NULL
pval <- NULL
zscore <- NULL
for(i in 1:length(ptPeaksTEsDNAPerChrom)) {
  noOfFeaturesi <- print(length(TEsDNAGR[[i]]))
  noOfFeatures <- c(noOfFeatures, noOfFeaturesi)
  expectedi <- print(round(mean(ptPeaksTEsDNAPerChrom[[i]]$numOverlaps$permuted)))
  expected <- c(expected, expectedi)
  observedi <- print(ptPeaksTEsDNAPerChrom[[i]]$numOverlaps$observed)
  observed <- c(observed, observedi)
  pvali <- print(round(ptPeaksTEsDNAPerChrom[[i]]$numOverlaps$pval, 4))
  pval <- c(pval, pvali)
  zscorei <- print(round(ptPeaksTEsDNAPerChrom[[i]]$numOverlaps$zscore, 4))
  zscore <- c(zscore, zscorei)
}
ptPeaksTEsDNAPerChromDataFrame <- cbind(noOfFeatures, expected, observed, pval, zscore)
write.table(ptPeaksTEsDNAPerChromDataFrame,
            file = paste0(outDir,
                          "permTest_kss_SPO11_RPI34_RPI35_rangerPeaks_vs_TEsDNA_DataFrame.txt"),
            sep = "\t", row.names = F)

# plot graphical summaries of results
for(i in 1:length(ptPeaksTEsDNAPerChrom)) {
  pdf(file = paste0(DNAplotDir, DNAfamNames[i], "_permTest_nperm10000_kss_SPO11_RPI34_RPI35_rangerPeaks_perChrom.pdf"), width = 10, height = 7)
  plot(ptPeaksTEsDNAPerChrom[[i]], main = paste0("kss_SPO11_RPI34_RPI35_rangerPeaks vs ", DNAfamNames[i]), xlab = "Number of overlaps", ylab = "Density")
  dev.off()

  # Using the localZScore() function, evaluate whether the association between peaks and other is highly dependent on their exact position
  lz_1kb <- localZScore(pt = ptPeaksTEsDNAPerChrom[[i]], A = peaksGR, B = TEsDNAGR[[i]],
                        window = 1000, step = 50, count.once = TRUE)
  lz_10kb <- localZScore(pt = ptPeaksTEsDNAPerChrom[[i]], A = peaksGR, B = TEsDNAGR[[i]],
                         window = 10000, step = 500, count.once = TRUE)
  lz_custom <- localZScore(pt = ptPeaksTEsDNAPerChrom[[i]], A = peaksGR, B = TEsDNAGR[[i]],
                           window = 10*mean(width(peaksGR)), step = mean(width(peaksGR))/2, count.once = TRUE)
  win <- as.character(round((10*mean(width(peaksGR)))/1000))
  step <- as.character(round(mean(width(peaksGR))/2))
  pdf(file = paste0(DNAplotDir, DNAfamNames[i], "_localZscore_permTest_nperm10000_kss_SPO11_RPI34_RPI35_rangerPeaks_w1kb_s50bp_w10kb_s500bp_w", win ,"kb_s", step, "bp_perChrom.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(lz_1kb, main = paste0("kss_SPO11_RPI34_RPI35_rangerPeaks vs ", DNAfamNames[i], " (1-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("kss_SPO11_RPI34_RPI35_rangerPeaks vs ", DNAfamNames[i], " (1-kb shift)"))
  plot(lz_10kb, main = paste0("kss_SPO11_RPI34_RPI35_rangerPeaks vs ", DNAfamNames[i], " (10-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("kss_SPO11_RPI34_RPI35_rangerPeaks vs ", DNAfamNames[i], " (10-kb shift)"))
  plot(lz_custom, main = paste0("kss_SPO11_RPI34_RPI35_rangerPeaks vs ", DNAfamNames[i], " (~", win, "-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("kss_SPO11_RPI34_RPI35_rangerPeaks vs ", DNAfamNames[i], " (~", win, "-kb shift)"))
  dev.off()
}


### RNA TEs

TEsRNAGR <- lapply(seq_along(RNAfamNames), function(x) {
  TEsRNA <- read.table(file = paste0(RNAdir, "TAIR10_Buisine_TEs_strand_tab_ann_", RNAfamNames[x], ".txt"), header = T)
  GRanges(seqnames = TEsRNA$chr, ranges = IRanges(start = TEsRNA$start, end = TEsRNA$end), strand = "*")
})

# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions in B as in A
set.seed(123)
ptPeaksTEsRNAPerChrom <- lapply(seq_along(TEsRNAGR), function(x) {
  permTest(A = peaksGR, B = TEsRNAGR[[x]], genome = genome,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE, per.chromosome = TRUE,
           evaluate.function = numOverlaps, count.once = TRUE,
           ntimes = 10000, mc.set.seed = FALSE, mc.cores = 47)
})

for(i in 1:length(ptPeaksTEsRNAPerChrom)) {
  assign(paste0(RNAfamNames[i]), ptPeaksTEsRNAPerChrom[[i]])
}
save(ptPeaksTEsRNAPerChrom,
     file = paste0(outDir,
                   "permTest_kss_SPO11_RPI34_RPI35_rangerPeaks_vs_TEsRNA.RData"))

# Summarise results in a table
noOfFeatures <- NULL
expected <- NULL
observed <- NULL
pval <- NULL
zscore <- NULL
for(i in 1:length(ptPeaksTEsRNAPerChrom)) {
  noOfFeaturesi <- print(length(TEsRNAGR[[i]]))
  noOfFeatures <- c(noOfFeatures, noOfFeaturesi)
  expectedi <- print(round(mean(ptPeaksTEsRNAPerChrom[[i]]$numOverlaps$permuted)))
  expected <- c(expected, expectedi)
  observedi <- print(ptPeaksTEsRNAPerChrom[[i]]$numOverlaps$observed)
  observed <- c(observed, observedi)
  pvali <- print(round(ptPeaksTEsRNAPerChrom[[i]]$numOverlaps$pval, 4))
  pval <- c(pval, pvali)
  zscorei <- print(round(ptPeaksTEsRNAPerChrom[[i]]$numOverlaps$zscore, 4))
  zscore <- c(zscore, zscorei)
}
ptPeaksTEsRNAPerChromDataFrame <- cbind(noOfFeatures, expected, observed, pval, zscore)
write.table(ptPeaksTEsRNAPerChromDataFrame,
            file = paste0(outDir,
                          "permTest_kss_SPO11_RPI34_RPI35_rangerPeaks_vs_TEsRNA_DataFrame.txt"),
            sep = "\t", row.names = F)

# plot graphical summaries of results
for(i in 1:length(ptPeaksTEsRNAPerChrom)) {
  pdf(file = paste0(RNAplotDir, RNAfamNames[i], "_permTest_nperm10000_kss_SPO11_RPI34_RPI35_rangerPeaks_perChrom.pdf"), width = 10, height = 7)
  plot(ptPeaksTEsRNAPerChrom[[i]], main = paste0("kss_SPO11_RPI34_RPI35_rangerPeaks vs ", RNAfamNames[i]), xlab = "Number of overlaps", ylab = "Density")
  dev.off()

  # Using the localZScore() function, evaluate whether the association between peaks and other is highly dependent on their exact position
  lz_1kb <- localZScore(pt = ptPeaksTEsRNAPerChrom[[i]], A = peaksGR, B = TEsRNAGR[[i]],
                        window = 1000, step = 50, count.once = TRUE)
  lz_10kb <- localZScore(pt = ptPeaksTEsRNAPerChrom[[i]], A = peaksGR, B = TEsRNAGR[[i]],
                         window = 10000, step = 500, count.once = TRUE)
  lz_custom <- localZScore(pt = ptPeaksTEsRNAPerChrom[[i]], A = peaksGR, B = TEsRNAGR[[i]],
                           window = 10*mean(width(peaksGR)), step = mean(width(peaksGR))/2, count.once = TRUE)
  win <- as.character(round((10*mean(width(peaksGR)))/1000))
  step <- as.character(round(mean(width(peaksGR))/2))
  pdf(file = paste0(RNAplotDir, RNAfamNames[i], "_localZscore_permTest_nperm10000_kss_SPO11_RPI34_RPI35_rangerPeaks_w1kb_s50bp_w10kb_s500bp_w", win ,"kb_s", step, "bp_perChrom.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(lz_1kb, main = paste0("kss_SPO11_RPI34_RPI35_rangerPeaks vs ", RNAfamNames[i], " (1-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("kss_SPO11_RPI34_RPI35_rangerPeaks vs ", RNAfamNames[i], " (1-kb shift)"))
  plot(lz_10kb, main = paste0("kss_SPO11_RPI34_RPI35_rangerPeaks vs ", RNAfamNames[i], " (10-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("kss_SPO11_RPI34_RPI35_rangerPeaks vs ", RNAfamNames[i], " (10-kb shift)"))
  plot(lz_custom, main = paste0("kss_SPO11_RPI34_RPI35_rangerPeaks vs ", RNAfamNames[i], " (~", win, "-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("kss_SPO11_RPI34_RPI35_rangerPeaks vs ", RNAfamNames[i], " (~", win, "-kb shift)"))
  dev.off()
}

