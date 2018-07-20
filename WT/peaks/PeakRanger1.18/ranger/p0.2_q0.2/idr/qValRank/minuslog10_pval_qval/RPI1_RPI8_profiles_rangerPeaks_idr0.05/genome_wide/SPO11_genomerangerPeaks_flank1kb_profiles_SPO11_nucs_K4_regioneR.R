####################################################################################################################################################
# Calculate and plot mean SPO11-1-oligos, nucleosomes and H3K4me3 coverage profiles and base frequency profiles around SPO11-1-oligos peaks        #
####################################################################################################################################################

library(segmentSeq)
library(parallel)
library(EnrichedHeatmap)
library(genomation)
library(regioneR)
#library(circlize)
#library(RColorBrewer)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
seqlevels(genome) <- sub("Chr", "", seqlevels(genome))

inDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/"
matDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/matrices/"
plotDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/plots/"
histDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/hist/"

load(paste0(inDir, "armrangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData"))
load(paste0(inDir, "perirangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData"))
peaksGR <- sort(c(armrangerPeaksGRmerge, perirangerPeaksGRmerge))
armrangerPeaksGRmerge <- NULL
perirangerPeaksGRmerge <- NULL
strand(peaksGR) <- "*"
seqlevels(peaksGR) <- sub("Chr", "", seqlevels(peaksGR))
#peaksGR <- peaksGR[width(peaksGR) <= 5000]
ranLocGR <- randomizeRegions(peaksGR, genome = genome, per.chromosome = TRUE, allow.overlaps = TRUE)

print("print(length(peaksGR))")
print(length(peaksGR))
print("print(mean(width(peaksGR)))")
print(mean(width(peaksGR)))
print("print(median(width(peaksGR)))")
print(median(width(peaksGR)))
print("print(range(width(peaksGR)))")
print(range(width(peaksGR)))

pdf(paste0(histDir, "hist_SPO11_genomerangerPeaksGR_RPI1_RPI8_idr0.05_width_brks250.pdf"))
par(mfrow = c(1,1))
# Plot histogram of SPO11 peak widths
hist(width(peaksGR), breaks = 250, col = "black", ann = FALSE)
mtext(side = c(1, 2), line = c(3.5, 2.5), text = c(expression(atop("Hotspot width", "(SPO11-1 RPI1 and RPI8; IDR <= 0.05)")), expression("Hotspots")))
dev.off()

SPO11oligos <- system("ls /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab.bed", intern = T)
H3K4me3 <- system("ls /projects/ajt200/BAM_masters/H3K4me3/WT/coverage/log2ChIPinput/log2wtH3K4me3ChIPwtH3K9me2input_norm_allchrs_coverage_coord_tab.bed", intern = T)
NUC <- system("ls /projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_trim51_input/log2ChIPinput/log2wtNucNakedDNA_norm_allchrs_coverage_coord_tab.bed", intern = T)

# Create lists/vectors of path and library names
libPaths <- list(SPO11oligos, H3K4me3, NUC)
libNames <- c("SPO11oligos", "H3K4me3", "NUC")

# Read in coverage files as GRanges objects and assign to library names
grTmp <- mclapply(seq_along(libPaths), function(x) {
  readGeneric(libPaths[[x]], meta.col = list(coverage = 4))
}, mc.cores = 3, mc.preschedule = F)
for(i in 1:length(grTmp)) {
  seqlevels(grTmp[[i]]) <- sub("Chr", "", seqlevels(grTmp[[i]]))
  assign(paste0(libNames[i]), grTmp[[i]])
}

# Create GRangesList object containing per base coverage for each library
grl <- GRangesList("SPO11oligos" = SPO11oligos, "H3K4me3" = H3K4me3, "NUC" = NUC)

rm(grTmp, SPO11oligos, H3K4me3, NUC)
gc()

winSize <- 1
targetSize <- mean(width(peaksGR))
flankSize <- 1000

# Function to create coverage matrices for target loci and random loci (incl. flanking regions)
## and to calculate mean levels per window
covMatrix <- function(signal, target, ranLoc, x) {
  #target loci
  set.seed(2840)
  mat1 <- normalizeToMatrix(signal, target, value_column = "coverage",
                            extend = flankSize, mean_mode = "absolute", w = winSize,
                            empty_value = 0, smooth = FALSE,
                            include_target = TRUE, target_ratio = targetSize/(targetSize+(flankSize*2)))
  print(mat1)
  print(length(mat1))
  #print("mat1 rows = ")
  #print(length(mat1)/(targetSize+(flankSize*2)))
  
  #mat1_failed_rows <- attr(mat1, "failed_rows")
  #print("mat1 failed rows = ")
  #print(length(mat1_failed_rows))
  #mat1 <- mat1[-mat1_failed_rows,]
  #print(mat1)
  #print(length(mat1))
  #print("mat1 rows less failed rows = ")
  #print(length(mat1)/(targetSize+(flankSize*2)))
  
  #numberPeaks <- print(length(mat1)/(targetSize+(flankSize*2)))
  numberPeaks <- length(peaksGR)
  mat1_DF <- data.frame(mat1)
  mat1_DF_colMeans <- as.vector(colMeans(mat1_DF))

##save(mat1, file = outFile[[x]][[1]])
  #save(mat1_failed_rows, file = outFailedRows[[x]][[1]])
  save(numberPeaks, file = outNoLoc[[x]][[1]])
##write.table(mat1_DF, file = outDF[[x]][[1]])
  write.table(mat1_DF_colMeans, file = outDFCM[[x]][[1]])

  #random loci
  set.seed(8472)
  mat2 <- normalizeToMatrix(signal, ranLoc, value_column = "coverage",
                            extend = flankSize, mean_mode = "absolute", w = winSize,
                            empty_value = 0, smooth = FALSE,
                            include_target = TRUE, target_ratio = targetSize/(targetSize+(flankSize*2)))
  print(mat2)
  print(length(mat2))
  #print("mat2 rows = ")
  #print(length(mat2)/(targetSize+(flankSize*2)))
  #mat2_failed_rows <- attr(mat2, "failed_rows")
  #print("mat2 failed rows = ")
  #print(length(mat2_failed_rows))
  #mat2 <- mat2[-mat2_failed_rows,]
  #print(mat2)
  #print(length(mat2))
  #print("mat2 rows less failed rows = ")
  #print(length(mat2)/(targetSize+(flankSize*2)))
 
  #numberRanLoc <- print(length(mat2)/(targetSize+(flankSize*2)))
  numberRanLoc <- length(ranLocGR)
  mat2_DF <- data.frame(mat2)
  mat2_DF_colMeans <- as.vector(colMeans(mat2_DF))

##save(mat2, file = outFile[[x]][[2]])
  #save(mat2_failed_rows, file = outFailedRows[[x]][[2]])
  save(numberRanLoc, file = outNoLoc[[x]][[2]])
##write.table(mat2_DF, file = outDF[[x]][[2]])
  write.table(mat2_DF_colMeans, file = outDFCM[[x]][[2]])
}

# Run covMatrix() on log2(ChIP/input) coverage GRanges objects to obtain matrices containing
## normalised methylation values around SPO11-1-oligos peaks
#outFile <- lapply(seq_along(libNames), function(x)
#  list(paste0(matDir, libNames[[x]],
#              "_norm_cov_SPO11peaks_mat1_target_and_flank.RData"),
#       paste0(matDir, libNames[[x]],
#              "_norm_cov_ranLoc_mat2_target_and_flank.RData")))
#outFailedRows <- lapply(seq_along(libNames), function(x)
#  list(paste0(matDir, libNames[[x]],
#              "_norm_cov_SPO11peaks_mat1_target_and_flank_failed_rows.RData"),
#       paste0(matDir, libNames[[x]],
#              "_norm_cov_ranLoc_mat2_target_and_flank_failed_rows.RData")))
outNoLoc <- lapply(seq_along(libNames), function(x)
  list(paste0(matDir, libNames[[x]],
              "_norm_cov_SPO11peaks_mat1_target_and_flank_numberPeaks.RData"),
       paste0(matDir, libNames[[x]],
              "_norm_cov_ranLoc_mat2_target_and_flank_numberRanLoc.RData")))
#outDF <- lapply(seq_along(libNames), function(x)
#  list(paste0(matDir, libNames[[x]],
#              "_norm_cov_SPO11peaks_mat1_target_and_flank_dataframe.txt"),
#       paste0(matDir, libNames[[x]],
#              "_norm_cov_ranLoc_mat2_target_and_flank_dataframe.txt")))
outDFCM <- lapply(seq_along(libNames), function(x)
  list(paste0(matDir, libNames[[x]],
              "_norm_cov_SPO11peaks_mat1_target_and_flank_dataframe_colMeans.txt"),
       paste0(matDir, libNames[[x]],
              "_norm_cov_ranLoc_mat2_target_and_flank_dataframe_colMeans.txt")))

mclapply(seq_along(grl), function(x) {
  covMatrix(grl[[x]], peaksGR, ranLocGR, x)
}, mc.cores = 3, mc.preschedule = F)

# Plot SPO11oligos coverage and other coverage profiles around SPO11oligos peaks and random loci

SPO11oligosDat <- read.table(file = paste0(matDir, "SPO11oligos_norm_cov_SPO11peaks_mat1_target_and_flank_dataframe_colMeans.txt"))
SPO11oligosRanDat <- read.table(file = paste0(matDir, "SPO11oligos_norm_cov_ranLoc_mat2_target_and_flank_dataframe_colMeans.txt"))

## Function to overlay SPO11oligos with other log2-transformed coverage profiles in separate plots
SPO11oligosvsOthersPlot <- function(xplot, otherName, otherYlabel) {
  load(file = paste0(matDir, otherName, "_norm_cov_SPO11peaks_mat1_target_and_flank_numberPeaks.RData"))
  load(file = paste0(matDir, otherName, "_norm_cov_ranLoc_mat2_target_and_flank_numberRanLoc.RData"))
  otherDat <- read.table(file = paste0(matDir, otherName, "_norm_cov_SPO11peaks_mat1_target_and_flank_dataframe_colMeans.txt"))
  otherRanDat <- read.table(file = paste0(matDir, otherName, "_norm_cov_ranLoc_mat2_target_and_flank_dataframe_colMeans.txt"))

  plot(xplot, SPO11oligosDat[,1], ylim = c(min(SPO11oligosDat[,1], SPO11oligosRanDat[,1]), max(SPO11oligosDat[,1], SPO11oligosRanDat[,1])), type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(SPO11oligosDat[,1], SPO11oligosRanDat[,1])))
  mtext(side = 2, line = 2, cex = 0.8, text = "SPO11-1", col = mycols[1])
  #mtext(side = 3, line = 1, text = paste0("Euchromatic SPO11-1 peaks (n = ", print(numberPeaks), ")"))
  par(new = T)
  plot(xplot, otherDat[,1], ylim = c(min(otherDat[,1], otherRanDat[,1]), max(otherDat[,1], otherRanDat[,1])), type = "l", lwd = 1.5, col = mycols[2], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(otherDat[,1], otherRanDat[,1])))
  #mtext(side = 4, line = 2, cex = 0.8, text = otherYlabel, col = mycols[2])
  axis(side = 1, at = c(0, flankSize, length(SPO11oligosDat[,1])-flankSize, length(SPO11oligosDat[,1])), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(0, flankSize, length(SPO11oligosDat[,1])-flankSize, length(SPO11oligosDat[,1])), text = c("-1 kb", "Start", "End", "+1 kb")) 
  #mtext(side = 1, line = c(3.5), cex = 0.8, text = expression(atop(paste("Position relative to"), "SPO11-1 peak (bp)")))
  abline(v = c(flankSize, length(SPO11oligosDat[,1])-flankSize), lty = 3)
  box(lwd = 1.5)

  plot(xplot, SPO11oligosRanDat[,1], ylim = c(min(SPO11oligosDat[,1], SPO11oligosRanDat[,1]), max(SPO11oligosDat[,1], SPO11oligosRanDat[,1])), type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 2, at = pretty(c(SPO11oligosDat[,1], SPO11oligosRanDat[,1])))
  #mtext(side = 2, line = 2, cex = 0.8, text = "SPO11-1", col = mycols[1])
  #mtext(side = 3, line = 1, text = paste0("Euchromatic random loci (n = ", print(numberRanLoc), ")"))
  par(new = T)
  plot(xplot, otherRanDat[,1], ylim = c(min(otherDat[,1], otherRanDat[,1]), max(otherDat[,1], otherRanDat[,1])), type = "l", lwd = 1.5, col = mycols[2], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(otherDat[,1], otherRanDat[,1])))
  mtext(side = 4, line = 2, cex = 0.8, text = otherYlabel, col = mycols[2])
  axis(side = 1, at = c(0, flankSize, length(SPO11oligosRanDat[,1])-flankSize, length(SPO11oligosRanDat[,1])), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(0, flankSize, length(SPO11oligosRanDat[,1])-flankSize, length(SPO11oligosRanDat[,1])), text = c("-1 kb", "Start", "End", "+1 kb"))
  #mtext(side = 1, line = c(3.5), cex = 0.8, text = expression(atop(paste("Position relative to"), "random locus (bp)")))
  abline(v = c(flankSize, length(SPO11oligosDat[,1])-flankSize), lty = 3)
  box(lwd = 1.5)
}

# Include base composition plots (A, T, G, C separately and (A+T)/2, (G+C)/2)
baseCompDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/base_composition/"

genome.coords <- list(read.table(file = paste0(baseCompDir, "peaks.a.txt"))[,1],
                   read.table(file = paste0(baseCompDir, "peaks.t.txt"))[,1],
                   read.table(file = paste0(baseCompDir, "peaks.g.txt"))[,1],
                   read.table(file = paste0(baseCompDir, "peaks.c.txt"))[,1])
genome.ran.coords <- list(read.table(file = paste0(baseCompDir, "ranLoc.a.txt"))[,1],
                       read.table(file = paste0(baseCompDir, "ranLoc.t.txt"))[,1],
                       read.table(file = paste0(baseCompDir, "ranLoc.g.txt"))[,1],
                       read.table(file = paste0(baseCompDir, "ranLoc.c.txt"))[,1])

baseCompPlot <- function(xplot, coords, ran.coords, ylim, mycols) {
  # peaks
  plot(xplot, coords[[1]], col = mycols[1], lwd = 1.5, type = "l", ylim = ylim, ann = F, xaxt = "n")
  mtext(side = 2, line = 2, cex = 0.8, text = "Base relative frequency")
  lines(xplot, coords[[2]], col = mycols[2], lwd = 1.5)
  lines(xplot, coords[[3]], col = mycols[3], lwd = 1.5)
  lines(xplot, coords[[4]], col = mycols[4], lwd = 1.5)
  axis(side = 1, at = c(-1000, 0, 1000), labels = c("-1 kb", "Midpoint", "+1 kb"))
  axis(side = 2, at = pretty(c(coords[[1]], coords[[2]], coords[[3]], coords[[4]], ran.coords[[1]], ran.coords[[2]], ran.coords[[3]], ran.coords[[4]])))
  abline(v = 0, lty = 3)
  #abline(v = c(-100, 99), lty = 3)
  legend("right",
         legend = c("A", "T", "G", "C"),
         col = mycols,
         text.col = mycols,
         ncol = 1, cex = 0.8, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
  # ranLoc
  plot(xplot, ran.coords[[1]], col = mycols[1], lwd = 1.5, type = "l", ylim = ylim, ann = F, xaxt = "n")
  lines(xplot, ran.coords[[2]], col = mycols[2], lwd = 1.5)
  lines(xplot, ran.coords[[3]], col = mycols[3], lwd = 1.5)
  lines(xplot, ran.coords[[4]], col = mycols[4], lwd = 1.5)
  axis(side = 1, at = c(-1000, 0, 1000), labels = c("-1 kb", "Midpoint", "+1 kb"))
  axis(side = 2, at = pretty(c(coords[[1]], coords[[2]], coords[[3]], coords[[4]], ran.coords[[1]], ran.coords[[2]], ran.coords[[3]], ran.coords[[4]])))
  abline(v = 0, lty = 3)
  #abline(v = c(-100, 99), lty = 3)
  legend("right",
         legend = c("A", "T", "G", "C"),
         col = mycols,
         text.col = mycols,
         ncol = 1, cex = 0.8, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
}

# Plot (A+T) and (G+C) around SPO11-1-oligos peaks
genome.at.coords <- (genome.coords[[1]]+genome.coords[[2]])
genome.gc.coords <- (genome.coords[[3]]+genome.coords[[4]])
genome.at.ran.coords <- (genome.ran.coords[[1]]+genome.ran.coords[[2]])
genome.gc.ran.coords <- (genome.ran.coords[[3]]+genome.ran.coords[[4]])

mergeBaseCompPlot <- function(xplot, at.coords, gc.coords, at.ran.coords, gc.ran.coords, ylim, mycols) {
  # peaks
  plot(xplot, at.coords, col = mycols[2], lwd = 1.5, type = "l", ylim = ylim, ann = F, xaxt = "n")
  mtext(side = 2, line = 2, cex = 0.8, text = "Base relative frequency")
  lines(xplot, gc.coords, col = mycols[1], lwd = 1.5)
  axis(side = 1, at = c(-1000, 0, 1000), labels = c("-1 kb", "Midpoint", "+1 kb"))
  axis(side = 2, at = pretty(c(at.coords, gc.coords, at.ran.coords, gc.ran.coords)))
  abline(v = 0, lty = 3)
  #abline(v = c(-100, 99), lty = 3)
  legend("right",
         legend = c("AT", "GC"),
         col = c(mycols[2], mycols[1]),
         text.col = c(mycols[2], mycols[1]),
         ncol = 1, cex = 0.8, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
  # ranLoc
  plot(xplot, at.ran.coords, col = mycols[2], lwd = 1.5, type = "l", ylim = ylim, ann = F, xaxt = "n")
  lines(xplot, gc.ran.coords, col = mycols[1], lwd = 1.5)
  axis(side = 1, at = c(-1000, 0, 1000), labels = c("-1 kb", "Midpoint", "+1 kb"))
  axis(side = 2, at = pretty(c(at.coords, gc.coords, at.ran.coords, gc.ran.coords)))
  abline(v = 0, lty = 3)
  #abline(v = c(-100, 99), lty = 3)
  legend("right",
         legend = c("AT", "GC"),
         col = c(mycols[2], mycols[1]),
         text.col = c(mycols[2], mycols[1]),
         ncol = 1, cex = 0.8, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
}


pdf(paste0(plotDir, "SPO11_vs_nucs_H3K4me3_AT_GC_at_SPO11_RPI1_RPI8_idr0.05_genomerangerPeaks_NTM_regioneR_v100118.pdf"), height = 7.5, width = 6)
par(mfrow = c(3, 2))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))
xplot <- seq(1, length(SPO11oligosDat[,1]), by = 1)
xplot2 <- seq(-1000, 999, by = 1)
mycols <- c("red", "blue")
mycols2 <- c("deepskyblue", "midnightblue", "red", "tomato4")
genome.ylim2 <- c(min(genome.at.coords, genome.gc.coords, genome.at.ran.coords, genome.gc.ran.coords),
               max(genome.at.coords, genome.gc.coords, genome.at.ran.coords, genome.gc.ran.coords))

SPO11oligosvsOthersPlot(xplot, "NUC", "Nucleosomes")
SPO11oligosvsOthersPlot(xplot, "H3K4me3", "H3K4me3")
mergeBaseCompPlot(xplot = xplot2, at.coords = genome.at.coords, gc.coords = genome.gc.coords, at.ran.coords = genome.at.ran.coords, gc.ran.coords = genome.gc.ran.coords, ylim = genome.ylim2, mycols = mycols)

dev.off()

sessionInfo()
