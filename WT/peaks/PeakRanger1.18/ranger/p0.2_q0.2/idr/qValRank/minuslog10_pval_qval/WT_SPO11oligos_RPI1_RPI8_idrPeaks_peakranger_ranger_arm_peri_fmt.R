# Separate arm and pericentromeric IDR-validated peaks
# Convert to GRanges objects (with overlapping peaks unmerged or merged)

library(GenomicRanges)
library(dplyr)

# Definitions
peakDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/"
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

rangerPeaks <- read.table(file = paste0(peakDir, "WT_SPO11oligos_RPI1_RPI8_idrValues_qValRank_idr0.05.narrowPeak"))
rangerPeaks <- rangerPeaks[,1:3]
colnames(rangerPeaks) <- c("chr", "start0based", "end")
rangerPeaks <- as.data.frame(cbind(rangerPeaks$chr, rangerPeaks$start0based+1, rangerPeaks$end))
colnames(rangerPeaks) <- c("chr", "start", "end")

rangerPeaks.arms <- NULL
rangerPeaks.peri <- NULL
for(i in 1:5) {
  chr.rangerPeaks <- rangerPeaks[rangerPeaks$chr == i,]
  chr.rangerPeaks.arms <- chr.rangerPeaks %>%
    filter(end < pericenStart[i] | start > pericenEnd[i])
  chr.rangerPeaks.peri <- chr.rangerPeaks %>%
    filter(end > pericenStart[i] & start < pericenEnd[i])
  rangerPeaks.arms <- rbind(rangerPeaks.arms, chr.rangerPeaks.arms)
  rangerPeaks.peri <- rbind(rangerPeaks.peri, chr.rangerPeaks.peri)
}
print("Unfiltered arm peaks:")
print(dim(rangerPeaks.arms)[[1]])
print("Unfiltered pericentromeric and centromeric peaks:")
print(dim(rangerPeaks.peri)[[1]])
print("Unfiltered peaks:")
print(dim(rangerPeaks.arms)[[1]] + dim(rangerPeaks.peri)[[1]])

write.table(rangerPeaks.arms, file = paste0(peakDir, "rangerPeaks.arms.RPI1_RPI8_idr0.05.txt"))
write.table(rangerPeaks.peri, file = paste0(peakDir, "rangerPeaks.peri.RPI1_RPI8_idr0.05.txt"))

### Create GRanges objects
armrangerPeaksGR <- sort(GRanges(seqnames = paste0("Chr", rangerPeaks.arms$chr), ranges = IRanges(start = rangerPeaks.arms$start, end = rangerPeaks.arms$end), strand = "+"))
save(armrangerPeaksGR, file = paste0(peakDir, "armrangerPeaksGR_RPI1_RPI8_idr0.05_noMinWidth.RData"))
armrangerPeaksGRmerge <- reduce(armrangerPeaksGR)
save(armrangerPeaksGRmerge, file = paste0(peakDir, "armrangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData"))
# remove peaks < 50 bp
armrangerPeaksGR <- armrangerPeaksGR[width(armrangerPeaksGR) >= 50]
save(armrangerPeaksGR, file = paste0(peakDir, "armrangerPeaksGR_RPI1_RPI8_idr0.05.RData"))

perirangerPeaksGR <- sort(GRanges(seqnames = paste0("Chr", rangerPeaks.peri$chr), ranges = IRanges(start = rangerPeaks.peri$start, end = rangerPeaks.peri$end), strand = "+"))
save(perirangerPeaksGR, file = paste0(peakDir, "perirangerPeaksGR_RPI1_RPI8_idr0.05_noMinWidth.RData"))
perirangerPeaksGRmerge <- reduce(perirangerPeaksGR)
save(perirangerPeaksGRmerge, file = paste0(peakDir, "perirangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData"))
## remove peaks < 50 bp
perirangerPeaksGR <- perirangerPeaksGR[width(perirangerPeaksGR) >= 50]
save(perirangerPeaksGR, file = paste0(peakDir, "perirangerPeaksGR_RPI1_RPI8_idr0.05.RData"))

### Merge overlapping peaks
armrangerPeaksGRmerge <- reduce(armrangerPeaksGR)
perirangerPeaksGRmerge <- reduce(perirangerPeaksGR)
#armrangerPeaksGRmerge <- armrangerPeaksGRmerge[start(armrangerPeaksGRmerge) > 0]
save(armrangerPeaksGRmerge, file = paste0(peakDir, "armrangerPeaksGRmerge_RPI1_RPI8_idr0.05.RData"))
save(perirangerPeaksGRmerge, file = paste0(peakDir, "perirangerPeaksGRmerge_RPI1_RPI8_idr0.05.RData"))

sessionInfo()

