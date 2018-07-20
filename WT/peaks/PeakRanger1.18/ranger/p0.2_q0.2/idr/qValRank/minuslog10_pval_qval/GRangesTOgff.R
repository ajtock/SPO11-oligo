# Convert GRanges objects to gff files

library(GenomicRanges)

inDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/"

load(paste0(inDir, "armrangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData"))
load(paste0(inDir, "perirangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData"))

armPeaks <- data.frame(armrangerPeaksGRmerge)
periPeaks <- data.frame(perirangerPeaksGRmerge)

armPeaksgff <- cbind(armPeaks[,1], rep(".", length(armPeaks[,1])), rep("SPO11_1_oligo_hotspot"), armPeaks[,2], armPeaks[,3], rep(".", length(armPeaks[,1])), rep(".", length(armPeaks[,1])), rep(".", length(armPeaks[,1])), rep(".", length(armPeaks[,1])))
colnames(armPeaksgff) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
write.table(armPeaksgff, file = paste0(inDir, "SPO11_1_oligo_RPI1_RPI8_idr0.05_ranger_armPeaks.gff"), row.names = F, col.names = F, quote = F, sep = "\t")

periPeaksgff <- cbind(periPeaks[,1], rep(".", length(periPeaks[,1])), rep("SPO11_1_oligo_hotspot"), periPeaks[,2], periPeaks[,3], rep(".", length(periPeaks[,1])), rep(".", length(periPeaks[,1])), rep(".", length(periPeaks[,1])), rep(".", length(periPeaks[,1])))
colnames(periPeaksgff) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
write.table(periPeaksgff, file = paste0(inDir, "SPO11_1_oligo_RPI1_RPI8_idr0.05_ranger_periPeaks.gff"), row.names = F, col.names = F, quote = F, sep = "\t")

