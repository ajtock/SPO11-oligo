# Use permutation test function in regioneR to determine if peaks overlap with features of interest (e.g., genes, TEs) more than expected
library(regioneR)

chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
#seqlevels(genome) <- sub("Chr", "", seqlevels(genome))

inDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/"
outDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/regioneR/"
plotDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/regioneR/plots/"

load(paste0(inDir, "armrangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData"))
load(paste0(inDir, "perirangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData"))
peaksGR <- sort(c(armrangerPeaksGRmerge, perirangerPeaksGRmerge))
armrangerPeaksGRmerge <- NULL
perirangerPeaksGRmerge <- NULL
strand(peaksGR) <- "*"
#seqlevels(peaksGR) <- sub("Chr", "", seqlevels(peaksGR))
#peaksGR <- peaksGR[width(peaksGR) <= 5000]
ranLocGR <- randomizeRegions(peaksGR, genome = genome, per.chromosome = TRUE, allow.overlaps = TRUE)

print(length(peaksGR))
#[1] 5914


load("/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/armPeaksSH99GRmerge.RData")
load("/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/periPeaksSH99GRmerge.RData")
nucleRnucsGR <- sort(c(armPeaksSH99GRmerge, periPeaksSH99GRmerge))
armPeaksSH99GRmerge <- NULL
periPeaksSH99GRmerge <- NULL
print(length(nucleRnucsGR))
#[1] 57734

load("/projects/ajt200/BAM_masters/nucleosomes/WT/peaks/PeakRanger1.18/ranger/nakedDNA_untrimmed_input_p0.05_q0.05_l147/armrangerPeaksGRmerge_WT_nuc_p0.05_q0.05_noMinWidth.RData")
load("/projects/ajt200/BAM_masters/nucleosomes/WT/peaks/PeakRanger1.18/ranger/nakedDNA_untrimmed_input_p0.05_q0.05_l147/perirangerPeaksGRmerge_WT_nuc_p0.05_q0.05_noMinWidth.RData")
rangernucsGR <- sort(c(armrangerPeaksGRmerge, perirangerPeaksGRmerge))
armrangerPeaksGRmerge <- NULL
perirangerPeaksGRmerge <- NULL
rangernucsGR <- rangernucsGR[width(rangernucsGR) <= 500]
print(length(rangernucsGR))
#[1] 55008

load("/projects/ajt200/GBS_CO/HS_CU_080617/wt/COsGRcoords.RData")
COsGR <- COsGRcoords
print(length(COsGR))
#[1] 3320

genes <- read.table(file = "/projects/ajt200/TAIR10/representative_genes/representative_genes_uniq_fmt_strand.txt", header = T)
genesGR <- GRanges(seqnames = genes$chr, ranges = IRanges(start = genes$start, end = genes$end), strand = "*")
seqlevels(genesGR) <- sub("", "Chr", seqlevels(genesGR))
print(length(genesGR))
#[1] 27204

genesGRprom <- GRanges(seqnames = genes$chr, ranges = IRanges(start = genes$start, end = genes$end), strand = genes$strand)
print(length(genesGRprom))
#[1] 27204
promotersGR <- promoters(genesGRprom, upstream = 500, downstream = 0)
seqlevels(promotersGR) <- sub("", "Chr", seqlevels(promotersGR))
strand(promotersGR) <- "*"
print(length(promotersGR))
#[1] 27204

source("/projects/ajt200/Rfunctions/downstream_TTS.r")
genesGRterm <- GRanges(seqnames = genes$chr, ranges = IRanges(start = genes$start, end = genes$end), strand = genes$strand)
print(length(genesGRterm))
#[1] 27204
terminatorsGR <- ttsPlus(genesGRterm, upstream = -1, downstream = 500)
seqlevels(terminatorsGR) <- sub("", "Chr", seqlevels(terminatorsGR))
strand(terminatorsGR) <- "*"
print(length(terminatorsGR))
#[1] 27204

TEs <- read.table(file = "/projects/ajt200/TAIR10/TAIR10_Buisine_TEs_strand_tab_ann.txt", header=T)
TEsGR <- GRanges(seqnames = TEs$Chr, ranges = IRanges(start = TEs$start, end = TEs$end), strand = "*")
print(length(TEsGR))
#[1] 31189

otherNames <- c("COs", "nucleRnucs", "rangernucs",
                 "genes", "promoters", "terminators", "TEs")

grl <- GRangesList("COsGR" = COsGR, "nucleRnucsGR" = nucleRnucsGR, "rangernucsGR" = rangernucsGR,
                   "genesGR" = genesGR, "promotersGR" = promotersGR, "terminatorsGR" = terminatorsGR, "TEsGR" = TEsGR)

numOverlaps(A = peaksGR, B = grl[[1]], count.once=TRUE)
#[1] 388
numOverlaps(A = peaksGR, B = grl[[1]], count.once=FALSE)
#[1] 418

# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions in B as in A
set.seed(123)
ptPeaksOtherPerChrom <- lapply(seq_along(grl), function(x) {
  permTest(A = peaksGR, B = grl[[x]], genome = genome,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE, per.chromosome = TRUE,
           evaluate.function = numOverlaps, count.once = TRUE,
           ntimes = 10000, mc.set.seed = FALSE, mc.cores = 47)
})
 
for(i in 1:length(ptPeaksOtherPerChrom)) {
  assign(paste0(otherNames[i]), ptPeaksOtherPerChrom[[i]])
}
save(ptPeaksOtherPerChrom, file = paste0(outDir, "pt_SPO11_RPI1_RPI8_idr0.05_peaks_COs_nucleRnucs_rangernucs_genes_promoters_terminators_TEs_noMinWidth.RData"))

for(i in 1:length(ptPeaksOtherPerChrom)) {
  pdf(file = paste0(plotDir, otherNames[i], "_permTest_nperm10000_SPO11_RPI1_RPI8_idr0.05_genomerangerPeaks_perChrom.pdf"), width = 10, height = 7) 
  plot(ptPeaksOtherPerChrom[[i]], main = paste0("SPO11_RPI1_RPI8_idr0.05_genomerangerPeaks vs ", otherNames[i], " genome wide"), xlab = "Number of overlaps", ylab = "Proportion")
  dev.off()

  # Using the localZScore() function, evaluate whether the association between peaks and other is highly dependent on their exact position
  lz_1kb <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = grl[[i]],
                        window = 1000, step = 50, count.once = TRUE)
  lz_10kb <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = grl[[i]],
                         window = 10000, step = 500, count.once = TRUE)
  lz_custom <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = grl[[i]],
                           window = 10*mean(width(peaksGR)), step = mean(width(peaksGR))/2, count.once = TRUE)

  pdf(file = paste0(plotDir, otherNames[i], "_localZscore_permTest_nperm10000_SPO11_RPI1_RPI8_idr0.05_genomerangerPeaks_w1kb_s50bp_w10kb_s500bp_w8kb_s412bp_perChrom.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(lz_1kb, main = paste0("SPO11_RPI1_RPI8_idr0.05_genomerangerPeaks vs ", otherNames[i], " genome wide (1-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("SPO11_RPI1_RPI8_idr0.05_genomerangerPeaks vs ", otherNames[i], " genome wide (1-kb shift)"))
  plot(lz_10kb, main = paste0("SPO11_RPI1_RPI8_idr0.05_genomerangerPeaks vs ", otherNames[i], " genome wide (10-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("SPO11_RPI1_RPI8_idr0.05_genomerangerPeaks vs ", otherNames[i], " genome wide (10-kb shift)"))
  plot(lz_custom, main = paste0("SPO11_RPI1_RPI8_idr0.05_genomerangerPeaks vs ", otherNames[i], " genome wide (~8-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("SPO11_RPI1_RPI8_idr0.05_genomerangerPeaks vs ", otherNames[i], " genome wide (~8-kb shift)"))
  dev.off()
}

## Use below to access the various parameters of the ptPeaksOtherPerChrom object (mean "permuted" (expected) overlaps, "observed" overlaps, "pval", "zscore")
#load(file = paste0(outDir, "pt_SPO11_RPI1_RPI8_idr0.05_peaks_COs_nucleRnucs_rangernucs_genes_promoters_terminators_TEs_noMinWidth.RData"))
#for(i in 1:length(ptPeaksOtherPerChrom)) {
#  print(mean(ptPeaksOtherPerChrom[[i]]$numOverlaps$permuted))
#}


