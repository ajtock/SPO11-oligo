#!/applications/R/R-3.3.2/bin/Rscript

# Use permutation test function in regioneR to determine if peaks overlap features
# of interest (e.g., other peaks, genes, TEs) more or less than expected by chance

# Usage on hydrogen node7:
# csmit -m 24G -c 24 "Rscript permTest_SPO11_RPI1_RPI8_rangerPeaks_vs_others.R"

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

# REC8 peaks
load("/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1_armrangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData")
load("/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep1_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData")
REC8_HA_Rep1GR <- sort(c(armrangerPeaksGRmergedOverlaps, perirangerPeaksGRmergedOverlaps))
armrangerPeaksGRmergedOverlaps <- NULL
perirangerPeaksGRmergedOverlaps <- NULL
strand(REC8_HA_Rep1GR) <- "*"
print("***********peaks***********")
print(REC8_HA_Rep1GR)
print(length(REC8_HA_Rep1GR))
#[1] 87738

load("/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep2_armrangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData")
load("/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_HA_Rep2_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData")
REC8_HA_Rep2GR <- sort(c(armrangerPeaksGRmergedOverlaps, perirangerPeaksGRmergedOverlaps))
armrangerPeaksGRmergedOverlaps <- NULL
perirangerPeaksGRmergedOverlaps <- NULL
strand(REC8_HA_Rep2GR) <- "*"
print("***********peaks***********")
print(REC8_HA_Rep2GR)
print(length(REC8_HA_Rep2GR))
#[1] 82900

load("/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_MYC_Rep1_armrangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData")
load("/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/REC8_MYC_Rep1_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData")
REC8_MYC_Rep1GR <- sort(c(armrangerPeaksGRmergedOverlaps, perirangerPeaksGRmergedOverlaps))
armrangerPeaksGRmergedOverlaps <- NULL
perirangerPeaksGRmergedOverlaps <- NULL
strand(REC8_MYC_Rep1GR) <- "*"
print("***********peaks***********")
print(REC8_MYC_Rep1GR)
print(length(REC8_MYC_Rep1GR))
#[1] 60883

load("/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/kss_REC8_HA_Rep1_armrangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData")
load("/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/kss_REC8_HA_Rep1_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData")
kss_REC8_HA_Rep1GR <- sort(c(armrangerPeaksGRmergedOverlaps, perirangerPeaksGRmergedOverlaps))
armrangerPeaksGRmergedOverlaps <- NULL
perirangerPeaksGRmergedOverlaps <- NULL
strand(kss_REC8_HA_Rep1GR) <- "*"
print("***********peaks***********")
print(kss_REC8_HA_Rep1GR)
print(length(kss_REC8_HA_Rep1GR))
#[1] 103616

# Others
load("/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/armPeaksSH99GRmerge.RData")
load("/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/periPeaksSH99GRmerge.RData")
nucleRnucsGR <- sort(c(armPeaksSH99GRmerge, periPeaksSH99GRmerge))
armPeaksSH99GRmerge <- NULL
periPeaksSH99GRmerge <- NULL
strand(nucleRnucsGR) <- "*"
print(length(nucleRnucsGR))
#[1] 57734

load("/projects/ajt200/BAM_masters/nucleosomes/WT/peaks/PeakRanger1.18/ranger/nakedDNA_untrimmed_input_p0.05_q0.05_l147/armrangerPeaksGRmerge_WT_nuc_p0.05_q0.05_noMinWidth.RData")
load("/projects/ajt200/BAM_masters/nucleosomes/WT/peaks/PeakRanger1.18/ranger/nakedDNA_untrimmed_input_p0.05_q0.05_l147/perirangerPeaksGRmerge_WT_nuc_p0.05_q0.05_noMinWidth.RData")
rangernucsGR <- sort(c(armrangerPeaksGRmerge, perirangerPeaksGRmerge))
armrangerPeaksGRmerge <- NULL
perirangerPeaksGRmerge <- NULL
strand(rangernucsGR) <- "*"
rangernucsGR <- rangernucsGR[width(rangernucsGR) <= 500]
print(length(rangernucsGR))
#[1] 55008

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

load("/projects/ajt200/BAM_masters/SPO11_ChIP/Xiaohui_BAM_and_coverage_files/WT/peaks/PeakRanger1.18/ranger/REC8_MYC_Rep2_input_p0.05_q0.05/armrangerPeaksGRmerge_WT_SPO11_ChIP4_p0.05_q0.05_noMinWidth.RData")
load("/projects/ajt200/BAM_masters/SPO11_ChIP/Xiaohui_BAM_and_coverage_files/WT/peaks/PeakRanger1.18/ranger/REC8_MYC_Rep2_input_p0.05_q0.05/perirangerPeaksGRmerge_WT_SPO11_ChIP4_p0.05_q0.05_noMinWidth.RData")
SPO11_ChIP4GR <- sort(c(armrangerPeaksGRmerge, perirangerPeaksGRmerge))
armrangerPeaksGRmerge <- NULL
perirangerPeaksGRmerge <- NULL
strand(SPO11_ChIP4GR) <- "*"
print(length(SPO11_ChIP4GR))
#[1] 25825

load("/projects/ajt200/BAM_masters/SPO11_ChIP/Xiaohui_BAM_and_coverage_files/WT/peaks/PeakRanger1.18/ranger/REC8_MYC_Rep2_input_p0.05_q0.05/armrangerPeaksGRmerge_WT_SPO11_ChIP13_p0.05_q0.05_noMinWidth.RData")
load("/projects/ajt200/BAM_masters/SPO11_ChIP/Xiaohui_BAM_and_coverage_files/WT/peaks/PeakRanger1.18/ranger/REC8_MYC_Rep2_input_p0.05_q0.05/perirangerPeaksGRmerge_WT_SPO11_ChIP13_p0.05_q0.05_noMinWidth.RData")
SPO11_ChIP13GR <- sort(c(armrangerPeaksGRmerge, perirangerPeaksGRmerge))
armrangerPeaksGRmerge <- NULL
perirangerPeaksGRmerge <- NULL
strand(SPO11_ChIP13GR) <- "*"
print(length(SPO11_ChIP13GR))
#[1] 20370

load("/projects/ajt200/BAM_masters/H3K4me3/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/armrangerPeaksGRmerge_WT_H3K4me3_ChIP14_WT_H3K4me3_ChIP15_idr0.05_noMinWidth.RData")
load("/projects/ajt200/BAM_masters/H3K4me3/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/perirangerPeaksGRmerge_WT_H3K4me3_ChIP14_WT_H3K4me3_ChIP15_idr0.05_noMinWidth.RData")
H3K4me3GR <- sort(c(armrangerPeaksGRmerge, perirangerPeaksGRmerge))
armrangerPeaksGRmerge <- NULL
perirangerPeaksGRmerge <- NULL
strand(H3K4me3GR) <- "*"
print(length(H3K4me3GR))
#[1] 13951

load("/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/ranger/p0.05_q0.05/WT_H3K9me2_ChIP_armrangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth.RData")
load("/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/ranger/p0.05_q0.05/WT_H3K9me2_ChIP_perirangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth.RData")
H3K9me2GR <- sort(c(armrangerPeaksGRmergedOverlaps, perirangerPeaksGRmergedOverlaps))
armrangerPeaksGRmergedOverlaps <- NULL
perirangerPeaksGRmergedOverlaps <- NULL 
strand(H3K9me2GR) <- "*"
print(length(H3K9me2GR))
#[1] 20289

load("/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/BCP/p0.05/WT_H3K9me2_ChIP_armbcpPeaksGRmergedOverlaps_p0.05_noMinWidth.RData")
load("/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/BCP/p0.05/WT_H3K9me2_ChIP_peribcpPeaksGRmergedOverlaps_p0.05_noMinWidth.RData")
H3K9me2GRbcp <- sort(c(armbcpPeaksGRmergedOverlaps, peribcpPeaksGRmergedOverlaps))
armbcpPeaksGRmergedOverlaps <- NULL
peribcpPeaksGRmergedOverlaps <- NULL
strand(H3K9me2GRbcp) <- "*"
print(length(H3K9me2GRbcp))
#[1] 479

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
promotersGR <- promoters(genesGRprom, upstream = 500, downstream = 0)
seqlevels(promotersGR) <- sub("", "Chr", seqlevels(promotersGR))
strand(promotersGR) <- "*"
print(length(promotersGR))
#[1] 27204

TSSdownstream500GR <- promoters(genesGRprom, upstream = 0, downstream = 500)
seqlevels(TSSdownstream500GR) <- sub("", "Chr", seqlevels(TSSdownstream500GR))
strand(TSSdownstream500GR) <- "*"
print(length(TSSdownstream500GR))
#[1] 27204

source("/projects/ajt200/Rfunctions/TTSplus.R")
genesGRterm <- GRanges(seqnames = genes$chr, ranges = IRanges(start = genes$start, end = genes$end), strand = genes$strand)
terminatorsGR <- TTSplus(genesGRterm, upstream = -1, downstream = 500)
seqlevels(terminatorsGR) <- sub("", "Chr", seqlevels(terminatorsGR))
strand(terminatorsGR) <- "*"
print(length(terminatorsGR))
#[1] 27204

TTSupstream500GR <- TTSplus(genesGRterm, upstream = 499, downstream = 0)
seqlevels(TTSupstream500GR) <- sub("", "Chr", seqlevels(TTSupstream500GR))
strand(TTSupstream500GR) <- "*"
print(length(TTSupstream500GR))
#[1] 27204

# Import exons as GRanges object
genes_exons <- read.table("/projects/ajt200/TAIR10/all_exons.txt",
                          header = T)
exons <- genes_exons[genes_exons$ge.ex == "exon",]
levels(exons$strand) <- c("-", "*", "+")
exonsGR <- GRanges(seqnames = exons$chr,
                   ranges = IRanges(start = exons$start, end = exons$end),
                   strand = "*")
print(length(exonsGR))
#[1] 152155

# Import introns tables and convert to GRanges object
intronsPlus <- read.table("/projects/ajt200/TAIR10/all_plus_introns.txt", header = T)
intronsMinus <- read.table("/projects/ajt200/TAIR10/all_minus_introns.txt", header = T)
intronsPlusGR <- GRanges(seqnames = paste0("Chr", intronsPlus$chr),
                         ranges = (IRanges(start = intronsPlus$all.intron.starts+1,
                                           end = intronsPlus$all.intron.stops-1)),
                         strand = "*")
intronsMinusGR <- GRanges(seqnames = paste0("Chr", intronsMinus$chr),
                          ranges = (IRanges(start = intronsMinus$all.intron.stops+1,
                                            end = intronsMinus$all.intron.starts-1)),
                          strand = "*")
intronsGR <- sort(append(intronsPlusGR, intronsMinusGR), by = ~ seqnames + start + end)
print(length(intronsGR))
#[1] 119808

TEs <- read.table(file = "/projects/ajt200/TAIR10/TAIR10_Buisine_TEs_strand_tab_ann.txt", header=T)
TEsGR <- GRanges(seqnames = TEs$Chr, ranges = IRanges(start = TEs$start, end = TEs$end), strand = "*")
print(length(TEsGR))
#[1] 31189

otherNames <- c("REC8_HA_Rep2GR", "REC8_MYC_Rep1GR", "kss_REC8_HA_Rep1GR",
                "nucleRnucsGR", "rangernucsGR",
                "REC8_HA_Rep1GR", "SPO11_ChIP4GR", "SPO11_ChIP13GR",
                "H3K4me3GR", "H3K9me2GR", "H3K9me2GRbcp", "COsGR",
                "genesGR", "promotersGR", "terminatorsGR",
                "TSSdownstream500GR", "TTSupstream500GR",
                "exonsGR", "intronsGR", "TEsGR")

grl <- GRangesList("REC8_HA_Rep2GR" = REC8_HA_Rep2GR, "REC8_MYC_Rep1GR" = REC8_MYC_Rep1GR, "kss_REC8_HA_Rep1GR" = kss_REC8_HA_Rep1GR,
                   "nucleRnucsGR" = nucleRnucsGR, "rangernucsGR" = rangernucsGR,
                   "REC8_HA_Rep1GR" = REC8_HA_Rep1GR, "SPO11_ChIP4GR" = SPO11_ChIP4GR, "SPO11_ChIP13GR" = SPO11_ChIP13GR,
                   "H3K4me3GR" = H3K4me3GR, "H3K9me2GR" = H3K9me2GR, "H3K9me2GRbcp" = H3K9me2GRbcp, "COsGR" = COsGR,
                   "genesGR" = genesGR, "promotersGR" = promotersGR, "terminatorsGR" = terminatorsGR,
                   "TSSdownstream500GR" = TSSdownstream500GR, "TTSupstream500GR" = TTSupstream500GR,
                   "exonsGR" = exonsGR, "intronsGR" = intronsGR, "TEsGR" = TEsGR)

# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions in B as in A
set.seed(123)
ptPeaksOtherPerChrom <- lapply(seq_along(grl), function(x) {
  permTest(A = peaksGR, B = grl[[x]], genome = genome,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE, per.chromosome = TRUE,
           evaluate.function = numOverlaps, count.once = TRUE,
           ntimes = 10000, mc.set.seed = FALSE, mc.cores = 24)
})

for(i in 1:length(ptPeaksOtherPerChrom)) {
  assign(paste0(otherNames[i]), ptPeaksOtherPerChrom[[i]])
}
save(ptPeaksOtherPerChrom,
     file = paste0(outDir,
                   "permTest_SPO11_RPI1_RPI8_rangerPeaks_vs_others.RData"))

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
                          "permTest_SPO11_RPI1_RPI8_rangerPeaks_vs_others_DataFrame.txt"),
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


