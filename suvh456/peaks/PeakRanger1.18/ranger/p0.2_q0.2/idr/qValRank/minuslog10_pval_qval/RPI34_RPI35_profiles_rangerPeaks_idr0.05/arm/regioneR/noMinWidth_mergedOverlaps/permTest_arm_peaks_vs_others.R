#!/applications/R/R-3.3.2/bin/Rscript

# Use permutation test function in regioneR to determine
# if peaks overlap features of interest more or less
# than expected by chance

# Usage on hydrogen node7:
# csmit -m 100G -c 32 "/applications/R/R-3.3.2/bin/Rscript permTest_arm_peaks_vs_others.R REC8_MYC_Rep1_ChIP_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData REC8_HA_Rep2_ChIP_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData 10000"

#peakFile1 <- "REC8_MYC_Rep1_ChIP_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData"
#peakFile2 <- "REC8_HA_Rep2_ChIP_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData"
#perms <- 10000

args <- commandArgs(trailingOnly = T)
peakFile1 <- args[1]
peakFile2 <- args[2]
perms <- as.numeric(args[3])

library(regioneR)

# Genomic definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
mask <- toGRanges(data.frame(chrs, pericenStart, pericenEnd))

inDir <- "/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep1_input_p0.001_q0.01/"
outDir <- "./"

plotDir <- "./histograms/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Test peaks
load("/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/armrangerPeaksGRmerge_RPI34_RPI35_idr0.05_noMinWidth.RData")
peaksGR <- armrangerPeaksGRmerge
armrangerPeaksGRmerge <- NULL
strand(peaksGR) <- "*"
print(length(peaksGR))

# Others
load(paste0(inDir,
            peakFile2))
REC8_HA_Rep2GR <- rangerPeaksGR_arm_mergedOverlaps
rangerPeaksGR_arm_mergedOverlaps <- NULL
strand(REC8_HA_Rep2GR) <- "*"
print("***********REC8_HA_Rep2_peaks***********")
print(REC8_HA_Rep2GR)
print(length(REC8_HA_Rep2GR))

load("/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/armPeaksSH99GRmerge.RData")
nucleRnucsGR <- armPeaksSH99GRmerge
armPeaksSH99GRmerge <- NULL
strand(nucleRnucsGR) <- "*"
print(length(nucleRnucsGR))
#[1] 37529

load(paste0(inDir,
            peakFile1))
REC8_MYC_Rep1GR <- rangerPeaksGR_arm_mergedOverlaps
rangerPeaksGR_arm_mergedOverlaps <- NULL
strand(REC8_MYC_Rep1GR) <- "*"
print("***********REC8_MYC_Rep1_peaks***********")
print(REC8_MYC_Rep1GR)
print(length(REC8_MYC_Rep1GR))

load("/projects/ajt200/BAM_masters/SPO11_ChIP/WT/peaks/PeakRanger1.18/ranger/MYC_Rep1_input_p0.001_q0.01/WT_SPO11_ChIP4_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData")
SPO11_ChIP4GR <- rangerPeaksGR_arm_mergedOverlaps
rangerPeaksGR_arm_mergedOverlaps <- NULL
strand(SPO11_ChIP4GR) <- "*"
print(length(SPO11_ChIP4GR))
#[1] 3093

load("/projects/ajt200/BAM_masters/SPO11_ChIP/WT/peaks/PeakRanger1.18/ranger/MYC_Rep1_input_p0.001_q0.01/WT_SPO11_ChIP13_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData")
SPO11_ChIP13GR <- rangerPeaksGR_arm_mergedOverlaps
rangerPeaksGR_arm_mergedOverlaps <- NULL
strand(SPO11_ChIP13GR) <- "*"
print(length(SPO11_ChIP13GR))
#[1] 1189 

load("/projects/ajt200/BAM_masters/H3K4me3/WT/peaks/PeakRanger1.18/ranger/H3K9me2_input_p0.05_q0.05/WT_H3K4me3_ChIP14_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData")
H3K4me3_ChIP14GR <- rangerPeaksGR_arm_mergedOverlaps
rangerPeaksGR_arm_mergedOverlaps <- NULL
strand(H3K4me3_ChIP14GR) <- "*"
print(length(H3K4me3_ChIP14GR))
#[1] 17866 

load("/projects/ajt200/BAM_masters/H3K4me3/WT/peaks/PeakRanger1.18/ranger/H3K9me2_input_p0.05_q0.05/WT_H3K4me3_ChIP15_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData")
H3K4me3_ChIP15GR <- rangerPeaksGR_arm_mergedOverlaps
rangerPeaksGR_arm_mergedOverlaps <- NULL
strand(H3K4me3_ChIP15GR) <- "*"
print(length(H3K4me3_ChIP15GR))
#[1] 17101

load("/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/ranger/REC8_MYC_Rep1_input_p0.05_q0.05/WT_H3K9me2_ChIP_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData")
H3K9me2GR <- rangerPeaksGR_arm_mergedOverlaps
rangerPeaksGR_arm_mergedOverlaps <- NULL
strand(H3K9me2GR) <- "*"
print(length(H3K9me2GR))
#[1] 327

load("/home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me1/peaks/PeakRanger1.18/ranger/H3K9me2_input_p0.05_q0.05/WT_H3K4me1_ChIP_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData")
H3K4me1GR <- rangerPeaksGR_arm_mergedOverlaps
rangerPeaksGR_arm_mergedOverlaps <- NULL
strand(H3K4me1GR) <- "*"
print(length(H3K4me1GR))
#[1] 9618 

load("/home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K4me2/peaks/PeakRanger1.18/ranger/H3K9me2_input_p0.05_q0.05/WT_H3K4me2_ChIP_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData")
H3K4me2GR <- rangerPeaksGR_arm_mergedOverlaps
rangerPeaksGR_arm_mergedOverlaps <- NULL
strand(H3K4me2GR) <- "*"
print(length(H3K4me2GR))
#[1] 14230

load("/home/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/H3K27me1/peaks/PeakRanger1.18/ranger/H3K9me2_input_p0.001_q0.01/WT_H3K27me1_ChIP_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData")
H3K27me1GR <- rangerPeaksGR_arm_mergedOverlaps
rangerPeaksGR_arm_mergedOverlaps <- NULL
strand(H3K27me1GR) <- "*"
print(length(H3K27me1GR))
#[1] 13090 

load("/home/ajt200/analysis/H3K27me3_bud_UWMadison_2015/peaks/PeakRanger1.18/ranger/H3K9me2_input_p0.05_q0.05/H3K27me3_ChIP_SRR1509478_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData")
H3K27me3GR <- rangerPeaksGR_arm_mergedOverlaps
rangerPeaksGR_arm_mergedOverlaps <- NULL
strand(H3K27me3GR) <- "*"
print(length(H3K27me3GR))
#[1] 14920

load("/projects/ajt200/GBS_CO/HS_CU_080617/wt/COsGRcoords.RData")
COsGR <- COsGRcoords
print(length(COsGR))
# Remove COs located within pericentromeric regions
maskCOsOverlaps <- findOverlaps(query = mask,
                                subject = COsGR,
                                ignore.strand = TRUE,
                                select = "all")
COsGR <- COsGR[-subjectHits(maskCOsOverlaps)]
print(length(COsGR))

genes <- read.table("/projects/ajt200/TAIR10/representative_genes/representative_genes_uniq_fmt_strand.txt",
                    header = T)
genesGR <- GRanges(seqnames = genes$chr,
                   ranges = IRanges(start = genes$start,
                                    end = genes$end),
                                    strand = "*")
seqlevels(genesGR) <- sub("", "Chr", seqlevels(genesGR))
print(length(genesGR))
#[1] 27204
# Remove genes located within pericentromeric regions
maskgenesOverlaps <- findOverlaps(query = mask,
                                  subject = genesGR,
                                  ignore.strand = TRUE,
                                  select = "all")
genesGR <- genesGR[-subjectHits(maskgenesOverlaps)]
print(length(genesGR))
#[1] 23473

genesGRprom <- GRanges(seqnames = genes$chr,
                       ranges = IRanges(start = genes$start,
                                        end = genes$end),
                       strand = genes$strand)
promotersGR <- promoters(genesGRprom, upstream = 500, downstream = 0)
seqlevels(promotersGR) <- sub("", "Chr", seqlevels(promotersGR))
strand(promotersGR) <- "*"
print(length(promotersGR))
#[1] 27204
# Remove promoters located within pericentromeric regions
maskpromotersOverlaps <- findOverlaps(query = mask,
                                      subject = promotersGR,
                                      ignore.strand = TRUE,
                                      select = "all")
promotersGR <- promotersGR[-subjectHits(maskpromotersOverlaps)]
print(length(promotersGR))
#[1] 23475

TSSdownstream500GR <- promoters(genesGRprom, upstream = 0, downstream = 500)
seqlevels(TSSdownstream500GR) <- sub("", "Chr", seqlevels(TSSdownstream500GR))
strand(TSSdownstream500GR) <- "*"
print(length(TSSdownstream500GR))
#[1] 27204
# Remove TSSdownstream500 located within pericentromeric regions
maskTSSdownstreamOverlaps <- findOverlaps(query = mask,
                                          subject = TSSdownstream500GR,
                                          ignore.strand = TRUE,
                                          select = "all")
TSSdownstream500GR <- TSSdownstream500GR[-subjectHits(maskTSSdownstreamOverlaps)]
print(length(TSSdownstream500GR))
#[1] 23475

source("/projects/ajt200/Rfunctions/TTSplus.R")
genesGRterm <- GRanges(seqnames = genes$chr,
                       ranges = IRanges(start = genes$start,
                                        end = genes$end),
                       strand = genes$strand)
terminatorsGR <- TTSplus(genesGRterm, upstream = -1, downstream = 500)
seqlevels(terminatorsGR) <- sub("", "Chr", seqlevels(terminatorsGR))
strand(terminatorsGR) <- "*"
print(length(terminatorsGR))
#[1] 27204
# Remove terminators located within pericentromeric regions
maskterminatorsOverlaps <- findOverlaps(query = mask,
                                        subject = terminatorsGR,
                                        ignore.strand = TRUE,
                                        select = "all")
terminatorsGR <- terminatorsGR[-subjectHits(maskterminatorsOverlaps)]
print(length(terminatorsGR))
#[1] 23478

TTSupstream500GR <- TTSplus(genesGRterm, upstream = 499, downstream = 0)
seqlevels(TTSupstream500GR) <- sub("", "Chr", seqlevels(TTSupstream500GR))
strand(TTSupstream500GR) <- "*"
print(length(TTSupstream500GR))
#[1] 27204
# Remove TTSupstream500 located within pericentromeric regions
maskTTSupstreamOverlaps <- findOverlaps(query = mask,
                                        subject = TTSupstream500GR,
                                        ignore.strand = TRUE,
                                        select = "all")
TTSupstream500GR <- TTSupstream500GR[-subjectHits(maskTTSupstreamOverlaps)]
print(length(TTSupstream500GR))
#[1] 23478

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
# Remove exons located within pericentromeric regions
maskexonsOverlaps <- findOverlaps(query = mask,
                                  subject = exonsGR,
                                  ignore.strand = TRUE,
                                  select = "all")
exonsGR <- exonsGR[-subjectHits(maskexonsOverlaps)]
print(length(exonsGR))
#[1] 129044

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
# Remove introns located within pericentromeric regions
maskintronsOverlaps <- findOverlaps(query = mask,
                                    subject = intronsGR,
                                    ignore.strand = TRUE,
                                    select = "all")
intronsGR <- intronsGR[-subjectHits(maskintronsOverlaps)]
print(length(intronsGR))
#[1] 103930

TEs <- read.table(file = "/projects/ajt200/TAIR10/TAIR10_Buisine_TEs_strand_tab_ann.txt", header=T)
TEsGR <- GRanges(seqnames = TEs$Chr, ranges = IRanges(start = TEs$start, end = TEs$end), strand = "*")
print(length(TEsGR))
#[1] 31189
# Remove TEs located within pericentromeric regions
maskTEsOverlaps <- findOverlaps(query = mask, 
                                subject = TEsGR,
                                ignore.strand = TRUE,
                                select = "all")
TEsGR <- TEsGR[-subjectHits(maskTEsOverlaps)]
print(length(TEsGR))
#[1] 13016

kss_hypoCHG_DMRs <- read.table("/home/ajt200/BS_Seq/Stroud_2013/DMRs/suvh456_hypoCHG_DMR_vs3reps_min4filter_mg200.bed", header = F)
kss_hypoCHG_DMRsGR <- sort(GRanges(seqnames = kss_hypoCHG_DMRs[,1],
                                   ranges = IRanges(start = kss_hypoCHG_DMRs[,2],
                                                    end = kss_hypoCHG_DMRs[,3]),
                                   strand = "*"))
seqlevels(kss_hypoCHG_DMRsGR) <- sub("chr", "Chr", seqlevels(kss_hypoCHG_DMRsGR))
mask_kss_hypoCHG_DMRsOverlaps <- findOverlaps(query = mask,
                                              subject = kss_hypoCHG_DMRsGR,
                                              ignore.strand = TRUE,
                                              select = "all")
kss_hypoCHG_DMRsGR <- kss_hypoCHG_DMRsGR[-subjectHits(mask_kss_hypoCHG_DMRsOverlaps)]
print(length(kss_hypoCHG_DMRsGR))
#[1] 3136

kss_hypoCHH_DMRs <- read.table("/home/ajt200/BS_Seq/Stroud_2013/DMRs/suvh456_hypoCHH_DMR_vs3reps_min4filter_mg200.bed", header = F)
kss_hypoCHH_DMRsGR <- sort(GRanges(seqnames = kss_hypoCHH_DMRs[,1],
                                   ranges = IRanges(start = kss_hypoCHH_DMRs[,2],
                                                    end = kss_hypoCHH_DMRs[,3]),
                                   strand = "*"))
seqlevels(kss_hypoCHH_DMRsGR) <- sub("chr", "Chr", seqlevels(kss_hypoCHH_DMRsGR))
mask_kss_hypoCHH_DMRsOverlaps <- findOverlaps(query = mask,
                                              subject = kss_hypoCHH_DMRsGR,
                                              ignore.strand = TRUE,
                                              select = "all")
kss_hypoCHH_DMRsGR <- kss_hypoCHH_DMRsGR[-subjectHits(mask_kss_hypoCHH_DMRsOverlaps)]
print(length(kss_hypoCHH_DMRsGR))
#[1] 1918

otherNames <- c("REC8_HA_Rep2GR",
                "nucleRnucsGR",
                "REC8_MYC_Rep1GR",
                "SPO11_ChIP4GR",
                "SPO11_ChIP13GR",
                "H3K4me3_ChIP14GR",
                "H3K4me3_ChIP15GR",
                "H3K4me1GR",
                "H3K4me2GR",
                "H3K27me1GR",
                "H3K27me3GR",
                "H3K9me2GR",
                "COsGR",
                "genesGR",
                "promotersGR",
                "terminatorsGR",
                "TSSdownstream500GR",
                "TTSupstream500GR",
                "exonsGR",
                "intronsGR",
                "TEsGR",
                "kss_hypoCHG_DMRsGR",
                "kss_hypoCHH_DMRsGR")

grl <- c("REC8_HA_Rep2GR" = REC8_HA_Rep2GR,
         "nucleRnucsGR" = nucleRnucsGR,
         "REC8_MYC_Rep1GR" = REC8_MYC_Rep1GR,
         "SPO11_ChIP4GR" = SPO11_ChIP4GR,
         "SPO11_ChIP13GR" = SPO11_ChIP13GR,
         "H3K4me3_ChIP14GR" = H3K4me3_ChIP14GR,
         "H3K4me3_ChIP15GR" = H3K4me3_ChIP15GR,
         "H3K4me1GR" = H3K4me1GR,
         "H3K4me2GR" = H3K4me2GR,
         "H3K27me1GR" = H3K27me1GR,
         "H3K27me3GR" = H3K27me3GR,
         "H3K9me2GR" = H3K9me2GR,
         "COsGR" = COsGR,
         "genesGR" = genesGR,
         "promotersGR" = promotersGR,
         "terminatorsGR" = terminatorsGR,
         "TSSdownstream500GR" = TSSdownstream500GR,
         "TTSupstream500GR" = TTSupstream500GR,
         "exonsGR" = exonsGR,
         "intronsGR" = intronsGR,
         "TEsGR" = TEsGR,
         "kss_hypoCHG_DMRsGR" = kss_hypoCHG_DMRsGR,
         "kss_hypoCHH_DMRsGR" = kss_hypoCHH_DMRsGR)

# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions as in A
set.seed(38402)
ptPeaksOtherPerChrom <- lapply(seq_along(grl), function(x) {
  permTest(A = peaksGR,
           B = grl[[x]],
           genome = genome,
           mask = mask,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE,
           per.chromosome = TRUE,
           evaluate.function = numOverlaps,
           count.once = TRUE,
           ntimes = perms,
           mc.set.seed = FALSE,
           mc.cores = detectCores())
})

for(i in 1:length(ptPeaksOtherPerChrom)) {
  assign(paste0(otherNames[i]), ptPeaksOtherPerChrom[[i]])
}
save(ptPeaksOtherPerChrom,
     file = paste0(outDir,
                   "permTest_", as.character(perms), "perms_kss_SPO11oligo_hotspots_arm_vs_others.RData"))

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
                          "permTest_", as.character(perms), "perms_kss_SPO11oligo_hotspots_arm_vs_others_DataFrame.txt"),
            sep = "\t", row.names = F)

# plot graphical summaries of results
for(i in 1:length(ptPeaksOtherPerChrom)) {
  pdf(file = paste0(plotDir, otherNames[i], "_permTest_", as.character(perms), "perms_kss_SPO11oligo_hotspots_arm_perChrom.pdf"), width = 10, height = 7)
  plot(ptPeaksOtherPerChrom[[i]], main = paste0("kss SPO11oligo hotspots arm vs ", otherNames[i]), xlab = "Number of overlaps", ylab = "Density")
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
  pdf(file = paste0(plotDir, otherNames[i], "_localZscore_permTest_", as.character(perms), "perms_kss_SPO11oligo_hotspots_arm_w1kb_s50bp_w10kb_s500bp_w", win ,"kb_s", step, "bp_perChrom.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(lz_1kb, main = paste0("kss SPO11oligo hotspots arm vs ", otherNames[i], " (1-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("kss SPO11oligo hotspots arm vs ", otherNames[i], " (1-kb shift)"))
  plot(lz_10kb, main = paste0("kss SPO11oligo hotspots arm vs ", otherNames[i], " (10-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("kss SPO11oligo hotspots arm vs ", otherNames[i], " (10-kb shift)"))
  plot(lz_custom, main = paste0("kss SPO11oligo hotspots arm vs ", otherNames[i], " (~", win, "-kb shift)"))
  mtext(side = 3, at = 2, text = paste0("kss SPO11oligo hotspots arm vs ", otherNames[i], " (~", win, "-kb shift)"))
  dev.off()
}

