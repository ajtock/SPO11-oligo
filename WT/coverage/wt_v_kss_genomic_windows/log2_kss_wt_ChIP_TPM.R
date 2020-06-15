#!/applications/R/R-3.5.0/bin/Rscript

# contact: ajt200@cam.ac.uk (Andy Tock)

# Calculate ChIP-seq TPM values and mean DNA methylation values
# in genomic windows, and identify windows in which log2(kss/wt)
# SPO11-1-oligos levels exceed log2FC

# Usage on hydrogen node7:
# /scripts/csmit -m 20G -c 1 "Rscript ./log2_kss_wt_ChIP_TPM.R 1000 1 SPO11oligos"

#winSize <- 1000
#L2FCthreshold <- 1
#L2FCfactor <- "SPO11oligos"
args <- commandArgs(trailingOnly = T)
winSize <- as.numeric(args[1])
L2FCthreshold <- as.numeric(args[2])
L2FCfactor <- args[3]

outDir <- paste0(L2FCfactor, "/")
gainDir <- paste0(outDir, "gain/")
lossDir <- paste0(outDir, "loss/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", gainDir, " ] || mkdir ", gainDir))
system(paste0("[ -d ", lossDir, " ] || mkdir ", lossDir))

library(GenomicAlignments)
library(dplyr)
library(ggplot2)

# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

# Define windows as GRanges object
windowsGR <- GRanges()
for(i in 1:length(chrs)) {
  seqWindows <- seq(1, chrLens[i], by = winSize)
  windowsIR <- IRanges(start = seqWindows,
                       width = winSize)
  windowsIR <- windowsIR[-length(windowsIR)]
  windowsIR <- append(windowsIR,
                      IRanges(start = seqWindows[length(seqWindows)],
                              end = chrLens[i]))
  chrWindowsGR <- GRanges(seqnames = chrs[i],
                          ranges = windowsIR,
                          strand = "*")
  print(chrWindowsGR)
  windowsGR <- append(windowsGR, chrWindowsGR)
}

# Load BAMs and create GAlignment (single-end reads) and
# GAlignmentPairs (paired-end reads) objects
# Single-end
wt_SPO11oligos_Rep1 <- readGAlignments("reads/wt_SPO11oligos_Rep1.bam")
wt_SPO11oligos_Rep2 <- readGAlignments("reads/wt_SPO11oligos_Rep2.bam")
wt_SPO11oligos_Rep3 <- readGAlignments("reads/wt_SPO11oligos_Rep3.bam")
kss_SPO11oligos_Rep1 <- readGAlignments("reads/kss_SPO11oligos_Rep1.bam")
kss_SPO11oligos_Rep2 <- readGAlignments("reads/kss_SPO11oligos_Rep2.bam")
gDNA_Rep1 <- readGAlignments("reads/wt_gDNA_Rep1_R1.bam")
# Paired-end
wt_REC8_HA_Rep1 <- readGAlignmentPairs("reads/wt_REC8_HA_Rep1_ChIP.bam")
wt_REC8_HA_Rep2 <- readGAlignmentPairs("reads/wt_REC8_HA_Rep2_ChIP.bam")
wt_REC8_Myc_Rep1 <- readGAlignmentPairs("reads/wt_REC8_Myc_Rep1_ChIP.bam")
kss_REC8_HA_Rep1 <- readGAlignmentPairs("reads/kss_REC8_HA_Rep1_ChIP.bam")
kss_REC8_HA_Rep2 <- readGAlignmentPairs("reads/kss_REC8_HA_Rep2_ChIP.bam")
wt_H3K9me2_Rep1 <- readGAlignmentPairs("reads/wt_H3K9me2_Rep1_ChIP.bam")
kss_H3K9me2_Rep1 <- readGAlignmentPairs("reads/kss_H3K9me2_Rep1_ChIP.bam")
input_Myc_Rep1 <- readGAlignmentPairs("reads/wt_REC8_Myc_Rep1_input.bam")

# Convert into GRanges
# SPO11oligos
wt_SPO11oligos_Rep1GR <- GRanges(wt_SPO11oligos_Rep1[seqnames(wt_SPO11oligos_Rep1) %in%
                                                     sub("Chr", "", chrs)])
seqlevels(wt_SPO11oligos_Rep1GR) <- sub("", "Chr", seqlevels(wt_SPO11oligos_Rep1GR))
wt_SPO11oligos_Rep2GR <- GRanges(wt_SPO11oligos_Rep2[seqnames(wt_SPO11oligos_Rep2) %in%
                                                     sub("Chr", "", chrs)])
seqlevels(wt_SPO11oligos_Rep2GR) <- sub("", "Chr", seqlevels(wt_SPO11oligos_Rep2GR))
wt_SPO11oligos_Rep3GR <- GRanges(wt_SPO11oligos_Rep3[seqnames(wt_SPO11oligos_Rep3) %in%
                                                     sub("Chr", "", chrs)])
seqlevels(wt_SPO11oligos_Rep3GR) <- sub("", "Chr", seqlevels(wt_SPO11oligos_Rep3GR))
kss_SPO11oligos_Rep1GR <- GRanges(kss_SPO11oligos_Rep1[seqnames(kss_SPO11oligos_Rep1) %in%
                                                       sub("Chr", "", chrs)])
seqlevels(kss_SPO11oligos_Rep1GR) <- sub("", "Chr", seqlevels(kss_SPO11oligos_Rep1GR))
kss_SPO11oligos_Rep2GR <- GRanges(kss_SPO11oligos_Rep2[seqnames(kss_SPO11oligos_Rep2) %in%
                                                       sub("Chr", "", chrs)])
seqlevels(kss_SPO11oligos_Rep2GR) <- sub("", "Chr", seqlevels(kss_SPO11oligos_Rep2GR))

# REC8
wt_REC8_HA_Rep1GR <- GRanges(wt_REC8_HA_Rep1[seqnames(wt_REC8_HA_Rep1) %in%
                                             sub("Chr", "", chrs)])
seqlevels(wt_REC8_HA_Rep1GR) <- sub("", "Chr", seqlevels(wt_REC8_HA_Rep1GR))
wt_REC8_HA_Rep2GR <- GRanges(wt_REC8_HA_Rep2[seqnames(wt_REC8_HA_Rep2) %in%
                                             sub("Chr", "", chrs)])
seqlevels(wt_REC8_HA_Rep2GR) <- sub("", "Chr", seqlevels(wt_REC8_HA_Rep2GR))
wt_REC8_Myc_Rep1GR <- GRanges(wt_REC8_Myc_Rep1[seqnames(wt_REC8_Myc_Rep1) %in%
                                               sub("Chr", "", chrs)])
seqlevels(wt_REC8_Myc_Rep1GR) <- sub("", "Chr", seqlevels(wt_REC8_Myc_Rep1GR))
kss_REC8_HA_Rep1GR <- GRanges(kss_REC8_HA_Rep1[seqnames(kss_REC8_HA_Rep1) %in%
                                               sub("Chr", "", chrs)])
seqlevels(kss_REC8_HA_Rep1GR) <- sub("", "Chr", seqlevels(kss_REC8_HA_Rep1GR))
kss_REC8_HA_Rep2GR <- GRanges(kss_REC8_HA_Rep2[seqnames(kss_REC8_HA_Rep2) %in%
                                               sub("Chr", "", chrs)])
seqlevels(kss_REC8_HA_Rep2GR) <- sub("", "Chr", seqlevels(kss_REC8_HA_Rep2GR))

# H3K9me2
wt_H3K9me2_Rep1GR <- GRanges(wt_H3K9me2_Rep1[seqnames(wt_H3K9me2_Rep1) %in%
                                             sub("Chr", "", chrs)])
seqlevels(wt_H3K9me2_Rep1GR) <- sub("", "Chr", seqlevels(wt_H3K9me2_Rep1GR))
kss_H3K9me2_Rep1GR <- GRanges(kss_H3K9me2_Rep1[seqnames(kss_H3K9me2_Rep1) %in%
                                               sub("Chr", "", chrs)])
seqlevels(kss_H3K9me2_Rep1GR) <- sub("", "Chr", seqlevels(kss_H3K9me2_Rep1GR))

# controls
input_Rep1GR <- GRanges(input_Rep1[seqnames(input_Rep1) %in%
                                   sub("Chr", "", chrs)])
seqlevels(input_Rep1GR) <- sub("", "Chr", seqlevels(input_Rep1GR))
gDNA_Rep1GR <- GRanges(gDNA_Rep1[seqnames(gDNA_Rep1) %in%
                                 sub("Chr", "", chrs)])
seqlevels(gDNA_Rep1GR) <- sub("", "Chr", seqlevels(gDNA_Rep1GR))


# Function to calculate feature TPM values for a given library
featureCovCalc <- function(features, reads, control) {
  # ChIP reads
  feature_reads <- countOverlaps(query = features,
                                 subject = reads,
                                 type = "any",
                                 ignore.strand = TRUE)
  feature_RPK <- feature_reads/(width(features)/1e3)
  RPKM_scaling_factor <- sum(feature_RPK)/1e6
  feature_TPM <- feature_RPK/RPKM_scaling_factor
  # control reads  
  feature_reads_control <- countOverlaps(query = features,
                                         subject = control,
                                         type = "any",
                                         ignore.strand = TRUE)
  feature_RPK_control <- feature_reads_control/(width(features)/1e3)
  RPKM_scaling_factor_control <- sum(feature_RPK_control)/1e6
  feature_TPM_control <- feature_RPK_control/RPKM_scaling_factor_control
  # Subtract control TPM from ChIP TPM and replace values < 0 with 0
  feature_TPM_ChIP_less_control <- feature_TPM-feature_TPM_control
  feature_TPM_ChIP_less_control[which(feature_TPM_ChIP_less_control < 0)] <- 0
  # Combine in dataframe
  data.frame(reads = as.integer(feature_reads),
             reads_control = as.integer(feature_reads_control),
             ChIP_TPM =  as.numeric(feature_TPM),
             control_TPM = as.numeric(feature_TPM_control),
             TPM = as.numeric(feature_TPM_ChIP_less_control),
             stringsAsFactors = F)
}

# Apply featureCovCalc() to each library,
# calculate mean TPM for each window (or feature),
# and calculate log2(kss/wt) TPM for each window
# SPO11oligos
wt_SPO11oligos_Rep1_winCov <- featureCovCalc(features = windowsGR,
                                             reads = wt_SPO11oligos_Rep1GR,
                                             control = gDNA_Rep1GR)
wt_SPO11oligos_Rep2_winCov <- featureCovCalc(features = windowsGR,
                                             reads = wt_SPO11oligos_Rep2GR,
                                             control = gDNA_Rep1GR)
wt_SPO11oligos_Rep3_winCov <- featureCovCalc(features = windowsGR,
                                             reads = wt_SPO11oligos_Rep3GR,
                                             control = gDNA_Rep1GR)
kss_SPO11oligos_Rep1_winCov <- featureCovCalc(features = windowsGR,
                                              reads = kss_SPO11oligos_Rep1GR,
                                              control = gDNA_Rep1GR)
kss_SPO11oligos_Rep2_winCov <- featureCovCalc(features = windowsGR,
                                              reads = kss_SPO11oligos_Rep2GR,
                                              control = gDNA_Rep1GR)

SPO11oligos_mean_TPM <- sapply(seq_along(windowsGR), function(x) {
  mean(c(wt_SPO11oligos_Rep1_winCov$TPM[x],
         wt_SPO11oligos_Rep2_winCov$TPM[x],
         wt_SPO11oligos_Rep3_winCov$TPM[x],
         kss_SPO11oligos_Rep1_winCov$TPM[x],
         kss_SPO11oligos_Rep2_winCov$TPM[x]))
})

log2_kss_wt_SPO11oligos_TPM <- sapply(seq_along(windowsGR), function(x) {
  log2(
    ( mean(c(wt_SPO11oligos_Rep1_winCov$TPM[x],
             wt_SPO11oligos_Rep2_winCov$TPM[x],
             wt_SPO11oligos_Rep3_winCov$TPM[x])) + 1 ) /
    ( mean(c(kss_SPO11oligos_Rep1_winCov$TPM[x],
             kss_SPO11oligos_Rep2_winCov$TPM[x])) + 1 )
  )
})

# REC8
wt_REC8_HA_Rep1_winCov <- featureCovCalc(features = windowsGR,
                                         reads = wt_REC8_HA_Rep1GR,
                                         control = input_Rep1GR)
wt_REC8_HA_Rep2_winCov <- featureCovCalc(features = windowsGR,
                                         reads = wt_REC8_HA_Rep2GR,
                                         control = input_Rep1GR)
wt_REC8_Myc_Rep1_winCov <- featureCovCalc(features = windowsGR,
                                          reads = wt_REC8_Myc_Rep1GR,
                                          control = input_Rep1GR)
kss_REC8_HA_Rep1_winCov <- featureCovCalc(features = windowsGR,
                                          reads = kss_REC8_HA_Rep1GR,
                                          control = input_Rep1GR)
kss_REC8_HA_Rep2_winCov <- featureCovCalc(features = windowsGR,
                                          reads = kss_REC8_HA_Rep2GR,
                                          control = input_Rep1GR)

REC8_HA_mean_TPM <- sapply(seq_along(windowsGR), function(x) {
  mean(c(wt_REC8_HA_Rep1_winCov$TPM[x],
         wt_REC8_HA_Rep2_winCov$TPM[x],
         kss_REC8_HA_Rep1_winCov$TPM[x],
         kss_REC8_HA_Rep2_winCov$TPM[x]))
})

log2_kss_wt_REC8_HA_TPM <- sapply(seq_along(windowsGR), function(x) {
  log2(
    ( mean(c(wt_REC8_HA_Rep1_winCov$TPM[x],
             wt_REC8_HA_Rep2_winCov$TPM[x])) + 1 ) /
    ( mean(c(kss_REC8_HA_Rep1_winCov$TPM[x],
             kss_REC8_HA_Rep2_winCov$TPM[x])) + 1 )
  )
})

# H3K9me2
wt_H3K9me2_Rep1_winCov <- featureCovCalc(features = windowsGR,
                                         reads = wt_H3K9me2_Rep1GR,
                                         control = input_Rep1GR)
kss_H3K9me2_Rep1_winCov <- featureCovCalc(features = windowsGR,
                                          reads = kss_H3K9me2_Rep1GR,
                                          control = input_Rep1GR)

H3K9me2_mean_TPM <- sapply(seq_along(windowsGR), function(x) {
  mean(c(wt_H3K9me2_Rep1_winCov$TPM[x],
         kss_H3K9me2_Rep1_winCov$TPM[x]))
})

log2_kss_wt_H3K9me2_TPM <- sapply(seq_along(windowsGR), function(x) {
  log2( 
    ( wt_H3K9me2_Rep1_winCov$TPM[x] + 1 ) /
    ( kss_H3K9me2_Rep1_winCov$TPM[x] + 1 )
  )
})

# Combine wt and kss TPM values for each library in one dataframe
featuresDF <- data.frame(windowsGR,
                         wt_SPO11oligos_Rep1_TPM = wt_SPO11oligos_Rep1_winCov$TPM,
                         wt_SPO11oligos_Rep2_TPM = wt_SPO11oligos_Rep2_winCov$TPM,
                         wt_SPO11oligos_Rep3_TPM = wt_SPO11oligos_Rep3_winCov$TPM,
                         kss_SPO11oligos_Rep1_TPM = kss_SPO11oligos_Rep1_winCov$TPM,
                         kss_SPO11oligos_Rep2_TPM = kss_SPO11oligos_Rep2_winCov$TPM,
                         SPO11oligos_mean_TPM = SPO11oligos_mean_TPM,
                         log2_kss_wt_SPO11oligos_TPM = log2_kss_wt_SPO11oligos_TPM,
                         wt_REC8_HA_Rep1_TPM = wt_REC8_HA_Rep1_winCov$TPM,
                         wt_REC8_HA_Rep2_TPM = wt_REC8_HA_Rep2_winCov$TPM,
                         wt_REC8_Myc_Rep1_TPM = wt_REC8_Myc_Rep1_winCov$TPM,
                         kss_REC8_HA_Rep1_TPM = kss_REC8_HA_Rep1_winCov$TPM,
                         kss_REC8_HA_Rep2_TPM = kss_REC8_HA_Rep2_winCov$TPM,
                         REC8_HA_mean_TPM = REC8_HA_mean_TPM,
                         log2_kss_wt_REC8_HA_TPM = log2_kss_wt_REC8_HA_TPM,
                         wt_H3K9me2_Rep1_TPM = wt_H3K9me2_Rep1_winCov$TPM,
                         kss_H3K9me2_Rep1_TPM = kss_H3K9me2_Rep1_winCov$TPM,
                         H3K9me2_mean_TPM = H3K9me2_mean_TPM,
                         log2_kss_wt_H3K9me2_TPM = log2_kss_wt_H3K9me2_TPM,
                         stringsAsFactors = F)
                       
# Select features with greatest gain and loss of L2FCfactor in kss
featuresDF_kssGain <- featuresDF[featuresDF[,grepl(paste0("log2_kss_wt_", L2FCfactor),
                                                   colnames(featuresDF))] >= L2FCthreshold,]
featuresDF_kssLoss <- featuresDF[featuresDF[,grepl(paste0("log2_kss_wt_", L2FCfactor),
                                                   colnames(featuresDF))] <= -L2FCthreshold,]
write.table(featuresDF_kssGain,
            file = paste0(outDir,
                          winSize/1e3, "kb_genomic_windows_with_log2_kss_wt_", L2FCfactor,
                          "_TPM_moreThanOrEqualToPlus", as.character(L2FCthreshold),
                          ".tsv"),
            quote = F, sep = "\t", col.names = T, row.names = F)
write.table(featuresDF_kssLoss,
            file = paste0(outDir,
                          winSize/1e3, "kb_genomic_windows_with_log2_kss_wt_", L2FCfactor,
                          "_TPM_lessThanOrEqualToMinus", as.character(L2FCthreshold),
                          ".tsv"),
            quote = F, sep = "\t", col.names = T, row.names = F)


cov_columns <- which(grepl("Rep", colnames(featuresDF_kssGain)))
cov_columnNames <- colnames(featuresDF_kssGain)[cov_columns]
cov_columnNamesPlot <- gsub("_", " ", cov_columnNames)
cov_columnNamesPlot <- sub(" TPM", "", cov_columnNamesPlot)
cov_columnNamesPlot <- sub("REC8 ", "REC8-", cov_columnNamesPlot)
featuresDF_kssGain_mean <- sapply(cov_columns, function(x) {
  mean(featuresDF_kssGain[,colnames(featuresDF_kssGain) == colnames(featuresDF_kssGain)[x]])
})
featuresDF_kssGain_SD <- sapply(cov_columns, function(x) {
  sd(featuresDF_kssGain[,colnames(featuresDF_kssGain) == colnames(featuresDF_kssGain)[x]])
})
featuresDF_kssGain_SEM <- sapply(seq_along(featuresDF_kssGain_SD), function(x) {
  featuresDF_kssGain_SD[x] / sqrt( (dim(featuresDF_kssGain)[1] - 1) )
})
featuresDF_kssGain_CIlower <- sapply(seq_along(featuresDF_kssGain_mean), function(x) {
  featuresDF_kssGain_mean[x] -
    ( qt(0.975, df = dim(featuresDF_kssGain)[1] - 1) *
      featuresDF_kssGain_SEM[x] )
})
featuresDF_kssGain_CIupper <- sapply(seq_along(featuresDF_kssGain_mean), function(x) {
  featuresDF_kssGain_mean[x] +
    ( qt(0.975, df = dim(featuresDF_kssGain)[1] - 1) *
      featuresDF_kssGain_SEM[x] )
})

featuresDF_quantileSEM <- sapply(1:quantiles, function(k) {
  featuresDF_quantileSD[k] / sqrt( (dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k),])[1] - 1) )
})
featuresDF_quantileCIlower <- sapply(1:quantiles, function(k) {
  featuresDF_quantileMean[k] -
    ( qt(0.975, df = dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k),])[1]-1 ) *
      featuresDF_quantileSEM[k] )
})
featuresDF_quantileCIupper <- sapply(1:quantiles, function(k) {
  featuresDF_quantileMean[k] +
    ( qt(0.975, df = dim(featuresDF[featuresDF$quantile == paste0("Quantile ", k),])[1]-1 ) *
      featuresDF_quantileSEM[k] )
})
featuresDF_summary_stats <- data.frame(quantile = paste0("Quantile ", 1:quantiles),
                                       Mean = featuresDF_quantileMean,
                                       SD = featuresDF_quantileSD,
                                       SEM = featuresDF_quantileSEM,
                                       CIlower = featuresDF_quantileCIlower,
                                       CIupper = featuresDF_quantileCIupper,
                                       stringsAsFactors = F)

# Make MA plot           
feature_meanTPM_log2TPM_df <- as_data_frame(cbind(mean_TPM = featuresDF[,grepl(paste0(L2FCfactor, "_mean_TPM"),
                                                                               colnames(featuresDF))],
                                                  log2_TPM = featuresDF[,grepl(paste0("log2_kss_wt_", L2FCfactor),
                                                                               colnames(featuresDF))]))
# MA plotting function
MAplotFun <- function(L2FCthreshold, SCMvalues) {
 ggplot(data = feature_meanTPM_log2TPM_df,
        mapping = aes(x = mean_TPM,
                      y = log2_TPM)) +
 geom_point(aes(colour = cut(log2_TPM,
                             c(-Inf, -L2FCthreshold, L2FCthreshold-0.001, Inf))),
            shape = 20) +
 scale_color_manual(name = expression("Log"[2] * "(fold change)"),
                    values = SCMvalues,
                    labels = c(bquote("<= -" * .(L2FCthreshold)),
                               bquote("-" * .(L2FCthreshold) * " < change < " * .(L2FCthreshold)),
                               bquote(">= " * .(L2FCthreshold)))) +
 geom_hline(yintercept = c(-L2FCthreshold, L2FCthreshold),
            linetype = "dashed",
            colour = "black",
            size = 1.0) +
 xlab("Mean TPM") +
 ylab(bquote("Log"[2] * "(kss " *
             .(L2FCfactor) * "/wt " * .(L2FCfactor) *
             ") TPM")) +
 xlim(0, max(feature_meanTPM_log2TPM_df$mean_TPM)) +
 ylim(-(max(abs(min(feature_meanTPM_log2TPM_df$log2_TPM)),
            max(feature_meanTPM_log2TPM_df$log2_TPM))),
      max(abs(min(feature_meanTPM_log2TPM_df$log2_TPM)),
          max(feature_meanTPM_log2TPM_df$log2_TPM))) +
 theme_bw() +
 theme(axis.ticks.length = unit(5, "pt"),
       axis.ticks.x = element_line(size = 1, colour = "black"),
       axis.ticks.y = element_line(size = 1, colour = "black"),
       axis.text.x = element_text(size = 14, colour = "black"),
       axis.text.y = element_text(size = 14, colour = "black"),
       axis.title = element_text(size = 14, colour = "black"),
       legend.text = element_text(size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black"),
       panel.border = element_rect(size = 1.5, colour = "black"),
       panel.grid = element_blank(),
       panel.background = element_blank(),
       plot.margin = unit(c(0.3, 0.9, 0.9, 0.3), "cm"),
       plot.title = element_text(hjust = 0.5, size = 14)) +
 ggtitle(bquote("MA plot of log"[2] * "(kss " *
                .(L2FCfactor) * "/wt " * .(L2FCfactor) *
                ") TPM in " * .(winSize/1e3) * "-kb windows"))
}

MAplot <- MAplotFun(L2FCthreshold = L2FCthreshold,
                    SCMvalues = c("dodgerblue2", "grey40", "orangered"))
ggsave(MAplot,
       file = paste0(outDir,
                     winSize/1e3, "kb_genomic_windows_with_log2_kss_wt_", L2FCfactor,
                     "_TPM_L2FC", as.character(L2FCthreshold), "_MAplot.pdf"),
       height = 6, width = 10)


