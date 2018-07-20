# Generate SPO11-1-oligo read counts for each peak and ranLoc (with equivalent width distribution to IDR-validated SPO11-1-oligos peaks)
# and plot read counts against locus width

library(GenomicAlignments)
library(ShortRead)
library(rtracklayer)
library(regioneR)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
seqlevels(genome) <- sub("Chr", "", seqlevels(genome))

inDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/"
RangedDataDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/peak_read_counts/"
outDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/peak_read_counts/RPKM/"

load(paste0(inDir, "armrangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData"))
load(paste0(inDir, "perirangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData"))
peaksGR <- sort(c(armrangerPeaksGRmerge, perirangerPeaksGRmerge))
armrangerPeaksGRmerge <- NULL
perirangerPeaksGRmerge <- NULL
seqlevels(peaksGR) <- sub("Chr", "", seqlevels(peaksGR))
#peaksGR <- peaksGR[width(peaksGR) <= 5000]
ranLocGR <- randomizeRegions(peaksGR, genome = genome, per.chromosome = TRUE, allow.overlaps = TRUE)

strand(peaksGR) <- "*"
strand(ranLocGR) <- "*"

print("print(length(peaksGR))")
print(length(peaksGR))
#[1] 5914
print("print(mean(width(peaksGR)))")
print(mean(width(peaksGR)))
#[1] 823.36
print("print(median(width(peaksGR)))")
print(median(width(peaksGR)))
#[1] 729
print("print(range(width(peaksGR)))")
print(range(width(peaksGR)))
#[1]  124 3306

# Load lib_ranged object 
load("/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/RPI1_RPI8_profiles_rangerPeaks_idr0.05/arms/peak_read_counts/WT_SPO11-oligo_RPI1_RangedData.RData")

# Calculate library size
chr_size <- NULL
for(i in 1:5) {
  print(i)
  chr_lib_ranged <- lib_ranged[i]
  chr_size <- c(chr_size, length(space(chr_lib_ranged)))
}
lib_size <- sum(chr_size)
# Calculate "per million" scaling factor
RPM_scaling_factor <- lib_size/1e+06

# Calculate RPM and RPKM for each hotspot and ranLoc
lib_rangedGR <- as(lib_ranged, "GRanges")

peak_reads <- countOverlaps(peaksGR, lib_rangedGR)
peak_RPM <- peak_reads/RPM_scaling_factor
peak_RPKM <- peak_RPM/(width(peaksGR)/1e+03)
peak_RPMplus1 <- peak_RPM+1
peak_RPKMplus1 <- peak_RPKM+1

ranLoc_reads <- countOverlaps(ranLocGR, lib_rangedGR)
ranLoc_RPM <- ranLoc_reads/RPM_scaling_factor
ranLoc_RPKM <- ranLoc_RPM/(width(ranLocGR)/1e+03)
ranLoc_RPMplus1 <- ranLoc_RPM+1
ranLoc_RPKMplus1 <- ranLoc_RPKM+1

# Calcualte TPM (transcripts per kilobase per million) for each hotspot and ranLoc
peak_RPK <- peak_reads/(width(peaksGR)/1e+03)
RPKPM_scaling_factor <- sum(peak_RPK)/1e+06
peak_TPM <- peak_RPK/RPKPM_scaling_factor

ranLoc_RPK <- ranLoc_reads/(width(ranLocGR)/1e+03)
RPKPM_scaling_factor <- sum(ranLoc_RPK)/1e+06
ranLoc_TPM <- ranLoc_RPK/RPKPM_scaling_factor

# Plot genome ranLoc width histogram, and SPO11-1-oligo RPKM+1 vs ranLoc width and cumulative fraction of random loci
pdf(file = paste0(outDir, "random_locus_width_hist_and_RPKMplus1_vs_random_locus_width_and_ecdf.pdf"), height = 4, width = 12)
par(mfrow = c(1, 3), mar =  c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
hist(width(ranLocGR), breaks = 250, col = "grey60", border = NA, lwd = 2,
     xlab = "Random locus width (bp)", ylab = "Random loci", main = "", cex.lab = 2, cex.axis = 2)
abline(v = mean(width(ranLocGR)), col = "black", lty = 2, lwd = 1)

plot(x = ranLoc_RPKMplus1, y = width(ranLocGR), pch = ".", log = "xy",
     xlim = c(1, 100), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "SPO11-1-oligo count (RPKM+1)", ylab = "Random locus width (bp)", cex.lab = 2)
axis(side = 1, at = c(1:10, seq(20, 100, by = 10)), labels = c("1", rep("", 8), "10", rep("", 8), expression("10"^"2")), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = c(seq(10, 90, by = 10), seq(100, 900, by = 100), seq(1000, 10000, by = 1000)), labels = c("10", rep("", 8), expression("10"^"2"), rep("", 8), expression("10"^"3"), rep("", 8), expression("10"^"4")), lwd.tick = 2, cex.axis = 2)
abline(h = mean(width(ranLocGR)), col = "black", lty = 2, lwd = 1)
box(lwd = 2)

plot(ecdf(ranLoc_RPKMplus1), log = "x", xlim = c(1, 100), xaxt = "n", yaxt = "n", pch = 20,
     xlab = "SPO11-1-oligo count (RPKM+1)", ylab = "Cumulative fraction of random loci", main = "", cex.lab = 2)
axis(side = 1, at = c(1:10, seq(20, 100, by = 10)), labels = c("1", rep("", 8), "10", rep("", 8), expression("10"^"2")), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = seq(0, 1, by = 0.25), labels = c("0", "", "0.5", "", "1"), lwd.tick = 2, cex.axis = 2)
box(lwd = 2)
dev.off()

# Plot SPO11-1-oligo genome hotspot width histogram, and SPO11-1-oligo RPKM vs hotspot width and cumulative fraction of loci
pdf(file = paste0(outDir, "random_locus_width_hist_and_RPKM_vs_random_locus_width_and_ecdf.pdf"), height = 4, width = 12)
par(mfrow = c(1, 3), mar =  c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
hist(width(ranLocGR), breaks = 250, col = "grey60", border = NA, lwd = 2,
     xlab = "Random locus width (bp)", ylab = "Random loci", main = "", cex.lab = 2, cex.axis = 2)
abline(v = mean(width(ranLocGR)), col = "black", lty = 2, lwd = 1)

plot(x = ranLoc_RPKM, y = width(ranLocGR), pch = ".", log = "xy",
     xlim = c(1, 100), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "SPO11-1-oligo count (RPKM)", ylab = "Random locus width (bp)", cex.lab = 2)
axis(side = 1, at = c(1:10, seq(20, 100, by = 10)), labels = c("1", rep("", 8), "10", rep("", 8), expression("10"^"2")), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = c(seq(10, 90, by = 10), seq(100, 900, by = 100), seq(1000, 10000, by = 1000)), labels = c("10", rep("", 8), expression("10"^"2"), rep("", 8), expression("10"^"3"), rep("", 8), expression("10"^"4")), lwd.tick = 2, cex.axis = 2)
abline(h = mean(width(ranLocGR)), col = "black", lty = 2, lwd = 1)
box(lwd = 2)

plot(ecdf(ranLoc_RPKM), log = "x", xlim = c(1, 100), xaxt = "n", yaxt = "n", pch = 20,
     xlab = "SPO11-1-oligo count (RPKM)", ylab = "Cumulative fraction of random loci", main = "", cex.lab = 2)
axis(side = 1, at = c(1:10, seq(20, 100, by = 10)), labels = c("1", rep("", 8), "10", rep("", 8), expression("10"^"2")), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = seq(0, 1, by = 0.25), labels = c("0", "", "0.5", "", "1"), lwd.tick = 2, cex.axis = 2)
box(lwd = 2)
dev.off()

# Plot SPO11-1-oligo genome hotspot width histogram, and SPO11-1-oligo TPM+1 vs hotspot width and cumulative fraction of loci
pdf(file = paste0(outDir, "random_locus_width_hist_and_TPMplus1_vs_random_locus_width_and_ecdf.pdf"), height = 4, width = 12)
par(mfrow = c(1, 3), mar =  c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
hist(width(ranLocGR), breaks = 250, col = "grey60", border = NA, lwd = 2,
     xlab = "Random locus width (bp)", ylab = "Random loci", main = "", cex.lab = 2, cex.axis = 2)
abline(v = mean(width(ranLocGR)), col = "black", lty = 2, lwd = 1)

plot(x = ranLoc_TPMplus1, y = width(ranLocGR), pch = ".", log = "xy",
     xlim = c(1, 100), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "SPO11-1-oligo count (TPM+1)", ylab = "Random locus width (bp)", cex.lab = 2)
axis(side = 1, at = c(1:10, seq(20, 100, by = 10)), labels = c("1", rep("", 8), "10", rep("", 8), expression("10"^"2")), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = c(seq(10, 90, by = 10), seq(100, 900, by = 100), seq(1000, 10000, by = 1000)), labels = c("10", rep("", 8), expression("10"^"2"), rep("", 8), expression("10"^"3"), rep("", 8), expression("10"^"4")), lwd.tick = 2, cex.axis = 2)
abline(h = mean(width(ranLocGR)), col = "black", lty = 2, lwd = 1)
box(lwd = 2)

plot(ecdf(ranLoc_TPMplus1), log = "x", xlim = c(1, 100), xaxt = "n", yaxt = "n", pch = 20,
     xlab = "SPO11-1-oligo count (TPM+1)", ylab = "Cumulative fraction of random loci", main = "", cex.lab = 2)
axis(side = 1, at = c(1:10, seq(20, 100, by = 10)), labels = c("1", rep("", 8), "10", rep("", 8), expression("10"^"2")), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = seq(0, 1, by = 0.25), labels = c("0", "", "0.5", "", "1"), lwd.tick = 2, cex.axis = 2)
box(lwd = 2)
dev.off()

# Plot SPO11-1-oligo genome hotspot width histogram, and SPO11-1-oligo TPM vs hotspot width and cumulative fraction of loci
pdf(file = paste0(outDir, "random_locus_width_hist_and_TPM_vs_random_locus_width_and_ecdf.pdf"), height = 4, width = 12)
par(mfrow = c(1, 3), mar =  c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
hist(width(ranLocGR), breaks = 250, col = "grey60", border = NA, lwd = 2,
     xlab = "Random locus width (bp)", ylab = "Random loci", main = "", cex.lab = 2, cex.axis = 2)
abline(v = mean(width(ranLocGR)), col = "black", lty = 2, lwd = 1)

plot(x = ranLoc_TPM, y = width(ranLocGR), pch = ".", log = "xy",
     xlim = c(1, 100), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "SPO11-1-oligo count (TPM)", ylab = "Random locus width (bp)", cex.lab = 2)
axis(side = 1, at = c(1:10, seq(20, 100, by = 10)), labels = c("1", rep("", 8), "10", rep("", 8), expression("10"^"2")), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = c(seq(10, 90, by = 10), seq(100, 900, by = 100), seq(1000, 10000, by = 1000)), labels = c("10", rep("", 8), expression("10"^"2"), rep("", 8), expression("10"^"3"), rep("", 8), expression("10"^"4")), lwd.tick = 2, cex.axis = 2)
abline(h = mean(width(ranLocGR)), col = "black", lty = 2, lwd = 1)
box(lwd = 2)

plot(ecdf(ranLoc_TPM), log = "x", xlim = c(1, 100), xaxt = "n", yaxt = "n", pch = 20,
     xlab = "SPO11-1-oligo count (TPM)", ylab = "Cumulative fraction of random loci", main = "", cex.lab = 2)
axis(side = 1, at = c(1:10, seq(20, 100, by = 10)), labels = c("1", rep("", 8), "10", rep("", 8), expression("10"^"2")), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = seq(0, 1, by = 0.25), labels = c("0", "", "0.5", "", "1"), lwd.tick = 2, cex.axis = 2)
box(lwd = 2)
dev.off()


# Plot genome hotspots width histogram, and SPO11-1-oligo RPKM+1 or RPM+1 vs loci width and cumulative fraction of loci (hotspots and random)
pdf(file = paste0(outDir, "Hotspot_width_hist_and_RPKMplus1_or_RPMplus1_vs_locus_width_and_ecdf.pdf"), height = 8, width = 12)
par(mfrow = c(2, 3), mar =  c(6, 6, 2, 2), mgp = c(4, 1.5, 0))
# RPKM+1
hist(width(peaksGR), breaks = 250, col = "grey60", border = NA, lwd = 2,
     xlab = "Hotspot width (bp)", ylab = "Hotspots", main = "", cex.lab = 2, cex.axis = 2)
abline(v = mean(width(peaksGR)), col = "black", lty = 2, lwd = 1)

plot(x = ranLoc_RPKMplus1, y = width(ranLocGR), pch = ".", log = "xy",
     xlim = c(1, 100), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = "red")
axis(side = 1, at = c(1:10, seq(20, 100, by = 10)), labels = c("1", rep("", 8), "10", rep("", 8), expression("10"^"2")), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = c(seq(10, 90, by = 10), seq(100, 900, by = 100), seq(1000, 10000, by = 1000)), labels = c("10", rep("", 8), expression("10"^"2"), rep("", 8), expression("10"^"3"), rep("", 8), expression("10"^"4")), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(x = peak_RPKMplus1, y = width(peaksGR), pch = ".", log = "xy",
     xlim = c(1, 100), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "SPO11-1-oligo count (RPKM+1)", ylab = "Locus width (bp)", cex.lab = 2)
abline(h = mean(width(peaksGR)), col = "black", lty = 2, lwd = 1)
box(lwd = 2)

plot(ecdf(ranLoc_RPKMplus1), log = "x", xlim = c(1, 100), xaxt = "n", yaxt = "n", pch = 20,
     xlab = "", ylab = "", main = "", col = "red")
axis(side = 1, at = c(1:10, seq(20, 100, by = 10)), labels = c("1", rep("", 8), "10", rep("", 8), expression("10"^"2")), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = seq(0, 1, by = 0.25), labels = c("0", "", "0.5", "", "1"), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(ecdf(peak_RPKMplus1), log = "x", xlim = c(1, 100), xaxt = "n", yaxt = "n", pch = 20,   
     xlab = "SPO11-1-oligo count (RPKM+1)", ylab = "Cumulative fraction of loci", main = "", cex.lab = 2)
box(lwd = 2)

# RPM+1
hist(width(peaksGR), breaks = 250, col = "grey60", border = NA, lwd = 2,
     xlab = "Hotspot width (bp)", ylab = "Hotspots", main = "", cex.lab = 2, cex.axis = 2)
abline(v = mean(width(peaksGR)), col = "black", lty = 2, lwd = 1)

plot(x = ranLoc_RPMplus1, y = width(ranLocGR), pch = ".", log = "xy",
     xlim = c(1, 100), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", col = "red")
axis(side = 1, at = c(1:10, seq(20, 100, by = 10)), labels = c("1", rep("", 8), "10", rep("", 8), expression("10"^"2")), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = c(seq(10, 90, by = 10), seq(100, 900, by = 100), seq(1000, 10000, by = 1000)), labels = c("10", rep("", 8), expression("10"^"2"), rep("", 8), expression("10"^"3"), rep("", 8), expression("10"^"4")), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(x = peak_RPMplus1, y = width(peaksGR), pch = ".", log = "xy",
     xlim = c(1, 100), ylim = c(10, 10000), xaxt = "n", yaxt = "n",
     xlab = "SPO11-1-oligo count (RPM+1)", ylab = "Locus width (bp)", cex.lab = 2)
abline(h = mean(width(peaksGR)), col = "black", lty = 2, lwd = 1)
box(lwd = 2)

plot(ecdf(ranLoc_RPMplus1), log = "x", xlim = c(1, 100), xaxt = "n", yaxt = "n", pch = 20,
     xlab = "", ylab = "", main = "", col = "red")
axis(side = 1, at = c(1:10, seq(20, 100, by = 10)), labels = c("1", rep("", 8), "10", rep("", 8), expression("10"^"2")), lwd.tick = 2, cex.axis = 2)
axis(side = 2, at = seq(0, 1, by = 0.25), labels = c("0", "", "0.5", "", "1"), lwd.tick = 2, cex.axis = 2)
par(new = T)
plot(ecdf(peak_RPMplus1), log = "x", xlim = c(1, 100), xaxt = "n", yaxt = "n", pch = 20,
     xlab = "SPO11-1-oligo count (RPM+1)", ylab = "Cumulative fraction of loci", main = "", cex.lab = 2)
box(lwd = 2)

dev.off()


