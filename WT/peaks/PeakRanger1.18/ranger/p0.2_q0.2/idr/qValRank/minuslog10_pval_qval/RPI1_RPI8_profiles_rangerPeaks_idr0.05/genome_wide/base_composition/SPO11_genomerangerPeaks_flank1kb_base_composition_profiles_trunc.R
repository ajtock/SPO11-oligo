######################################################################################
# sequence composition around SPO11-1-oligo hotspots                                 #
######################################################################################

library(Biostrings)
library(segmentSeq)
library(BSgenome)
library(BSgenome.Athaliana.TAIR.TAIR9)
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
outDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/base_composition/"
plotDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/RPI1_RPI8_profiles_rangerPeaks_idr0.05/genome_wide/base_composition/plots/"

load(paste0(inDir, "armrangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData"))
load(paste0(inDir, "perirangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData"))
peaksGR <- sort(c(armrangerPeaksGRmerge, perirangerPeaksGRmerge))
armrangerPeaksGRmerge <- NULL
perirangerPeaksGRmerge <- NULL
strand(peaksGR) <- "*"
seqlevels(peaksGR) <- sub("Chr", "", seqlevels(peaksGR))
#peaksGR <- peaksGR[width(peaksGR) <= 5000]
ranLocGR <- randomizeRegions(peaksGR, genome = genome, per.chromosome = TRUE, allow.overlaps = TRUE)

print(length(peaksGR))
#[1] 5914

peaksGR <- GRanges(seqnames = seqnames(peaksGR), ranges = IRanges(start = start(peaksGR), end = end(peaksGR)), strand = "+",
                   midpoint = round((start(peaksGR))+((end(peaksGR)-start(peaksGR))/2)))
ranLocGR <- GRanges(seqnames = seqnames(ranLocGR), ranges = IRanges(start = start(ranLocGR), end = end(ranLocGR)), strand = "+",
                    midpoint = round((start(ranLocGR))+((end(ranLocGR)-start(ranLocGR))/2)))

chr1 <- Athaliana$Chr1
chr2 <- Athaliana$Chr2
chr3 <- Athaliana$Chr3
chr4 <- Athaliana$Chr4
chr5 <- Athaliana$Chr5

a.coords <- rep(0,times=2000)
t.coords <- rep(0,times=2000)
g.coords <- rep(0,times=2000)
c.coords <- rep(0,times=2000)
ran.a.coords <- rep(0,times=2000)
ran.t.coords <- rep(0,times=2000)
ran.g.coords <- rep(0,times=2000)
ran.c.coords <- rep(0,times=2000)
peak.chr.tots <- NULL
ranLoc.chr.tots <- NULL
for(i in 1:5) {
        print(i)	
	chr.seq <- Athaliana[[i]]
        chrPeaksGR <- peaksGR[seqnames(peaksGR) == i]
        peakStart <- chrPeaksGR$midpoint-1000
        peakEnd <- chrPeaksGR$midpoint+999
	print(i)
	peak.chr.tots <- c(peak.chr.tots, length(chrPeaksGR))
	print(i)
        chrRanLocGR <- ranLocGR[seqnames(ranLocGR) == i]
        ranLocStart <- chrRanLocGR$midpoint-1000
        ranLocEnd <- chrRanLocGR$midpoint+999
        ranLoc.chr.tots <- c(ranLoc.chr.tots, length(chrRanLocGR))
	print(i)
	chr.acoords <- rep(0,times=2000)
	chr.tcoords <- rep(0,times=2000)
	chr.gcoords <- rep(0,times=2000)
	chr.ccoords <- rep(0,times=2000)
	for(j in 1:length(peakStart)) {
		acoords <- rep(0,times=2000)
		tcoords <- rep(0,times=2000)
		gcoords <- rep(0,times=2000)
		ccoords <- rep(0,times=2000)
		print(j)
		sel.seq <- unlist(strsplit(as.character(chr.seq[peakStart[j]:peakEnd[j]]),split=""))
		acoords[which(sel.seq=="A")] <- 1
		tcoords[which(sel.seq=="T")] <- 1
		gcoords[which(sel.seq=="G")] <- 1
		ccoords[which(sel.seq=="C")] <- 1	   
		chr.acoords <- chr.acoords+acoords
		chr.tcoords <- chr.tcoords+tcoords
		chr.gcoords <- chr.gcoords+gcoords
		chr.ccoords <- chr.ccoords+ccoords
	}
	a.coords <- a.coords+chr.acoords
	t.coords <- t.coords+chr.tcoords
	g.coords <- g.coords+chr.gcoords
	c.coords <- c.coords+chr.ccoords
	print(i)
        chr.ran.acoords <- rep(0,times=2000)
        chr.ran.tcoords <- rep(0,times=2000)
        chr.ran.gcoords <- rep(0,times=2000)
        chr.ran.ccoords <- rep(0,times=2000)
        for(j in 1:length(peakStart)) {
                acoords <- rep(0,times=2000)
                tcoords <- rep(0,times=2000)
                gcoords <- rep(0,times=2000)
                ccoords <- rep(0,times=2000)
                print(j)
                sel.seq <- unlist(strsplit(as.character(chr.seq[ranLocStart[j]:ranLocEnd[j]]),split=""))
                acoords[which(sel.seq=="A")] <- 1
                tcoords[which(sel.seq=="T")] <- 1
                gcoords[which(sel.seq=="G")] <- 1
                ccoords[which(sel.seq=="C")] <- 1
                chr.ran.acoords <- chr.ran.acoords+acoords
                chr.ran.tcoords <- chr.ran.tcoords+tcoords
                chr.ran.gcoords <- chr.ran.gcoords+gcoords
                chr.ran.ccoords <- chr.ran.ccoords+ccoords
        }
        ran.a.coords <- ran.a.coords+chr.ran.acoords
        ran.t.coords <- ran.t.coords+chr.ran.tcoords
        ran.g.coords <- ran.g.coords+chr.ran.gcoords
        ran.c.coords <- ran.c.coords+chr.ran.ccoords
}
peak.tot <- sum(peak.chr.tots)
print(peak.tot)
ranLoc.tot <- sum(ranLoc.chr.tots)
print(ranLoc.tot)
a.coords <- a.coords/peak.tot
t.coords <- t.coords/peak.tot
g.coords <- g.coords/peak.tot
c.coords <- c.coords/peak.tot
ran.a.coords <- ran.a.coords/ranLoc.tot
ran.t.coords <- ran.t.coords/ranLoc.tot
ran.g.coords <- ran.g.coords/ranLoc.tot
ran.c.coords <- ran.c.coords/ranLoc.tot

write.table(a.coords, file = paste0(outDir, "peaks.a.txt"))
write.table(t.coords, file = paste0(outDir, "peaks.t.txt"))
write.table(g.coords, file = paste0(outDir, "peaks.g.txt"))
write.table(c.coords, file = paste0(outDir, "peaks.c.txt"))
write.table(ran.a.coords, file = paste0(outDir, "ranLoc.a.txt"))
write.table(ran.t.coords, file = paste0(outDir, "ranLoc.t.txt"))
write.table(ran.g.coords, file = paste0(outDir, "ranLoc.g.txt"))
write.table(ran.c.coords, file = paste0(outDir, "ranLoc.c.txt"))

genome.coords <- list(read.table(file = paste0(outDir, "peaks.a.txt"))[,1],
                   read.table(file = paste0(outDir, "peaks.t.txt"))[,1],
                   read.table(file = paste0(outDir, "peaks.g.txt"))[,1],
                   read.table(file = paste0(outDir, "peaks.c.txt"))[,1])
genome.ran.coords <- list(read.table(file = paste0(outDir, "ranLoc.a.txt"))[,1],
                       read.table(file = paste0(outDir, "ranLoc.t.txt"))[,1],
                       read.table(file = paste0(outDir, "ranLoc.g.txt"))[,1],
                       read.table(file = paste0(outDir, "ranLoc.c.txt"))[,1])

baseCompPlot <- function(xplot, coords, ran.coords, ylim, mycols) {
  # peaks
  plot(xplot, coords[[1]], col = mycols[1], lwd = 1.5, type = "l", ylim = ylim, main = "SPO11-1 peaks (RPI1 & RPI8; IDR < 0.05)", xlab = "", ylab = "Base relative frequency", xaxt = "n")
  lines(xplot, coords[[2]], col = mycols[2], lwd = 1.5)
  lines(xplot, coords[[3]], col = mycols[3], lwd = 1.5)
  lines(xplot, coords[[4]], col = mycols[4], lwd = 1.5)
  axis(side = 1, at = c(-1000, 0, 1000), labels = c("-1 kb", "0", "+1 kb"))
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
  plot(xplot, ran.coords[[1]], col = mycols[1], lwd = 1.5, type = "l", ylim = ylim, main = "Random loci", xlab = "", ylab = "Base relative frequency", xaxt = "n")
  lines(xplot, ran.coords[[2]], col = mycols[2], lwd = 1.5)
  lines(xplot, ran.coords[[3]], col = mycols[3], lwd = 1.5)
  lines(xplot, ran.coords[[4]], col = mycols[4], lwd = 1.5)
  axis(side = 1, at = c(-1000, 0, 1000), labels = c("-1 kb", "0", "+1 kb"))
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
       
pdf(paste0(plotDir, "SPO11genomerangerPeaks_idr0.05_base_comp_v100118.pdf"))
xplot <- seq(-1000, 999, by = 1)
mycols <- c("deepskyblue", "midnightblue", "red", "tomato4")
genome.ylim <- c(min(genome.coords[[1]], genome.coords[[2]], genome.coords[[3]], genome.coords[[4]], genome.ran.coords[[1]], genome.ran.coords[[2]], genome.ran.coords[[3]], genome.ran.coords[[4]]), 
              max(genome.coords[[1]], genome.coords[[2]], genome.coords[[3]], genome.coords[[4]], genome.ran.coords[[1]], genome.ran.coords[[2]], genome.ran.coords[[3]], genome.ran.coords[[4]]))
par(mfrow = c(2, 2))
par(mar = c(5, 4, 4, 4))
par(mgp = c(3, 0.75, 0)) 
baseCompPlot(xplot = xplot, coords = genome.coords, ran.coords = genome.ran.coords, ylim = genome.ylim, mycols = mycols)
dev.off()


# Plot (A+T) and (G+C) around SPO11-1-oligos peaks
genome.at.coords <- (genome.coords[[1]]+genome.coords[[2]])
genome.gc.coords <- (genome.coords[[3]]+genome.coords[[4]])
genome.at.ran.coords <- (genome.ran.coords[[1]]+genome.ran.coords[[2]])
genome.gc.ran.coords <- (genome.ran.coords[[3]]+genome.ran.coords[[4]])

mergeBaseCompPlot <- function(xplot, at.coords, gc.coords, at.ran.coords, gc.ran.coords, ylim, mycols) {
  # peaks
  plot(xplot, at.coords, col = mycols[1], lwd = 1.5, type = "l", ylim = ylim, main = "SPO11-1 peaks (RPI1 & RPI8; IDR < 0.05)", xlab = "", ylab = "Base relative frequency", xaxt = "n")
  lines(xplot, gc.coords, col = mycols[2], lwd = 1.5)
  axis(side = 1, at = c(-1000, 0, 1000), labels = c("-1 kb", "0", "+1 kb"))
  axis(side = 2, at = pretty(c(at.coords, gc.coords, at.ran.coords, gc.ran.coords)))
  abline(v = 0, lty = 3)
  #abline(v = c(-100, 99), lty = 3)
  legend("right",
         legend = c("AT", "GC"),
         col = c(mycols[1], mycols[2]),
         text.col = c(mycols[1], mycols[2]),
         ncol = 1, cex = 0.8, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
  # ranLoc
  plot(xplot, at.ran.coords, col = mycols[1], lwd = 1.5, type = "l", ylim = ylim, main = "Random loci", xlab = "", ylab = "Base relative frequency", xaxt = "n")
  lines(xplot, gc.ran.coords, col = mycols[2], lwd = 1.5)
  axis(side = 1, at = c(-1000, 0, 1000), labels = c("-1 kb", "0", "+1 kb"))
  axis(side = 2, at = pretty(c(at.coords, gc.coords, at.ran.coords, gc.ran.coords)))
  abline(v = 0, lty = 3)
  #abline(v = c(-100, 99), lty = 3)
  legend("right",
         legend = c("AT", "GC"),
         col = c(mycols[1], mycols[2]),
         text.col = c(mycols[1], mycols[2]),
         ncol = 1, cex = 0.8, lwd = 1.5, bty = "n")
  box(lwd = 1.5)
}

pdf(paste0(plotDir, "SPO11genomerangerPeaks_idr0.05_AT_GC_merge_base_comp_v100118.pdf"))
xplot <- seq(-1000, 999, by = 1)
mycols <- c("blue", "red")
genome.ylim <- c(min(genome.at.coords, genome.gc.coords, genome.at.ran.coords, genome.gc.ran.coords),
              max(genome.at.coords, genome.gc.coords, genome.at.ran.coords, genome.gc.ran.coords))
par(mfrow = c(2, 2))
par(mar = c(5, 4, 4, 4))
par(mgp = c(3, 0.75, 0))
mergeBaseCompPlot(xplot = xplot, at.coords = genome.at.coords, gc.coords = genome.gc.coords, at.ran.coords = genome.at.ran.coords, gc.ran.coords = genome.gc.ran.coords, ylim = genome.ylim, mycols = mycols)
dev.off()


sessionInfo()

