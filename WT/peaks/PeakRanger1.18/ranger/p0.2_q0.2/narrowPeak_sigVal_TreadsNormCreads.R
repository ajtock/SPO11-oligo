# Adjust "signalValue" (region treatment reads; column 7) by region control reads
# signalValue = region treatment reads/region control reads
# Use this signalValue with caution - unknown whether these read counts are normalised by library size by PeakRanger
# May be preferrable to use qValue (FDR) or pValue in downstream analyses, including IDR calculation

peakDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/SPO11_peaks/PeakRanger1.18/ranger/p0.2_q0.2/"

for(i in c(1, 3, 8)) {
  peaks <- read.table(file = paste0(peakDir, "WT_SPO11oligos_RPI", i, "_peaks_peakranger_ranger_p0.2_q0.2_Treads_Creads.narrowPeak"))
  colnames(peaks) <- c("chr", "start0based", "end", "name", "score", "strand",
                       "signalVal", "pVal", "qVal", "summit0based", "treads", "creads")
  peaks <- cbind(peaks[,1:6], peaks[,11]/peaks[,12], peaks[,8:10])
  colnames(peaks) <- c("chr", "start0based", "end", "name", "score", "strand",
                       "TreadsNormCreads", "pVal", "qVal", "summit0based")
  head(peaks)
  write.table(peaks, file = paste0(peakDir, "WT_SPO11oligos_RPI", i, "_peaks_peakranger_ranger_p0.2_q0.2_TreadsNormCreads.narrowPeak"),
              col.names = F, row.names = F, sep = "\t", quote = F)
}

