# Adjust "signalValue" (region treatment reads; column 7) by region width
# signalValue = region treatment reads/region width
# Use this signalValue with caution - unknown whether these read counts are normalised by library size by PeakRanger
# May be preferrable to use qValue (FDR) or pValue in downstream analyses, including IDR calculation

peakDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/SPO11_peaks/PeakRanger1.18/ranger/p0.2_q0.2/"

for(i in c(1, 3, 8)) {
  peaks <- read.table(file = paste0(peakDir, "WT_SPO11oligos_RPI", i, "_peaks_peakranger_ranger_p0.2_q0.2.narrowPeak"))
  peak_widths <- (peaks[,3]-peaks[,2])+1
  TreadsNormWidth <- peaks[,7]/peak_widths
  peaks_TreadsNormWidth <- cbind(peaks[,1:6], TreadsNormWidth, peaks[,8:10])
  colnames(peaks_TreadsNormWidth) <- c("chr", "start0based", "end", "name", "score", "strand",
                                       "TreadsNormWidth", "pVal", "qVal", "summit0based")
  head(peaks_TreadsNormWidth)
  write.table(peaks_TreadsNormWidth, file = paste0(peakDir, "WT_SPO11oligos_RPI", i, "_peaks_peakranger_ranger_p0.2_q0.2_TreadsNormWidth.narrowPeak"),
              col.names = F, row.names = F, sep = "\t", quote = F)
}

