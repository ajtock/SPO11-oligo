####################################################################
# Calculate log2 ratio of ChIP and input coverage values           #
# and write to 1-based BED-like file suitable for EnrichedHeatmap  #
####################################################################

library(doParallel)
registerDoParallel(cores = 2)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

inDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/"
outDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/"
inNames <- c("WT_SPO11-oligo_RPI3_nakedDNA_R1", "WT_SPO11-oligo_RPI8_nakedDNA_R1")
outNames <- c("log2wtSPO11oligoRPI3NakedDNA", "log2wtSPO11oligoRPI8NakedDNA")

foreach(i = 1:length(inNames), .combine = 'c') %dopar% {
  ChIP_input <- read.table(file = paste0(inDir, inNames[i], "_norm_allchrs_coverage_coord_tab.bed"))
  print(head(ChIP_input))
  ChIP_input_offset <- cbind(ChIP_input, ChIP_input$V4+1, ChIP_input$V5+1)  
  ChIP_input_offset <- ChIP_input_offset[,-4:-5]
  colnames(ChIP_input_offset) <- c("V1", "V2", "V3", "V4", "V5")
  print(head(ChIP_input_offset))
  norm <- log2(ChIP_input_offset$V4/ChIP_input_offset$V5)
  norm <- (norm-mean(norm, na.rm = T))/sd(norm, na.rm = T)
  #norm[is.na(norm)] <- 0
  #norm[is.infinite(norm)] <- 0
  ChIPnormInput <- cbind(ChIP_input_offset, norm)
  log2ChIPinput <- ChIPnormInput[,-4:-5] 
  write.table(log2ChIPinput, file = paste0(outDir, outNames[i], "_norm_allchrs_coverage_coord_tab.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  rm(ChIPnormInput)
}
gc()

