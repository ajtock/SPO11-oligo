###########################################################################
# Plot heatmaps of log2-transformed coverage levels at and flanking genes #
###########################################################################

###########################################################################################################################################
# Plot SPO11-1, MNase and H3K4me3 heatmaps with rows (features; e.g., genes) ordered according to TSS–TTS coverage levels of each dataset #
# to enable comparison of samples                                                                                                         #
###########################################################################################################################################

library(parallel)
library(doParallel)
library(EnrichedHeatmap)
#library(genomation)
library(circlize)
library(RColorBrewer)

plotDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/geneHeatmaps/analysis_04_SPO11_nuc_H3K4me3_sixcols/plots/"
matDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/geneHeatmaps/analysis_01_SPO11_nuc_H3K4me3/matrices/"
rowOrderDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/geneHeatmaps/analysis_01_SPO11_nuc_H3K4me3/row_order/"

# import genes as GRanges object
#representive_genes_uniq <- system("ls /projects/ajt200/TAIR10/representative_genes/representative_genes_uniq_fmt_strand.txt", intern=T)
#genes <- readGeneric(representive_genes_uniq, header = T, strand = 4, meta.col = list(gene_model = 5))
#print("***********genes***********")
#print(genes)
#tss <- promoters(genes, upstream = 0, downstream = 1)
#print(tss)

#WT_SPO11_oligos_RPI1 <- system("ls /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2WTSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab.bed", intern = T)
#WT_MNase <- system("ls /projects/ajt200/BAM_masters/nucleosomes/WT/coverage/log2ChIPinput/log2WTnucNakedDNA_norm_allchrs_coverage_coord_tab.bed", intern = T)
#WT_H3K4me3_ChIP12 <- system("ls /projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/WT_H3K4me3_ChIP12_log2ChIPinput_norm_allchrs_coverage_coord_tab.bed", intern = T)

## create lists/vectors of path and library names
#libPaths <- list(WT_SPO11_oligos_RPI1, WT_MNase, WT_H3K4me3_ChIP12)
libNames <- c("WT_SPO11_oligos_RPI1", "WT_MNase", "WT_H3K4me3_ChIP12") 
orderNames <- c("WT_SPO11_oligos_RPI1", "WT_MNase", "WT_H3K4me3_ChIP12")

# import coverage files as GRanges objects and assign to library names
#grTmp <- mclapply(seq_along(libPaths), function(x) {
#  readGeneric(libPaths[[x]], meta.col = list(coverage = 4))
#}, mc.cores = 3, mc.preschedule = F)
#for(i in 1:length(grTmp)) {
#  seqlevels(grTmp[[i]]) <- sub("Chr", "", seqlevels(grTmp[[i]]))
#  assign(paste0(libNames[i]), grTmp[[i]])
#}

#print("***WT_SPO11_oligos_RPI1***")
#print(WT_SPO11_oligos_RPI1)
#print("***WT_MNase***")
#print(WT_MNase)
#print("***WT_H3K4me3_ChIP12***")
#print(WT_H3K4me3_ChIP12)

#grl <- GRangesList("WT_SPO11_oligos_RPI1" = WT_SPO11_oligos_RPI1, "WT_MNase" = WT_MNase, "WT_H3K4me3_ChIP12" = WT_H3K4me3_ChIP12)
#orderGrl <- GRangesList("WT_SPO11_oligos_RPI1" = WT_SPO11_oligos_RPI1, "WT_MNase" = WT_MNase, "WT_H3K4me3_ChIP12" = WT_H3K4me3_ChIP12)

## define window size for calculation of average coverage values within windows, and define lower and upper coverage value quantiles to be trimmed
#w <- 20
#trim <- c(0, 0.01)

# function to create coverage matrices and heatmaps for genes
## and flanking regions for each dataset, and to extract row order
#geneHeatmap <- function(signal, target, x) {
#  set.seed(7539)
#  mat1_smoothed <- normalizeToMatrix(signal, target, value_column = "coverage",
#                                     extend = 0, mean_mode = "w0", w = w, trim = trim,
#                                     empty_value = 0, smooth = TRUE,
#                                     include_target = TRUE, target_ratio = 1)
#  print(mat1_smoothed)
#  print(length(mat1_smoothed))
#  save(mat1_smoothed, file = outFile[[x]])
#
#  rich8to6equal <- c("#000041", "#0000CB", "#0081FF", "#FDEE02", "#FFAB00", "#FF3300")
#  col_fun1 <- colorRamp2(quantile(mat1_smoothed, c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)), rich8to6equal)
#  ht1 <- print(EnrichedHeatmap(mat1_smoothed, col = col_fun1, name = "Coverage",
#                               column_title = paste0(orderNames[x], " around genes"),
#                               axis_name = c("TSS", "TTS"), border = FALSE))
#  ht1 <- draw(ht1)
#  row_order_ht1 <- row_order(ht1)
#  save(row_order_ht1, file = rowOrderFiles[[x]])
#}

## run geneHeatmap() on each dataset to obtain heatmap row orders
### for subsequent ordering of genes for each dataset by each dataset
#outFile <- lapply(seq_along(orderNames), function(x) {
#             paste0(matDir, orderNames[x],
#                    "_norm_cov_genes_mat1_smoothed_target.RData")
#           })
#rowOrderFiles <- lapply(seq_along(orderNames), function(x) {
#                   paste0(rowOrderDir, orderNames[x],
#                          "_norm_cov_genes_mat1_smoothed_row_order.RData")
#                 })

#mclapply(seq_along(orderGrl), function(x) {
#  geneHeatmap(orderGrl[[x]], genes, x)
#}, mc.cores = 3, mc.preschedule = F)

# function to create coverage matrices and heatmaps for genes
## and flanking regions for each dataset, with heatmaps organised by above-defined row order
geneHeatmapOrder <- function(h, i) {
  #set.seed(7539)
  #mat2_smoothed <- normalizeToMatrix(signal, target, value_column = "coverage",
  #                                   extend = 2000, mean_mode = "w0", w = w, trim = trim,
  #                                   empty_value = 0, smooth = TRUE,
  #                                   include_target = TRUE, target_ratio = 0.3333333)
  #print(mat2_smoothed)
  #print(length(mat2_smoothed))
  #save(mat2_smoothed, file = outFile[[h]][[i]])
  load(file = outFile[[h]][[i]])
  
  rich8to6equal <- c("#000041", "#0000CB", "#0081FF", "#FDEE02", "#FFAB00", "#FF3300")
  col_fun1 <- colorRamp2(quantile(mat2_smoothed, c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)), rich8to6equal)
  load(rowOrderFiles[[h]])
  pdf(plotFiles[[h]][[i]])
  ht1 <- print(EnrichedHeatmap(mat2_smoothed, row_order = row_order_ht1[[1]],
                               col = col_fun1, name = "Coverage",
                               column_title = paste0(libNames[i], " around genes, by ", orderNames[h]),
                               axis_name = c("-2 kb", "TSS", "TTS", "2 kb"), border=FALSE))
  dev.off()
  ht1 <- draw(ht1)
  row_order_ht1 <- row_order(ht1)
  #save(row_order_ht1, file = rowOrderFilesConstr[[h]][[i]])
} 

# run geneHeatmapOrder() on each dataset, with genes order by each dataset for each dataset
outFile <- lapply(seq_along(orderNames), function(h) {
             lapply(seq_along(libNames), function(i) {
               paste0(matDir, libNames[i], "_by_", orderNames[h],
                      "_norm_cov_genes_mat2_smoothed_target_and_flank.RData")
             })
           })
rowOrderFiles <- lapply(seq_along(orderNames), function(h) {
                   paste0(rowOrderDir, orderNames[h],
                          "_norm_cov_genes_mat1_smoothed_row_order.RData")
                 })
#rowOrderFilesConstr <- lapply(seq_along(orderNames), function(h) {
#                         lapply(seq_along(grl), function(i) {
#                            paste0(rowOrderDir, libNames[i], "_by_", orderNames[h],
#                                   "_norm_cov_genes_mat2_smoothed_row_order.RData")
#                         })
#                       })
plotFiles <- lapply(seq_along(orderNames), function(h) {
               lapply(seq_along(libNames), function(i) {
                 paste0(plotDir, libNames[i], "_by_", orderNames[h],
                        "_norm_cov_genes_heatmap_trim0.00_0.01_smoothed.pdf")
               })
             })
             
lapply(seq_along(orderNames), function(h) {
  mclapply(seq_along(libNames), function(i) {
    geneHeatmapOrder(h, i)
  }, mc.cores = 3, mc.preschedule = F)
})


sessionInfo() 
