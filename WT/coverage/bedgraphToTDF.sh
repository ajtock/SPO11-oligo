#######################################################################################
# Convert concatenated bedgraph files to TDF format for visualisation in IGV          #
#######################################################################################

for i in WT_SPO11-oligo_RPI1 WT_SPO11-oligo_RPI3 WT_SPO11-oligo_RPI8
  do
    igvtools toTDF ${i}_norm_allchrs_coverage.bedgraph ${i}_norm_allchrs_coverage.bedgraph.tdf /projects/ajt200/TAIR10/tair10.chrom.sizes
  done

#lib.names <- c("REC8_ChIP", "MSH4_ChIP")
#out.dir <- c("/projects/ajt200/ajt_Chris_Sasha_REC8_MSH4/coverage_xiaohui/merged1/")
#bed.files <- list(paste0(out.dir, lib.names[1], "_norm_allchrs_coverage_merged1.bedgraph"),
#                  paste0(out.dir, lib.names[2], "_norm_allchrs_coverage_merged1.bedgraph"))
#print(bed.files)

#system("igvtools toTDF bed.files[[1]] bed.files[[1]].tdf /projects/ajt200/tair10.chrom.sizes")
#system("igvtools toTDF bed.files[[2]] bed.files[[2]].tdf /projects/ajt200/tair10.chrom.sizes")

