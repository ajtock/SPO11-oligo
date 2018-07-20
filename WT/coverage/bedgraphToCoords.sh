#######################################################################################
# Convert bedgraph files to bed files containing coordinates and normalised coverage  #
#######################################################################################

for i in WT_SPO11-oligo_RPI1 WT_SPO11-oligo_RPI3 WT_SPO11-oligo_RPI8
  do
    awk '{print $1, $3, $3, $4}' ${i}_norm_allchrs_coverage.bedgraph > ${i}_norm_allchrs_coverage_coord.bed
    sed 's/ /\t/g' ${i}_norm_allchrs_coverage_coord.bed > ${i}_norm_allchrs_coverage_coord_tab.bed
  done

