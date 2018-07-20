#!/bin/bash

######################
# Reformat bed files #
######################

for i in log2SPO11oligoNakedDNA
  do
  awk '{print $1, $2, $4}' ${i}_norm_allchrs_coverage_coord_tab.bed > ${i}_norm_allchrs_coverage_coord.bed
  sed 's/ /\t/g' ${i}_norm_allchrs_coverage_coord.bed > ${i}_norm_allchrs_coverage_coord_tab_cp.bed
  rm ${i}_norm_allchrs_coverage_coord.bed
done

