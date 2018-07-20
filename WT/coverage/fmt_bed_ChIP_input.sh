#!/bin/bash

# Append "input" file columns onto "ChIP" file columns and reformat bed-like file to include coverage for both and to remove redundant columns

for i in WT_SPO11-oligo_meanAllReps
do
  paste -d "\t" ${i}_norm_allchrs_coverage_coord_tab.bed /projects/ajt200/BAM_masters/WT_nakedDNA/single_trim50/R1/coverage/WT_nakedDNA_R1_norm_allchrs_coverage_coord_tab.bed > WT_SPO11-oligo_meanAllReps_nakedDNA_R1_norm_allchrs_coverage_coord_tab_tmp.bed
done

for i in WT_SPO11-oligo_meanAllReps_nakedDNA_R1
do
  awk '{print $1, $2, $3, $4, $8}' ${i}_norm_allchrs_coverage_coord_tab_tmp.bed > ${i}_norm_allchrs_coverage_coord.bed
  rm ${i}_norm_allchrs_coverage_coord_tab_tmp.bed
  sed 's/ /\t/g' ${i}_norm_allchrs_coverage_coord.bed > ${i}_norm_allchrs_coverage_coord_tab.bed
  rm ${i}_norm_allchrs_coverage_coord.bed
done

