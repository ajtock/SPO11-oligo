#!/bin/bash

## PeakRanger v1.18
## Use peakranger ranger to call narrow peaks in each SPO11-1-oligos replicate, using naked gDNA as an input control sample
## From the PeakRanger manual at http://ranger.sourceforge.net/manual1.18.html :

### The algorithm

## ranger uses a staged algorithm to discover enriched regions and the summits within them. In the first step, PeakRanger implements a FDR based adapative thresholding algorithm, which was originally proposed by PeakSeq. PeakRanger uses this thresholder to find regions with enriched reads that exceed expects. After that, PeakRanger searches for summits in these regions. The summit-search algorithm first looks for the location with largest number of reads. It then searchs for sub-summits with the sensitivity, the delta -r, specified by the user. Smaller -r will generate more summits.The coverage profiles are smoothed and padded before calling summits. The smoothing grade varies with -b. Higher smoothing bandwidth results less false summits at the cost of degraded summit accuracy .To measure the significance of the enriched regions, PeakRanger uses binormial distribution to model the relative enrichment of sample over control. A p value is generated as a result. Users can thus select highly significant peaks by using a smaller -p.

bamDir=/projects/ajt200/BAM_masters/SPO11-oligo/WT/SPO11_peaks

for i in 1 3 8
do
  peakranger ranger -d $bamDir/WT_SPO11-oligo_RPI${i}_k10_bt2_mapped_lowmiss_unique_both_sort_rmdup_new_sort.bam \
                    -c $bamDir/Col_trim50_R1_k10_bt2_mapped_lowmiss_unique_both_sort_rmdup.bam \
                    --format bam -o WT_SPO11oligos_RPI${i}_peaks_peakranger_ranger_p0.2_q0.2 \
                    -p 0.2 -q 0.2 -l 200 --pad -t 24 --verbose 
done

