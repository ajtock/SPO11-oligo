#######################################################################################
# Concatenate chromosome bedgraph files and add new first line "track type=bedGraph"  #
#######################################################################################

lib.names <- c("WT_SPO11-oligo_RPI1", "WT_SPO11-oligo_RPI3", "WT_SPO11-oligo_RPI8")
chr.index <- c(1:5)
out.dir <- c("/projects/ajt200/BAM_masters/SPO11-oligo/coverage/")
bed.files <- list(paste0(out.dir, lib.names[1], "_norm_Chr", chr.index, "_coverage.bedgraph"),
                  paste0(out.dir, lib.names[2], "_norm_Chr", chr.index, "_coverage.bedgraph"),
                  paste0(out.dir, lib.names[3], "_norm_Chr", chr.index, "_coverage.bedgraph"))
print(bed.files)

#WT_SPO11-oligo_RPI1
these.files <- bed.files[[1]]
myfiles1 <- read.table(these.files[1], header=F)
myfiles2 <- read.table(these.files[2], header=F)
myfiles3 <- read.table(these.files[3], header=F)
myfiles4 <- read.table(these.files[4], header=F)
myfiles5 <- read.table(these.files[5], header=F)

cat <- rbind(myfiles1, myfiles2, myfiles3, myfiles4, myfiles5)
cat("track type=bedGraph\n", file=paste0(out.dir, lib.names[1], "_norm_allchrs_coverage.bedgraph"))
write.table(cat, file=paste0(out.dir, lib.names[1], "_norm_allchrs_coverage.bedgraph"), append=T, row.names=F, col.names=F, quote=F)

#WT_SPO11-oligo_RPI3
these.files <- bed.files[[2]]
myfiles1 <- read.table(these.files[1], header=F)
myfiles2 <- read.table(these.files[2], header=F)
myfiles3 <- read.table(these.files[3], header=F)
myfiles4 <- read.table(these.files[4], header=F)
myfiles5 <- read.table(these.files[5], header=F)

cat <- rbind(myfiles1, myfiles2, myfiles3, myfiles4, myfiles5)
cat("track type=bedGraph\n", file=paste0(out.dir, lib.names[2], "_norm_allchrs_coverage.bedgraph"))
write.table(cat, file=paste0(out.dir, lib.names[2], "_norm_allchrs_coverage.bedgraph"), append=T, row.names=F, col.names=F, quote=F)

#WT_SPO11-oligo_RPI8
these.files <- bed.files[[3]]
myfiles1 <- read.table(these.files[1], header=F)
myfiles2 <- read.table(these.files[2], header=F)
myfiles3 <- read.table(these.files[3], header=F)
myfiles4 <- read.table(these.files[4], header=F)
myfiles5 <- read.table(these.files[5], header=F)

cat <- rbind(myfiles1, myfiles2, myfiles3, myfiles4, myfiles5)
cat("track type=bedGraph\n", file=paste0(out.dir, lib.names[3], "_norm_allchrs_coverage.bedgraph"))
write.table(cat, file=paste0(out.dir, lib.names[3], "_norm_allchrs_coverage.bedgraph"), append=T, row.names=F, col.names=F, quote=F)

