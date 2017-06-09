#!/usr/bin/env Rscript

# convert a table of counts of sequences in a set of samples from a 
# "mapfile" (long format) to "otu table" (wide format), where 
# rows = samples and columns = type of thing being counted

################################################################################
# ARGUMENTS:
# 1. infile
# 2. path for output file (to be written)
################################################################################

arguments <- commandArgs(TRUE)
infile <- arguments[1]
outfile <- arguments[2]

#-------------------------------------------------------------------------------
library(data.table) # fread, fwrite
library(gtools) # mixedsort/mixedorder
#-------------------------------------------------------------------------------


map.l <- fread(input = infile, stringsAsFactors = FALSE)
names(map.l) <- c("seq", "sample", "count")

map.w <- dcast.data.table(data = map.l, sample ~ seq, value.var = "count")

# compare times
# system.time(map.w <- dcast.data.table(data = map, sample ~ seq, value.var = "count"))
# system.time(map.w.mat <- as.matrix(map.w))

dt_na2zero <- function(DT) {
  # set NAs to zero in a data table
  for (i in seq_len(ncol(DT)))
    set(DT, which(is.na(DT[[i]])), i, 0)
}

dt_na2zero(map.w)

#-------------------------------------------------------------------------------
# order columns by otu name
map.w <- map.w[,c("sample", mixedsort(colnames(map.w)[-1])), with = FALSE]

fwrite(x = map.w, file = outfile, sep = ",")