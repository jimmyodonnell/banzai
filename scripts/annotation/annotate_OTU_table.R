#!/usr/bin/env Rscript

# Annotate the OTU table using the output from MEGAN

args <- commandArgs(TRUE)


if(length(args) != 2 ){
  stop("this script requires two arguments: \nthe paths to an OTU table and a meganout_mod.csv")
}

library(gtools)

OTU_table <- read.csv(args[1], row.names = 1)

megan_output <- read.csv(args[2], header = FALSE, stringsAsFactors = FALSE)

if (length(grep("size=", megan_output[,2])) == 0 ) {
  stop("the megan input does not contain expected fields")
}
if (ncol(megan_output) != 3 ) {
  stop("the megan input does not contain three columns")
}

# in which row of the MEGAN output does each OTU occur?
OTU_in_megan <- match(rownames(OTU_table), megan_output[,1])

name_for_OTU <- megan_output[OTU_in_megan,3]

name_for_OTU[is.na(name_for_OTU)] <- "Not annotated"

# for each column (sample), sum the rows (taxa) that belong to the same OTU
Taxon_table <- aggregate(OTU_table, list(name_for_OTU), FUN = sum)

# Make the rownames the values stored in the new first colum
rownames(Taxon_table) <- Taxon_table[,1]

# remove that column
Taxon_table <- Taxon_table[,-1]

# sort by abundance
Taxon_table <- Taxon_table[order(rowSums(Taxon_table), decreasing = TRUE),]

# write output to CSV file
write.csv(Taxon_table, "taxon_table.csv")

