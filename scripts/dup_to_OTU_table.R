#!/usr/bin/env Rscript

# Map counts of unique sequences to corresponding OTUs 
# (previously named collapse_dup_by_OTU)

################################################################################
# Arguments: 
#-------------------------------------------------------------------------------
# 1: path to file mapping unique sequences ("dup") to samples, with 3 columns: 
#   (col 1) unique sequence id, (col 2) sample id, (col 3) count
# 
# 2: dup to otu table; two columns with headers "Query" and "Match"
#   each row contains sequence headers of query and corresponding match from 
#   fasta files of OTU clustering process
# 
# 3: otu table path (to be written)
################################################################################

arguments <- commandArgs(TRUE)

#-------------------------------------------------------------------------------
library(gtools) # mixedsort
library(data.table)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Read in unique sequence to sample file
dups <- fread(arguments[1], header = FALSE, stringsAsFactors = FALSE)
colnames(dups) <- c("uniq", "sample", "count")

# as of NextSeq run in July 2015, dereplication had to be majorly overhauled, resulting in transposed duplicate table.
# dups <- t(dups)


# Read in dups to OTUs files
dups_to_OTUs <- fread(arguments[2], header=TRUE, stringsAsFactor=FALSE)
colnames(dups_to_OTUs) <- c("uniq", "otu")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# merge the two tables based on the unique sequence ID; 
dt1 <- merge(dups, dups_to_OTUs)

# sum the counts by sample and otu
otu_map <- dt1[,list(count = sum(count)), by = list(otu, sample)]

# OTU_table <- aggregate(dups, list(OTUs), FUN = sum)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# tests
# any(duplicated(otu_map[otu == "DUP_1",sample]))
# any(duplicated(dt1[otu == "DUP_1",sample]))
#-------------------------------------------------------------------------------

# sort by otu name
otu_map <- otu_map[mixedorder(otu)]

# write output to tab separated "map" file
fwrite(x = otu_map, file = arguments[3], 
  quote = FALSE, sep = "\t", col.names = FALSE)
