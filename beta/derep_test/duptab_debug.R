#!/usr/bin/env Rscript

# old output in wide format
# note that singletons have been removed in this!
duptab_old <- as.matrix(read.csv("duplicate_table.csv", row.names = 1))

# new file in long format
df_long_file <- "derep.map"
df_long <- read.table(df_long_file, stringsAsFactors = FALSE)

# exclude singletons
non_singletons <- names(which(sapply(split(df_long[,3], df_long[,1]), sum) > 1))
df_long <- df_long[df_long[,1] %in% non_singletons,]

# check that they're gone
sum(sapply(split(df_long[,3], df_long[,1]),sum) < 2)

samples_all <- union(rownames(duptab_old), derep_long[,2])

# reshape to wide format
duptab_new <- reshape(df_long, idvar = "V2", 
                      timevar = "V1", direction = "wide")

# move column 1 (sample names) to row names 
rownames(duptab_new) <- duptab_new[,1]

# remove column 1 and make into matrix
duptab_new <- as.matrix(duptab_new[,-1])


################################################################################
# correct column names
colnames(duptab_new) <- gsub("V3.", replacement = "", colnames(duptab_new))

# change NA's to 0's
duptab_new[is.na(duptab_new)] <- 0

# exclude singletons (i.e. columns with only 1 occurrence across ALL samples!)
# derep_wide <- derep_wide[,colSums(derep_wide) > 1]

duptab_list <- list(old = duptab_old, new = duptab_new)
lapply(duptab_list, function(x) x[1:5,1:5])
lapply(duptab_list, sum)
lapply(duptab_list, dim)

# check column names
if(!identical(colnames(duptab_old), colnames(duptab_new))){
	"uh oh, the sequence IDs are different in the two matrices:"
	setdiff(colnames(duptab_old), colnames(duptab_new))
} else { "looks fine" }

# check row names
identical(rownames(duptab_old), rownames(duptab_new))
setdiff(rownames(duptab_old), rownames(duptab_new))
absent_from_new <- rownames(duptab_old)[!rownames(duptab_old) %in% rownames(duptab_new)]
absent_from_old <- rownames(duptab_new)[!rownames(duptab_new) %in% rownames(duptab_old)]
if(sum(duptab_old[absent_from_new,]) > 0){
	"yikes, some samples are missing from the new version that are > 0"
} else { "some samples are missing from the nrew version, but they have no reads" }

if(sum(duptab_new[absent_from_old,]) > 0){
	"yikes, some samples are missing from the original version that are > 0"
} else { "some samples are missing from the original version, but they have no reads" }

# the reason a sample was missing from the original duplicate_table is unintentional:
# that sample had 1 read after demultiplexing and primer removal, and therefore was a singleton, AND the only read in that sample.

# OK, let's move on and just deal with samples that have reads
duptab_new <- duptab_new[rowSums(duptab_new) > 0,]
duptab_old <- duptab_old[rowSums(duptab_old) > 0,]

identical(duptab_old, duptab_new)

duptab_list <- list(old = duptab_old, new = duptab_new)
lapply(duptab_list, function(x) x[1:5,1:5])
lapply(duptab_list, sum)
lapply(duptab_list, dim)

identical(colnames(duptab_old), colnames(duptab_new))

identical(rownames(duptab_old), rownames(duptab_new))
setdiff(rownames(duptab_old), rownames(duptab_new))

duptab_new <- duptab_new[rownames(duptab_old),]
identical(rownames(duptab_new), rownames(duptab_old))

identical(duptab_old, duptab_new) 

identical(rowSums(duptab_old), rowSums(duptab_new))
identical(colSums(duptab_old), colSums(duptab_new))

library(vegan)
dissim <- sapply(1:nrow(derep_orig), function(x) 
       vegdist(rbind(sort(derep_orig[x,]), sort(derep_wide[x,]))))

if(all(dissim == 0)){
	"The samples are all identical, the order of sequence ids is just jumbled."
}

# Whew. The problem was just that a set of sequences with identical abundance are given arbitrary output IDs in the two methods.
# e.g. if the 2nd and 3rd most abundant sequences both occur 22 times, 
# in method one, the sequence itself could be called DUP_2, while in the second set it could be DUP_3
# method 1 fasta:
# >DUP_2;size=22
# AAAAAAAAAAAAAAA
# >DUP_3;size=22
# TTTTTTTTTTTTTTT

# method 2 fasta:
# >DUP_2;size=22
# TTTTTTTTTTTTTTT
# >DUP_3;size=22
# AAAAAAAAAAAAAAA

# this causes lots of things to go wrong. rest easy.
