# translate duplicates to OTUs: was named collapse_dup_by_OTU

# Set working directory
# setwd("/Users/threeprime/Documents/GoogleDrive/Kelly_Lab_Big/Illumina_Data_Analyzed/16S/run_20150401/Analysis_20150403_1906/all_lib")
# arguments: 
# 1: duplicate table
# 2: dup to otu table
# 3: otu table path
# 4: concatenated directory (obsolete?))

# automated script:
arguments <- commandArgs(TRUE)
setwd(arguments[4])

# load gtools to use function mixedsort
library(gtools)

# Read in duplicates files
# read.csv("dups.csv", row.names = 1)
dups <- read.csv(arguments[1], row.names = 1)

# Read in dups to OTUs files
dups_to_OTUs <- read.csv(arguments[2], header=TRUE, stringsAsFactor=FALSE)

OTUs <- dups_to_OTUs$Match[
  match(
    rownames(dups),
    dups_to_OTUs$Query
  )
]

# for each column (sample), sum the rows (duplicates) that belong to the same OTU
OTU_table <- aggregate(dups, list(OTUs), FUN = sum)

# Make the rownames the values stored in the new first colum
rownames(OTU_table) <- OTU_table[,1]

# remove that column
OTU_table <- OTU_table[,-1]

# sort by rowname
OTU_table <- OTU_table[mixedsort(rownames(OTU_table)),]

# write output to CSV file
write.csv(x = OTU_table, file = arguments[3])





# Most of this was an artifact of the weirdness of the old usearch output (uc format)
# unique(c(dup_to_OTU[,1], dup_to_OTU[,2]))
# dups_to_collapse <- split(dup_to_OTU[,1], dup_to_OTU[,2])
# dups_to_collapse <- lapply(dups_to_collapse, as.character)
# dups_to_collapse <- mapply(c, as.list(names(dups_to_collapse)), dups_to_collapse)

# dups_to_collapse <- sapply(dups_to_collapse, unique)
# dups_to_collapse <- sapply(dups_to_collapse, sort)
# 
# 
# no_clusters <- dups[!rownames(dups) %in% unique(unlist(dups_to_collapse)),]
# 
# consolidated_dups <- list()
# for(i in 1:length(dups_to_collapse)){
# 	consolidated_dups[[i]] <- colSums(dups[dups_to_collapse[[i]],])
# }
# consolidated_dups <- do.call(rbind, consolidated_dups)
# rownames(consolidated_dups) <- sapply(dups_to_collapse, function(x) {x[[1]]})
# 
# sort(as.numeric(gsub("DUP_", "", rownames(ALL_CLUSTERS))))
# ALL_CLUSTERS <- rbind(as.data.frame(consolidated_dups), no_clusters)
# order(rownames(ALL_CLUSTERS))
# write.csv(ALL_CLUSTERS, "all_clusters.csv")
# ALL_CLUSTERS <- read.csv("")

# confirm there are no duplicate otus:
# which(duplicated(rownames(ALL_CLUSTERS)))


# this one is weird, among many others
# dups_to_collapse[2416]

# write.table(dups_to_collapse, "dups_to_collapse.txt")
# dups_to_collapse_vec <- unlist(dups_to_collapse)
# 
# dupped <- dups_to_collapse_vec[duplicated(dups_to_collapse_vec)]


# dups_removed <- dup_to_OTU[-which(dup_to_OTU[,2] %in% dupped),]
# dups_removed_tmp <- lapply(split(dups_removed[,1], dups_removed[,2]), as.character)
# dups_removed <- mapply(c, as.list(names(dups_removed_tmp)), dups_removed_tmp)
# which(duplicated(unlist(dups_removed)))


# X <- stack(setNames(dups_to_collapse, seq_along(dups_to_collapse)))
# TAB <- table(X)
# TAB.mat <- as.matrix(TAB)
# dup_rows <- TAB.mat[which(rowSums(TAB.mat) > 1),]
# dup_cols <- dup_rows[,which(colSums(dup_rows) > 0)]
# # dup_cols[,which(colSums(dup_cols) > 1)]
# identical(sort(rownames(dup_cols)),sort(dupped))
# # edit(dup_cols)
# dups_to_collapse[as.numeric(colnames(dup_cols))]
# 
# TRASH:
# Read in files of chimaeras vs "not chimaeras"
# chimaeras <- read.table("chimaeras.txt")
# not_chimaeras <- read.table("not_chimaeras.txt")
# tail(chimaeras)

# dups_chimaeras <- dups[as.character(chimaeras[,1]),]
# dups_no_chimaeras <- dups[as.character(not_chimaeras[,1]),]


# WHAT THE FUCK IS GOING ON HERE???

# sapply(consolidated_otus, duplicated)

# consolidated_otus <- list()
# for(i in 1:nrow(dup_cols)){
	# consolidated_otus[[i]] <- Reduce(union, dups_to_collapse[as.numeric(names(which(dup_cols[i,] > 0)))])
# }
# consol <- stack(setNames(consolidated_otus, seq_along(consolidated_otus)))
# consol.mat <- as.matrix(table(consol))
# colSums(consol.mat)

# duplicated(unlist(consolidated_otus))

# consolidated_2 <- list()
# for(i in 1:nrow(consol.mat)){
	# consolidated_2[[i]] <- Reduce(union, consolidated_otus[as.numeric(names(which(consol.mat[i,] > 0)))])
# }
# consol2 <- stack(setNames(consolidated_2, seq_along(consolidated_2)))
# consol2.mat <- as.matrix(table(consol2))
# colSums(consol2.mat)

# consolidated_3 <- list()
# for(i in 1:nrow(consol.mat)){
	# consolidated_3[[i]] <- Reduce(union, consolidated_2[as.numeric(names(which(consol2.mat[i,] > 0)))])
# }

# Reduce(union, dups_to_collapse[as.numeric(names(which(dup_cols[35,] > 0)))])
# Reduce(intersect, dups_to_collapse)