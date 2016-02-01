# Annotate the OTU table using the output from MEGAN

setwd("/Users/threeprime/Desktop")

library(gtools)

OTU_table <- read.csv("/Users/threeprime/temp_big_upload/Analysis_20150512_0457/all_lib/OTU_table.csv", row.names = 1)

megan_output <- read.csv("/Users/threeprime/temp_big_upload/Analysis_20150512_0457/all_lib/megan/mod/meganout_Phylum_mod.csv", header = FALSE, stringsAsFactors = FALSE)

OTUs_not_in_megan <- paste(setdiff(rownames(OTU_table), megan_output[,1]), collapse = " ")

OTU_to_taxon <- megan_output[match(rownames(OTU_table), megan_output[,1]),3]

OTU_to_taxon[which(is.na(OTU_to_taxon))] <- OTUs_not_in_megan #"NoHits"

# for each column (sample), sum the rows (taxa) that belong to the same OTU
Taxon_table <- aggregate(OTU_table, list(OTU_to_taxon), FUN = sum)

# Make the rownames the values stored in the new first colum
rownames(Taxon_table) <- Taxon_table[,1]

# remove that column
Taxon_table <- Taxon_table[,-1]

# sort by rowname
OTU_table <- OTU_table[mixedsort(rownames(OTU_table)),]

# write output to CSV file
write.csv(Taxon_table, "taxon_table_phylum.csv")

