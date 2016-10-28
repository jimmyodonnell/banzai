#!/usr/bin/env Rscript

arguments   <- commandArgs(TRUE)

df_long_file <- arguments[1]  # path to the long format file
df_wide_file <- arguments[2]  # path to the wide format file (to write)
rm_single    <- arguments[3]  # from remove_singletons=[YES|NO]

df_long <- read.table(df_long_file, stringsAsFactors = FALSE)

# exclude singletons
single_msg <- "not"
if( rm_single == "YES"){
  non_singletons <- names(which(sapply(split(df_long[,3], df_long[,1]), sum) > 1))
  df_long <- df_long[df_long[,1] %in% non_singletons,]
  single_msg <- ""
}

# check that they're gone
# sum(sapply(split(df_long[,3], df_long[,1]),sum) < 2)

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

write.csv(x = duptab_new, file = df_wide_file, quote = FALSE)

duptab_dim <- dim(duptab_new)
exit_msg <- list(
  c("Contingency table written to file:", df_wide_file),
  paste("It contains counts of", duptab_dim[2], "unique sequences in", 
  duptab_dim[1], "samples."),
  c("Samples with 0 reads do not appear in this table."),
  paste("Singletons were", single_msg, "removed.")
  )

invisible(lapply(exit_msg, writeLines))
