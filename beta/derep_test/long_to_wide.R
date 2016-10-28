setwd("~/banzai_out_20161028_0508/all_lib")

# note that singletons have been removed in this!
derep_orig <- read.csv("duplicate_table.csv", row.names = 1)
derep_orig[1:5,1:5]

derep_long <- read.table("derep.map", stringsAsFactors = FALSE)
derep_long <- derep_long[derep_long[,3] > 1,]

derep_wide <- reshape(derep_long, idvar = "V2", timevar = "V1", direction = "wide")
derep_wide[1:5,1:5]
rownames(derep_wide) <- derep_wide[,1]
derep_wide <- derep_wide[,-1]

derep_wide[1:5,1:5]
colnames(derep_wide) <- gsub("V3.", replacement = "", colnames(derep_wide))
derep_wide[1:5,1:5]
derep_wide[is.na(derep_wide)] <- 0
derep_wide[1:5,1:5]

dim(derep_wide)
dim(derep_orig)
boxplot(sapply(list(derep_orig, derep_wide), rowSums))
boxplot(sapply(list(derep_orig, derep_wide), colSums))

rownames(derep_orig) %in% rownames(derep_wide)
rownames(derep_wide) %in% rownames(derep_orig)

setdiff(colnames(derep_orig), colnames(derep_wide))
samples_missing <- setdiff(rownames(derep_orig), rownames(derep_wide))
rowSums(derep_orig[samples_missing,])
# this explaing the difference in rows!

derep_orig <- derep_orig[rowSums(derep_orig) > 1,]
dim(derep_orig)
dim(derep_wide)

# can't explain the mising columns yet...

plot(
sort(as.matrix(derep_orig), decreasing = TRUE)[1:100],
sort(as.matrix(derep_wide), decreasing = TRUE)[1:100]
)

sort(as.matrix(derep_orig), decreasing = TRUE)[1:100],
sort(as.matrix(derep_wide), decreasing = TRUE)[1:100]

table(as.matrix(derep_orig))
table(as.matrix(derep_wide))

hist(as.matrix(derep_orig), breaks = 31, col = hsv(1,1,1,0.2))
hist(as.matrix(derep_wide), breaks = 31, col = hsv(0.6,1,1,0.2), add = TRUE)
# lots more zeros in the original

plot(sort(colSums(derep_orig), decreasing = TRUE)[1:200], log = "y")
points(sort(colSums(derep_wide), decreasing = TRUE), log = "y", col = "red")


derep_wide <- derep_wide[rownames(derep_orig), colnames(derep_orig)]
class(derep_orig)
class(derep_wide)
identical(derep_orig, derep_wide)
dim(derep_wide)
dim(derep_orig)