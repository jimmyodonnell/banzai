# install.packages("taxize")

# todo: replace do.call(rbind, ...) with data.table::rbindlist(...)
# from http://stackoverflow.com/questions/14972998/how-to-avoid-renaming-of-rows-when-using-rbind-inside-do-call
library(taxize)

blast_results_file_path  <- "/Users/threeprime/Documents/GoogleDrive/Kelly_Lab/Projects/Lemonade/Data/blast_20151125_1530/blast_results_all.txt"
blast_results <- read.table(
	file = blast_results_file_path, 
	sep = "\t", 
	stringsAsFactors = FALSE, 
	quote = NULL, 
	comment = ''
	)

# table columns order: output_format="6 qseqid sallseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"

query_col=1
evalue_col=11
bitscore_col=12
title_col=13
gi_col=2
taxid_col="taxid_all" # change this to a column number if it was returned in the blast output

# get gi numbers
gi_all <- do.call(c, lapply(strsplit(blast_results[,gi_col], split = "|", fixed = TRUE), "[", 2))
blast_results <- cbind.data.frame(blast_results, gi_all, stringsAsFactors = FALSE)
gi_unique <- as.character(unique(gi_all))

time_start <- Sys.time()
# get taxid from gi number. This could be avoided by having the taxid given back by blastn
taxids <- vector(mode = "character")
for(i in 1:length(gi_unique)){
	taxids[i] <- genbank2uid(id = gi_unique[i])[1]
}


gi_taxid <- data.frame(
	gi = gi_unique, 
	taxid = taxids, 
	stringsAsFactors = FALSE
	)
	

write.table(
	x = gi_taxid, 
	file = "gi_taxid_20160202.txt", 
	quote = FALSE, 
	row.names = FALSE
)

# put taxon ids onto blast results
taxid_all <- gi_taxid$taxid[match(gi_all, gi_taxid$gi)]
blast_results <- cbind.data.frame(blast_results, taxid_all, stringsAsFactors = FALSE)


# get taxonomic hierarchy from taxon ids
taxid_uniq <- unique(taxids)
classifications <- classification(x = taxid_uniq, db = "ncbi")
save(classifications, file = "classifications20160202.RData")

# extract the names only (exclude rank name, e.g. "Genus")
names_only <- lapply(classifications, "[[", 1)

# what is the group that is common to all results?
common_ancestor <- Reduce(intersect, names_only)

# Reduce(intersect, names_only[c("239049", "219410")])
# taxids_to_collapse <- c("239049", "219410")
# Reduce(intersect, names_only[taxids_to_collapse])

blast_queries <- split(blast_results, blast_results[, query_col])

identical(unique(blast_results[, query_col]), names(blast_queries))

besthits <- function(x){
	x[x[, evalue_col] == min(x[, evalue_col]),]
}

##############################################################################################
# Holy macaroni, that's it! Come back and clean this up!
least_common_ancestor <- list()
for(i in 1:length(blast_queries)){
	besthit_gis <- as.character(besthits(blast_queries[[i]])[,"gi_all"])
	besthit_taxids <- unique(gi_taxid[match(besthit_gis, gi_taxid[,"gi"]),"taxid"])
	least_common_ancestor[[i]] <- Reduce(intersect, names_only[besthit_taxids])
}
# Holy macaroni, that's it! Come back and clean this up!
##############################################################################################
time_end <- Sys.time()
# 1.771811 hours for length(least_common_ancestor) == 1601, mostly over network getting taxid

# extract only the rows corresponding to the lowest e-value (best hit)
blast_queries_best <- lapply(X = blast_queries, FUN = besthits)

# assess the number of taxonomic ties per e-value
blast_queries[[1]][ , "taxid_all"]

taxid_per_query <- split(blast_results[ , taxid_col], blast_results[, query_col])
taxid_per_query_best <- lapply(blast_queries_best, "[[", taxid_col)

taxon_hit_index <- function(x)
{
	# this function takes a character vector of taxa (names or id numbers)
	# and calculates the ratio of unique taxa per hit
	# varies between 0 and 1, lower is better.
	length(unique(x))/length(x)
}

taxon_hit_index(blast_queries[[1]][ , "taxid_all"])

taxon_hit_all <- sapply(X = taxid_per_query, FUN = taxon_hit_index)
taxon_hit_best <- sapply(X = taxid_per_query_best, FUN = taxon_hit_index)

hit_index_evalue <- list(taxon_hit_all, taxon_hit_best)

pdf(file = "tax_hit_index.pdf")
boxplot(
	x = hit_index_evalue, 
	ylab = "N taxa / N hits", 
	names = c("all e-values", "only best e-value")
)
stripchart(
	x = hit_index_evalue, 
	vertical = TRUE, 
	add = TRUE, 
	method = "jitter",  
	jitter = 0.2, 
	pch = 19, 
	col = hsv(h = 0, s = 1, v = 1, alpha = 0.2), 
	cex = 0.8
)
dev.off()
for(i in 1:length(hit_index_evalue)){
	points(
		x = jitter(rep(i, length(hit_index_evalue[[i]])), factor = 7), 
		y = hit_index_evalue[[i]], 
		pch = 19, 
		col = hsv(h = 0, s = 0, v = 0, alpha = 0.5), 
		cex = 0.8
	)
}

query_seq <- unique(blast_results[ , query_col ])
# or if using factor this might work levels(blast_results[ , query_col ])

# what is the lowest evalue or highest bitscore per query sequence
best_evalue <- lapply(split(blast_results[ , evalue_col ], blast_results[,query_col]), min)
best_bitscore <- lapply(split(blast_results[ , bitscore_col ], blast_results[,query_col]), max)


max(blast_results[,evalue_col])


best_tax <- function(x){
  best_matches <- which(x[ , bitscore_col] == max(x[ , bitscore_col]))
  return(x[best_matches , title_col])
}

best_gi <- function(x, gi_col = 2){
  best_matches <- which(x[ , bitscore_col] == max(x[ , bitscore_col]))
  return(x[best_matches , gi_col])
}

stitles <- lapply(
              split(
                  blast_results, 
                  blast_results[ , query_col] 
                ), 
              best_tax
              )


#

splitter <- function(x){
  Reduce(paste, strsplit(x, split = " ")[[1]][1:2])
}
splitter(stitles[[2]])


Reduce(paste, strsplit(x, split = " ")[[1]][1:2])

# THIS WORKS!!!

do.call(rbind, lapply(stitles[[2]], function(x) Reduce(paste, strsplit(x, split = " ")[[1]][1:2])))
# THIS WORKS!!!

lapply(stitles[[2]], splitter)


unique(as.vector(do.call(rbind, lapply(stitles[[2]], function(x) Reduce(paste, strsplit(x, split = " ")[[1]][1:2])))))
table(as.vector(do.call(rbind, lapply(stitles[[2]], function(x) Reduce(paste, strsplit(x, split = " ")[[1]][1:2])))))

# Holy shit, finally.
taxa_hits <- lapply(
  stitles,
  function(x){
    as.vector(
      do.call(
        rbind, 
        lapply(
          x, function(x){
              Reduce(paste, strsplit(x, split = " ")[[1]][1:2])
            }
          )
        )
    )
  }
)

taxa_unique <- sapply(taxa_hits, unique)
sapply(taxa_hits, table)


# Put this on hold until PICKUP

#############################################
# as one function
tax_list <- function(x){
  best_matches <- which(x[ , bitscore_col] == max(x[ , bitscore_col]))
  titles_sub <- x[best_matches , title_col]
  
}
##############################################





# PICKUP
classification("Quietula y-cauda", db = "itis")



# TAXIZE
########
taxa_unique

# this requires a network connection, and could take a while
tax_hier_ncbi <- list()
for(i in 1:length(taxa_unique)){
  tax_hier_ncbi[[i]] <- classification(taxa_unique[[i]], db = "ncbi")
}

save(tax_hier_ncbi, file = "tax_hier_ncbi_primerbias20151205.RData")

tax_hier_collapser <- function(x)
{
  consol <- do.call(rbind, x)
  return(unique(consol[duplicated(consol), ]))
}
# 2,4,5,7,8
tax_hier_collapser(tax_hier_ncbi[[1]])
tax_hier_ncbi[[9]]

# collapse (intersect taxonomic hierarc)
tax_hier_intersect <- list()
for(i in 1:length(tax_hier_ncbi)){
  tax_hier_intersect[[i]] <- Reduce(intersect, 
         sapply(tax_hier_ncbi[[i]][!is.na(tax_hier_ncbi[[i]])], "[", 1, simplify = TRUE)
  )
}
names(tax_hier_intersect) <- names(taxa_unique)
# STILL a problem with 9! (because of "uncultered organism)
tax_hier_intersect

######## THIS IS IT
best_hit <- sapply(tax_hier_intersect, function(x) x[length(x)])
######## THIS IS IT
names(best_hit)
identical(names(best_evalue), names(best_hit))

FINAL_TABLE <-   cbind(query = names(best_hit), best_hit, best_evalue, best_bitscore)
rownames(FINAL_TABLE) <- NULL
write.csv(
  x = FINAL_TABLE, 
  file = "blast_hit_summary.csv", 
  row.names = FALSE
  )
#
#
#
#
#


Reduce(intersect, sapply(tax_hier_ncbi[[15]], "[", 1, simplify = TRUE))

sapply(tax_hier_ncbi[[15]], length)

which.min(sapply(tax_hier_ncbi[[15]], length))
tax_hier_ncbi[[15]][!is.na(tax_hier_ncbi[[15]])]

lapply(tax_hier_ncbi, is.na)


# probably better to get ids first:
uids <- get_uid(c("Chironomus riparius", "Chaetopteryx"))


# GRAVEYARD
classification("Amphiprion", db = "itis")
specieslist <- c("Abies procera","Pinus contorta")
highertax_itis <- classification(specieslist, db = 'itis')
highertax_ncbi <- classification(specieslist, db = 'ncbi')
highertax_ncbi

clowns <- c("Amphiprion sandaracinos","Premnas biaculeatus", "Cymatogaster aggregata")
clownclass <- classification(clowns, db = "itis")
clown_consol <- do.call(rbind, clownclass)
rownames(clown_consol) <- NULL


unique(clown_consol[duplicated(clown_consol), ])

duplicated(clown_consol)
merge(
  clownclass[[1]],
  clownclass[[2]], 
  by = "name"
)
do.call(merge, clownclass)

Reduce(merge, clownclass)



Reduce(intersect, clownclass)
Reduce(intersect, 
       sapply(clownclass, "[", 1, simplify = TRUE)
       )
clownclass

lapply(highertax_ncbi, "[", c("name", "rank"))

consolidated <- do.call(rbind, highertax_ncbi)

consolidated[duplicated(consolidated), ]



# OTHER PLOTS
# plot bitscore against evalue
plot(
	x = blast_results[, evalue_col], 
	y = blast_results[, bitscore_col], 
	# log = "x", 
	xlab = "evalue", 
	ylab = "bitscore", 
	pch = 21, 
	col = "black", 
	bg = rgb(red = 0, green = 0, blue = 0, alpha = 0.1)
	)

