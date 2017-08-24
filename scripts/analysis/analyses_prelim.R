#!/usr/bin/Rscript

# analyze the data from a collection of runs
# Input:
# 1. (PDF) File path to the PDF that will be generated
# 2. (CSV) File path to a table of counts of things (DNA sequences, clusters, duplicates, OTUs, taxa), where columns are samples, and rows are things counted (DNA sequences, OTUs, taxa).
# 3. (CSV) File path to a sequencing pool metadata spreadsheet
# 4. (string) Name of the column of the metadata file containing library names
# 5. (string) Name of the column of the metadata file containing tag sequences
# 6. (string) Name of the column of the metadata file containing sample names
# 7. (string) Name of the column of the metadata file containing sample types

# If you want to run this interactively, just comment out the following line, un-comment the 8 lines after that, and make all of those the appropriate file paths/arguments (see numbers above).

arguments<-commandArgs(TRUE)
arguments <- c(
  "/Users/threeprime/Desktop/debug.pdf",
  "/Users/threeprime/temp_big/20150717_nextseq/Analysis_20151019_1918/all_lib/OTUs_swarm/OTU_table.csv",
  "/Users/threeprime/temp_big/20150717_nextseq/SEQUENCING_POOL_20150618.csv",
  "library",
  "tag_sequence",
  "sample_name",
  "sample_type")


pdf_file			<- arguments[1]
otu_file			<- arguments[2]
metadata_file	<- arguments[3]
col_library		<- arguments[4]
col_tagseq		<- arguments[5]
col_samplename	<- arguments[6]
col_sampletype	<- arguments[7]



# load required packages
library(vegan)

# initialize PDF
pdf(file = pdf_file)

par_orig <- par()

################################################################################
# READ IN THE DATA (and do basic cleaning)

# OTU TABLE
DATA <- read.csv(otu_file, row.names = 1)

# for clustered data, replace "DUP" with "OTU"
# rownames(DATA) <- gsub("DUP", "OTU", rownames(DATA))

# transpose OTU data to the appropriate orientation (samples are rows, OTUs are columns)
DATA <- t(as.matrix(DATA))

# order the duplicates/OTUs by decreasing order of abundance
DATA <- DATA[,order(colSums(DATA), decreasing = TRUE)]

# check matrix is the expected dimensions
# dim(DATA)


# METADATA
# Read in metadata spreadsheet from labwork, which contains columns of sample names and corresponding tag sequences
# for original formatting see "/Users/threeprime/Documents/GoogleDrive/Data_Illumina/16S/run_20141113_time_series/sample_data.csv"
metadata <- read.csv(metadata_file)

# eliminate any columns that aren't going to be used
metadata <- metadata[,c(col_library, col_tagseq, col_samplename, col_sampletype)]






################################################################################
################################################################################
################################################################################
# check that the two data files can be linked up

# The rownames of the OTU table should ID/differentiate distinct samples from the sequencer
# This information was glommed together from the metadata by banzai, so we should be able to go back and piece things back together

# make a new column in the metadata that will link the OTU file and the sequencing pool
# NOTE: The following line must directly correspond to the naming scheme used by banzai!
# originally was a separate vector: sample_id <- paste("lib", metadata[,col_library], "tag", metadata[,col_tagseq], sep = "_")
metadata <- cbind(
				metadata, 
				sample_id = paste(
					"lib", 
					metadata[,col_library], 
					"tag", 
					metadata[,col_tagseq], 
					sep = "_"
					)
				)


# Check for differences between the sample IDs in metadata and OTU table.
if(length(setdiff(rownames(DATA), metadata$sample_id)) == 0 ){
	print("The sample IDs in the metadata and OTU table match up -- great jorb!")
} else {
	print("The sample IDs in the metadata and OTU table do not match up.")
	# Which has the 
	maxboth <- max(length(rownames(DATA)), length(metadata$sample_id))
	
	# how many of the OTU table rownames are contained in the "sample_id" column
	inboth <- sum(rownames(DATA) %in% metadata$sample_id)	
}


! 
rownames(DATA) %in% metadata$sample_id
metadata$sample_id %in% rownames(DATA)



if( inboth == 0){
	stop("Sample identifiers in metadata and OTU table do not correspond.")
} else if(inboth == maxboth){
	
} else if(){
	
} else {
	print("something else")
}

paste(nrow(DATA), "rows found in OTU table,", nrow(metadata), "rows found in metadata")
nrow(DATA) 
if(TRUE){
	stop(nrow(DATA), " rows found in OTU table, ", nrow(metadata), " rows found in metadata, ")
}
print()

length(metadata$sample_id)

# link the metadata file to the OTU file
tag_to_sequencing_data <- match(rownames(DATA), metadata$sample_id)
tag_to_samplename <- metadata[tag_to_sequencing_data, arguments[6]]

################################################################################
################################################################################
################################################################################










# plot the number of reads per sample binned by sample origin (environmental samples and controls)
data_by_sample <- split(rowSums(DATA), tag_to_samplename)
stripchart(
	data_by_sample, 
	las = 2, 
	# cex.axis = 0.4, 
	pch = 1, 
	vertical = TRUE, 
	xaxt = "n", 
	main = "Reads per sample"
	)
	
axis(
	side = 1, 
	at = 1:length(data_by_sample), 
	labels = names(data_by_sample), 
	cex.axis = 0.4, 
	las = 2
	)

SD_reads <- c(
	median(rowSums(DATA)) - sd(rowSums(DATA)), 
	median(rowSums(DATA)) + sd(rowSums(DATA))
	)
	
polygon(
	x = c(0, length(data_by_sample)+1, length(data_by_sample)+1, 0), 
	y = rep(SD_reads, each = 2), 
	col = "#BEBEBE50", 
	border = NA
	)
	
abline(h = median(rowSums(DATA)), lty = 2)

# plot total sum of reads per dup/OTU across samples
par() <- par_orig
plot(
	colSums(DATA), 
	# log = "y", 
	pch = 20, 
	cex = 0.5, 
	main = "Total reads per OTU across samples", 
	xlab = "OTU rank", 
	ylab = "Number of reads"
	)

# create matrix of tag sequence and library number
TAG_LIB <- strsplit(rownames(DATA), "_")
TAG_LIB <- do.call(rbind, TAG_LIB)[,c(4,2)]
colnames(TAG_LIB) <- c("Tag", "Lib")

# bind tag sequence and library number to OTU table; must be a dataframe to store numbers/text
DATA.df <- cbind(TAG_LIB, as.data.frame(DATA))


# Order the OTU data the same as the sequencing pool sample data
DATA.df <- DATA.df[
					match(
					interaction(metadata[c("tag_sequence", "library")]), 
					interaction(DATA.df[c("Tag", "Lib")])),
				]


# make a vector of sample names in the right order (corresponding to the tag sequence from OTU table)
# sample_name <- metadata$sample_name[match(interaction(DATA.df[c("Tag", "Lib")]), interaction(metadata[c("tag_sequence", "library")]))]
# note I think it is safer to reorder the OTU data to match up with the sequencing pool first,
# thus, just use
# metadata$sample_name

# for some reason there was an empty level, which caused problems down the road
# sample_name <- droplevels(sample_name)

# SAMPLE TYPE
# Extract the data from the sequencing sample pool spreadsheet
# metadata$sample_type
# make a vector of sample types (environmental, tissue, filter blank)
# first just make them all environmental
# sample_type <- rep("environ", nrow(DATA.df))
# switch the tissue and filter blanks on the basis of the sample name
# sample_type[which(sample_name == "TILAPIA")] <- "tissue"
# sample_type[which(sample_name == "DIH20-20140709")] <- "filter_blank"

# add those to the data frame and check it out
DATA.df <- cbind(
	sample_name = metadata[, arguments[6]], 
	sample_type = metadata[,arguments[7]], 
	DATA.df
	)
	
# DATA.df[,1:5]

# Incorporate sample name into rownames of DATA (otu table stored as matrix)
# rownames(DATA) <- paste(sample_name, rownames(DATA), sep = "_")

# order by sample name
DATA.df <- DATA.df[order(DATA.df$sample_name),]

# plot number of reads per sample:
# reads_per_sample <- rowSums(DATA.df[,5:ncol(DATA.df)])
# plot(sort(reads_per_sample), main = "reads per sample", xaxt = "n", ann=FALSE)
# axis(side = 1, at = 1:length(reads_per_sample), labels = names(sort(reads_per_sample)), las = 2, cex.axis = 0.5)
# legend("topleft", legend = , col = )


# community analyses
sample_shan <- as.numeric(diversity(DATA.df[,5:ncol(DATA.df)])); plot(sample_shan)
sample_rich <- as.numeric(rowSums(DATA.df[,5:ncol(DATA.df)] > 0)); plot(sample_rich)
sample_reads <- as.numeric(rowSums(DATA.df[,5:ncol(DATA.df)]))
# plot(sample_reads) # no need to plot this; plotted above, but this is in order as above

# PLOT RICHNESS AGAINST READS
plot(sample_reads, sample_rich, xlab = "Reads per Sample", ylab = "Total Number of OTUs")
lm_rich_reads <- lm(sample_rich~sample_reads)
# summary(lm_rich_reads) # MAKE THIS PRINT

split_shan <- split(sample_shan, paste(DATA.df[, "sample_name"], DATA.df[, "Tag"], DATA.df[, "Lib"]))
split_rich <- split(sample_rich, paste(DATA.df[, "sample_name"], DATA.df[, "Tag"], DATA.df[, "Lib"]))
split_reads <- split(sample_reads, paste(DATA.df[, "sample_name"], DATA.df[, "Tag"], DATA.df[, "Lib"]))

# BOXPLOT SHANNON INDEX
par(mar=c(4,6,1,1), cex.axis=0.5)
boxplot(
	split_shan,
	horizontal=TRUE,
	las=1,
	cex.names=0.2,
	border = as.numeric(DATA.df[,2]),
	xlab = "Shannon Index"
	# border=c(1, rep(2:4, each = 2 , len = length(split_shan)-1))
)
abline(h = which(!duplicated(rev(DATA.df[,"sample_name"]))) -0.5, lty = 2)

# PLOT CLUSTER RICHNESS
# pdf(file = "cluster_richness.pdf", width = 10, height = 7)
par(mar=c(4,6,1,1), cex.axis=0.5)
stripchart(
	x = rev(split_rich),
	vertical = FALSE,
	method = "jitter",
	pch = 21, col = "aquamarine4", bg = "aquamarine",
	las = 1,
	xlab = "Number of Clusters"
	# xlim = c(0,3000)
	# group.names=as.character(reps)
)
abline(h = which(!duplicated(rev(DATA.df[,"sample_name"]))) -0.5, lty = 2, col = "lightgrey")

# dev.off()



# PLOT TOTAL READS
# pdf(file = "total_reads.pdf", width = 10, height = 7)
par(mar=c(4,6,1,1), cex.axis=0.5)
stripchart(
	x = split_reads,
	vertical = FALSE,
	# method = "jitter",
	pch = 21, col = "aquamarine4", bg = "aquamarine",
	xlab = "Reads per Sample",
	# xlim = c(0,3000)
	# group.names=as.character(reps)
	las = 1
)
abline(h = which(!duplicated(rev(DATA.df[,"sample_name"]))) -0.5, lty = 2, col = "lightgrey")
# dev.off()


# PLOT SHANNON INDEX
# pdf(file = "shannon_index.pdf", width = 10, height = 7)
par(mar=c(4,6,1,1), cex.axis=0.5)
stripchart(
	x = split_shan,
	vertical = FALSE,
	method = "jitter",
	pch = 21, col = "aquamarine4", bg = "aquamarine",
	las = 1,
	xlab = "Shannon Index"
	# xlim = c(0,3000)
	# group.names=as.character(reps)
)
abline(h = which(!duplicated(rev(DATA.df[,"sample_name"]))) -0.5, lty = 2, col = "lightgrey")
# dev.off()


##########################################################################################################
# Correct richness for number of reads
# library(vegan)
# S <- specnumber(BCI) ## rowSums(BCI > 0) does the same...
# plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
# abline(0, 1)
# rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)

# Note the removal of the controls!
DATA.env <- droplevels(DATA.df[DATA.df$sample_type=="environmental",])
rownames(DATA.env) <- paste(DATA.env$sample_name, DATA.env$Tag, DATA.env$Lib)
DATA.mat <- as.matrix(DATA.env[,5:ncol(DATA.env)])
# DATA.env[,1:6]

# Calculate richness for the environmental samples
rich_env <- as.numeric(rowSums(DATA.env[,5:ncol(DATA.df)] > 0))


# Calculate the total reads per sample; use the smallest as the value to rarefy samples against
# (If using controls, remember the filtration blanks are the smallest three, and have an order of magnitude fewer reads)
min_reads <- min(rowSums(DATA.mat))
RARE <- rarefy(DATA.mat, min_reads)
# str(RARE)
rare_rich <- split(as.numeric(RARE), paste(DATA.env[, "sample_name"], DATA.env[, "Tag"], DATA.env[, "Lib"]))

# PLOT CLUSTER RICHNESS
# pdf(file = "rarefied_richness.pdf", width = 10, height = 7)
par(mar=c(4,6,1,1), cex.axis=0.5)
stripchart(
	x = rev(rare_rich),
	vertical = FALSE,
	method = "jitter",
	pch = 21, col = "aquamarine4", bg = "aquamarine",
	las = 1,
	xlab = "Rarefied cluster richness"
	# xlim = c(0,3000)
	# group.names=as.character(reps)
)
abline(h = which(!duplicated(rev(DATA.env[,"sample_name"]))) -0.5, lty = 2, col = "lightgrey")
# dev.off()


# Plot rarefaction curves
rarecurve(DATA.mat, step = 10000, sample = min_reads, col = "blue", cex = 0.6, label=FALSE, main = "Rarefaction Curves")
points(rowSums(DATA.mat), rich_env, pch=20, cex = 0.5, col = DATA.env$sample_name)

# make a vector of tag levels
tag_levels <- unname(do.call(c, lapply(split(DATA.env[,"Tag"], DATA.env[,"sample_name"]), function(x) as.numeric(droplevels(x)))))

rich_df <- data.frame(DATA.env[,1:4], tag_levels, rarefied_richness=as.numeric(RARE))

# aov_rich <- aov(rarefied_richness ~ sample_name + tag_levels*Lib, data = rich_df)
# aov_rich <- aov(rarefied_richness ~ sample_name + Lib*tag_levels, data = rich_df)
# summary(aov_rich)
# as.matrix(aov_rich)
# capture.output(summary(aov_rich),file="/Users/threeprime/test.txt")


# mean(rich_df$rarefied_richness)
# sd(rich_df$rarefied_richness)
# plot(rich_df$rarefied_richness, col=rich_df$sample_name)



# plot proportion of reads for a sampling of OTUs
reads_per_tag <- rowSums(DATA.df[,5:ncol(DATA.df)])
par(mfrow=c(2,2))
for(i in 1:4){
	plot(DATA.df[,4+i]/reads_per_tag, ylab = "Proportion of reads", xaxt = "n", xlab="", pch=19, cex = 0.9, main=paste("Cluster ", i))
	abline(v = which(!duplicated(DATA.df[,"sample_name"])) -0.5, lty = 2, col = "lightgrey")
	axis(1, at=1:nrow(DATA.df), labels = DATA.df$sample_name, las = 2, cex.axis = 0.8)
}


########################################################################################################
# COMMUNITY LEVEL INDICES
########################################################################################################
# Compute distance metric between each replicate (library) for each tag (25 in the case of the time series data)
## takes a while...
## note 20150326; changed "DATA.df" to "DATA.env" to run this on ONLY the environmental samples.
BC_within_tag <- sapply(split(DATA.env[,5:ncol(DATA.env)], DATA.env$Tag), vegdist)

# only two replicates per tag, so mean is unnecessary
# BC_within_tag_mean <- apply(BC_within_tag, 2, mean)

## Compute distance metric between each replicate (library) for each tag (25 in the case of the time series data -- 22 excluding controls)
BC_between_tag <- sapply(split(DATA.env[,5:ncol(DATA.env)], DATA.env$sample_name), vegdist, simplify = FALSE)
## drop samples with fewer than 15 points ()
# BC_between_tag <- BC_between_tag[-which(sapply(BC_between_tag, length) < 15)]
## convert to standard (i.e. not distance) matrix
BC_between_tag <- lapply(BC_between_tag, as.matrix)

lower_quadrant_mean <- function(x){
	mean_x <- mean(x[2:4, 1:2])
	mean_x
}


BC_between_tag_mean <- sapply(BC_between_tag, lower_quadrant_mean)

########################################################################################################
# GENERATE A BOXPLOT
########################################################################################################
# pdf(file = "BC_dis_boxplot.pdf")
# boxplot(
	# x = list(BC_within_tag, BC_between_tag_mean),
	# ylim = c(0,1),
	# ylab = "Mean Bray-Curtis Dissimilarity",
	# xaxt = "n"
# )
# axis(
	# side = 1,
	# at = c(1,2),
	# tick = FALSE,
	# line = 1,
	# labels = c(
		# paste("within tags,\namong libraries\nN = ", length(BC_within_tag), sep = ""),
		# paste("among tags,\namong libraries\nN = ", length(BC_between_tag_mean), sep = "")
	# )
# )
# dev.off()
########################################################################################################



########################################################################################################
# GENERATE A STRIPCHART
########################################################################################################
# pdf(file = "BC_dis_stripchart.pdf")
par(mfrow=c(1,1))
stripchart(
	x = list(BC_within_tag, BC_between_tag_mean),
	vertical = TRUE,
	method = "jitter",
	jitter = 0.2,
	pch = 1, col = "black", #bg = "grey",
	cex = 0.8,
	# las = 2,
	ylab = "Mean Bray-Curtis Dissimilarity",
	ylim = c(0,1),
	xaxt = "n",
	group.names = c("within tags,\namong libraries", "between tags,\nbetween libraries")
)
axis(
	side = 1,
	at = c(1,2),
	tick = FALSE,
	line = 1,
	labels = c(
		paste("within tags,\namong libraries\nN = ", length(BC_within_tag), sep = ""),
		paste("among tags,\namong libraries\nN = ", length(BC_between_tag_mean), sep = "")
	)
)
# dev.off()
########################################################################################################






########################################################################################################
# OTU specific effects
########################################################################################################

# First, consider which OTUs comprise the bulk of the sequence data across samples, to do this,
## How many reads were just environmental samples?
reads_total_env <- sum(rowSums(DATA.env[,5:ncol(DATA.env)]))

## What proportion of those reads were from each OTU?
propOTU_total_env <- colSums(DATA.env[,5:ncol(DATA.env)])/reads_total_env

## What proportion of the reads are comprised of the ten most abundant OTUs?
top_ten_cumsum <- sum(sort(propOTU_total_env, decreasing = TRUE)[1:10])
top10OTU_env <- names(sort(propOTU_total_env, decreasing = TRUE)[1:10])


# Make a new data frame identical to DATA.df except that instead of raw read numbers, numbers are the proportional contribution of each OTU
DATA.prop <- data.frame(DATA.df[1:4], DATA.df[,5:ncol(DATA.df)] / rowSums(DATA.df[,5:ncol(DATA.df)]))

# to confirm this worked, make sure columns sum to 1!
# rowSums(DATA.prop[,5:ncol(DATA.prop)])

# consider only the environmental samples, and only the top ten most abundant OTUs:
DATA.prop_env <- DATA.prop[DATA.prop$sample_type == "environmental", c(colnames(DATA.prop)[1:4],  top10OTU_env)]

# make a vector of sample tags (now we have tag_levels vector)
# sample_tag <- rep(as.factor(c(1,1,1,2,2,2)), 11)

## insert into data frame
DATA.prop_env <- data.frame(DATA.prop_env[,1:4] , tag_levels, DATA.prop_env[,5:ncol(DATA.prop_env)])

########################################################################################################
# PLOT ABUNDANCE OF EACH OTU IN EACH ENVIRONMENTAL SAMPLE BY TAG
# pdf(file = "OTU_abundance_by_tag.pdf", width = 5, height = 5)
LINES_POS <- which(!duplicated(DATA.prop_env$sample_name))/2
for(i in top10OTU_env){ # for just OTU_1 use "[1]"
	plot(DATA.prop_env[ DATA.prop_env["tag_levels"] == 1 , i], cex = 0.8, xaxt='n', xlab = "Environmental Sample Number", ylim = c(0, max(DATA.prop_env[, i ])), ylab = "proportional abundance", main = i)
	points(DATA.prop_env[ DATA.prop_env["tag_levels"] == 2 , i], pch = 4, cex = 0.8)
	abline(v = LINES_POS, lty = 2, col = "grey")
	axis(side = 1, at = LINES_POS + 1.5, labels = seq(1:length(which(!duplicated(DATA.prop_env$sample_name)))))
	legend("topright", legend = c("tag 1", "tag 2"), pch = c(1,4), bg = "white")
}
dev.off()
