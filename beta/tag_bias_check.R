# analyze time series data

# initialize PDF
# pdf(file = "analysis_output.pdf")
# pdf(file = "/Users/threeprime/Documents/Dropbox/Kelly_Lab/")

# READ IN THE DATA
# unclustered:
# DATA <- read.csv("~/Documents/Data/IlluminaData/16S/run_20141113_time_series/all_libraries/dups.csv")
# clustered:
# DATA <- read.csv("/Users/threeprime/Documents/GoogleDrive/Data_Illumina/16S/run_20141113_time_series/all_libraries/all_clusters.csv", row.names = 1)
DATA <- read.csv("/Users/threeprime/Documents/GoogleDrive/Data_Illumina/16S/test/Analysis_20150403_0109/all_lib/dups.csv", row.names = 1)

# for clustered data, replace "DUP" with "OTU"
# rownames(DATA) <- gsub("DUP", "OTU", rownames(DATA))

# transpose to the appropriate orientation (samples are rows, OTUs are columns)
DATA <- t(as.matrix(DATA))

# order the duplicates/OTUs by decreasing order of abundance
DATA <- DATA[,order(colSums(DATA), decreasing = TRUE)]

# check matrix is the expected dimensions
dim(DATA)

# view number of reads per sample:
reads_per_sample <- rowSums(DATA)
plot(reads_per_sample, main = "reads per sample")

# view total sum of reads per dup/OTU across samples
plot(colSums(DATA), pch = 20, cex = 0.5, main = "total reads per OTU")

# create matrix of tag sequence and library number
TAG_LIB <- strsplit(rownames(DATA), "_")
TAG_LIB <- do.call(rbind, TAG_LIB)[,c(3,1)]
colnames(TAG_LIB) <- c("Tag", "Lib")

# bind tag sequence and library number to OTU table; must be a dataframe to store numbers/text
DATA.df <- cbind(TAG_LIB, as.data.frame(DATA))

# Read in spreadsheet from labwork, which contains columns of sample names and corresponding tag sequences
# for original formatting see "/Users/threeprime/Documents/GoogleDrive/Data_Illumina/16S/run_20141113_time_series/sample_data.csv"
SAMPLES <- read.csv("/Users/threeprime/Documents/GoogleDrive/Data_Illumina/16S/run_20150401/20150317_sequencing_pool.csv")

# Order the OTU data the same as the sequencing pool sample data

DATA.df <- DATA.df[match(interaction(SAMPLES[c("tag_sequence", "library")]), interaction(DATA.df[c("Tag", "Lib")])),]


# rows after row 25 are garbage, delete them.
# SAMPLES <- SAMPLES[1:25,]

# make a vector of sample names in the right order (corresponding to the tag sequence from OTU table)
# sample_name <- SAMPLES$sample_name[match(interaction(DATA.df[c("Tag", "Lib")]), interaction(SAMPLES[c("tag_sequence", "library")]))]
# note I think it is safer to reorder the OTU data to match up with the sequencing pool first,
# thus, just use
SAMPLES$sample_name

# for some reason there was an empty level, which caused problems down the road
# sample_name <- droplevels(sample_name)

# SAMPLE TYPE
# Extract the data from the sequencing sample pool spreadsheet
SAMPLES$sample_type
# make a vector of sample types (environmental, tissue, filter blank)
# first just make them all environmental
# sample_type <- rep("environ", nrow(DATA.df))
# switch the tissue and filter blanks on the basis of the sample name
# sample_type[which(sample_name == "TILAPIA")] <- "tissue"
# sample_type[which(sample_name == "DIH20-20140709")] <- "filter_blank"

# add those to the data frame and check it out
DATA.df <- cbind(sample_name = SAMPLES$sample_name, sample_type = SAMPLES$sample_type, DATA.df)
DATA.df[,1:5]

# Order by sample name
DATA.df <- DATA.df[order(DATA.df$sample_name),]

# Incorporate sample name into rownames of DATA (otu table stored as matrix)
# rownames(DATA) <- paste(sample_name, rownames(DATA), sep = "_")

# order by sample name
DATA.df <- DATA.df[order(DATA.df$sample_name),]

# load vegan to do some community analyses
library(vegan)
sample_shan <- as.numeric(diversity(DATA.df[,5:ncol(DATA.df)])); plot(sample_shan)
sample_rich <- as.numeric(rowSums(DATA.df[,5:ncol(DATA.df)] > 0)); plot(sample_rich)
sample_reads <- as.numeric(rowSums(DATA.df[,5:ncol(DATA.df)])); plot(sample_reads)

# PLOT RICHNESS AGAINST READS
plot(sample_reads, sample_rich, xlab="Reads per Sample", ylab = "Total Number of Clusters")
lm_rich_reads <- lm(sample_rich~sample_reads)
summary(lm_rich_reads)

split_shan <- split(sample_shan, paste(DATA.df[, "sample_name"], DATA.df[, "Tag"], DATA.df[, "Lib"]))
split_rich <- split(sample_rich, paste(DATA.df[, "sample_name"], DATA.df[, "Tag"], DATA.df[, "Lib"]))
split_reads <- split(sample_reads, paste(DATA.df[, "sample_name"], DATA.df[, "Tag"], DATA.df[, "Lib"]))

# BOXPLOT SHANNON INDEX
par(mar=c(4,6,1,1), cex.axis=0.5)
boxplot(
	rev(split_shan),
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
dev.off()


##########################################################################################################
# Correct richness for number of reads
library(vegan)
# S <- specnumber(BCI) ## rowSums(BCI > 0) does the same...
# plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
# abline(0, 1)
# rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)

# Note the removal of the controls!
DATA.env <- DATA.df[DATA.df$sample_type=="environmental",]
DATA.env <- droplevels(DATA.env)
rownames(DATA.env) <- paste(DATA.env$sample_name, DATA.env$Tag, DATA.env$Lib)
DATA.mat <- as.matrix(DATA.env[,5:ncol(DATA.env)])
DATA.env[,1:6]

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

aov_rich <- aov(rarefied_richness ~ sample_name + tag_levels*Lib, data = rich_df)
aov_rich <- aov(rarefied_richness ~ sample_name + Lib*tag_levels, data = rich_df)
summary(aov_rich)
as.matrix(aov_rich)
capture.output(summary(aov_rich),file="/Users/threeprime/test.txt")


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
rowSums(DATA.prop[,5:ncol(DATA.prop)])

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
# dev.off()
########################################################################################################
# greatest difference for OTU 1 is in Sample 9


########################################################################################################
# difference in mean proportional abundance among tags
#### For each replicate, calculate the difference between mean proportion using each tag for each of the top 10 OTUs
DIFF_OTU <- vector(mode = "list", length = length(top10OTU_env))
names(DIFF_OTU) <- top10OTU_env
MEAN_OTU <- vector(mode = "list", length = length(top10OTU_env))
names(MEAN_OTU) <- top10OTU_env
for(i in top10OTU_env){
	TEMP	<-	split(DATA.prop_env[c("tag_levels", i)], DATA.prop_env$sample_name, drop = TRUE)
	MEAN_OTU[i][[1]]	<-	lapply(TEMP, function(x) sapply(split(x[,2], x[,1]), mean))
	DIFF_OTU[i][[1]]	<-	sapply(TEMP2, diff)
}

#which sample had the greatest between-tag difference in mean abundance for OTU 1?
ENV_SAMPLE_MAX_DIFF_OTU1 <- which.max(sapply(MEAN_OTU$OTU_1, diff))[[1]]
# what were those values?
MEAN_OTU["OTU_1"][[1]][ENV_SAMPLE_MAX_DIFF_OTU1]

# Plot the difference in mean proportaional abundance among tags
# pdf(file = "tag_diff_by_OTU.pdf", width = 14, height = 6)
par(mar = c(4,4,1,1), mfrow = c(2,5))
for(i in 1:length(DIFF_OTU)){
	plot(DIFF_OTU[i][[1]], xlab = "environmental sample", ylab = "difference in mean proportional abundance among tags")
}
# dev.off()
########################################################################################################

# Could do some T tests...
# T_TESTER <- function(x){
	# for(i in unique(x$sample_name)) {
		# DAT <- x[which(x$sample_name == i),]
		# t.test(DAT$OTU_1, DAT$Tag)
		# print(DAT)
	# }
# }
# T_TESTER(DATA.prop_env)




DATA.prop_10 <- DATA.prop[,1:14]
diff(range(TEMP[12][[1]][4:6]))
var(TEMP[12][[1]][1:3])


colnames(DATA.prop)[5:14]

OTU01 <- lapply(
	split(
		data.frame(tag = DATA.df$Tag, prop = DATA.df[,5] / rowSums(DATA.df[,5:ncol(DATA.df)])),
		paste(DATA.df$sample_name)	
	),
	function(x) split(x[,2], droplevels(x[,1]))
)

mean(OTU01[2][[1]][1][[1]]) - mean(OTU01[2][[1]][2][[1]])



split(
	DATA.df[,5] / rowSums(DATA.df[,5:ncol(DATA.df)]), 
	paste(DATA.df$sample_name, DATA.df$Tag)
	)


# Made it to here 20150323



########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################



# split(DATA.df[,:10], DATA.df$sample_name)
# # Compute distance metric (Bray-Curtis) for each sample, within and among tags
# ## Data must first be in matrix format
# DATA.mat <- as.matrix(DATA.df[,5:ncol(DATA.df)])
# rownames(DATA.mat) <- DATA.df$sample_name


########################################################################################################
### Plot a dendrogram
## make distance matrix. Don't actually need to compute ALL values, but this only takes ~2s for 27973 OTUs and 75 samples
DATA.dist <- vegdist(DATA.df[,5:ncol(DATA.df)], method = "bray", diag=TRUE, upper = TRUE)
# or compute on DATA.mat

par(mar=c(1,1,1,1))
plot(
	hclust(DATA.dist),
	axes=FALSE,
	# main="Cluster Dendrogram of Clustered Time Series Samples (Bray-Curtis)",
	main = "", xlab="",
	ylab="",
	sub="")
mtext("Cluster Dendrogram of Time Series Samples (Bray-Curtis on Clustered sequences)", side = 2)

## to color the dendrogram:
### first, make the cluster output a dendrogram object
DENDRO <- as.dendrogram(hclust(DATA.dist))

### get leaf labels
LEAVES <- levels(as.factor(labels(DENDRO)))

### function to make colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

labelColors <- gg_color_hue(length(LEAVES))
names(labelColors) <- levels(as.factor(labels(DATA.dist)))

### function to assign label colors
colLab <- function(n) {
    if (is.leaf(n)) {
        a <- attributes(n)
        labCol <- labelColors[a$label][[1]]
        attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    n
}

clusDendro <- dendrapply(DENDRO, colLab)
plot(clusDendro)







# load field data, calculate time since tides
TIDES <- read.table("~/Documents/Projects/eDNA_Puget_Sound/Data/tides_shilshole_201406.txt", header=TRUE)
FIELD_DATA <- read.csv("~/Desktop/field_data.csv")
FIELD_DATA <- cbind(FIELD_DATA, "DateTime" = as.POSIXct(paste(FIELD_DATA$date_collected, FIELD_DATA$time_collected)))

# If elapsed time since first sample isn't already included in data frame:
difftime(FIELD_DATA$DateTime, FIELD_DATA$DateTime[1], units = "hours")
TIDES$Date <- as.Date(TIDES$Date)
TIDES <- TIDES[TIDES$Date > "2014/06/24",]
TIDES <- cbind(TIDES, "DateTime"=as.POSIXct(paste(TIDES$Date, TIDES$Time)))

# plot sampling points in time wrt Tide
plot(TIDES$DateTime[3:19], TIDES$Pred[3:19], type="b", ylab="Tidal Height (ft)", xlab="Time", main = "Sample Times WRT Tide")
abline(v = FIELD_DATA$DateTime, col="red")

nearest_tides <- data.frame(nearestL=numeric(), nearest_H=numeric())
for(i in 1:length(FIELD_DATA$DateTime)){
	time_to_H <- difftime(FIELD_DATA$DateTime[i], TIDES[TIDES$High.Low == "H",]$DateTime)
	time_to_L <- difftime(FIELD_DATA$DateTime[i], TIDES[TIDES$High.Low == "L",]$DateTime)
	nearest_L <- time_to_L[which(abs(time_to_L) == min(abs(time_to_L)))]
	nearest_H <- time_to_H[which(abs(time_to_H) == min(abs(time_to_H)))]
	nearest_tides[i,1] <- nearest_L
	nearest_tides[i,2] <- nearest_H
}










# individual duplicate plots
# for(i in 1:40){
	# pdf(file= paste("DUP_",i,".pdf",sep=""), width=9, height=7)
	# par(mar=c(4,10,5,2))
	# boxplot(split(DATA.df[,i+3], DATA.df$sample_name), horizontal=TRUE, las=1)
	# dev.off()
# }

# This is probably bad practice: it creates duplicate rownames. Not sure why R allows it:
# rownames(DATA) <- DATA.df$sample_name

barplot(DATA/rowSums(DATA), horiz=TRUE, las=1)
barplot(t(DATA),horiz=TRUE, las=1)


N_READS <- colSums(DATA)
boxplot(colSums(DATA), ylab="number of reads per sample")
mtext("N = 25")


head(DATA)
DATA.df10<-DATA.df[,1:13]
bu <- DATA.df10
rownames(DATA.df10)<-DATA.df10$sample_name

DATA.df10 <- DATA.df10[with(DATA.df10, order(sample_name)),]

rev(rowSums(DATA.df10[4:13]))
DATA.mat <- as.matrix(DATA.df10[,4:10])
rowSums(DATA.mat)

pdf(file="time_series_community6.pdf", width=10, height=7)
barplot(t(DATA.mat/rowSums(DATA.mat)), horiz=TRUE, names=DATA.df10$sample_name,
cex.names = 0.45, las=1,
col=c('#3168FF', '#00CCFF', '#98CBF8', '#CCFCCC', '#00FD03', '#97CB00', '#FFFD00', '#FE9900', '#FB02FE', '#FF0200'))
dev.off()



library(reshape2)
DATA.melt <- melt(DATA.mat)

colnames(DATA.melt) <- c("Sample", "Taxon", "Reads")

mycols <- c('#3168FF', '#00CCFF', '#98CBF8', '#CCFCCC', '#00FD03', '#97CB00', '#FFFD00', '#FE9900', '#FB02FE', '#FF0200') # Set milieu colors

ggplot(DATA.melt) +
	geom_bar(aes(factor(Sample), fill=fill(Taxon), position='fill')


	+
	scale_fill_manual(values = mycols)
