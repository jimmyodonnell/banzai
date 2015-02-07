# Assess the output from the analysis pipeline

# In which directory are the demultiplexed data folders?
ANALYSIS_DIRECTORY <- "/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141030_2020/demultiplexed"

# Identify a file containing sample metadata.
# This must be a csv file with a header row that contains at least two columns, named "Tag_Sequence" and "Sample_Type"
# column names must be EXACT, including capitalization and whitespace.
SAMPLE_METADATA <- read.csv("~/Documents/Projects/eDNA_Hopkins/Data/12S_samples/12S_tagged_run_metadata.csv")

# assign the number of tags in that folder to the object N_TAGS
TAG_FOLDERS <- system(paste("ls ", ANALYSIS_DIRECTORY, sep=""), intern=TRUE)
N_TAGS <- length(TAG_FOLDERS)

# Read in the CSV files called "meganout_mod.csv" from each of those folders, and make them into a list of data frames.
# Storing them as a list makes it easy to analyze them in loops or using the apply functions.
TAGS <- list()
for (i in TAG_FOLDERS){
	FILE <- paste(ANALYSIS_DIRECTORY, "/", i, "/meganout_mod.csv", sep = "")
	tryCatch({
		TAGS[[i]] <- read.csv(FILE, header = FALSE, col.names=c("ClusterID","N_Reads","Taxon"))
	}, error=function(e){cat("ERROR :",conditionMessage(e), i, "\n")})
}

# create a vector of tag names and apply them to the list of data frames.
# TAG_NAMES <- vector()
# for (i in 1:N_TAGS){
# 	TAG_NAMES[i] <- paste("tag_", i, sep = "")
# }
# names(TAGS) <- TAG_NAMES

# eliminate ones with no reads (this might include those that ERRONEOUSLY have no reads for now!)
EMPTY_TAGS <- names(TAGS[which(sapply(TAGS, nrow) == 0)])
TAGS[which(sapply(TAGS, nrow) == 0)] <- NULL

N_READS <- list()
for (i in 1:length(TAGS)){
	N_READS[[i]] <- sort(sapply(split(TAGS[[i]]$N_Reads, TAGS[[i]]$Taxon), sum))
}

names(N_READS) <- gsub("tag_", "", names(TAGS))

reads_per_tag <- sapply(N_READS, sum)

READS <- reads_per_tag[match(SAMPLE_METADATA$Tag_Sequence, names(N_READS))]
SAMPLE_METADATA <- cbind(SAMPLE_METADATA, READS)

sample_size <- paste("N=",table(SAMPLE_METADATA$Sample_Type), sep="")

boxplot(SAMPLE_METADATA$READS~SAMPLE_METADATA$Sample_Type, names=paste(levels(SAMPLE_METADATA$Sample_Type), "\n", sample_size), ylab="Number of Reads", xlab="Sample Type", main="Reads per tagged 12S primer set")

reads_by_type <- reads_per_tag


barplot(reads_per_tag, names=c(SampleTypes), main="Reads per tagged primer set", las=2)

# these are all of the taxon groups
TAXA <- sort(unique(do.call(c, sapply(N_READS, names))))

# Build master dataframe
DATA <- matrix(data=NA, nrow = length(N_READS), ncol=length(TAXA), dimnames=list(names(TAGS), TAXA))

for (i in 1:length(N_READS)){
	DATA[i,] <- N_READS[[i]][match(TAXA, names(N_READS[[i]]))]
}

# Replace NAs with zeros
DATA[is.na(DATA)] <- 0

# write the data frame
write.csv(DATA, file = "Reads_per_tag_by_OTU.csv")

# Examples:
# Boxplot of number of reads by taxon
par(mar=c(5,10,1,1))
boxplot(DATA, horizontal=TRUE, las=2)

# Boxplot of number of Scombridae reads by tag (or sample)
barplot(DATA[,"Scombridae"], main = "Number of Scombridae reads per tag")
barplot(DATA[,"Euteleostomi"], main = "Number of Eutelostomi reads per tag")
