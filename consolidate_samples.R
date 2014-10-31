# Assess the output from the analysis pipeline

# In which directory are the demultiplexed data folders?
ANALYSIS_DIRECTORY <- "~/Documents/Data/IlluminaData/16S/20141020/Analysis_20141023_1328/Reanalysis_20141028_1351"

# assign the number of tags in that folder to the object N_TAGS
N_TAGS <- length(system(paste("ls ", ANALYSIS_DIRECTORY, sep=""), intern=TRUE))

# Read in the CSV files called "meganout_mod.csv" from each of those folders, and make them into a list of data frames.
# Storing them as a list makes it easy to analyze them in loops or using the apply functions.
TAGS <- list()
for (i in 1:N_TAGS){
	FILE <- paste(ANALYSIS_DIRECTORY, "/tag_", i, "/meganout_mod.csv", sep = "")
	tryCatch({
		TAGS[[i]] <- read.csv(FILE, header = FALSE, col.names=c("ClusterID","N_Reads","Taxon"))
	}, error=function(e){cat("ERROR :",conditionMessage(e), i, "\n")})
}

# create a vector of tag names and apply them to the list of data frames.
TAG_NAMES <- vector()
for (i in 1:N_TAGS){
	TAG_NAMES[i] <- paste("tag_", i, sep = "")
}
names(TAGS) <- TAG_NAMES

# eliminate ones with no reads (this might include those that ERRONEOUSLY have no reads for now!)
TAGS[which(sapply(TAGS, nrow) == 0)] <- NULL

N_READS <- list()
for (i in 1:length(TAGS)){
	N_READS[[i]] <- sort(sapply(split(TAGS[[i]]$N_Reads, TAGS[[i]]$Taxon), sum))
}


reads_per_tag <- vector()
for(i in 1:60){
	reads_per_tag[[i]] <- sum(all_tags[[i]][,2])
}

barplot(reads_per_tag, names=seq(1:length(reads_per_tag)), main="Reads per tagged primer set")

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
