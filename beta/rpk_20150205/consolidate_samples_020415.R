#!usr/bin/Rscript
# Assess the output from the analysis pipeline

setwd(commandArgs(TRUE)[1])
ANALYSIS_DIRECTORY=commandArgs(TRUE)[1]
TAXLEVEL=commandArgs(TRUE)[2]

current_directory <- getwd()





# In which directory are the demultiplexed data folders?
#ANALYSIS_DIRECTORY <- "/Users/rpk/Analysis/PE2x150_JCpoolC_jessieport093014/Analysis_20150202_1517/demultiplexed"

# assign the number of tags in that folder to the object N_TAGS; note you need to move the params file out of the folder before counting...
TAGLIST=system(paste("cd ", ANALYSIS_DIRECTORY, "; ls -d */"), intern=T)
TAGLIST=gsub("/", "", TAGLIST)
N_TAGS <- length(TAGLIST)

# Read in the CSV files called "meganout_mod.csv" from each of those folders, and make them into a list of data frames.
# Storing them as a list makes it easy to analyze them in loops or using the apply functions.
TAGS <- list()
for (i in TAGLIST){
FILE <- paste(ANALYSIS_DIRECTORY, "/", i, "/meganout_",TAXLEVEL,"_mod.csv", sep = "")
tryCatch({
TAGS[[i]] <- read.csv(FILE, header = FALSE, col.names=c("ClusterID","N_Reads","Taxon"))
}, error=function(e){cat("ERROR :",conditionMessage(e), i, "\n")})
}

# create a vector of tag names and apply them to the list of data frames.
# TAG_NAMES <- vector()
# for (i in 1:N_TAGS){
# TAG_NAMES[i] <- paste("tag_", i, sep = "")
# }
# names(TAGS) <- TAG_NAMES

# eliminate ones with no reads (this might include those that ERRONEOUSLY have no reads for now!)
#TAGS[which(sapply(TAGS, nrow) == 0)] <- NULL

N_READS <- list()

for (i in 1:length(TAGS)){
	tryCatch({
		N_READS[[i]] <- sort(sapply(split(TAGS[[i]]$N_Reads, TAGS[[i]]$Taxon), sum))
		}, error=function(e){cat("Note :",conditionMessage(e), i, "\n")})
}

# these are all of the taxon groups
TAXA <- sort(unique(do.call(c, sapply(N_READS, names))))

# Build master dataframe
DATA <- matrix(data=NA, nrow = length(N_READS), ncol=length(TAXA), dimnames=list(names(TAGS), TAXA))

for (i in 1:length(N_READS)){
	if (is.null(N_READS[i][[1]])==FALSE)
	DATA[i,] <- N_READS[[i]][match(TAXA, names(N_READS[[i]]))]
}

# Replace NAs with zeros
DATA[is.na(DATA)] <- 0

write.table(DATA[,names(sort(colSums(DATA), decreasing=T)[1:10])], "/tmp/temp.txt", sep="\t")

# write the data frame
write.csv(DATA, file = paste("~/Desktop/Reads_per_tag_by_OTU_12S_",TAXLEVEL, format(Sys.time(), "%Y%m%d_%X"), ".csv", sep=""))

# samples=read.csv("~/Desktop/20141008_sequencing_sample.csv")
# rownames(DATA)=paste(rownames(DATA), samples$Sample, sep="_")

# #Normalized by NREADS
# DATANORM=DATA[,-which(colnames(DATA)=="No hits")]
# DATANORM=DATANORM[,-which(colnames(DATANORM)=="Not assigned")]
# DATANORM = sweep(DATANORM, 1, rowSums(DATANORM), FUN="/")

# # Examples:
# # Boxplot of number of reads by taxon
# par(mar=c(5,10,1,1))

# # Boxplot of number of Scombridae reads by tag (or sample)
# barplot(DATA[,""], main = "Number of Balanidae reads per tag")


# PCB=DATANORM[c(3,4,7),]
# PCB=round(PCB, 4)
# PCB=PCB[,colSums(PCB)!=0]

# par(mfrow=c(3,1))
# barplot(PCB[1,], horiz=F, las=2)
# barplot(PCB[2,], horiz=F, las=2)
# barplot(PCB[3,], horiz=F, las=2)


# x=round(DATANORM,4)
# x=x[,colSums(x)>0.001]

# barplot(x[i,], las=2, main=rownames(x)[i])
# i=i+1
