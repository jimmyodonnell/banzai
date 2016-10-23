SAMPLES <- read.csv("~/Documents/Dropbox/Kelly_Lab/Lab Notebook/20141008_sequencing_sample.csv")
SAMPLES <- SAMPLES[1:8,]
DATA <- DATA[1:8,]
rownames(DATA) <- SAMPLES[,"Sample"]
# DATA <- cbind(SAMPLES["Sample"], DATA)

DATA <- DATA[,which(colSums(DATA)!=0)]
NoHits <- DATA[,"No hits"]
NotAssigned <- DATA[,"Not assigned"]
DATA.raw <- DATA
DATA <- DATA[,-which(colnames(DATA)=="No hits")]
DATA <- DATA[,-which(colnames(DATA)=="Not assigned")]

DATA <- DATA/rowSums(DATA)

par(mar=c(5,10,1,1))
barplot(rowSums(DATA), las=1, horiz=TRUE)

plot(SAMPLES[,"concentration_nguL"], rowSums(DATA))

library(vegan)
diversity(DATA)
barplot(diversity(DATA), horiz=TRUE, las=1)

COL="Clinus "; barplot(DATA[,COL]/rowSums(DATA), horiz=TRUE, las=1, xlab=paste("Proportion reads ", COL, sep=""))

for (COL in colnames(DATA)){
	barplot(DATA[,COL]/rowSums(DATA), horiz=TRUE, las=1, xlab=paste("Proportion reads ", COL, sep=""), )
}


