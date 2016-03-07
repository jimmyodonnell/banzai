setwd("~/Documents/Data/IlluminaData/16S/20141020/Analysis_20141023_1328/Reanalysis_20141028_1351")
tag1 <- read.csv("tag_1/meganout_mod.csv", header=FALSE)
tag8 <- read.csv("tag_8/meganout_mod.csv", header=FALSE)

DATA <- tag8


names(DATA)<-c("ClusterID","N_Reads","Taxon")

N_READS <- sort(sapply(split(DATA$N_Reads, DATA$Taxon), sum))
N_READS_HITS <- N_READS[names(N_READS) != "No hits"]
prop_no_hits <- N_READS[["No hits"]]/sum(DATA["N_Reads"])

par(mar=c(5,9,1,1))
barplot(N_READS_HITS, hor=TRUE, las=1, xlab="Number of Reads", cex.names=0.8)
text(0,0, paste("proportion no BLAST hits = ", round(prop_no_hits, digits=3), sep=""), adj = c(-0.2,0))
