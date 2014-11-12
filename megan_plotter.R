#!usr/bin/Rscript
# Assess the output from MEGAN

setwd(commandArgs(TRUE)[1])

DATA<-read.csv("meganout_mod.csv", header=FALSE)

names(DATA)<-c("ClusterID","N_Reads","Taxon")

N_READS <- sort(sapply(split(DATA$N_Reads, DATA$Taxon), sum))
N_READS_HITS <- N_READS[names(N_READS) != "No hits"]
prop_no_hits <- N_READS[["No hits"]]/sum(DATA["N_Reads"])

pdf(file="N_reads.pdf")
par(mar=c(5,9,1,1))
barplot(N_READS_HITS, hor=TRUE, las=1, xlab="Number of Reads", cex.names=0.8)
text(0,0, paste("proportion no BLAST hits = ", round(prop_no_hits, digits=3), sep=""), adj = c(-0.2,0))
dev.off()
