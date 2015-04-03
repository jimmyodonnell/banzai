#!/usr/bin/Rscript
# Assess the output from MEGAN

args<-commandArgs(TRUE)

setwd(args[1])

current_directory <- getwd()

filelist <- system(paste("find ", current_directory, " -name '*mod.csv'", sep=""), intern = TRUE)

for ( i in 1:length(filelist)){
	DATA<-read.csv(filelist[i], header=FALSE)
	names(DATA)<-c("ClusterID","N_Reads","Taxon")

	N_READS <- sort(sapply(split(DATA$N_Reads, DATA$Taxon), sum))
	write.csv(N_READS, file = paste("N_Reads" , basename(filelist[i]), sep = ""))
	N_READS_HITS <- N_READS[names(N_READS) != "No hits"]
	ifelse((sum(colnames(N_READS)=="No hits")==0), prop_no_hits <- 0 ,prop_no_hits <- N_READS[["No hits"]]/sum(DATA["N_Reads"]))

	pdf(file=paste("N_reads_", strsplit(basename(filelist[i]), "\\.")[[1]][1], ".pdf", sep = ""))
		par(mar=c(5,9,2,2))
		barplot(N_READS_HITS, hor=TRUE, las=1, xlab="Number of Reads", cex.names=0.8)
		text(0,0, paste("proportion no BLAST hits = ", round(prop_no_hits, digits=3), sep=""), adj = c(-0.2,0))
	dev.off()

}
