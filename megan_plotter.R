#!usr/bin/Rscript
# Assess the output from MEGAN

setwd(commandArgs(TRUE)[1])

current_directory <- getwd()

filelist <- system(paste("find ", current_directory, " -name '*mod.csv'", sep=""), intern = TRUE)

for ( i in 1:length(filelist)){
	DATA<-read.csv(filelist[i], header=FALSE)
	names(DATA)<-c("ClusterID","N_Reads","Taxon")
	
	strsplit(folderlist[1], split='/')
	
	INFILE <- strsplit(folderlist[1], split='/')
	INFILE <- INFILE[[1]][length(INFILE[[1]])]
	
	N_READS <- sort(sapply(split(DATA$N_Reads, DATA$Taxon), sum))
	write.csv(N.READS, file = paste("N_Reads" , i, ".csv", sep = ""))
	N_READS_HITS <- N_READS[names(N_READS) != "No hits"]
	prop_no_hits <- N_READS[["No hits"]]/sum(DATA["N_Reads"])
	
	pdf(file=paste("N_reads", i, ".pdf", sep = ""))
		par(mar=c(5,9,2,2))
		barplot(N_READS_HITS, hor=TRUE, las=1, xlab="Number of Reads", cex.names=0.8, main=)
		text(0,0, paste("proportion no BLAST hits = ", round(prop_no_hits, digits=3), sep=""), adj = c(-0.2,0))
	dev.off()

}

