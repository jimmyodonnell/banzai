grep_trimmed <- read.csv("~/Desktop/grep_trimmed.csv", stringsAsFactors=FALSE)
cutadapt_trimmed <-read.csv("~/Desktop/cutadapt_trimmed.csv")

total_reads <- 2282381
grep_trimmed <- grep_trimmed[2:nrow(grep_trimmed),]
sum(grep_trimmed[,2])/total_reads
prop <- grep_trimmed[,2]/total_reads
grep_trimmed <- cbind(grep_trimmed, prop)
sum(grep_trimmed[,3][1:4])
sum(grep_trimmed[,3][5:8])

