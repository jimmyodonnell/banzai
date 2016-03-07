# sequences per run and per library

mydata <- read.table("/Users/threeprime/Desktop/20150717/libraries/quality_reports/fastqc_counts.txt", sep = ",", stringsAsFactors = FALSE)

mydata <- cbind(as.data.frame(do.call(rbind,
        strsplit(x = mydata[,1], split = "_")
  )[,1:3]
)
, mydata[,2])

# add names
colnames(mydata) <- c("library_name", "library", "read", "sequences")

# Are the number of sequences in R1 and R2 always the same?
if(
  identical(
    mydata$sequences[mydata$read == "R1"],
    mydata$sequences[mydata$read == "R2"]
  )
) { 
  print("Same number of reads in both directions, consolidating. The following numbers are PER EACH READ DIRECTION")
  mydata <- mydata[mydata$read == "R1",]
} else { 
  print("Different number of reads in each direction! The following numbers are TOTAL reads!")
}


# how many sequences total?
total_seqs <- sum(mydata$sequences)

# add percentage column
mydata <- data.frame(mydata, percent = mydata$sequences*100/total_seqs)

# Which sequences were assignable to a library?
assignable_seqs <- mydata$sequences[1:14]

# how many sequences were in other study?
other_seqs <- sum(assignable_seqs[1:4])
other_seqs_percent <- other_seqs*100/sum(assignable_seqs)

# present study
present_seqs <- sum(assignable_seqs[5:14])
present_seqs_percent <- present_seqs*100/sum(assignable_seqs)

# jittering is random; to control it, set a seed number
pdf(file = "seqs_per_lib.pdf")
  
set.seed(1)
stripchart(
  ylim = c(0, max(mydata$sequences)),
  x = list(
    mydata$sequences[1:4], 
    mydata$sequences[5:14], 
    mydata$sequences[15]
  ), 
  vertical = TRUE,
  method = "jitter",
  jitter = 0.2,
  pch = 1, col = "black", #bg = "grey",
  cex = 0.8,
  group.names = c("other", "present", "undetermined"), 
  xlab = "Library", 
  ylab = "Number of Sequences", 
  main = "Sequences per library"
)
abline(h = sum(assignable_seqs)/length(assignable_seqs), lty = 3)
box()
dev.off()
# how many undetermined?
