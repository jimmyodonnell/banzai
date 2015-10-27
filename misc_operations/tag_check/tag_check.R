# plot the efficiency of a tagging protocol

infile <- "/Users/threeprime/Desktop/tag_check/hopkins_tags.txt"

tags <- read.table(file = infile)

pdf(file = "hopkins_tags.pdf")
plot(tags[1:100,1])

dev.off()
