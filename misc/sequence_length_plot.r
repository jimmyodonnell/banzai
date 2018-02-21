#!/usr/bin/env Rscript

# INFILE <- commandArgs(trailingOnly = FALSE)[1]
# INFILE <- '/Users/threeprime/Desktop/outfile1.txt'
INFILE <- 'AMCPMOO-101-1_S5_L001_R1_001.seqlen'

READ_LENGTH <- read.table(INFILE, header = TRUE)

head(READ_LENGTH)



MIN <- paste("min : ", min(READ_LENGTH[,'length']), sep = "")
MEAN <- paste("mean : ", 
  round(mean(with(READ_LENGTH, rep(length, times = count))), 
  digits = 2), sep = "")
MAX <- paste("max : ", max(READ_LENGTH[,'length']), sep = "")
MODE <- paste("mode : ", 
  READ_LENGTH[which.max(READ_LENGTH[,'count']),'length'], sep = "")


# pdf(file = "read_length_hist.pdf")
with(READ_LENGTH, {
  plot(length, count, las = 1, type = 'n', ylab = "Count", xlab = "Read Length")
  grid()
  points(length, count, las = 1, type = 'b')
})

# add text using function legend()
legend("topleft", legend = c(MIN, MEAN, MAX, MODE), bty = "b") #, lty = c(0,2,3,0)
# abline(h = c(mean(READ_LENGTH), Mode(READ_LENGTH)), lty = c(2,3))
# alternatively, add text using function text()
# PLOT_COORDS <- par('usr')
# PLOT_TEXT <- paste(MIN, MEAN, MAX, MODE, sep = "\n")
# text(x = PLOT_COORDS[2]*0.9, y = PLOT_COORDS[4]*0.9, labels = PLOT_TEXT)
# dev.off()
