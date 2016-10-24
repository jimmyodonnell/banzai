#!/bin/R

INFILE <- '/Users/threeprime/Desktop/outfile1.txt'

READ_LENGTH <- read.table(INFILE)[,1]

# Find the mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

MIN <- paste("min : ", min(READ_LENGTH), sep = "")
MEAN <- paste("mean : ", round(mean(READ_LENGTH), digits = 2), sep = "")
MAX <- paste("max : ", max(READ_LENGTH), sep = "")
MODE <- paste("mode : ", Mode(READ_LENGTH), sep = "")




# pdf(file = "read_length_hist.pdf")
hist(READ_LENGTH, breaks = 100)
# add text using function legend()
legend("topright", legend = c(MIN, MEAN, MAX, MODE), bty = "n", lty = c(0,2,3,0))
abline(v = c(mean(READ_LENGTH), Mode(READ_LENGTH)), lty = c(2,3))
# alternatively, add text using function text()
# PLOT_COORDS <- par('usr')
# PLOT_TEXT <- paste(MIN, MEAN, MAX, MODE, sep = "\n")
# text(x = PLOT_COORDS[2]*0.9, y = PLOT_COORDS[4]*0.9, labels = PLOT_TEXT)
# dev.off()
