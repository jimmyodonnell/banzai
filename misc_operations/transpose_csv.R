#!/usr/bin/Rscript

# Transpose a CSV file.

# takes one argument: the CSV file to transpose

# read in arguments from the command line
arguments <- commandArgs(TRUE)

# get the directory containing the file
my_dir <- dirname(arguments)

# get the file name and split it by the dot to get the extension
my_file <- strsplit(basename(arguments), split = '[.]')[[1]]

# get the base file name
my_filename <- my_file[1]

# make an output file name
output_file <- file.path(my_dir, paste(my_filename, "_t.csv", sep = ""))

# read in the file
my_csv <- read.csv(arguments, header = TRUE, row.names = 1)

# transpose it
my_csv_transposed <- t(my_csv)

# write the output
write.csv(x = my_csv_transposed, file = output_file, quote = FALSE)
