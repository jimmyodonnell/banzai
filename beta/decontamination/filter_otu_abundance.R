#!/usr/bin/env Rscript

# Remove OTUs lower than a given threshold

# Usage: (from terminal)
# Rscript this_script_name.R 0.005 /path/to/OTU_table.csv

# expects two arguments: 
# 1: the threshold percent below which to disregard OTUs
# 2: a file path to the OTU table to analyze

# read in a vector of arguments
arguments <- commandArgs(TRUE)


# From Nick Bokulich: 0.005% filtered on a per-sequencing run basis (as opposed to per-sample)
min_otu_percent <- as.numeric(arguments[1])
# min_otu_percent <- 0.005

min_otu_prop <- min_otu_percent / 100

# define input
infile <- arguments[2]
# infile <- "/Users/jimmy.odonnell/Downloads/Swarm_OTU_table.txt"

# read OTU table; expects rows=OTUs, cols=samples
otu_table <- as.matrix(read.csv(infile, row.names = 1))

# transpose OTU table
otu_table <- t(otu_table)

# assemble path for output directory
start_time <- gsub(
				pattern = ":|-| ", 
				replacement = "", 
				x = as.character(Sys.time())
				)

out_dir <- file.path(
				dirname(infile), 
				paste(
					"OTU_table_filtered", 
					start_time, 
					sep = "_"
					)
				)
system(paste("mkdir", out_dir))

# row/column operatio (output is vector)
# calculate total reads
total_reads <- sum(otu_table)

# calculate total OTUs
N_otus_in <- ncol(otu_table)

# calculate number of samples
N_samples <- nrow(otu_table)

# calculate number of OTUs needed to pass dataset wide filter
min_reads_dataset <- min_otu_prop * total_reads





# row/column operations (output is vector)
# calculate reads per OTU
reads_per_otu <- colSums(otu_table)

# calculate reads per sample
reads_per_sample <- rowSums(otu_table)
reads_per_sample_min <- min(reads_per_sample)
reads_per_sample_mean <- mean(reads_per_sample)
reads_per_sample_max <- max(reads_per_sample)

# calculate minimum reads per otu per sample
threshold_per_sample <- reads_per_sample * min_otu_prop
threshold_per_sample_min <- min(threshold_per_sample)
threshold_per_sample_mean <- mean(threshold_per_sample)
threshold_per_sample_max <- max(threshold_per_sample)





# Table-wide calculations (output is matrix)
# calculate proportional abundance of each OTU per sample
otu_prop_samp <- otu_table/reads_per_sample

# calculate proportional abundance of each OTU per dataset
otu_prop_whole <- otu_table/total_reads





#############################################################################################
# FILTER BASED ON ENTIRE DATASET
#############################################################################################
# Disregard OTUs whose overall abundance in the dataset is less than 0.005 % 
# For each OTU, is it greater than the minimum threshold for the overall dataset?
OTUs_pass_dataset <- reads_per_otu > min_reads_dataset

# remove OTUs that are removed from all samples:
otu_mat_filtered_dataset <- otu_table[, OTUs_pass_dataset ]

# otu_mat_filtered_whole <- otu_prop_whole > min_otu_prop
# otu_table_filtered_whole <- otu_table
# otu_table_filtered_whole[! otu_mat_filtered_whole] <- 0
# otu_table_stripped_whole <- otu_table_filtered_whole[, colSums(otu_table_filtered_whole) > 0]

# calculate how many OTUs were retained
N_otus_pass_dataset <- sum(OTUs_pass_dataset)
prop_otus_pass_dataset <- N_otus_pass_dataset/N_otus_in

# calculate how many reads were retained
reads_retained_dataset <- sum(otu_mat_filtered_dataset)
prop_retained_dataset <- reads_retained_dataset/total_reads

write.csv(
	x = otu_mat_filtered_dataset, 
	file = file.path(out_dir, "otu_table_filtered_per_run.csv"), 
	quote = FALSE
)




#############################################################################################
# FILTER ON PER-SAMPLE BASIS
#############################################################################################
# Disregard OTUs whose abundance in a given sample is less than 0.005 % of the reads in that sample
otu_filter_mat_samp <- otu_prop_samp > min_otu_prop
otu_mat_filtered_samp <- otu_table
otu_mat_filtered_samp[!otu_filter_mat_samp] <- 0

# remove OTUs that are removed from all samples:
otu_table_stripped_samp <- otu_mat_filtered_samp[, colSums(otu_mat_filtered_samp) > 0]

OTUs_pass_samp <- colnames(otu_table_stripped_samp)

# calculate how many OTUs were retained
N_otus_pass_samp <- ncol(otu_table_stripped_samp)
prop_otus_pass_samp <- N_otus_pass_samp/N_otus_in

# calculate how many reads were retained
reads_retained_sample <- sum(otu_table_stripped_samp)
prop_retained_sample <- reads_retained_sample/total_reads

write.csv(
	x = otu_table_stripped_samp, 
	file = file.path(out_dir, "otu_table_filtered_per_samp.csv"), 
	quote = FALSE
)


#############################################################################################
# link together stuff to write a log file:
# ex c("some thing" = "atta boy")
output_vec <- c(
	"OTU table file:" = infile, 
	"output directory:" = out_dir, 
	"minimum OTU abundance threshold (%):" = min_otu_percent, 
	"minimum OTU abundance threshold (prop):" = min_otu_prop, 
	"number of reads:" = total_reads, 
	"number of samples:" = N_samples, 
	"OTUs input:" = N_otus_in, 
	"reads per OTU (mean):" = mean(reads_per_otu), 
	"reads per sample (min):" = reads_per_sample_min,
	"reads per sample (mean):" = reads_per_sample_mean,
	"reads per sample (max):" = reads_per_sample_max,
	"minimum OTU abundance threshold per dataset:" = min_reads_dataset, 
	"minimum OTU abundance threshold per sample (min):" = threshold_per_sample_min, 
	"minimum OTU abundance threshold per sample (mean):" = threshold_per_sample_mean, 
	"minimum OTU abundance threshold per sample (max):" = threshold_per_sample_max, 
	"OTUs above threshold (N; dataset):" = N_otus_pass_dataset,
	"OTUs above threshold (prop; dataset):" = prop_otus_pass_dataset, 
	"reads above threshold (N; dataset):" = reads_retained_dataset, 
	"reads above threshold (prop; dataset):" = prop_retained_dataset, 
	"OTUs above threshold (N; per sample):" = N_otus_pass_samp,
	"OTUs above threshold (prop; per sample):" = prop_otus_pass_samp, 
	"reads above threshold (N; per sample):" = reads_retained_sample,
	"reads above threshold (prop; per sample):" = prop_retained_sample
)

# data.frame(variable_name = names(output_vec), 

output_df <- data.frame(value = output_vec)
write.table(
	x = output_df, 
	file = file.path(out_dir, "summary.txt"), 
	quote = FALSE, 
	col.names = FALSE
	)


# plots
create_plots <- TRUE

if( create_plots ){
	# boxplot(min_otu_prop * reads_per_sample)
	
	pdf(file = file.path(out_dir, "plots.pdf"))
	plot(
		x = cumsum(sort(reads_per_otu/total_reads, decreasing = TRUE)), 
		type = "l", 
		lty = 1, 
		ylab = "cumulative proportion of total reads", 
		ylim = c(0,1), 
		log = "x", 
		xlab = "OTU ranked by abundance (log)"
		)
	abline(v = c(N_otus_pass_dataset, N_otus_pass_samp), lty = c(2, 4), col = c(1, 2))
	
	legend(
		"bottomright", 
		legend = c("dataset filter", "sample filter (mean)"), 
		lty = c(2, 4), 
		col = c(1,2), 
		bg = "white"
		)
	
	
	plot(
		x = sort(colSums(otu_table), 
		decreasing = TRUE), 
		log = "y", 
		col = rgb(0,0,0, alpha = 0.2), 
		pch = 20, 
		ylab = "number of reads", 
		xlab = "OTU ranked by abundance"
		)
	abline(h = c(1:10, 100, 1000), lty = 2, col = "lightgrey")
	abline(h = c(min_reads_dataset, threshold_per_sample_mean), lty = c(2, 4), col = c(2, 4))
	legend(
		"topright", 
		legend = c("dataset filter", "sample filter (mean)"), 
		lty = c(2, 4), 
		col = c(2, 4), 
		bg = "white"
		)
	
	dev.off()

}



#    GRAVEYARD

# test with fake data:
# fake_data <- matrix(
				# data = c(100, 50, 25, 10, 5, 2, 50, 25, 10, 5, 2, 0, 25, 10, 5, 2, 0, 0, 10, 5, 2, 0, 0, 0), 
				# ncol = 4, 
				# byrow = FALSE, 
				# dimnames = list(
								# c(paste("otu", seq(1, 6), sep = "")),
								# c(paste("sample", seq(1, 4), sep = ""))
								# )
				# )

# otu_table <- fake_data
# min_otu_percent <- 3


# text(
	# x = c(N_otus_pass_dataset, N_otus_pass_samp), 
	# y = c(0,0), 
	# labels = c("dataset filter", "sample filter"), 
	# pos = 4, 
	# cex = 1, 
	# srt = 90
	# )

