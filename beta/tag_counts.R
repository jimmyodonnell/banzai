
tag_count_file <- "/Users/jimmy.odonnell/Desktop/Analysis_20151013_1719/tag_count.txt"
tag_counts <- read.table(file = tag_count_file,
                         header = TRUE, 
                         sep = " "
                         )
head(tag_counts)

# All combinations of primary and secondary index are considered.
# omit combinations which have a very low number of reads relative to other samples
# samples with fewer than this proportion of the mean number of reads will be excluded
lower_percent_threshold <- 0.05


low_frequency_samples <- which(tag_counts$left_tagged < mean(tag_counts$left_tagged)*lower_percent_threshold)

low_frequency_data <- tag_counts[low_frequency_samples, ]

tag_counts <- tag_counts[-low_frequency_samples,]



plot(
  sort(tag_counts$left_tagged)
  )

plot(
  sort(
    (tag_counts$left_tagged - tag_counts$right_tagged) / tag_counts$left_tagged
    )
  )


tag_rate <- tag_counts[, "right_tagged"] / tag_counts[,"left_tagged"]
boxplot(tag_rate, ylim = c(0, 1))

mean(tag_rate)
sd(tag_rate)
range(tag_rate)
nrow(tag_counts)
