
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

# low_frequency_samples <- which(tag_counts$left_tagged < mean(tag_counts$left_tagged)*lower_percent_threshold)
# low_frequency_data <- tag_counts[low_frequency_samples, ]
# tag_counts <- tag_counts[-low_frequency_samples, ]

boxplot(split(tag_counts$left_tagged, tag_counts$library))

plot(sort(tag_counts$right_tagged_by_max_v))

plot(diff(sort(tag_counts$right_tagged_by_max_v)))
# BINGO!




left_tagged_by_max <- lapply(split(tag_counts$left_tagged, tag_counts$library), function(x) x/max(x))
right_tagged_by_max <- lapply(split(tag_counts$right_tagged, tag_counts$library), function(x) x/max(x))

k_membership_l <- unlist(lapply(left_tagged_by_max, function(x) kmeans(x, centers = 2)$cluster))
k_membership_r <- unlist(lapply(right_tagged_by_max, function(x) kmeans(x, centers = 2)$cluster))

left_tagged_by_max_v <- unlist(left_tagged_by_max)
right_tagged_by_max_v <- unlist(right_tagged_by_max)

tag_counts <- data.frame(tag_counts, left_tagged_by_max_v, right_tagged_by_max_v, k_membership_l, k_membership_r, stringsAsFactors = FALSE)
length(k_membership_l)
nrow(tag_counts)

split(right_tagged_by_max_v, c(tag_counts$library, k_membership_r))
membership <- 
tag_counts$right_tagged[k_membership_r == 1,]
split(, )

{
}

abundance_by_lib_k <- list()
for(i in 1:max(k_membership_r)){
	abundance_by_lib_k[[i]] <- 	lapply(split(tag_counts, tag_counts[,"library"]), function(x) x[x[,"k_membership_r"] == i, "right_tagged_by_max_v"])
}



stripchart(
	abundance_by_lib_k[[1]], 
	las = 1, 
	method = "jitter", 
	pch = 21, 
	col = 3, 
	ylab = "library"
	)
stripchart(
	abundance_by_lib_k[[2]], 
	las = 1, 
	ylab = "library", 
	method = "jitter", 
	pch = 21, 
	col = 4, 
	add = TRUE
	)



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
