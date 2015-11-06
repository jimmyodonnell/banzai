###################################################
#DECONTAMINATION STEP   #note -- needs tweaking; doesn't quite get all of the suspected contamination out, and I don't know why
###################################################
# Define otu table and metadata file paths
OTU_table_file <- "/Users/jimmy.odonnell/Desktop/Analysis_20151019_1918/all_lib/OTUs_swarm/OTU_table.csv"
metadata_file <- "/Users/jimmy.odonnell/Kelly_Lab_Big/Illumina_Data_Raw/20150717/libraries/kelly_lab/metadata_field_exp.csv"

# Specify project-related values
col_library		<- "library"
col_tagseq		<- "tag_sequence"
col_project <- "experiment"
name_project <- "PS_urban"
col_sampletype	<- "sample_type"
name_control <- "positive_control"# c("positive_control", "negative_control")


# read in OTU table
samples_all <- as.matrix(read.csv(OTU_table_file, row.names = 1))
# samples_all <- OTUs#read.csv("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/Data/Analysis_20150727_0004/OTU_analyses/Swarm_OTU_table.txt", row.names=1)

# read in metadata; set relevant column names
metadata <- read.csv(metadata_file)

# CONTROL=meta[meta$experiment=="control"&grepl("TILAPIA", meta$sample_name),]  #select relevant samples from the metadata
# CONTROL_tags=paste0("lib_", CONTROL$library,"_tag_", CONTROL$tag_sequence)  #generate tag-ids for samples from relevant experiment

# colnames(samples_all)
# colnames(metadata)

sample_id <- paste(
					"lib",
					metadata[,col_library],
					"tag",
					metadata[,col_tagseq],
					sep = "_"
					)

# calculate maximum read number per sample:
N_reads <- max(colSums(samples_all))

samples_non <- samples_all[, ! colnames(samples_all) %in% sample_id]
dim(samples_non)
colSums(samples_non)

samples_targ <- samples_all[, colnames(samples_all) %in% sample_id]
samples_all <- samples_targ

# get the IDs of the control samples
id_control <- sample_id[metadata[, col_sampletype] == name_control ]

# get the IDS of the samples from this project
id_project <- sample_id[metadata[, col_project] == name_project ]


# break data into control and experimental samples
samples_env <- samples_all[ , id_project]
# samples_env <- samples_all[,match(PSurban_tags,names(samples_all))]
samples_con <- samples_all[ , id_control]
# samples_con <- samples_all[,match(CONTROL_tags,names(samples_all))]  #here, just using positive (Tilapia) controls, rather than the DIH2O controls

# IDENTIFY CONTROL TAXON
# Which OTU is (must be a row name in the samples_all matrix)
# control_taxon <- "OTU_1"
control_taxon <- names(which.max(rowSums(samples_con)))

# make vector for coloring by sample type
color_by_type <- metadata[
						match(
							colnames(samples_all), 
							sample_id
							),
						col_sampletype
						]
# color_by_type <- as.character(color_by_type)

# color_by_type[which(color_by_type == "environmental")] <- "black"
# color_by_type[which(color_by_type == levels(metadata[, col_sampletype])[1])] <- "black"
# color_by_type[which(color_by_type == levels(metadata[, col_sampletype])[2])] <- "pink"
# color_by_type[which(color_by_type == levels(metadata[, col_sampletype])[3])] <- "red"
color_by_type[which(is.na(color_by_type))] <- "white"

plot(
	x = samples_all[control_taxon,]/colSums(samples_all), 
	col = as.numeric(color_by_type), 
	xlab = "Sample", 
	ylab = "Proportion reads positive control taxon"
	)

legend(
	x = 0, 
	y = 0.5, 
	legend = levels(color_by_type), 
	col = c(1,2,3), 
	# as.numeric(color_by_type), 
	pch = 1
	)


plot(colSums(samples_all), col = "grey")
points(
	x = samples_all[control_taxon,], 
	col = as.numeric(color_by_type)
	)
abline(v = 1:ncol(samples_all), col = rgb(0,0,0, alpha = 0.1))


################################################################################
# calculate two rates of contamination:
################################################################################
# 1. non-control taxa in control samples
contam_con <- samples_con[!(rownames(samples_all) %in% control_taxon),]

# 1a. max abundance of any single non-target taxon in control samples
contam_tax_max_raw <- max(contam_con)

# 1b. max proportion of any single non-target taxon in control samples
contam_tax_prop <- apply(X = contam_con, MARGIN = 2, FUN = max) / colSums(samples_con)
contam_tax_prop_max <- max(contam_tax_prop)

# 1c. max cumulative abundance of non-target taxa in control samples
contam_sum_raw <- colSums(contam_con)
contam_sum_raw_max <- max(contam_sum_raw)

# 1d. max cumulative proportion of non-target taxa in control samples 
contam_sum_prop <- colSums(contam_con) / colSums(samples_con)
contam_sum_prop_max <- max(contam_sum_prop)

# plot
# stripchart(contam_sum_prop, pch = 1, method = "jitter", jitter = 20)

# which taxa are contaminants?
contam_taxa_in_control <- names(rowSums(contam_con)[rowSums(contam_con) > 0])


################################################################################
# 2. control taxon in non-control samples
contam_env <- samples_env[control_taxon,]

# 2a. max abundance of control taxon in environmental samples
contam_env_raw_max <- max(contam_env)

# 2b. proportion of control taxon in environmental samples
contam_env_prop <- contam_env / colSums(samples_env)
contam_env_prop_max <- max(contam_env_prop)

# plot (figname control_taxon_in_env_samples )
stripchart(
	x = contam_env_prop, 
	xlim = c(0,1), 
	xlab = "proportion control taxon in environmental samples", 
	pch = 1, 
	method = "jitter"
	)
	
summary(contam_env_prop)

################################################################################


# THUS, you could run the script "filter_otu_abundance", using as input the number:
min_otu_percent <- contam_tax_prop_max*100
infile <- OTU_table_file
otu_table <- samples_env




# Jimmy Stopped here 2015 11 05
################################################################################








# CALCULATE MAXIMUM ABUNDANCE OF NON-TARGET OTUS IN CONTROL
contam_max <- max(samples_con[!(rownames(samples_all) %in% control_taxon),])
contam_sum_max <- max(colSums(samples_con[!(rownames(samples_all) %in% control_taxon),]))

contam_max_prop <- as.numeric(max(samples_con[!(rownames(samples_all) %in% control_taxon),])/colSums(samples_all)[which.max(sapply(samples_con[!(rownames(samples_all) %in% control_taxon),], max))])  #same thing, but as a proportion
contam_sum_max_prop <- max(colSums(samples_con[!(rownames(samples_all) %in% control_taxon),])/colSums(samples_con))


#Step 1: remove control taxon OTU from environmental samples in which it occurs
control_removed <- samples_all[(rownames(samples_all) %in% control_taxon),]
samples_no_control <- samples_all[!(rownames(samples_all) %in% control_taxon),]

#Step 2: calculate cumulative proportional abundance of least-frequent OTUs in each sample, and truncate at the threshold calculated above as contam_sum_max_prop

samples_props= scale(samples_no_control, center=FALSE, scale=colSums(samples_no_control))
samples_contam_prop_sum <- samples_props
smallest=list()
for(i in 1:ncol(samples_props)){
  	namevec=names(samples_props[order(samples_props[,i], decreasing = FALSE), i])
	smallest[[i]]=namevec[cumsum(samples_props[order(samples_props[,i], decreasing = FALSE), i])<contam_sum_max_prop]
  	}

for (i in 1:ncol(samples_contam_prop_sum)){
	samples_contam_prop_sum[smallest[[i]],i] <- 0
}

#sum(samples_props[,c(1)]>0) ; sum(samples_contam_prop_sum[,1]>0)  #view the change to an example sample

#Step 3: re-create counts from proportion data, given new colSums
newCountTotals=NA
for (i in 1:ncol(samples_no_control)){
	newCountTotals[i]= sum(samples_no_control[-which(row.names(samples_no_control)%in%smallest[[i]]),i])
}

DECONTAM=round(scale(samples_contam_prop_sum, center=FALSE, scale=1/newCountTotals), 0)
