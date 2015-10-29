###################################################
#DECONTAMINATION STEP   #note -- needs tweaking; doesn't quite get all of the suspected contamination out, and I don't know why
###################################################


samples_all <- as.matrix(read.csv("/Users/jimmy.odonnell/Downloads/Swarm_OTU_table.txt", row.names = 1))
# samples_all <- OTUs#read.csv("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/PS_urban_eDNA/Data/Analysis_20150727_0004/OTU_analyses/Swarm_OTU_table.txt", row.names=1)

# calculate maximum read number per sample:
N_reads <- max(colSums(samples_all))

metadata <- read.csv("/Users/jimmy.odonnell/Kelly_Lab_Big/Illumina_Data_Raw/20150717/libraries/kelly_lab/metadata_field_exp.csv")
col_library		<- "library"
col_tagseq		<- "tag_sequence"
col_project <- "experiment"
name_project <- "PS_urban"
col_sampletype	<- "sample_type"
name_control <- "positive_control"# c("positive_control", "negative_control")

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

# get the IDs of the control samples
id_control <- sample_id[metadata[, col_sampletype] == name_control ]

# get the IDS of the samples from this project
id_project <- sample_id[metadata[, col_project] == name_project ]





	samples_env <- samples_all[ , id_project]
	# samples_env <- samples_all[,match(PSurban_tags,names(samples_all))]
	samples_con <- samples_all[ , id_control]
	# samples_con <- samples_all[,match(CONTROL_tags,names(samples_all))]  #here, just using positive (Tilapia) controls, rather than the DIH2O controls

	# IDENTIFY CONTROL TAXON
	# Which OTU is (must be a row name in the samples_all matrix)
	# control_taxon <- "OTU_1"
	control_taxon <- names(which.max(rowSums(samples_con)))

	
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

