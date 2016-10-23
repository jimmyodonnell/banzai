#####################################################################################
#7) #####DESeq2 Analyses ################################################
#####################################################################################

library(DESeq2)

#Load DESeq2 Datasets
edna_deseq_count=read.csv("~/Desktop/eDNA_Project/12S_Hopkins_5May15_deseq_countData.csv", sep=",", header=FALSE, row.names=1) #countData datasheet (rows = taxa and columns = samples)
edna_deseq_col=read.csv("~/Desktop/eDNA_Project/12S_Hopkins_filtered_dataset_19Nov14_deseq_colData.csv", sep=",", header=TRUE) #colData datasheet (rows must correspond to columns in countData)

#Create DESeqDataSet object
edna_deseq=DESeqDataSetFromMatrix(edna_deseq_count, edna_deseq_col, design=~Habitat_bins)

#Estimate size factors when zeros are present in each row
fordeseq=counts(edna_deseq)
geoMeans = apply(fordeseq, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
edna_deseq= estimateSizeFactors(edna_deseq, geoMeans=geoMeans)

deseq=counts(edna_deseq, normalized=TRUE) #store normalized counts table
write.csv(deseq, "~/Downloads/deseq_scaled_5May15.csv") #Write out normalized counts


