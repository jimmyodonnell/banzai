# Remove spurious (false negative, false positive) sequences from environmental samples based on sequences in the control samples

# Controls: Samples derived from a single, known source (e.g. tissue).
# Assumption: Reads in the control sample that are NOT from the species of origin are 'spurious' - contamination in some form.

# There are several ways to deal with this:

# 1. Calculate the maximum abundance of any one *non-target* OTU in the controls, and
  # A. remove from the environmental samples any OTU with equal or lesser abundance.
  # B. remove from the environmental samples this amount from EVERY OTU.
# 2. Calculate the maximum abundance for *every* OTU in the controls, and remove this amount from any other samples in which those OTUs occur.
# 3. Calculate the sum of the abundance of non-target OTUs in the controls, sort the OTUs in the environmental samples by abundance, and remove from the environmental samples the OTUs whose sum is less than this value.
# 4. Set to zero in the environmental samples any OTU that is also found in the control samples 

# ALSO,
# 5. Remove any reads from the target OTU used in the control.

pdf(file = "remove_contamination.pdf")
# PROLOGUE
# either load your data:
# samples_all <- read.csv(file = "OTU_table.csv")
# isolate the columns from control and environmental samples like so:
# samples_env <- samples_all[,c(1,2,3)]
# samples_con <- samples_all[,c(4,5,6)]
# Where:
# 1, 2, and 3 are the numbers of the columns containing environmental samples, and
# 4, 5, and 6 are the numbers of the columns containing control samples

# OR
# Generate some data

N_OTU <- 10 # number of OTUs
N_env <- 3 # number of environmental samples
N_control <- 3 # number of control samples
N_reads <- 10000 # number of reads per sample

library(gtools) # rdirichlet

set.seed(2)
samples_env <- t(round(x = rdirichlet(n = N_env, alpha = rexp(N_OTU))*N_reads, digits = 0))
dimnames(samples_env) <- list(
  paste("OTU_", 1:N_OTU, sep = ""), 
  paste("E", 1:N_env, sep =""))
# barplot(samples_env)

probs <- c(0.9975, 5e-04, 4e-04, 3e-04, 2e-04, 1e-04)
probs_full <- c(probs, rep(x = 0, times = N_OTU - length(probs)))

samples_con <- replicate(N_control,
  table(
    factor(
      sample(
        x = 1:N_OTU, 
        size = N_reads, 
        replace = TRUE, 
        prob = probs_full
        ), 
      levels = 1:N_OTU
      )
    )
  )
dimnames(samples_con) <- list(
  paste("OTU_", 1:N_OTU, sep = ""), 
  paste("C", 1:N_control, sep ="")
)

# barplot(samples_con)

samples_all <- cbind(samples_env, samples_con)
barplot(samples_all, ylab = "number of reads", main = "Full Dataset")


################################################################################
# IDENTIFY CONTROL TAXON
# Which OTU is (must be a row name in the samples_all matrix)
# control_taxon <- "OTU_1"
control_taxon <- names(which.max(rowSums(samples_con)))


################################################################################
# CALCULATE MAXIMUM ABUNDANCE OF NON-TARGET OTUS IN CONTROL
contam_max <- max(samples_con[!(rownames(samples_all) %in% control_taxon),])
contam_sum_max <- max(colSums(samples_con[!(rownames(samples_all) %in% control_taxon),]))


################################################################################
# 1A. remove from the environmental samples any OTU with equal or lesser abundance.
# Set to zero any OTU occurrence below the maximum contaminating OTU abundance
samples_threshold <- samples_all
samples_threshold[samples_threshold <= contam_max] <- 0
barplot(samples_threshold, ylab = "number of reads", main = "Remove OTUs rarer than most abundant contamination")

################################################################################
# 1B. remove from the environmental samples this amount from EVERY OTU.
samples_contam_removed <- samples_all - contam_max
barplot(samples_contam_removed, ylab = "number of reads", main = "Every OTU trimmed of contaminant proportion")

################################################################################
# 2. Calculate the maximum abundance for *every* OTU in the controls, and remove this amount from this OTU in any other samples in which it occurs.
contam_max_OTU <- apply(X = samples_con, MARGIN = 1, FUN = max)
samples_OTU_subtract <- samples_all - contam_max_OTU
samples_OTU_subtract[samples_OTU_subtract < 0] <- 0
barplot(samples_OTU_subtract, ylab = "number of reads", main = "OTU-specific contamination removed")


################################################################################
# 3. Calculate the sum of the abundance of non-target OTUs in the controls, sort the OTUs in the environmental samples by abundance, and remove from the environmental samples the OTUs whose sum is less than this value.
# Set to zero the OTUs comprising the proportion of each sample less than or equal to the proportion of contamination in the controls

samples_contam_sum <- samples_all
for(i in 1:ncol(samples_all)){
  smallest <- names(which(
      cumsum(
        samples_all[
          order(samples_all[,i], decreasing = FALSE), i
        ]) <= contam_sum_max
      )
  )
  samples_contam_sum[smallest,i] <- 0
}
barplot(samples_contam_sum, ylab = "number of reads", main = "Remove tail equal to proportion of contamination in control")


################################################################################
# 4. Set to zero in the environmental samples any OTU that is also found in the control samples 
OTUs_contam <- names(which(rowSums(samples_con) > 0))
samples_contam_OTU <- samples_all
samples_contam_OTU[OTUs_contam,] <- 0
barplot(samples_contam_OTU, ylab = "number of reads", main = "Remove all OTUs found in controls")


################################################################################
# 5. Remove any reads from the target OTU used in the control.
# REMOVE CONTROL TAXON FROM ALL SAMPLES
control_removed <- samples_all[(rownames(samples_all) %in% control_taxon),]
samples_no_control <- samples_all[!(rownames(samples_all) %in% control_taxon),]
barplot(samples_no_control, ylab = "number of reads", main = "Remove control taxon")

dev.off()