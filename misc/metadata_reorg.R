#!/usr/bin/env R

getwd()

metadata_file <- "metadata.csv"

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

metadata <- metadata[with(metadata, order(sample_name, library, tag_number)),]

temp <- split(metadata$sample_name, metadata$sample_name)

sample_id <- unlist(
  sapply(temp,
    function(x) paste(x, 1:length(x), sep = ".")
  )
)

# check I haven't screwed up the order
temp2 <- sapply(strsplit(sample_id, split = "\\."), "[", 1)
names(temp2) <- NULL
identical(temp2, metadata$sample_name)


metadata$sample_type[metadata$sample_type == "positive_control"] <- "tissue"

dna_source <- metadata$sample_type
dna_source[dna_source == "positive_control"] <- "tissue"
dna_source[dna_source == "environmental"] <- "filter-env"
dna_source[dna_source == "negative_control"] <- "filter-dih2o"
dna_source

# desired fields:
metadata_new <- data.frame(

lab_order       = metadata$order,

dna_id          = metadata$sample_name,
dna_source      = dna_source,

sample_id       = sample_id,

pri_index_name  = metadata$library,
pri_index_seq   = metadata$primary_index_seq,

sec_index_start = rep(4, nrow(metadata)),
sec_index_name  = metadata$tag_number,
sec_index_seq   = metadata$tag_sequence,

index_combo     = metadata$library_tag_combo,

insert_size     = rep(180, nrow(metadata)),

locus           = metadata$locus,
primerF_seq     = metadata$primer_sequence_F,
primerR_seq     = metadata$primer_sequence_R,

file1           = metadata$file1,
file2           = metadata$file2

)

write.csv(
  metadata_new, 
  file = "metadata_new.csv", 
  row.names = FALSE, 
  quote = FALSE
)

# "order"
# "sample_name",
# "library_tag_combo",
# "library",
# "primary_index_seq",
# "sample_type",
# "locus",
# "primer_sequence_F",
# "primer_sequence_R",
# "tag_number",
# "tag_sequence",
# "file1",
# "file2"


# exclude: "date_PCR2",
