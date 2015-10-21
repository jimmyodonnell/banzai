Sequencing Metadata
=====================

Metadata must be stored as a CSV file that meets the following requirements:  

- Unicode text encoding
- Unix ("LF") newlines/line breaks
- There are no SPACES, COLONS, COMMAS, or PARENTHESES

You can check that your file doesn’t contain any weirdness by opening it with a plain text editor like TextWrangler (in TextWrangler, the bottom bar should display the current text encoding and newline/linebreaks).

*NOTE*: Every time you open the file in Excel for Mac (even if you don’t change anything), it will change newline characters to Classic Mac (CR). Every time you open the file in LibreOffice Calc, it will add quotes.



### Columns that should be included in sequencing metadata
- sample_id : think of this as the sample with respect to sequencing, not environmental
- locus
- library_index
- library_name
- primer1_name
- primer1_sequence_core
- primer1_index_name
- primer1_index_sequence
- primer1_sequence_full
- primer2_name
- primer2_sequence_core
- primer2_index_name
- primer2_index_sequence
- interprimer_size_min
- interprimer_size_max
- interprimer_size_mean
- env_sample_name
- DNA_extraction_name
- PCR_name
- project : useful for combining multiple projects per run (e.g. grid, spatial)

environmental data:

- env_sample_name
- date_collected
- time_collected
- latitude
- longitude
- ?tidal_height
- ?salinity
- ?temp
- ?habitat
- ?distance_from_shore


from an existing file:
order
sample_id
library_tag_combo
lib_primer_index
library
tag_sequence
sample_name
experiment
dist_from_shore
transect
PCR_replicate
sample_type

move these fields to a "library_prep_worksheet"
vol_product_unclean_uL
vol_clean_input
vol_clean_elute_uL
vol_clean_total_uL
vol_removed_gel_uL
vol_removed_qubit_uL
vol_removed_BA
vol_remain_uL
conc_DNA_sample_ng/uL
mass_DNA_total_ng
fragment_size_BA
conc_DNA_BA_nguL
mass_DNA_pooled_ng
vol_pooled_undiluted_uL
notes

