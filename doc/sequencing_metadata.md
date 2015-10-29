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
- library_index (primary_index_seq)
- library_name (primary_index_name)
- primer1_name
- primer1_sequence_core
- primer1_index_name (secondary_index_1_name)
- primer1_index_sequence (secondary_index_1_seq)
- primer1_sequence_full
- primer2_name
- primer2_sequence_core
- primer2_index_name (secondary_index_2_name)
- primer2_index_sequence (secondary_index_2_name)
- primer2_sequence_full
- interprimer_size_mean (in bp; min/max?)
- frag_size_sequenced (size in bp of fragments when they went on the sequencer)

- env_sample_name
- DNA_extraction_name
- PCR_name

- sample_type : ( control_pos | control_neg | environmental )

- project : useful for combining multiple projects per run (e.g. grid, spatial)
- locus
- order

Don't include environmental data (e.g. temp, salinity, lat, lon). Put this in another file and cross-reference using sequencing sample_id.

Don't include lab data (e.g. mass_DNA_total_ng, vol_pooled_undiluted_uL, conc_DNA_sample_ng/uL). Put these in another file and cross-reference using sequencing sample_id.

Don't include a "notes" column, especially before the last column. This is likely to contain spaces, and if those come before another required field, it could cause some terribleness.
