#!/usr/bin/env bash

# written by Reiko at MBARI

# convert an OTU table, with taxonomic assignments, to an hdf5 biom file and the second bit of code to attach the biom metadata to the hdf5 file

biom convert \
    -i "${ANALYSIS_DIR}"/all_lib/OTUs_swarm/OTU_table.txt \
    -o "${ANALYSIS_DIR}"/all_lib/BOG18S_hdf5.biom \
    --table-type="OTU table" \
    --to-hdf5

biom add-metadata \
    -i "${ANALYSIS_DIR}"/all_lib/BOG18S_hdf5.biom \
    -o "${ANALYSIS_DIR}"/all_lib/BOG18S_hdf5_md.biom \
    --sample-metadata-fp /home/mbonteam/MBARI/reiko/scripts/BOG18S_biom_metadata_all.txt

# Adding the taxonomic assignments with this code.
biom add-metadata \
    -i "${ANALYSIS_DIR}"/all_lib/BOG18S_hdf5_md.biom \
    -o "${ANALYSIS_DIR}"/all_lib/BOG18S_hdf5_md_obs.biom \
    --observation-metadata-fp "${DIR}"/BOG18S_obs_md.txt \
    --int-fields size \
    --sc-separated taxonomy

# where the file BOG18S_obs_md.txt looks like:
#OTU_ID    taxonomy
# OTU_14311    Bacteria;Verrucomicrobia;Opitutae;Puniceicoccales;Puniceicoccaceae;Pelagicoccus;Pelagicoccus
# OTU_75047    Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Verrucomicrobiaceae;Rubritalea;Rubritalea
# OTU_101324    Bacteria;Cyanobacteria;Cyanobacteria;Cyanobacteria;Cyanobacteria;Cyanobacteria;Cyanobacteria
#
# or a Phinch version:
# #OTU_ID    taxonomy
# OTU_14311    k__Bacteria;p__Verrucomicrobia;c__Opitutae;o__Puniceicoccales;f__Puniceicoccaceae;g__Pelagicoccus;s__Pelagicoccus
# OTU_75047    k__Bacteria;p__Verrucomicrobia;c__Verrucomicrobiae;o__Verrucomicrobiales;f__Verrucomicrobiaceae;g__Rubritalea;s__Rubritalea
# OTU_101324    k__Bacteria;p__Cyanobacteria;c__Cyanobacteria;o__Cyanobacteria;f__Cyanobacteria;g__Cyanobacteria;s__Cyanobacteria
# OTU_101349    k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium
