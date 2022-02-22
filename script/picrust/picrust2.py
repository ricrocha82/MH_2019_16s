#!/usr/bin/env python

# get the data (see picrust2_data_prep.R)

# run on command line
# <https://github.com/picrust/picrust2/wiki/Full-pipeline-script>
# https://github.com/picrust/picrust2/wiki/PICRUSt2-Tutorial-(v2.4.1)

# Run process
conda activate picrust2

# go to the directory where the files are and run the script

# The entire PICRUSt2 pipeline can be run using a single script, called picrust2_pipeline.py
# This script will run each of the 4 key steps:
# (1) sequence placement, 
# (2) hidden-state prediction of genomes, 
# (3) metagenome prediction, 
# (4) pathway-level predictions.
picrust2_pipeline.py -s seq_bac.fasta -i otutab_bac.txt \
-o picrust2_out -p 1

# Result annotation

add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                    -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
                    -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o pathways_out/path_abun_unstrat_descrip.tsv.gz
                    
# make a table with pathways -> MetaCyc website (https://metacyc.org - smart table)
                    
# STAMP (https://github.com/picrust/picrust2/wiki/STAMP-example)
conda deactivate
conda install -c bioconda stamp

stamp # or STAMP
