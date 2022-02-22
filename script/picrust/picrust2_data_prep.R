folder_path <- getwd()
library(tidyverse)

#--------------------
# Converting files
#--------------------
# code to convert csv to fasta
csv = read.csv('seq_arc.csv', header = T)
fa = character(2 * nrow(csv))
fa[c(TRUE, FALSE)] = sprintf("> %s", csv$id)
fa[c(FALSE, TRUE)] = csv$seq
writeLines(fa, "seq_arc.fasta")

csv = read.csv('seq_bac.csv', header = T)
fa = character(2 * nrow(csv))
fa[c(TRUE, FALSE)] = sprintf("> %s", csv$id)
fa[c(FALSE, TRUE)] = csv$seq
writeLines(fa, "seq_bac.fasta")

## or Biostrings package
library(Biostrings)
csv = read.csv('bac_2019.csv', header = T)
seq = csv$seq
names(seq) = csv$id
dna = DNAStringSet(seq)
writeXStringSet(dna, "seq_bac.fasta")

csv = read.csv('arc_2019.csv', header = T)
seq = csv$seq
names(seq) = csv$id
dna = DNAStringSet(seq)
writeXStringSet(dna, "seq_arc.fasta")

# otutab for picrsut2
# make on excel

# working the MetaCyc smart table to be used for picrust downstream analysis (see picrust2.py)
#tax_metacyc <- read.delim(file.choose(), stringsAsFactor = FALSE)
tax_metacyc <- read_delim(paste0(folder_path,"/data/picrust2_out/pathway_id_ontology.txt"), delim = "\t")
tax_metacyc <- tax_metacyc %>% rename('ontology' = 2)
tax_metacyc <- tax_metacyc %>% separate(ontology, into = letters[1:16], sep = " // ")
tax_metacyc <- tax_metacyc %>% select(-a, -b)
list.onto <- tax_metacyc %>% map(~ count(tibble(x = .), x))
# merge with picrust output (i.e. description)
otu_fun <- read_csv(paste0(folder_path,"/data/picrust2_out/path_abun_unstrat_descrip.csv")) %>% rename(pathway_id = pathway)
tax_metacyc <- otu_fun %>% select(pathway_id, description) %>% full_join(tax_metacyc)
# merge with Eric's table
eric_otu_fun <- read_csv(paste0(folder_path,"/data/picrust2_out/TAX_METACYC.csv")) %>% rename(pathway_id = pathway)
tax_metacyc <- tax_metacyc %>% full_join(eric_otu_fun)
write_csv(tax_metacyc, paste0(folder_path,"/data/picrust2_out/pathway_id_ontology.csv"))
# work on the table on excel
tax_metacyc <- read_csv(paste0(folder_path,"/data/picrust2_out/tax_metacyc.csv"))
tax_metacyc <- tax_metacyc %>% 
  mutate(across(c(pathway, ontology), ~str_replace_all(.," ", "-"))) %>%
  mutate(across(c(pathway, ontology), str_to_lower)) 
# check duplicates
t <- tax_metacyc %>% distinct(ontology) %>% arrange(ontology)
t <-tax_metacyc %>% distinct(pathway) %>% arrange(pathway)
# change names
tax_metacyc <- tax_metacyc %>%
  mutate(pathway = case_when(pathway == 'aromatic-compounds-biosynthesis' ~ 'aromatic-compound-biosynthesis',
                             pathway == 'carbohydrates-biosynthesis' ~ 'carbohydrate-biosynthesis', 
                             pathway == 'carbohydrates-degradation' ~ 'carbohydrate-degradation', 
                             pathway == 'sugar-and-acids-degradation' ~ 'sugar-acid-degradation',
                             pathway == 'sugars-and-polysaccharides-degradation' ~ 'carbohydrate-degradation',
                             TRUE ~ as.character(pathway)
  )) %>%
  mutate(ontology = case_when(ontology == 'adenosine-deoxy-denovo-bbiosynthesis' ~ 'adenosine-deoxy-de-novo-biosynthesis', 
                              ontology == 'aromatic-compounds-degradation' ~'aromatic-compound-degradation',
                              ontology == 'biosynthesis' ~ 'amine-polyamine-biosynthesis',
                              ontology == 'carbohydrates-biosynthesis' ~ 'carbohydrate-biosynthesis', 
                              ontology == 'cofactor,-carrier,-and-vitamin-biosynthesis' ~ 'adenosylcobamide-biosynthesis',
                              ontology == 'fatty-acid-biosynthesis' ~ 'fatty-acid-and-lipid-biosynthesis',
                              ontology == 'fatty-acid-degradation' ~ 'fatty-acid-and-lipid-degradation',
                              ontology == 'glycolysis' ~ 'glycolysis-variants',
                              ontology == 'polysaccharides-de' ~ 'polysaccharide-degradation',
                              ontology == 'reductive-tca-cycle-ii_autotrophic-co2-fixation' ~ 'autotrophic-co2-fixation_reductive-tca-cycles',
                              ontology == 'siderophores-biosynthesis' ~ 'siderophore-and-metallophore-biosynthesis',
                              ontology == 'sugar-derivatives' ~ 'sugar-derivatives-degradation',
                              ontology == 'super-pathways' ~ 'super-pathway',
                              ontology == 'superpathways' ~ 'super-pathway',
                              ontology == 'unsaturated-fatty-acids-biosynthesis' ~ 'unsaturated-fatty-acid-biosynthesis',
                              TRUE ~ as.character(ontology)
  ))

write_csv(tax_metacyc, paste0(folder_path,"/data/picrust2_out/tax_metacyc.csv"))
