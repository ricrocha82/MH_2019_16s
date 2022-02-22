
library(dada2)
library(tidyverse)
packageVersion(c("dada2","tidyverse"))
# https://benjjneb.github.io/dada2/index.html
# Assign taxonomy - DADA2
# FASTA to tabular = https://usegalaxy.org/

folder_path <- "~/R/git_hub/MH_2019/data/"
#folder_path <- '~/Documents/R/phd/2019/data/'
#setwd(folder_path)
seqtab.nochim.bac <- readRDS(paste0(folder_path,"sequences_fasta_16S/dada2/seqtab.nochim.RDS"))


doParallel::registerDoParallel()

# Taxonomic reference data @ <https://benjjneb.github.io/dada2/training.html>
# 16S ----
refFasta.16s <- '~/R/git_hub/MH_2019/silva_nr99_v138.1_train_set.fa.gz' # updated Mar 10 ,2021


# 2019 data ----


# BACTERIA -----
taxa_ASV_bac_MH2019 <- assignTaxonomy(seqtab.nochim.bac, refFasta.16s, minBoot = 50, tryRC = FALSE,
                                      outputBootstraps = FALSE, 
                                      taxLevels = c("Kingdom", "Phylum", "Class",
                                                    "Order", "Family", "Genus", "Species"), 
                                      multithread = T, verbose = FALSE)

write.csv(taxa_ASV_bac_MH2019, '/data/tax_bac_2019.csv')


# making data frames ----
# table with samples id and extraction id
extract_id <- read_csv(paste0(folder_path,"/data/extraction_ID.csv")) 
extract_id <- extract_id %>% arrange(ID_extraction)
extract_id <- extract_id %>% 
  mutate(sample_extr = str_c("RRPS-MH-", ID_extraction)) %>%
  select(sample_extr, Samples)
sample.names <- extract_id %>% pull(Samples)
sample.extr <- extract_id %>% pull(sample_extr)

# to upload sequences information to NCBI
# meta_date_collection <- read_csv(paste0(folder_path,"/samples_info.csv")) %>% select(sample_id, date_collected)
# meta <- read_csv(paste0(folder_path,"/meta.csv")) 
# meta_extr <- extract_id %>% 
#                 rename(sample_id = Samples) %>% 
#                 left_join(meta) %>% 
#                 left_join(meta_date_collection) %>%
#                 mutate(lat_long = str_c(lat, long, sep = "_"))
# meta_extr %>% write_csv(paste0(folder_path,"/data/extraction_meta_seq_ID.csv")) 

 
## Bac ----

asv_bac <- t(seqtab.nochim.bac) %>% as.data.frame() %>% 
  rownames_to_column("seq")  %>% as_tibble() %>%
  rename_with(~ sample.names[which(sample.extr == .x)], .cols = sample.extr) # 25,505 ASVs

asv_bac %>% 
   mutate(total = rowSums(across(where(is.numeric)))) %>%
 # rowwise() %>% mutate(total = sum(c_across(where(is.numeric)))) %>% 
  filter(total ==1) # 4,209 ASVs
asv_bac %>%  mutate(total = reduce(select(., where(is.numeric)), `+`)) %>%
  filter(total ==2) # 3,161 ASVs

write_csv(asv_bac, paste0(folder_path,'ASV_bac_2019.csv'))
tax_bac <- read_csv(paste0(getwd(),'/data/tax_bac_2019.csv')) %>% rename(seq = 1)

tax_bac %>% nrow()
tax_bac %>% group_by(Kingdom) %>% count(Phylum) %>% View()
asv_bac <- asv_bac %>% left_join(tax_bac)
asv_bac %>% nrow() # 25505
asv_bac %>% group_by(Kingdom) %>% count(Phylum) %>% View() # checking
# filter non-bacteria
asv_bac <- asv_bac %>% filter(Kingdom == "Bacteria")
asv_bac %>% count(Kingdom)
asv_bac %>% count(Phylum, sort = T)
# check NAs at the Class level
asv_bac %>% filter(is.na(Class)) %>% count(Phylum) %>% View
# remove Phylum NAs
asv_bac <- asv_bac %>% drop_na(Phylum)
asv_bac <- asv_bac %>%
  mutate(id = str_c('b_asv_',1:nrow(.))) %>%
  relocate(id) 

write_csv(asv_bac, paste0(folder_path,'bac_2019.csv')) 






