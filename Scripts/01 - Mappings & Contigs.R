# LIBRARIES ----

require(pacman)
pacman::p_load(tidyverse, janitor, here, phylotools, assertr, 
               readxl, patchwork, reshape2, furrr)

# This script can take a considerable amount of time to run and require a lot of 
# memory depending on hardware.

# Set up `furrr` for `purrr` optimization 
future::plan(multisession)

# READ MAPPINGS ----

# N.B. Sample 42 did not have unique names and needed to process separately before
# adding to the main data frame.

mapping_42 <- future_map_dfr(list.files(pattern = ".xlsx",
                                        path = here("Raw/Mappings/42"),
                                        full.names = TRUE),
                             ~read_excel(.) %>% clean_names()) %>% 
  mutate(across(name, ~ gsub("bowtie2_merged", "42_bowtie2", name)))

mappings <- future_map_dfr(list.files(pattern = ".xlsx",
                                      path = here("Raw/Mappings"),
                                      full.names = TRUE),
                           ~read_excel(.) %>% clean_names()) %>% 
  rbind(mapping_42) %>% 
  select(1, 3, 6) %>% 
  mutate(name = gsub("_mapping", "", name)) %>% 
  mutate(name = gsub("_bowtie2_1_contig_", "-", name)) %>% 
  mutate(name = gsub("_bowtie2_contig_", "-", name)) %>% 
  rename("id" = "name",
         "coverage" = "average_coverage") %>% 
  separate_wider_delim(cols = id, delim = "-",
                       names = c("sample_number", "contig")) %>% 
  verify(n_distinct(sample_number) == 45) %>% 
  unite(col = "id", sample_number:contig, sep = "-")

# READ CONTIGS ----

# *****The following commented block should be run on the first time only.*****
# *****THIS BLOCK REQUIRES AT LEAST 16GB OF RAM*****

# Compile all contigs from 45 fasta files into a data frame


# contigs_42 <- future_map_dfr(list.files(pattern = "*assembly",
#                                         path = here("Raw/Contigs/42"),
#                                         full.names = TRUE),
#                              ~read.fasta(.)) %>%
#   mutate(across(seq.name, ~ gsub("bowtie2_", "42_bowtie2_", seq.name)))
# 
# 
# contigs_raw <- future_map_dfr(list.files(pattern = "*assembly",
#                                  path = here("Raw/Contigs"),
#                                  full.names = TRUE),
#                       ~read.fasta(.)) %>%
#   rbind(contigs_42)
# 
# contigs_raw %>%
#   write.csv(here("Data/contigs_raw.csv"))


# Attach mapping data to contig list

contigs_raw <- read_csv(here("Data/contigs_raw.csv")) %>% 
  select(-1)

contigs <- contigs_raw %>% 
  separate_wider_delim(seq.name, delim = ": ",
                       names = c("name", "coverage")) %>% 
  separate_wider_delim(name, delim = " ",
                       names = c("id", "junk"), too_many = "merge") %>% 
  select(-"junk", -"coverage") %>% 
  mutate(across(id, ~ gsub("_bowtie2_1_contig_", "-", id))) %>% 
  mutate(across(id, ~ gsub("_bowtie2_contig_", "-", id))) %>% 
  mutate(across(id, ~gsub("_bowtie2_merged_contig_", "-", id))) %>% 
  left_join(mappings, by = "id")

# SAVE DATA ----

mappings %>% 
  write.csv(here("Data/Mappings.csv"))

contigs %>% 
  write.csv(here("Data/Contigs.csv"))

