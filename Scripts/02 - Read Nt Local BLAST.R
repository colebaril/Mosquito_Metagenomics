# LIBRARIES ----

require(pacman)
pacman::p_load(tidyverse, janitor, here, phylotools, assertr, 
               readxl, patchwork, reshape2, furrr)

# Set up `furrr` for `purrr` optimization 
future::plan(multisession)

# Load Contigs Dataset

contigs <- read_csv(here("Data/Contigs.csv"))


# READ BLAST ----

# N.B. Sample 42 did not have unique names and needed to process separately before
# adding to the main data frame.

sample_42_no_name <- future_map_dfr(list.files(pattern = ".xlsx",
                                               path = here("Raw/BLAST_2023/viruses/42"),
                                               full.names = TRUE),
                                    ~read_excel(.) %>% clean_names()) %>% 
  mutate(across(query, ~ gsub("bowtie2_merged", "42_bowtie2_merged", query)))


virus_master_2023 <- future_map_dfr(list.files(pattern = ".xlsx",
                                               path = here("Raw/BLAST_2023/Viruses"),
                                               full.names = TRUE),
                                    ~read_excel(.) %>% clean_names()) %>% 
  rbind(sample_42_no_name) %>% 
  mutate(result_desc = description_identity_percent,
         accession = accession_identity_percent) %>% 
  select(-contains("description"), -contains("accession_"), -contains("positive"), 
         -contains("bit"), -contains("number_of")) %>% 
  drop_na(lowest_e_value) %>% 
  filter(lowest_e_value <=1*10^-100) %>%
  mutate(across(query, ~ gsub("1_contig_", "", query))) %>% 
  mutate(across(query, ~ gsub("contig_", "", query))) %>% 
  mutate(across(query, ~ gsub("_bowtie2_", "-", query))) %>% 
  mutate(across(query, ~gsub("_bowtie2_merged_contig_", "-", query))) %>% 
  mutate(across(query, ~gsub("merged_", "", query))) %>% 
  left_join(contigs, by = c("query" = "id")) %>% rename("contig_sequence" = "seq.text") %>% 
  mutate(contig_length = nchar(contig_sequence)) %>% 
  
  relocate(contig_sequence, contig_length, coverage, .after = query) %>% 
  separate_wider_delim(query, delim = "-",
                       names = c("sample_number", "contig_number")) %>% 
  verify(n_distinct(sample_number) == 45) %>% 
  mutate(coverage = as.numeric(coverage)) %>% 
  filter(coverage >= 10) %>% 
  filter(contig_length >= 250) %>% 
  verify(contig_length >= 250) %>% 
  verify(coverage >= 10) %>% 
  verify(lowest_e_value <= 1*10^-100) %>% 
  mutate(novel_flag = case_when(greatest_identity_percent < 90 ~ "Presumptive Novel", TRUE ~ "Not Novel"))

# WRITE FASTA FOR tBLASTx

virus_master_2023_contigs <- virus_master_2023 %>% 
  select(sample_number, contig_number, contig_sequence) %>% 
  unite(col = "id", sample_number, contig_number, sep = "-") %>% 
  rename(seq.name = id,
         seq.text = contig_sequence)

dat2fasta(virus_master_2023_contigs, outfile = "Data/virus_master_2023_contigs.fasta")