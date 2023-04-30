# This script takes in viral reads hits as well as all contig sequences. Contig
# sequences are joined together, resulting in contig sequences only for contigs
# that had viral hits. This enables me to save time by only re-blasting hits. 
# I needed to re-BLAST hits because new more recent relatively were deposited to
# GenBank which could impact BLAST results.

# A link to the contig sequences saved in a file called Raw/Contigs in my working directory
# can be found in the README.md file as it is too large for GitHub.

require(pacman)
p_load(tidyverse, here, phylotools, janitor, readxl, openxlsx)

# VIRUS HITS DF ----

virus_df <- read_csv(here("Raw/Mosquito_Meta - Sheet1.csv")) %>% 
  clean_names() %>% 
  separate_wider_delim(query, delim = "_bowtie2_", 
                       names = c("sample_number", "contig")) %>% 
  mutate(clc_desc = description_hsp_length) %>% 
  separate(description_hsp_length, sep = "(?<=.\\d)\\s(?=\\w)",
           into = c("accession", "name"), extra = "merge") %>% 
  mutate(across(contig, ~ gsub("1_contig_", "", contig))) %>% 
  mutate(across(contig, ~ gsub("contig_", "", contig)))

virus_df_join <- virus_df %>% 
  select(sample_number, contig, ) %>% 
  unite(col = "id", sample_number, contig, sep = "-")



# CONTIGS ----

contigs_raw <- map_df(list.files(pattern = "*assembly",
                                 path = here("Raw/Contigs"),
                                 full.names = TRUE),
                      ~read.fasta(.)
)

contigs <- contigs_raw %>% 
  separate_wider_delim(seq.name, delim = ": ",
                       names = c("name", "coverage")) %>% 
  separate_wider_delim(name, delim = " ",
                       names = c("id", "junk"), too_many = "merge") %>% 
  select(-"junk", -"coverage") %>% 
  mutate(across(id, ~ gsub("_bowtie2_1_contig_", "-", id))) %>% 
  mutate(across(id, ~ gsub("_bowtie2_contig_", "-", id)))


# JOIN CONTIGS ----

contig_joined <- virus_df_join %>% 
  left_join(contigs, by = "id") %>% 
  rename(seq.name = id)

# CREATE FASTA ---- 

dat2fasta(contig_joined, outfile = "mosquito_contigs.fasta")

# 2023 BLAST UPLOAD ----

blastn_2023 <- map_df(list.files(pattern = ".xlsx",
                                 path = here("Raw/Blastn_2023"),
                                 full.names = TRUE),
                      ~read_excel(.) %>% clean_names() %>% 
                        select(query, greatest_identity_percent, 
                               accession_hsp_length, description_hsp_length) %>% 
                      rename(id = query,
                             greatest_identity_percent_2023 = greatest_identity_percent,
                             accession_hsp_length_2023 = accession_hsp_length,
                             name_2023 = description_hsp_length,
                             accession_2023 = accession_hsp_length) %>% 
                        mutate(across(name_2023, ~ gsub("MAG: ", "", name_2023)))
                      )

# 2020 BLAST ----

virus_df_2020_blast <- virus_df %>% 
  select(sample_number, contig, name, accession, greatest_identity_percent) %>% 
  unite(col = "id", sample_number, contig, sep = "-") %>% 
  rename(name_2020 = name,
         accession_2020 = accession,
         greatest_identity_percent_2020 = greatest_identity_percent) %>% 
  mutate(across(accession_2020, ~ gsub("\\.1", "", accession_2020)))

# COMPARE 2020/2023 BLAST RESULT ----

blast_2020_2023_comparison <- virus_df_2020_blast %>% 
  left_join(blastn_2023, by = "id") %>% 
  mutate(name_flag = case_when(name_2020 != name_2023 ~ "Different", TRUE ~ "Same")) %>% 
  mutate(accession_flag = case_when(accession_2020 != accession_2023 ~ "Different", TRUE ~ "Same")) %>% 
  select(-greatest_identity_percent_2020, -greatest_identity_percent_2023) %>% 
  relocate(name_2023, .before = accession_2020) %>% 
  left_join(contigs, by = "id") %>% 
  rename("contig_sequence" = "seq.text") %>% 
  separate_wider_delim(id, delim = "-",
                       names = c("sample_number", "contig"))

# In many cases the name changed but accession remained the same. Every time the 
# accession changed, it's a different virus. 

differences_accession <- blast_2020_2023_comparison %>% 
  group_by(accession_flag) %>% 
  summarise(n = n())

wb <- createWorkbook(creator = "Cole Baril")

addWorksheet(wb, "Summary")
addWorksheet(wb, "Count")

writeData(wb, sheet = "Summary", blast_2020_2023_comparison)
writeData(wb, sheet = "Count", differences_accession)

saveWorkbook(wb, file = "2020-2023_Blast_Differences.xlsx", overwrite = TRUE)




