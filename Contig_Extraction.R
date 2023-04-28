# This sceipt takes in viral reads hits as well as all contig sequences. Contig
# sequences are joined together, resulting in contig sequences only for contigs
# that had viral hits. This enables me to save time by only re-blasting hits. 
# I needed to re-BLAST hits because new more recent relatively were deposited to
# GenBank which could impact BLAST results.

# A link to the contig sequences saved in a file called Raw/Contigs in my working directory
# can be found in the README.md file as it is too large for GitHub.

require(pacman)
p_load(tidyverse, here, phylotools, janitor)

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
  select(sample_number, contig) %>% 
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



