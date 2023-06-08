require(pacman)
pacman::p_load(tidyverse, janitor, here, DT, gt, phylotools, assertr, readxl, patchwork, reshape2, ggtext, openxlsx)

# Compile all contigs from 45 fasta files into a data frame
# ====== WARNING: THIS BLOCK REQUIRES AT LEAST 16GB OF RAM 
#contigs_raw <- map_df(list.files(pattern = "*assembly",
#                                 path = here("Raw/Contigs"),
#                                 full.names = TRUE),
#                      ~read.fasta(.)
#)
#contigs_raw %>% 
#  write.csv(here("Data/contigs_raw.csv"))

contigs_raw <- read_csv(here("Data/contigs_raw.csv")) %>% 
  select(-1)

contigs <- contigs_raw %>% 
  separate_wider_delim(seq.name, delim = ": ",
                       names = c("name", "coverage")) %>% 
  separate_wider_delim(name, delim = " ",
                       names = c("id", "junk"), too_many = "merge") %>% 
  select(-"junk") %>% 
  mutate(across(id, ~ gsub("_bowtie2_1_contig_", "-", id))) %>% 
  mutate(across(id, ~ gsub("_bowtie2_contig_", "-", id)))

virus_master_2023 <- map_df(list.files(pattern = ".xlsx",
                                       path = here("Raw/BLAST_2023/Viruses"),
                                       full.names = TRUE),
                            ~read_excel(.) %>% clean_names()) %>% 
  mutate(result_desc = description_identity_percent,
         accession = accession_identity_percent) %>% 
  select(-contains("description"), -contains("accession_"), -contains("positive"), 
         -contains("bit"), -contains("number_of")) %>% 
  drop_na(lowest_e_value) %>% 
  filter(lowest_e_value <=1*10^-100) %>%
  mutate(across(query, ~ gsub("1_contig_", "", query))) %>% 
  mutate(across(query, ~ gsub("contig_", "", query))) %>% 
  mutate(across(query, ~ gsub("_bowtie2_", "-", query))) %>% 
  left_join(contigs, by = c("query" = "id")) %>% rename("contig_sequence" = "seq.text") %>% 
  mutate(contig_length = nchar(contig_sequence)) %>% 
  #left_join(virus_df_coverage_stats, by = c("query" = "id")) %>% 
  # rename("read_count" = "count",
  #"contig_length" = "length",
  # "coverage_depth" = "coverage") %>% 
  relocate(contig_sequence, contig_length, coverage, .after = query) %>% 
  separate_wider_delim(query, delim = "-",
                       names = c("sample_number", "contig_number")) %>% 
  mutate(coverage = as.numeric(coverage)) %>% 
  filter(coverage >= 10) %>% 
  filter(contig_length >= 250) %>% 
  verify(contig_length >= 250) %>% 
  verify(coverage >= 10) %>% 
  verify(lowest_e_value <= 1*10^-100) %>% 
  mutate(novel_flag = case_when(greatest_identity_percent < 90 ~ "Presumptive Novel", TRUE ~ "Not Novel"))

virus_master_2023 <- virus_master_2023 %>% 
  # COLLECTION YEAR
  mutate(collection_year = case_when(sample_number %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                          "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                                                          "20", "24", "36", "32") ~ "2020", 
                                     TRUE ~ "2021")) %>% relocate(collection_year, .after = contig_sequence) %>% 
  # COLLECTION LOCATION
  mutate(location_pool = case_when(sample_number %in% c("1", "10", "13", "19", "20", "29", "30", "31") ~ "Brandon",
                                   sample_number %in% c("2", "14", "21", "22", "23", "24", "40", "45", "36") ~ "Shoal Lake",
                                   sample_number %in% c("3", "15") ~ "Virden & Souris",
                                   sample_number %in% c("17", "18", "8") ~ "WPG Insect Contol",
                                   sample_number %in% c("7", "37", "38", "39", "44") ~ "Cypress River",
                                   sample_number %in% c("12", "6") ~ "Cypress River & Carberry",
                                   sample_number %in% c("5", "35") ~ "Killarney",
                                   sample_number %in% c("11") ~ "Shoal Lake, Virden & Souris",
                                   sample_number %in% c("25", "26") ~ "Souris",
                                   sample_number %in% c("27", "28") ~ "Virden",
                                   sample_number %in% c("32") ~ "Multiple",
                                   sample_number %in% c("42", "43") ~ "Carberry",
                                   sample_number %in% c("41") ~ "Newdale",
                                   sample_number %in% c("4", "34", "33") ~ "Boissevain",
                                   sample_number %in% c("9") ~ "Western MB",
                                   sample_number %in% c("16") ~ "Boissevain & Killarney",
  )) %>% relocate(location_pool, .after = collection_year) %>% 
  # MOSQUITO SPECIES
  mutate(mosquito_species = case_when(sample_number %in% c("1", "2", "3", "4", "5", "6", "8", 
                                                           "20", "22", "25", "27", "29", 
                                                           "33", "35", "36", "37", "41", "42", "45") ~ "Aedes vexans",
                                      sample_number %in% c("31") ~ "Aedes canadensis",
                                      sample_number %in% c("10", "11", "12", "17", "19", "23", "26", 
                                                           "30", "34", "38", "43" ) ~ "Culex tarsalis",
                                      sample_number %in% c("13", "14", "15", "16", "40") ~ "Ochlerotatus dorsalis",
                                      sample_number %in% c("24") ~ "Ochlerotatus flavescens",
                                      sample_number %in% c("32") ~ "Ochlerotatus triseriatus",
                                      sample_number %in% c("21") ~ "Anopheles earlei",
                                      sample_number %in% c("7", "9", "18", "28", "39", "44") ~ "Coquillettidia perturbans"
  )) %>% relocate(mosquito_species, .after = location_pool) %>% 
  # SIMPLIFIED VIRUS NAME
  mutate(virus_name = case_when(str_detect(result_desc, regex("Soybean thrips iflavirus 4", ignore_case = TRUE)) ~ "Soybean thrips iflavirus 4",
                                str_detect(result_desc, regex("Soybean thrips dicistrovirus", ignore_case = TRUE)) ~ "Soybean thrips dicistrovirus",
                                str_detect(result_desc, regex("Black queen cell virus isolate", ignore_case = TRUE)) ~ "Black queen cell virus",
                                str_detect(result_desc, regex("Hubei noda-like virus 12", ignore_case = TRUE)) ~ "Hubei noda-like virus 12",
                                str_detect(result_desc, regex("Mekrijarvi Negevirus", ignore_case = TRUE)) ~ "Mekrijarvi Negevirus",
                                str_detect(result_desc, regex("Cordoba virus", ignore_case = TRUE)) ~ "Cordoba virus",
                                str_detect(result_desc, regex("Hubei macula-like virus 3", ignore_case = TRUE)) ~ "Hubei macula-like virus 3",
                                str_detect(result_desc, regex("Culex densovirus", ignore_case = TRUE)) ~ "Culex densovirus",
                                str_detect(result_desc, regex("Gouley virus", ignore_case = TRUE)) ~ "Gouley virus",
                                str_detect(result_desc, regex("Hubei virga-like virus 2", ignore_case = TRUE)) ~ "Hubei virga-like virus 2",
                                str_detect(result_desc, regex("Placeda virus", ignore_case = TRUE)) ~ "Placeda virus",
                                str_detect(result_desc, regex("Culex narnavirus 1", ignore_case = TRUE)) ~ "Culex narnavirus 1",
                                str_detect(result_desc, regex("Marma virus", ignore_case = TRUE)) ~ "Marma virus",
                                str_detect(result_desc, regex("Hubei mosquito virus 4", ignore_case = TRUE)) ~ "Hubei mosquito virus 4",
                                str_detect(result_desc, regex("Culex iflavilike virus 4", ignore_case = TRUE)) ~ "Culex iflavilike virus 4",
                                str_detect(result_desc, regex("Culex Iflavi-like virus 4", ignore_case = TRUE)) ~ "Culex Iflavi-like virus 4",
                                str_detect(result_desc, regex("Culex iflavilike virus 3", ignore_case = TRUE)) ~ "Culex iflavilike virus 3",
                                str_detect(result_desc, regex("Culex iflavi-like virus 3", ignore_case = TRUE)) ~ "Culex iflavilike virus 3",
                                str_detect(result_desc, regex("Partiti-like culex mosquito virus", ignore_case = TRUE)) ~ 
                                  "Partiti-like culex mosquito virus",
                                str_detect(result_desc, regex("Elisy virus", ignore_case = TRUE)) ~ "Elisy virus",
                                str_detect(result_desc, regex("Canya virus", ignore_case = TRUE)) ~ "Canya virus",
                                str_detect(result_desc, regex("Manitoba virus", ignore_case = TRUE)) ~ "Manitoba virus",
                                str_detect(result_desc, regex("Flanders hapavirus", ignore_case = TRUE))  ~ "Flanders hapavirus",
                                str_detect(result_desc, regex("Hapavirus flanders", ignore_case = TRUE)) ~ "Flanders hapavirus",
                                str_detect(result_desc, regex("Culex rhabdovirus", ignore_case = TRUE)) ~ "Culex rhabdovirus",
                                str_detect(result_desc, regex("Merida virus", ignore_case = TRUE)) ~ "Merida virus",
                                str_detect(result_desc, regex("Merida-like virus", ignore_case = TRUE)) ~ "Merida virus",
                                str_detect(result_desc, regex("Flanders hapavirus", ignore_case = TRUE)) ~ "Flanders hapavirus",
                                str_detect(result_desc, regex("Wuhan mosquito virus 6", ignore_case = TRUE)) ~ "Wuhan mosquito virus 6",
                                str_detect(result_desc, regex("Chuvirus Mos8Chu0", ignore_case = TRUE)) ~ "Chuvirus Mos8Chu0",
                                str_detect(result_desc, regex("Flanders hapavirus", ignore_case = TRUE)) ~ "Flanders hapavirus",
                                str_detect(result_desc, regex("Flanders virus", ignore_case = TRUE)) ~ "Flanders hapavirus",
                                str_detect(result_desc, regex("Culex bunyavirus 2", ignore_case = TRUE)) ~ "Culex bunyavirus 2",
                                str_detect(result_desc, regex("Riverside virus 1", ignore_case = TRUE)) ~ "Manitoba Rhabdovirus 1",
                                str_detect(result_desc, regex("Culex rhabdo-like virus", ignore_case = TRUE)) & 
                                  novel_flag == "Presumptive Novel" ~ "Manitoba Rhabdovirus 1",
                                str_detect(result_desc, regex("Yongsan picorna-like virus 1", ignore_case = TRUE)) ~ 
                                  "Manitoba picorna-like virus 1",
                                str_detect(result_desc, regex("Thrace picorna-like virus 1", ignore_case = TRUE)) ~ 
                                  "Manitoba picorna-like virus 1",
                                str_detect(result_desc, regex("Nor picorna-like virus", ignore_case = TRUE)) ~ 
                                  "Manitoba picorna-like virus 1",
                                str_detect(result_desc, regex("Bro virus", ignore_case = TRUE)) ~ 
                                  "Manitoba mononega-like virus 1",
                                str_detect(result_desc, regex("Joensuu anphevirus", ignore_case = TRUE)) ~ 
                                  "Manitoba mononega-like virus 1",
                                str_detect(result_desc, regex("Budalangi Iflavi-like virus", ignore_case = TRUE)) ~ 
                                  "Manitoba iflavi-like virus 1",
                                str_detect(result_desc, regex("Hypsignathus monstrosus tombus-like virus", ignore_case = TRUE)) ~ 
                                  "Manitoba tombus-like virus 1",
                                str_detect(result_desc, regex("Tiger mosquito bi-segmented tombus-like virus", ignore_case = TRUE)) ~ 
                                  "Manitoba tombus-like virus 1",
                                str_detect(result_desc, regex("Culicine-associated Z virus", ignore_case = TRUE)) ~ "Ballard Lake virus",
                                str_detect(result_desc, regex("Port Bolivar virus ", ignore_case = TRUE)) ~ "Ballard Lake virus",
                                str_detect(result_desc, regex("Espirito Santo virus", ignore_case = TRUE)) ~ "Ballard Lake virus",
                                str_detect(result_desc, regex("Culicine-associated Z virus", ignore_case = TRUE)) ~ "Ballard Lake virus",
                                str_detect(result_desc, regex("Partitivirus-like Culex mosquito virus", ignore_case = TRUE)) ~ 
                                  "Partitivirus-like Culex mosquito virus",
                                str_detect(result_desc, regex("Snelk virus", ignore_case = TRUE)) ~ "Snelk virus",
                                str_detect(result_desc, regex("Cafluga virus", ignore_case = TRUE)) ~ "Cafluga virus",
                                str_detect(result_desc, regex("chuvirus", ignore_case = TRUE)) ~ "Chuvirus",
                                str_detect(result_desc, regex("Ballard Lake virus", ignore_case = TRUE)) ~ "Ballard Lake virus",
                                str_detect(result_desc, regex("Hanko iflavirus 1", ignore_case = TRUE)) ~ "Hanko iflavirus 1",
                                str_detect(result_desc, regex("Astopletus virus", ignore_case = TRUE)) ~ "Astopletus virus",
                                str_detect(result_desc, regex("Des Moines River virus", ignore_case = TRUE)) ~ "Des Moines River virus",
                                str_detect(result_desc, regex("Big Cypress virus", ignore_case = TRUE)) ~ "Manitoba mononega-like virus 2",
                                str_detect(result_desc, regex("Atrato picorna-like virus", ignore_case = TRUE)) ~ "Manitoba picorna-like virus 2",
                                str_detect(result_desc, regex("Pedersore iflavirus", ignore_case = TRUE)) ~ "Manitoba iflavirus 1",
                                str_detect(result_desc, regex("Hubei mosquito virus 1", ignore_case = TRUE)) ~ "Manitoba toti-like virus 1",
                                str_detect(result_desc, regex("Dicistroviridae sp. isolate 78", ignore_case = TRUE)) ~ "Manitoba picorna-like virus 3",
                                str_detect(result_desc, regex("Culex originated Tymoviridae", ignore_case = TRUE)) ~ "Manitoba tymo-like virus 1",
                                str_detect(result_desc, regex("Black queen", ignore_case = TRUE)) ~ "Black Queen Cell Virus",
                                str_detect(result_desc, regex("Culex rhabdo-like", ignore_case = TRUE)) ~ "Culex Rhabdo-like virus"
                                
  )) %>% 
  # VIRAL GENOME
  mutate(genome = case_when(str_detect(virus_name, regex("Soybean thrips iflavirus 4", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Soybean thrips dicistrovirus", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Black queen cell virus isolate", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Hubei noda-like virus 12", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Mekrijarvi Negevirus", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Cordoba virus", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Hubei macula-like virus 3", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Culex densovirus", ignore_case = TRUE)) ~ "ssDNA",
                            str_detect(virus_name, regex("Gouley virus", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Hubei virga-like virus 2", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Placeda virus", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Culex narnavirus 1", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Marma virus", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Hubei mosquito virus 4", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Culex iflavilike virus 4", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Culex Iflavi-like virus 4", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Culex iflavilike virus 3", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Culex iflavi-like virus 3", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Partiti-like culex mosquito virus", ignore_case = TRUE)) ~ "dsRNA",
                            str_detect(virus_name, regex("Elisy virus", ignore_case = TRUE)) ~ "-ssRNA",
                            str_detect(virus_name, regex("Canya virus", ignore_case = TRUE)) ~ "-ssRNA",
                            str_detect(virus_name, regex("Manitoba virus", ignore_case = TRUE)) ~ "-ssRNA",
                            str_detect(virus_name, regex("Flanders hapavirus", ignore_case = TRUE))  ~ "-ssRNA",
                            str_detect(virus_name, regex("Hapavirus flanders", ignore_case = TRUE)) ~ "-ssRNA",
                            str_detect(virus_name, regex("Culex rhabdovirus", ignore_case = TRUE)) ~ "-ssRNA",
                            str_detect(virus_name, regex("Merida virus", ignore_case = TRUE)) ~ "-ssRNA",
                            str_detect(virus_name, regex("Merida-like virus", ignore_case = TRUE)) ~ "-ssRNA",
                            str_detect(virus_name, regex("Flanders hapavirus", ignore_case = TRUE)) ~ "-ssRNA",
                            str_detect(virus_name, regex("Wuhan mosquito virus 6", ignore_case = TRUE)) ~ "-ssRNA",
                            str_detect(virus_name, regex("Chuvirus Mos8Chu0", ignore_case = TRUE)) ~ "-ssRNA",
                            str_detect(virus_name, regex("Flanders hapavirus", ignore_case = TRUE)) ~ "-ssRNA",
                            str_detect(virus_name, regex("Flanders virus", ignore_case = TRUE)) ~ "-ssRNA",
                            str_detect(virus_name, regex("Culex bunyavirus 2", ignore_case = TRUE)) ~ "-ssRNA",
                            str_detect(virus_name, regex("Riverside virus 1", ignore_case = TRUE)) ~ "-ssRNA",
                            str_detect(virus_name, regex("Culex rhabdo-like virus", ignore_case = TRUE)) & 
                              novel_flag == "Presumptive Novel" ~ "-ssRNA",
                            str_detect(virus_name, regex("Yongsan picorna-like virus 1", ignore_case = TRUE)) ~ 
                              "+ssRNA",
                            str_detect(virus_name, regex("Thrace picorna-like virus 1", ignore_case = TRUE)) ~ 
                              "+ssRNA",
                            str_detect(virus_name, regex("Nor picorna-like virus", ignore_case = TRUE)) ~ 
                              "+ssRNA",
                            str_detect(virus_name, regex("Manitoba mononega-like virus 2", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Manitoba mononega-like virus 1", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Budalangi Iflavi-like virus", ignore_case = TRUE)) ~ 
                              "+ssRNA",
                            str_detect(virus_name, regex("Hypsignathus monstrosus tombus-like virus", ignore_case = TRUE)) ~ 
                              "+ssRNA",
                            str_detect(virus_name, regex("Tiger mosquito bi-segmented tombus-like virus", ignore_case = TRUE)) ~ 
                              "+ssRNA",
                            str_detect(virus_name, regex("Culicine-associated Z virus", ignore_case = TRUE)) ~ "dsRNA",
                            str_detect(virus_name, regex("Port Bolivar virus ", ignore_case = TRUE)) ~ "dsRNA",
                            str_detect(virus_name, regex("Espirito Santo virus", ignore_case = TRUE)) ~ "dsRNA",
                            str_detect(virus_name, regex("Culicine-associated Z virus", ignore_case = TRUE)) ~ "dsRNA",
                            str_detect(virus_name, regex("Partitivirus-like Culex mosquito virus", ignore_case = TRUE)) ~ 
                              "dsRNA",
                            str_detect(virus_name, regex("Snelk virus", ignore_case = TRUE)) ~ "dsRNA",
                            str_detect(virus_name, regex("Cafluga virus", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("chuvirus", ignore_case = TRUE)) ~ "-ssRNA",
                            str_detect(virus_name, regex("Ballard Lake virus", ignore_case = TRUE)) ~ "dsRNA",
                            str_detect(virus_name, regex("Hanko iflavirus 1", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Astopletus virus", ignore_case = TRUE)) ~ "-ssRNA",
                            str_detect(virus_name, regex("Des Moines River virus", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Manitoba mononega-like virus 2", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Manitoba picorna-like virus 2", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Manitoba iflavi-like virus 1", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Manitoba toti-like virus 1", ignore_case = TRUE)) ~ "dsRNA",
                            str_detect(virus_name, regex("Manitoba picorna-like virus 3", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Manitoba tymo-like virus 1", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Black Queen Cell Virus", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Manitoba tombus-like virus 1", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Manitoba picorna-like virus 1", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Manitoba rhabdovirus 1", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Manitoba iflavirus 1", ignore_case = TRUE)) ~ "+ssRNA",
                            str_detect(virus_name, regex("Culex Rhabdo-like virus", ignore_case = TRUE)) ~ "-ssRNA",
                            
                            
                            
  )) %>% 
  # VIRAL FAMILY
  mutate(viral_family = case_when(str_detect(virus_name, regex("Soybean thrips iflavirus 4", ignore_case = TRUE)) ~ "Iflaviridae",
                            str_detect(virus_name, regex("Soybean thrips dicistrovirus", ignore_case = TRUE)) ~ "Dicistroviridae",
                            str_detect(virus_name, regex("Black queen cell virus isolate", ignore_case = TRUE)) ~ "Dicistroviridae",
                            str_detect(virus_name, regex("Hubei noda-like virus 12", ignore_case = TRUE)) ~ "Nodaviridae",
                            str_detect(virus_name, regex("Mekrijarvi Negevirus", ignore_case = TRUE)) ~ "Negevirus",
                            str_detect(virus_name, regex("Cordoba virus", ignore_case = TRUE)) ~ "Negevirus",
                            str_detect(virus_name, regex("Hubei macula-like virus 3", ignore_case = TRUE)) ~ "Tymoviridae",
                            str_detect(virus_name, regex("Culex densovirus", ignore_case = TRUE)) ~ "Parvoviridae",
                            str_detect(virus_name, regex("Gouley virus", ignore_case = TRUE)) ~ "Totiviridae",
                            str_detect(virus_name, regex("Hubei virga-like virus 2", ignore_case = TRUE)) ~ "Virgaviridae",
                            str_detect(virus_name, regex("Placeda virus", ignore_case = TRUE)) ~ "Flaviviridae",
                            str_detect(virus_name, regex("Culex narnavirus 1", ignore_case = TRUE)) ~ "Narnaviridae",
                            str_detect(virus_name, regex("Marma virus", ignore_case = TRUE)) ~ "Luteoviridae",
                            str_detect(virus_name, regex("Hubei mosquito virus 4", ignore_case = TRUE)) ~ "Tombusviridae",
                            str_detect(virus_name, regex("Culex iflavilike virus 4", ignore_case = TRUE)) ~ "Iflaviridae",
                            str_detect(virus_name, regex("Culex Iflavi-like virus 4", ignore_case = TRUE)) ~ "Iflaviridae",
                            str_detect(virus_name, regex("Culex iflavilike virus 3", ignore_case = TRUE)) ~ "Iflaviridae",
                            str_detect(virus_name, regex("Culex iflavi-like virus 3", ignore_case = TRUE)) ~ "Iflaviridae",
                            str_detect(virus_name, regex("Partiti-like culex mosquito virus", ignore_case = TRUE)) ~ "Partitiviridae",
                            str_detect(virus_name, regex("Elisy virus", ignore_case = TRUE)) ~ "Rhabdoviridae",
                            str_detect(virus_name, regex("Canya virus", ignore_case = TRUE)) ~ "Rhabdoviridae",
                            str_detect(virus_name, regex("Manitoba virus", ignore_case = TRUE)) ~ "Rhabdoviridae",
                            str_detect(virus_name, regex("Flanders hapavirus", ignore_case = TRUE))  ~ "Flanders hapavirus",
                            str_detect(virus_name, regex("Hapavirus flanders", ignore_case = TRUE)) ~ "Flanders hapavirus",
                            str_detect(virus_name, regex("Culex rhabdovirus", ignore_case = TRUE)) ~ "Culex rhabdovirus",
                            str_detect(virus_name, regex("Merida virus", ignore_case = TRUE)) ~ "Rhabdoviridae",
                            str_detect(virus_name, regex("Merida-like virus", ignore_case = TRUE)) ~ "Rhabdoviridae",
                            str_detect(virus_name, regex("Flanders hapavirus", ignore_case = TRUE)) ~ "Rhabdoviridae",
                            str_detect(virus_name, regex("Wuhan mosquito virus 6", ignore_case = TRUE)) ~ "Orthomyxoviridae",
                            str_detect(virus_name, regex("Chuvirus Mos8Chu0", ignore_case = TRUE)) ~ "Chuviridae",
                            str_detect(virus_name, regex("Flanders hapavirus", ignore_case = TRUE)) ~ "Rhabdoviridae",
                            str_detect(virus_name, regex("Flanders virus", ignore_case = TRUE)) ~ "Rhabdoviridae",
                            str_detect(virus_name, regex("Culex bunyavirus 2", ignore_case = TRUE)) ~ "Peribunyaviridae",
                            str_detect(virus_name, regex("Riverside virus 1", ignore_case = TRUE)) ~ "Rhabdoviridae",
                            str_detect(virus_name, regex("Culex rhabdo-like virus", ignore_case = TRUE)) & 
                              novel_flag == "Presumptive Novel" ~ "Rhabdoviridae",
                            str_detect(virus_name, regex("Yongsan picorna-like virus 1", ignore_case = TRUE)) ~ 
                              "Iflaviridae",
                            str_detect(virus_name, regex("Thrace picorna-like virus 1", ignore_case = TRUE)) ~ 
                              "Iflaviridae",
                            str_detect(virus_name, regex("Nor picorna-like virus", ignore_case = TRUE)) ~ 
                              "+ssRNA",
                            str_detect(virus_name, regex("Manitoba mononega-like virus 2", ignore_case = TRUE)) ~ "Negevirus",
                            str_detect(virus_name, regex("Manitoba mononega-like virus 1", ignore_case = TRUE)) ~ "Negevirus",
                            str_detect(virus_name, regex("Budalangi Iflavi-like virus", ignore_case = TRUE)) ~ 
                              "Iflaviridae",
                            str_detect(virus_name, regex("Hypsignathus monstrosus tombus-like virus", ignore_case = TRUE)) ~ 
                              "Tombusviridae",
                            str_detect(virus_name, regex("Tiger mosquito bi-segmented tombus-like virus", ignore_case = TRUE)) ~ 
                              "Tombusviridae",
                            str_detect(virus_name, regex("Culicine-associated Z virus", ignore_case = TRUE)) ~ "Birnaviridae",
                            str_detect(virus_name, regex("Port Bolivar virus ", ignore_case = TRUE)) ~ "Birnaviridae",
                            str_detect(virus_name, regex("Espirito Santo virus", ignore_case = TRUE)) ~ "Birnaviridae",
                            str_detect(virus_name, regex("Culicine-associated Z virus", ignore_case = TRUE)) ~ "Birnaviridae",
                            str_detect(virus_name, regex("Partitivirus-like Culex mosquito virus", ignore_case = TRUE)) ~ 
                              "PArtitiviridae",
                            str_detect(virus_name, regex("Snelk virus", ignore_case = TRUE)) ~ "Totiviridae",
                            str_detect(virus_name, regex("Cafluga virus", ignore_case = TRUE)) ~ "Cafluga virus",
                            str_detect(virus_name, regex("chuvirus", ignore_case = TRUE)) ~ "Chuviridae",
                            str_detect(virus_name, regex("Ballard Lake virus", ignore_case = TRUE)) ~ "Birnaviridae",
                            str_detect(virus_name, regex("Hanko iflavirus 1", ignore_case = TRUE)) ~ "Iflaviridae",
                            str_detect(virus_name, regex("Astopletus virus", ignore_case = TRUE)) ~ "Orthomyxoviridae",
                            str_detect(virus_name, regex("Des Moines River virus", ignore_case = TRUE)) ~ "Tombusviridae",
                            str_detect(virus_name, regex("Manitoba mononega-like virus 2", ignore_case = TRUE)) ~ "Negevirus",
                            str_detect(virus_name, regex("Manitoba picorna-like virus 2", ignore_case = TRUE)) ~ "Iflaviridae",
                            str_detect(virus_name, regex("Manitoba iflavi-like virus 1", ignore_case = TRUE)) ~ "Iflaviridae",
                            str_detect(virus_name, regex("Manitoba toti-like virus 1", ignore_case = TRUE)) ~ "Totiviridae",
                            str_detect(virus_name, regex("Manitoba picorna-like virus 3", ignore_case = TRUE)) ~ "Iflaviridae",
                            str_detect(virus_name, regex("Manitoba tymo-like virus 1", ignore_case = TRUE)) ~ "Tymoviridae",
                            str_detect(virus_name, regex("Black Queen Cell Virus", ignore_case = TRUE)) ~ "Dicistroviridae",
                            str_detect(virus_name, regex("Manitoba tombus-like virus 1", ignore_case = TRUE)) ~ "Tombusviridae",
                            str_detect(virus_name, regex("Manitoba picorna-like virus 1", ignore_case = TRUE)) ~ "Iflaviridae",
                            str_detect(virus_name, regex("Manitoba rhabdovirus 1", ignore_case = TRUE)) ~ "Rhabdoviridae",
                            str_detect(virus_name, regex("Manitoba iflavirus 1", ignore_case = TRUE)) ~ "Iflaviridae",
                            str_detect(virus_name, regex("Culex Rhabdo-like virus", ignore_case = TRUE)) ~ "Rhabdoviridae",
                            
                            
                            
  )) %>% 
  
  filter(!grepl("Homo sapiens", result_desc)) %>% 
  filter(!grepl("Human gamma", result_desc)) %>% 
  filter(!grepl("Human end", result_desc)) %>% 
  filter(!grepl("Cowpox", result_desc))

virus_master_2023 %>% 
write.csv(here("Data/virus_master_2023.csv"))

# Missing Values identifier 

missing <- virus_master_2023 %>% 
  filter(is.na(virus_name)) %>% 
  distinct(result_desc, .keep_all = TRUE)

missing_blast <- missing %>% 
  select(sample_number, contig_number, contig_sequence) %>% 
  unite(col = "id", sample_number, contig_number, sep = "-") %>% 
  rename(seq.name = id,
         seq.text = contig_sequence)

dat2fasta(missing_blast, outfile = "missing_blast.fasta")


# NOVEL ----

novel_viruses <- virus_master_2023 %>% 
  filter(grepl("Manitoba", virus_name)) %>% 
  filter(virus_name != "Manitoba virus")

novel_viruses_max <- virus_master_2023 %>% 
  filter(grepl("Manitoba", virus_name)) %>% 
  filter(virus_name != "Manitoba virus") %>% 
  group_by(virus_name) %>% 
  filter(contig_length == max(contig_length))

novel_list <- list(novel_viruses, novel_viruses_max)

write.xlsx(novel_list, file = here("Data/novel_viruses.xlsx"))