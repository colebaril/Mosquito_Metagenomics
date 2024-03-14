require(pacman)
p_load(httr, jsonlite, here, tidyverse, openxlsx, countrycode, 
       maps, sf, gt, assertthat, igraph, ggraph, ggmap, 
       shadowtext, rgeos, rworldmap, ggflags, png, patchwork, tiff)

# Read accession numbers from CSV
accession_numbers <- read.csv(here("Accession/accession_numbers.csv"), stringsAsFactors = FALSE)$AccessionNumber

# Iterate through each accession number and fetch metadata - run once then comment
# for (accession_number in accession_numbers) {
#   url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id=", accession_number, "&retmode=json")
#   
#   response <- GET(url)
#   
#   if (status_code(response) == 200) {
#     metadata_json <- content(response, "text")
#     metadata <- fromJSON(metadata_json)
#     
#     filename <- paste0("Accession/", accession_number, "_metadata.json")
#     write(jsonlite::toJSON(metadata), file = filename)
#     
#     cat("Metadata saved to", filename, "\n")
#   } else {
#     cat("Failed to fetch metadata for", accession_number, "\n")
#   }
#   
#   Sys.sleep(1) # Rate limit once per second 
# }

# Initialize lists to store extracted values
extracted_subnames <- list()
extracted_subtypes <- list()

json_files <- list.files(pattern = ".json", path = here("Accession/"), full.names = TRUE)

for (file_path in json_files) {
  json_data <- fromJSON(file_path, flatten = TRUE)  # Load JSON data
  
  # Extract subname and subtype values
  subtype_value <- json_data$result[[names(json_data$result)[2]]]$subtype
  subname_value <- json_data$result[[names(json_data$result)[2]]]$subname
  
  # Store the values in the lists
  extracted_subtypes[[file_path]] <- subtype_value
  extracted_subnames[[file_path]] <- subname_value
}

# Create a list to store data frames
data_frames <- list()

for (i in seq_along(json_files)) {
  subname <- extracted_subnames[[json_files[i]]]
  subtype <- extracted_subtypes[[json_files[i]]]
  
  # Split values and names
  subname_values <- unlist(strsplit(subname, "\\|"))
  subtype_values <- unlist(strsplit(subtype, "\\|"))
  
  # Combine into a data frame
  data_frames[[i]] <- data.frame(subtype = subtype_values, subname = subname_values, file_name = basename(json_files[i]))
}

# Combine all data frames into one
combined_df <- bind_rows(data_frames)

# Pivot the data
pivot_df <- pivot_wider(combined_df, names_from = subtype, values_from = subname)

df <- pivot_df %>% 
  rename(Accession = file_name) %>% 
  mutate(Accession = gsub("_metadata.json", "", Accession)) %>% 
  mutate(host = coalesce(host, isolation_source)) %>% 
  mutate(host = coalesce(host, cell_line)) %>% 
  select(-"segment":-last_col()) %>% 
  select(-"strain":-"isolate", -"isolation_source") %>% 
  separate_wider_delim(cols = host, delim = " ", names = c("genus_accession", "species_accession"),
                       too_few = "align_start", too_many = "merge") %>% 
  mutate(across(c("species_accession"),
                ~replace(., str_detect(., "sp.|sp|pool"), NA))) %>% 
  mutate(organism = case_when(genus_accession %in% c("Aedes", "Culex", "Anopheles", "Culiseta", "Ochlerotatus", "mosquito", 
                                                     "Mosquitoes", "mosquitoes", "Mosquito", "culicine", "Culicidae") ~ "Mosquito",
                              genus_accession %in% c("Grus", "Branta", "Ramphocelus") ~ "Bird",
                              genus_accession %in% c("Apis") ~ "Bee",
                              genus_accession %in% c("Neohydatothrips") ~ "Thrips",
                              genus_accession %in% c("spiders") ~ "Spiders",
                              genus_accession %in% c("Homalodisca") ~ "Hemiptera",
                              TRUE ~ NA)) %>% 
  relocate(organism, .after = species_accession) %>% 
  mutate(across(c("country"),
                ~replace(., str_detect(., "USA"), "USA"))) %>% 
  mutate(across(c("country"),
                ~replace(., str_detect(., "Brazil"), "Brazil"))) %>% 
  mutate(across(c("country"),
                ~replace(., str_detect(., "Australia"), "Australia"))) %>% 
  mutate(across(c("country"),
                ~replace(., str_detect(., "China"), "China"))) %>% 
  mutate(across(c("country"),
                ~replace(., str_detect(., "Canada"), "Canada"))) %>% 
  mutate(across(c("country"),
                ~replace(., str_detect(., "South Korea"), "South Korea"))) %>%
  mutate(across(c("country"),
                ~replace(., str_detect(., "Russia"), "Russia"))) %>% 
  mutate(across(c("country"),
                ~replace(., str_detect(., "Finland"), "Finland"))) 

tblastx <- read_csv("C:/Projects/Mosquito Metagenomics/Data/tblastx_master.csv") %>% 
  left_join(df, by = c("accession" = "Accession")) %>% 
  mutate(organism_simple = case_when(organism == "Mosquito" ~ "Mosquito",
                                     TRUE ~ "Non-Mosquito")) %>% 
  mutate(country_simple = case_when(country %in% c("USA", "Canada", "Mexico") ~ "North America",
                                    TRUE ~ "Outside America")) %>% 
  mutate(virus_name = str_replace(virus_name, "Manitoba virus", "Manitoba_virus")) %>% 
  mutate(novel_flag = case_when(str_detect(virus_name, "Manitoba ") ~ "Novel",
                                TRUE ~ "Not Novel"))

# gt Table for summary 
tblastx %>% 
  distinct(virus_name, country, .keep_all = TRUE) %>% 
  drop_na(country) %>% 
  mutate(country = gsub("USA", "United States of America", country)) %>% 
  mutate(country = gsub("Serbia", "Republic of Serbia", country)) %>% 
  # left_join(centroid_df, by = c("country" = "name")) %>% 
  select(virus_name, novel_flag, organism, country) %>% 
  summarise(n = n(), .by = c(novel_flag, country)) %>% 
  pivot_wider(names_from = novel_flag, values_from = n) %>% 
  mutate(country_code = countrycode(country, "country.name", "iso2c")) %>% 
  relocate(country_code, .before = country) %>% 
  arrange(country) %>% 
  gt() %>% 
  fmt_flag(., columns = country_code) %>% 
  cols_label(country = "Country",
             country_code = "") %>% 
  cols_align(., align = c("center"), columns = c("Not Novel", "Novel")) %>% 
  gtsave(., here("country_virus_table.html"))