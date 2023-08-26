require(pacman)
p_load(httr, jsonlite, here, tidyverse, openxlsx, countrycode, 
       maps, sf, gt, assertthat, igraph, ggraph, ggmap, 
       shadowtext, rgeos, rworldmap, ggflags)

# Scrape Data ----

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

tblastx %>% 
  distinct(virus_name, accession, .keep_all = TRUE) %>% 
  group_by(novel_flag, country_simple) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  gt() %>% 
  cols_label(novel_flag = "Novel",
             country_simple = "Region")

mosquito_summary <- tblastx %>% 
  separate_wider_delim(cols = c("mosquito_species"), delim = " ", names = c("genus", "species")) %>% 
  mutate(genus_consensus = case_when(genus == genus_accession ~ "Same Genus",
                                     TRUE ~ "Different Genus")) %>% 
  # distinct(virus_name, accession, .keep_all = TRUE) %>%
  group_by(novel_flag, genus_consensus) %>% 
  summarise(n = n())


# World Map Fill ----

mosq_world <- tblastx %>% 
  distinct(virus_name, accession, country, .keep_all = TRUE) %>% 
  mutate(virus_name = str_replace(virus_name, "Manitoba virus", "Manitoba_virus")) %>% 
  drop_na(country) %>% 
  # mutate(iso = countrycode(country, origin = "country.name", destination = "iso2c")) %>% 
  group_by(country) %>%
  summarise(value = n())
# rename(country = iso)



# NON NOVEL
  
mosq_world <- tblastx %>% 
  distinct(virus_name, accession, country, .keep_all = TRUE) %>% 
  mutate(virus_name = str_replace(virus_name, "Manitoba virus", "Manitoba_virus")) %>% 
  filter(., !grepl("Manitoba ", virus_name)) %>% 
  drop_na(country) %>% 
  # mutate(iso = countrycode(country, origin = "country.name", destination = "iso2c")) %>% 
  group_by(country) %>% 
  summarise(value = n())
  # rename(country = iso)


# World map data (without Antarctica)
world <- map_data("world")
world <- subset(world, region != "Antarctica")

# Create the plot
non_novel_countries <- world %>% 
  merge(mosq_world, by.x = "region", by.y = "country", all.x = T) %>% 
  arrange(group, order) %>% 
  ggplot(aes(x = long, y = lat, group = group, fill = value)) + geom_polygon() +
  theme_void() +
  scale_fill_viridis_c("Number of Distinct \nVirus-Accession Combinations") +
  theme(legend.position = c(0.2, 0.2))

# NOVEL

mosq_world <- tblastx %>% 
  distinct(virus_name, country, .keep_all = TRUE) %>% 
  mutate(virus_name = str_replace(virus_name, "Manitoba virus", "Manitoba_virus")) %>% 
  filter(., grepl("Manitoba ", virus_name)) %>% 
  drop_na(country) %>% 
  # mutate(iso = countrycode(country, origin = "country.name", destination = "iso2c")) %>% 
  group_by(country) %>% 
  summarise(value = n())
# rename(country = iso)


# World map data (without Antarctica)
world <- map_data("world")
world <- subset(world, region != "Antarctica")

# Create the plot
novel_countries <- world %>% 
  merge(mosq_world, by.x = "region", by.y = "country", all.x = T) %>% 
  arrange(group, order) %>% 
  ggplot(aes(x = long, y = lat, group = group, fill = value)) + geom_polygon() +
  theme_void() +
  scale_fill_viridis_c("Number of Distinct \nVirus-Accession Combinations") +
  theme(legend.position = c(0.2, 0.2)) 

# SPECIES ----

tblastx %>% 
  select(virus_name, mosquito_species, country, genus, species) %>% 
  summarise(n_contigs = n(), .by = c(virus_name, mosquito_species, country, genus, species)) %>% 
  arrange(virus_name) %>% 
  mutate(organism = case_when(genus %in% c("Aedes", "Culex", "Anopheles", "Culiseta", "Ochlerotatus", "mosquito", 
                                           "Mosquitoes", "mosquitoes", "Mosquito", "culicine", "Culicidae") ~ "Mosquito",
                              genus %in% c("Grus", "Branta", "Ramphocelus") ~ "Bird",
                              genus %in% c("Apis") ~ "Bee",
                              genus %in% c("Neohydatothrips") ~ "Thrips",
                              genus %in% c("spiders") ~ "Spiders",
                              genus %in% c("Homalodisca") ~ "Hemiptera",
                              TRUE ~ NA)) %>% 
  relocate(organism, .after = species) %>% 
  relocate(n_contigs, .after = mosquito_species) %>% 
  rename("Mosquito Species (Canada)" = mosquito_species,
         "Country" = country,
         "Genus" = genus,
         "Species" = species,
         "Organism" = organism,
         "Number of Contigs" = n_contigs) %>% 
  ungroup() %>% 
  gt() %>% 
  tab_spanner(label = "Closest BLAST Relative", columns = c(Country, Genus, Species, Organism))



df %>% 
  write.xlsx(here("hit_metadata.xlsx"))

# Network ----

# Get centroid coordinates for each country
wmap <- getMap(resolution="high")
centroids <- gCentroid(wmap, byid=TRUE)
centroid_df <- as.data.frame(centroids) %>% 
  rownames_to_column(., "name") %>% 
  rename(lon = x) %>% rename(lat = y)

edges <- tblastx %>% 
  distinct(virus_name, accession, .keep_all = TRUE) %>% 
  drop_na(country) %>% 
  mutate(country = gsub("USA", "United States of America", country)) %>% 
  mutate(country = gsub("Serbia", "Republic of Serbia", country)) %>% 
  left_join(centroid_df, by = c("country" = "name")) %>% 
  select(novel_flag, organism, country, lat, lon) %>% 
  rename(name = country) %>% 
  group_by(name) %>% 
  mutate(id = cur_group_id()) %>% ungroup()
 

nodes <- edges %>% 
  distinct(name, .keep_all = TRUE) %>% 
  select(-novel_flag, -organism) %>% 
  relocate(., name, .after = id) %>% 
  relocate(., id, .before = lat) %>% 
  mutate(iso2_lower = tolower(countrycode(name, origin = "country.name", destination = "iso2c")))

edges <- edges %>% 
  select(-name) %>% 
  select(-lat, -lon) %>% 
  mutate(to = 4) %>% 
  rename(from = id) %>% 
  rename(category = novel_flag) %>% 
  mutate(category = as.factor(category)) %>% 
  mutate(weight = runif(length(to))) %>% 
  ungroup() %>% 
  relocate(., category, .after = weight) %>% 
  mutate_at(c(1), ~replace(., is.na(.), "No Data")) %>% 
  relocate(., organism, .after = category)
# Remove Canada 
edges <- edges[edges$from != edges$to, ]




g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)

edges_for_plot <- edges %>%
  inner_join(nodes %>% select(id, lon, lat), by = c('from' = 'id')) %>%
  rename(x = lon, y = lat) %>%
  inner_join(nodes %>% select(id, lon, lat), by = c('to' = 'id')) %>%
  rename(xend = lon, yend = lat)

assert_that(nrow(edges_for_plot) == nrow(edges))

nodes$weight = degree(g)

maptheme <- theme(panel.grid = element_blank()) +
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(legend.position = c(0.8, 0.1)) +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_rect(fill = "cadetblue")) +
  theme(plot.margin = unit(c(0, 0, 0.5, 0), 'cm')) +
  theme(legend.key.size = unit(0.5, 'cm')) +
  theme(legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

country_shapes <- geom_polygon(aes(x = long, y = lat, group = group),
                               data = map_data('world'),
                               fill = "#CECECE", color = "#515151",
                               size = 0.15)
mapcoords <- coord_fixed(xlim = c(-150, 180), ylim = c(-55, 80))


na_flags <- nodes %>% 
  filter(iso2_lower %in% c("us", "ca"))

non_na_flags <- nodes %>% 
  filter(!iso2_lower %in% c("us", "ca", "fi"))

fin <- nodes %>% 
  filter(iso2_lower %in% c("fi"))

network_map <- ggplot(nodes) + country_shapes +
  geom_curve(aes(x = x, y = y, xend = xend, yend = yend,     
                 color = category, size = weight),
             data = edges_for_plot, curvature = 0.25,
             alpha = 0.5) +
  scale_colour_manual("Noveltly", values = c("Not Novel" = "#440154", "Novel" = "#fde725")) +
  scale_size_continuous(guide = FALSE, range = c(0.25, 2)) + 
  geom_point(aes(x = lon, y = lat, size = weight),           
             shape = 21, fill = 'white',
             color = 'black', stroke = 0.5) +
  scale_size_continuous(guide = FALSE, range = c(1, 14.5)) +  
  ggflags::geom_flag(data = non_na_flags, aes(x = lon+3, y = lat, country = iso2_lower)) +
  ggflags::geom_flag(data = na_flags, aes(x = lon-6, y = lat, country = iso2_lower)) +
  ggflags::geom_flag(data = fin, aes(x = lon+3.5, y = lat, country = iso2_lower)) +
  # geom_shadowtext(aes(x = lon, y = lat, label = name),             
  #           hjust = 0, nudge_x = 0, nudge_y = 4,
  #           size = 4, color = "white", fontface = "bold") +
  mapcoords + maptheme +
  guides(colour = guide_legend(override.aes = list(linewidth=3)))

tiff(here("network_map.tiff"), units="in", width=20, height=10, res=1000)
network_map
dev.off()


