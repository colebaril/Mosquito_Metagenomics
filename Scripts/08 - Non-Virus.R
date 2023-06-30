require(pacman)
pacman::p_load(tidyverse, readxl, janitor, here, ggtext, shadowtext, patchwork)

# FUNGI ----

# DATA 

fungi <- read_excel(here("Raw/Meta2.xlsx"),
                        sheet = "Fungi") %>% 
  row_to_names(1) %>% clean_names() %>% 
  mutate(sample_number = as.integer(sample_number)) %>% 
  mutate(nt_reads = as.numeric(nt_reads)) %>% 
  select(-"na", -"depth", -"nt_rpm", -"nr_rpm":-last_col()) %>% 
  drop_na(sample_number)

# plot data
fungi_reads_summary <- fungi %>% 
  group_by(mosquito_species, fungal_division_or_phylum) %>% 
  summarise(total_reads = sum(nt_reads)) %>% 
  rename("division" = fungal_division_or_phylum) %>% 
  mutate(group = "Fungi")

distinct_organism_summary <- fungi %>% 
  group_by(mosquito_species, fungal_division_or_phylum) %>% 
  summarise(total_unique_organisms = n_distinct(name)) %>% 
  rename("division" = fungal_division_or_phylum) %>% 
  mutate(group = "Fungi")
# End 

total_organism_summary <- fungi %>% 
  group_by(fungal_division_or_phylum) %>% 
  summarise(total_unique_organisms = n_distinct(name))

fungi_sample <- fungi %>% 
  group_by(mosquito_species, fungal_division_or_phylum) %>% 
  summarise(total_unique_organisms = n_distinct(name))

fungi %>% 
  group_by(fungal_division_or_phylum) %>% 
  summarise(n = n_distinct(name)) %>% 
  adorn_totals()

# PARASITES ----

# DATA

parasites <- read_excel(here("Raw/Meta2.xlsx"),
                    sheet = "Parasites") %>% 
  row_to_names(1) %>% clean_names() %>% 
  mutate(sample_number = as.integer(sample_number)) %>% 
  mutate(nt_reads = as.numeric(nt_reads)) %>% 
  select(-"na", -"na_2", -"coverage_depth", -"nt_rpm", -"nr_rpm":-last_col()) %>% 
  drop_na(sample_number)

parasites_reads_summary <- parasites %>% 
  group_by(mosquito_species, phylum) %>% 
  summarise(total_reads = sum(nt_reads)) %>% 
  rename("division" = phylum) %>% 
  mutate(group = "Parasites")

distinct_parasites_summary <- parasites %>% 
  group_by(mosquito_species, phylum) %>% 
  summarise(total_unique_organisms = n_distinct(name)) %>% 
  rename("division" = phylum) %>% 
  mutate(group = "Parasites")

total_parasites_summary <- parasites %>% 
  group_by(mosquito_species, phylum) %>% 
  summarise(total_unique_organisms = n())

# PLANTS ----

plants <- read_excel(here("Raw/Meta2.xlsx"),
                        sheet = "Plant") %>% 
  row_to_names(1) %>% clean_names() %>% 
  mutate(sample_number = as.integer(sample_number)) %>% 
  mutate(nt_reads = as.numeric(nt_reads)) %>% 
  select(-"na", -"na_2", -"nt_rpm", -"nr_rpm":-last_col()) %>% 
  drop_na(sample_number)

plant_reads_summary <- plants %>% 
  group_by(mosquito_species, order) %>% 
  summarise(total_reads = sum(nt_reads)) %>% 
  rename("division" = order) %>% 
  mutate(group = "Plants")

distinct_plants_summary <- plants %>% 
  group_by(mosquito_species, order) %>% 
  summarise(total_unique_organisms = n_distinct(name)) %>% 
  rename("division" = order) %>% 
  mutate(group = "Plants")

total_plants_summary <- plants %>% 
  group_by(mosquito_species, order) %>% 
  summarise(total_unique_organisms = n())

# CHORDATES ----

vert <- read_excel(here("Raw/Meta2.xlsx"),
                     sheet = "Vertebrate") %>% 
  row_to_names(1) %>% clean_names() %>% 
  mutate(sample_number = as.integer(sample_number)) %>% 
  mutate(nt_reads = as.numeric(nt_reads)) %>% 
  select(-"na", -"na_2", -"nt_rpm", -"nr_rpm":-last_col()) %>% 
  drop_na(sample_number)

vert_reads_summary <- vert %>% 
  group_by(mosquito_species, type) %>% 
  summarise(total_reads = sum(nt_reads)) %>% 
  rename("division" = type) %>% 
  mutate(group = "Vertebrates")

distinct_vert_summary <- vert %>% 
  group_by(mosquito_species, type) %>% 
  summarise(total_unique_organisms = n_distinct(name)) %>% 
  rename("division" = type) %>% 
  mutate(group = "Vertebrates")

total_vert_summary <- vert %>% 
  group_by(mosquito_species, type) %>% 
  summarise(total_unique_organisms = n())


# PLOTS ----

combined_distinct_organism <- rbind(distinct_organism_summary, distinct_parasites_summary, 
                                    distinct_plants_summary, distinct_vert_summary)

combined_reads <- rbind(fungi_reads_summary, parasites_reads_summary, plant_reads_summary,
                        vert_reads_summary)

distinct_organisms <- combined_distinct_organism %>% 
  mutate(mosquito_species = gsub("An. earlei", "Anopheles earlei", mosquito_species),
         mosquito_species = gsub("Cx. tarsalis", "Culex tarsalis", mosquito_species)) %>% 
  ggplot(aes(x = reorder(division, -total_unique_organisms, sum), 
             y = total_unique_organisms, fill = mosquito_species, label = total_unique_organisms)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d("Host Species", option = "viridis") +
  geom_shadowtext(size = 4, fontface = "bold", position = position_dodge(width = .9)) +
  labs(x = "",
       y = "Total Distinct Organisms") +
  facet_wrap(~group, scales = "free") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        legend.position = c(0.9, 0.8),
        legend.background = element_rect(linetype = 2, linewidth = 0.5, colour = 1))

tiff(here("Figures/non-virus_distinct_organisms.tiff"), units="in", width=16, height=10, res=300)
distinct_organisms 
dev.off()


combined_reads_non_virus <- combined_reads %>% 
  mutate(mosquito_species = gsub("An. earlei", "Anopheles earlei", mosquito_species),
         mosquito_species = gsub("Cx. tarsalis", "Culex tarsalis", mosquito_species)) %>% 
  ggplot(aes(x = reorder(division, -total_reads, sum), 
             y = log10(total_reads), fill = mosquito_species)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d("Host Species", option = "viridis") +
  labs(x = "",
       y = "Total Reads (Log 10)") +
  facet_wrap(~group, scales = "free") +
  theme_bw() +
  theme(legend.text = element_text(face = "italic"),
        legend.position = c(0.93, 0.89),
        legend.background = element_rect(linetype = 2, linewidth = 0.5, colour = 1))

tiff(here("Figures/non-virus_reads.tiff"), units="in", width=16, height=10, res=300)
combined_reads_non_virus
dev.off()



# TABLES ----




