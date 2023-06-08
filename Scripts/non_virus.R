require(pacman)
pacman::p_load(tidyverse, readxl, janitor, here, ggtext, shadowtext)

# FUNGI ----

# DATA 

fungi <- read_excel(here("Raw/Meta2.xlsx"),
                        sheet = "Fungi") %>% 
  row_to_names(1) %>% clean_names() %>% 
  mutate(sample_number = as.integer(sample_number)) %>% 
  mutate(nt_reads = as.numeric(nt_reads)) %>% 
  select(-"na", -"depth", -"nt_rpm", -"nr_rpm":-last_col()) %>% 
  drop_na(sample_number)

fungi_reads_summary <- fungi %>% 
  group_by(mosquito_species, fungal_division_or_phylum) %>% 
  summarise(total_reads = sum(nt_reads))

distinct_organism_summary <- fungi %>% 
  group_by(mosquito_species, fungal_division_or_phylum) %>% 
  summarise(total_unique_organisms = n_distinct(name))

total_organism_summary <- fungi %>% 
  group_by(fungal_division_or_phylum) %>% 
  summarise(total_unique_organisms = n_distinct(name))

fungi_sample <- fungi %>% 
  group_by(mosquito_species, fungal_division_or_phylum) %>% 
  summarise(total_unique_organisms = n_distinct(name))

# PLOTS

distinct_fungi_organism_plot <- distinct_organism_summary %>% 
  ggplot(aes(x = reorder(fungal_division_or_phylum, -total_unique_organisms, sum), 
             y = total_unique_organisms, fill = mosquito_species, label = total_unique_organisms)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d("Host Species", labels = c("*Ae. vexans*", "*An. earlei*", "*Cq. perturbans*", "*Cx. tarsalis*", "*Oc. dorsalis*")) +
  geom_shadowtext(size = 4, fontface = "bold", position = position_dodge(width = 0.9)) +
  labs(x = "Fungal Taxon",
       y = "Total Distinct Organisms") +
  theme_bw() +
  theme(legend.text = element_markdown(),
        legend.position = c(0.9, 0.3),
        legend.background = element_rect(linetype = 2, linewidth = 0.5, colour = 1))

tiff(here("Figures/distinct_fungi_organisms.tiff"), units="in", width=14, height=7, res=300)
distinct_fungi_organism_plot
dev.off()

fungi_reads_plot <- fungi_reads_summary %>% 
  mutate(total_reads_log = round(log10(total_reads), digits = 2)) %>% 
  ggplot(aes(x = reorder(fungal_division_or_phylum, -total_reads, sum), 
             y = total_reads_log, fill = mosquito_species)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d("Host Species", labels = c("*Ae. vexans*", "*An. earlei*", "*Cq. perturbans*", "*Cx. tarsalis*", "*Oc. dorsalis*")) +
  labs(x = "Fungal Taxon",
       y = "Total Reads (Log 10)") +
  theme_bw() +
  theme(legend.text = element_markdown(),
        legend.position = c(0.9, 0.7),
        legend.background = element_rect(linetype = 2, linewidth = 0.5, colour = 1))

tiff(here("Figures/fungi_reads_plot.tiff"), units="in", width=14, height=7, res=300)
fungi_reads_plot
dev.off()

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
  summarise(total_reads = sum(nt_reads))

distinct_parasites_summary <- parasites %>% 
  group_by(mosquito_species, phylum) %>% 
  summarise(total_unique_organisms = n_distinct(name))

total_parasites_summary <- parasites %>% 
  group_by(mosquito_species, phylum) %>% 
  summarise(total_unique_organisms = n())


# PLOT

distinct_parasite_organisms_plot <- distinct_parasites_summary %>% 
  ggplot(aes(x = reorder(phylum, -total_unique_organisms, sum), 
             y = total_unique_organisms, fill = mosquito_species, label = total_unique_organisms)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d("Host Species", option = "cividis", labels = c("*Ae. vexans*", "*Cq. perturbans*", "*Cx. tarsalis*", "*Oc. dorsalis*", "*Oc. flavescens*")) +
  geom_shadowtext(size = 4, fontface = "bold", position = position_dodge(width = .9)) +
  labs(x = "Parasite Taxon",
       y = "Total Distinct Organisms") +
  theme_bw() +
  theme(legend.text = element_markdown(),
        legend.position = c(0.9, 0.3),
        legend.background = element_rect(linetype = 2, linewidth = 0.5, colour = 1))

tiff(here("Figures/distinct_parasite_organisms.tiff"), units="in", width=14, height=7, res=300)
distinct_parasite_organisms_plot
dev.off()

parasite_reads_plot <- parasites_reads_summary %>% 
  mutate(total_reads_log = round(log10(total_reads), digits = 2)) %>% 
  ggplot(aes(x = reorder(phylum, -total_reads, sum), 
             y = total_reads_log, fill = mosquito_species)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d("Host Species", option = "cividis",labels = c("*Ae. vexans*", "*Cq. perturbans*", "*Cx. tarsalis*", "*Oc. dorsalis*", "*Oc. flavescens*")) +
  labs(x = "Parasite Taxon",
       y = "Total Reads (Log 10)") +
  theme_bw() +
  theme(legend.text = element_markdown(),
        legend.position = c(0.9, 0.7),
        legend.background = element_rect(linetype = 2, linewidth = 0.5, colour = 1))

tiff(here("Figures/parasites_reads_plot.tiff"), units="in", width=14, height=7, res=300)
parasite_reads_plot
dev.off()

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
  summarise(total_reads = sum(nt_reads))

distinct_plants_summary <- plants %>% 
  group_by(mosquito_species, order) %>% 
  summarise(total_unique_organisms = n_distinct(name))

total_plants_summary <- plants %>% 
  group_by(mosquito_species, order) %>% 
  summarise(total_unique_organisms = n())

# PLOT

distinct_plant_organisms_plot <- distinct_plants_summary %>% 
  ggplot(aes(x = reorder(order, -total_unique_organisms, sum), 
             y = total_unique_organisms, fill = mosquito_species, label = total_unique_organisms)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d("Host Species", option = "turbo", labels = c("*Ae. vexans*", "*An. earlei*", "*Cq. perturbans*", "*Cx. tarsalis*")) +
  geom_shadowtext(size = 4, fontface = "bold", position = position_dodge(width = .9)) +
  labs(x = "Plant Taxon",
       y = "Total Distinct Organisms") +
  theme_bw() +
  theme(legend.text = element_markdown(),
        legend.position = c(0.9, 0.5),
        legend.background = element_rect(linetype = 2, linewidth = 0.5, colour = 1))

tiff(here("Figures/distinct_plant_organisms.tiff"), units="in", width=14, height=7, res=300)
distinct_plant_organisms_plot
dev.off()

plant_reads_plot <- plant_reads_summary %>% 
  mutate(total_reads_log = round(log10(total_reads), digits = 2)) %>% 
  ggplot(aes(x = reorder(order, -total_reads, sum), 
             y = total_reads_log, fill = mosquito_species)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d("Host Species", option = "turbo",labels = c("*Ae. vexans*", "*An. earlei*", "*Cq. perturbans*", "*Cx. tarsalis*")) +
  labs(x = "Plant Taxon",
       y = "Total Reads (Log 10)") +
  theme_bw() +
  theme(legend.text = element_markdown(),
        legend.position = c(0.9, 0.7),
        legend.background = element_rect(linetype = 2, linewidth = 0.5, colour = 1))

tiff(here("Figures/plant_reads_plot.tiff"), units="in", width=14, height=7, res=300)
plant_reads_plot
dev.off()

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
  summarise(total_reads = sum(nt_reads))

distinct_vert_summary <- vert %>% 
  group_by(mosquito_species, type) %>% 
  summarise(total_unique_organisms = n_distinct(name))

total_vert_summary <- vert %>% 
  group_by(mosquito_species, type) %>% 
  summarise(total_unique_organisms = n())

# PLOT

distinct_vert_organisms_plot <- distinct_vert_summary %>% 
  ggplot(aes(x = reorder(type, -total_unique_organisms, sum), 
             y = total_unique_organisms, fill = mosquito_species, label = total_unique_organisms)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d("Host Species", option = "turbo", labels = c("*An. earlei*", "*Cx. tarsalis*")) +
  geom_shadowtext(size = 4, fontface = "bold", position = position_dodge(width = .9)) +
  labs(x = "Vertebrate Taxon",
       y = "Total Distinct Organisms") +
  theme_bw() +
  theme(legend.text = element_markdown(),
        legend.position = c(0.7, 0.5),
        legend.background = element_rect(linetype = 2, linewidth = 0.5, colour = 1))

tiff(here("Figures/distinct_vertebrate_organisms.tiff"), units="in", width=14, height=7, res=300)
distinct_vert_organisms_plot
dev.off()

vert_reads_plot <- vert_reads_summary %>% 
  mutate(total_reads_log = round(log10(total_reads), digits = 2)) %>% 
  ggplot(aes(x = reorder(type, -total_reads, sum), 
             y = total_reads_log, fill = mosquito_species)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d("Host Species", option = "turbo",labels = c("*An. earlei*", "*Cx. tarsalis*")) +
  labs(x = "Plant Taxon",
       y = "Total Reads (Log 10)") +
  theme_bw() +
  theme(legend.text = element_markdown(),
        legend.position = c(0.7, 0.9),
        legend.background = element_rect(linetype = 2, linewidth = 0.5, colour = 1))

tiff(here("Figures/vert_reads_plot.tiff"), units="in", width=14, height=7, res=300)
vert_reads_plot
dev.off()


