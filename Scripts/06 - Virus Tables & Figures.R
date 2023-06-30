# LIBRARIES ----
require(pacman)
p_load(tidyverse, janitor, here, gt, forcats)

# READ DATA ----
virus_master_2023 <- read_csv(here("Data/tblastx_master.csv"))


# SUMMARY TABLE ----


library_summary <- virus_master_2023 %>%
  group_by(mosquito_species, sample_number) %>%
  summarise(n = n_distinct(mosquito_species)) %>%
  group_by(mosquito_species) %>%
  summarise(n = n()) %>%
  adorn_totals()



years <- virus_master_2023 %>% 
  group_by(collection_year, virus_name) %>% 
  summarise(count = n_distinct(virus_name, collection_year)) %>% 
  pivot_wider(names_from = collection_year, values_from = count, values_fill = 0)


virus_lineage <- virus_master_2023 %>% 
  group_by(virus_name, viral_family, genome) %>% 
  summarise(n = n_distinct(virus_name, viral_family, genome)) %>% 
  select(-"n")

virus_master_2023 %>% 
  group_by(virus_name) %>% 
  summarise(n_contigs = n(),
            mean_cov = mean(coverage),
            min_cov = min(coverage),
            max_cov = max(coverage),
            mean_pid = mean(greatest_identity_percent),
            med_pid = median(greatest_identity_percent),
            min_pid = min(greatest_identity_percent),
            max_pid = max(greatest_identity_percent),
            total_reads = sum(total_read_count)) %>% 
  left_join(virus_lineage, by = "virus_name") %>% 
  group_by(genome, viral_family) %>% 
  arrange(genome, viral_family) %>% 
  relocate(total_reads, .after = "n_contigs") %>% 
  mutate(across(.cols = "virus_name",
                ~case_when(grepl("Manitoba", virus_name) & !grepl("Manitoba virus", virus_name) ~ 
                             paste0(., "*"),
                           TRUE ~ .))) %>% 
  gt() %>% 
  fmt_number(columns = mean_cov:max_pid,
             decimals = 2) %>% 
  fmt_number(columns = total_reads,
             sep_mark = ",",
             decimals = 0) %>% 
  tab_spanner(label = "Coverage Depth", columns = c(mean_cov, min_cov, max_cov)) %>% 
  tab_spanner(label = "aa Percent Identity", columns = c(mean_pid, min_pid, max_pid, med_pid)) %>% 
  cols_label(
    n_contigs = "Contigs",
    total_reads = "Reads",
    mean_cov = "Mean", min_cov = "Min", max_cov = "Max",
    mean_pid = "Mean", min_pid = "Min", max_pid = "Max", med_pid = "Median",
    virus_name = "Virus") %>% 
  tab_style(
    style = list(cell_fill(color = "grey"),
                 cell_text(weight = "bold")),
    locations = cells_row_groups(groups = everything())
  ) %>% 
  tab_style(
    style = cell_borders(
      sides = "left",
      weight = px(2),
      color = "grey"),
    locations = cells_body(
      columns = c(mean_cov, mean_pid, "n_contigs")
    )
  ) %>% 
  data_color(
    columns = mean_pid:max_pid,
    palette = "viridis"
    
  ) %>% 
  tab_footnote("* Putatively novel virus.") %>% 
  gtsave(filename = "Figures/reads_contigs_summary.html")

# Detected in mosquitoes table 
virus_master_2023 %>% 
  group_by(mosquito_species, virus_name) %>% 
  summarise(count = n_distinct(virus_name, sample_number)) %>% 
  pivot_wider(names_from = mosquito_species, values_from = count, values_fill = 0) %>% 
  mutate("Aedes canadensis \n Total: 1" = (`Aedes canadensis` / 1) * 100,
         "Aedes vexans \n Total: 19" = (`Aedes vexans` / 19) * 100,
         "Anopheles earlei \n Total: 1" = (`Anopheles earlei` / 1) * 100,
         "Coquillettidia perturbans \n Total: 6" = (`Coquillettidia perturbans` / 6) * 100,
         "Culex tarsalis \n Total: 11" = (`Culex tarsalis` / 11) * 100,
         "Ochlerotatus dorsalis \n Total: 5" = (`Ochlerotatus dorsalis` / 5) * 100,
         "Ochlerotatus flavescens \n Total: 1" = (`Ochlerotatus flavescens` / 1) * 100,
         "Ochlerotatus triseriatus \n Total: 1" = (`Ochlerotatus triseriatus` / 1) * 100,
  ) %>% select(-`Aedes canadensis`:-`Ochlerotatus triseriatus`) %>% 
  left_join(virus_lineage, by = "virus_name") %>% 
   left_join(years, by = "virus_name") %>% 
   relocate("2020", .after = "virus_name") %>% relocate("2021", .after = "2020") %>% 
   mutate(across(.cols = c("2020", "2021"),
                 ~case_when(. == "1" ~ TRUE,
                            . == "0" ~ FALSE))) %>% 
   mutate(across(.cols = "virus_name",
                 ~case_when(grepl("Manitoba", virus_name) & !grepl("Manitoba virus", virus_name) ~ 
                              paste0(., "*"),
                            TRUE ~ .))) %>% 
  group_by(genome, viral_family) %>% 
  arrange(genome, viral_family) %>% 
  gt() %>% 
  fmt_number(columns = "Aedes canadensis \n Total: 1":"Ochlerotatus triseriatus \n Total: 1",
             decimals = 2) %>% 
  tab_spanner(label = "Percent of Libraries Detected by Species", 
              columns = c("Aedes canadensis \n Total: 1":"Ochlerotatus triseriatus \n Total: 1")) %>% 
   tab_spanner(label = html("Years Detected"),
               columns = c("2020":"2021")) %>% 
  cols_label(
    "Aedes canadensis \n Total: 1" = html("<em>Aedes canadensis</em> <br> Total: 1"),
    "Aedes vexans \n Total: 19" = html("<em>Aedes vexans</em> <br> Total: 19"),
    "Anopheles earlei \n Total: 1" = html("<em>Anopheles earlei</em> <br> Total: 1"),
    "Coquillettidia perturbans \n Total: 6" = html("<em>Coquillettidia perturbans</em> <br> Total: 6"),
    "Culex tarsalis \n Total: 11" = html("<em>Culex tarsalis</em> <br> Total: 11"),
    "Ochlerotatus dorsalis \n Total: 5" = html("<em>Ochlerotatus dorsalis</em> <br> Total: 5"),
    "Ochlerotatus flavescens \n Total: 1" = html("<em>Ochlerotatus flavescens</em> <br> Total: 1"),
    "Ochlerotatus triseriatus \n Total: 1" = html("<em>Ochlerotatus triseriatus</em> <br> Total: 1"),
    virus_name = "Virus") %>% 
   cols_align(align = "center",
              columns = "2020":"Ochlerotatus triseriatus \n Total: 1") %>% 
   tab_style(
     style = list(cell_fill(color = "grey"),
                  cell_text(weight = "bold")),
     locations = cells_row_groups(groups = everything())) %>% 
   data_color(columns = `Aedes canadensis \n Total: 1`:last_col(),
              palette = "viridis", 
              direction = "row") %>% 
   tab_style(
     style = cell_borders(
       sides = "left",
       weight = px(2),
       color = "black"),
     locations = cells_body(
       columns = c("Aedes canadensis \n Total: 1"))) %>% 
  tab_footnote("* Putatively novel virus") %>% 
  gtsave(filename = "Figures/virus_detections_by_species_summary.html")

# Organism figs ----

virus_master_2023 %>% 
  distinct(virus_name, mosquito_species, .keep_all = TRUE) %>% 
  group_by(viral_family, genome, mosquito_species) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = reorder(factor(viral_family), -n), y = n, fill = genome)) +
  geom_col() +
  theme_bw(base_size = 14) +
  scale_fill_viridis_d("Genome") +
  scale_y_continuous(name = "Number of Viruses", breaks = c(1,5,10)) + 
  facet_wrap(~mosquito_species, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 14),
        strip.text = element_text(face = "bold.italic"),
        panel.border = element_rect(fill = NA, color = "gray90"),
        axis.title.y = element_blank(),
        legend.position = c(0.8, 0.1)) +
  labs(x = "")



  
  




# SEQ FIGS ----

# Read count & contigs
virus_master_2023 %>% 
  mutate(reads_log10 = log10(total_read_count)) %>% 
  drop_na(genome) %>% 
  ggplot(aes(x = contig_length, y = reads_log10)) +
  geom_point(aes(colour = genome)) +
  scale_colour_viridis_d("Genome") +
  theme_bw() +
  labs(x = "Contig Length (nt)",
       y = "Reads (Log 10)")

# Read count & contigs (mosquito species)

virus_master_2023 %>% 
  mutate(reads_log10 = log10(total_read_count)) %>% 
  drop_na(genome) %>% 
  ggplot(aes(x = contig_length, y = reads_log10)) +
  geom_point(aes(colour = mosquito_species)) +
  scale_colour_viridis_d("Species") +
  theme_bw() +
  labs(x = "Contig Length (nt)",
       y = "Reads (Log 10)")

# Contig length and coverage
virus_master_2023 %>% 
  mutate(coverage_log10 = log10(coverage)) %>% 
  ggplot(aes(x = contig_length, y = coverage_log10)) +
  geom_point(aes(colour = genome)) +
  scale_colour_viridis_d() +
  theme_bw()

# Coverage reads

virus_master_2023 %>% 
  mutate(coverage_log10 = log10(coverage)) %>% 
  mutate(reads_log10 = log10(total_read_count)) %>% 
  ggplot(aes(x = coverage_log10, y = reads_log10)) +
  geom_point(aes(colour = genome)) +
  scale_colour_viridis_d("Genome") +
  theme_bw() +
  labs(x = "Coverage Depth (Log 10)",
       y = "Reads (Log 10)")

# Contig Length and % ID
virus_master_2023 %>% 
  ggplot(aes(x = contig_length, y = greatest_identity_percent)) +
  geom_point(aes(colour = genome)) +
  scale_colour_viridis_d() +
  theme_bw()

# VIRUSES ----

# Overall


# Family

viruses <- virus_master_2023 %>% 
  distinct(mosquito_species, virus_name, viral_family) %>% 
  group_by(mosquito_species, viral_family) %>% 
  summarise(n = n()) %>% 
  adorn_totals()

