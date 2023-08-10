# LIBRARIES ----
require(pacman)
p_load(tidyverse, janitor, here, gt, forcats, ggVennDiagram, tiff, openxlsx)

# READ DATA ----
virus_master_2023 <- read_csv(here("Data/tblastx_master.csv"))

# READS SUMMARIES ----

# Reads per species/genus

species_reads <- virus_master_2023 %>% 
  group_by(mosquito_species) %>% 
  summarise(sum_viral_reads = sum(total_read_count)) %>% 
  adorn_totals()

genus_reads <- virus_master_2023 %>% 
  separate_wider_delim(., cols = mosquito_species, delim = " ", names = c("genus", "species")) %>% 
  group_by(genus) %>% 
  summarise(sum_viral_reads = sum(total_read_count)) %>% 
  adorn_totals()

# Reads per genome / viral family 

genome_reads <- virus_master_2023 %>% 
  group_by(genome) %>% 
  summarise(sum_viral_reads = sum(total_read_count)) %>% 
  adorn_totals()

viral_family_reads <- virus_master_2023 %>% 
  group_by(viral_family) %>% 
  summarise(sum_viral_reads = sum(total_read_count)) %>% 
  adorn_totals()

# Save 

reads_list <- list(species_reads, genus_reads, genome_reads, viral_family_reads)

write.xlsx(reads_list, file = here("Data/viral_reads_summary.xlsx"))


# SUMMARY TABLE ----

# BLAST SUMMARY ----

# NON NOVEL

library_summary <- virus_master_2023 %>%
  filter(novel_flag != "Presumptive Novel") %>% 
  group_by(mosquito_species, sample_number) %>%
  summarise(n = n_distinct(mosquito_species)) %>%
  group_by(mosquito_species) %>%
  summarise(n = n()) %>%
  adorn_totals()



years <- virus_master_2023 %>% 
  filter(novel_flag != "Presumptive Novel") %>%
  group_by(collection_year, virus_name) %>% 
  summarise(count = n_distinct(virus_name, collection_year)) %>% 
  pivot_wider(names_from = collection_year, values_from = count, values_fill = 0)


virus_lineage <- virus_master_2023 %>% 
  filter(novel_flag != "Presumptive Novel") %>%
  group_by(virus_name, viral_family, genome) %>% 
  summarise(n = n_distinct(virus_name, viral_family, genome)) %>% 
  select(-"n")

virus_master_2023 %>% 
  filter(novel_flag != "Presumptive Novel") %>%
  group_by(virus_name) %>% 
  summarise(n_contigs = n(),
            mean_cov = mean(coverage),
            min_cov = min(coverage),
            max_cov = max(coverage),
            mean_pid = mean(greatest_identity_percent),
            med_pid = median(greatest_identity_percent),
            min_pid = min(greatest_identity_percent),
            max_pid = max(greatest_identity_percent),
            total_reads = sum(total_read_count),
            longest_contig = max(contig_length)) %>% 
  left_join(virus_lineage, by = "virus_name") %>% 
  group_by(genome, viral_family) %>% 
  arrange(genome, viral_family) %>% 
  relocate(total_reads, .after = "n_contigs") %>% 
  relocate(longest_contig, .after = "n_contigs") %>% 
  # mutate(across(.cols = "virus_name",
  #               ~case_when(grepl("Manitoba", virus_name) & !grepl("Manitoba virus", virus_name) ~ 
  #                            paste0(., "*"),
  #                          TRUE ~ .))) %>% 
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
    longest_contig = "Longest Contig (nt)",
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
    
  )  %>% 
  gtsave(filename = "Figures/NON-NOVEL_reads_contigs_summary.html")

# NOVEL

library_summary <- virus_master_2023 %>%
  filter(novel_flag == "Presumptive Novel") %>% 
  group_by(mosquito_species, sample_number) %>%
  summarise(n = n_distinct(mosquito_species)) %>%
  group_by(mosquito_species) %>%
  summarise(n = n()) %>%
  adorn_totals()



years <- virus_master_2023 %>% 
  filter(novel_flag == "Presumptive Novel") %>% 
  group_by(collection_year, virus_name) %>% 
  summarise(count = n_distinct(virus_name, collection_year)) %>% 
  pivot_wider(names_from = collection_year, values_from = count, values_fill = 0)


virus_lineage <- virus_master_2023 %>% 
  filter(novel_flag == "Presumptive Novel") %>% 
  group_by(virus_name, viral_family, genome) %>% 
  summarise(n = n_distinct(virus_name, viral_family, genome)) %>% 
  select(-"n")

virus_master_2023 %>% 
  filter(novel_flag == "Presumptive Novel") %>% 
  group_by(virus_name) %>% 
  summarise(n_contigs = n(),
            mean_cov = mean(coverage),
            min_cov = min(coverage),
            max_cov = max(coverage),
            mean_pid = mean(greatest_identity_percent),
            med_pid = median(greatest_identity_percent),
            min_pid = min(greatest_identity_percent),
            max_pid = max(greatest_identity_percent),
            total_reads = sum(total_read_count),
            longest_contig = max(contig_length)) %>% 
  left_join(virus_lineage, by = "virus_name") %>% 
  group_by(genome, viral_family) %>% 
  arrange(genome, viral_family) %>% 
  relocate(total_reads, .after = "n_contigs") %>% 
  relocate(longest_contig, .after = "n_contigs") %>% 
  # mutate(across(.cols = "virus_name",
  #               ~case_when(grepl("Manitoba", virus_name) & !grepl("Manitoba virus", virus_name) ~ 
  #                            paste0(., "*"),
  #                          TRUE ~ .))) %>% 
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
    longest_contig = "Longest Contig (nt)",
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
    
  )  %>% 
  gtsave(filename = "Figures/NOVEL_reads_contigs_summary.html")

# Detected in mosquitoes table ----

# NOT NOVEL

library_summary <- virus_master_2023 %>%
  filter(novel_flag != "Presumptive Novel") %>% 
  group_by(mosquito_species, sample_number) %>%
  summarise(n = n_distinct(mosquito_species)) %>%
  group_by(mosquito_species) %>%
  summarise(n = n()) %>%
  adorn_totals()

years <- virus_master_2023 %>% 
  filter(novel_flag != "Presumptive Novel") %>%
  group_by(collection_year, virus_name) %>% 
  summarise(count = n_distinct(virus_name, collection_year)) %>% 
  pivot_wider(names_from = collection_year, values_from = count, values_fill = 0)


virus_lineage <- virus_master_2023 %>% 
  filter(novel_flag != "Presumptive Novel") %>%
  group_by(virus_name, viral_family, genome) %>% 
  summarise(n = n_distinct(virus_name, viral_family, genome)) %>% 
  select(-"n")

virus_master_2023 %>% 
  filter(novel_flag != "Presumptive Novel") %>% 
  group_by(mosquito_species, virus_name) %>% 
  summarise(count = n_distinct(virus_name, sample_number)) %>% 
  pivot_wider(names_from = mosquito_species, values_from = count, values_fill = 0) %>% 
  mutate("Aedes canadensis \n Total: 1" = (`Aedes canadensis` / 1) * 100,
         "Aedes vexans \n Total: 19" = (`Aedes vexans` / 19) * 100,
         # "Anopheles earlei \n Total: 1" = (`Anopheles earlei` / 1) * 100,
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
   # mutate(across(.cols = "virus_name",
   #               ~case_when(grepl("Manitoba", virus_name) & !grepl("Manitoba virus", virus_name) ~ 
   #                            paste0(., "*"),
   #                          TRUE ~ .))) %>% 
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
    # "Anopheles earlei \n Total: 1" = html("<em>Anopheles earlei</em> <br> Total: 1"),
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
  gtsave(filename = "Figures/NON-NOVEL_virus_detections_by_species_summary.html")

# novel

library_summary <- virus_master_2023 %>%
  filter(novel_flag == "Presumptive Novel") %>% 
  group_by(mosquito_species, sample_number) %>%
  summarise(n = n_distinct(mosquito_species)) %>%
  group_by(mosquito_species) %>%
  summarise(n = n()) %>%
  adorn_totals()

years <- virus_master_2023 %>% 
  filter(novel_flag == "Presumptive Novel") %>% 
  group_by(collection_year, virus_name) %>% 
  summarise(count = n_distinct(virus_name, collection_year)) %>% 
  pivot_wider(names_from = collection_year, values_from = count, values_fill = 0)


virus_lineage <- virus_master_2023 %>% 
  filter(novel_flag == "Presumptive Novel") %>% 
  group_by(virus_name, viral_family, genome) %>% 
  summarise(n = n_distinct(virus_name, viral_family, genome)) %>% 
  select(-"n")

virus_master_2023 %>% 
  filter(novel_flag == "Presumptive Novel") %>% 
  group_by(mosquito_species, virus_name) %>% 
  summarise(count = n_distinct(virus_name, sample_number)) %>% 
  pivot_wider(names_from = mosquito_species, values_from = count, values_fill = 0) %>% 
  mutate("Aedes canadensis \n Total: 1" = (`Aedes canadensis` / 1) * 100,
         "Aedes vexans \n Total: 19" = (`Aedes vexans` / 19) * 100,
         "Anopheles earlei \n Total: 1" = (`Anopheles earlei` / 1) * 100,
         "Coquillettidia perturbans \n Total: 6" = (`Coquillettidia perturbans` / 6) * 100,
         "Culex tarsalis \n Total: 11" = (`Culex tarsalis` / 11) * 100,
         "Ochlerotatus dorsalis \n Total: 5" = (`Ochlerotatus dorsalis` / 5) * 100,
         # "Ochlerotatus flavescens \n Total: 1" = (`Ochlerotatus flavescens` / 1) * 100,
         "Ochlerotatus triseriatus \n Total: 1" = (`Ochlerotatus triseriatus` / 1) * 100,
  ) %>% select(-`Aedes canadensis`:-`Ochlerotatus triseriatus`) %>% 
  left_join(virus_lineage, by = "virus_name") %>% 
  left_join(years, by = "virus_name") %>% 
  relocate("2020", .after = "virus_name") %>% relocate("2021", .after = "2020") %>% 
  mutate(across(.cols = c("2020", "2021"),
                ~case_when(. == "1" ~ TRUE,
                           . == "0" ~ FALSE))) %>% 
  # mutate(across(.cols = "virus_name",
  #               ~case_when(grepl("Manitoba", virus_name) & !grepl("Manitoba virus", virus_name) ~ 
  #                            paste0(., "*"),
  #                          TRUE ~ .))) %>% 
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
    # "Ochlerotatus flavescens \n Total: 1" = html("<em>Ochlerotatus flavescens</em> <br> Total: 1"),
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
  gtsave(filename = "Figures/NOVEL_virus_detections_by_species_summary.html")

# Organism figs ----

virus_master_2023 %>%
  filter(novel_flag == "Presumptive Novel") %>% 
  distinct(virus_name, mosquito_species, .keep_all = TRUE) %>%
  group_by(viral_family, genome, mosquito_species) %>%
  summarise(n = n()) %>%
  ungroup() %>% group_by(mosquito_species) %>%
  summarise(sum = sum(n))

virus_org <- virus_master_2023 %>% 
  distinct(virus_name, mosquito_species, .keep_all = TRUE) %>% 
  group_by(viral_family, genome, mosquito_species) %>% 
  summarise(n = n()) %>% 
  arrange(genome, desc(n)) %>% 
  ggplot(aes(x = fct_inorder(viral_family), y = n, fill = genome)) +
  geom_col() +
  theme_bw(base_size = 14) +
  scale_fill_viridis_d("Genome") +
  scale_y_continuous(name = "Number of Distinct Viruses", breaks = c(0, 2, 4, 6, 8)) + 
  facet_wrap(~mosquito_species, scales = "free_x", labeller = labeller(
    mosquito_species =
      c(
        "Aedes canadensis" = "<strong>*Aedes canadensis* (n = 5 distinct viruses, 1 library)</strong>",
        "Aedes vexans" = "<strong>*Aedes vexans* (n = 28 distinct viruses, 19 libraries)</strong>",
        "Anopheles earlei" = "<strong>*Anopheles earlei* (n = 1 distinct virus, 1 library)</strong>",
        "Coquillettidia perturbans" = "<strong>*Coquillettidia perturbans* (n = 13 distinct viruses, 6 libraries)</strong>",
        "Culex tarsalis" = "<strong>*Culex tarsalis* (n = 34 distinct viruses, 11 libraries)</strong>",
        "Ochlerotatus dorsalis" = "<strong>*Ochlerotatus dorsalis* (n = 15 distinct viruses, 5 libraries)</strong>",
        "Ochlerotatus flavescens" = "<strong>*Ochlerotatus flavescens* (n = 2 distinct viruses, 1 library)</strong>",
        "Ochlerotatus triseriatus" = "<strong>*Ochlerotatus triseriatus* (n = 2 distinct viruses, 1 library)</strong>"
      )
  )) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 14),
        strip.text = ggtext::element_markdown(),
        panel.border = element_rect(fill = NA, color = "gray90"),
        legend.position = c(0.8, 0.1)) +
  labs(x = "") 

tiff(here("virus_organism_plot.tiff"), units="in", width=16, height=10, res=300)
virus_org
dev.off()

# Non novel

virus_org <- virus_master_2023 %>% 
  filter(novel_flag == "Not Novel") %>% 
  distinct(virus_name, mosquito_species, .keep_all = TRUE) %>% 
  group_by(viral_family, genome, mosquito_species) %>% 
  summarise(n = n()) %>% 
  arrange(genome, desc(n)) %>% 
  ggplot(aes(x = fct_inorder(viral_family), y = n, fill = genome)) +
  geom_col() +
  theme_bw(base_size = 14) +
  scale_fill_viridis_d("Genome") +
  scale_y_continuous(name = "Number of Distinct Viruses", breaks = c(0, 2, 4, 6, 8)) + 
  facet_wrap(~mosquito_species, scales = "free_x", labeller = labeller(
    mosquito_species =
      c(
        "Aedes canadensis" = "<strong>*Aedes canadensis* (n = 1 distinct viruses, 1 library)</strong>",
        "Aedes vexans" = "<strong>*Aedes vexans* (n = 23 distinct viruses, 19 libraries)</strong>",
        "Coquillettidia perturbans" = "<strong>*Coquillettidia perturbans* (n = 12 distinct viruses, 6 libraries)</strong>",
        "Culex tarsalis" = "<strong>*Culex tarsalis* (n = 31 distinct viruses, 11 libraries)</strong>",
        "Ochlerotatus dorsalis" = "<strong>*Ochlerotatus dorsalis* (n = 11 distinct viruses, 5 libraries)</strong>",
        "Ochlerotatus flavescens" = "<strong>*Ochlerotatus flavescens* (n = 2 distinct viruses, 1 library)</strong>",
        "Ochlerotatus triseriatus" = "<strong>*Ochlerotatus triseriatus* (n = 1 distinct viruses, 1 library)</strong>"
      )
  )) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 14),
        strip.text = ggtext::element_markdown(),
        panel.border = element_rect(fill = NA, color = "gray90"),
        legend.position = c(0.8, 0.1)) +
  labs(x = "") 

tiff(here("NON-NOVEL_virus_organism_plot.tiff"), units="in", width=16, height=10, res=300)
virus_org
dev.off()

# Novel

virus_org <- virus_master_2023 %>% 
  filter(novel_flag == "Presumptive Novel") %>% 
  distinct(virus_name, mosquito_species, .keep_all = TRUE) %>% 
  group_by(viral_family, genome, mosquito_species) %>% 
  summarise(n = n()) %>% 
  arrange(genome, desc(n)) %>% 
  ggplot(aes(x = fct_inorder(viral_family), y = n, fill = genome)) +
  geom_col() +
  theme_bw(base_size = 14) +
  scale_fill_viridis_d("Genome") +
  scale_y_continuous(name = "Number of Distinct Viruses", breaks = c(0, 2, 4, 6, 8)) + 
  facet_wrap(~mosquito_species, scales = "free_x", labeller = labeller(
    mosquito_species =
      c(
        "Aedes canadensis" = "<strong>*Aedes canadensis* (n = 4 distinct viruses, 1 library)</strong>",
        "Aedes vexans" = "<strong>*Aedes vexans* (n = 6 distinct viruses, 19 libraries)</strong>",
        "Anopheles earlei" = "<strong>*Anopheles earlei* (n = 1 distinct virus, 1 library)</strong>",
        "Coquillettidia perturbans" = "<strong>*Coquillettidia perturbans* (n = 1 distinct viruses, 6 libraries)</strong>",
        "Culex tarsalis" = "<strong>*Culex tarsalis* (n = 4 distinct viruses, 11 libraries)</strong>",
        "Ochlerotatus dorsalis" = "<strong>*Ochlerotatus dorsalis* (n = 4 distinct viruses, 5 libraries)</strong>",
        "Ochlerotatus triseriatus" = "<strong>*Ochlerotatus triseriatus* (n = 1 distinct viruses, 1 library)</strong>"
      )
  )) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 14),
        strip.text = ggtext::element_markdown(),
        panel.border = element_rect(fill = NA, color = "gray90"),
        legend.position = c(0.8, 0.1)) +
  labs(x = "") 

tiff(here("NOVEL_virus_organism_plot.tiff"), units="in", width=16, height=10, res=300)
virus_org
dev.off()

virus_master_2023 %>%
  group_by(viral_family, genome, mosquito_species) %>% 
  summarise(sum_reads = sum(total_read_count)) %>% 
  ungroup() %>% group_by(mosquito_species) %>%
  summarise(sum = sum(sum_reads))

virus_reads <- virus_master_2023 %>% 
  group_by(viral_family, genome, mosquito_species) %>% 
  summarise(sum_reads = sum(total_read_count)) %>% 
  arrange(genome, desc(sum_reads)) %>% 
  mutate(log_reads = log10(sum_reads)) %>% 
  ggplot(aes(x = fct_inorder(viral_family), y = log_reads, fill = genome)) +
  geom_col() +
  theme_bw(base_size = 14) +
  scale_fill_viridis_d("Genome") +
  scale_y_continuous(name = "Total Viral Reads (Log 10)") + 
  facet_wrap(~mosquito_species, scales = "free_x", labeller = labeller(
    mosquito_species =
      c(
        "Aedes canadensis" = "<strong>*Aedes canadensis* (n = 12,031 viral reads, 1 library)</strong>",
        "Aedes vexans" = "<strong>*Aedes vexans* (n = 1,212,834 viral reads, 19 libraries)</strong>",
        "Anopheles earlei" = "<strong>*Anopheles earlei* (n = 391 viral reads, 1 library)</strong>",
        "Coquillettidia perturbans" = "<strong>*Coquillettidia perturbans* (n = 167,718 viral reads, 6 libraries)</strong>",
        "Culex tarsalis" = "<strong>*Culex tarsalis* (n = 945,820 viral reads, 11 libraries)</strong>",
        "Ochlerotatus dorsalis" = "<strong>*Ochlerotatus dorsalis* (n = 1,718,651 viral reads, 5 libraries)</strong>",
        "Ochlerotatus flavescens" = "<strong>*Ochlerotatus flavescens* (n = 1,576,622 viral reads, 1 library)</strong>",
        "Ochlerotatus triseriatus" = "<strong>*Ochlerotatus triseriatus* (n = 3,019 viral reads, 1 library)</strong>"
      )
  )) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 14),
        strip.text = ggtext::element_markdown(),
        panel.border = element_rect(fill = NA, color = "gray90"),
        legend.position = c(0.8, 0.1)) +
  labs(x = "") 

tiff(here("virus_reads_plot.tiff"), units="in", width=16, height=10, res=300)
virus_reads
dev.off()


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

total_family <- virus_master_2023 %>% 
  distinct(virus_name, viral_family, .keep_all = TRUE) %>% 
  group_by(genome, viral_family) %>% 
  summarise(n = n()) %>% 
  arrange(genome, n) %>%
  ggplot(aes(x = fct_inorder(viral_family), y = n)) +
  geom_col(aes(fill = genome)) +
  scale_fill_viridis_d("Genome") +
  scale_y_continuous(limits = c(0, 14), breaks = c(0, 2, 4, 6, 8, 10, 12, 14), expand = c(0,0)) +
  coord_flip() + 
  theme_bw() +
  labs(x = "Viral Family",
       y = "Number of Distinct Viruses Detected")

tiff(here("all_viruses_family_plot.tiff"), units="in", width=16, height=10, res=300)
total_family
dev.off()

non_novel_family <- virus_master_2023 %>% 
  filter(novel_flag == "Not Novel") %>% 
  distinct(virus_name, viral_family, .keep_all = TRUE) %>% 
  group_by(genome, viral_family) %>% 
  summarise(n = n()) %>% 
  arrange(genome, n) %>%
  ggplot(aes(x = fct_inorder(viral_family), y = n)) +
  geom_col(aes(fill = genome)) +
  scale_fill_viridis_d("Genome") +
  scale_y_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12), expand = c(0,0)) +
  coord_flip() + 
  theme_bw() +
  labs(x = "Viral Family",
       y = "Number of Distinct Viruses Detected")

tiff(here("non_novel_family_plot.tiff"), units="in", width=16, height=10, res=300)
non_novel_family
dev.off()

novel_family <- virus_master_2023 %>% 
  filter(novel_flag != "Not Novel") %>% 
  distinct(virus_name, viral_family, .keep_all = TRUE) %>% 
  group_by(genome, viral_family) %>% 
  summarise(n = n()) %>% 
  arrange(genome, n) %>%
  ggplot(aes(x = fct_inorder(viral_family), y = n)) +
  geom_col(aes(fill = genome)) +
  scale_fill_viridis_d("Genome") +
  scale_y_continuous(limits = c(0, 4), breaks = c(0, 1, 2, 3, 4), expand = c(0,0)) +
  coord_flip() + 
  theme_bw() +
  labs(x = "Viral Family",
       y = "Number of Distinct Viruses Detected")

tiff(here("novel_family_plot.tiff"), units="in", width=16, height=10, res=300)
novel_family
dev.off()

 

virus_master_2023 %>% 
  distinct(virus_name, viral_family, .keep_all = TRUE) %>% 
  group_by(genome, viral_family) %>% 
  summarise(n = n()) %>% 
  arrange(genome, desc(n)) %>% 
  adorn_totals() %>% 
  gt() %>% 
  cols_label(
    genome = "Genome",
    viral_family = "Family")

# Novel viruses ----


virus_master_2023 %>% 
  distinct(virus_name, viral_family, .keep_all = TRUE) %>% 
  filter(str_detect(virus_name, "Manitoba ")) %>% 
  filter(!str_detect(virus_name, "Manitoba virus")) %>% 
  group_by(genome, viral_family) %>% 
  summarise(n = n())

# VENN DIAGRAM ----

venn <- virus_master_2023 %>% 
  separate_wider_delim(cols = mosquito_species, names = c("genus", "species"), delim = " ") %>% 
  distinct(virus_name, genus, .keep_all = TRUE) %>% 
  group_by(virus_name, genus) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  complete(virus_name, genus, fill = list(n = 0)) %>% 
  pivot_wider(names_from = "genus", values_from = "n") %>% 
  mutate(across(c("Aedes", "Anopheles", "Coquillettidia", "Culex", "Ochlerotatus"),
                ~ as.character(.))) %>%

  mutate(across(c("Aedes", "Anopheles", "Coquillettidia", "Culex", "Ochlerotatus"), 
                ~ case_when(. == "1" ~ virus_name,
                            TRUE ~ NA_character_))) %>% 
  select(-1)
  

aedes_venn <- venn %>%
  select("Aedes") %>%
  drop_na()

culex_venn <- venn %>%
  select("Culex") %>%
  drop_na()

coquillettidia_venn <- venn %>%
  select("Coquillettidia") %>%
  drop_na()

anopheles_venn <- venn %>%
  select("Anopheles") %>%
  drop_na()

ochlerotatus_venn <- venn %>%
  select("Ochlerotatus") %>%
  drop_na()

venn_map <- map(c(aedes_venn, culex_venn, coquillettidia_venn, anopheles_venn, ochlerotatus_venn), unique)


venn_data <- Venn(venn_map)

data <- process_data(venn_data) 
  

venn_diagram <- ggplot() +
  geom_sf(aes(fill = count), data = venn_region(data)) +
  # 2. set edge layer
  geom_sf(aes(color = id), data = venn_setedge(data), show.legend = FALSE) +
  # 3. set label layer
  geom_sf_text(aes(label = name), fontface = "bold.italic", data = venn_setlabel(data)) +
  # 4. region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void() +
  scale_fill_viridis_c("Number of Viruses")


tiff(here("virus_venn_diagram.tiff"), units="in", width=16, height=10, res=300)
venn_diagram
dev.off()

# Non novel

venn <- virus_master_2023 %>% 
  filter(novel_flag == "Not Novel") %>% 
  separate_wider_delim(cols = mosquito_species, names = c("genus", "species"), delim = " ") %>% 
  distinct(virus_name, genus, .keep_all = TRUE) %>% 
  group_by(virus_name, genus) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  complete(virus_name, genus, fill = list(n = 0)) %>% 
  pivot_wider(names_from = "genus", values_from = "n") %>% 
  mutate(across(c("Aedes", "Coquillettidia", "Culex", "Ochlerotatus"),
                ~ as.character(.))) %>%
  
  mutate(across(c("Aedes","Coquillettidia", "Culex", "Ochlerotatus"), 
                ~ case_when(. == "1" ~ virus_name,
                            TRUE ~ NA_character_))) %>% 
  select(-1)


aedes_venn <- venn %>%
  select("Aedes") %>%
  drop_na()

culex_venn <- venn %>%
  select("Culex") %>%
  drop_na()

coquillettidia_venn <- venn %>%
  select("Coquillettidia") %>%
  drop_na()

# anopheles_venn <- venn %>%
#   select("Anopheles") %>%
#   drop_na()

ochlerotatus_venn <- venn %>%
  select("Ochlerotatus") %>%
  drop_na()

venn_map <- map(c(aedes_venn, culex_venn, coquillettidia_venn, ochlerotatus_venn), unique)


venn_data <- Venn(venn_map)

data <- process_data(venn_data) 


venn_diagram <- ggplot() +
  geom_sf(aes(fill = count), data = venn_region(data)) +
  # 2. set edge layer
  geom_sf(aes(color = id), data = venn_setedge(data), show.legend = FALSE) +
  # 3. set label layer
  geom_sf_text(aes(label = name), fontface = "bold.italic", data = venn_setlabel(data)) +
  # 4. region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void() +
  scale_fill_viridis_c("Number of Viruses")


tiff(here("NON-NOVEL_virus_venn_diagram.tiff"), units="in", width=16, height=10, res=300)
venn_diagram
dev.off()

# NOVEL venn 

venn <- virus_master_2023 %>% 
  filter(novel_flag == "Presumptive Novel") %>% 
  separate_wider_delim(cols = mosquito_species, names = c("genus", "species"), delim = " ") %>% 
  distinct(virus_name, genus, .keep_all = TRUE) %>% 
  group_by(virus_name, genus) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  complete(virus_name, genus, fill = list(n = 0)) %>% 
  pivot_wider(names_from = "genus", values_from = "n") %>% 
  mutate(across(c("Aedes", "Anopheles", "Coquillettidia", "Culex", "Ochlerotatus"),
                ~ as.character(.))) %>%
  
  mutate(across(c("Aedes", "Anopheles", "Coquillettidia", "Culex", "Ochlerotatus"), 
                ~ case_when(. == "1" ~ virus_name,
                            TRUE ~ NA_character_))) %>% 
  select(-1)


aedes_venn <- venn %>%
  select("Aedes") %>%
  drop_na()

culex_venn <- venn %>%
  select("Culex") %>%
  drop_na()

coquillettidia_venn <- venn %>%
  select("Coquillettidia") %>%
  drop_na()

anopheles_venn <- venn %>%
  select("Anopheles") %>%
  drop_na()

ochlerotatus_venn <- venn %>%
  select("Ochlerotatus") %>%
  drop_na()

venn_map <- map(c(aedes_venn, culex_venn, coquillettidia_venn, ochlerotatus_venn), unique)


venn_data <- Venn(venn_map)

data <- process_data(venn_data) 


venn_diagram <- ggplot() +
  geom_sf(aes(fill = count), data = venn_region(data)) +
  # 2. set edge layer
  geom_sf(aes(color = id), data = venn_setedge(data), show.legend = FALSE) +
  # 3. set label layer
  geom_sf_text(aes(label = name), fontface = "bold.italic", data = venn_setlabel(data)) +
  # 4. region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void() +
  scale_fill_viridis_c("Number of Viruses")


tiff(here("NOVEL_virus_venn_diagram.tiff"), units="in", width=16, height=10, res=300)
venn_diagram
dev.off()





