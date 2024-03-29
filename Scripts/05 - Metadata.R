require(pacman)
pacman::p_load(tidyverse, janitor, here, ggtext)

# LOAD AND CLEAN ----
metadata <- read_csv(here("Raw/Mosquito sequencing metadata.csv")) %>% 
  clean_names() %>% 
  separate(sample_id, 
           into = c("id", "id_"),
           sep = "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])") %>% select(-"id_") %>% 
  mutate(species = case_when(species == "Cx. Tarsalis" ~ "Cx. tarsalis",
                             species == "Cq. Perturbans" ~ "Cq. perturbans",
                             TRUE ~ as_factor(species))) %>% 
  # Converting reads to log form - easier to interpret. 
  mutate(total_reads_log = log10(total_reads)) %>% 
  relocate(total_reads_log, .after = total_reads) %>% 
  mutate(number_reads_passing_qc_log = log10(number_reads_passing_qc)) %>% 
  relocate(number_reads_passing_qc_log, .after = number_reads_passing_qc) %>% 
  mutate(number_reads_passing_host_filters_log = log10(number_reads_passing_host_filters)) %>% 
  relocate(number_reads_passing_host_filters_log, .after = number_reads_passing_host_filters)

# Save data

write.csv(metadata, here("Data/metadata_cleaned.csv"))

# PLOT ----
metadata_plot <- metadata %>% 
  ggplot(aes(x = total_reads_log, y = reorder(id, total_reads_log))) + 
  geom_segment(aes(x = number_reads_passing_host_filters_log, xend = total_reads_log, yend = id),
               col = "black", linewidth = 1) +
  geom_point(aes(colour = species, shape = year, size = 1.2)) +
  geom_point(aes(x = number_reads_passing_qc_log, colour = species, shape = year, size = 1.2)) +
  geom_point(aes(x = number_reads_passing_host_filters_log, colour = species, shape = year, size = 1.2)) +
  
  geom_text(aes(label = number_in_pool), hjust = 1.5) +
  theme_bw() +
  scale_x_reverse() +
  scale_colour_viridis_d("Species", labels = c("*Ae. canadensis*", "*Ae. vexans*", "*An. earlei*", "*Cq. perturbans*", 
                                               "*Cx. tarsalis*", "*Oc. dorsalis*", "*Oc. flavescens*", "*Oc. triseriatus*")) +
  scale_shape_manual("Year", values = c(15, 17, 18, 19)) +
  theme(panel.grid.major.x  = element_blank(),
        panel.grid.minor.x  = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_markdown(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.95, 0.3),
        legend.background = element_rect(linetype = 2, linewidth = 0.5, colour = 1)) +
  guides(colour = guide_legend(override.aes = list(size = 5)),
         shape = guide_legend(override.aes = list(size = 5)),
         size = "none") +
  labs(x = "Total Reads (Log 10)",
       y = "")



# SAVE PLOT ----

# Saving as a .tiff at 300 dpi. 

tiff(here("metadata_plot.tiff"), units="in", width=16, height=10, res=300)
metadata_plot
dev.off()

# PLOT ----
host_filters_plot <- metadata %>% 
  ggplot(aes(x = number_reads_passing_qc_log, y = reorder(id, number_reads_passing_qc_log))) + 
  geom_segment(aes(x = number_reads_passing_host_filters_log, xend = number_reads_passing_qc_log, yend = id),
               col = "black", linewidth = 1) +
  geom_point(aes(colour = species, shape = year, size = 1.2)) +
  geom_point(aes(x = number_reads_passing_host_filters_log, colour = species, shape = year, size = 1.2)) +
  
  geom_text(aes(label = number_in_pool), hjust = -0.75) +
  theme_bw() +
  scale_colour_viridis_d("Species", labels = c("*Ae. canadensis*", "*Ae. vexans*", "*An. earlei*", "*Cq. perturbans*", 
                                               "*Cx. tarsalis*", "*Oc. dorsalis*", "*Oc. flavescens*", "*Oc. triseriatus*")) +
  scale_shape_manual("Year", values = c(15, 17, 18, 19)) +
  theme(panel.grid.major.x  = element_blank(),
        panel.grid.minor.x  = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_markdown(),
        legend.position = c(0.95, 0.3),
        legend.background = element_rect(linetype = 2, linewidth = 0.5, colour = 1)) +
  guides(colour = guide_legend(override.aes = list(size = 5)),
         shape = guide_legend(override.aes = list(size = 5)),
         size = "none") +
  labs(x = "Total Reads (Log 10)",
       y = "Sample ID",
       title = "Reads remaining after host filtering (Mosquito)",
       subtitle = "Number of specimens represented by each sample displayed to the right")


# SAVE PLOT ----

# Saving as a .tiff at 300 dpi. 

tiff(here("Figures/host_filters_plot.tiff"), units="in", width=16, height=10, res=300)
host_filters_plot
dev.off()