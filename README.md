# Mosquito Metagenomics

Code and workflow for analyzing, summarizing, visualizing and tabulating mosquito metagenomics CLC Genomics Workbench BLAST results files.

# Tables

[Click here to view the HTML table](Tables/known_virus_detections_by_species_summary.html)








# Figures


<figure>
  <img src="https://github.com/colebaril/Mosquito_Metagenomics/assets/110275137/f1762075-8f07-425b-b6d3-61f5aabb91a6" alt="Combined_Virus_Venn_Diagram_Mar2024">
  <figcaption>Figure 1: Venn Diagram depicting the number of unique viruses identified in each mosquito genus separated by known viruses (A) which have been previously reported and novel viruses (B) which we newly identified.</figcaption>
</figure>

<details>
  <summary>Click to view code</summary>

```r
require(pacman)
p_load(tidyverse, janitor, here, gt, forcats, tiff, openxlsx, ggVennDiagram, ggpattern, svglite)

virus_master_2023 <- read_csv(here("Data/tblastx_master.csv"))

# KNOWN VIRUSES

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

  mutate(across(c("Aedes", "Coquillettidia", "Culex", "Ochlerotatus"), 
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

data_range_known <- venn_regionedge(process_data(venn_data))


venn_diagram_known <- ggplot() +
  # 1. region count layer
  geom_polygon(aes(X, Y, fill = count, group = id), 
               data = venn_regionedge(data)) +
  # 2. set edge layer
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(data), 
            show.legend = FALSE) +
  # 3. set label layer
  geom_text(aes(X, Y, label = name), 
            data = venn_setlabel(data),
            fontface = "bold.italic",
            size = 5,
            nudge_x = -0.01) +
  # 4. region label layer
  geom_label(aes(X, Y, label = count), 
             data = venn_regionlabel(data),
             fontface = "bold",
             size = 7) +
  coord_equal() +
  theme_void() +
  theme(
    # legend.title = element_text(size = 15),
    # legend.text = element_text(size = 13),
    plot.title = element_text(size = 20, face = "bold")
  ) +
  # scale_fill_viridis_c("Number of Viruses") +
  labs(title = "A) Known Viruses")


tiff(here("virus_venn_diagram.tiff"), units="in", width=18, height=10, res=300)
venn_diagram
dev.off()


#  NOVEL VIRUSES ----


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


aedes_venn_novel <- venn %>%
  select("Aedes") %>%
  drop_na()

culex_venn_novel <- venn %>%
  select("Culex") %>%
  drop_na()

coquillettidia_venn_novel <- venn %>%
  select("Coquillettidia") %>%
  drop_na()

anopheles_venn_novel <- venn %>%
  select("Anopheles") %>%
  drop_na()

ochlerotatus_venn_novel <- venn %>%
  select("Ochlerotatus") %>%
  drop_na()

venn_map <- map(c(aedes_venn_novel, culex_venn_novel, anopheles_venn_novel, 
                  coquillettidia_venn_novel, ochlerotatus_venn_novel), unique)


venn_data <- Venn(venn_map)

data <- process_data(venn_data)

data_range_novel <- venn_regionedge(process_data(venn_data))


venn_diagram_novel <- ggplot() +
  # 1. region count layer
  geom_polygon(aes(X, Y, fill = count, group = id), 
               data = venn_regionedge(data)) +
  # 2. set edge layer
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(data), 
            show.legend = FALSE) +
  # 3. set label layer
  geom_text(aes(X, Y, label = name), 
            data = venn_setlabel(data),
            fontface = "bold.italic",
            size = 5) +
  # 4. region label layer
  geom_label(aes(X, Y, label = count), 
             data = venn_regionlabel(data),
             fontface = "bold",
             size = 7) +
  coord_equal() +
  theme_void() +
  theme(
    # legend.title = element_text(size = 15),
    # legend.text = element_text(size = 13),
    plot.title = element_text(size = 20, face = "bold")
  ) +
  # scale_fill_viridis_c("Number of Viruses") +
  labs(title = "B) Novel Viruses")

# Combine venn diagrams 

library(patchwork)

(venn_diagram_known + venn_diagram_novel) +
  plot_layout(guides = "collect") &
  scale_fill_viridis_c("Number of Viruses", limits = range(c(data_range_novel$count, data_range_known$count))) &
  theme(
    legend.title = element_text(size = 15, vjust = 0.75),
    legend.text = element_text(size = 13),
    legend.position = "bottom"
  ) 

ggsave("Combined_Virus_Venn_Diagram_Mar2024.png", plot = last_plot(), width=18, height=10)
```
</details>

<figure>
  <img src="https://github.com/colebaril/Mosquito_Metagenomics/assets/110275137/c11fec09-5867-4147-8fa8-25070553b89b" alt="Viruses Plot">
  <figcaption>Figure 2: Faceted bar plot displaying the number of unique viruses identified (y-axis) in each mosquito species. The total number of viruses identified in a mosquito species and the number of sequencing libraries representing each mosquito species is displayed at the top of each facet. Viral family is shown on the x-axis and bars are coloured based on genome type of the family. The number of novel viruses is indicated by hash marks.</figcaption>
</figure>
<details>
  <summary>Click to view code</summary>

```r
require(pacman)
p_load(tidyverse, janitor, here, gt, forcats, tiff, openxlsx, ggVennDiagram, ggpattern, svglite)

virus_master_2023 <- read_csv(here("Data/tblastx_master.csv"))

virus_master_2023 %>%
  filter(novel_flag == "Presumptive Novel") %>% 
  distinct(virus_name, mosquito_species, .keep_all = TRUE) %>%
  group_by(viral_family, genome, mosquito_species) %>%
  summarise(n = n()) %>%
  ungroup() %>% group_by(mosquito_species) %>%
  summarise(sum = sum(n))

p_load(ggpattern, svglite)

virus_org <- virus_master_2023 %>% 
  distinct(virus_name, mosquito_species, .keep_all = TRUE) %>% 
  group_by(viral_family, genome, mosquito_species, novel_flag) %>% 
  summarise(n = n()) %>% 
  arrange(genome, desc(n)) %>% 
  ggplot(aes(x = fct_inorder(viral_family), y = n, fill = genome, pattern = novel_flag)) +
  geom_col_pattern(pattern_fill = "white",
                   pattern_fill2 = "white",
                   pattern_colour = "white",
                   colour = "white",
                   pattern_size = 0.5,
                   pattern_alpha = 0.5,
                   show.legend = TRUE,
                   pattern_key_scale_factor = .5) +
  theme_bw(base_size = 14) +
  scale_fill_viridis_d("Genome", guide = guide_legend(override.aes = list(pattern = "none"))) +
  scale_pattern_manual("", values = c("Not Novel" = "none", "Presumptive Novel" = "stripe"), 
                       labels = c("Known Virus", "Novel Virus")) +
  scale_y_continuous(name = "Number of Distinct Viruses", breaks = c(0, 2, 4, 6, 8)) + 
  facet_wrap(~ mosquito_species, ncol = 2, scales = "free_x", labeller = labeller(
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
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        strip.text = ggtext::element_markdown(size = 18),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18),
        legend.position = "bottom",
        legend.key.size = unit(1.5, 'cm')
        ) +
        # panel.border = element_rect(fill = NA, color = "gray90")) +
        # legend.position = c(0.8, 0.1)) +
  labs(x = "")

tiff(here("virus_organism_plot.tiff"), units="in", width=16, height=16, res=300)
virus_org
dev.off()

ggsave(here("virus_organism_plot.png"), plot = virus_org, dpi = 300, units = "in", width = 16, height = 16)
```
</details>




# Links to Files

Some files generated are too large to store on GitHub (e.g., contig sequences, publication-quality images). Large files can be accessed from DropBox folders below:

[Dropbox link to raw assemblies separated by sample](https://www.dropbox.com/s/m194auk7oxlpnwa/Contigs.zip?dl=0)

[Dropbox link to raw assemblies concatenated into a .csv](https://www.dropbox.com/s/ppe83l8lbvdz262/contigs_raw.zip?dl=0)

[Dropbox link to publication quality figures](https://www.dropbox.com/sh/2lw74xdap2ofxtb/AAACTMIU-ipvd0GgKVGxuzjma?dl=0)

Raw RNA sequencing reads can be retrieved from the NCBI short sequence read archive under the SRA accession number PRJNA793247.


# Re-BLAST 2023

Due to the significant changes to GenBank we are re-BLASTing results. To speed things up, I am curating custom BLAST databases using sequences downloaded from GenBank's Nucleotide database using Entrez Queries to eliminate sequences that are extremely unlikely to be present in our samples that make up a large portion of the Nucleotide database. BLAST is being run locally using CLC Genomics Workbench. I am BLASTing by taxon categories (e.g., viruses, fungi, parasites) as we are limited by computational availability. For example, eliminating the top few viruses (unlikely to be found in mosquitoes) in the nucleotide database reduces the number of sequences by 90%. 

## Viruses

**Entrez Query:**

Viruses [ORGN] NOT Coronavirus [ORGN] NOT Human immunodeficiency virus 1 [ORGN] NOT Influenza A virus [ORGN] NOT Hepacivirus C [ORGN] NOT Hepatitis B virus [ORGN] NOT Influenza B virus [ORGN] NOT Rotavirus A [ORGN] NOT Norwalk virus [ORGN] NOT Simian immunodeficiency virus [ORGN] 

## Re-BLAST 2023 Steps

1. Concatenate all read mapping data and contigs: `01 - Mappings & Contigs.R`
2. Re-BLAST all contig sequences against local NCBI virus blast (curated via above Entrez Query) using CLC
3. Read local BLAST results and filter results: `02 - Read Nt Local BLAST.R`
4. BLAST all contig sequences passing filters at NCBI against all organisms to eliminate false positives using CLC
5. Read BLAST at NCBI results and filter out non-virus sequences: `03 - BLAST at NCBI.R`
6. tBLASTx on contig sequences passing filters from step using local database on CLC to determine aa identites
7. Read tBLASTx results and filter: `04 - tBLASTx.R`
8. Scripts 05 - 09 are to produce tables, figures and summaries for publications and supervisors

# Non-Virus Sequences

Non-virus sequences (e.g., Fungi, parasites, bacteria, plants, vertebrates) were considerably lower in overall number as well as quality. Furthermore, for these reasons, it was difficult to discern species from the recovered sequences. Some of the reasons for this is that compared to viruses where often we recovered the near complete genome, for fungi, parasites, bacteria, plants and vertebrates we recovered mostly rRNA, mitochondrial sequences, or in the case of plants, chloroplast sequences. Therefore, we decided to conduct analyses for non-virus sequences recovered at a higher level (e.g., Family, Genus) rather than species. This still gives us a good idea about what organisms are harboured by mosquitoes. 

# Methods

## Host and Quality Filtering

Chan Zuckerberg ID Metagenomic Pipeline v6.8 (Chan Zuckerberg Biohub; CZID), an open-sourced cloud-based bioinformatics platform (https://czid.org/) was used for quality control and host filtration of reads as well as de novo assembly and taxonomic binning as described by Batson et al., (2021) and Kalantar et al., (2020). The CZID pipeline employs STAR and Bowtie2 to perform host filtration (human and mosquito), Trimmomatic for adapter trimming, Price Seq for removal of low-quality reads, LZW for the removal of low complexity reads and CZIDdedup for duplicate read identification.

## *De Novo* Assembly

The host and quality filtered reads were allowed to continue through the CZID pipeline, which involves de novo assembly with SPADES using default settings and only the assembly module. After assembly, reads are mapped back to contigs with Bowtie2. The host and quality filtered reads from CZID (the Bowtie2 output) were downloaded and assembled with the CLC Genomics Workbench version 20 assembler with a minimum contig length of 250 nt, mismatch cost of 2, insertion cost of 3, deletion cost of 3, length fraction of 0.7 and a similarity fraction of 0.95. Contigs were subject to BLASTn and BLASTp searches on the NCBI nt and nr databases, respectively. The BLAST results were very similar between CZID and CLC, thus we opted to use CLC Genomics Workbench version 20 for subsequent analyses.

## BLAST

Assembled contigs were subject to BLASTn searches on the NCBI non-redundant nucleotide and protein databases, respectively and contigs were assigned to taxa based on BLAST results. Positive contigs from BLASTn were subject to tBLASTx to identify amino acid identities. We were looking for non-mosquito sequences of viral, bacterial, parasitic, fungal and plant origin and discarded sequences that were of mosquito origin.


To begin, BLASTn search results were filtered by E-value (≤1x 10<sup>-100</sup>) and contig length (≥250). Contigs with a match length of ≥250 nt, ≥90% for both nt and aa sequence similarity were classified as hits. Coverage depth was analyzed on a per-sequence basis: for viruses, we used a minimum threshold of 10X coverage depth, and we were less strict for protozoan, fungal, bacterial, plant and chordate coverage. For virus hits, contigs meeting these criteria were subject to tBLASTx searches to identify amino acid identities.

Contigs with a percent nt identity <90% were flagged as potentially novel. While the ICTV sets specific standards for different viral taxa for percent identity to claim a novel virus, many of the viruses we recovered are unclassified beyond the order or family classification, which can make selecting an identity threshold difficult. Therefore, we selected ≤90% as the cut-off because Kalantar et al., (2020) suggests that <90% nt identity is a good general threshold.
