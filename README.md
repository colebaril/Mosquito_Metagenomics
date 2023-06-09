# Mosquito Metagenomics

Code and workflow for analyzing, summarizing, and tabulating mosquito metagenomics CLC Genomics Workbench BLAST results files.

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
