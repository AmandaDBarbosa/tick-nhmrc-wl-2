pkgs <- c("qiime2R", "phyloseq", "tidyverse", "ampvis2", "ampvis2extras", 
          "ggpubr", "agricolae", "plotly", "viridis", "cowplot", "MicrobeR", 
          "microbiome", "reshape", "decontam", "data.table", "ape", "DESeq2", 
          "vegan", "microbiomeutilities", "knitr", "tibble", "dplyr", 
          "patchwork", "Biostrings")
lapply(pkgs, require, character.only = TRUE)

install.packages("devtools")
devtools::install_github("microsud/microbiomeutilities")

library(readr)
library(tibble)

#Importing ASV Count Data

dada2_ASVs <- read_tsv("~/Desktop/feature-table.tsv", skip = 1)
dada2_ASVs_lab = column_to_rownames(dada2_ASVs, var = "#OTU ID")

#Importing Taxonomy Data

taxonomy <- read_tsv("~/Desktop/taxonomy.tsv")
taxonomy_lab = column_to_rownames(taxonomy, var = "Feature ID")

#Add ASV sequences by reading in FASTA file
#Conversion required from qza to fasta (DONE FOR LOD)
#mkdir Rep-seqs
#qiime tools export \
#--input-path rep_seq.qza \
#--output-path Rep-seqs/

# read sequence file
rep.seqs <- Biostrings::readDNAStringSet("~/Desktop/skin-dna-sequences.fasta", format = "fasta")

# important sequence metadata- *Important: make sure sample names and matched AGRF-ID are listed in the spreadsheet in the same order as "dada2_ASVs". 
library(readr)
metadata <- read_tsv("~/Desktop/metadata-skin.tsv")


#metadata <- read_csv("~/Desktop/NGS-Bioinformatics-Jan2023/Limit of Detection/sample-metadata-new-decont-corrected.csv",
#col_types = cols(
#SampleType = col_factor(levels =c("categorical","Sample", "Control", "Zymo")),
#SampleCategory = col_factor(levels = c("categorical", "Blood", 
#"EBlank", "PC","NC")),
#gDNAID = col_factor(levels = c("categorical","pw","mtb","ctp1795-1","ctp1795-2","mdb68","k372","hc7","ra",
#"pw_2","mtb_2","ctp1795-1_2","ctp1795-2_2","mdb68_2","k372_2","hc7_2","ra_2",
#"pw_3","mtb_3","ctp1795-1_3","ctp1795-2_3","mdb68_3","k372_3","hc7_3","ra_3",
#"co1","co2","co3","co4","co5","co6","co7","co8",
#"co1_2","co2_2","co3_2","co4_2","co5_2","co6_2","co7_2","co8_2",
#"co1_3","co2_3","co3_3","co4_3","co5_3","co6_3","co7_3","co8_3","co","co_2","co_3",
#"ad5","ad5_2","ad5_3","ajm","ajm_2","ajm_3","donor","donor_2","donor_3",
#"EBC1","EBC2","EBC3","EBC4","EBC5","PC","NTC"))))

library(phyloseq)

# Remove second row which contains column data for QIIME2 format
metadata <- metadata[-c(1), ]
metadata_lab = column_to_rownames(metadata, var = "lab-id")
sampledata = sample_data(data.frame(metadata_lab))
sampledata

#To read %<% command below
library(dplyr)

#Renaming sample-ids for consistency and matching to ENA accession numbers.
OldName <- colnames(dada2_ASVs_lab)
NewName <- rownames(metadata_lab)
dada2_ASVs_lab_renamed <- dada2_ASVs_lab %>% rename_at(all_of(OldName), ~ NewName)

# ***Check if AGRF IDs were replaced by matched sample Id's from metadata table.

# Make OTU matrix
otumat <- as.matrix(dada2_ASVs_lab_renamed)
# Make Tax matrix
taxmat <- as.matrix(taxonomy_lab)

#Create phyloseq object
class(otumat)
class(taxmat)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX)
physeq
sample_names(physeq)

#Merge sequence and metadata to create a final phyloseq object
#OLD:
#ps_raw_bact = merge_phyloseq(physeq, metadata_lab, rep.seqs)
#NEW:
phyloseq_object_all = phyloseq(otu_table(OTU), tax_table(TAX), sample_data(metadata_lab))

#Visualize tables using Nice.Table

Nice.Table(ps_raw_bact@sam_data)

#Decontamin

#Install Bioconductor package for R >4.2

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("decontam")

library(decontam)

#Make plot of library size of Samples vs Controls

library(ggplot2)

df <- as.data.frame(sample_data(phyloseq_object_all)) # Put sample_data into a ggplot-friendly data.frame

asv_table <- otu_table(phyloseq_object_all)
sample_data <- data.frame(sample_data(phyloseq_object_all))

df$LibrarySize <- sample_sums(phyloseq_object_all)

df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
libQC <- ggplot(data=df, aes(x=Index, y=LibrarySize, color="SampleType")) + geom_point()
ggsave("libQC.pdf", plot = libQC, path = "~/Desktop", width = 10, height = 10, units = "cm")

#Make html plot with plotly

library(plotly)
libQCplotly <- ggplotly(libQC)
htmlwidgets::saveWidget(libQCplotly, "~/Desktop/libQCplotly.html")

#Identify contaminating ASVs - define control samples and threshold (e.g. 0.05)

library(decontam)

df <- as.data.frame(sample_data(phyloseq_object_all)) # Put sample_data into a ggplot-friendly data.frame
sample_data(phyloseq_object_all)$is.neg <- sample_data(phyloseq_object_all)$"SampleType" == "control"
contamdf.prev <- isContaminant(phyloseq_object_all, method="prevalence", neg="is.neg", threshold = 0.05)

#Identify contaminants

table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))


con_ASVs <- contamdf.prev %>% 
  filter(contaminant == "TRUE")
con_ASVs  <- rownames(con_ASVs)

#Take a look at the number of times several of these taxa were observed in negative controls and positive samples.

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(phyloseq_object_all, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(phyloseq_object_all)$"SampleType" == "control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$"SampleType" == "sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

# Make plot and save as pdf
deconplot <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
ggsave("deconplot.pdf", plot = deconplot, path = "~/Desktop", width = 10, height = 10, units = "cm")


#Make distribution plot of reads using microbiomeutilities

library(BiocManager)
BiocManager::install("microbiome")
library(microbiome)

#install.packages("devtools")
library(devtools)
devtools::install_github("microsud/microbiomeutilities")
library(microbiomeutilities)

library(knitr)
library(tibble)
library(dplyr)

distrib <- plot_read_distribution(phyloseq_object_all, groups = "site", 
                                  plot.type = "density") + xlab("Reads per sample") + ylab("Density")
distrib <- distrib + geom_density(alpha = 0.5, fill = "grey")
ggsave("distrib.pdf", plot = distrib, path ="~/Desktop", width = 15, height = 10, units = "cm")

hist(contamdf.prev$p, n=100)

#Reminder: ps_raw_bact file in Siobhon's script = phyloseq_object_all in my script
#If no contaminants identified using this method (apparently common if analysing low biomass samples), consider removing taxa found in the EBs and NCs from clinical samples)
#i.e. explore how much overlap there is between the taxa found in the extraction control and the true samples, and then perhaps remove most if not all of the taxa present in the negative control.
#Findings: all taxa in EB and NC are present in samples in higher numbers, so will not remove them manually. I would do if for instance there were taxa with 10 reads in clinical sample A, which had 50000 reads in the negative control.
#I will however inspect all ASVs present in EBs and NC and if they are environmental (i.e. not pathogens), I will create a separate file phyloseq_no_EBNCcont (ASVs in EBs, NCs) for additional analysis and compare at the end. Code below:

badTaxa = c("12cf958aac38df753bd74152c9de7b19")
allTaxa = taxa_names(phyloseq_object_all)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
phyloseq_no_cont = prune_taxa(allTaxa, phyloseq_object_all)
# new phyloseq object with just the taxa I kept= phyloseq_no_EBNC

#Save R data for phyloseq object - saving “raw data” (phyloseq_object_all) and “decontaminated data” (phyloseq_no_EBNCcont)
save(phyloseq_object_all, file = "~/Desktop/phyloseq_object_all.RData")
save(phyloseq_no_cont, file = "~/Desktop/phyloseq_no_cont.RData")

#If after inspection want to remove additional or all taxa from EBs and NCs, use "badTaxa" again and then the code below to convert key phyloseq objects into excel files. These will be created in the "home" directory:

phyloseq_no_cont@tax_table
otu_df <- phyloseq_no_cont@otu_table
write.csv (otu_df, "otu_df")

tax_df <- phyloseq_no_cont@tax_table
write.csv(tax_df, "tax_df")

sam_df <- phyloseq_no_cont@sam_data
write.csv(sam_df,"sam_df.csv")

#NEXT STEPS: 3 AND 4 (http://siobhonlegan.com/wildlife-bacteria/phyloseq.html#3_Load_data_and_subset)
#NEXT STEPS: Microbiome visualisation and stats (https://microsud.github.io/microbiomeutilities/articles/microbiomeutilities.html#abundance-prevalence-relationship and http://siobhonlegan.com/wildlife-bacteria/microbiome-viz.html#3_Abundance-Prevalence_relationship)

#Add Brucelaceae in the taxa of interest group
#Do rarefied barplot/table? Rarefy phyloseq_object_all
