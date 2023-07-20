### Making inputs for picrust

#load libraries
library(phyloseq)
library(tidyverse)
#install.packages("phylotools", dependencies = TRUE)
library(phylotools)


#load sequence table
sequence_table <- readRDS("/home/lgschaer/old/MWarke_temp/seqtab.rds")
sequence_table[1:5,1:5]
dim(sequence_table)

pre_fna_summary <- sequence_table %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "seq.text") %>%
  mutate(seq.name = paste("sp", 1:67730, sep="")) %>%
  select(c("seq.name", "seq.text"))
head(pre_fna_summary)
dim(pre_fna_summary)


phylotools::dat2fasta(pre_fna_summary, outfile = "/home/lgschaer/old/MWarke_temp/arsenic.fna")


#### ATTEMPT 2: Using only T7 SoilR samples

ps <- readRDS("/home/lgschaer/old/MWarke/phyloseq_output/phyloseq_notrarefied_09132022.rds") %>%
  subset_samples(Timepoint == "T7" & Sample.type == "soilR")
ps

write_rds(ps, "/home/lgschaer/old/MWarke_temp/unrarefied_T7soilR_phyloseq_07192023")

head(sample_data(ps))
unique(sample_data(ps)$Sample.type)
unique(sample_data(ps)$Timepoint)

#extract the subsetted sequence table
sequence_table <- otu_table(ps)
sequence_table[1:5,1:5]
dim(sequence_table)

#we will use the original sequence table to get the sequences we need
og_sequence_table <- readRDS("/home/lgschaer/old/MWarke_temp/seqtab.rds") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "seq.text") %>%
  mutate(seq.name = paste("sp", 1:67730, sep="")) %>% #dimensions of OG sequence table
  select(c("seq.name", "seq.text"))
head(og_sequence_table)
dim(og_sequence_table)

pre_fna_summary <- sequence_table %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "seq.name") %>%
  select(c("seq.name")) %>%
  left_join(og_sequence_table, by = "seq.name")
head(pre_fna_summary)
dim(pre_fna_summary)


phylotools::dat2fasta(pre_fna_summary, outfile = "/home/lgschaer/old/MWarke_temp/arsenic.fna")
