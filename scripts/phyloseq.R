#-----MAKING A PHYLOSEQ OBJECT-----#

library(ggpubr)
library(tidyverse)
library(phyloseq)
library(csv)
library(pwr)
library(FSA)
library(tsnemicrobiota)
library(vegan)
library(ape)

# Load data into R for phyloseq analysis
# We will need (1) sample/meta data, (2) sequence table, and (3) taxa table
# (1) is a csv file with information about each sample (media, carbon, conditions, enrichment details, etc)
# (2) and (3) are output files from dada2.


#### TROUBLESHOOTING MISSING FILES

path <- "/home/lgschaer/old/MWarke/New_Fastqs_05182022/Manas/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

list_of_files <- list.files(path) %>%
  as_tibble() %>%
  mutate(
    File_Name = value
  ) %>%
  separate(value, into = c("SampleID", "File_Extension"), sep = "_S")
head(list_of_files)

sdata <- as.csv("/home/lgschaer/old/MWarke/samples_metadata_soilarsenic_miceobiome_06222022.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
head(sdata)

sdataA <- sdata %>%
  rownames_to_column(var = "SampleID") %>%
  mutate(
    SampleID = gsub("_", "-", SampleID)
  ) %>%
  left_join(list_of_files, by = "SampleID") %>%
  filter(is.na(File_Extension))
View(sdataA)

#### START PHYLOSEQ

#load sample data
sdata <- as.csv("/home/lgschaer/old/MWarke/samples_metadata_soilarsenic_miceobiome_06222022.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
head(sdata)
dim(sdata)

# reformat metadata (may or may not be necessary depending on how it is formatted)
sdata2 <- sdata %>% 
  rownames_to_column(var = "SampleID") %>%
  mutate(
    SampleID = gsub("_", "-", SampleID),
    SampleID = gsub(" ppm-3", "-3", SampleID),
    rownames = SampleID,
    Arsenic.concentration = as.character(Arsenic.concentration)
  ) %>%
  column_to_rownames(var = "rownames")
head(sdata2)                                                       #view data to make sure everything is OK

#load sequence table
sequence_table <- readRDS("/home/lgschaer/old/MWarke/dada2_output/seqtab.rds")
dim(sequence_table)
colnames(sequence_table) <- NULL                                   #remove column names
sequence_table[1:5, 1:5]

#make nonzero subset to remove all columns with zero taxa counts
sequence_table <- as.matrix(sequence_table)                                #change to matrix format
m <- (colSums(sequence_table, na.rm=TRUE) != 0)                        #T if colSum is not 0, F otherwise
nonzero <- sequence_table[, m]                                         #all the non-zero columns
nonzero[1:5,1:5]

#load taxa table
taxa_table <- readRDS("/home/lgschaer/old/MWarke/dada2_output/taxa.rds")
taxa_table <- as.matrix(taxa_table)                                #change to matrix format
taxa_table[1:5,1:5]                                                #view to make sure everything looks good

#make phyloseq object
samdata = sample_data(sdata2)                                      #define sample data
colnames(nonzero) <- NULL                                          #remove column names from "nonzero"
seqtab = otu_table(nonzero, taxa_are_rows = FALSE)                 #define sequence table
taxtab = tax_table(taxa_table)                                     #define taxa table
rownames(taxtab) <- NULL                                           #remove rownames from taxa table


sample_names(samdata)
sample_names(seqtab)

phyloseq_object_all = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata))
phyloseq_object_all


justbacteria <- phyloseq_object_all %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "Mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  ) 
justbacteria

any(sample_sums(justbacteria) == 0)
sum(sample_sums(justbacteria) == 0)

any(sample_sums(justbacteria) < 1000)
sum(sample_sums(justbacteria) < 1000)


#normalize data
#Delete samples with a mean of less than 1000
samplesover1000_all <- subset_samples(justbacteria, sample_sums(justbacteria) > 1000)
min(sample_sums(samplesover1000_all))

#Check if there are OTUs with no counts, if so how many?
any(taxa_sums(samplesover1000_all) == 0)
sum(taxa_sums(samplesover1000_all) == 0)

#Prune OTUs with no counts 
prune_samplesover1000_all <- prune_taxa(taxa_sums(samplesover1000_all) > 0, samplesover1000_all)
any(taxa_sums(prune_samplesover1000_all) == 0)

#make sure seed is set the same each time, set to 81 here
rarefy_samplesover1000_all <- rarefy_even_depth(prune_samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(prune_samplesover1000_all)))

#number rarefied at:
min(sample_sums(prune_samplesover1000_all))
justbacteria2 <- rarefy_samplesover1000_all
justbacteria2

#no_filt_metadata <- sample_data(phyloseq_object_all) %>%
 #as_tibble() 
#head(no_filt_metadata)

#filt_metadata <- sample_data(justbacteria2) %>%
 #as_tibble() 
#head(filt_metadata)

#only_mitochondria_samples <- no_filt_metadata %>%
 #anti_join(filt_metadata)
#head(only_mitochondria_samples)

#dim(only_mitochondria_samples)
#362-246

#write_csv(only_mitochondria_samples, "/home/lgschaer/old/MWarke/phyloseq_output/samples_with_only_mitochondria.csv")


#head(sample_data(test))


#sample counts
sample_counts <- sample_data(justbacteria) %>%
  group_by(Timepoint, Fern.type, Arsenic.concentration, Sample.type) %>%
  mutate(Count = 1) %>%
  summarise(SumCount = sum(Count))
View(sample_counts)

as.csv(sample_counts, "/home/lgschaer/old/MWarke/phyloseq_output/sample_counts_06222022.csv")

head(sample_data(justbacteria))

justbacteria_filt <- justbacteria %>%
  subset_samples(Sample.type != "shoot" & (Fern.type == "PE" | Fern.type == "PV"))
justbacteria_filt
justbacteria

#saveRDS file
saveRDS(justbacteria, file = "/home/lgschaer/old/MWarke/phyloseq_output/phyloseq_rarefied_nochloroplasts_07202022.rds")
saveRDS(justbacteria_filt, file = "/home/lgschaer/old/MWarke/phyloseq_output/phyloseq_rarefied_nochloroplasts_noPE_noT7_07202022.rds")

#-----ALPHA DIVERISTY-----#

#Violin plot of alpha diversity Observed OTUs and Shannon Diversity (with color)
colors <- c("lightblue", "purple", "firebrick", "orange", "cyan")
greens <- c("greenyellow", "chartreuse3", "chartreuse4")
greens2 <- c("greenyellow", "chartreuse", "chartreuse3", "chartreuse4")

colors <- c("firebrick", "dodgerblue")

PlotA <- justbacteria_filt %>%                                                     #phyloseq object
  plot_richness(
    x = "Fern.type",                                                #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  geom_violin(aes(fill = Fern.type), show.legend = FALSE)+          #make violin plot, set fill aes to sampletype
  geom_boxplot(width=0.1) +                                          #add boxplot, set width
  theme_classic()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  scale_x_discrete(limits = c("PV", "PE"))+
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = colors, breaks=c("PV", "PE"))+
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position
PlotA

sample_type_colors <- c("lightblue", "green", "firebrick", "tan4")

justbacteria_filt %>%                                                     #phyloseq object
  plot_richness(
    x = "Timepoint",                                                #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  geom_boxplot(show.legend = FALSE) +                                          #add boxplot, set width
  geom_jitter(aes(shape = Sample.type, fill = Fern.type), color = "black", size = 5)+
  scale_shape_manual(values = c(21, 23, 24)) +
  theme_classic()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  scale_x_discrete(limits = c("T0", "T6", "T7", "T12"))+
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size =  20, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 18))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = colors)+                         #set fill colors
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position



justbacteria_filt %>%                                                     #phyloseq object
  plot_richness(
    x = "Sample.type",                                                #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  geom_boxplot(show.legend = FALSE) +                                          #add boxplot, set width
  geom_jitter(aes(fill = Fern.type), color = "black", size = 5, shape = 21)+
  #scale_shape_manual(values = c(21, 23, 24)) +
  theme_classic()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  #scale_x_discrete(limits = c("T0", "T6", "T7", "T12"))+
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size =  20, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 18))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = colors)+                         #set fill colors
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position


blues <- c("lightblue", "dodgerblue", "darkblue")

justbacteria_filt %>%                                                     #phyloseq object
  plot_richness(
    x = "Timepoint",                                                #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  geom_boxplot(show.legend = FALSE) +                                          #add boxplot, set width
  geom_jitter(aes(fill = Arsenic.concentration, shape = Fern.type), color = "black", size = 5)+
  scale_shape_manual(values = c(21, 23, 24))+
  theme_classic()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  scale_x_discrete(limits = c("T0", "T6", "T7", "T12"))+
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size =  20, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 18))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = blues)+                         #set fill colors
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position


#-----ALPHA DIVERISTY STATISTICS-----#

#add alpha diversity data to a data frame
richness <- justbacteria_filt %>%
  estimate_richness(measures = c("Observed", "Shannon")) %>%           #specify which measures
  rownames_to_column(var = "SampleID") %>%                             #add column name to SampleID column
  as_tibble() %>%
  mutate(SampleID = gsub("\\.", "-", SampleID))
head(richness)

head(sdata2)

alphadiv <- richness %>%
  left_join(sdata2, by = "SampleID")
head(alphadiv)

#alpha diversity summary statistics

alphadiv_summary <- alphadiv %>%
  group_by(Fern.type, Sample.type) %>%
  summarise(
    minShannon = min(Shannon),
    meanShannon = mean(Shannon),
    maxShannon = max(Shannon),
    sdShannon = sd(Shannon),
    minObserved = min(Observed),
    meanObserved = mean(Observed),
    maxObserved = max(Observed),
    sdObserved = sd(Observed)
  )
View(alphadiv_summary)

#Kruskal-Wallis Test
set.seed(81)

##BY FERN TYPE
#Observed
kruskal.test(Observed ~ Fern.type, data = alphadiv) 

#Shannon
kruskal.test(Shannon ~ Fern.type, data = alphadiv)


##BY SAMPLE TYPE
head(alphadiv)

#Observed
kruskal.test(Observed ~ Sample.type, data = alphadiv) 

#Shannon
kruskal.test(Shannon ~ Sample.type, data = alphadiv) 


#Dunn test (post hoc)

##Observed
dunnObserved <- dunnTest(Observed ~ Sample.type,
                  data=alphadiv,
                  method="bh")
dunnObserved


##Shannon
dunnShannon <- dunnTest(Shannon ~ Sample.type,
                  data=alphadiv,
                  method="bh")
dunnShannon


##BY ARSENIC CONCENTRATION
head(alphadiv)

#Observed
kruskal.test(Observed ~ Arsenic.concentration, data = alphadiv) 

#Shannon
kruskal.test(Shannon ~ Arsenic.concentration, data = alphadiv) 



##BY TIME
head(alphadiv)

#Observed
kruskal.test(Observed ~ Timepoint, data = alphadiv) 

#Shannon
kruskal.test(Shannon ~ Timepoint, data = alphadiv) 


#Dunn test (post hoc)

##Observed
dunnObserved <- dunnTest(Observed ~ Timepoint,
                         data=alphadiv,
                         method="bh")
dunnObserved


##Shannon
dunnShannon <- dunnTest(Shannon ~ Timepoint,
                        data=alphadiv,
                        method="bh")
dunnShannon

#-----BETA DIVERSITY-----#

#PCOA Plot

#ordination
all_pcoa <- ordinate(
  physeq = justbacteria_filt, 
  method = "PCoA", 
  distance = "bray"
)

colors3 <- c("green", "blue", "lightblue")

#plot
PlotB <- plot_ordination(
  physeq = justbacteria_filt,                                                          #phyloseq object
  ordination = all_pcoa)+                                                #ordination
  geom_point(aes(fill = Sample.type), size = 6, shape = 21) +                         #sets fill color to sampletype
  scale_fill_manual(values = sample_type_colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank())+                                      #removes legend title
    #legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.position = "bottom",
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
PlotB

#plot
plot_ordination(
  physeq = justbacteria_filt,                                                          #phyloseq object
  ordination = all_pcoa)+                                                #ordination
  geom_point(aes(shape = Sample.type, fill = Timepoint), size = 6) +                         #sets fill color to sampletype
  #scale_fill_manual() +
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_fill_manual(values = greens, breaks=c("T0","T6","T12"))+
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank())+                                      #removes legend title
  #legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.position = "bottom",
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command


#-----BETA DIVERSTIY STATISTICS-----#

#PERMANOVA
set.seed(81)
head(sample_data(justbacteria2))

#looking at differences between enrichments
#subset phyloseq object
E2_L1 <- subset_samples(justbacteria2, Enrichment %in% c("Emma2", "Laura1"))
E2_L2 <- subset_samples(justbacteria2, Enrichment %in% c("Emma2", "Laura2"))
L1_L2 <- subset_samples(justbacteria2, Enrichment %in% c("Laura1", "Laura2"))

# Calculate bray curtis distance matrix, all samples
E2_L1bray1 <- phyloseq::distance(E2_L1, method = "bray")
E2_L2bray2 <- phyloseq::distance(E2_L2, method = "bray")
L1_L2bray3 <- phyloseq::distance(L1_L2, method = "bray")

# make a data frame from the sample_data, all samples
E2_L1sam1 <- data.frame(sample_data(E2_L1))
E2_L2sam2 <- data.frame(sample_data(E2_L2))
L1_L2sam3 <- data.frame(sample_data(L1_L2))

# Adonis test, all samples
adonis(E2_L1bray1 ~ Enrichment, data = E2_L1sam1)
adonis(E2_L2bray2 ~ Enrichment, data = E2_L2sam2)
adonis(L1_L2bray3 ~ Enrichment, data = L1_L2sam3)

#looking at differences between media/carbon type within each enrichment
#subset phyloseq object
E2<-subset_samples(justbacteria2, Enrichment == "Emma2")
L1<-subset_samples(justbacteria2, Enrichment == "Laura1")
L2<-subset_samples(justbacteria2, Enrichment == "Laura2")

e2one <- subset_samples(E2, Media_Carbon %in% c("ARW_DCPET", "AMW_DCPET"))
e2two <- subset_samples(E2, Media_Carbon %in% c("ARW_DCPET", "BH_DCPET"))
e2three <- subset_samples(E2, Media_Carbon %in% c("ARW_DCPET", "BH_TPA"))
e2four <- subset_samples(E2, Media_Carbon %in% c("AMW_DCPET", "BH_DCPET"))
e2five <- subset_samples(E2, Media_Carbon %in% c("AMW_DCPET", "BH_TPA"))
e2six <- subset_samples(E2, Media_Carbon %in% c("BH_DCPET", "BH_TPA"))

l1one <- subset_samples(L1, Media_Carbon %in% c("ARW_DCPET", "AMW_DCPET"))
l1two <- subset_samples(L1, Media_Carbon %in% c("ARW_DCPET", "BH_DCPET"))
l1three <- subset_samples(L1, Media_Carbon %in% c("ARW_DCPET", "BH_TPA"))
l1four <- subset_samples(L1, Media_Carbon %in% c("AMW_DCPET", "BH_DCPET"))
l1five <- subset_samples(L1, Media_Carbon %in% c("AMW_DCPET", "BH_TPA"))
l1six <- subset_samples(L1, Media_Carbon %in% c("BH_DCPET", "BH_TPA"))

l2one <- subset_samples(L2, Media_Carbon %in% c("ARW_DCPET", "AMW_DCPET"))
l2two <- subset_samples(L2, Media_Carbon %in% c("ARW_DCPET", "BH_DCPET"))
l2three <- subset_samples(L2, Media_Carbon %in% c("ARW_DCPET", "BH_TPA"))
l2four <- subset_samples(L2, Media_Carbon %in% c("AMW_DCPET", "BH_DCPET"))
l2five <- subset_samples(L2, Media_Carbon %in% c("AMW_DCPET", "BH_TPA"))
l2six <- subset_samples(L2, Media_Carbon %in% c("BH_DCPET", "BH_TPA"))


# Calculate bray curtis distance matrix, all samples
E2bray1 <- phyloseq::distance(e2one, method = "bray")
E2bray2 <- phyloseq::distance(e2two, method = "bray")
E2bray3 <- phyloseq::distance(e2three, method = "bray")
E2bray4 <- phyloseq::distance(e2four, method = "bray")
E2bray5 <- phyloseq::distance(e2five, method = "bray")
E2bray6 <- phyloseq::distance(e2six, method = "bray")

L1bray1 <- phyloseq::distance(l1one, method = "bray")
L1bray2 <- phyloseq::distance(l1two, method = "bray")
L1bray3 <- phyloseq::distance(l1three, method = "bray")
L1bray4 <- phyloseq::distance(l1four, method = "bray")
L1bray5 <- phyloseq::distance(l1five, method = "bray")
L1bray6 <- phyloseq::distance(l1six, method = "bray")

L2bray1 <- phyloseq::distance(l2one, method = "bray")
L2bray2 <- phyloseq::distance(l2two, method = "bray")
L2bray3 <- phyloseq::distance(l2three, method = "bray")
L2bray4 <- phyloseq::distance(l2four, method = "bray")
L2bray5 <- phyloseq::distance(l2five, method = "bray")
L2bray6 <- phyloseq::distance(l2six, method = "bray")

# make a data frame from the sample_data, all samples
E2sam1 <- data.frame(sample_data(e2one))
E2sam2 <- data.frame(sample_data(e2two))
E2sam3 <- data.frame(sample_data(e2three))
E2sam4 <- data.frame(sample_data(e2four))
E2sam5 <- data.frame(sample_data(e2five))
E2sam6 <- data.frame(sample_data(e2six))

L1sam1 <- data.frame(sample_data(l1one))
L1sam2 <- data.frame(sample_data(l1two))
L1sam3 <- data.frame(sample_data(l1three))
L1sam4 <- data.frame(sample_data(l1four))
L1sam5 <- data.frame(sample_data(l1five))
L1sam6 <- data.frame(sample_data(l1six))

L2sam1 <- data.frame(sample_data(l2one))
L2sam2 <- data.frame(sample_data(l2two))
L2sam3 <- data.frame(sample_data(l2three))
L2sam4 <- data.frame(sample_data(l2four))
L2sam5 <- data.frame(sample_data(l2five))
L2sam6 <- data.frame(sample_data(l2six))

# Adonis test, all samples
adonis(E2bray1 ~ Media_Carbon, data = E2sam1)
adonis(E2bray2 ~ Media_Carbon, data = E2sam2)
adonis(E2bray3 ~ Media_Carbon, data = E2sam3)
adonis(E2bray4 ~ Media_Carbon, data = E2sam4)
adonis(E2bray5 ~ Media_Carbon, data = E2sam5)
adonis(E2bray6 ~ Media_Carbon, data = E2sam6)

adonis(L1bray1 ~ Media_Carbon, data = L1sam1)
adonis(L1bray2 ~ Media_Carbon, data = L1sam2)
adonis(L1bray3 ~ Media_Carbon, data = L1sam3)
adonis(L1bray4 ~ Media_Carbon, data = L1sam4)
adonis(L1bray5 ~ Media_Carbon, data = L1sam5)
adonis(L1bray6 ~ Media_Carbon, data = L1sam6)

adonis(L2bray1 ~ Media_Carbon, data = L2sam1)
adonis(L2bray2 ~ Media_Carbon, data = L2sam2)
adonis(L2bray3 ~ Media_Carbon, data = L2sam3)
adonis(L2bray4 ~ Media_Carbon, data = L2sam4)
adonis(L2bray5 ~ Media_Carbon, data = L2sam5)
adonis(L2bray6 ~ Media_Carbon, data = L2sam6)



#EXPLORING TAXA

sample_sums(justbacteria)
sample_sums(test)
taxa_sums(test)

subset_taxa(test, Family == "Mitochondria")
test

#Summarize abundance of each class
genusabundance <- justbacteria %>%
  tax_glom(taxrank = "Family") %>%                      # agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Family) 
head(genusabundance)


genusabundance_test <- genusabundance %>%
  filter(Fern.type == "PV" & Timepoint == "T7" & Arsenic.concentration == 0 & Sample.type == "shoot") %>%
  filter(Abundance != 0)
View(genusabundance_test)

#save color palatte
colors10 <- c(
  "black",   "darkcyan",     "orchid1",   "green",       "blue",   
  "grey47",  "cyan",    "coral1",     "yellow",    "darkgreen",   "palegoldenrod",    
  "grey77",  "darkblue",     "orange",    "red",         "mediumpurple1", "tan4",   "purple4",
  "dodgerblue",    "white", "firebrick", "yellowgreen", "magenta", "blue", "green", "red", "orchid", "lightblue"
) 

head(genusabundance)
#Select and summarize necessary variables
all <- genusabundance %>%
  select(Phylum, Class, Family, Sample, Abundance, Timepoint, Fern.type, Arsenic.concentration, Sample.type, Replicates) %>%
  filter(Abundance != 0) %>%
  mutate(
    Phylum = as.character(Phylum),
    Class = as.character(Class),
    Family = as.character(Family)#,
    #Genus = as.character(Genus)
    )
head(all)

phylum <- all %>%
  dplyr::group_by(Sample, Timepoint, Fern.type, Sample.type, Arsenic.concentration, Phylum)%>%
  summarise(Abundance = sum(Abundance)/n()) %>%
  mutate(Phylum.1p = ifelse(Abundance < 0.01, "<1%", Phylum)) #%>%
  #mutate(Fern.type = ifelse(is.na(Fern.type), "Control", Fern.type))
View(phylum)

range(phylum$Abundance)

class <- all %>%
  group_by(Timepoint, Fern.type, Arsenic.concentration, Sample.type, Class)%>%
  summarise(Abundance = sum(Abundance)/n())  %>%
  mutate(Class.1p = ifelse(Abundance < 0.01, "<1%", Class)) %>%
  filter(!is.na(Sample.type)) %>%
  mutate(Fern.type = ifelse(is.na(Fern.type), "Control", Fern.type)) #%>%
  #mutate(Class = ifelse(Abundance < 0.01, "<1%", Class))
head(class)

family <- all %>%
  group_by(Timepoint, Fern.type, Sample.type, Family)%>%
  summarise(Abundance = sum(Abundance)/n()) %>%
  mutate(Family.5p = ifelse(Abundance < 0.05, "<5%", Family))  %>%
  filter(!is.na(Sample.type)) %>%
  mutate(Fern.type = ifelse(is.na(Fern.type), "Control", Fern.type))
View(family)

#genus <- all %>%
#  group_by(Timepoint, Fern.type, Sample.type, Genus)%>%
 # summarise(Abundance = sum(Abundance)/n()) %>%
  #mutate(Genus.13p = ifelse(Abundance < 0.13, "<13%", Genus))  %>%
#  filter(!is.na(Sample.type)) %>%
 # mutate(Fern.type = ifelse(is.na(Fern.type), "Control", Fern.type))
#head(genus)


#MAKING A TAXA PLOT 

#save color palette
colors9 <- c(
  "black",      "yellow",    "olivedrab",  "darkcyan",       "coral3",   
"cyan",         "blue",      "grey47",    "tomato",       "darkgreen",   "palegoldenrod",    
  "grey77",     "darkblue",     "orange",           "mediumpurple1", "tan4",   "purple4",
  "dodgerblue", "red",    "yellowgreen", "magenta3",     "white",  "orchid1", "lightblue",   "black",      "yellow",    "olivedrab",  "darkcyan",       "coral3",   
"cyan",         "blue",      "grey47",    "tomato",       "darkgreen",   "palegoldenrod",    
"grey77",     "darkblue",     "orange",           "mediumpurple1", "tan4",   "purple4",
"dodgerblue", "red",    "yellowgreen", "magenta3",     "white",  "orchid1", "lightblue"
)  


length(colors10)

#plot

dim(phylum)

phy <- ggplot(phylum)+
  geom_col(mapping = aes(x = Timepoint, y = Abundance, fill = Phylum), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(cols = vars(Fern.type, Sample.type), rows = vars(Arsenic.concentration), scales = "free")+
  scale_fill_manual(values = colors9) +
  xlab(NULL)+
  theme_minimal()+
  scale_x_discrete(limits = c("T0", "T6", "T7", "T12"))+
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, angle = 0, vjust = 0.5, hjust = 0.5, color = "black"),
        axis.title.y = element_text(size = 15),
        axis.title.y.right = element_text(size = 15),
        axis.text.y.right = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "right",
        #plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 15, face = "bold", angle = 0),
        strip.text.y = element_text(size = 15, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 20))
phy

head(phylum)

ggplot(phylum)+
  geom_col(mapping = aes(x = Sample, y = Abundance, fill = Phylum), color = "black", position = "fill", show.legend = TRUE)+
  #facet_grid(cols = vars(Fern.type, Sample.type), rows = vars(Arsenic.concentration), scales = "free")+
  scale_fill_manual(values = colors9) +
  xlab(NULL)+
  theme_minimal()+
 # scale_x_discrete(limits = c("T0", "T6", "T7", "T12"))+
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.title.y = element_text(size = 15),
        axis.title.y.right = element_text(size = 15),
        axis.text.y.right = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "bottom",
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 15, face = "bold", angle = 0),
        strip.text.y = element_text(size = 15, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 20))

cla <- ggplot(class)+
  geom_col(mapping = aes(x = Timepoint, y = Abundance, fill = Class.1p), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(cols = vars(Fern.type, Sample.type), rows = vars(Arsenic.concentration), scales = "free")+
  scale_fill_manual(values = colors9) +
  scale_y_continuous(name = "Proportion of Community", sec.axis = sec_axis(trans=~.*0, name="Arsenic Concentration (mg/kg)"))+
  xlab(NULL)+
  theme_minimal()+
  scale_x_discrete(limits = c("T0", "T6", "T7", "T12"))+
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, angle = 0, vjust = 0.5, hjust = 0.5, color = "black"),
        axis.title.y = element_text(size = 15),
        axis.title.y.right = element_text(size = 15),
        axis.text.y.right = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "right",
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 15, face = "bold", angle = 0),
        strip.text.y = element_text(size = 15, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 20))
cla

#summaries

sample_counts <- sample_data(justbacteria2) %>%
  group_by(Enrichment) %>%
  mutate(Count = 1)%>%
  summarise(SumCount = sum(Count))# %>%
#mutate(countpercategory = mean(SumCount))
sample_counts


phy.summary <- phylum  %>%
  mutate(totalSum = sum(Abundance),
         RelAb = Abundance/totalSum) %>%
  filter(Phylum =="Actinobacteriota" | Phylum == "Proteobacteria") %>%
  ungroup() %>%
  dplyr::group_by(Enrichment, Phylum) %>%
  summarise(
    maxRelAb = max(RelAb),
    minRelAb = min(RelAb)
    ) %>%
  unique()
View(phy.summary)

class.summary <- class  %>%
  mutate(totalSum = sum(Abundance),
         RelAb = Abundance/totalSum) %>%
  ungroup() %>%
  dplyr::group_by(Enrichment, Class) %>%
  summarise(
    maxRelAb = format(max(RelAb), scientific = F, digits = 3),
    minRelAb = format(min(RelAb), scientific = F, digits = 3)
  ) 
View(class.summary)


family.genus.summary <- all  %>%
  group_by(Enrichment, Family, Genus)%>%
  summarise(maxRelAb = format(max(Abundance), scientific = FALSE, digits = 3),
            meanRelAb = format(mean(Abundance), scientific = FALSE, digits = 3),
            minRelAb = format(min(Abundance), scientific = FALSE, digits = 3)) 
View(family.genus.summary)
write_csv(family.genus.summary, "/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/03172022_DCPET_Full/phyloseq_output/max_mean_min_relative_abundance_family_genus_level.csv")


#genus by replicate taxa plot
genus2 <- all %>%
  group_by(Enrichment, Media_Carbon, Genus, Replicate)%>%
  summarise(Abundance = sum(Abundance)) %>%
  mutate(Genus.2p = ifelse(Abundance < 0.02, "<2%", Genus))#%>%
  #arrange(desc(Abundance, Sample))                           #Arrange with descending abundances
head(genus2)

ggplot(genus2)+
  geom_col(mapping = aes(x = Replicate, y = Abundance, fill = Genus.2p), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(rows = vars(Enrichment), cols = vars(Media_Carbon))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors10) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        title = element_text(size = 25))

# A closer look at Laura2 at the genus level

L2_genus <- all %>%
    mutate(Genus.2p = ifelse(Abundance < 0.02, "<2%", Genus)) %>%
    filter(Enrichment == "Laura2")
head(L2_genus)

ggplot(L2_genus)+
  geom_col(mapping = aes(x = Replicate, y = Abundance, fill = Genus.2p), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(cols = vars(Media_Carbon))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors10) +
  scale_x_discrete(labels = Media_Carbon_Labels)+
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5, color = "black"),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 12, face = "bold", angle = 0),
        #  legend.title = element_blank(),
        title = element_text(size = 18))
