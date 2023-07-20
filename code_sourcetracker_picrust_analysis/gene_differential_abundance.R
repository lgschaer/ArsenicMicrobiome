#load packages
library(tidyverse)
library(csv)
library(phyloseq)
library(DESeq2) #differential abundance
library("KEGGREST") #getting KO annotations



#load phyloseq object
ps2 <- read_rds("/home/lgschaer/old/MWarke_temp/unrarefied_T7soilR_phyloseq_07192023")
ps2
#View(otu_table(ps)[1:5,1:5])
#ps2 <- subset_samples(ps, Timepoint == "T7")
unique(sample_data(ps2)$Fern.type)


#load sample data and format for phyloseq
sdata <- as_tibble(sample_data(ps2)) %>%
  mutate(SampleID = str_replace_all(SampleID, "_", "-"),
         SampleID2 = SampleID) %>%
  column_to_rownames(var = "SampleID2") %>%
  filter(!is.na(Fern.type))
head(sdata)
#View(sdata)
unique(sdata$Fern.type)
dim(sdata)


#load picrust THIS WILL BE THE SEQUENCE TABLE
pred_mg <- read_tsv("/home/lgschaer/old/MWarke_temp/picrust2_out_07192023/KO_metagenome_out/pred_metagenome_unstrat.tsv") 
#pred_mg <- read_table("/home/lgschaer/old/MWarke_temp/picrust2_out_07192023/pathways_out/path_abun_unstrat.tsv")
pred_mg[1:10,1:10]
colnames(pred_mg)

pmg2 <- pred_mg %>%
  column_to_rownames(var = "function") %>%
  t()
pmg2[1:10,1:10]
rownames(pmg2)


# getting a list of gene names to match KO numbers THIS WILL BE THE TAXA TABLE
KO <- pred_mg %>%
  column_to_rownames(var = "function") %>%
  rownames_to_column(var = "KO") %>%
  select(KO)
head(KO)



#make phyloseq object
samdata = sample_data(sdata)                                      #define sample data
colnames(pmg2) <- NULL                                          #remove column names from "nonzero"
seqtab = otu_table(pmg2, taxa_are_rows = FALSE)                 #define sequence table
taxtab = tax_table(KO)                                     #define taxa table
rownames(taxtab) <- NULL                                           #remove rownames from taxa table


sample_names(samdata)
sample_names(seqtab)

psGene = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata)) 

#View(sample_data(phyloseq_object_all))
unique(sample_data(psGene)$Fern.type)

psGene <- subset_samples(psGene, sample_sums(psGene) > 0)
psGene <- subset_samples(psGene, taxa_sums(psGene) > 0)
psGene

####### Now we are ready for DESeq2

#saveRDS file
saveRDS(psGene, file = "/home/lgschaer/old/MWarke_temp/picrust_genes_ps_format_07132023.rds")
head(sample_data(psGene))

#add a pseudo count of one to the OTU table
#justbacteria@otu_table <- as.matrix(justbacteria@otu_table)+1


#View(sample_data(phyloseq_object_all))
unique(sample_data(psGene)$Fern.type)

ds.ps <- phyloseq_to_deseq2(psGene, ~ Fern.type)

ds.ps$Fern.type<-relevel(ds.ps$Fern.type, "PE")

# Run DESeq2 analysis (all taxa at once!)
dds_ps <- DESeq(ds.ps)

# Investigate results
resultsNames(dds_ps)

resultsNames(dds_ps)

colnames(dds_ps)

mcols(dds_ps,use.names=TRUE)

dds_ps$SampleID
dds_ps@colData
dds_ps@assays

#colnames(mcols(dds_ps))

# Put DESeq output into data frames
res <- as.data.frame(results(dds_ps, contrast=c("Fern.type","PV","PE"))) %>% mutate(Comparison="PE vs. PV") %>% rownames_to_column(var = "taxon")
head(res)

# Add taxonomy to DESeq output
res_tax <- as.data.frame(tax_table(phyloseq_object_all)) %>% rownames_to_column(var = "taxon") %>% full_join(res) %>%
  mutate(KO = ta1) %>%
  select(-ta1)
head(res_tax)

##check dimensions
dim(res)
dim(res_tax)

# Join everything together
enriched <- res_tax %>% 
  filter(!is.na(padj)) %>%
  mutate(
    threshold = ifelse(padj <= 0.01 & abs(log2FoldChange) >= 2, "Enriched", "Not_Enriched"),
    Where_Enriched = ifelse(log2FoldChange < 0, "PV", "PE")#,
    #Enriched_Genus = ifelse(is.na(Enriched_Genus), "Unclassified Genus", Enriched_Genus)
  ) %>%
  filter(threshold == "Enriched")
head(enriched)
#View(enriched)

# Are any ASVs enriched?
any(enriched$threshold == "Enriched")
sum(enriched$threshold == "Enriched")
sum(enriched$Where_Enriched == "PE")
sum(enriched$Where_Enriched == "PV")

# matching up the differential abundance significant genes to names

#KO <- c("K00001", "K00002", "K00003", "K00004", "K00005")
KO <- as.vector(enriched$KO)
length(KO)

#I think keggList() has a limit of 100 rows at a time, run it 3x to get all the genes

## if length(KO) is > 100
#KO1 <- KO[1:100] %>% keggList()# %>% as.data.frame()
#KO2 <- KO[101:200] %>% keggList() %>% as.data.frame()
#KO3 <- KO[201:279] %>% keggList() %>% as.data.frame()

## if length(KO) is < 100
KO <- KO %>% keggList() %>% as.data.frame()

KO_details <- KO   #merge all KO numbers together
dim(KO_details)
colnames(KO_details) <- c("details")
KO_details2 <- KO_details %>%
  rownames_to_column(var = "KO") %>%
  separate(details, into = c("Gene_Abbreviation", "Gene_Name"), sep = "; ") %>%
  separate(Gene_Name, into = c("Gene_Name", "EC_Number"), sep = " \\[") %>%
  separate(EC_Number, into = c("EC_Number", NA), sep = "\\]") %>%
  mutate(Gene_Prefix = str_extract(Gene_Abbreviation, "^..."),
         Variable = "x")
  #separate(Gene_Abbreviation, into = c("Gene_Prefix", "NA"), sep = "\\]")
head(KO_details2)
#View(KO_details2)

enriched_w_KO_details <- enriched %>%
  left_join(KO_details2, by = "KO")
head(enriched_w_KO_details)
dim(enriched)
dim(KO_details2)
dim(enriched_w_KO_details)

length(unique(enriched_w_KO_details$KO))
range(enriched_w_KO_details$log2FoldChange)

#View(enriched_w_KO_details)

# Save a csv of the results
write.csv(enriched,"/home/lgschaer/old/MWarke_temp/out_enriched_genes_KO_soilR.csv")
write.csv(KO_details2,"/home/lgschaer/old/MWarke_temp/KO_details_soilR.csv")
write.csv(enriched_w_KO_details,"/home/lgschaer/old/MWarke_temp/enriched_genes_w_KO_details_soilR.csv")

colFunc <- colorRampPalette(c("purple4", "orange"))
color_list <- colFunc(100)
length(color_list)

colnames(enriched_w_KO_details)
range(enriched_w_KO_details$log2FoldChange)


dim(enriched_w_KO_details)

#testTable <- enriched_w_KO_details[1:20,1:15]
head(enriched_w_KO_details)

ggplot(enriched_w_KO_details, aes(x=Variable, y = KO))+
  geom_tile(mapping = aes(fill = log2FoldChange), color = "black", show.legend = TRUE)+
  #facet_nested(cols = vars(), rows = vars(Pathway), space = "free", scales = "free")+
  #ylab("Proportion of Community") +
  scale_fill_gradientn(colors = color_list, breaks = c(-6.5, 5.361385), labels=c("PV", "PE")) +
  xlab(NULL)+
  ylab(NULL)+
  labs(fill=expression(paste(Log[2], " Fold Change")))+
  theme_linedraw()+ 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 18, face = "bold", angle = 0),
        strip.text.y = element_text(size = 18, face = "bold", angle = 0))+
  guides(fill=guide_colorbar(barwidth = 7, direction = "horizontal",  title.vjust = 1, title.hjust = 1, 
                             title.theme = element_text(face = "bold", size = 22), label.vjust = 0.5))#+
  #guides(fill=guide_legend(title=expression("Log[2]FoldChange")))



cl <- colors(distinct = TRUE)
set.seed(15887) # to set random generator seed
colors11 <- sample(cl, length(unique(enriched$Enriched_Genus)))

#colors11 <- c(
# "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
#"magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
#  "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
# "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
#"magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","gray63","white"
#)

shapes <- c(
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24, 
  21, 22, 23, 24, 
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  22, 21, 21
)
length(shapes)

#colors11 <- c(
# "orangered",      "purple",        "green",           "cyan",          "orange",        "khaki4",             "mediumslateblue",
#"mediumpurple1",  "darkmagenta",   "darkgreen",       "wheat2",        "yellow",        "lawngreen",          "plum",  
#  "royalblue",      "magenta",       "mediumseagreen",  "palegoldenrod", "grey47",        "chocolate4",         "darkorange3",        
# "lightblue",      "firebrick",     "yellowgreen",     "turquoise3",    "purple4",       "blue",               "red",            
#"lightcyan",       "coral1",       "cyan",            "goldenrod",     "yellowgreen",   "turquoise3",    "purple4",       "blue",               "red",            
#  "lightcyan",       "coral1",       "cyan",            "goldenrod",     "black",         "white"   
#) 
#length(colors11)

#colors11 <- c(
# "palegoldenrod","palegoldenrod","palegoldenrod","palegoldenrod",
#"orange",  "orange",  "orange",  "orange",  
#  "firebrick","firebrick","firebrick","firebrick",
# "pink","lightpink","lightpink","lightpink",
#"purple","purple","purple","purple",
#  "darkblue","darkblue","darkblue","darkblue",
# "darkcyan", "darkcyan", "darkcyan", "darkcyan", 
#"lightblue","lightblue","lightblue","lightblue",
#  "olivedrab2", "olivedrab2", "olivedrab2", "olivedrab2", 
# "darkgreen", "darkgreen", "darkgreen", "darkgreen", 
#"grey77","grey77","white"
#)

#volcano plot:
ggplot(data=enriched, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(fill=Enriched_Genus), shape = 21, color = "black", size=6) +
  facet_grid(cols = vars(Comparison))+
  #scale_fill_manual(values=colors11) +
  #scale_shape_manual(values=shapes) +
  labs(x = "log2 fold change", 
       y = "-log10 p-value") +
  theme_classic(base_size = 14)+
  geom_hline(yintercept = -log10(0.001), colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = 2, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = -2, colour="#990000", linetype="dashed")+ 
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 25, face = "bold", angle = 0),
        legend.position = "bottom",
        title = element_text(size = 18))#+
#guides(fill = guide_legend(override.aes = list(shape = shapes)))

