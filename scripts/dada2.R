#packages used
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", force = TRUE)
library(dada2)
library(phyloseq)
library(csv)
library(tidyverse)


#-----MAKING ASV AND TAXA TABLES-----#

#path to fastq files
path <- "/home/lgschaer/old/MWarke/New_Fastqs_05182022/Manas/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
length(list.files(path))/4

list_of_files <- list.files(path) %>%
  as.data.frame()
head(list_of_files)

#forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

length(fnFs)

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1)
head(sample.names)

#visualize quality profiles
plotQualityProfile(fnFs[150:153])           #forward reads
plotQualityProfile(fnRs[150:153])           #reverse reads

#place filtered files in "filtered", a subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))

#standard filtering parameters: 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = 40, trimRight = 70,
                     maxN=0, maxEE=Inf, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE


head(out)

out2 <- out %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  filter(reads.in<1000)

dim(out2)

ggplot(out2, aes(x=SampleID, y=reads.in))+
  geom_col()+
  theme(axis.text.x = element_text(angle=90))

read.stats <- out %>%
  as_tibble() %>%
  summarise(
    max.in = max(reads.in),
    min.in = min(reads.in),
    mean.in = mean(reads.in),
    median.in = median(reads.in),
    sum.in = sum(reads.in),
    count.in = sum(ifelse(reads.in > 0, 1, 0)),
    under1000.in = sum(ifelse(reads.in < 1000, 1, 0)),
    max.out = max(reads.out),
    min.out = min(reads.out),
    mean.out = mean(reads.out),
    median.out = median(reads.out),
    sum.out = sum(reads.out),
    count.out = sum(ifelse(reads.out > 0, 1, 0)),
    under1000.out = sum(ifelse(reads.out < 1000, 1, 0))
  )
View(t(read.stats))

#save csvs
write.csv(out, "/home/lgschaer/old/MWarke/dada2_output/out.csv")
write.csv(read.stats, "/home/lgschaer/old/MWarke/dada2_output/read.stats.csv")

#learn about error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#plot errors
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

#name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#sample inferance
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#inspecting dada-class object
dadaFs[[1]]
dadaRs[[1]]

#merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate = TRUE, verbose=TRUE)

#inspect the merger data.frame from the first sample
head(mergers[[1]])

#construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

View(track)


#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/home/lgschaer/old/Silva_Versions/silva_nr_v138_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


#save sequence table and taxa table
saveRDS(seqtab.nochim, "/home/lgschaer/old/MWarke/dada2_output/seqtab.rds")                  #sequence table
saveRDS(taxa.print, "/home/lgschaer/old/MWarke/dada2_output/taxa.rds")          #taxa table


