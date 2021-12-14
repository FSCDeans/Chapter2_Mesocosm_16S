#Enable packages
library(dada2)
library(phyloseq)  
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(plyr)
library(dplyr)
library(scales)
library(reshape)
library(reshape2)
library(RColorBrewer)
library(Rmisc)
library(grid)
library(microbiome)
library(graphics)
library(tidyr)
library('ggedit')
library(ggrepel)



#Update path inside R
path <- "/SCRATCH/Sven_Data/DADA2/FD_2017/Samples/"

#List the files we have
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


#Inspect read quality
pdf("Plotquality.pdf")
plotQualityProfile(fnFs[1])
plotQualityProfile(fnRs[1]) #Reverse is worse quality - common in illumina sequencing
dev.off()

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Filter and Trim - need to match IDs as sometimes not all reverse sequences pass quality checks and therefore the forward and reverse fastq files sequencing order do not match
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE, matchIDs=TRUE) # On Windows set multithread=FALSE
head(out)

#Issues during previous step can occur - usually due to a replicate sample or it randomly assigning a different name, I have no idea why it does it, but you can check for it using:
#any(duplicated(c(fnFs, fnRs)))
#any(duplicated(c(filtFs, filtRs)))
#Both should be false, if they are and previous step still fails it means that a different name has been assigned to files and you'll have to manually check file names and look for inconsistencies.
#This could also be because barcodes during demultiplexing were mismatched (try deleting one of the reverse options when re-running) and thus truncation cut off got rid of all reads.


#Learn the error rate
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#Sanity check
pdf("Error_rates.pdf")
plotErrors(errF, nominalQ=TRUE)
plotErrors(errF, nominalQ=TRUE)
dev.off()

#If some sequencing does not pass the quality checks their files are not generated and the dada function fails, so we need to ensure that all inputs exist
exists <- file.exists(filtFs)

#Dereplicate only fastq files that exist
derepFs <- derepFastq(filtFs[exists], verbose=TRUE)
derepRs <- derepFastq(filtRs[exists], verbose=TRUE)

#Only retain existing fastq files
names(derepFs) <- sample.names[exists]
names(derepRs) <- sample.names[exists]

##Dereplicate the fastq file and carry out sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE) #Maybe look at pooling in the future, see how it affects ASV creation?
dadaRs <- dada(derepRs, err=errR, multithread=TRUE) #Maybe look at pooling in the future, see how it affects ASV creation?

dadaFs[[1]]


##Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


##Construct a sequence Table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


##Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)


##Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, "Track_reads_through_pipeline.csv")


##Assign Taxonomy
#Only goes down to family
taxa <- assignTaxonomy(seqtab.nochim, paste(path, "taxa_DB/silva_nr_v132_train_set.fa.gz", sep = ""), multithread=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#Based on exact matching and goes down to species if possible
taxa <- addSpecies(taxa, paste(path, "taxa_DB/silva_species_assignment_v132.fa.gz", sep = ""))

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

write.csv(taxa.print, "Taxa.csv")



##Convert to phyloseq object
map_file <- paste(path, "Mapping_file.txt", sep = "")
bmsd <- import_qiime_sample_data(map_file)

samples.out <- rownames(seqtab.nochim)

#You may have to manipulate this a little to add .fastq to the end of the mapping file names, but that was due to the names and should not be a problem if names are carefully chosen.
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(bmsd[exists]), 
               tax_table(taxa))

#Store DNA sequences of ASVs in the phyloseq refseq slot, renaming taxa to a short string
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

saveRDS(ps, file=paste(path, "Phyloseq_object.rds", sep = ""))
