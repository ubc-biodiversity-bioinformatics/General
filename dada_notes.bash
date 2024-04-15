######################
# dada notes - for phylogeny of environmental DNA
# ITS bands, CO1, 16S, 12S - need large divergences between groups, e.g. 
# Mar 2024
# Evan Morien
# https://github.com/morien/dada2-tutorial-16S
######################

# clone the github repo

git clone https://github.com/morien/dada2-tutorial-16S.git

module load StdEnv/2020 r/4.3.2 python/3.11.5

#--------------------------------------
#Create a new virtual environment and install Cutadapt into it:

virtualenv cutadapt-venv
cutadapt-venv/bin/pip install cutadapt

#Cutadapt is now available as cutadapt-venv/bin/cutadapt:
#cutadapt-venv/bin/cutadapt --version

source cutadapt-venv/bin/activate
cutadapt --version

#---------------------------------
# install packages

R

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.11")

install.packages("dada2")
install.packages("phyloseq")
install.packages("tidyverse")
install.packages("reshape2")
install.packages("stringr")
install.packages("data.table")
install.packages("ShortRead")
install.packages("Biostrings")
install.packages("seqinr")
install.packages("DECIPHER")

#-----------------------------------

#tutorial for for processing 16S MiSeq using dada2
#author: Evan Morien
#last modified: March 20th 2024
#working directory: /path/to/your/working_directory

####Libraries####
library(dada2)
library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(ShortRead)
library(Biostrings)
library(seqinr)
library(DECIPHER)

####Environment Setup####
theme_set(theme_bw())
setwd("~/projects/dada2_tutorial_16S")

####File Path Setup####
#this is so dada2 can quickly iterate through all the R1 and R2 files in your read set
path <- "~/projects/dada2_tutorial_16S/raw_data/" # CHANGE ME to the directory containing the fastq fi
les
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE)) #change the pattern to m
atch all your R1 files
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #change the delimiter in quotes and the
number at the end of this command to decide how to split up the file name, and which element to extrac
t for a unique sample name
#check your sample names
sample.names

####fastq Quality Plots####
#if you have hundreds of samples in your dataset, it's easiest to just look at a few (10-50) fastq fil
es
numplot <- 5 #CHANGE ME to the number of samples you want to include in the quality plot ( can be anyw
here from 2 to length(fnFs) )
a <- sample(fnFs, numplot) #randomly select a set of N samples
b <- which(fnFs %in% a) #identify the indices of those samples
plotfnFs <- fnFs[b] #use the indices to create two lists of corresponding forward and reverse files to
 plot
plotfnRs <- fnRs[b]
pdf("quality_plots.dada2.R1s.pdf", width = 16, height = 9) # define plot width and height. completely
up to user.
  plotQualityProfile(plotfnFs) #this plots the quality profiles for each sample
dev.off()

pdf("quality_plots.dada2.R2s.pdf", width = 16, height = 9) # define plot width and height. completely
up to user.
  plotQualityProfile(plotfnRs)
dev.off()

####Primer Removal####
####identify primers####
FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer sequence
REV <- "CCGYCAATTYMTTTRAGTTT"  ## CHANGE ME to your reverse primer sequence


allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, trimLeft = c(0,0), trimRight = c(0,0), truncLen=c(27
0,170), maxN = 0, multithread = 32, compress = TRUE, matchIDs=TRUE)


primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
index <- 6 #this is the index of the file we want to check for primers, within the lists "fn*s.filtN",
 it can be any number from 1 to N, where N is the number of samples you are processing
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[index]]), #the index of the
sample you'd like to use for this test is used here (your first sample may be a blank/control and not
have many sequences in it, be mindful of this)
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[index]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[index]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[index]]))

####OPTIONAL!!!!####
#REV <- REV.orients[["Reverse"]] #IMPORTANT!!! change orientation ONLY IF NECESSARY. see the dada2 ITS
_workflow guide section "Identify Primers" for details

#### primer removal ####
cutadapt <- "/usr/local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs.filtN))
fnRs.cut <- file.path(path.cut, basename(fnRs.filtN))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

#Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "-j", 32,# -n 2 required to remove FWD and R
EV from reads #-j sets no. threads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#sanity check, should report zero for all orientations and read sets
index <- 6 #this is the index of the file we want to check for primers, within the lists "fn*s.cut", i
t can be any number from 1 to N, where N is the number of samples you are processing
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[index]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[index]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[index]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[index]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1", full.names = TRUE)) #remember to change this so it
matches ALL your file names!
cutRs <- sort(list.files(path.cut, pattern = "R2", full.names = TRUE)) #remember to change this so it
matches ALL your file names!

####filter and trim reads####
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

####trim & filter####
#filter and trim command. dada2 can canonically handle lots of errors, I am typically permissive in th
e maxEE parameter set here, in order to retain the maximum number of reads possible. error correction
steps built into the dada2 pipeline have no trouble handling data with this many expected errors.
#it is best, after primer removal, to not truncate with 18s data, or with data from any region in whic
h the length is broadly variable. you may exclude organisms that have a shorter insert than the trunca
tion length (definitely possible, good example is giardia). defining a minimum sequence length is best
.
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(0,0), minLen = c(150,120),
                     maxN=c(0,0), maxEE=c(6,8), truncQ=c(2,2), rm.phix=TRUE, matchIDs=TRUE,
                     compress=TRUE, multithread=36)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
write.table(retained, "retained_reads.filterAndTrim_step.txt", sep="\t", row.names=TRUE, col.names=TRU
E, quote=FALSE)

####learn error rates####
#the next three sections (learn error rates, dereplication, sample inference) are the core of dada2's
sequence processing pipeline. read the dada2 paper and their online documentation (linked at top of th
is guide) for more information on how these steps work
errF <- learnErrors(filtFs, multithread=32)
errR <- learnErrors(filtRs, multithread=32)

pdf("error_rates.dada2.R1s.pdf", width = 10, height = 10) # define plot width and height. completely u
p to user.
  plotErrors(errF, nominalQ=TRUE) #assess this graph. it shows the error rates observed in your datase
t. strange or unexpected shapes in the plot should be considered before moving on.
dev.off()
pdf("error_rates.dada2.R2s.pdf", width = 10, height = 10) # define plot width and height. completely u
p to user.
  plotErrors(errR, nominalQ=TRUE) #assess this graph. it shows the error rates observed in your datase
t. strange or unexpected shapes in the plot should be considered before moving on.
dev.off()

#sometimes, samples with low sequence count don't produce any filtered output, so we have to do someth
ing extra here to get rid of their entries in filtFs and filtRs
filtFs <- sort(list.files("~/projects/dada2_tutorial_16S/raw_data/cutadapt/filtered/", pattern = "_R1_
001.fastq.gz", full.names = TRUE)) #remember to change the pattern so it matches ALL your file names!
filtRs <- sort(list.files("~/projects/dada2_tutorial_16S/raw_data/cutadapt/filtered/", pattern = "_R2_
001.fastq.gz", full.names = TRUE)) #remember to change the pattern so it matches ALL your file names!

####dereplication####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# rewrite names of derepFs derepRs so that they match what's in sample.names (just take the sampleID f
rom the first element, split on _)
names(derepFs) <- sapply(strsplit(names(derepFs), "_"), `[`, 1)
names(derepRs) <- sapply(strsplit(names(derepRs), "_"), `[`, 1)

####sample inference####
dadaFs <- dada(derepFs, err=errF, multithread=32)
dadaRs <- dada(derepRs, err=errR, multithread=32)

dadaFs[[1]]
dadaRs[[1]]

####OPTIONAL: remove low-sequence samples before merging####
#a "subscript out of bounds" error at the next step (merging) may indicate that you aren't merging any
 reads in one or more samples.
#NB, NOT getting this error doesn't necessarily mean that all of your samples end up with more than 0
merged reads, as i found out while processing a large 18s dataset. your guess is as good as mine as to
 why this error does or does not appear, but filtering out the samples that cause it is necessary for
completion of the pipeline.
#samples_to_keep <- as.numeric(out[,"reads.out"]) > 100 #example of simple method used above after the
 filter and trim step. if you already did this but still got an error when merging, try the steps belo
w
getN <- function(x) sum(getUniques(x)) #keeping track of read retention, number of unique sequences af
ter ASV inference
track <- cbind(sapply(derepFs, getN), sapply(derepRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN
))
samples_to_keep <- track[,4] > 1 #your threshold. try different ones to get the lowest one that will w
ork. #this method accounts for dereplication/ASVs left after inference
samples_to_remove <- names(samples_to_keep)[which(samples_to_keep == FALSE)] #record names of samples
you have the option of removing
write.table(names(which(samples_to_keep == TRUE)), "samples_retained.txt", row.names=FALSE, quote=F, s
ep="\n")
write.table(setdiff(sample.names, names(which(samples_to_keep == TRUE))), "samples_removed.txt", row.n
ames=FALSE, quote=F, sep="\n")

####merge paired reads####
#OPTION 2: modify command when removing low-sequence samples
mergers <- mergePairs(dadaFs[samples_to_keep], derepFs[samples_to_keep], dadaRs[samples_to_keep], dere
pRs[samples_to_keep], verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

####construct sequence table####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

####View Sequence Length Distribution Post-Merging####
#most useful with merged data. this plot will not show you much for forward reads only, which should h
ave a nearly uniform length distribution.
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab)))) #tabulate sequence length distri
bution
pdf("length_histogram.merged_reads.pdf", width = 10, height = 8) # define plot width and height. compl
etely up to user.
  plot(x=length.histogram[,1], y=length.histogram[,2],
       xlab = "length (bp)",
       ylab = "frequency (# observed ASVs)",
       main = "ASV Length Histogram") #view length distribution plot
dev.off()

####remove low-count singleton ASVs####
otus <- otu_table(t(seqtab), taxa_are_rows = TRUE)

#generate counts of sample per ASV
otu_pres_abs <- otus
otu_pres_abs[otu_pres_abs >= 1] <- 1 #creating a presence/absence table
otu_pres_abs_rowsums <- rowSums(otu_pres_abs) #counts of sample per ASV

#IF you want to filter out rare variants (low-read-count singleton ASVs) you can use phyloseq's "trans
form_sample_counts" to create a relative abundance table, and then filter your ASVs by choosing a thre
shold of relative abundance: otus_rel_ab = transform_sample_counts(otus, function(x) x/sum(x))
dim(seqtab) # sanity check
dim(otus) # (this should be the same as last command, but the dimensions reversed)

otus_rel_ab <- transform_sample_counts(otus, function(x) x/sum(x)) #create relative abundance table
df <- as.data.frame(unclass(otus_rel_ab)) #convert to plain data frame
df[is.na(df)] <- 0 #if there are samples with no merged reads in them, and they passed the merge step
(a possiblity, converting to a relative abundance table produes all NaNs for that sample. these need t
o be set to zero so we can do the calculations in the next steps.)
otus_rel_ab.rowsums <- rowSums(df) #compute row sums (sum of relative abundances per ASV. for those on
ly present in one sample, this is a value we can use to filter them for relative abundance on a per-sa
mple basis)
a <- which(as.data.frame(otu_pres_abs_rowsums) == 1) #which ASVs are only present in one sample
b <- which(otus_rel_ab.rowsums <= 0.001) #here is where you set your relative abundance threshold #whi
ch ASVs pass our filter for relative abundance
removed <- length(intersect(a,b)) #how many of our singleton ASVs fail on this filter
rows_to_remove <- intersect(a,b) #A also in B (we remove singleton ASVs that have a lower relative abu
ndance value than our threshold)
if (removed > 0) {
otus_filt <- otus[-rows_to_remove,] #filter OTU table we created earlier
} else {
otus_filt <- otus #nothing to filter
}
dim(otus_filt) #how many ASVs did you retain?
seqtab.nosingletons <- t(as.matrix(unclass(otus_filt))) #convert filtered OTU table back to a sequence
 table matrix to continue with dada2 pipeline

#Start ASV report
cat("dimensions of unfiltered ASV table:", dim(otus),file="ASV_report.txt",sep="\t",append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)
cat("# singleton ASVs Removed:", length(intersect(a,b)),file="ASV_report.txt",sep="\t",append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)
cat("dimensions of ASV table after singleton removal:", dim(otus_filt),file="ASV_report.txt",sep="\t",
append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)

####remove chimeras####
#here we remove "bimeras" or chimeras with two sources. look at "method" to decide which type of pooli
ng you'd like to use when judging each sequence as chimeric or non-chimeric
seqtab.nosingletons.nochim <- removeBimeraDenovo(seqtab.nosingletons, method="pooled", multithread=36,
 verbose=TRUE) #this step can take a few minutes to a few hours, depending on the size of your dataset
cat("dimensions of ASV table after chimera removal", dim(seqtab.nosingletons.nochim),file="ASV_report.
txt",sep="\t",append=TRUE)
cat("", file="ASV_report.txt", sep="\n", append=TRUE)
sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons) #proportion of nonchimeras #it should be rela
tively high after filtering out your singletons/low-count ASVs, even if you lose a lot of ASVs, the nu
mber of reads lost should be quite low
cat("proportion of non-chimeric to chimeric reads (0-1):", sum(seqtab.nosingletons.nochim)/sum(seqtab.
nosingletons),file="ASV_report.txt",sep="\t",append=TRUE)

####track read retention through steps####
#first remove samples in the "out" R object that aren't in samples_to_keep
row.names(out) <- as.character(t(as.data.frame(strsplit(row.names(out), "_")))[,1]) #make sure that th
e output you're assigning to the row.names of "out" matches your sample names before proceeding. selec
t the correct delimiter (_ in exampe) and element (1 in example) to extract a unique name for each sam
ple that corresponds to any metadata you will work with in downstream analyses
getN <- function(x) sum(getUniques(x))
track <- cbind(out[which(row.names(out) %in% names(samples_to_keep)[which(samples_to_keep == TRUE)]),]
, sapply(dadaFs[samples_to_keep], getN), sapply(dadaRs[samples_to_keep], getN), sapply(mergers, getN),
 rowSums(seqtab.nosingletons), rowSums(seqtab.nosingletons.nochim))
# If processing only a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with
getN(dadaFs)
track <- cbind(track, 100-track[,6]/track[,5]*100, 100-track[,7]/track[,6]*100, track[,7]/track[,1]*10
0)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nosingletons", "nochime
ras", "percent_singletons", "percent_chimeras", "percent_retained_of_total")

####save output from sequnce table construction steps####
write.table(data.frame("row_names"=rownames(track),track),"read_retention.16S_merged.txt", row.names=F
ALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(seqtab.nosingletons.nochim),seqtab.nosingletons.nochim),"s
equence_table.16S.merged.txt", row.names=FALSE, quote=F, sep="\t")


# #OPTIONAL: read in sequence table on cluster, for taxonomy assignment
# #code to read in before taxonomy assignment, if done on a separate machine
# seqtab.nosingletons.nochim <- fread("sequence_table.16S.merged.txt", sep="\t", header=T, colClasses
= c("row_names"="character"), data.table=FALSE)
# row.names(seqtab.nosingletons.nochim) <- seqtab.nosingletons.nochim[,1] #set row names
# seqtab.nosingletons.nochim <- seqtab.nosingletons.nochim[,-1] #remove column with row names in it
# seqtab.nosingletons.nochim <- as.matrix(seqtab.nosingletons.nochim) #cast the object as a matrix
# mode(seqtab.nosingletons.nochim) <- "numeric"


####taxonomy assignment with DECIPHER and GTDB####
### Taxonomy
## creating a DNAStringSet object from the ASVs
seq_16S <- DNAStringSet(getSequences(seqtab.nosingletons.nochim))

## downloading DECIPHER database #choose here between GTDB, SILVA, Contax for 16S. DECIPHER devs recom
mend GTDB. SILVA still has the broadest taxonomic content, but GTDB is the newest published reference
taxonomy. see here for details: http://www2.decipher.codes/ClassifyOrganismsFAQ.html
download.file("http://www2.decipher.codes/Classification/TrainingSets/GTDB_r214-mod_April2023.RData",
"GTDB_r214-mod_April2023.RData")
# loading ref taxonomy object
load("GTDB_r214-mod_April2023.RData")

#taxonomic inference #trying with lower threshold
tax_info_16S <- IdTaxa(seq_16S, trainingSet, strand = "both", threshold=40, processors = 38) #this is
faster than other methods we've used before. with 38 processors, it only took about 15 min


## making and writing out standard output files:
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs_16S <- colnames(seqtab.nosingletons.nochim)

asv_headers_16S <- vector(dim(seqtab.nosingletons.nochim)[2], mode = "character")
for (i in 1:dim(seqtab.nosingletons.nochim)[2]) {
  asv_headers_16S[i] <- paste(">ASV", i, sep = "")
}

# tax table:
    # creating vector of desired ranks
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
    # creating table of taxonomy and setting any that are unclassified as "NA"
tax_tab_16S <- t(sapply(tax_info_16S, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
}))


colnames(tax_tab_16S) <- ranks
row.names(tax_tab_16S) <- NULL
tax_tab_16S <- data.frame("ASV_ID" = sub(">", "", asv_headers_16S), tax_tab_16S, check.names = FALSE)

write.table(tax_tab_16S, "taxonomy_table.16S.DECIPHER.GTDB.threshold_40.txt", sep = "\t", quote = F, r
ow.names = FALSE)


#### save sequences and do taxonomy assignment with blast ####
##### replace the long ASV names (the actual sequences) with human-readable names####
#save the new names and sequences as a .fasta file in your project working directory, and save a table
 that shows the mapping of sequences to new ASV names
my_otu_table <- t(as.data.frame(seqtab.nosingletons.nochim)) #transposed (OTUs are rows) data frame. u
nclassing the otu_table() output avoids type/class errors later on
ASV.seq <- as.character(unclass(row.names(my_otu_table))) #store sequences in character vector
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') #create new names
write.table(cbind(ASV.num, ASV.seq), "sequence_ASVname_mapping.txt", sep="\t", quote=F, row.names=F, c
ol.names=F)
write.fasta(sequences=as.list(ASV.seq), names=ASV.num, "16S_ASV_sequences.fasta") #save sequences with
 new names in fasta format
#IMPORTANT: sanity checks
colnames(seqtab.nosingletons.nochim) == ASV.seq #only proceed if this tests as true for all elements

#rename your ASVs in the taxonomy table and sequence table objects
colnames(seqtab.nosingletons.nochim) <- ASV.num

#re-save sequence and taxonomy tables with updated names
write.table(data.frame("row_names"=rownames(seqtab.nosingletons.nochim),seqtab.nosingletons.nochim),"s
equence_table.16S.merged.w_ASV_names.txt", row.names=FALSE, quote=F, sep="\t")

####taxonomy assignment with SILVA v138 database (dada2 devs formatted)####
colnames(seqtab.nosingletons.nochim) <- ASV.seq #change column names back to sequences

taxa_boot <- assignTaxonomy(seqtab.nosingletons.nochim, "/data/taxonomyDBs/silva_for_dada2/v138_by_dad
a2team/silva_nr99_v138.1_train_set.fa.gz", multithread=38, outputBootstraps = TRUE)
taxa <- addSpecies(taxa_boot[[1]], "/data/taxonomyDBs/silva_for_dada2/v138_by_dada2team/silva_species_
assignment_v138.1.fa.gz") #here we operate on the first element of the list taxa_boot, the character m
atrix. the second element is a numeric matrix of confidence values up to genus equivalent taxonomic ra
nk

row.names(taxa) == ASV.seq #only proceed if this tests as true for all elements
row.names(taxa_boot$tax) == ASV.seq #only proceed if this tests as true for all elements
row.names(taxa_boot$boot) == ASV.seq #only proceed if this tests as true for all elements

#assign new ASV names
row.names(taxa) <- ASV.num
row.names(taxa_boot$boot) <- ASV.num

#re-save sequence and taxonomy tables with updated names
write.table(data.frame("row_names"=rownames(taxa),taxa),"taxonomy_table.16S.RDP_SILVA138_dada2.txt", r
ow.names=FALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(taxa_boot$boot),taxa_boot$boot),"taxonomy_table.16S.RDP_SI
LVA138_dada2.bootstraps.txt", row.names=FALSE, quote=F, sep="\t")


####taxonomy assignment with NCBI NT blast db and LCA####
#this is example code, you need a local copy of NCBI's NT database formatted for use with blast, which
 you can download from NCBI's FTP website
#depends on: blastn, galaxy-tool-BLAST (find on github), and galaxy-tool-LCA (find on github)
#now assign taxonomy with blast NT database at 96% similarity threshold using both 'LCA + besthit' and
 'LCA only' parameters
mkdir blast_96_sim_LCA_besthit
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 96 -qcov_hsp_pe
rc 50 -db /data/taxonomyDBs/NCBI_NT/2023-06-08/nt -outfmt '6 qseqid stitle sacc staxid pident qcovs ev
alue bitscore' -query 16S_ASV_sequences.fasta  -out blast_96_sim_LCA_besthit/16S_ASV_sequences.blast.o
ut
python2 ~/programs/galaxy-tool-BLAST/blastn_add_taxonomy_lite.py -i blast_96_sim_LCA_besthit/16S_ASV_s
equences.blast.out -t /data/taxonomyDBs/NCBI_taxonomy/2023-06-08/rankedlineage.dmp -m /data/taxonomyDB
s/NCBI_taxonomy/2023-06-08/merged.dmp -o blast_96_sim_LCA_besthit/taxonomy
cat <(head -n 1 ~/programs/galaxy-tool-lca/example/example.tabular) taxonomy_16S_ASV_sequences.blast.o
ut > tmp
python2 ~/programs/galaxy-tool-lca/lca.py -i tmp -o blast_96_sim_LCA_besthit/taxonomy_table.16S.NCBI_N
T.96sim.LCA+besthit.txt -b 100 -id 96 -cov 50 -t best_hit -tid 98 -tcov 80 -fh environmental,unidentif
ied,kingdom -flh unclassified

mkdir blast_96_sim_LCA_only
python2 ~/programs/galaxy-tool-lca/lca.py -i tmp -o blast_96_sim_LCA_only/taxonomy_table.16S.NCBI_NT.9
6sim.LCA_only.txt -b 100 -id 96 -cov 50 -t only_lca -fh environmental,unidentified -flh unclassified

#cleanup
rm blast_96_sim_LCA_besthit/16S_ASV_sequences.blast.out #remove blast output without taxonomy
rm taxonomy_16S_ASV_sequences.blast.out #remove redundant file
mv tmp blast_96_sim_LCA_besthit/16S_ASV_sequences.blast.out #replace with taxonomy added blast output






