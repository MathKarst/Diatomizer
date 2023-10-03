##### Packages activation #####
library("dada2")
library("stringr")
##### Setting of the fastq path directory #####
path<-"A"

list.files(path)


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq.gz and SAMPLENAME_R2_001.fastq.gz


fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

##Alternative formats

fnFs <- sort(list.files(path, pattern="_L001_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2.fastq.gz", full.names = TRUE))

fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq.gz", full.names = TRUE))

fnFs <- sort(list.files(path, pattern=".R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern=".R2.fastq.gz", full.names = TRUE))



# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


##### Plot quality #####

plotQualityProfile(fnRs[1:8])
plotQualityProfile(fnFs[1:8])

#this plot may be useful to check the quality of your reads, the quality of the reverse read usually decline sooner than the forward quality  


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

##### Trimming & Filtering #####
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200),trimLeft =c(21,27),   ##c(27,22) for European(diat.barcode) Primers and c(21,27) for UK primers
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, #These argument values should work for most MiSeq Runs but can be changed , for more info check the dada2 website
                     compress=TRUE, multithread=16, verbose=TRUE) # On Windows set multithread=FALSE

#ratio of filtered seq
out<-as.data.frame(out)

out$ratio<-(out[,2]/out[,1])*100

head(out)
mean(out$ratio)

##### Errors rate learning #####

errF <- learnErrors(filtFs, multithread= T, verbose=TRUE, randomize = TRUE)

errR <- learnErrors(filtRs, multithread=F,  verbose=TRUE, randomize = TRUE)

plotErrors(errF, nominalQ=TRUE)

plotErrors(errR, nominalQ=TRUE)

##### Dereplication #####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

##### Denoising using the DADA2 algorythm ##### 
dadaFs <- dada(derepFs, err=errF, multithread=F,verbose=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=FALSE,verbose = TRUE)

dadaFs[[1]]

##### Merging paired reads  ####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


seqtab <- makeSequenceTable(mergers) # ASV table 
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

##### Chimera removal #####

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=F, verbose=TRUE) # If not using windows Multithreading could be used but might crash

dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

write.csv(seqtab.nochim,file=str_glue("ASV_table_nochime_{basename(path)}.csv"))# export ASV table cleaned of chimera

#write.csv(seqtab.nochim,file="ASV_table_nochim_INRA_Run2.csv")# export ASV table cleaned of chimera
#write.csv(seqtab.nochim,file="ASV_table_nochim_B2B6R_check.csv") # export ASV table cleaned of chimera
#write.csv(seqtab.nochim,file="ASV_table_nochim_INRA_both_runs.csv") # export ASV table cleaned of chimera
#write.csv(seqtab.nochim,file="ASV_table_nochime_RiverMeso.csv")

###### Tracking file #####

getN <- function(x) sum(getUniques(x))
track <- cbind(out[1:2], sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv2(track,file=str_glue("TrackFile_{basename(path)}.csv"))


##### Taxonomic assignment #####

## There we use 3 different classifiers : One created with the UK GoldStandard (corrected), one  from diat.barcode (rsyst) and one custom one based on diat.barcode with the addition of non diatoms taxa.



#taxa_GS_corrected <- assignTaxonomy(seqtab.nochim, "GoldStandard_UK_diatoms_rbcL_dada2.fasta",  taxLevels = c("Class","Genus", "Species","ID","clone"),outputBootstraps = FALSE, verbose = TRUE, multithread=FALSE, minBoot=60)
#write.csv(taxa_GS_corrected ,file=str_glue("taxa_{basename(path)}_GS_corrected_2019.csv"))

#taxa_GS_corrected_dash <- assignTaxonomy(seqtab.nochim, "corrected_barcodeUK_dada2_dash.fasta",  taxLevels = c("Class","Genus", "Species","ID","clone"),outputBootstraps = FALSE, verbose = TRUE, multithread=FALSE, minBoot=60)
#write.csv(taxa_GS_corrected_dash ,file=str_glue("taxa_{basename(path)}_GS_corrected_dash.csv"))

taxa_GS_corrected_2021_11 <- assignTaxonomy(seqtab.nochim, "corrected_barcodeUK_dada2_2021_11_15_1.fasta",  taxLevels = c("Class","Genus", "Species","ID","clone"),outputBootstraps = FALSE, verbose = F, multithread=FALSE, minBoot=60)
write.csv(taxa_GS_corrected_2021_11 ,file=str_glue("taxa_{basename(path)}_GS_corrected_2021_11.csv"))


#taxa_GS_2017_errors<- assignTaxonomy(seqtab.nochim, "FINAL_FERA2017_REFERENCE_POURRAVE.fasta",  taxLevels = c("Genus", "Species"),outputBootstraps = FALSE, verbose=TRUE,multithread = FALSE, minBoot=60)
#write.csv(taxa_GS_2017_errors ,file=str_glue("taxa_{basename(path)}_GS2017_errors.csv"))


taxa_diatbarcode <- assignTaxonomy(seqtab.nochim, "Rsyst__1401seqs_312bp_taxonomy_CLASSIFIER_DADA2.fasta", taxLevels = c("Domain", "Kingdom","infraKingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","Clone"), multithread=F,outputBootstraps = FALSE, verbose = TRUE, minBoot=60)
write.csv(taxa_diatbarcode ,file=str_glue("taxa_{basename(path)}_diatbarcode.csv"))


taxa_custom2019 <- assignTaxonomy(seqtab.nochim, "Rsyst__1401seqs_312bp_taxonomy_CLASSIFIER_DADA2_10feb2020_version.fasta", taxLevels = c("Domain", "Kingdom","infraKingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","Clone"),outputBootstraps = FALSE, verbose = TRUE, multithread=F, minBoot=60)    #On Windows set multithread=FALSE
write.csv(taxa_custom2019 ,file=str_glue("taxa_{basename(path)}_custom2019.csv"))


taxa_custom2022 <- assignTaxonomy(seqtab.nochim,"Diat.barcode_CLASSIFIER_DADA2_08_03_2022_version2.fasta", taxLevels = c("Domain", "Kingdom","infraKingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","Clone"),outputBootstraps = FALSE, verbose = TRUE, multithread=F, minBoot=60)
write.csv(taxa_custom2022 ,file=str_glue("taxa_{basename(path)}_custom2022_2.csv"))


bbtaxa.print <- taxa_XXXXX # Removing sequence rownames for display only, change taxa_XXXX by the name of the taxa file you want to display

rownames(taxa.print) <- NULL
head(taxa.print)

write.csv(taxa.print,file="taxa_{basename(path)}.csv")


##Step not related to the official pipeline but to the fragilaria paper with MARIA

SpeciesF2<-assignSpecies(seqtab.nochim, "Fragilaria_UK_fishing.fasta")
write.csv2(SpeciesF2 ,file=str_glue("AssignSPecies_{basename(path)}.csv")) 


taxaF <- assignTaxonomy(seqtab.nochim, "Fragilaria_UK.fasta", taxLevels=c("Species","cloneID"),outputBootstraps = FALSE, verbose = TRUE, multithread=FALSE)
write.csv(taxaF ,file=str_glue("taxa_{basename(path)}_Fragilaria.csv"))


##alternatively : Without chimera removing, facultative step in order to see the effect of chimera removal on the assignment  
#taxa1 <- assignTaxonomy(seqtab, "GoldStandard_UK_diatoms_rbcL_dada2.fasta", taxLevels = c("Domain", "Kingdom","infraKingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","Clone"), multithread=16)

#taxa.print1 <- taxa1 # Removing sequence rownames for display only
#rownames(taxa.print1) <- NULL
#head(taxa.print1)

##### Graphical representations using Phyloseq #####

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

#Palette 
library(ggsci)
library(RColorBrewer)
getPalette= colorRampPalette (c(pal_simpsons("springfield")(16),pal_futurama("planetexpress")(12),brewer.pal(9,"Set1") ,pal_simpsons("springfield")(16)))


##### Phyloseq object creation #####
 
#sample file made of sample names, facultative for most phyloseq object but needed for Krona plot later
s_data<-data.frame(row.names=sample.names,SampleID=sample.names,Pool="1") 

##chose one of them
taxtable<-taxa_GS_corrected_dash
taxtable<-taxa_diatbarcode
taxtable<-taxa_custom2021
basename(taxa_custom)

ps <- phyloseq(sample_data(s_data),otu_table(seqtab.nochim, taxa_are_rows=FALSE),tax_table(taxa_custom2022)) #change taxa with the taxonomic assignment you prefer

ps

#rarefying step (facultative)
ps.rare<-rarefy_even_depth(ps, sample.size = 6000)
write.csv(tax_table(ps.rare),str_glue("{basename(path)}fragi_rare_taxa.csv"))

# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(ps.rare), "matrix")
# transpose if necessary
if(taxa_are_rows(ps.rare)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

write.csv(OTU1,file=str_glue("{basename(path)}OTU1.csv"))
write.csv(OTUdf,file=str_glue("{basename(path)}OTUdf.csv"))



# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, title="Bray NMDS")

top50 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:50]

allnames <- names(sort(taxa_sums(ps), decreasing=TRUE))
sort(taxa_sums(ps))

ps.top50 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top50 <- prune_taxa(top50, ps.top50)

ps.all <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.all <- prune_taxa(allnames, ps.all)



## Assigned taxa : greyscale, Unassigned taxa : red 
plot_bar(ps.all,fill="Species")+scale_fill_grey()+scale_color_grey()+ 
  geom_bar(aes(color= Species,fill=Species), stat="identity", position="stack") +
  theme(legend.position ="bottom")## Species


ggsave(str_glue("{basename(path)}_Species_custom2022.jpeg"),width=60,height=60,units="cm")#export jpeg
ggsave(str_glue("{basename(path)}_Species_custom2022.pdf"),width=60,height=60,units="cm")#export pdf



plot_bar(ps.all,fill="Genus")+scale_fill_grey()+scale_color_grey()+ 
  geom_bar(aes(color= Genus,fill=Genus), stat="identity", position="stack") +
  theme(legend.position ="bottom")## Genus

ggsave(str_glue("{basename(path)}_Genus_custom2022.jpeg"),width=60,height=60,units="cm")#export jpeg
ggsave(str_glue("{basename(path)}_Genus_custom2022.pdf"),width=60,height=60,units="cm")#export pdf


plot_bar(ps.all,fill="Class")+scale_fill_grey()+scale_color_grey()+ 
  geom_bar(aes(color= Class,fill=Class), stat="identity", position="stack") +
  theme(legend.position ="bottom")## Class

ggsave(str_glue("{basename(path)}_Genus_custom2022.jpeg"),width=60,height=60,units="cm")#export jpeg
ggsave(str_glue("{basename(path)}_Genus_custom2022.pdf"),width=60,height=60,units="cm")#export pdf



##alternative palette
plot_bar(ps.all,fill = "Genus")+scale_fill_manual(values=getPalette(XXXXXXXX)) ## change the get_palette() value to the number of genus in your data
## you can change "Genus" for "Species" or "Phylum", or any other taxonomic rank present in your taxa file.

plot_bar(ps.top50,fill = "Genus")+scale_fill_manual(values=getPalette(40))

##### Miscellaneous #####
wh0 = genefilter_sample(ps, filterfun_sample(function(x) x > 5))
GP1 = prune_taxa(wh0, ps)
GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))

plot_bar(GP1,fill = "Genus")+scale_fill_manual(values=getPalette(80))

#####Species richness

phylorichness<-plot_richness(ps, measures=c("Shannon", "Simpson"))
write.csv(phylorichness$data,str_glue("{basename(path)}_richness.jpeg"))


###### KRONA graphs #####

library("psadd")

plot_krona(ps.all,output=str_glue("krona_{basename(path)}_taxa_custom2022_pool"),variable="Pool")


plot_krona(ps.all,output=str_glue("krona_{basename(path)}_taxacustom2022_all"),variable="SampleID")


