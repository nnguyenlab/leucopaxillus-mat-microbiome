#load libraries
library(phyloseq)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

setwd("/Volumes/GoogleDrive/Shared drives/Research/active_projects/mat_bacteria/Leucopaxillus_mat/analyses/R analyses")

#Color Brewer Pallette
#RColorBrewer::display.brewer.all()
#Manual color scales
PhylaColors <- c("Acidobacteriota" = "#1c9e77", "Actinobacteriota" = "#756fb4", "Bacteroidota" = "#f45c85", 
                 "Chloroflexi" = "#436fb6", "Cyanobacteria" = "seagreen", "Desulfobacterota" = "gold1", 
                 "Firmicutes" = "#d671aa", "Gemmatimonadota" = "#ed4131", "Myxococcota" = "#c2e871", 
                 "Patescibacteria" = "deepskyblue", "Planctomycetota" = "#7a81ff","Proteobacteria" = "orange", "Verrucomicrobiota" = "#66a621")

#Import a phyloseq object (biom should have taxonomy attached), start with UNRAREFIED otu table
leucopax_bact <- import_biom("table-final-with-taxonomy.biom")#, treefilename="Leucopaxillus-bact-tree.nwk", refseqfilename="seqs-leucopax-bact.fasta")

#Rename columns 
colnames(tax_table(leucopax_bact)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Save otu table (optional) -- this will save a simple OTU table without taxonomy
#write.csv(as.data.frame(otu_table(leucopax_bact)),"otu_leucopax_bact_initial.csv")

#Read in metadata; standard QIIME2 metadata table can be used here
metadata <- read.csv("metadata-leucopax-bact.csv")

#Create another row for "SampleID", create a dataframe for the metadata
row.names(metadata) <- metadata$X.SampleID
metadata$SampleID <- metadata$X.SampleID
sample_data(leucopax_bact) <- metadata
#View(data.frame(sample_data(leucopax_bact)))

#De-QIIME-ify the taxa table -- this will separate taxonomic ranks into each separate columns. This _0 table is necessary downstream.
tax_table_leucopax_bact_0 <- as.data.frame(tax_table(leucopax_bact))

#OPTIONAL export taxa table
# write.csv(tax_table_leucopax_bact_0, "tax_table_leucopax_bact_0.csv")

#Make a copy of the table
tax_table_leucopax_bact <- tax_table_leucopax_bact_0

#Renaming the taxonomy if not standard (this may not be necessary depending on the taxonomic database used)
#tax_table_leucopax_bact <- data.frame(lapply(tax_table_leucopax_bact, function(x) {gsub("Acidobacteriota", "Acidobacteria", x)}))
#tax_table_leucopax_bact <- data.frame(lapply(tax_table_leucopax_bact, function(x) {gsub("Actinobacteriota", "Actinobacteria", x)}))
#tax_table_leucopax_bact <- data.frame(lapply(tax_table_leucopax_bact, function(x) {gsub("Armatimonadota", "Armatimonadetes", x)}))
#tax_table_leucopax_bact <- data.frame(lapply(tax_table_leucopax_bact, function(x) {gsub("Bacteroidota", "Bacteroidetes", x)}))
#tax_table_leucopax_bact <- data.frame(lapply(tax_table_leucopax_bact, function(x) {gsub("Gemmatimonadota", "Gemmatimonadetes", x)}))
#tax_table_leucopax_bact <- data.frame(lapply(tax_table_leucopax_bact, function(x) {gsub("Halobacterota", "Halobacteria", x)}))
#tax_table_leucopax_bact <- data.frame(lapply(tax_table_leucopax_bact, function(x) {gsub("Planctomycetota", "Planctomycetes", x)}))
#tax_table_leucopax_bact <- data.frame(lapply(tax_table_leucopax_bact, function(x) {gsub("Verrucomicrobiota", "Verrucomicrobia", x)}))

#Remove taxonomic notations
tax_table_leucopax_bact$Kingdom<- gsub("d__", "", tax_table_leucopax_bact$Kingdom)#sometimes "k" is replaced by "d" so make sure that is is properly removed.
tax_table_leucopax_bact$Phylum <- gsub("p__", "", tax_table_leucopax_bact$Phylum)
tax_table_leucopax_bact$Class <- gsub("c__", "", tax_table_leucopax_bact$Class)
tax_table_leucopax_bact$Order <- gsub("o__", "", tax_table_leucopax_bact$Order)
tax_table_leucopax_bact$Family <- gsub("f__", "", tax_table_leucopax_bact$Family)
tax_table_leucopax_bact$Genus <- gsub("g__", "", tax_table_leucopax_bact$Genus)
tax_table_leucopax_bact$Species <- gsub("s__", "", tax_table_leucopax_bact$Species)
#tax_table_leucopax_bact$Gen_Fam <- paste(tax_table_leucopax_bact$Genus, " (", tax_table_leucopax_bact$Family,")",sep="")

#View(tax_table_leucopax_bact_0)
#View(tax_table_leucopax_bact)

row.names(tax_table_leucopax_bact) <- row.names(tax_table_leucopax_bact_0)
tax_table(leucopax_bact) <- as.matrix(tax_table_leucopax_bact)
# View(data.frame(tax_table(leucopax_bact)))

#subsetting your datasets (often it will requires a slow narrowing down of each category until you get the samples you want)
#can also use for subsetting taxa: leucopax_bact_sub0 = subset_taxa(leucopax_bact, Kingdom=="Bacteria")
#leucopax_bact_treatment = subset_samples(leucopax_bact, SampleType == "mat")#includes the following column category
#leucopax_bact_treatment = subset_samples(leucopax_bact_treatment, Treatment != "Control")#excludes the following column category
#leucopax_bact = subset_samples(leucopax_bact, SampleID != "02.2.5.1.16S.a")#exclude the following sample based on ID

##leucopax_bact <- leucopax_bact_treatment

# View(data.frame(otu_table(leucopax_bact_treatment)))

# Check rarefaction of the data
# rarecurve(t(otu_table(leucopax_bact)), step=50, cex=0.5)

#removing samples that didn't work -- low read counts (this would best be done in QIIME)
#leucopax_bact_0 <- prune_samples(sample_sums(leucopax_bact) >= 100, leucopax_bact)#if done here, pass this object onto downstream workflow instead of leucopax_bact
# #View(data.frame(sample_data(leucopax_bact_0)))
# #View(data.frame(otu_table(leucopax_bact_0)))
# rarecurve(t(otu_table(leucopax_bact_0)), step=100, cex=0.5)
# leucopax_bact_0_otu <- data.frame(otu_table(leucopax_bact_0))
# write.csv(leucopax_bact_0_otu,"leucopax_bact_0_otu.csv")

##########################################
###DESeq Comparison [soil (left) vs mat (right), both species together]
##########################################
###Note: This approach requires you to make a new phyloseq object for each comparison, which is shown below (using "subset_samples")

#Subsetting your dataset to make various comparisons
##leucopax_bact_sp <- subset_samples(leucopax_bact, Time.Point == "1")

#remove empty cells due to subsetting; empty cells can cause issues later
#can use this to filter out lower abundance taxa (e.g. x > 10)
leucopax_bact_filtered <- filter_taxa(leucopax_bact, function(x) sum(x) > 100, TRUE)

#An error may occur later on when running DESeq because all OTUs have one 0. Adding a pseudocount of +1 will solve this issue. If error does not occur, skip this step.
#The following code extracts 3 objects from the leucopax_bact phyloseq object, adds pseudocount +1 to the otu_table and then put everything back together into the original phyloseq object.
leucopax_bact_filtered <- phyloseq((otu_table(leucopax_bact_filtered)+1), sample_data(leucopax_bact_filtered), tax_table(leucopax_bact_filtered))

#Convert phyloseq to DESeq Data Set object (dds)
leucopax_bact_dds <- phyloseq_to_deseq2(leucopax_bact_filtered, ~SampleType)

#Determine which level should the dataset be set as the REFERENCE sample class
leucopax_bact_dds$SampleType <- relevel(leucopax_bact_dds$SampleType, "soil")

#Perform the contrast using standard and recognized parameters and tests
leucopax_bact_dds = DESeq(leucopax_bact_dds, test="Wald", fitType="parametric")

#Check to see if the comparision conditions are present
resultsNames(leucopax_bact_dds)

#Performing the final calculations and extracting the points
leucopax_bact_dds_results = results(leucopax_bact_dds, cooksCutoff = TRUE)

###Contrast reports results such that positive fold change means the first level is enriched with a specific taxa, and a negative fold change means the second level is enriched with a specific taxa. For instance, a positive log fold change in the results below indicates enrichment in "soil", and a negative fold change indicates enrichment in "mat".
#Reduce over estimation of fold changes in graphical format (makes the dots closer to better show on a figure)
#The last item of the contrast list should be the reference
#This dataset did not need shrinking
#leucopax_bact_dds_results = lfcShrink(dds=leucopax_bact_dds, contrast = c("SampleType","Soil","Mat"), res=leucopax_bact_dds_results, type="apeglm") #This did not work for me
#lfcShrink(dds = leucopax_bact_dds, coef = 2, type = "apeglm") #This worked

#choosing an alpha of 0.05 I feel is pretty conservative, especially because DESeq is already conservative, but I still typically go with it.
alpha = 0.05

#Extract information from the DESeq object. The objects designated as "sigtab" have a p-value < or equal to the alpha set above. The objects names "notsig" are the results that have p-values > the alpha. These latter results can provide insight into common OTUs/ASVs.
#finding the differential abundant for each ASVs -- between the treatments
sig_table = leucopax_bact_dds_results[which(leucopax_bact_dds_results$padj <= alpha), ]

#Bind taxa names to tables of significant taxa
sig_table = cbind(as(sig_table, "data.frame"), as(tax_table(leucopax_bact)[rownames(sig_table), ], "matrix"))

#View(sig_table_C1)
write.csv(sig_table, "sig_table.csv")

#Find ASVs that are not significant; ASVs that are common among the treatments
#notsig_table_C1 = leucopax_bact_C1dds_results[which(leucopax_bact_C1dds_results$padj > alpha), ]

#Bind taxa names to tables of not significant taxa
#notsig_table_C1 = cbind(as(notsig_table_C1, "data.frame"), as(tax_table(leucopax_bact_C1)[rownames(notsig_table_C1), ], "matrix"))
#head(notsig_table_C1)
# write.csv(notsig_table_C1, "notsig_table_C1.csv")

#This will provide a list of the phyla, so you can make sure all of your phyla are in the colors in "scale_colour_manual" below
sig_table_Phyla <- unique(sig_table$Phylum)
sig_table_Phyla

#Remove anything that do not have a family or phylum taxonomy
sig_table_sub <- subset(sig_table, Family!="N/A" & Phylum != "N/A" & Family!="" & Phylum !="" & Phylum !="RCP2-54")

#Rename items in the Genus column for formatting. Each dataset will vary but this list below is a standard list; add to it as needed.
sig_table_sub$Genus <- gsub("N/A", "unidentified bacterium", sig_table_sub$Genus)
sig_table_sub$Genus <- gsub("_", " ", sig_table_sub$Genus)
sig_table_sub$Genus <- gsub("uncultured", "unidentified bacterium", sig_table_sub$Genus)
sig_table_sub$Genus <- gsub("Burkholderia-Caballeronia-Paraburkholderia", "Para/Burkholderia-Caballeronia", sig_table_sub$Genus)
sig_table_sub$Genus <- gsub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Allo/Neo/Para-Rhizobium", sig_table_sub$Genus)

#Plot the logfold changes
sig_table_subp <- ggplot(sig_table_sub, aes(x=log2FoldChange, y=reorder(Genus,desc(Genus)), color=Phylum, alpha=0.9)) + 
  geom_point(size=2, stroke = 0.8) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.25, size="fit")) + 
  theme(legend.position="none") + 
  geom_vline(xintercept = 0, linetype = "solid", color = "black") + 
  ggtitle("Non-mat                                   Mat") + 
  theme(plot.title = element_text(hjust = 0.5, size=11))

#Plot and beautify by faceting based on phyla using manual colors
Diff_Abund_Plot <- sig_table_subp + 
  facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + 
  scale_y_discrete(position = "right") + 
  theme(strip.text.y.left = element_text(angle = 0)) + 
  theme(axis.text.y=element_text(size=6, face="italic")) + 
  theme(axis.title.y=element_blank()) + 
  theme(panel.spacing = unit(0.3, "lines")) + 
  theme(text = element_text(size = 8)) + 
  scale_colour_manual(values = PhylaColors) #+scale_color_brewer(palette="Dark2")
Diff_Abund_Plot

#Save a PDF of your file. Final edits can be made in Inkscape.
ggsave("Figure 3 - DESeq differential abundance.pdf", device="pdf", width=8, height=9.5)

