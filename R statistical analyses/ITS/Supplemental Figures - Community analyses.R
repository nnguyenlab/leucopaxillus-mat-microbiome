setwd("/Volumes/GoogleDrive/Shared drives/Research/active_projects/mat_bacteria/Leucopaxillus_mat/analyses/R analyses/ITS")

library(tidyverse)
library(vegan)
library(phyloseq)
library(patchwork)
library(pairwiseAdonis)
library(RColorBrewer)
#display.brewer.all()
#library(agricolae)
#library(qiime2R)

##===================================================================================================
##                                  Read in and wrangle data 
##===================================================================================================

#Read in rarefied OTU table
otu_table <- read.table(file="rarefied-table-with-taxonomy.tsv", header=TRUE, sep='\t', check.names=FALSE, row.names=1) %>% 
  select(-taxonomy) %>% 
  t() %>%
  as_tibble()

# #Read in taxonomy
# taxa_table <- read.table(file="taxonomy.tsv", header=TRUE, sep='\t', check.names=FALSE) %>% 
#   rename(Feature.ID="Feature ID") %>% 
#   parse_taxonomy() %>% 
#   as.matrix()

#Read in diversity metrics
Observed_otus <- read.table(file="observed-otus.tsv", header=TRUE)
Shannon_div <- read.table(file="shannon-diversity.tsv", header=TRUE)
#Faith_PD <- read.table(file="faith-PD.tsv", header=TRUE)

#Read in metadata file and merge Obseverd OTUs and Shannon Diversity with the metadata
metadata_table <- read.table(file="metadata-leucopax.txt", header=TRUE) %>% 
  left_join(Observed_otus, by="SampleID", row.names="SampleID") %>% 
  left_join(Shannon_div, by="SampleID") %>% 
  #left_join(Faith_PD, by="SampleID") %>% 
  drop_na() %>% 
  as_tibble()

#Format data tables for phyloseq
otu <- otu_table(otu_table, taxa_are_rows=FALSE)
# taxa <- tax_table(taxa_table)
metadata <- sample_data(metadata_table)
# tree<-read_tree("Leucopaxillus-bact-tree.nwk")

#Checking for consistent OTU names
#taxa_names(taxa)
#taxa_names(otu)

#Merge data tables into a phyloseq object
merged_tables <- phyloseq(otu, metadata)#, tree, taxa)


##===================================================================================================
##                                  Build alpha diversity Plots
##===================================================================================================

observed_asv_plot_fungi <- ggplot(metadata, aes(x=HostSampleType, y=observed_features, fill=HostSampleType)) +
  #stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_boxplot(outlier.color="gray", show.legend = FALSE) +
  theme_classic() +
  scale_fill_manual(values=c("#B3DE68", "#8DD3C7", "#FDB462", "#FB8072"), name="Sample type", 
                    labels=c("L. albissimus mat", "L. albissimus non-mat", "L. gentianeus mat", "L. gentianeus non-mat")) +
  #scale_fill_brewer(palette="Dark2") +
  scale_x_discrete(limit=c("L.gentianeus_mat", "L.gentianeus_soil", "L.albissimus_mat", "L.albissimus_soil"), 
                   labels = c("L. gentianeus \n mat", "L. gentianeus \n non-mat", "L. albissimus \n mat", "L. albissimus \n non-mat")) +
  labs(x="", y = "Observed ASVs", tag="A") +
  theme(axis.text.x = element_text(face = "italic")) +
  theme(plot.tag = element_text(size=18))
  
  
faiths_pd_plot_fungi <- ggplot(metadata, aes(x=HostSampleType, y=shannon_entropy, fill=HostSampleType)) + 
  geom_boxplot(outlier.color="gray") +
  theme_classic() +
  scale_fill_manual(values=c("#B3DE68", "#8DD3C7", "#FDB462", "#FB8072"), name="Sample type", 
                    labels=c("L. albissimus mat", "L. albissimus non-mat", "L. gentianeus mat", "L. gentianeus non-mat")) +
  #scale_fill_brewer(palette="Dark2") +
  scale_x_discrete(limit=c("L.gentianeus_mat", "L.gentianeus_soil", "L.albissimus_mat", "L.albissimus_soil"), 
                   labels = c("L. gentianeus \n mat", "L. gentianeus \n non-mat", "L. albissimus \n mat", "L. albissimus \n non-mat")) +
  labs(x="", y = "Shannon's Diveristy", tag="B") +
  theme(axis.text.x = element_text(face = "italic")) +
  theme(plot.tag = element_text(size=18))

#Plot together
observed_asv_plot_fungi + faiths_pd_plot_fungi

#Save the plot
ggsave("FigureS1-alphadiv.pdf", device="pdf", width=10.5, height=5)

##===================================================================================================
##                                  Testing ANOVA assumptions; perform comparisons
##===================================================================================================
# Check for normality of data using the Shapiro-Wilk test -- if p-value not significant = data is normal => OK for ANOVA
shapiro.test(metadata$observed_features)#S
#shapiro.test(metadata$faith_pd)
shapiro.test(metadata$shannon_entropy)#S

# Test for homogeneity of variances using the Bartlett test -- if p-value not significant = data is homogeneous => OK for ANOVA
bartlett.test(metadata$observed_features~metadata$HostSampleType)#S
#bartlett.test(metadata$faith_pd~metadata$HostSampleType)
bartlett.test(metadata$shannon_entropy~metadata$HostSampleType)#S

# ANOVA assumptions failed; use the non-parametric Kruskal-Wallis test
kruskal_asv <- kruskal.test(observed_features ~ SampleType, data = metadata_table)
kruskal_shannon <- kruskal.test(shannon_entropy ~ SampleType, data = metadata_table)

# perform test for HostSampleType (group by fungal species)
kruskal_asv <- kruskal.test(observed_features ~ HostSampleType, data = metadata_table)
kruskal_shannon <- kruskal.test(shannon_entropy ~ HostSampleType, data = metadata_table)

# Pairwise Kruskal-Wallis comparisons
wilcox_asv <- pairwise.wilcox.test(metadata_table$observed_features, metadata_table$HostSampleType, p.adjust.method="BH", paired=FALSE)
wilcox_shannon <- pairwise.wilcox.test(metadata_table$shannon_entropy, metadata_table$HostSampleType, p.adjust.method="BH", paired=FALSE)

##===================================================================================================
##                                  Build PCoA plots 
##===================================================================================================

#Perform the ordination
pcoa_bray <- ordinate(merged_tables, method="PCoA", distance="bray", weighted=TRUE)

#Plot PCoA
cpoa_plot <- plot_ordination(merged_tables, pcoa_bray, type="samples", color="SampleType", shape="HostSpecies")

#Beautify the plot
cpoa_plot +
  theme_bw() +
  geom_point(size=2, show.legend=FALSE) +
  stat_ellipse() +
  #ggtitle("Fungi at depth") +
  theme(text=element_text(size=10)) +
  scale_color_manual(values = c("mat"="darkorange", "soil" = "#7a81ff"), name="Sample type", labels=c("Mat", "Non-mat")) +
  scale_shape_manual(values = c(16, 8), name="Host fungus", labels=c("L. albissimus", "L. gentianeus"))

#Save the plot
ggsave("Figure 2 - ordination.pdf", device="pdf", width=4.3, height=3)


##===================================================================================================
##                                Making distance matrices
##===================================================================================================
#Extract the metadata table from phyloseq object
#metadata_table <- data.frame(sample_data(merged_tables))

#Bray-Curtis distance matrix
braydist <- phyloseq::distance(merged_tables, method="bray", weighted=TRUE)

##===================================================================================================
##                                Examining dispersion of data
##===================================================================================================
#Non-significant results from betadisper means that the groups are homogeneous and thus the results from adonis are trustworthy.
#Significant results from betadisper and significant results from adonis means that group differences could be due to within group variation, and interpretation from adonis is cautioned.

#Measure dispersion of Bray-Curtis distances
betadisper(braydist, metadata_table$HostSampleType) %>% 
  permutest(pair=TRUE)
#Significant result -- meaning that interpretation should be cautioned.

##===================================================================================================
##                                  PERMANOVA to compare groups
##===================================================================================================
#Run PERMANOVA with adonis2 with Bray-Curtis distance
set.seed(2022)
adonis2(braydist ~ SampleType*HostSpecies, data=metadata_table, permutations=999)

set.seed(2022)
adonis2(braydist ~ HostSampleType, data=metadata_table, permutations=999)

#Blocked by species
set.seed(2022)
adonis2(braydist ~ SampleType, strata=metadata_table$HostSpecies, data=metadata_table, permutations=999)

#All results are significant but the samples clustered very strongly so we will go ahead and use this data anyway.

##===================================================================================================
##                       Performing a pair-wise PERMANOVA
##===================================================================================================
#Pair-wise adonis
set.seed(2002)
pairwise.adonis2(braydist ~ HostSampleType, data = metadata_table, permutations=999)