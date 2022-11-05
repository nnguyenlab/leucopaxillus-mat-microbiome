setwd("/Volumes/GoogleDrive/Shared drives/Research/active_projects/mat_bacteria/Leucopaxillus_mat/analyses/R analyses/16S")

library(tidyverse)
library(vegan)
library(phyloseq)
library(multcompView)
library(patchwork)
library(pairwiseAdonis)
library(RColorBrewer)
#library(qiime2R)
#library(agricolae)
#display.brewer.all()
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
Faith_PD <- read.table(file="faith-PD.tsv", header=TRUE)

#Read in metadata file and merge Obseverd OTUs and Shannon Diversity with the metadata
metadata_table <- read.table(file="metadata-leucopax-bact.tsv", header=TRUE) %>% 
  left_join(Observed_otus, by="SampleID", row.names="SampleID") %>% 
  left_join(Shannon_div, by="SampleID") %>% 
  left_join(Faith_PD, by="SampleID") %>% 
  drop_na() %>% 
  as_tibble()

#Format data tables for phyloseq
otu <- otu_table(otu_table, taxa_are_rows=FALSE)
# taxa <- tax_table(taxa_table)
metadata <- sample_data(metadata_table)
tree<-read_tree("Leucopaxillus-bact-tree.nwk")

#Checking for consistent OTU names
#taxa_names(taxa)
#taxa_names(otu)

#Merge data tables into a phyloseq object
merged_tables <- phyloseq(otu, metadata, tree)#, taxa)


##===================================================================================================
##                                  Build alpha diversity Plots
##===================================================================================================

observed_asv_plot <- ggplot(metadata, aes(x=HostSampleType, y=observed_features, fill=HostSampleType)) +
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
  
  
faiths_pd_plot <- ggplot(metadata, aes(x=HostSampleType, y=faith_pd, fill=HostSampleType)) + 
  geom_boxplot(outlier.color="gray") +
  theme_classic() +
  scale_fill_manual(values=c("#B3DE68", "#8DD3C7", "#FDB462", "#FB8072"), name="Sample type", 
                    labels=c("L. albissimus mat", "L. albissimus non-mat", "L. gentianeus mat", "L. gentianeus non-mat")) +
  #scale_fill_brewer(palette="Dark2") +
  scale_x_discrete(limit=c("L.gentianeus_mat", "L.gentianeus_soil", "L.albissimus_mat", "L.albissimus_soil"), 
                   labels = c("L. gentianeus \n mat", "L. gentianeus \n non-mat", "L. albissimus \n mat", "L. albissimus \n non-mat")) +
  labs(x="", y = "Faith's Phylogenetic Diveristy", tag="B") +
  theme(axis.text.x = element_text(face = "italic")) +
  theme(plot.tag = element_text(size=18))

#Plot together
observed_asv_plot + faiths_pd_plot

#Save the plot
ggsave("Figure1-alphadiv.pdf", device="pdf", width=10.5, height=5)

##===================================================================================================
##                                  Testing ANOVA assumptions; perform ANOVA
##===================================================================================================
# Using the Shapiro Test -- if p-value not significant = data is normal = OK for ANOVA
shapiro.test(metadata$observed_features)#NS
shapiro.test(metadata$faith_pd)#NS
#shapiro.test(metadata$shannon_entropy)#S

# Test for homogeneity of variances -- if p-value not significant = data is homogeneous = OK for ANOVA
bartlett.test(metadata$observed_features~metadata$HostSampleType)#NS
bartlett.test(metadata$faith_pd~metadata$HostSampleType)#NS
#bartlett.test(metadata$shannon_entropy~metadata$HostSampleType)#S

# Will use Observed ASV and Faith's PD for manuscript; ANOVA assumptions OK
anova_observed_asv <- aov(metadata$observed_features ~ metadata$HostSampleType)
summary(anova_observed_asv)
tukey_observed_asv <- TukeyHSD(anova_observed_asv, ordered="TRUE")
tukey_observed_asv

anova_faiths_pd <- aov(metadata$faith_pd ~ metadata$HostSampleType)
summary(anova_faiths_pd)
tukey_faiths_pd <- TukeyHSD(anova_faiths_pd, ordered="TRUE")
tukey_faiths_pd

# Extract letters associated with significance with the "multcompView" package
observed_letters <- multcompLetters4(anova_observed_asv, tukey_observed_asv)
as.data.frame.list(observed_letters$`metadata$HostSampleType`)

faiths_letters <- multcompLetters4(anova_faiths_pd, tukey_faiths_pd)
as.data.frame.list(faiths_letters$`metadata$HostSampleType`)

##===================================================================================================
##                                  Build PCoA plots 
##===================================================================================================

#Perform the ordination
pcoa_unifrac <- ordinate(merged_tables, method="PCoA", distance="unifrac", weighted=FALSE)

#Plot PCoA
cpoa_plot <- plot_ordination(merged_tables, pcoa_unifrac, type="samples", color="SampleType", shape="HostSpecies")

#Beautify the plot
cpoa_plot +
  theme_bw() +
  geom_point(size=2, show.legend=FALSE) +
  stat_ellipse(level=0.7) +
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
metadata_table <- data.frame(sample_data(merged_tables))

#Bray-Curtis distance matrix
braydist <- phyloseq::distance(merged_tables, method="bray", weighted=TRUE)

#Unifrac distance matrix
unifracdist <- phyloseq::distance(merged_tables, method="unifrac", weighted=FALSE)

#Weighted Unifrac distance matrix
unifracdist_w <- phyloseq::distance(merged_tables, method="unifrac", weighted=TRUE)


##===================================================================================================
##                                Examining dispersion of data
##===================================================================================================
#Non-significant results from betadisper means that the groups are homogeneous and thus the results from adonis are trustworthy.
#Significant results from betadisper and significant results from adonis means that group differences could be due to within group variation, and interpretation from adonis is cautioned.

#Measure dispersion of Bray-Curtis distances
betadisper(braydist, metadata_table$HostSampleType) %>% 
  permutest(pair=FALSE)
#Non-Significant result from Bray-curtis distances

#Measure dispersion of unweighted Unifrac distances
betadisper(unifracdist, metadata_table$HostSampleType) %>% 
  permutest(pair=FALSE)
#Non-Significant result from unweighted matrix

#Measure dispersion of weighted Unifrac distances
betadisper(unifracdist_w, metadata_table$HostSampleType) %>% 
  permutest(pair=FALSE)
#Non-Significant result from weighted matrix

##===================================================================================================
##                                  PERMANOVA to compare groups
##===================================================================================================
#Run PERMANOVA with adonis2 with Bray-Curtis distance; not used in manuscript
# set.seed(2022)
# adonis2(braydist ~ SampleType*HostSpecies, data=metadata_table, permutations=999)
# 
# set.seed(2022)
# adonis2(braydist ~ HostSampleType, data=metadata_table, permutations=999)
# 
# #Blocked by species
# set.seed(2022)
# adonis2(braydist ~ SampleType, strata=metadata_table$HostSpecies, data=metadata_table, permutations=999)


#Run PERMANOVA with adonis2 with Unifrac distance; used in manuscript
set.seed(2022)
adonis2(unifracdist ~ SampleType, data=metadata_table, permutations=999)

set.seed(2022)
adonis2(unifracdist ~ SampleType*HostSpecies, data=metadata_table, permutations=999)

set.seed(2022)
adonis2(unifracdist ~ HostSampleType, data=metadata_table, permutations=999)

#Blocked by species
set.seed(2022)
adonis2(unifracdist ~ SampleType, strata=metadata_table$HostSpecies, data=metadata_table, permutations=999)


##===================================================================================================
##                       Performing a pair-wise PERMANOVA
##===================================================================================================
#Pair-wise adonis
set.seed(2002)
pairwise.adonis2(unifracdist ~HostSampleType, data = metadata_table, permutations=999)

# set.seed(2002)
# pairwise.adonis2(braydist ~ HostSampleType, data = metadata_table, permutations=999)