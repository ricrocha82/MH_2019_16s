install.packages("BiocManager")
install.packages("ggplot2")
install.packages("plyr")
install.packages("ape")
install.packages("picante")
install.packages("tidyverse")

library(BiocManager)

#BiocManager::install("phyloseq")
#BiocManager::install("DESeq2")

library(ggplot2)
library(phyloseq)
library(DESeq2)
library(plyr)
library(ape)
library(picante)
library(tidyverse)
library(rjson)




### Reproducing the categorize by function (level 3) functionality in plain-text tables.
### Doing this because adding a column of KEGG Pathways to a table and then converting
### that table to BIOM is difficult.

categorize_by_function_l3 <- function(in_ko, kegg_brite_mapping) {
  # Function to create identical output as categorize_by_function.py script,
  # but with R objects instead of BIOM objects in Python.
  # Input KO table is assumed to have rownames as KOs and sample names as columns.

  out_pathway <- data.frame(matrix(NA, nrow=0, ncol=(ncol(in_ko) + 1)))

  colnames(out_pathway) <- c("pathway", colnames(in_ko))

  for(ko in rownames(in_ko)) {
    
    # Skip KO if not in KEGG BRITE mapping df
    # (this occurs with newer KOs that weren't present in PICRUSt1).
    if(! ko %in% rownames(kegg_brite_mapping)) {
      next
    }
    
    pathway_list <- strsplit(kegg_brite_mapping[ko, "metadata_KEGG_Pathways"], "\\|")[[1]]
    
    for(pathway in pathway_list) {
      
      pathway <- strsplit(pathway, ";")[[1]][3]
      
      new_row <- data.frame(matrix(c(NA, as.numeric(in_ko[ko,])), nrow=1, ncol=ncol(out_pathway)))
      colnames(new_row) <- colnames(out_pathway)
      new_row$pathway <- pathway
      out_pathway <- rbind(out_pathway, new_row)
    }
    
  }
  
  out_pathway = data.frame(aggregate(. ~ pathway, data = out_pathway, FUN=sum))
  
  rownames(out_pathway) <- out_pathway$pathway
  
  out_pathway <- out_pathway[, -which(colnames(out_pathway) == "pathway")]
  
  if(length(which(rowSums(out_pathway) == 0)) > 0) {
    out_pathway <- out_pathway[-which(rowSums(out_pathway) == 0), ]
  }
  
  return(out_pathway)
  
}


### Example commands:
### Read in BRITE hierarchy per KO.
# kegg_brite_map <- read.table("/path/to/picrust1_KO_BRITE_map.tsv",
#                              header=TRUE, sep="\t", quote = "", stringsAsFactors = FALSE, comment.char="", row.names=1)
#
# 
### When reading in tab-delimited file of KO predictions (PICRUSt2 output):
# test_ko <- read.table("/path/to/test_ko.tsv", header=TRUE, sep="\t", row.names=1)
#
#
### Alternatively, when reading in legacy TSV BIOM file (PICRUSt1 output): 
### test_ko <- read.table("/path/to/test_ko.tsv",
###                       header=TRUE, sep="\t", row.names=1, skip=1, comment.char="")
### if(length(which(colnames(test_ko) == "KEGG_Pathways")) > 0)) {
###     test_ko <- test_ko[, -which(colnames(test_ko) == "KEGG_Pathways")]
### }
#
#
#
#
### Run function to categorize all KOs by level 3 in BRITE hierarchy.
# test_ko_L3 <- categorize_by_function_l3(test_ko, kegg_brite_map)
# test_ko_L3_sorted <- test_ko_L3[rownames(orig_ko_L3), ]
#
#
### Commands that could be used to compare the KO levels from this function with the actual output of categorize_by_function.py:
# orig_ko_L3 <- read.table("/path/to/test_ko_L3.tsv",
#                          header=TRUE, sep="\t", row.names=1, skip=1, comment.char="", quote="")
# 
# orig_ko_L3 <- orig_ko_L3[, -which(colnames(orig_ko_L3) == "KEGG_Pathways")]
# 
# orig_ko_L3 <- orig_ko_L3[-which(rowSums(orig_ko_L3) == 0),]
#
#
### The below command will be True when the output is exactly the same.
# identical(test_ko_L3_sorted, orig_ko_L3)



kegg_brite_map <- read.table("C:/Users/rae02e/Google Drive/FSP Environomics/Derwent trips/Derwent_16S_data/ASV Derwent/FSPEnvironomics_BaselineDataset/PICRust2/picrust1_KO_BRITE_map.tsv",
                             header=TRUE, sep="\t", quote = "", stringsAsFactors = FALSE,
                             comment.char="", row.names=1)

test_ko <- read.table("C:/Users/rae02e/Google Drive/FSP Environomics/Derwent trips/Derwent_16S_data/ASV Derwent/FSPEnvironomics_BaselineDataset/PICRust2/PICRUST2_Derw_KO_.tsv", header=TRUE, sep="\t", row.names=1)

test_ko_L3 <- categorize_by_function_l3(test_ko, kegg_brite_map)
#test_ko_L3_sorted <- test_ko_L3[rownames(orig_ko_L3), ]

write.csv(test_ko_L3, file="C:/Users/rae02e/Google Drive/FSP Environomics/Derwent trips/Derwent_16S_data/ASV Derwent/FSPEnvironomics_BaselineDataset/PICRust2/KEGG_functions_Derwent.csv")  


##############

library(BiocManager)
library(microbiome)
#BiocManager::install("microbiome")

#otu_mat <- read.csv("Derwent_B16S_ASV_16Scopynumbernormalised_July2020_ASV_nosequences.csv", header=TRUE)
#otu_mat <- otu_mat %>% select (-ASV)
#otu_mat<-round(otu_mat, digits = 0)


###MetaCyc
otu_fun_MetaCyc <- "C:/Users/rae02e/Google Drive/FSP Environomics/Derwent trips/Derwent_16S_data/ASV Derwent/FSPEnvironomics_BaselineDataset/PICRust2/Pathways_out_Derwent_path_abun_unstrat.csv"
tax_mat_MetaCyc <- "C:/Users/rae02e/Google Drive/FSP Environomics/Derwent trips/Derwent_16S_data/ASV Derwent/FSPEnvironomics_BaselineDataset/PICRust2/TAX_METACYC.csv"
samples_df <- "C:/Users/rae02e/Google Drive/FSP Environomics/Derwent trips/Derwent_16S_data/ASV Derwent/FSPEnvironomics_BaselineDataset/PICRust2/Meta_DERWENT_ASV.csv"
file.exists(otu_fun_MetaCyc)
file.exists(tax_mat_MetaCyc)
file.exists(samples_df)

carbom<-read_csv2phyloseq(
  otu.file = otu_fun_MetaCyc,
  taxonomy.file = tax_mat_MetaCyc,
  metadata.file = samples_df,
  sep = ","
)

###KEGG/KO
otu_fun_KO <- "C:/Users/rae02e/Google Drive/FSP Environomics/Derwent trips/Derwent_16S_data/ASV Derwent/FSPEnvironomics_BaselineDataset/PICRust2/Picrust2_Derwent_97OTU/KO_Derwent_pred_metagenome_unstrat.csv"
tax_mat_KO <- "C:/Users/rae02e/Google Drive/FSP Environomics/Derwent trips/Derwent_16S_data/ASV Derwent/FSPEnvironomics_BaselineDataset/PICRust2/TAX_KO.csv"
samples_df <- "C:/Users/rae02e/Google Drive/FSP Environomics/Derwent trips/Derwent_16S_data/ASV Derwent/FSPEnvironomics_BaselineDataset/PICRust2/Meta_DERWENT_ASV.csv"
file.exists(otu_fun_KO)
file.exists(tax_mat_KO)
file.exists(samples_df)

carbom<-read_csv2phyloseq(
  otu.file = otu_fun_KO,
  taxonomy.file = tax_mat_KO,
  metadata.file = samples_df,
  sep = ","
)

###KEGG/KO
#otu_fun_KEGG <- "C:/Users/rae02e/Google Drive/FSP Environomics/Derwent trips/Derwent_16S_data/ASV Derwent/FSPEnvironomics_BaselineDataset/PICRust2/Pathways_out_Derwent_path_abun_unstrat.csv"
#tax_mat_KEGG <- "C:/Users/rae02e/Google Drive/FSP Environomics/Derwent trips/Derwent_16S_data/ASV Derwent/FSPEnvironomics_BaselineDataset/PICRust2/TAX_METACYC.csv"
#samples_df <- "C:/Users/rae02e/Google Drive/FSP Environomics/Derwent trips/Derwent_16S_data/ASV Derwent/FSPEnvironomics_BaselineDataset/PICRust2/Meta_DERWENT_ASV.csv"
#file.exists(otu_fun_KEGG)
#file.exists(tax_mat_KEGG)
#file.exists(samples_df)

carbom<-read_csv2phyloseq(
  otu.file = otu_fun_KEGG,
  taxonomy.file = tax_mat_KEGG,
  metadata.file = samples_df,
  sep = ","
)


sample_sums(carbom)

# have a look at the read distribution
#Make a data frame with a column for the read counts of each sample##### 
sample_sum_df <- data.frame(sum = sample_sums(carbom))
# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 500000) +
  ggtitle("Distribution of sample sequencing depth - 16S ASV data") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

####### Rarefy
physeq1 = filter_taxa(carbom, function(x) sum(x) > 10, TRUE)
carbom_rar <-rarefy_even_depth(physeq1, sample.size = min(sample_sums(physeq1)),
                               rngseed = 123, replace = TRUE, trimOTUs = TRUE, verbose=TRUE)
sample_sums(carbom_rar)
taxa_names(carbom_rar)

###
library(devtools)
devtools::install_github("gauravsk/ranacapa")

ggrare <- function(physeq_object, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  # Get sample data
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  # Add, any custom-supplied plot-mapped variables
  if ( length(color) > 1 ) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  
  if ( length(label) > 1 ) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size",
                                           y = ".S",
                                           group = "Sample",
                                           color = color))
  
  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")
  
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x",
                                                    y = "y",
                                                    label = label,
                                                    color = color),
                                size = 4, hjust = 0)
  }
  
  p <- p + ggplot2::geom_line()
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}
p <- ggrare(carbom, step = 10, color = "Month_year", se = FALSE)


#########################
sample_data(carbom)$Site = factor(sample_data(carbom)$Site, levels=c("B1", "B3", "RBN", "G2","E", "KB","NTB5", "NTB13", "U2","U4", "U7", "U12"))
sample_data(carbom)$Estuary = factor(sample_data(carbom)$Estuary, levels=c("Lower-Estuary", "Mid-Estuary", "Zinc-Factory", "Upper-Estuary"))

my.cols <- c("darkorchid","orange", "blue1",
             "green" ,"black", 
             "cyan", "red", "gray",
             "goldenrod", "saddlebrown",
             "gold1", "slateblue", "chartreuse" )
carbom_hell<-transform(carbom, transform = "hellinger", target = "OTU")
#logt=transform_sample_counts(carbom_rar, function(x) sqrt(x))
out.pcoa.logt <- ordinate(carbom_hell, method = "CCA", distance = "bray", formula=~Site)
evals <- out.pcoa.logt$values$Eigenvalues
p3<-plot_ordination(carbom_hell, out.pcoa.logt, type = "samples", 
                    color = "Site")
p3 + scale_colour_manual(values=my.cols) + theme_bw()+ geom_point(size=3)+ggtitle("CCA Derwent MetaCyc") +
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=16, vjust = 0.3),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#,shape = "Month_year"
#########################

sample_data(carbom_rar)$Month_year <- factor(sample_data(carbom_rar)$Month_year, 
                          levels=c("May '18", "July '18", "August '18", "September '18", "October '18" ,
                                    "November '18","January '19","February '19",
                                   "March '19","April '19","May '19","June '19","July '19","August '19",
                                   "September '19","October '19","November '19",
                                   "December '19","February '20"))
                          
sample_data(carbom_rar)$Month = factor(sample_data(carbom_rar)$Month, levels=c("January","February","March","April","May","June","July","August","September","October","November","December"),
                                                                      labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

sample_data(carbom_rar)$Site = factor(sample_data(carbom_rar)$Site, levels=c("B1", "B3", "RBN", "G2","E", "KB","NTB5", "NTB13", "U2","U4", "U7", "U12"))
                                                                    
p<-plot_richness(carbom_rar,x="Site", color="Month",
                 measures=c("Chao1", "Shannon"))
p + theme_bw()+geom_point(shape = 1,size = 1.7,colour = "black")+
  scale_color_manual(values=my.cols)+ggtitle("Alpha diversity MetaCyc") +
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=16, vjust = 0.3),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #scale_x_continuous(name="Latitude", breaks=seq(-70, 5, 10))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13,face="bold"))+ theme(axis.text.x = element_text(angle = 90))+
  theme(panel.grid.major = element_blank())+stat_smooth(method = "lm", formula = y ~ x, se = FALSE)

####Shannon only
p<-plot_richness(carbom_rar,x="Latitude..decimal.degrees.", color="Water.mass",
                 measures=c("Chao1"))
p + theme_bw()+geom_point(shape = 1,size = 1.7,colour = "black")+
  scale_color_manual(values=c("red", "blue", "darkgreen", "orange"))+ggtitle("Alpha diveristy_3depths_MetaCyc_97_20k") +
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=16, vjust = 0.3),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(name="Latitude", breaks=seq(-70, 5, 10))+
  #scale_y_continuous(name="Alpha Diversity Measure", limits=c(5, 5.5),breaks=seq(5,5.5,0.5))+
  scale_y_continuous(name="Alpha Diversity Measure", limits=c(280, 380),breaks=seq(280,380,40))+
  #scale_y_continuous(name="Alpha Diversity Measure", limits=c(5, 5.5),breaks=seq(5,5.5,0.2))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13,face="bold"))+
  theme(panel.grid.major = element_blank())+stat_smooth(method = "lm", formula = y ~ x, se = FALSE)

# Calculate richness
plot_richness(carbom_rar,x="Latitude..decimal.degrees.", color="Water.mass",
              measures=c("Chao1"))
results = estimate_richness(carbom_rar, measures = 'Chao1')
d = sample_data(carbom_rar)

write.csv(results, "MetaCyc_Richness_results_97_MLD_rar_20613.csv")
Richness<-read.csv("MetaCyc_Richness_results_97_MLD_rar_20613.csv", header=TRUE)

compare_means(Chao1 ~ Water.mass,  data = Richness, p.adjust.method = "bonferroni")
library(dplyr)

ddply(Richness, ~Water.mass, plyr:::summarise, mean = mean(Chao1), sd = sd(Chao1),
      max = max(Chao1), min = min(Chao1))


#######
install.packages("ggpubr")
library(ggpubr)
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")

carbom_2 = transform_sample_counts(carbom_rar, function(x) x / sum(x) )

carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Amino-Acid-Biosynthesis")) #1
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Nucleotide-Biosynthesis")) #2
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Energy-Metabolism"))       #3
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Lipid-Biosynthesis"))      #4
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Carbohydrates-Biosynthesis")) #5
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Cell-Structure/ Cell-Wall-Biosynthesis")) #6
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Cofactor-Biosynthesis"))   #7
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Secondary-Metabolite-Biosynthesis")) #8
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Vitamin-Biosynthesis"))     #9
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Fermentation"))       #10

carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Sulfur-Metabolism")) #IndiSpecies
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Aromatic-compound-degradation"))
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Carbohydrates-Degradation"))
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Sugars-And-Acids-Degradation"))
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Amino-Acid-Degradation"))
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Sugars-And-Polysaccharides-Degradation"))
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Denitrification"))
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Alcohol-Degradation "))
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Phosphorus-Compounds"))
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("CO2-Fixation"))
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Siderophores-Biosynthesis"))
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Methanogenesis"))
carbom_1_CO2 = subset_taxa(carbom_2, Pathway==c("Antibiotic-Biosynthesis"))

#carbom_1_CO2 = subset_taxa(carbom_2, Ontology...parents.of.class==c("Antibiotic-Resistance "))

glom <- tax_glom(carbom_1_CO2, taxrank = 'Pathway')
#glom1 <- tax_glom(carbom_1_CO2, taxrank = 'Ontology...parents.of.class')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$Pathway <- as.character(data_glom$Pathway) #convert to character
#data_glom$Water.mass <- factor(data_glom$Water.mass, levels = c("Southern Ocean","STF","SPSG","PED"))
#p<-ggplot(data_glom, aes(x=Latitude..decimal.degrees. , y=Abundance, fill=Water.mass))
abundance_percent<-data_glom$Abundance*100
p<-ggplot(data_glom, aes(x=Site , y=abundance_percent, fill=Month))
p <- p +geom_point(shape = 21,size = 2,colour = "black")+
  scale_fill_manual(values=my.cols) +
  #p <- p + geom_smooth(method = "lm", se = FALSE, 
  #                     formula = my.formula, 
  #                     colour = "red")+
  ggtitle("Methanogenesis") +
  #scale_x_continuous(name="Latitude", limits=c(-70, 10), breaks=seq(-70,5,10)) +
  scale_y_continuous(name="Relative Abundance (%)", labels = scales::number_format(accuracy = 0.001, decimal.mark = '.'))+
  #limits=c(0,0.0002), breaks=seq(0,0.02,0.0005))+
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=16, vjust = 0.3),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p 
#+geom_smooth(method = "lm")
compare_means(Abundance ~ Water.mass,  data = data_glom, p.adjust.method = "bonferroni")
library(dplyr)
ddply(data_glom, ~Water.mass, plyr:::summarise, mean = mean(Abundance), sd = sd(Abundance),
      max = max(Abundance), min = min(Abundance))

#CO2:  stat_compare_means(comparisons = my_comparisons,label.y = c(0.0159,0.017, 0.0151, 0.0145)
#Fermentation:   stat_compare_means(comparisons = my_comparisons,label.y = c(0.0263,0.0215,0.0205)
#Sulfur-Metabolism:   stat_compare_means(comparisons = my_comparisons,label.y = c(0.00215,0.0013,0.0015, 0.0011)
#Lipid-Biosynthesis:  stat_compare_means(comparisons = my_comparisons,label.y = c(0.10,0.097,0.099, 0.094)
#Polyamine-Biosynthesis: stat_compare_means(comparisons = my_comparisons,label.y = c(0.0055,0.0094, 0.01)
#Carbohydrates-Biosynthesis:   stat_compare_means(comparisons = my_comparisons,label.y = c(0.051,0.0517, 0.048)
#Degradation_SFC:   stat_compare_means(comparisons = my_comparisons,label.y = c(0.0285,0.026, 0.0268,0.024)
#Cofactor-Biosynthesis_SFC:    stat_compare_means(comparisons = my_comparisons,label.y = c(0.124,0.139, 0.137)
#Cell-Structure-Biosynthesis_SFC:    stat_compare_means(comparisons = my_comparisons,label.y = c(0.0155,0.019, 0.0195)
#Cell-Wall-Biosynthesis:   stat_compare_means(comparisons = my_comparisons,label.y = c(0.0129,0.0143, 0.0146)
#Nucleotide-Biosynthesis:  stat_compare_means(comparisons = my_comparisons,label.y = c(0.176,0.177, 0.174)

######

################################
#Deseq2 https://github.com/Microbial-Ecology/DESeq2-for-Metagenome/blob/master/DESeq2_Metagenome_Script.R
##########################
library(DESeq2)
#cut alpha level at 0.05. We might want to change this later to higher level of significance if we get heaps of otu responding

library(DESeq2)
#cut alpha level at 0.05. We might want to change this later to higher level of significance if we get heaps of otu responding

alpha = 0.05

#alpha can be changed to 0.005  or 0.01 if necessary

alpha = 0.005
alpha = 0.01
alpha = 0.001

#subset first otherwise you get NA error!!!!
Lower_Zinc = subset_samples(carbom, Estuary == "Lower-Estuary" | Estuary =="Zinc-Factory")
Lower_Zinc <- tax_glom(Lower_Zinc, taxrank="Pathway")

#set up model Lithifying vs Nonlithifying
#sets it up so that the first factor will be on the bottom, so > 0 means it is higher in Nonlithifying < 0 means higher in Lithifying
sample_data(Lower_Zinc)$Estuary <- factor(sample_data(Lower_Zinc)$Estuary, levels = c("Lower-Estuary", "Zinc-Factory"))

#1st option look at lith Vs nonlith -> don't use, use option 2 which variance stabilses the data
#dds<-phyloseq_to_deseq2(Lower_Zinc,~Estuary)
#<-DESeq(dds,test="Wald",fitType = "parametric")
#com<-results(dds)


#####2nd option https://github.com/joey711/phyloseq/issues/283 # GetVarianceStabilizedData
diagdds<-phyloseq_to_deseq2(Lower_Zinc,~Estuary)
diagdds = estimateSizeFactors(diagdds)
diagdds = estimateDispersions(diagdds)
diagvst = getVarianceStabilizedData(diagdds)
dim(diagvst)

# Save the untransformed data as a separate variable so you can go back to
# it
Lower_Zinc0 = Lower_Zinc
otu_table(Lower_Zinc) <- otu_table(diagvst, taxa_are_rows = TRUE)

dds<-DESeq(diagdds,test="Wald",fitType = "parametric")
com<-results(dds)

#pull out model parameters for plotting
sigtab = com[which(com$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Lower_Zinc)[rownames(sigtab), ], "matrix"))

test=as.data.frame(com)
#we can write out the above as a csv if we want to join and save it with other data later on.

# order the graph (Genus). You can change this to other taxonomic levels too
x = tapply(sigtab$log2FoldChange, sigtab$Pathway, function(x) max(x))
x = sort(x, TRUE)
sigtab$Pathway = factor(as.character(sigtab$Pathway), levels=names(x))

#graph
DA<-ggplot(sigtab, aes(y=Pathway, x=log2FoldChange)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 1)+geom_point(aes(fill=Pathway), size=3.5,color="black",shape=21) + 
  theme(legend.key = element_rect(fill="white"), panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"), panel.background = element_rect(fill="white",colour="grey50"), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

DA<-DA+theme(panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank())
DA


##########