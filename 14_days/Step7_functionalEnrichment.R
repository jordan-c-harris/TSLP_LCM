# Introduction to this script ----
# for the purposes of this script we'll want several data objects generated in previous scripts, including:
# 1) your normalized filtered expression data, in the form of a data matrix with symbols as rownames.
# 2) your study design file
# 3) your contrast matrix that lays out the pairwise comparisons you're interested in testing
# 4) Individual signatures or 'collections' of signatures to test for enrichment in your data.
# These signatures can be downloaded from gene signature databases such as MSigDB
# Signatures can also be custom made based on your interests.
# Signatures can also be pulled from R/Bioconductor as described below

# Load packages ----
library(tidyverse)
library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots

# Carry out GO enrichment using gProfiler2 ----
# use topTable result to pick the top genes for carrying out a Gene Ontology (GO) enrichment analysis
myTopHits <- topTable(ebFit, adjust.method ="BH", coef=1, number=50, sort.by="logFC")
# use the 'gost' function from the gprofiler2 package to run GO enrichment analysis
gost.res <- gost(rownames(myTopHits), organism = "mmusculus", correction_method = "fdr")
# produce an interactive manhattan plot of enriched GO terms
gostplot(gost.res, interactive = T, capped = T) #set interactive=FALSE to get plot for publications
# produce a publication quality static manhattan plot with specific GO terms highlighted
# rerun the above gostplot function with 'interactive=F' and save to an object 'mygostplot'
publish_gostplot(
  mygostplot, #your static gostplot from above
  highlight_terms = c("GO:0034987"),
  filename = NULL,
  width = NA,
  height = NA)

#you can also generate a table of your gost results
publish_gosttable(
  gost.res,
  highlight_terms = NULL,
  use_colors = TRUE,
  show_columns = c("source", "term_name", "term_size", "intersection_size"),
  filename = NULL,
  ggplot=TRUE)
# now repeat the above steps using only genes from a single module from the step 6 script, by using `rownames(myModule)`
# what is value in breaking up DEGs into modules for functional enrichment analysis?

# Perform GSEA using clusterProfiler ----
# there are a few ways to get msigDB collections into R
# option1: download directly from msigdb and load from your computer
# can use the 'read.gmt' function from clusterProfiler package to create a dataframe,
# alternatively, you can read in using 'getGmt' function from GSEABase package if you need a GeneSetCollection object
c2cp <- read.gmt("/Users/danielbeiting/Dropbox/MSigDB/c2.cp.v7.1.symbols.gmt")

# option2: use the msigdb package to access up-to-date collections
# this option has the additional advantage of providing access to species-specific collections
# are also retrieved as tibbles
msigdbr_species()
mm_gsea <- msigdbr(species = "Mus musculus") #gets all collections/signatures with human gene IDs
#take a look at the categories and subcategories of signatures available to you
mm_gsea %>%
  dplyr::distinct(gs_cat, gs_subcat) %>%
  dplyr::arrange(gs_cat, gs_subcat)

# choose a specific msigdb collection/subcollection
# since msigdbr returns a tibble, we'll use dplyr to do a bit of wrangling
mm_gsea_c2 <- msigdbr(species = "Mus musculus", # change depending on species your data came from
                      category = "C2") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature

# Now that you have your msigdb collections ready, prepare your data
# grab the dataframe you made in step3 script
# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC)
# construct a named vector
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=mm_gsea_c2, verbose=FALSE)
myGSEA.df <- as_tibble(myGSEA.res@result)

# view results as an interactive table
datatable(myGSEA.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Signatures enriched with TSLP AAV treatment',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:10), digits=2)
# create enrichment plots using the enrichplot package
gsea_inflam <- gseaplot2(myGSEA.res,
          geneSetID = c(919), #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          title = "Inflammatory gene sets enriched with TSLP AAV treatment") #can also turn off this title

## Initiate writing to PDF file
pdf("gsea_inflammatory_gene_sets.pdf", height = 5, width = 8)

## Create a graphical object g here
gsea_inflam

## Stop writing to the PDF file
dev.off()

gsea_lipid <- gseaplot2(myGSEA.res,
                         geneSetID = c(692), #can choose multiple signatures to overlay in this plot
                         pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
                         title = "Fatty acid metabolism gene set decreased with 14 day TSLP AAV treatment") #can also turn off this title

## Initiate writing to PDF file
pdf("gsea_lipid_gene_sets.pdf", height = 5, width = 9)

## Create a graphical object g here
gsea_lipid

## Stop writing to the PDF file
dev.off()

gsea_death <- gseaplot2(myGSEA.res,
                        geneSetID = c(635), #can choose multiple signatures to overlay in this plot
                        pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
                        title = "Programmed cell death gene set enriched with 14 day TSLP AAV treatment") #can also turn off this title

## Initiate writing to PDF file
pdf("gsea_celldeath_gene_set.pdf", height = 5, width = 8)

## Create a graphical object g here
gsea_death

## Stop writing to the PDF file
dev.off()

# add a variable to this result that matches enrichment direction with phenotype
myGSEA.df <- myGSEA.df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "TSLP",
    NES < 0 ~ "Control"))
myGSEA.df_desc <- myGSEA.df[order(-myGSEA.df$NES),]

# create 'bubble plot' to summarize y signatures across x phenotypes
TSLP_bubble <- ggplot(myGSEA.df_desc[1:20,], aes(x=phenotype, y=ID)) +
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
 # scale_color_gradient(low="blue", high="red") +
  theme_bw()

## Initiate writing to PDF file
pdf("gsea_TSLP_bubbleplot.pdf", height = 8, width = 10)

## Create a graphical object g here
TSLP_bubble

## Stop writing to the PDF file
dev.off()
#
# # Competitive GSEA using CAMERA----
# # for competitive tests the null hypothesis is that genes in the set are, at most, as often differentially expressed as genes outside the set
# # first let's create a few signatures to test in our enrichment analysis
# myTopHits_50 <- topTable(ebFit, adjust.method ="BH", coef=1, number=50, sort.by="logFC")
# mySig <- rownames(myTopHits_50) #if your own data is from mice, wrap this in 'toupper()' to make gene symbols all caps
# mySig2 <- sample((rownames(v.DEGList.filtered.norm$E)), size = 50, replace = FALSE)
# collection <- list(real = mySig, fake = mySig2)
# # now test for enrichment using CAMERA
# camera.res <- camera(v.DEGList.filtered.norm$E, collection, design, contrast.matrix[,1])
# camera.df <- as_tibble(camera.res, rownames = "setName")
# camera.df
#
# # Self-contained GSEA using ROAST----
# # remember that for self-contained the null hypothesis is that no genes in the set are differentially expressed
# mroast(v.DEGList.filtered.norm$E, collection, design, contrast=1) #mroast adjusts for multiple testing
#
# # now repeat with an actual gene set collection
# # camera requires collections to be presented as a list, rather than a tibble, so we must read in our signatures using the 'getGmt' function
# broadSet.C2.ALL <- getGmt("c2.all.v7.1.symbols.gmt", geneIdType=SymbolIdentifier())
# #extract as a list
# broadSet.C2.ALL <- geneIds(broadSet.C2.ALL)
# camera.res <- camera(v.DEGList.filtered.norm$E, broadSet.C2.ALL, design, contrast.matrix[,1])
# camera.df <- as_tibble(camera.res, rownames = "setName")
# camera.df
#
# # filter based on FDR and display as interactive table
# camera.df <- filter(camera.df, FDR<=0.05)
#
# datatable(camera.df,
#           extensions = c('KeyTable', "FixedHeader"),
#           caption = 'Signatures enriched in leishmaniasis',
#           options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
#   formatRound(columns=c(2,4,5), digits=2)
#
# #as before, add a variable that maps up/down regulated pathways with phenotype
# camera.df <- camera.df %>%
#   mutate(phenotype = case_when(
#     Direction == "Up" ~ "disease",
#     Direction == "Down" ~ "healthy"))
#
# #easy to filter this list based on names of signatures using 'str_detect'
# #here is an example of filtering to return anything that has 'CD8' or 'CYTOTOX' in the name of the signature
# camera.df.sub <- camera.df %>%
#   dplyr::filter(str_detect(setName, "CD8|CYTOTOX"))
#
# # graph camera results as bubble chart
# ggplot(camera.df[1:25,], aes(x=phenotype, y=setName)) +
#   geom_point(aes(size=NGenes, color = Direction, alpha=-log10(FDR))) +
#   theme_bw()

# Single sample GSEA using the GSVA package----
# the GSVA package offers a different way of approaching functional enrichment analysis.
# A few comments about the approach:
# In contrast to most GSE methods, GSVA performs a change in coordinate systems,
# transforming the data from a gene by sample matrix to a gene set (signature) by sample matrix.
# this allows for the evaluation of pathway enrichment for each sample.
# the method is both non-parametric and unsupervised
# bypasses the conventional approach of explicitly modeling phenotypes within enrichment scoring algorithms.
# focus is therefore placed on the RELATIVE enrichment of pathways across the sample space rather than the absolute enrichment with respect to a phenotype.
# however, with data with a moderate to small sample size (< 30), other GSE methods that explicitly include the phenotype in their model are more likely to provide greater statistical power to detect functional enrichment.
#
# # be aware that if you choose a large MsigDB file here, this step may take a while
# GSVA.res.C2CP <- gsva(v.DEGList.filtered.norm$E, #your data
#                       broadSet.C2.ALL, #signatures
#                       min.sz=5, max.sz=500, #criteria for filtering gene sets
#                       mx.diff=FALSE,
#                       method="gsva") #options for method are "gsva", ssgsea', "zscore" or "plage"
#
# # Apply linear model to GSVA result
# # now using Limma to find significantly enriched gene sets in the same way you did to find diffGenes
# # this means you'll be using topTable, decideTests, etc
# # note that you need to reference your design and contrast matrix here
# fit.C2CP <- lmFit(GSVA.res.C2CP, design)
# ebFit.C2CP <- eBayes(fit.C2CP)
#
# # use topTable and decideTests functions to identify the differentially enriched gene sets
# topPaths.C2CP <- topTable(ebFit.C2CP, adjust ="BH", coef=1, number=50, sort.by="logFC")
# res.C2CP <- decideTests(ebFit.C2CP, method="global", adjust.method="BH", p.value=0.05, lfc=0.58)
# # the summary of the decideTests result shows how many sets were enriched in induced and repressed genes in all sample types
# summary(res.C2CP)
#
# # pull out the gene sets that are differentially enriched between groups
# diffSets.C2CP <- GSVA.res.C2CP[res.C2CP[,1] !=0,]
# head(diffSets.C2CP)
# dim(diffSets.C2CP)
#
# # make a heatmap of differentially enriched gene sets
# hr.C2CP <- hclust(as.dist(1-cor(t(diffSets.C2CP), method="pearson")), method="complete") #cluster rows by pearson correlation
# hc.C2CP <- hclust(as.dist(1-cor(diffSets.C2CP, method="spearman")), method="complete") #cluster columns by spearman correlation
#
# # Cut the resulting tree and create color vector for clusters.  Vary the cut height to give more or fewer clusters, or you the 'k' argument to force n number of clusters
# mycl.C2CP <- cutree(hr.C2CP, k=2)
# mycolhc.C2CP <- rainbow(length(unique(mycl.C2CP)), start=0.1, end=0.9)
# mycolhc.C2CP <- mycolhc.C2CP[as.vector(mycl.C2CP)]
#
# # assign your favorite heatmap color scheme. Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]). Type demo.col(20) to see more color schemes.
# myheatcol <- colorRampPalette(colors=c("yellow","white","blue"))(100)
# # plot the hclust results as a heatmap
# heatmap.2(diffSets.C2CP,
#           Rowv=as.dendrogram(hr.C2CP),
#           Colv=as.dendrogram(hc.C2CP),
#           col=myheatcol, scale="row",
#           density.info="none", trace="none",
#           cexRow=0.9, cexCol=1, margins=c(10,14)) # Creates heatmap for entire data set where the obtained clusters are indicated in the color bar.
# # just as we did for genes, we can also make an interactive heatmap for pathways
# # you can edit what is shown in this heatmap, just as you did for your gene level heatmap earlier in the course

# the essentials ----
library(tidyverse)
library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots
gost.res_up <- gost(rownames(myModule_up), organism = "mmusculus", correction_method = "fdr")
mygostplot_up <- gostplot(gost.res_up, interactive = T, capped = T) #set interactive=FALSE to get plot for publications
# rerun the above gostplot function with 'interactive=F' and save to an object 'mygostplot'
publish_gostplot(
  mygostplot_up, #your static gostplot from above
  highlight_terms = c("GO:0002252", "GO:0006954", "GO:0019221", "GO:0042110",
                      "GO:0046649", "GO:0098542", "GO:0032103", "GO:0042330",
                      "GO:0016477", "GO:0012501", "GO:0031347", "GO:0042098",
                      "GO:0002440", "GO:0009617"),
  filename = "GO_plot_upreg.pdf",
  width = 8,
  height = 10)

#you can also generate a table of your gost results
publish_gosttable(
  gost.res_up,
  highlight_terms = NULL,
  use_colors = TRUE,
  show_columns = c("source", "term_name", "term_size", "intersection_size"),
  filename = "GO_table_upreg.png",
  ggplot=TRUE)

gost.res_down <- gost(rownames(myModule_down), organism = "mmusculus", correction_method = "fdr")
mygostplot_down <- gostplot(gost.res_down, interactive = F, capped = F) #set interactive=FALSE to get plot for publications

# rerun the above gostplot function with 'interactive=F' and save to an object 'mygostplot'
publish_gostplot(
  mygostplot_down, #your static gostplot from above
  highlight_terms = c("GO:0008610", "GO:0006629", "GO:0006694", "GO:0006631", "GO:0006695",
                      "GO:0046513", ""),
  filename = "GO_plot_downreg.pdf",
  width = NA,
  height = NA)

#you can also generate a table of your gost results
publish_gosttable(
  gost.res_down,
  highlight_terms = NULL,
  use_colors = TRUE,
  show_columns = c("source", "term_name", "term_size", "intersection_size"),
  filename = "GO_table_downreg.png",
  ggplot=TRUE)


hs_gsea_c2 <- msigdbr(species = "Mus musculus", # change depending on species your data came from
                      category = "C2") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature

# Now that you have your msigdb collections ready, prepare your data
# grab the dataframe you made in step3 script
# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC)
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE)
myGSEA.df <- as_tibble(myGSEA.res@result)

# view results as an interactive table
datatable(myGSEA.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:10), digits=2)
# create enrichment plots using the enrichplot package
gseaplot2(myGSEA.res,
          geneSetID = 47, #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          title = myGSEA.res$Description[47]) #can also turn off this title

# add a variable to this result that matches enrichment direction with phenotype
myGSEA.df <- myGSEA.df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "disease",
    NES < 0 ~ "healthy"))

# create 'bubble plot' to summarize y signatures across x phenotypes
ggplot(myGSEA.df[1:20,], aes(x=phenotype, y=ID)) +
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()


######## Perform GO analysis again, but on gene lists with p < 0.05 (instead of 0.01)
gost.res_up_05 <- gost(rownames(myModule_up_05), organism = "mmusculus", correction_method = "fdr")
mygostplot_up_05 <- gostplot(gost.res_up_05, interactive = T, capped = T) #set interactive=FALSE to get plot for publications
# rerun the above gostplot function with 'interactive=F' and save to an object 'mygostplot'
publish_gostplot(
  mygostplot_up, #your static gostplot from above
  highlight_terms = c("GO:0002252", "GO:0006954", "GO:0019221", "GO:0042110",
                      "GO:0046649", "GO:0098542", "GO:0032103", "GO:0042330",
                      "GO:0016477", "GO:0012501", "GO:0031347", "GO:0042098",
                      "GO:0002440", "GO:0009617"),
  filename = "GO_plot_upreg.pdf",
  width = 8,
  height = 10)

#you can also generate a table of your gost results
publish_gosttable(
  gost.res_up,
  highlight_terms = NULL,
  use_colors = TRUE,
  show_columns = c("source", "term_name", "term_size", "intersection_size"),
  filename = "GO_table_upreg.png",
  ggplot=TRUE)

gost.res_down_05 <- gost(rownames(myModule_down_05), organism = "mmusculus", correction_method = "fdr")
mygostplot_down_05 <- gostplot(gost.res_down_05, interactive = T, capped = F) #set interactive=FALSE to get plot for publications

# rerun the above gostplot function with 'interactive=F' and save to an object 'mygostplot'
publish_gostplot(
  mygostplot_down, #your static gostplot from above
  highlight_terms = c("GO:0008610", "GO:0006629", "GO:0006694", "GO:0006631", "GO:0006695",
                      "GO:0046513", ""),
  filename = "GO_plot_downreg.pdf",
  width = NA,
  height = NA)

#you can also generate a table of your gost results
publish_gosttable(
  gost.res_down,
  highlight_terms = NULL,
  use_colors = TRUE,
  show_columns = c("source", "term_name", "term_size", "intersection_size"),
  filename = "GO_table_downreg.png",
  ggplot=TRUE)



######### UP GO BAR PLOT

gost.res_up_data <- gost.res_up_05$result
gost.res_up_data_filter <- gost.res_up_data[c(132, 134, 148, 180, 253, 286, 292,
                                              329, 354, 496),]
gost.res_up_data_filter$group <- c("Cell Death", "Cell Death", "Cell Death", "Cell Death", 'Cell Death', 'Cell Death', 'Cell Death',
                                     "Cell Death","Cell Death","Cell Death")

GO_bar <- ggplot(data = gost.res_up_data_filter, aes(x=reorder(term_name, -p_value), y=-log10(p_value),
                                                       fill = group)) +
  geom_bar(stat="identity", color = "black") +
  geom_text(aes(label=intersection_size), vjust = 0.5,  hjust = -0.5 ,size = 4) +
  scale_fill_manual(values = "#A1A2A2") +
  theme_classic() +
  coord_flip() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title.x = element_text(size = 12)) +
  xlab("Gene Ontology Term") +
  labs(fill = "Gene Group")
GO_bar

pdf("GO_bar_plot_upreg.pdf", height = 4, width = 13)

## Create a graphical object g here
ggplot(data = gost.res_up_data_filter, aes(x=reorder(term_name, -p_value), y=-log10(p_value),
                                             fill = group)) +
  geom_bar(stat="identity", color = "black") +
  geom_text(aes(label=intersection_size), vjust = 0.5,  hjust = -0.5 ,size = 4) +
  scale_fill_manual(values = "#A1A2A2") +
  theme_classic() +
  coord_flip() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title.x = element_text(size = 12),
        legend.position = c(.8,.2)) +
  xlab("Gene Ontology Term") +
  labs(fill = "Gene Group") +
  ggtitle("Upregulated Gene Ontology Terms with 14 day TSLP treatment")
## Stop writing to the PDF file
dev.off()


#### DOWN GO BAR PLOT

gost.res_down_data <- gost.res_down_05$result
gost.res_down_data_filter <- gost.res_down_data[c(11, 13, 14, 17, 24, 31, 40, 42, 44, 48),]
gost.res_down_data_filter$group <- c("Lipid Metabolism", "Lipid Metabolism", "Lipid Metabolism", "Lipid Metabolism",
                                     "Lipid Metabolism","Lipid Metabolism","Lipid Metabolism", "Lipid Metabolism",
                                     "Lipid Metabolism", "Lipid Metabolism")

GO_bar <- ggplot(data = gost.res_down_data_filter, aes(x=reorder(term_name, -p_value), y=-log10(p_value),
                                                       fill = group)) +
  geom_bar(stat="identity", color = "black") +
  geom_text(aes(label=intersection_size), vjust = 0.5,  hjust = -0.5 ,size = 4) +
  theme_classic() +
  coord_flip() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title.x = element_text(size = 12)) +
  xlab("Gene Ontology Term") +
  labs(fill = "Gene Group")
GO_bar

pdf("GO_bar_plot_downreg.pdf", height = 4, width = 11)

## Create a graphical object g here
ggplot(data = gost.res_down_data_filter, aes(x=reorder(term_name, -p_value), y=-log10(p_value),
                                             fill = group)) +
  geom_bar(stat="identity", color = "black") +
  scale_fill_manual(values = "gold") +
  geom_text(aes(label=intersection_size), vjust = 0.5,  hjust = -0.5 ,size = 4) +
  theme_classic() +
  coord_flip() +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title.x = element_text(size = 12),
        legend.position = c(.8,.2)) +
  xlab("Gene Ontology Term") +
  labs(fill = "Gene Group") +
  ggtitle("Downregulated Gene Ontology Terms with 14 day TSLP treatment")
## Stop writing to the PDF file
dev.off()
