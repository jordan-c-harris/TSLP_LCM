# Introduction to this script -----------
# this script walks thorough techniques for data exploration and expands on last week's data wrangling theme
# we'll also continue to create publication-quality graphics
# This script starts with your filtered and normalized abundance data from the Step 2 script.

# Load packages ------
library(tidyverse) # you're familiar with this fromt the past two lectures
library(DT) # for making interactive tables
library(plotly) # for making interactive plots
library(gt) # A layered 'grammar of tables' - think ggplot, but for tables

# Identify variables of interest in study design file ----
targets
group <- targets$group
group <- factor(group)
group <- factor(group, levels = c("Control", "TSLP"))

# Prepare your data -------
# for this part of the class you'll use your normalized and filtered data in log2 cpm
# make sure you have this object already in your work environment
# if you don't, go back to the Step2 script and generate it
log2.cpm.filtered.norm.df

# Hierarchical clustering ---------------
#hierarchical clustering can only work on a data matrix, not a data frame
#try using filtered and unfiltered data...how does this change the result?
#try other distance methods (e.g. switch from 'maximum' to 'euclidean')...how does this change the result?
distance <- dist(t(log2.cpm.filtered.norm), method = "euclidean") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(clusters, labels=sampleLabels)

## Initiate writing to PDF file
pdf("cluster_dendrogram.pdf", height = 5, width = 5)

## Create a graphical object g here
plot(clusters, labels=sampleLabels)

## Stop writing to the PDF file
dev.off()

# Principal component analysis (PCA) -------------
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
#look at the PCA result (pca.res) that you just created
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
pca.res$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
pca.res$x # 'x' shows you how much each sample influenced each PC (called 'loadings')
#note that these have a magnitude and a direction (this is the basis for making a PCA plot)
screeplot(pca.res) # A screeplot is a standard way to view eigenvalues for each PCA
pc.var<-pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
pc.per

# Visualize your PCA result ------------------
#lets first plot any two PCs against each other
#We know how much each sample contributes to each PC (loadings), so let's plot
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = group, fill = group) +
  geom_point(size=3, shape = 21, color = "black") +
  scale_fill_manual(values = c("blue", "red")) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) +
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="Control vs 7 day TSLP treatment") +
  coord_fixed() +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))
pca.plot

## Initiate writing to PDF file
pdf("PCA_plot.pdf", height = 5, width = 5)

## Create a graphical object g here
pca.plot

## Stop writing to the PDF file
dev.off()

# Let's discuss and iteratively refine the PCA code and plot from above
# First, take note of the fact that we can use information from our PCA analysis to label our axes
# Remember that PCA is unsupervised, so knows nothing about group assignment (healthy vs disease)
# But *we* know, and so we can use this knowledge to enhance the plot.  Add a 'color=group' mapping to the aes of the plot above
# Can we figure out the identity of the outlier?  We have already provided samplelabel mapping in aes, so just uncomment the 'geom_label()'
# Uncomment 'coord_fixed()' to apply the correct aspect ratio
# Uncomment 'stat_ellipse()' to see how you can circle clusters on the PCA
# How would this PCA look if you used raw counts (myCounts) instead of log2 CPM?
# What are the disadvantages of looking at a PCA result using such a simple XY plot?

# Create a PCA 'small multiples' chart ----
# this is another way to view PCA laodings to understand impact of each sample on each pricipal component
pca.res.df <- pca.res$x[,1:4] %>% # note that this is the first time you've seen the 'pipe' operator from the magrittr package
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group)

pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)

pca.pivot$sample <- factor(pca.pivot$sample, levels = sampleLabels)

sm <- ggplot(pca.pivot) +
  aes(x=forcats::fct_rev(sample),  y=loadings, fill=group) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot") +
  theme_bw() +
  coord_flip()
sm
## Initiate writing to PDF file
pdf("smallmultiples_plot.pdf", height = 5, width = 5)

## Create a graphical object g here
sm

## Stop writing to the PDF file
dev.off()

# Use dplyr 'verbs' to modify our dataframe ----
# use dplyr 'mutate' function to add new columns based on existing data
mydata.df <- log2.cpm.filtered.norm.df %>%
  mutate(Control.AVG = (Control1 + Control2 + Control3)/3,
         TSLP.AVG = (TSLP4 + TSLP5 + TSLP6)/3,
         #now make columns comparing each of the averages above that you're interested in
         LogFC = (TSLP.AVG - Control.AVG)) %>%
  mutate_if(is.numeric, round, 2)

#now look at this modified data table
mydata.df

# Use dplyr 'arrange' and 'select' to sort your dataframe based on any variable
# first, we'll use dplyr "arrange" function to sort rows based on the values in a column of interest
# then we'll display 'select' only the columns we're interested in seeing
mydata.sort <- mydata.df %>%
  dplyr::arrange(desc(LogFC)) %>%
  dplyr::select(geneID, LogFC)


# # the essentials ----
# library(tidyverse)
# library(DT)
# library(gt)
# library(plotly)
#
# group <- targets$group
# group <- factor(group)
#
# pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
# pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
# pc.per <- round(pc.var/sum(pc.var)*100, 1)
# pca.res.df <- as_tibble(pca.res$x)
# pca.plot <- ggplot(pca.res.df) +
#   aes(x=PC1, y=PC2, label=sampleLabels, color = group) +
#   geom_point(size=4) +
#   #stat_ellipse() +
#   xlab(paste0("PC1 (",pc.per[1],"%",")")) +
#   ylab(paste0("PC2 (",pc.per[2],"%",")")) +
#   labs(title="Sebaceous gland vs Bulk skin PCA plot"#,
#  #      caption=paste0("produced on ", Sys.time())) +
#   , color = "Group") +
#   coord_fixed() +
#   theme_bw()
# pca.plot
#
# ## Initiate writing to PDF file
# pdf("PCA_plot.pdf", height = 5, width = 5)
#
# ## Create a graphical object g here
# pca.plot
#
# ## Stop writing to the PDF file
# dev.off()
#
# ggplotly(pca.plot)
#
# mydata.df <- mutate(log2.cpm.filtered.norm.df,
#                     SG.AVG = (SG1 + SG2 + SG3)/3,
#                     Skin.AVG = (Skin4 + Skin5 + Skin6 + Skin7)/4,
#                     #now make columns comparing each of the averages above that you're interested in
#                     LogFC = (SG.AVG - Skin.AVG)) %>% #note that this is the first time you've seen the 'pipe' operator
#   mutate_if(is.numeric, round, 2)
#
# datatable(mydata.df[,c(1,8:10)],
#           extensions = c('KeyTable', "FixedHeader"),
#           filter = 'top',
#           options = list(keys = TRUE,
#                          searchHighlight = TRUE,
#                          pageLength = 10,
#                          #dom = "Blfrtip",
#                          #buttons = c("copy", "csv", "excel"),
#                          lengthMenu = c("10", "25", "50", "100")))
