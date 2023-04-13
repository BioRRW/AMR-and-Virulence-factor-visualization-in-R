### Reed Woyda
# 11/17/2020
# Script for generating metadata-labeled heatmaps with presence/absence data
library(ggplot2)


library(dplyr)
library(viridis)
library(ggrepel)

library(gplots)
library(pheatmap)
library(dendextend)
library("fastcluster")
library("dendextend")
library("ggplot2")
library("ggdendro")
library("IntClust")
library("ggtree")
library("factoextra")
library(NbClust)
library(cluster)
library(mclust)
library("FactoMineR")
library("factoextra")
library("corrplot")


### Read in data in format: rownames = isolate, colnames = genes
# Ideally place metadata preceeding the presence/absence data
data.df <- as.data.frame(read.csv("Example-AMR-data.csv"))

# make factors, factors
# These are only for the metadata and result in categorical output on the heatmap
data.df$Sample.Date <- as.Date(data.df$Sample.Date, tryFormats = c("%m/%d/%Y"))
data.df$Broiler.House <- as.factor(data.df$Broiler.House)
data.df$Section <- as.factor(data.df$Section)
data.df$Area <- as.factor(data.df$Area)
data.df$Serovar <- as.factor(data.df$Serovar)
data.df$age <- as.factor(data.df$age)
data.df$Phase <- as.factor(data.df$Phase)
data.df$Flock.Cohort <- as.factor(data.df$Flock.Cohort)
data.df$House.temp. <- as.factor(data.df$House.temp.)
data.df$Moisture <- label_percent()(data.df$Moisture)
data.df$pH <- as.factor(data.df$pH)

# set row names of dataframe as the isolate names
# can alternatively be done upon read.csv() as well
rownames(data.df) <- data.df$ID.1

# simple heatmap to get an idea
# data.df[,19:length(colnames(data.df))]: select the data columns you wish to be included in the heatmap
pheatmap(data.df[,19:length(colnames(data.df))])

### Hierarchical clustering
# Fist calc distance
# data.df[,19:length(colnames(data.df))]: select the data columns you wish to be included in the heatmap
distance.mat <- vegdist(data.df[,19:length(colnames(data.df))], method = "jaccard", binary = TRUE)

# "complete" for abundance
# "ward.D" for pres/absence
hclust.results <- hclust(as.dist(distance.mat), method = "average", members = NULL)
hclust.results.den <- as.dendrogram(hclust.results)


# Determine optimal number of clusters
# data.df[,19:length(colnames(data.df))]: select the data columns you wish to be included in the heatmap
fviz_nbclust(data.df[,19:length(colnames(data.df))], kmeans, method='silhouette', k.max = 4)

# Use that optimal number of clusters to cut the tree
my_clusters <- data.frame(cutree(tree = as.dendrogram(hclust.results), k = 4))
colnames(my_clusters)[1] <- "Cluster"

# Plot the dendrogram
as.dendrogram(hclust.results) %>%
  #color_branches(k=19) %>%
  plot(horiz = TRUE)

# View Clusters
print(my_clusters)

# Add other metadata to annotate rows with
# my_clusters <- my_clusters[-1]: this removes the clusters, identified above, from the displayed metadata. Can remove if want to retain.
my_clusters <- my_clusters[-1]
#my_clusters$Section <- data.df$section
#my_clusters$Area <- data.df$area
#my_clusters$Temperature <- data.df$temp.house.cat
#my_clusters$Sample.Month <- data.df$Sample.date.cat
my_clusters$Serovar <- data.df$Serovar
#my_clusters$Phase <- data.df$Phase
my_clusters$Serovar <- data.df$Serovar
my_clusters$House <- data.df$Broiler.House
my_clusters$'Flock Cohort' <- data.df$Flock.Cohort
#my_clusters$pH <- data.df$ph.cat
#my_clusters$MLST <- data.df$MLST

# Add custom colors for metadata
my_colors = list(
  Sample.Month = c(August = "red", March = "green", June = "blue", May = "purple"),
  House = c('1' = "red", '2' = "green", '3' = "blue", '4' = "purple"),
  MLST = c('464' = "purple", '48' = "blue", '-' = "grey"),
  'Flock Cohort' = c('1' = "mediumpurple1", '2' = "deepskyblue1", '3' = "springgreen3"),
  Serovar = c(Senftenberg = "Blue",
              Enteritidis = "Red",
              Kentucky = "Green")
  )

# Generate heatmap to check
# cutree_rows = 4: must be changed to matche the optimally identified clusters above
# data.df[,19:length(colnames(data.df))]: select the data columns you wish to be included in the heatmap
# show_colnames = TRUE: this displays the column names (sample isolates in this case) on the heatmap
#    This is set to FALSE for the output to PDF
pheatmap(data.df[,19:length(colnames(data.df))], annotation_row = my_clusters[], annotation_colors = my_colors, cutree_rows = 4, show_colnames = TRUE)


# Generate same heatmap and output to PDF
# cutree_rows = 4: must be changed to matche the optimally identified clusters above
# data.df[,19:length(colnames(data.df))]: select the data columns you wish to be included in the heatmap
# show_colnames = TRUE: this displays the column names (sample isolates in this case) on the heatmap
#    This is set to FALSE for the output to PDF
pdf(file = "Heatmap.pdf", height = 7.8, width = 7)
pheatmap(data.df[,19:length(colnames(data.df))], annotation_row = my_clusters[], annotation_colors = my_colors, cutree_rows = 4, show_colnames = FALSE)
dev.off()

