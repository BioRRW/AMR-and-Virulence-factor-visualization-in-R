### Reed Woyda
# 11/17/2020
# Script for doing PCA with contingency tables

# Load Libraries
# *must be installed prior to loading* tools->install packages"
library(viridis)
library(ggrepel)
library(scales)
library("fastcluster")
library("ggplot2")
library("factoextra")
library("FactoMineR")
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
rownames(data.df) <- data.df$isolate.id

# Initial view of CA to get idea
res1 <- CA(data.df[19:length(colnames(data.df))], ncp = 5, graph = TRUE)

print(res1)
get_eigenvalue(res1)
fviz_eig(res1)

# Output plot of dimension contributions
get_ca_row(res1)

#Various other metrics to look at
get_ca_col(res1)
fviz_ca_row(res1)
fviz_ca_col(res1)
fviz_ca_biplot(res1)


# repel= TRUE to avoid text overlapping (slow if many points)
# Output plot with labels to get idea
fviz_ca_biplot(res1, repel = TRUE)

# Color by cos2 values: quality on the factor map
# If a row item is well represented by two dimensions, the sum of the cos2 is closed to one.
# variables with low cos2 values will be colored in “white”
# variables with mid cos2 values will be colored in “blue”
# variables with high cos2 values will be colored in red
fviz_ca_row(res1, col.row = "cos2",
            gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
            repel = TRUE)

# +++++++++++++++++++ Alternative view: +++++++++++++++++++++++++++++
# Change the transparency by cos2 values
fviz_ca_row(res1, alpha.row="cos2")

# Example colors "green", "blue","red","black","aquamarine","chartreuse","chocolate1","darkgoldenrod1","darkmagenta","darkorange4","firebrick4","darksalmon","darkslateblue")

# Changing label on points
# Increasing this value will increase the number of labeled points on the plot
options(ggrepel.max.overlaps = 25)


# Output PDF file
# NOTE: legend needs to be edited externally 
# col.row = data.df$Serovar: "$Serovar" can be changed to any included metadata column of choice to color the data points
pdf(file = "Example-CA.pdf", height = 8, width = 8)
fviz_ca_row(res1, col.row = data.df$Serovar, title ="Correspondence analysis", repel = TRUE, addEllipses=FALSE, ellipse.level=0.95, label="all")
dev.off()