##################################CLANS_Figure3B##################################

# Load necessary libraries
library(ggpubr)
library(factoextra)
library(ggrepel)

# Read the data
df <- read.table("clans_2d_DEF.csv", sep="\t", header=TRUE)

# Elbow method
fviz_nbclust(df[, c(2,3)], kmeans, method = "wss") +
  labs(subtitle = "Elbow method")
ggsave("clans.redo.test.elbow.pdf")

# Silhouette method
fviz_nbclust(df[, c(2,3)], kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette method")
ggsave("redo.test.silhouette.pdf")

# Gap statistic method
fviz_nbclust(df[, c(2,3)], kmeans, nstart = 25, method = "gap_stat", nboot = 100) +
  labs(subtitle = "Gap statistic method")
ggsave("redo.test.gap.pdf")

# K-means clustering
res.km <- kmeans(df[, c(2,3)], 5, nstart = 25)

# Add cluster information to the dataframe
newdf <- df
newdf$clusters <- as.character(res.km$cluster)

# Save the clustering results to a CSV file
write.table(newdf, file="clustering_results_2.csv", sep="\t", row.names=FALSE, quote=FALSE)

# Filter out 'Other' from the dataframe
testnewdf <- newdf[grep("Other", newdf$Proteins, invert=TRUE),]

# Create a vector of specific protein names to label
proteins_to_label <- c("NVE26090", "NVE21851", "NVE23160", "NVE12853", "NVE16599", "NVE16598", "NVE4706", "manIFIH1_HUMAN_7-99", "man_IFIH1_HUMAN_115-200", "MAVS_HUMAN", "DDX58_HUMAN_99-190", "DDX58_HUMAN_1-93")

# Filter the dataframe to include only the specific proteins to label
label_df <- newdf[newdf$Proteins %in% proteins_to_label, ]

# Create the plot
#ggplot() +
#  geom_point(data = newdf, mapping = aes(x = x, y = y, colour = Phyla)) +
#  geom_label_repel(data = label_df, colour = "black",
#                   aes(x = x, y = y, label = Proteins, fill = Phyla))


phyla_colors <- brewer.pal(6, "Dark2") 

# Define custom colors for the Phyla groups
phyla_colors <- c(
  "Cephalochordata" = "#CC79A7",  # orange - 5
  "Cnidaria" = "#D55E00",  # vermilion -3 
  "Echinodermata" = "#009E73",  # bluish Green -4
  "Protostomia" = "#56B4E9",  # sky blue -2 
  "Urochordata" = "#F0E442",  # Yellow -6
  "Vertebrates" = "#0072B2"   # Blue -1
)

# Create the plot with custom colors
ggplot() +
  geom_point(data = newdf, mapping = aes(x = x, y = y, colour = Phyla)) +
  geom_label_repel(data = label_df, colour = "black",
                   aes(x = x, y = y, label = Proteins, fill = Phyla)) +
  scale_colour_manual(values = phyla_colors) +  # Apply custom colors
  theme_minimal()


# Save the plot
ggsave("redo.5.cluster.mavs.ggpoint.pdf")

######################################################################################