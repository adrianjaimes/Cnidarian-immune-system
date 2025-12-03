#######################Figure 4f_Heatmap ###################################


#######################colorblind people###################################

library(pheatmap)
library(ComplexHeatmap)
library(circlize)  # For colorRamp2
library(viridis)   # For colorblind-friendly color palettes

df <- read.csv("ImmuneGenes_56.csv")
rownames(df) <- df[, 1]  # Set the first column as row names
df <- df[, -1]           # Remove the first column from the data frame

df <- as.matrix(df)

# Define custom color breaks and colors (adjust these values based on your data range and needs)
breaks <- c(-4, -2, 0, 2, 4, 6)

# Use the viridis color palette, which is colorblind-friendly
colors <- viridis(length(breaks), option = "D")

# Create a custom color mapping
color_mapping <- colorRamp2(breaks, colors)

# Generate the heatmap with the custom color scale
Heatmap(df, 
        name = "LFC",  # Customize the heatmap title
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        col = color_mapping)