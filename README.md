Cnidarian Immune System – Analysis Scripts

R code used to generate the analyses and figures for the manuscript:

DEGs.R

Pipeline for differential gene expression using DESeq2 and edgeR, PCA, VST/rlog visualization, correlation heatmaps, and GO-term enrichment.

Heatmap.R

Generates a colorblind-friendly heatmap of immune gene expression using ComplexHeatmap, viridis, and circlize (Figure 4F).

Correlations_SingleCell.R

Computes and visualizes the correlation between RLRb–CARDIB expression across metacells and experimental conditions using Seurat data (Figure 2C–F).

CLANS.R

Used to generate 2D CLANS clustering of CARD-like proteins and highlight key antiviral signaling proteins (Figure 3B).
