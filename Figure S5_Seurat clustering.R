# Run on Seurat 3.1.1

# Load necessary libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(raster)

# Specify data file path
data_file_path <- ".../Data/Figure S5/rawdata/MPOA_female_summary_result.csv"

# Load the data frame storing all cell information
MPN_HCR_df <- read.csv(file = data_file_path)

# Subset the MPOA_pbmc to MPOM/L (MPN)
MPN_HCR_df <- MPN_HCR_df %>% 
  filter(Analysis.major.region == 'MPN')

# Define a vector of gene names of interest
genes <- c("Fos", "Vglut2", "Vgat", "Esr1", "Gal", "Calcr", "Nts", "Prlr");


# Initialize empty vectors to store column names
subset_columns_raw_copies <- c()
subset_columns_cell_intensity <- c()
subset_columns_regressed_copies <- c()
subset_columns_positive_cells <- c()

# Generate column names for each gene
for(gene in genes){
  subset_columns_raw_copies <- c(subset_columns_raw_copies, sprintf("%s.Copies", gene))
  subset_columns_cell_intensity <- c(subset_columns_cell_intensity, sprintf("%s.Cell.Intensity", gene))
  subset_columns_regressed_copies <- c(subset_columns_regressed_copies, sprintf("%s.regressed_copies", gene))
  subset_columns_positive_cells <- c(subset_columns_positive_cells, sprintf("%s.positive_cells", gene))
}

# Generate a sequence of row names
row_names <- seq(from = 1, to = nrow(MPN_HCR_df), by = 1)

# Assign new row names to the dataframe
row.names(MPN_HCR_df) <- row_names

# Define the metadata columns
meta_columns <- c("Analysis.Region", "Analysis.major.region", "Hemisphere",
                  "ID", "section", "Sex", "Stim", "Sex_Stim")

# Subset data based on the generated column names
raw_copies <- MPN_HCR_df[subset_columns_raw_copies]
cell_intensity <- MPN_HCR_df[subset_columns_cell_intensity]
regressed_copies <- MPN_HCR_df[subset_columns_regressed_copies]
positive_cells <- MPN_HCR_df[subset_columns_positive_cells]

# Create a metadata dataframe
meta <- MPN_HCR_df[c(meta_columns, subset_columns_cell_intensity)]

# Cutoff the negative cells using boolean matrix
for(i in seq(from = 1, to = 8, by = 1)){
  # If a cell is not positive (as indicated by positive_cells), set its value to 0 in regressed_copies and raw_copies
  regressed_copies[!as.logical(positive_cells[,i]), i] <- 0
  raw_copies[!as.logical(positive_cells[,i]), i] <- 0
}

# Remove cells that have no gene expression at all
# Calculate the maximum expression value for each row (cell) across all genes
max_expr <- apply(raw_copies, 1, max)

# Create a logical vector indicating which cells have a maximum expression value greater than 0
high_expr_rows <- max_expr > 0

# Subset the dataframes based on the cells with non-zero maximum expression
meta <- meta[high_expr_rows, ]
regressed_copies <- regressed_copies[high_expr_rows, ]


# Transpose the 'regressed_copies' dataframe and convert it back to a dataframe
regressed_copies.t <- transpose(regressed_copies)
regressed_copies.t <- as.data.frame(t(as.matrix(regressed_copies.t)))

# Remove the Fos data for the clustering
remove_Fos <- regressed_copies[, 2:8]

# Transpose the 'remove_Fos' dataframe and convert it back to a dataframe
remove_Fos.t <- transpose(remove_Fos)
remove_Fos.t <- as.data.frame(t(as.matrix(remove_Fos.t)))


# Create a Seurat object
MPOA_pbmc = CreateSeuratObject(remove_Fos.t, project = "MPOA_HCR", assay = "RNA",
                               min.cells = 0, min.features = 1, names.field = 1,
                               names.delim = "_", meta.data = meta
)


# Define all genes as the rownames of the Seurat object
all.genes <- rownames(MPOA_pbmc)

# Define a list of genes excluding 'Fos.Copies'
nonFos.genes <- all.genes[all.genes != "Fos.Copies"]

# Define the number of dimensions equal to the length of nonFos.genes
dims <- length(nonFos.genes)

# Scale the data in the Seurat object
MPOA_pbmc <- ScaleData(
  object = MPOA_pbmc,
  features = all.genes,
  vars.to.regress = c('nCount_RNA', subset_columns_cell_intensity)
)

# Run PCA on the data in the Seurat object
MPOA_pbmc <- RunPCA(
  object = MPOA_pbmc,
  features = nonFos.genes,
  approx = FALSE
)

# Visualize the loadings of the PCA
VizDimLoadings(
  object = MPOA_pbmc,
  dims = 1:dims,
  reduction = "pca"
)

# Create a dimensionality reduction plot
DimPlot(
  object = MPOA_pbmc,
  reduction = "pca"
)


# Find nearest neighbors in the PCA space
MPOA_pbmc <- FindNeighbors(MPOA_pbmc, dims = 1:dims)

# Find clusters based on the nearest neighbors graph
MPOA_pbmc <- FindClusters(MPOA_pbmc, resolution = 0.08)  # for regressed copies

# Run UMAP on the Seurat object
MPOA_pbmc <- RunUMAP(MPOA_pbmc, dims = 1:dims, n.neighbors = 100, min.dist = 0.05)

# Define the features to be used for plots
features <- rownames(MPOA_pbmc)

# Plot the UMAP visualization
DimPlot(MPOA_pbmc, reduction = "umap")

# Create a DotPlot for the features
DotPlot(MPOA_pbmc, features = features) + 
  scale_x_discrete(labels = rlist::list.reverse(features))  # reverse the order of labels

# Create a FeaturePlot based on 'umap' reduction
FeaturePlot(
  object = MPOA_pbmc,
  reduction = 'umap',
  features = subset_columns_cell_intensity,
  max.cutoff = 20
)

# Create dimension plots grouped by different variables
DimPlot(MPOA_pbmc, reduction = "umap", group.by = "Analysis.major.region")
DimPlot(MPOA_pbmc, reduction = "umap", group.by = "Stim")
DimPlot(MPOA_pbmc, reduction = "umap", group.by = "ID")

# Create a FeaturePlot for 'Fos.regressed-copies'
FeaturePlot(
  object = MPOA_pbmc,
  reduction = 'umap',
  features = "Fos.regressed-copies",
  max.cutoff = 100
)

# Create a combined plot with DimPlot, FeaturePlot and VlnPlot
DimPlot(MPOA_pbmc, reduction = "umap") + 
  FeaturePlot(
    object = MPOA_pbmc,
    reduction = 'umap',
    features = "Fos.regressed-copies"
  ) + 
  VlnPlot(
    object = MPOA_pbmc,
    features = "Fos.regressed-copies", 
    split.by = "Sex_Stim"
  )



# Load the ggplot2 package
library(ggplot2)

# Set the folder path
folder_path <- ".../Data/Figure S5/Seurat figures"

# Set the file names for the PDF, PNG, and SVG files
pdf_file <- paste0(folder_path, "FeaturePlot.pdf")
png_file <- paste0(folder_path, "FeaturePlot.png")
svg_file <- paste0(folder_path, "FeaturePlot.svg")

# Generate a FeaturePlot
FeaturePlot <- FeaturePlot(MPOA_pbmc,reduction = 'umap',features = subset_columns_cell_intensity,max.cutoff = 20,)

# Save the FeaturePlot as a PDF
pdf(pdf_file)
print(FeaturePlot)
dev.off()

# Save the FeaturePlot as a PNG
png(png_file)
print(FeaturePlot)
dev.off()

# Save the FeaturePlot as a SVG
svg(svg_file)
print(FeaturePlot)
dev.off()


# Set the file names for the PDF, PNG, and SVG files
pdf_file <- paste0(folder_path, "DimPlot.pdf")
png_file <- paste0(folder_path, "DimPlot.png")
svg_file <- paste0(folder_path, "DimPlot.svg")

# Generate a DimPlot
DimPlot <- DotPlot(MPOA_pbmc, features = features,col.max = 0.5) + scale_x_discrete(labels = rlist::list.reverse(genes))

# Save the DimPlot as a PDF
pdf(pdf_file)
print(DimPlot)
dev.off()

# Save the DimPlot as a PNG
png(png_file)
print(DimPlot)
dev.off()

# Save the DimPlot as a SVG
svg(svg_file)
print(DimPlot)
dev.off()

# Set the file names for the PDF, PNG, and SVG files
pdf_file <- paste0(folder_path, "DimPlot2.pdf")
png_file <- paste0(folder_path, "DimPlot2.png")
svg_file <- paste0(folder_path, "DimPlot2.svg")

# Generate a DimPlot
DimPlot <- DimPlot(MPOA_pbmc, reduction = "umap")


# Save the DimPlot as a PDF
pdf(pdf_file)
print(DimPlot)
dev.off()

# Save the DimPlot as a PNG
png(png_file)
print(DimPlot)
dev.off()

# Save the DimPlot as a SVG
svg(svg_file)
print(DimPlot)
dev.off()

# Load the ggplot2 package
library(ggplot2)


# Basic plot of clusters by Sex_Stim
plot1 <- ggplot(MPOA_pbmc@meta.data, aes(x=seurat_clusters, fill=Sex_Stim)) + geom_bar()


ggsave(paste0(folder_path, "plot1_clusters_by_Sex_Stim.pdf"), plot1)
ggsave(paste0(folder_path, "plot1_clusters_by_Sex_Stim.png"), plot1)

# Plot as proportion or percentage of cluster
plot2 <- ggplot(MPOA_pbmc@meta.data, aes(x=seurat_clusters, fill=Sex_Stim)) + geom_bar(position = "fill")
ggsave(paste0(folder_path, "plot2_clusters_by_Sex_Stim_proportion.pdf"), plot2)
ggsave(paste0(folder_path, "plot2_clusters_by_Sex_Stim_proportion.png"), plot2)

# Basic plot of clusters by Analysis.major.region
plot3 <- ggplot(MPOA_pbmc@meta.data, aes(x=seurat_clusters, fill=Analysis.major.region)) + geom_bar()
ggsave(paste0(folder_path, "plot3_clusters_by_Analysis_major_region.pdf"), plot3)
ggsave(paste0(folder_path, "plot3_clusters_by_Analysis_major_region.png"), plot3)

# Plot as proportion or percentage of cluster
plot4 <- ggplot(MPOA_pbmc@meta.data, aes(x=seurat_clusters, fill=Analysis.major.region)) + geom_bar(position = "fill")
ggsave(paste0(folder_path, "plot4_clusters_by_Analysis_major_region_proportion.pdf"), plot4)
ggsave(paste0(folder_path, "plot4_clusters_by_Analysis_major_region_proportion.png"), plot4)

# Basic plot of clusters by section
plot5 <- ggplot(MPOA_pbmc@meta.data, aes(x=seurat_clusters, fill=section)) + geom_bar()
ggsave(paste0(folder_path, "plot5_clusters_by_section.pdf"), plot5)
ggsave(paste0(folder_path, "plot5_clusters_by_section.png"), plot5)

# Plot as proportion or percentage of cluster
plot6 <- ggplot(MPOA_pbmc@meta.data, aes(x=seurat_clusters, fill=section)) + geom_bar(position = "fill")
ggsave(paste0(folder_path, "plot6_clusters_by_section_proportion.pdf"), plot6)
ggsave(paste0(folder_path, "plot6_clusters_by_section_proportion.png"), plot6)



'
#basic plot of clusters by Sex_Stim
ggplot(MPOA_pbmc@meta.data, aes(x=seurat_clusters, fill=Sex_Stim)) + geom_bar()

#plot as proportion or percentage of cluster
ggplot(MPOA_pbmc@meta.data, aes(x=seurat_clusters, fill=Sex_Stim)) + geom_bar(position = "fill")

#basic plot of clusters by Sex_Stim
ggplot(MPOA_pbmc@meta.data, aes(x=seurat_clusters, fill=Analysis.major.region)) + geom_bar()

#plot as proportion or percentage of cluster
ggplot(MPOA_pbmc@meta.data, aes(x=seurat_clusters, fill=Analysis.major.region)) + geom_bar(position = "fill")

#basic plot of clusters by Sex_Stim
ggplot(MPOA_pbmc@meta.data, aes(x=seurat_clusters, fill=section)) + geom_bar()

#plot as proportion or percentage of cluster
ggplot(MPOA_pbmc@meta.data, aes(x=seurat_clusters, fill=section)) + geom_bar(position = "fill")
'

regressed_copies_csv_path = ".../Data/Figure S3/Seurat figures/MPOA_female_summary_regressed_Seurat_MPN.csv"
# save csvs 
write.csv(regressed_copies, regressed_copies_csv_path,)

raw_copies_csv_path = ".../Data/Figure S3/Seurat figures/MPOA_female_summary_raw_Seurat_MPN.csv"

# save csvs 
write.csv(raw_copies, raw_copies_csv_path, )
#[,MPOA_pbmc$orig.ident]

seurat_clusters_path = ".../Data/Figure S3/Seurat figures/MPOA_female_summary_seurat_clusters_Seurat_MPN.csv"
# Save the cluster assignments as a CSV file
write.csv(x = MPOA_pbmc$seurat_clusters, file = seurat_clusters_path, row.names = TRUE)


meta_path = ".../Data/Figure S3/Seurat figures/MPOA_female_summary_meta_Seurat_MPN.csv"
# Save the cluster assignments as a CSV file
write.csv(x = meta, file = meta_path, row.names = TRUE)

# Save the Seurat object
object_path = ".../Data/Figure S3/Seurat figures/MPOA_female_summary_object_Seurat_MPN.rds"

saveRDS(MPOA_pbmc, file = "MPOA_pbmc.rds")
MPOA_pbmc = readRDS( file = "MPOA_pbmc.rds")
