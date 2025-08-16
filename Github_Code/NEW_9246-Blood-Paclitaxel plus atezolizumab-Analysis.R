# -------------------- Uninstalled Packages

renv::init()
remove.packages("Matrix")
install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-4.tar.gz", repos = NULL, type = "source")
install.packages("~/Seurat_4.3.0.tar.gz", repos = NULL, type = "source")
install.packages(c("Matrix", "Rcpp", "ggplot2"))
install.packages("remotes")
remotes::install_version("Seurat", version = "4.3.0", repos = "http://cran.us.r-project.org")
packageVersion("Seurat") # 4.3.0

renv::init()
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", Version = "3.18")
BiocManager::install(version = "3.18")
BiocManager::version()
BiocManager::install("GenomicAlignment")
packageVersion("GenomicAlignment")
#install.packages(c("pkg1", "pkg2"))
install.packages("Matrix", version ="1.6.4")
#remotes::install_version("Seurat", version = "4.3.0")

#install.packages("GenomicAlignments", lib = "C:/Users/HC/AppData/Local/R/win-library/4.3")
#install.packages("C:/seurat-4.3.0.zip", repos = NULL, type = "win.binary")
install.packages("Seurat", version = "4.3.0")
remotes::install_version("SeuratObject", version = "4.1.3", repos = "http://cran.us.r-project.org")
# ----------------------- Installed Packages

# Install dependencies for AnnotationHub
BiocManager::install(c("BiocGenerics", "BiocFileCache", "RSQLite", "AnnotationDbi", "S4Vectors", "interactiveDisplayBase"))
# Install AnnotationHub from source
install.packages("https://bioconductor.org/packages/3.18/bioc/src/contrib/AnnotationHub_3.10.1.tar.gz", repos = NULL, type = "source")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", Version = "3.12")
BiocManager::install("SingleCellExperiment", version = "1.24.0")
BiocManager::install("rtracklayer", lib = "C:/Users/HC/AppData/Local/R/win-library/4.3", force = TRUE)
BiocManager::install("SingleR")
BiocManager::install("ensembldb", Version = "2.26.0") # --------- packageVersion("ensembldb")= 2.26.1
BiocManager::install("celldex")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("rjson")
library(rjson)
BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
library(ensembldb)

BiocManager::install("ensembldb")

BiocManager::install("ensembldb")
install.packages("devtools")
install.packages("patchwork")
install.packages("Rcpp")
install.packages("conflicted")
install.packages("DropletUtils")
install.packages("Rcmdr", dependencies = TRUE)
install.packages("umap")
install.packages("ExperimentHub")

remotes::install_version("scales", version = "1.3.0")
remotes::install_version("cowplot", version = "1.1.3")
remotes::install_version("RCurl", version = "1.98.1.13")
install.packages("devtools")
devtools::install_version("Matrix", version = "1.6.4")
devtools::install_version("tidyverse", version = "2.0.0")
BiocManager::install("DropletUtils")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("multtest")

install.packages("metap")

# -------------------------- Packages Version
packageVersion("Seurat") # 4.3.0
packageVersion("ggplot2")
packageVersion("tidyr")
packageVersion("dplyr")
packageVersion("AnnotationHub")
packageVersion("SingleCellExperiment")
packageVersion("BiocFileCache")
packageVersion("SingleR")
packageVersion("patchwork")
packageVersion("rtracklayer")
packageVersion("devtools")
packageVersion("Rcmdr")
packageVersion("Matrix")
packageVersion("SeuratObject")


# --------------- Necessary Packages Version 
packageVersion("celldex") # 1.12.0
packageVersion("GenomicAlignment")
packageVersion("ensembldb") # 2.26.0
remotes::install_version("tidyverse", version = "2.0.0")

# ----------------- Load Library
#library(SingleCellExperiment)
library(Matrix)
library(ggplot2)
library(Rcpp)
library(celldex)
library(Seurat)
library(conflicted)
library(dplyr)
library(BiocManager)
library(AnnotationHub)
library(ExperimentHub)
library(DropletUtils)
library(SingleR)
library(umap)
#detach("package:devtools", unload = TRUE)
#library(devtools)
library(Matrix)

library(tidyverse)
library(scales)
library(cowplot)
library(RCurl)

library(patchwork)
library(tidyr)
library(ensembldb)
library(tibble)
library(metap)
library(tibble)
library(tibble)  
library(readr)
library(readxl)
library(tibble)  
library(metap)
#install.packages("purrr")
library(purrr)
#install.packages("stringr")
library(stringr)
library(readr)

# Install the package
#install.packages("remotes")
#install.packages("renv")
#renv::init()
#remotes::install_version("Seurat", version = "4.3.0")
#------------------------------------------------------------------------------------
# --------------------- Read Samples GSE169246
# --------------------- Read Samples GSE169246------Post_P001_b---paclitaxel plus atezolizumab
sample1_path1 <- "D:\\Thesis\\Proposal\\Thesis_Data\\GSE169246\\files\\Post_P001_b"
list.files(sample1_path1)
data1 <- Read10X(data.dir = sample1_path1)
head(data1)
seurat_JCML1 <- CreateSeuratObject(counts = data1, min.features = 100,project="Post_P001_b")
seurat_JCML1$orig.ident <- "Post_P001_b"
seurat_JCML1$treatment_conditions <- "Post"
seurat_sample1 <- seurat_JCML1
view(seurat_JCML1)
unique(seurat_JCML1@meta.data$orig.ident)
# --------------------- Read Samples GSE169246------Post_P002_b---paclitaxel plus atezolizumab
sample1_path1 <- "D:\\Thesis\\Proposal\\Thesis_Data\\GSE169246\\files\\Post_P002_b"
list.files(sample1_path1)
data2 <- Read10X(data.dir = sample1_path1)
head(data2)
seurat_JCML2 <- CreateSeuratObject(counts = data2, min.features = 100,project="Post_P002_b")
seurat_JCML2$orig.ident <- "Post_P002_b"
seurat_JCML2$treatment_conditions <- "Post"
seurat_sample2 <- seurat_JCML2
view(seurat_JCML2)
unique(seurat_JCML2@meta.data$orig.ident)

# --------------------- Read Samples GSE169246------Post_P004_b---paclitaxel plus atezolizumab
sample1_path1 <- "D:\\Thesis\\Proposal\\Thesis_Data\\GSE169246\\files\\Post_P004_b"
list.files(sample1_path1)
data3 <- Read10X(data.dir = sample1_path1)
head(data3)
seurat_JCML3 <- CreateSeuratObject(counts = data3, min.features = 100,project="Post_P004_b")
seurat_JCML3$orig.ident <- "Post_P004_b"
seurat_JCML3$treatment_conditions <- "Post"
seurat_sample3 <- seurat_JCML3
view(seurat_JCML3)
unique(seurat_JCML3@meta.data$orig.ident)

# --------------------- Read Samples GSE169246------Post_P005_b---paclitaxel plus atezolizumab
sample1_path1 <- "D:\\Thesis\\Proposal\\Thesis_Data\\GSE169246\\files\\Post_P005_b"
list.files(sample1_path1)
data4 <- Read10X(data.dir = sample1_path1)
head(data4)
seurat_JCML4 <- CreateSeuratObject(counts = data4, min.features = 100,project="Post_P005_b")
seurat_JCML4$orig.ident <- "Post_P005_b"
seurat_JCML4$treatment_conditions <- "Post"
seurat_sample4 <- seurat_JCML4
view(seurat_JCML4)
unique(seurat_JCML4@meta.data$orig.ident)

# --------------------- Read Samples GSE169246------Post_P007_b---paclitaxel plus atezolizumab
sample1_path1 <- "D:\\Thesis\\Proposal\\Thesis_Data\\GSE169246\\files\\Post_P007_b"
list.files(sample1_path1)
data5 <- Read10X(data.dir = sample1_path1)
head(data5)
seurat_JCML5 <- CreateSeuratObject(counts = data5, min.features = 100,project="Post_P007_b")
seurat_JCML5$orig.ident <- "Post_P007_b"
seurat_JCML5$treatment_conditions <- "Post"
seurat_sample5 <- seurat_JCML5
view(seurat_JCML5)
unique(seurat_JCML5@meta.data$orig.ident)
# --------------------- Read Samples GSE169246------Post_P019_b---paclitaxel plus atezolizumab
sample1_path1 <- "D:\\Thesis\\Proposal\\Thesis_Data\\GSE169246\\files\\Post_P019_b"
list.files(sample1_path1)
data6 <- Read10X(data.dir = sample1_path1)
head(data6)
seurat_JCML6 <- CreateSeuratObject(counts = data6, min.features = 100,project="Post_P019_b")
seurat_JCML6$orig.ident <- "Post_P019_b"
seurat_JCML6$treatment_conditions <- "Post"
seurat_sample6 <- seurat_JCML6
view(seurat_JCML6)
unique(seurat_JCML6@meta.data$orig.ident)
# --------------------- Read Samples GSE169246------Pre_P001_b---paclitaxel plus atezolizumab
sample1_path1 <- "D:\\Thesis\\Proposal\\Thesis_Data\\GSE169246\\files\\Pre_P001_b"
list.files(sample1_path1)
data7 <- Read10X(data.dir = sample1_path1)
head(data7)
seurat_JCML7 <- CreateSeuratObject(counts = data7, min.features = 100,project="Pre_P001_b")
seurat_JCML7$orig.ident <- "Pre_P001_b"
seurat_JCML7$treatment_conditions <- "Pre"
seurat_sample7 <- seurat_JCML7
view(seurat_JCML7)
unique(seurat_JCML7@meta.data$orig.ident)

# --------------------- Read Samples GSE169246------Pre_P002_b---paclitaxel plus atezolizumab
sample1_path1 <- "D:\\Thesis\\Proposal\\Thesis_Data\\GSE169246\\files\\Pre_P002_b"
list.files(sample1_path1)
data8 <- Read10X(data.dir = sample1_path1)
head(data8)
seurat_JCML8 <- CreateSeuratObject(counts = data8, min.features = 100,project="Pre_P002_b")
seurat_JCML8$orig.ident <- "Pre_P002_b"
seurat_JCML8$treatment_conditions <- "Pre"
seurat_sample8 <- seurat_JCML8
view(seurat_JCML8)
unique(seurat_JCML8@meta.data$orig.ident)

# --------------------- Read Samples GSE169246------Pre_P004_b---paclitaxel plus atezolizumab
sample1_path1 <- "D:\\Thesis\\Proposal\\Thesis_Data\\GSE169246\\files\\Pre_P004_b"
list.files(sample1_path1)
data9 <- Read10X(data.dir = sample1_path1)
head(data9)
seurat_JCML9 <- CreateSeuratObject(counts = data9, min.features = 100,project="Pre_P004_b")
seurat_JCML9$orig.ident <- "Pre_P004_b"
seurat_JCML9$treatment_conditions <- "Pre"
seurat_sample9 <- seurat_JCML9
view(seurat_JCML9)
unique(seurat_JCML9@meta.data$orig.ident)

# --------------------- Read Samples GSE169246------Pre_P005_b---paclitaxel plus atezolizumab
sample1_path1 <- "D:\\Thesis\\Proposal\\Thesis_Data\\GSE169246\\files\\Pre_P005_b"
list.files(sample1_path1)
data10 <- Read10X(data.dir = sample1_path1)
head(data10)
seurat_JCML10 <- CreateSeuratObject(counts = data10, min.features = 100,project="Pre_P005_b")
seurat_JCML10$orig.ident <- "Pre_P005_b"
seurat_JCML10$treatment_conditions <- "Pre"
seurat_sample10 <- seurat_JCML10
view(seurat_JCML10)
unique(seurat_JCML10@meta.data$orig.ident)

# --------------------- Read Samples GSE169246------Pre_P007_b---paclitaxel plus atezolizumab
sample1_path1 <- "D:\\Thesis\\Proposal\\Thesis_Data\\GSE169246\\files\\Pre_P007_b"
list.files(sample1_path1)
data11 <- Read10X(data.dir = sample1_path1)
head(data11)
seurat_JCML11 <- CreateSeuratObject(counts = data11, min.features = 100,project="Pre_P007_b")
seurat_JCML11$orig.ident <- "Pre_P007_b"
seurat_JCML11$treatment_conditions <- "Pre"
seurat_sample11 <- seurat_JCML11
view(seurat_JCML11)
unique(seurat_JCML11@meta.data$orig.ident)
# --------------------- Read Samples GSE169246------Pre_P019_b---paclitaxel plus atezolizumab
sample1_path1 <- "D:\\Thesis\\Proposal\\Thesis_Data\\GSE169246\\files\\Pre_P019_b"
list.files(sample1_path1)
data12 <- Read10X(data.dir = sample1_path1)
head(data12)
seurat_JCML12 <- CreateSeuratObject(counts = data12, min.features = 100,project="Pre_P019_b")
seurat_JCML12$orig.ident <- "Pre_P019_b"
seurat_JCML12$treatment_conditions <- "Pre"
seurat_sample12 <- seurat_JCML12
view(seurat_JCML12)
unique(seurat_JCML12@meta.data$orig.ident)
# --------------------------------------------------Create a merged object GSE169246
merged_sample1 <- merge(x = seurat_sample1, 
                        y = c(seurat_sample2, seurat_sample3, seurat_sample4,
                              seurat_sample5, seurat_sample6, seurat_sample7,
                              seurat_sample8, seurat_sample9, seurat_sample10, 
                              seurat_sample11,seurat_sample12
                              ), 
                        add.cell.id = c("Post","Post","Post","Post","Post","Post",
                                        "Pre","Pre","Pre","Pre","Pre","Pre"))

# ------------------------------- NUMBER OF ALL CELLS
all_samples <- list(seurat_sample1,seurat_sample2, seurat_sample3, seurat_sample4, 
                    seurat_sample5, seurat_sample6, seurat_sample7, 
                    seurat_sample8, seurat_sample9,seurat_sample10, 
                    seurat_sample11,seurat_sample12)

cell_counts <- sapply(all_samples, function(x) ncol(x))
cell_counts
sum(cell_counts)

# ---------------------------------------------------First sample sample
merged_sample <- merged_sample1
head(merged_sample@meta.data)
tail(merged_sample@meta.data)
# Explore merged metadata
View(merged_sample@meta.data)

# Add number of genes per UMI for each cell to metadata
merged_sample$log10GenesPerUMI <- log10(merged_sample$nFeature_RNA) / log10(merged_sample$nCount_RNA)
View(merged_sample)
# Compute percent mito ratio
#NOTE: The pattern provided (“^MT-“) for human and "^mt-" for mouse works for human gene names. You may need to adjust depending on your organism of interest. 
#If you weren’t using gene names as the gene ID, then this function wouldn’t work. 
head(rownames(merged_sample))
grep("^MT-", rownames(merged_sample), value = TRUE)
mt_genes <- grep("^MT-", rownames(merged_sample), value = TRUE)
mt_data <- merged_sample[mt_genes, ]

#grep( "^mt-", rownames(merged_sample), value = T)
# Calculate the mitochondrial gene percentage
merged_sample$mitoRatio <- PercentageFeatureSet(object = merged_sample, pattern = "^MT-")
View(merged_sample)
head(merged_sample)
# Normalize the mitoRatio (if needed, depending on how the PercentageFeatureSet function works)
merged_sample$mitoRatio <- merged_sample$mitoRatio / 100
head(merged_sample)

# Create metadata dataframe
metadata <- merged_sample@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
head(metadata)
# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
head(metadata)

# Create sample column
# View the first few rows of metadata
head(merged_sample@meta.data)
merged_sample@meta.data <- metadata
head(metadata)
metadata$sample <- "Post"
metadata$sample[metadata$treatment_conditions == "Pre"] <- "Pre"
head(metadata)
tail(metadata)

# Add metadata back to Seurat object
merged_sample@meta.data <- metadata

# Create .RData object to load at any time
save(merged_sample, file="F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\RData1.RData")
load("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\RData1.RData")
ls()
head(merged_sample)

# Visualize the number of cell counts per sample
# Create directory if it doesn't exist
dir.create("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs", recursive = TRUE, showWarnings = FALSE)

# Create bar plot of cell counts per sample
metadata %>%
  ggplot(aes(x = sample, fill = sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Number of Cells per Condition")
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig1-Cell_Counts.png", width = 7, height = 9, bg = "transparent", dpi = 1200)
ls()

# Visualize the number UMIs/transcripts per cell
metadata %>%
  dplyr::filter(nUMI > 0) %>%
  ggplot(aes(color = sample, x = nUMI, fill = sample)) + 
  geom_density(alpha = 0.9) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("log 10 Cell Density") +
  geom_vline(xintercept = 15000)
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig2-UMIs.png", width = 7, height = 9, bg = "transparent", dpi = 1200)

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  dplyr::filter(nGene > 0) %>%
  ggplot(aes(color = sample, x = nGene, fill= sample)) + 
  geom_density(alpha = 0.9) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 600)
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig3-Genes_Per_Cell.png", width = 7, height = 9, bg = "transparent", dpi = 1200)

# Visualize the overall complexity of gene expression (novelty score)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.9) +
  theme_classic() +
  geom_vline(xintercept = 0.80)
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig4-Complexity.png", width = 7, height = 9, bg = "transparent", dpi = 1200)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  dplyr::filter(mitoRatio > 0) %>%
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.4) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.12)
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig5-Mitochondrial_Genes.png", width = 7, height = 9, bg = "transparent", dpi = 1200)

# Visualize the correlation between genes detected and number of UMIs
metadata %>% 
  dplyr::filter(nUMI > 0, nGene > 0) %>%
  ggplot(aes(x=nUMI, y=nGene)) + 
  geom_point(aes(color=mitoRatio)) +
  scale_colour_gradient(low = "Red", high = "black") +
  stat_smooth(method="loess", span=0.2, se=FALSE, aes(group=1)) +  
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 2000) +
  geom_hline(yintercept = 1000) +
  facet_wrap(~sample, scales = "free")
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig6-Correlation.png", width = 7, height = 9, bg = "transparent", dpi = 1200)

# Filter out low-quality cells using optimized thresholds
filtered_seurat <- subset(merged_sample, 
                          subset= (nUMI >= 2000) &  
                            (nGene >= 500) &  
                            (log10GenesPerUMI > 0.8) &  
                            (mitoRatio < 0.1))
filtered_seurat
ncol(merged_sample)
ncol(filtered_seurat)

merged_sample[["nCount_RNA"]] <- Matrix::colSums(merged_sample@assays$RNA@counts)
merged_sample[["nFeature_RNA"]] <- Matrix::colSums(merged_sample@assays$RNA@counts > 0)
merged_sample[["percent.mt"]] <- PercentageFeatureSet(merged_sample, pattern = "^MT-")

VlnPlot(merged_sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\VlnPlot.png", 
       plot = last_plot(), 
       width = 7, height = 9, 
       dpi = 300, bg = "transparent")

Layers(filtered_seurat)
counts <- GetAssayData(object = filtered_seurat, layer = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_sample <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Create .RData object to load at any time
save(filtered_sample, file="F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\sample_filtered.RData2.RData")

# Normalize the counts
sample_phase <- NormalizeData(filtered_sample)

# Load cell cycle markers
cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv") #------------Human
#cell_cycle_genes2<- read.csv(text = "cell_cycle_Homo_sapiens.csv")
cell_cycle_genes <- read.csv(text = cc_file)
# Connect to AnnotationHub
ah <- AnnotationHub()
#removeCache(ah)
# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)
# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

# Perform cell cycle scoring, # Score cells for cell cycle
sample_phase <- CellCycleScoring(sample_phase,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells 
View(sample_phase@meta.data) 

# Identify the most variable genes
sample_phase <- FindVariableFeatures(sample_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)
View(sample_phase@meta.data) 
length(VariableFeatures(sample_phase))

# Scale the counts
sample_phase <- ScaleData(sample_phase)

# Perform PCA
sample_phase <- RunPCA(sample_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(sample_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig7-PCA colored by cell cycle phase9246.png", 
       width = 7, height = 9, bg = "transparent", dpi = 1200)
options(future.globals.maxSize = 4000 * 1024^2)
save(sample_phase, file="F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\cell_cycle_HUM_RData3.RData")
merged_sample1 <- filtered_sample
# -----------------------------------------------------------------------------------
# Sample 1 => Split object by condition to perform cell cycle scoring and SCT on all samples

split_sample <- SplitObject(merged_sample1, split.by = "sample")
split_sample <- split_sample[c("Pre", "Post")]
for (i in names(split_sample)) {
  split_sample[[i]] <- NormalizeData(split_sample[[i]], verbose = TRUE)
  split_sample[[i]] <- CellCycleScoring(split_sample[[i]], g2m.features = g2m_genes, s.features = s_genes)
  split_sample[[i]] <- SCTransform(split_sample[[i]], variable.features.n = 2000, vars.to.regress = c("mitoRatio"))
}

# Sample 1 => Check which assays are stored in objects
split_sample1 <- split_sample
split_sample1$Pre@assays
split_sample1$Post@assays

memory.limit(size=100000)
# ------------------------------------------------------------
# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = c(split_sample1), 
                                            nfeatures = 3000) 
# Prepare the SCT list object for integration
split_sample <- PrepSCTIntegration(object.list = c(split_sample1), 
                                   anchor.features = integ_features)

# Find best buddies - can take a while to run  # 13/08/2025
integ_anchors <- FindIntegrationAnchors(object.list = split_sample, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
saveRDS(integ_anchors, "F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\integ_anchors.rds")
integ_anchors <- readRDS("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\integ_anchors.rds")
# Integrate across conditions
sample_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
# ------------------------------------------------------------------------------------
# Save integrated object
saveRDS(sample_integrated,"F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\integrated_sample.rds")
sample_integrated <- readRDS("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\integrated_sample.rds")

rm(list = c(paste0("seurat_JCML", 1:22), 
            paste0("seurat_sample", 1:22), 
            paste0("split_sample", 1:4), 
            "Post_sample", "Pre_sample"))

memory.limit(size=200000)

# Run PCA
sample_integrated <- RunPCA(object = sample_integrated)

# Plot PCA
PCAPlot(sample_integrated,
        split.by = "sample")
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig15-PCA-integrated_sample1.png", width=7, height=9, bg = "transparent", dpi = 1200)

# Run UMAP
sample_integrated <- RunUMAP(sample_integrated, 
                             dims = 1:40,
                             reduction = "pca")
# Plot UMAP   
DimPlot(sample_integrated) 
DimPlot(sample_integrated,
        split.by = "sample")
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig16-UMAP-integrated_sample.png", width=7, height=9, bg = "transparent", dpi = 1200)
#dev.off()

# Open a PNG graphics device
png("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig17-DimHeatmap-integrated_sample.png", 
    width = 7, height = 9, units = "in", res = 1200, bg = "transparent")

# Render the Heatmap (print it if needed to ensure it displays correctly)
DimHeatmap(sample_integrated, dims = 1:9, cells = 500, balanced = TRUE)
png("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig18-DimHeatmap-integrated_sample.png", 
    width = 7, height = 9, units = "in", res = 1200, bg = "transparent")
# Close the graphics device to save the image
dev.off()

# Printing out the most variable genes driving PCs
print(x = sample_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

# Plot the elbow plot
ElbowPlot(object = sample_integrated, 
          ndims = 40)
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig19-elbow.png", 
       width = 7, height = 9, bg = "white", dpi = 1200)

# Determine the K-nearest neighbor graph   
# Seurat uses a graph-based clustering approach, which embeds cells in a graph structure, using a K-nearest neighbor (KNN) graph (by default), with edges drawn between cells with similar gene expression patterns. Then, it attempts to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’ [Seurat - Guided Clustering Tutorial].
sample_integrated <- FindNeighbors(object = sample_integrated, 
                                   dims = 1:40)
# Determine the clusters for various resolutions
sample_integrated <- FindClusters(object = sample_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
sample_integrated@meta.data %>% 
  View()

# Assign identity of clusters
Idents(object = sample_integrated) <- "integrated_snn_res.0.8"

# Plot the UMAP
DimPlot(sample_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig20-integrated_uMAP1.png", width=7, height=9, bg = "transparent", dpi = 1200)

DimPlot(sample_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6, 
        split.by = "sample")
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig21-integrated_uMAP1_split.png", width=7, height=9, bg = "transparent", dpi = 1200)

# Assign identity of clusters
Idents(object = sample_integrated) <- "integrated_snn_res.0.4"

# Plot the UMAP
DimPlot(sample_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig22-Clus_plot.0.4.png", width = 12, height = 8, dpi = 600)
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig22-Clus_plot.0.4.pdf")

DimPlot(sample_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        split.by = "sample")
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig23-Clus_plot_2.0.4.png", width = 12, height = 8, dpi = 1200)

# Assign identity of clusters
Idents(object = sample_integrated) <- "integrated_snn_res.1.4"
# Plot the UMAP
DimPlot(sample_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig24-integrated_uMAP3.1.4.png", width=7, height=9, bg = "transparent", dpi = 1200)

DimPlot(sample_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6, 
        split.by = "sample")
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig25-integrated_uMAP3.1.4_split.png", width=7, height=9, bg = "transparent", dpi = 1200)
dev.off()

# Assign identity of clusters
Idents(object = sample_integrated) <- "integrated_snn_res.0.4"
# Plot the UMAP
DimPlot(sample_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig26-integrated_uMAP2.0.4.png", width=7, height=9, bg = "transparent", dpi = 1200)

DimPlot(sample_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6, 
        split.by = "sample")
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig27-integrated_uMAP2_split.0.4.png", width=7, height=9, bg = "transparent", dpi = 1200)
dev.off()

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(sample_integrated, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)
View(n_cells)
write.csv(n_cells, 
          "F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\n_cells_data.csv", 
          row.names = FALSE)
# UMAP of cells in each cluster by sample
DimPlot(sample_integrated, 
        label = TRUE, 
        split.by = "sample")  

ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig28-Clus2_plot.0.4.pdf")
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig28-Clus2_plot_2.0.4.png", width = 12, height = 8, dpi = 1200)

# Explore whether clusters segregate by cell cycle phase
DimPlot(sample_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig29-cell_cycle_phase.png", width = 12, height = 8, dpi = 1200)

DimPlot(sample_integrated,
        label = TRUE, 
        split.by = "sample")  + NoLegend()
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig30-cluster1.png", width = 12, height = 8, dpi = 1200)

# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")
FeaturePlot(sample_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig31-FeaturePlot.0.4.png", width = 12, height = 8, dpi = 1200)

# Determine metrics to plot present in seurat_integrated@meta.data # split by sample
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")
FeaturePlot(sample_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            split.by = "sample")
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig32-FeaturePlot-Sample.0.4.png", width = 12, height = 8, dpi = 1200)
columns <- c(paste0("PC_", 1:16),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(sample_integrated, 
                     vars = columns)
# Extract the UMAP coordinates for the first 10 cells
# sample_integrated@reductions$umap@cell.embeddings[1:10, 1:2]

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(sample_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
#install.packages("purrr")
library(purrr)
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig33-PCs_plot.pdf")
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig33-PCs_plot.png", width = 12, height = 8, dpi = 1200)

#We could look back at our genes driving this PC to get an idea of what the cell types might be:
# Examine PCA results
print(sample_integrated[["pca"]], dims = 1:16, nfeatures = 10)
DimPlot(object = sample_integrated, 
        reduction = "umap", 
        label = TRUE,) + NoLegend()
DimPlot(object = sample_integrated, 
        reduction = "umap", 
        label = TRUE,
        split.by = "sample") 
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig34-PCA results.png", width = 12, height = 8, dpi = 1200)

# Select the RNA counts slot to be the default assay
DefaultAssay(sample_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
sample_integrated <- NormalizeData(sample_integrated, verbose = FALSE)

pre_sample <- subset(sample_integrated, subset = sample == "Pre")
post_sample <- subset(sample_integrated, subset = sample == "Post")
features_to_plot = c("ANKRD30A", "AR", "CD24",	"CD44",	"CDH1",	"CK18",	"Cytokeratin-19",	"EMA",	"EPCAM",	"ERBB2",	"ESA",	"ESR1",	"KRT15",	"KRT23",	"KRT81",	"MKI67",	"MUC",	"MUC1",	"PGR",	"SYTL2")
pre_plot <- DotPlot(pre_sample, 
                    features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Pre-treatment")
post_plot <- DotPlot(post_sample, 
                     features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Post-treatment")
combined_plot <- pre_plot + post_plot
combined_plot
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig35-Luminal Epithelial Cell_combined.png", width = 12, height = 15, dpi = 1200)
memory.limit(size=200000)

# DotPlot for multiple markers # Cytotoxic T Cell
pre_sample <- subset(sample_integrated, subset = sample == "Pre")
post_sample <- subset(sample_integrated, subset = sample == "Post")
features_to_plot = c("CD2","CD25","CD28","CD3D","CD3E","CD3G","CD8","CD8A","GZMA","GZMB","LAPTM5","NKG7","PRF1","PTPRC","VIM")
pre_plot <- DotPlot(pre_sample, 
                    features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Pre-treatment")
post_plot <- DotPlot(post_sample, 
                     features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Post-treatment")
combined_plot <- pre_plot + post_plot
combined_plot
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig36-Cytotoxic T Cell_combined.png", width = 12, height = 15, dpi = 1200)

# DotPlot for multiple markers # Endothelial Cell
pre_sample <- subset(sample_integrated, subset = sample == "Pre")
post_sample <- subset(sample_integrated, subset = sample == "Post")
features_to_plot = c("CD31","CD34","CD36","CDH5","CLEC14A","EMCN","IL3RA","KDR","MCAM","PECAM1","PTPRB","SERPINH1","SOX18","VIM","VWF")
pre_plot <- DotPlot(pre_sample, 
                    features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Pre-treatment")
post_plot <- DotPlot(post_sample, 
                     features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Post-treatment")
combined_plot <- pre_plot + post_plot
combined_plot
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig37-Endothelial Cell_combined.png", width = 12, height = 15, dpi = 1200)

# DotPlot for multiple markers # Luminal Progenitor Cell
pre_sample <- subset(sample_integrated, subset = sample == "Pre")
post_sample <- subset(sample_integrated, subset = sample == "Post")
features_to_plot = c("ALDH1A3","c-Kit","CD24","CD44","CD49f","CK5","EpCAM","EZH2","GABRP","KIT","KRT19","PROM1","SLPI","TNFRSF11A")
pre_plot <- DotPlot(pre_sample, 
                    features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Pre-treatment")
post_plot <- DotPlot(post_sample, 
                     features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Post-treatment")
combined_plot <- pre_plot + post_plot
combined_plot
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig38-Luminal Progenitor Cell_combined.png", width = 12, height = 15, dpi = 1200)

# DotPlot for multiple markers # T Cell
pre_sample <- subset(sample_integrated, subset = sample == "Pre")
post_sample <- subset(sample_integrated, subset = sample == "Post")
features_to_plot = c("41BB","CCL5","CD2","CD28","CD3","CD3D","CD3E","CD3G","CD4","CD5","CD7","CD8","CD8A","CFD","ETS1","FOXP3","GZMK","ICOC","IFIT1","LAPTM5","PTPRC","VIM")
pre_plot <- DotPlot(pre_sample, 
                    features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Pre-treatment")
post_plot <- DotPlot(post_sample, 
                     features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Post-treatment")
combined_plot <- pre_plot + post_plot
combined_plot
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig39-T Cell_combined.png", width = 12, height = 15, dpi = 1200)

# DotPlot for multiple markers # Natural Killer Cell
pre_sample <- subset(sample_integrated, subset = sample == "Pre")
post_sample <- subset(sample_integrated, subset = sample == "Post")
features_to_plot = c("CCR7","CD3D","CD3G","CD56","CD8A","GNLY","HAVR2C","KIR2DL1","KIR2DL3","KIR2DL4","KIR2DS4","KIR3DL1","KIR3DL2","KIR3DL3","KLRB1","KLRC1","KLRD1","KLRF1","NCAM1","NCR1","NKG7","XCL1")
pre_plot <- DotPlot(pre_sample, 
                    features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Pre-treatment")
post_plot <- DotPlot(post_sample, 
                     features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Post-treatment")
combined_plot <- pre_plot + post_plot
combined_plot
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig40-Natural Killer Cell_combined.png", width = 12, height = 15, dpi = 1200)

# DotPlot for multiple markers # Fibroblast
pre_sample <- subset(sample_integrated, subset = sample == "Pre")
post_sample <- subset(sample_integrated, subset = sample == "Post")
features_to_plot = c("ACTA2","CK18","CK19","CK8","COL1A1","COL1A2","CXCL12","DCN","EpCAM","ITGB1","MMP11","PDGFRA","PDGFRB","PN1","POSTN","PTN","S100A4","SERPINH1","SPARC","TAGLN","THY1","VIM")
pre_plot <- DotPlot(pre_sample, 
                    features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Pre-treatment")
post_plot <- DotPlot(post_sample, 
                     features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Post-treatment")
combined_plot <- pre_plot + post_plot
combined_plot
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig41-Fibroblast_combined.png", width = 12, height = 15, dpi = 1200)

# DotPlot for multiple markers # Macrophages
pre_sample <- subset(sample_integrated, subset = sample == "Pre")
post_sample <- subset(sample_integrated, subset = sample == "Post")
features_to_plot = c("CD11b","CD11c","CD14","CD163","CD209","CD33","CD4","CD68","CD74","CD86","CDH5","CLEC5A","CSF1","CSF1R","FCGR3A","HLA-DPA1","IL1B","ITGAM","ITGAX","LAPTM5","LYZ","MSR1","PTPRC","VIM")
pre_plot <- DotPlot(pre_sample, 
                    features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Pre-treatment")
post_plot <- DotPlot(post_sample, 
                     features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Post-treatment")
combined_plot <- pre_plot + post_plot
combined_plot
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig42-Macrophage_combined.png", width = 12, height = 15, dpi = 1200)

# DotPlot for multiple markers # B Cell
pre_sample <- subset(sample_integrated, subset = sample == "Pre")
post_sample <- subset(sample_integrated, subset = sample == "Post")
features_to_plot = c("BANK1","CCDC50","CD14","CD19","CD20","CD27","CD38","CD40","CD74","CD79A","IGHD","IGHM","IGJ","IGKC","IGLL1","IGLL5","LAPTM5","MS4A1","PTPRC","SEL1L3","SSR4","TCF4","TMP1","VIM")
pre_plot <- DotPlot(pre_sample, 
                    features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Pre-treatment")
post_plot <- DotPlot(post_sample, 
                     features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Post-treatment")
combined_plot <- pre_plot + post_plot
combined_plot
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig43-B Cell_combined.png", width = 12, height = 15, dpi = 1200)

# DotPlot for multiple markers # Exhausted CD8+ T Cell
pre_sample <- subset(sample_integrated, subset = sample == "Pre")
post_sample <- subset(sample_integrated, subset = sample == "Post")
features_to_plot = c("ALOX5AP","BST2","CCL3","CD82","CXCL13","ENTPD1","GNLY","GZMB","HAVCR2","HLA-DQA2","HLA-DRA","HLA-DRB1","HLA-DRB5","IFI27","IFI44","IFI44L","IFI6","IFIT3","IFITM3","ISG15","KRT86","LY6E","MX1","OAS1","PRF1")
pre_plot <- DotPlot(pre_sample, 
                    features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Pre-treatment")
post_plot <- DotPlot(post_sample, 
                     features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Post-treatment")
combined_plot <- pre_plot + post_plot
combined_plot
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig44-Exhusted CD8+ T Cell_combined.png", width = 12, height = 15, dpi = 1200)

# DotPlot for multiple markers # Stem Cell
pre_sample <- subset(sample_integrated, subset = sample == "Pre")
post_sample <- subset(sample_integrated, subset = sample == "Post")
features_to_plot = c("ABCG2","ALDH1","ALDH1A1","ALDH1A2","ALDH1A3","ALDH1B1","ALDH2","CD133","CD24","CD29","CD44","CD49f","CD61","CD70","CD90","CDH5","CXCR4","EGFR","EpCAM","ID2","KLF4","LEPR","LGR4","LGR5","LGR6","MET","NANOG","Notch2","4-Oct","ProC-R","Procr")
pre_plot <- DotPlot(pre_sample, 
                    features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Pre-treatment")
post_plot <- DotPlot(post_sample, 
                     features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Post-treatment")
combined_plot <- pre_plot + post_plot
combined_plot
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig45-Stem Cell_combined.png", width = 12, height = 15, dpi = 1200)

# DotPlot for multiple markers # MesenchymalCell
pre_sample <- subset(sample_integrated, subset = sample == "Pre")
post_sample <- subset(sample_integrated, subset = sample == "Post")
features_to_plot = c("Cav-1","CD248","CD44","CD59","CD90","CDH2","collagens","CTSC","ECM1","EGFR","ETS1","FAK","FAP","FBLN5","FBN1","Fibronectin","fibronectins","FN1","MCAM","MMP2","MSN","N-cadherin","PDGFRB","PLAT","proteases","Rac1","Slug","SNAI1","Snail","SPARC","Src","TGFB1","TGFBR2","THY1","TIMP1","Twist","TWIST1","VIM","Vimentin","Zeb-1","ZEB1","ZEB2")
pre_plot <- DotPlot(pre_sample, 
                    features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Pre-treatment")
post_plot <- DotPlot(post_sample, 
                     features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Post-treatment")
combined_plot <- pre_plot + post_plot
combined_plot
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig46-MesenchymalCell_combined.png", width = 12, height = 15, dpi = 1200)

# DotPlot for multiple markers # EpithelialCell
pre_sample <- subset(sample_integrated, subset = sample == "Pre")
post_sample <- subset(sample_integrated, subset = sample == "Post")
features_to_plot = c("ACTA2","Aire","ALDH1A3","Avil","Ccl21a","CD117","CD166","CDH1","CDH1/E-cadherin","CK","CK14","CK18","CLDN12","CXADR","Cytokeratin","Cytokeratin-18","cytokeratins","Dclk1","desmoplakins","E-cadherin","EMA","EPCAM","ESR1","Fezf2","GCDFP15","Itga6","Ivl","Krt1","Krt10","KRT14","KRT18","KRT5","KRT8","Ly6a","MCAM","OCLN","P63","S100A9","SDC1","SMA","Spink5","ZEB1","ZO-1")
pre_plot <- DotPlot(pre_sample, 
                    features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Pre-treatment")
post_plot <- DotPlot(post_sample, 
                     features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Post-treatment")
combined_plot <- pre_plot + post_plot
combined_plot
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig47-EpithelialCell_combined.png", width = 12, height = 15, dpi = 1200)

# DotPlot for multiple markers # Cancer Stem Cell
pre_sample <- subset(sample_integrated, subset = sample == "Pre")
post_sample <- subset(sample_integrated, subset = sample == "Post")

features_to_plot <- c("ABCG2","ALCAM","ALDH","ALDH1","ALDH1A1","ALDH1A3","ALDH2","ALDH6A1","BMI1","c-Kit","CD10","CD117","CD133","CD15","CD166","CD200","CD24","CD29","CD34","CD44","CD49f","CD61","CD90","Connexin 43","CXCR4","EPCAM","ESA","ESR1","Feline sarcoma-related tyrosine kinase","GD2","GLI-1","GLI1","GTAT3","Her-2","HER2","LGR5","LSD1","MYC","Nanog","Nectin-4","Nestin","OCT3/4","4-Oct","P63","PGR","PODXL-1","POU5F1","PROCR","SOX2","SOX9","THY1","VEGF","VIM")
pre_plot <- DotPlot(pre_sample, 
                    features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Pre-treatment")
post_plot <- DotPlot(post_sample, 
                     features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Post-treatment")
combined_plot <- pre_plot + post_plot
combined_plot
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig48-Cancer_Stem_Cell_combined.png", plot = combined_plot, width = 24, height = 15, dpi = 1200)

# DotPlot for multiple markers # Cancer Cell
pre_sample <- subset(sample_integrated, subset = sample == "Pre")
post_sample <- subset(sample_integrated, subset = sample == "Post")

features_to_plot = c("ADH1B","AKR1A1","ALPI","ARHGAP10","B3GNT2","BCAM","BRCA1","BRCA2","C1GALT1C1","CA153","CAMK1D","CD36","CDH1","CDH5","CHST15","CK19","CPNE1","CPS-II","CRYBB2","CSC","CTNNA1","CYP3A4","CYP3A5","CYP3A7","DOCK9","EDNRA","EMT","ENG","EPCAM","ER","ERBB2","ESR1","FAM177A1","FASLG (soluble)","GATA3","GCDFP-15","GNA13","GOLM1","HER-2","HER2","IGFIR (soluble)","IL3RA","IR","IRF9","ISG15","ISLR2","JAG1","JAK1","KIN","KLRF1","KRT18","KRT19","KRT7","KRT8","LIFR (soluble)","MCCS","MET","MKI67","NOTCH1","OAS1","panCK","PBCS","PD-1","pKARMA","PRDM1","PSD","RIC8A","RSPO3","SELE","SELP","SEMA6A","STAT1","STAT2","STOM","SULF2","THSD1","TIE1 (soluble)","TLL1","TMPRSS11D","TNS2","TPST2","TRPS1","TTF-1","UKBGS","VCAM1","VEGFR2/KDR","VEGFR3/FLT4","VIM")
pre_plot <- DotPlot(pre_sample, 
                    features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Pre-treatment")
post_plot <- DotPlot(post_sample, 
                     features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Post-treatment")
combined_plot <- pre_plot + post_plot
combined_plot
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig49-Cancer Cell_combined.png", width = 12, height = 15, dpi = 1200)

#Identification of conserved markers in all conditions
DefaultAssay(sample_integrated) <- "RNA"

cluster0_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                   ident.1 = 0,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
View(annotations)
# Combine markers with gene descriptions 
# This is where rownames_to_column comes from
cluster0_ann_markers <- cluster0_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster0_ann_markers)

# comparing cluster 1 with other clusters 
cluster1_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                   ident.1 = 1,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster1_ann_markers <- cluster1_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster1_ann_markers)

# comparing cluster 2 with other clusters 
cluster2_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                   ident.1 = 2,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster2_ann_markers <- cluster2_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster2_ann_markers)

# comparing cluster 3 with other clusters 
cluster3_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                   ident.1 = 3,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster3_ann_markers <- cluster3_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster3_ann_markers)

# comparing cluster 4 with other clusters 
cluster4_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                   ident.1 = 4,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster4_ann_markers <- cluster4_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster4_ann_markers)

# comparing cluster 5 with other clusters 
cluster5_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                   ident.1 = 5,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster5_ann_markers <- cluster5_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster5_ann_markers)

# comparing cluster 6 with other clusters 
cluster6_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                   ident.1 = 6,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster6_ann_markers <- cluster6_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster6_ann_markers)

# comparing cluster 7 with other clusters 
cluster7_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                   ident.1 = 7,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster7_ann_markers <- cluster7_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster7_ann_markers)

# comparing cluster 8 with other clusters 
cluster8_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                   ident.1 = 8,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster8_ann_markers <- cluster8_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster8_ann_markers)

# comparing cluster 9 with other clusters 
cluster9_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                   ident.1 = 9,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster9_ann_markers <- cluster9_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster9_ann_markers)

# comparing cluster 10 with other clusters 
cluster10_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                    ident.1 = 10,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster10_ann_markers <- cluster10_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster10_ann_markers)

# comparing cluster 11 with other clusters 
cluster11_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                    ident.1 = 11,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster11_ann_markers <- cluster11_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster11_ann_markers)

# comparing cluster 12 with other clusters 
cluster12_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                    ident.1 = 12,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster12_ann_markers <- cluster12_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster12_ann_markers)

# comparing cluster 13 with other clusters 
cluster13_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                    ident.1 = 13,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster13_ann_markers <- cluster13_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster13_ann_markers)

# comparing cluster 14 with other clusters 
cluster14_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                    ident.1 = 14,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster14_ann_markers <- cluster14_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster14_ann_markers)

# comparing cluster 15 with other clusters 
cluster15_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                    ident.1 = 15,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster15_ann_markers <- cluster15_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster15_ann_markers)

# comparing cluster 16 with other clusters 
cluster16_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                    ident.1 = 16,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster16_ann_markers <- cluster16_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster16_ann_markers)

# comparing cluster 17 with other clusters 
cluster17_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                    ident.1 = 17,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster17_ann_markers <- cluster17_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster17_ann_markers)

# comparing cluster 18 with other clusters 
cluster18_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                    ident.1 = 18,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster18_ann_markers <- cluster18_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster18_ann_markers)

# comparing cluster 19 with other clusters 
cluster19_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                    ident.1 = 19,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster19_ann_markers <- cluster19_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster19_ann_markers)

# comparing cluster 20 with other clusters 
cluster20_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                    ident.1 = 20,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster20_ann_markers <- cluster20_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster20_ann_markers)

# comparing cluster 21 with other clusters 
cluster21_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                    ident.1 = 21,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster21_ann_markers <- cluster21_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster21_ann_markers)

# comparing cluster 22 with other clusters 
cluster22_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                    ident.1 = 22,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster22_ann_markers <- cluster22_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster22_ann_markers)

# comparing cluster 23 with other clusters 
cluster23_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                    ident.1 = 23,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster23_ann_markers <- cluster23_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster23_ann_markers)

# comparing cluster 24 with other clusters 
cluster24_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                    ident.1 = 24,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster24_ann_markers <- cluster24_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster24_ann_markers)

# comparing cluster 25 with other clusters 
cluster25_conserved_markers <- FindConservedMarkers(sample_integrated,
                                                    ident.1 = 25,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster25_ann_markers <- cluster25_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster25_ann_markers)


# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(sample_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       logfc.threshold = 0.25,
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}
conserved_markers <- map_dfr(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25), get_conserved)
saveRDS(conserved_markers, file = "F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\conserved_markers.rds")

#all_conserved_markers <- bind_rows(conserved_markers, cluster3_conserved_markers)
all_conserved_markers <- bind_rows(conserved_markers)

# Extract top 10 markers per cluster

top10 <- all_conserved_markers %>% 
  mutate(avg_fc = (Pre_avg_log2FC + Post_avg_log2FC) / 2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, wt = avg_fc)

colnames(all_conserved_markers)

# Visualize top 10 markers per cluster
View(top10)
# Ensure that top10 is a data frame and contains the data you want to export
#install.packages("readr")

write_csv(top10,"F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\TOP10.csv")
library(dplyr)

data_clusters <- read.csv("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\TOP10.csv")
output_file <- "F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\matched_results_combined.csv"

clusters <- data_clusters %>%
  group_by(cluster_id) %>%
  summarise(Genes = list(gene))

data2 <- read_excel("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Azimouth-Marker.xlsx",
                    sheet = "Azimouth-Marker", 
                    skip = 1)   
names(data2)  
data3 <- read_excel("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\CellMarker.xlsx")  
results <- list()
for (i in 1:nrow(clusters)) {
  cluster_id <- clusters$cluster_id[i]
  genes <- clusters$Genes[[i]]
  search_array <- genes
  matching_columns_data2 <- sapply(data2, function(column) {
    any(column %in% search_array)
  })
  matched_columns_data2 <- names(matching_columns_data2[matching_columns_data2 == TRUE])
  matching_columns_data3 <- sapply(data3, function(column) {
    any(column %in% search_array)
  })
  matched_columns_data3 <- names(matching_columns_data3[matching_columns_data3 == TRUE])
  matched_columns_combined <- c(
    if (length(matched_columns_data2) > 0) matched_columns_data2 else "None from data2",
    if (length(matched_columns_data3) > 0) matched_columns_data3 else "None from data3"
  )
  results[[i]] <- data.frame(
    ClusterID = cluster_id,
    Genes = paste(search_array, collapse = ", "),
    MatchedColumns = paste(matched_columns_combined, collapse = ", ")
  )
}

final_results <- do.call(rbind, results)
write.csv(final_results, output_file, row.names = FALSE)
print(final_results)
#-----------------------
# Rename all identities
sample_integrated <- RenameIdents(object = sample_integrated, 
                                  "0" = "T Cells",
                                  "1" = "Natural Killer Cells",
                                  "2" = "Natural Killer Cells",
                                  "3" = "T Cells",
                                  "4" = "Dendritic Cells",
                                  "5" = "Monocytes",
                                  "6" = "Dendritic Cells",
                                  "7" = "Natural Killer Cells",
                                  "8" = "Natural Killer Cells",
                                  "9" = "Monocytes",
                                  "10"= "T Cells",
                                  "11" = "Natural Killer Cells",
                                  "12" = "Monocytes",
                                  "13"= "T Cells",
                                  "14" = "Natural Killer Cells",
                                  "15" = "T Cells",
                                  "16"= "Natural Killer Cells",
                                  "17" = "Macrophages",
                                  "18" = "T Cells",
                                  "19" = "Natural Killer Cells",
                                  "20" = "Natural Killer Cells",
                                  "21" = "Natural Killer Cells",
                                  "22" = "T Cells",
                                  "23" = "Dendritic Cells",
                                  "24" = "Plasmablasts",
                                  "25" = "Natural Killer Cells"
                                  
)

saveRDS(sample_integrated, "F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\integrated_sample_Finalcluster.FinaSave.rds")
sample_integrated <- readRDS("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\integrated_sample_Finalcluster.FinaSave.rds")
cat("final cell count:", ncol(sample_integrated), "\n")
p <- DimPlot(
  object = sample_integrated,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  repel = TRUE,
  pt.size = 0.5
) +
  theme_classic(base_size = 16) +  
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.ticks = element_line(size = 0.3),
    panel.grid = element_blank(),                 
    panel.border = element_blank(),               
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(x = "UMAP 1", y = "UMAP 2")
print(p)
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig51-UMAP-CellTypes-Last.png",
       plot = p, width = 12, height = 15, dpi = 1200, bg = "white")

ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\Fig51-UMAP-CellTypes.pdf",
       plot = p, width = 12, height = 15, dpi = 1200, bg = "white")

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(sample_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       logfc.threshold = 0.25,
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    
    cbind(cluster_id = cluster, .)
}

Idents(sample_integrated) %>% unique()

Idents(sample_integrated) <- str_trim(Idents(sample_integrated))
Idents(sample_integrated) %>% unique()
conserved_markers <- map_dfr(c(
  "Natural Killer Cells",
  "T Cells",
  "Monocytes",
  "Macrophages",
  "Dendritic Cells",
  "Plasmablasts"), get_conserved)

# انتخاب ژن‌های مهم
top_pre <- conserved_markers %>%
  group_by(cluster_id) %>%
  slice_max(order_by = Pre_avg_log2FC, n = 5, with_ties = FALSE) %>%
  dplyr::select(cluster_id, gene, avg_log2FC = Pre_avg_log2FC) %>%
  mutate(condition = "Pre")

top_post <- conserved_markers %>%
  group_by(cluster_id) %>%
  slice_max(order_by = Post_avg_log2FC, n = 5, with_ties = FALSE) %>%
  dplyr::select(cluster_id, gene, avg_log2FC = Post_avg_log2FC) %>%
  mutate(condition = "Post")

# ادغام و مرتب‌سازی
top_combined <- bind_rows(top_pre, top_post) %>%
  mutate(gene = fct_reorder(gene, avg_log2FC))

# پلات نهایی
p <- ggplot(top_combined, aes(x = gene, y = avg_log2FC, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  coord_flip() +
  facet_wrap(~ cluster_id, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Pre" = "#66C2A5", "Post" = "#FC8D62")) +
  labs(title = "Top 5 Conserved Markers per Cluster",
       x = NULL, y = "avg_log2FC") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 8),
    axis.title = element_text(size = 10),
    strip.text = element_text(face = "bold", size = 9),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "top"
  )

ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\top_conserved_markers_clean.tiff",
       plot = p, width = 14, height = 10, dpi = 600, compression = "lzw", bg = "white")

ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\top_conserved_markers_clean.pdf",
       plot = p, width = 14, height = 10, dpi = 600)

#all_conserved_markers <- bind_rows(conserved_markers, cluster3_conserved_markers)
all_conserved_markers <- bind_rows(conserved_markers)

top10 <- all_conserved_markers %>% 
  mutate(avg_fc = (Pre_avg_log2FC + Post_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)
View(top10)

write_csv(top10,"F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\TOP10-Clusters-final-Blood-Pacli&Atezo-plus Azi.csv")

library(ggplot2)
library(readr)
library(dplyr)
top10 <- read_csv("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\TOP10-Clusters-final-Blood-Pacli&Atezo-plus Azi.csv")
top10 <- top10 %>%
  mutate(avg_fc = (Pre_avg_log2FC + Post_avg_log2FC) / 2)

p <- ggplot(top10, aes(x = avg_fc, y = reorder(gene, avg_fc), fill = cluster_id)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_wrap(~cluster_id, scales = "free_y", ncol = 2) +
  scale_fill_brewer(palette = "Paired") +
  labs(
    title = "Top 10 Marker Genes per Cluster",
    x = "Average log2 Fold Change",
    y = "Gene"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
print(p)
ggsave("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\top10_markers_per_cluster_ggplot.png",
       plot = p, width = 12, height = 10, dpi = 600)




