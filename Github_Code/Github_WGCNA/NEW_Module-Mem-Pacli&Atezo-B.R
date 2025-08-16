library(Seurat)
library(Matrix)
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
# #------------------------------read blood, paclitaxel plus atezolizumab-------------------------------------
 folder_names <- c("Post_P001_b","Post_P005_b","Post_P019_b","Pre_P001_b",
                   "Pre_P002_b","Pre_P004_b","Pre_P005_b","Pre_P007_b", 
                   "Pre_P019_b","Prog_P002_b","Prog_P004_b","Prog_P007_b")  

#------------------------------------------------------------------------------------
seurat_list <- list()
for (i in seq_along(folder_names)) {
  folder_name <- folder_names[i]
  sample_path <- file.path("F:\\WGCNA\\Important_modules_Blood", folder_name)
  data <- Read10X(data.dir = sample_path)
  seurat_object <- CreateSeuratObject(counts = data, min.features = 100, project = folder_name)
  seurat_object$orig.ident <- folder_name
  seurat_object$sample <- tolower(strsplit(folder_name, "_")[[1]][1])
  seurat_list[[folder_name]] <- seurat_object
}
#---------merge
merged_samples <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),   # Take the folder names
  project = "Blood_samples"
)

samples <- folder_names
post_samples <- samples[grepl("Post|Prog", samples)]
trait <- ifelse(grepl("Prog", post_samples), "2", "1")

rm(seurat_list)
rm(seurat_object)
rm(data)
#---------pseudo bulk
raw_counts <- GetAssayData(merged_samples, layer = "counts")
head(raw_counts)
dim(raw_counts)            # gene number* cell number
cell_counts <- table(merged_samples$orig.ident) #number of cell
cell_counts
pseudo_bulk <- t(rowsum(t(as.matrix(raw_counts)), group = merged_samples$orig.ident))
head(pseudo_bulk)
dim(pseudo_bulk)   # gene number* Sample number
mean_counts <- sweep(pseudo_bulk, 2, cell_counts, FUN = "/")
head(mean_counts)
dim(mean_counts)
pseudo_bulk<-mean_counts
rm(mean_counts)
rm(raw_counts)
rm(combined)
#----------split samples based on pre and post for construct two wgcna network
sample_names <- colnames(pseudo_bulk)

sample_group <- ifelse(grepl("^Pre", sample_names), "pre", "post")
table(sample_group)

pre_samples <-  colnames(pseudo_bulk)[sample_group == "pre"]
post_samples <- colnames(pseudo_bulk)[sample_group == "post"]

datExpr_pre <- as.data.frame(t(pseudo_bulk[, pre_samples]))
datExpr_post <- as.data.frame(t(pseudo_bulk[, post_samples]))
#----------
var_pre  <- apply(datExpr_pre, 2, var)
var_post <- apply(datExpr_post, 2, var)
common_zero_var_genes <- names(var_pre)[var_pre == 0 & var_post == 0]
datExpr_pre  <- datExpr_pre[, !colnames(datExpr_pre) %in% common_zero_var_genes]
View(datExpr_pre)
datExpr_post <- datExpr_post[, !colnames(datExpr_post) %in% common_zero_var_genes]
#---------------------------
data1 <- read.csv("F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Figs\\TOP10-Clusters-final-Blood-Pacli&Atezo-plus Azi.csv")
output_dir <- "F:\\WGCNA\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)  # ایجاد مسیر در صورت عدم وجود
gene_info <- data1[, c("gene")]
head(gene_info)
genes_to_use <- intersect(gene_info, colnames(datExpr_pre))
datExpr_pre1 <- datExpr_pre[, genes_to_use]

num_modules <- ncol(MEs_pr)
head(MEs_pr)
module_names <- colnames(MEs_pr)
module_names

# select common genes
genes_to_use <- intersect(gene_info, colnames(datExpr_pre))
datExpr_pre1 <- datExpr_pre[, genes_to_use]

common_samples <- intersect(rownames(MEs_pr), rownames(datExpr_pre1))
MEs_mat <- MEs_pr[common_samples, ]
expr_mat <- datExpr_pre1[common_samples, ]

# calculate module membership (kME)
module_membership <- cor(expr_mat, MEs_mat, use = "pairwise.complete.obs")
write.csv(module_membership, file = "F:\\WGCNA\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Results\\ModuleMembership_NewData.csv")
weighted_cor_matrix <- module_membership * module_importance  
weighted_cor_matrix
write.csv(weighted_cor_matrix, file = "F:\\WGCNA\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Results\\Weighted_ModuleCorrelation_Pre_vs_Post.csv")
#weighted_cor_matrix 
max_membership_per_gene <- apply(weighted_cor_matrix, 1, max)

# output data frame
max_membership_df <- data.frame(
  gene = rownames(weighted_cor_matrix),
  Max_Membership = max_membership_per_gene
)
write.csv(max_membership_df, file = "F:\\Blood-Patient-Paclitaxel-plus atezolizumab-Analysis\\Results\\Max_Module_Membership_Per_Gene.csv", row.names = FALSE)


