library(Seurat)
library(dplyr)
library(readr)
library(e1071)
library(caret)
library(rpart) 
library(ggplot2)
library(rpart.plot)
library(iml)
library(data.table)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(httr)
library(jsonlite)
library(readxl)
library(janitor)
library(reshape2)
library(randomForest)
library(reshape2)
library(xgboost)
library(Matrix)
library(openxlsx)
library(DALEX)
library(MLmetrics)
# ===============================================================
file_path <- "F:\\TOP-BIOMARKERS\\Blood_Pac&Atezo\\Top80_Biomarkers_WeightedScore.csv"
data <- read.csv(file_path, header = TRUE)
d <- data$gene
data_dir <- "D:/Thesis/Proposal/Thesis_Data/GSE169246/files"
names <- c("Pre_P001_b","Pre_P002_b","Pre_P004_b",
           "Pre_P005_b", "Pre_P007_b", "Pre_P019_b"
)  
sample_dirs <- list.dirs(data_dir, recursive = TRUE, full.names = TRUE)
sample_dirs <- sample_dirs[basename(sample_dirs) %in% names]
seurat_list <- list()
for (dir in sample_dirs) {
  sample_name <- basename(dir)
  counts <- Read10X(data.dir = dir)
  colnames(counts) <- paste(sample_name, colnames(counts), sep = "-")
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name)
  seurat_list[[sample_name]] <- seurat_obj@assays$RNA@counts
}
combined_counts <- do.call(cbind, seurat_list)
filtered_counts <- combined_counts[rownames(combined_counts) %in% d, ]
saveRDS(filtered_counts, "F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\filtered_counts_matrix.rds")
filtered_counts <- readRDS("F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\filtered_counts_matrix.rds")
cell_names <- colnames(combined_counts)
sample_ids <- sapply(strsplit(cell_names, "-"), `[`, 1)
unique_samples <- unique(sample_ids)
cell_labels <- data.frame(Cell = character(), Label = character(), stringsAsFactors = FALSE)
valid_labels <- c("progressive", "non-progressive")
for (sample in unique_samples) {
  repeat {
    label <- readline(prompt = paste0("Enter label for sample '", sample, "' (progressive / non-progressive): "))
    label <- trimws(tolower(label))  
    if (label %in% valid_labels) break
    cat("Invalid input. Please enter only 'progressive' or 'non-progressive'.\n")
  }
  cells <- cell_names[sample_ids == sample]
  sample_df <- data.frame(Cell = cells, Label = label, stringsAsFactors = FALSE)
  cell_labels <- rbind(cell_labels, sample_df)
}
saveRDS(cell_labels,
        "F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\cell_labels.rds")
cell_labels <- readRDS("F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\cell_labels.rds")
write.csv(cell_labels, "F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\cell_labels.csv", row.names = FALSE)
fwrite(cell_labels,
       file = "F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\cell_labels_fast.csv",
       row.names = FALSE)
# ====================================== GXBoots
filtered_counts <- readRDS("F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\filtered_counts_matrix.rds")
cell_labels <- readRDS("F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\cell_labels.rds")
common_cells <- intersect(colnames(filtered_counts), cell_labels$Cell)
filtered_counts <- filtered_counts[, common_cells]
cell_labels <- cell_labels[cell_labels$Cell %in% common_cells, ]
df <- as.data.frame(t(as.matrix(filtered_counts)))  
df$Cell <- rownames(df)
df <- merge(df, cell_labels, by = "Cell")
X_all <- df[, !(names(df) %in% c("Cell", "Label"))]  
y_all <- df$Label  
y_all <- gsub("-", "_", y_all)
preProc_all <- preProcess(X_all, method = c("center", "scale"))
X_all_scaled <- predict(preProc_all, X_all)
model_data <- data.frame(X_all_scaled, Label = factor(y_all, levels = c("non_progressive", "progressive")))
print(dim(model_data))
print(table(model_data$Label))
train_control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,  #Using twoClassSummary for ROC, Precision, and F1
  savePredictions = "final",
  verboseIter = TRUE
)
xgb_grid <- expand.grid(
  nrounds = 200,
  max_depth = c(3, 5, 7),
  eta = c(0.01, 0.05, 0.1),
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.8
)

# ============================ Training-XGBoost
xgb_model <- train(
  Label ~ .,
  data = model_data,
  method = "xgbTree",
  trControl = train_control,
  tuneGrid = xgb_grid,
  metric = "ROC"  
)
saveRDS(xgb_model, "F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\xgboost_cv200f.rds")
write.csv(xgb_model$results, "F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\xgboost_cv200f.csv", row.names = FALSE)
colnames(xgb_model$results)

# ========================= Extract the best model based on AUC
best_idx <- which.max(xgb_model$results$ROC)
best_row <- xgb_model$results[best_idx, ]
predictions <- predict(xgb_model, newdata = model_data)
cm <- confusionMatrix(predictions, model_data$Label)

# Extraction of Precision، F1-Score، Accuracy و Recall from diffusion matrix
precision <- cm$byClass["Pos Pred Value"]  # Precision
f1_score <- cm$byClass["F1"]  # F1-Score
accuracy <- cm$overall["Accuracy"]  # Accuracy
recall <- cm$byClass["Recall"]  # Recall
sink("F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\xgboost_best_model_summary.txt")
cat("========== XGBoost Best Model Summary ==========\n\n")
cat("Best Hyperparameters:\n")
print(xgb_model$bestTune)
cat(sprintf("\n Best ROC AUC: %.2f%%\n", best_row$ROC * 100))  # AUC
cat(sprintf(" Sensitivity: %.2f%%\n", best_row$Sens * 100))   # Sensitivity
cat(sprintf("️ Specificity: %.2f%%\n", best_row$Spec * 100))  
cat(sprintf(" Precision: %.2f%%\n", precision * 100))  # Precision
cat(sprintf(" Recall: %.2f%%\n", recall * 100))  # Recall
cat(sprintf(" F1-Score: %.2f%%\n", f1_score * 100))  # F1-Score
cat(sprintf(" Accuracy: %.2f%%\n", accuracy * 100))  # Accuracy
sink()
# ======================= Important Features
importance <- varImp(xgb_model, scale = FALSE)
top_features <- importance$importance
top_features <- top_features[order(-top_features$Overall), , drop = FALSE]
write.csv(top_features, "F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\xgboost_feature_importance100.csv")

png("F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\xgboost_feature_importance_plot100.png", width = 800, height = 600)
imp_df <- data.frame(Feature = rownames(top_features), Importance = top_features$Overall)
top_n <- head(imp_df, 20)
ggplot(top_n, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 Important Features in XGBoost", x = "Feature", y = "Importance") +
  theme_minimal()
dev.off()
top_features$Gene <- rownames(top_features)
save_gene_subset <- function(gene_list, filename) {
  df_subset <- model_data[, c(gene_list, "Label")]
  write.csv(df_subset, file = filename, row.names = FALSE)
}
top20 <- head(top_features$Gene, 20)
top50 <- head(top_features$Gene, 50)
top60 <- head(top_features$Gene, 60)
top70 <- head(top_features$Gene, 70)
top80 <- head(top_features$Gene, 80)
save_gene_subset(top20, "F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\data_top20-200f.csv")
save_gene_subset(top50, "F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\data_top50-200f.csv")
save_gene_subset(top60, "F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\data_top60-200f.csv")
save_gene_subset(top70, "F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\data_top70-200f.csv")
save_gene_subset(top80, "F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\data_top80-200f.csv")
d_mapped <- make.names(d)
initial200 <- model_data[, c(intersect(d_mapped, colnames(model_data)), "Label")]
write.csv(initial200, "F:\\ML\\ML_Blood_Pacli&Atezo\\GXBoots\\data_initial_79.csv", row.names = FALSE)
