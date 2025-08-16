library(readr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(extrafont)
library(clusterProfiler)
library(stringr)
# ===============================================
process_group <- function(group_name, deg_path, wgcna_path, go_path, kegg_path, cnv_genes = NULL) {
  deg <- read_csv(deg_path)
  wgcna <- read_csv(wgcna_path)
  go <- read_csv(go_path)
  kegg <- read_csv(kegg_path)
  deg_genes <- unique(deg$gene)
  wgcna_genes <- unique(wgcna$gene)
  go_genes <- unique(unlist(strsplit(go$geneID, "/")))
  kegg_genes <- unique(unlist(strsplit(kegg$geneID, "/")))
  cnv_genes <- if (!is.null(cnv_genes)) unique(cnv_genes$gene) else character(0)
  all_genes <- unique(c(deg_genes, wgcna_genes, go_genes, kegg_genes, cnv_genes))
  deg_info <- deg %>%
    dplyr::select(gene, cluster_id) %>%
    distinct()
  df <- data.frame(
    gene = all_genes,
    group = group_name,
    in_DEG = as.integer(all_genes %in% deg_genes),
    in_WGCNA = as.integer(all_genes %in% wgcna_genes),
    in_GO = as.integer(all_genes %in% go_genes),
    in_KEGG = as.integer(all_genes %in% kegg_genes),
    in_CNV = as.integer(all_genes %in% cnv_genes)
  )
  df <- df %>%
    left_join(wgcna, by = "gene") %>%
    left_join(deg_info, by = "gene") %>%
    mutate(evidence_score = in_DEG + in_WGCNA + in_GO + in_KEGG + in_CNV)
  return(df)
}
# ===============================================
g1 <- process_group("PBMC_Pac",
                    "F:/Blood-Patient-Paclitaxel-Analysis/Figs/TOP10-Clusters-final.Pacli.csv",
                    "F:/WGCNA/Blood-Patient-Paclitaxel-Analysis/Results/Max_Module_Membership_Per_Gene.csv",
                    "F:/Blood-Patient-Paclitaxel-Analysis/Figs/Enrich/GO_Enrichment.csv",
                    "F:/Blood-Patient-Paclitaxel-Analysis/Figs/Enrich/KEGG_Enrichment.csv"
)
all_biomarkers <- bind_rows(g1)
outdir <- "F:/TOP-BIOMARKERS/Blood_Pac"

write_csv(all_biomarkers, file.path(outdir, "All_Biomarkers_AllGroups.csv"))
# ===============================================
outdir <- "F:/TOP-BIOMARKERS/Blood_Pac"
write_csv(all_biomarkers, "F:/TOP-BIOMARKERS/Blood_Pac/All_Biomarkers_AllGroups.csv")
all_biomarkers <- read_csv(file.path(outdir, "All_Biomarkers_AllGroups.csv"))
all_biomarkers_weighted <- all_biomarkers %>%
  dplyr::filter(!is.na(Max_Membership)) %>% 
  mutate(weighted_score = 0.4 * evidence_score + 0.6 * Max_Membership)
top_gene_per_cluster <- all_biomarkers_weighted %>%
  group_by(gene) %>%
  slice_max(order_by = weighted_score, n = 1, with_ties = FALSE) %>%
  ungroup()
top80_weighted <- top_gene_per_cluster[order(-top_gene_per_cluster$weighted_score), ][1:80, ]
write_csv(top80_weighted, file.path(outdir, "Top80_Biomarkers_WeightedScore.csv"))
# =============================================== TOP 100 Heatmap Plot
outdir <- "F:/TOP-BIOMARKERS/Blood_Pac"
top80_file <- file.path(outdir, "Top80_Biomarkers_WeightedScore.csv")
df <- read_csv(top80_file)
cluster_scores <- df %>%
  group_by(cluster_id) %>%
  summarise(mean_score = mean(weighted_score, na.rm = TRUE)) %>%
  arrange(desc(mean_score))
p <- ggplot(cluster_scores, aes(x = reorder(cluster_id, mean_score), y = mean_score, fill = mean_score)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_text(aes(label = round(mean_score, 2)), hjust = -0.2, size = 3.5) +
  scale_fill_gradient(low = "#A6CEE3", high = "#08306B") +
  coord_flip() +
  labs(
    title = "Average Weighted Score per Cell Cluster",
    x = "Cell Type (Cluster)",
    y = "Mean Weighted Score",
    fill = "Score"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "none", 
    plot.margin = margin(20, 30, 20, 30)
  )
ggsave(file.path(outdir, "Cluster_Average_WeightedScore_ProfessionalPlot.png"),
       plot = p, width = 10, height = 7, dpi = 600, bg = "white")
outdir <- "F:/TOP-BIOMARKERS/Blood_Pac"
plot_path <- file.path(outdir, "Cluster_WeightedScore_Boxplot.png")
p <- ggplot(df, aes(x = reorder(cluster_id, weighted_score, FUN = median), y = weighted_score)) +
  geom_boxplot(fill = "#6BAED6", outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkblue", size = 1.5) +
  coord_flip() +
  labs(title = "Distribution of Weighted Scores per Cluster",
       x = "Cell Type (Cluster)", y = "Weighted Score") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )
ggsave(filename = plot_path,
       plot = p,
       width = 12,
       height = 8,
       dpi = 600,
       bg = "white")
# ======================== TOP 40 FOR CLEAN HEATMAP
outdir <- "F:/TOP-BIOMARKERS/Blood_Pac"
top_file <- file.path(outdir, "Top80_Biomarkers_WeightedScore.csv")
df <- readr::read_csv(top_file)
top40_df <- df[order(-as.numeric(df$weighted_score)), ][1:40, ]
mat <- top40_df %>%
  dplyr::select(gene, weighted_score) %>%
  column_to_rownames("gene")
annotation_row <- data.frame(CellType = top40_df$cluster_id)
rownames(annotation_row) <- top40_df$gene
png(file.path(outdir, "Top40_Biomarkers_Heatmap_With_CellTypes_Weighted.png"),
    width = 1200, height = 900)
pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "yellow", "darkblue"))(100),
         annotation_row = annotation_row,
         main = "Top 40 Biomarkers (Weighted Score) with Cell Types",
         fontsize_row = 10,
         fontsize_col = 10)
dev.off()
# ======================================
df <- readr::read_csv("F:/TOP-BIOMARKERS/Blood_Pac/Top80_Biomarkers_WeightedScore.csv")
df <- as.data.frame(df)
df_top20 <- df[order(-df$weighted_score), ][1:20, ]
df_top20$gene <- factor(df_top20$gene, levels = rev(df_top20$gene))
df_top20$cluster_id <- droplevels(factor(df_top20$cluster_id))
p <- ggplot(df_top20, aes(x = weighted_score, y = gene, fill = cluster_id)) +
  geom_col(width = 0.6) +
  labs(
    title = "Top 20 Biomarkers by Weighted Score",
    x = "Weighted Score",
    y = "Gene"
  ) +
  scale_fill_brewer(palette = "Set2", na.value = "gray") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
ggsave("F:/TOP-BIOMARKERS/Blood_Pac/Top20_Biomarkers_Barplot_Weighted_FIXED.png",
       plot = p, width = 9, height = 6, dpi = 600)
# ========== 
top_file <- "F:/TOP-BIOMARKERS/Blood_Pac/Top80_Biomarkers_WeightedScore.csv"
df <- readr::read_csv(top_file)
genes <- df %>%
  arrange(desc(weighted_score)) %>%
  slice_head(n = 40) %>%
  pull(gene)
# ========== 
id_map <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# ========== GO
go <- enrichGO(gene = id_map$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               keyType = "ENTREZID",
               pAdjustMethod = "BH",
               qvalueCutoff = 0.05)
# ========== 
go@result <- go@result %>%
  dplyr::filter(p.adjust < 0.01) %>%
  mutate(Description = stringr::str_wrap(Description, width = 45)) %>%
  head(10)
go@result$geneID <- sapply(go@result$geneID, function(x) {
  ids <- unlist(strsplit(x, "/"))
  symbols <- id_map$SYMBOL[match(ids, id_map$ENTREZID)]
  paste(symbols, collapse = "/")
})
# ========== KEGG
kegg <- enrichKEGG(gene = id_map$ENTREZID,
                   organism = "hsa",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
kegg@result <- kegg@result %>%
  dplyr::filter(p.adjust < 0.01) %>%
  dplyr::mutate(Description = stringr::str_wrap(Description, width = 45)) %>%
  head(10)
kegg@result$geneID <- sapply(kegg@result$geneID, function(x) {
  ids <- unlist(strsplit(x, "/"))
  symbols <- id_map$SYMBOL[match(ids, id_map$ENTREZID)]
  paste(symbols, collapse = "/")
})
# =====================
kegg <- enrichKEGG(gene = gene_list_entrez, organism = "hsa", pvalueCutoff = 0.05)
head(kegg@result)
write.csv(kegg@result, "F:/TOP-BIOMARKERS/Blood_Pac/KEGG_Results.csv", row.names = FALSE)
kegg_pathways_genes <- kegg@result %>%
  dplyr::filter(p.adjust < 0.1) %>%
  dplyr::select(Description, geneID) %>%
  dplyr::mutate(Genes = gsub("/", ", ", geneID)) %>%
  dplyr::select(-geneID)
write.csv(kegg_pathways_genes, "F:/TOP-BIOMARKERS/Blood_Pac/KEGG_Pathway_GeneList.csv", row.names = FALSE)
if (nrow(kegg@result) > 0) {
  ggsave("F:/TOP-BIOMARKERS/Blood_Pac/KEGG_barplot_top40_filtered.pdf",
         plot = barplot(kegg, showCategory = 10) + 
           ggtitle("Immune-related KEGG Pathways Enriched in Top 40 Biomarkers") +
           theme_classic(base_size = 14) +
           theme(
             plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
             axis.text = element_text(size = 14),
             axis.title = element_text(size = 14),
             plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
           ),
         width = 15, height = 10, bg = "white")
}
# ===================== 
library(patchwork)
library(enrichplot)
library(ggplot2)
p_go_body <- cnetplot(go,
                      showCategory = 10,
                      node_label = "all",
                      circular = TRUE,
                      color.params = list(edge = TRUE)) +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "") +
  theme(
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )
p_go <- wrap_elements(full = p_go_body) +
  plot_annotation(
    title = "Key GO Biological Processes Enriched in Top 40 Biomarkers",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20)
    )
  )
ggsave("F:/TOP-BIOMARKERS/Blood_Pac/Key_GO_Biological_Processes_Enriched_in_Top_40_Biomarkers.pdf",
       plot = p_go,  
       width = 15, height = 10, dpi = 600, bg = "white", device = cairo_pdf)



