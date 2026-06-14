# ===============================
# 0. CARTELLE DI LAVORO
# ===============================

dirData <- "data"
dirResults <- "results"

if (!dir.exists(dirData)) {
  dir.create(dirData)
}

if (!dir.exists(dirResults)) {
  dir.create(dirResults)
}

# ===============================
# 1. DOWNLOAD E IMPORTAZIONE DEI DATI DA GEO
# ===============================

library(GEOquery)

series <- "GSE32863"

getGEOSuppFiles(
  GEO = series,
  baseDir = dirData
)

file_expr <- list.files(
  path = file.path(dirData, series),
  pattern = "non-normalized.*\\.txt\\.gz$",
  full.names = TRUE,
  recursive = TRUE
)

expr <- read.delim(
  file_expr[1],
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

dim(expr)
head(expr[, 1:10])
colnames(expr)[1:10]

# ===============================
# 2. PRE-ELABORAZIONE
# ===============================

# Rimozione Detection Pval
gene_ids <- expr$ID_REF

expr_values <- expr[, !grepl("Detection", colnames(expr))]

rownames(expr_values) <- expr_values$ID_REF
expr_values$ID_REF <- NULL

dim(expr_values)
head(expr_values[, 1:6])

# Gruppi Tumor / Normal
sample_names <- colnames(expr_values)

group <- ifelse(
  grepl("_T$", sample_names),
  "Tumor",
  "Normal"
)

group <- factor(group, levels = c("Normal", "Tumor"))

table(group)

# Ricostruzione coppie matched
patient_id <- gsub("_(N|T)$", "", sample_names)

pair_table <- table(patient_id, group)

complete_pairs <- rownames(pair_table)[
  pair_table[, "Normal"] == 1 &
    pair_table[, "Tumor"] == 1
]

length(complete_pairs)

keep_samples <- patient_id %in% complete_pairs

expr_matched <- expr_values[, keep_samples]

group_matched <- group[keep_samples]
patient_matched <- patient_id[keep_samples]

table(group_matched)
length(unique(patient_matched))

# Matrice numerica
expr_matrix <- as.matrix(expr_matched)

mode(expr_matrix) <- "numeric"

# Trasformazione log2
expr_log2 <- log2(expr_matrix + 1)

# Calcolo IQR
iqr_values <- apply(expr_log2, 1, IQR)

# Grafico IQR
hist(
  iqr_values,
  breaks = 50,
  main = "Distribuzione IQR",
  xlab = "IQR",
  col = "gray"
)

# Soglia 10° percentile
iqr_threshold <- quantile(iqr_values, 0.10)

iqr_threshold

# Grafico con soglia
hist(
  iqr_values,
  breaks = 50,
  main = "Distribuzione IQR",
  xlab = "IQR",
  col = "gray"
)

abline(
  v = iqr_threshold,
  col = "red",
  lwd = 2
)

# Filtraggio
expr_filtered <- expr_log2[iqr_values > iqr_threshold, ]

dim(expr_log2)
dim(expr_filtered)

# Export dati filtrati
write.table(
  expr_filtered,
  file = file.path(dirResults, "GSE32863_expr_log2_IQRfiltered.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

metadata <- data.frame(
  Sample = colnames(expr_filtered),
  Patient = patient_matched,
  Group = group_matched
)

write.table(
  metadata,
  file = file.path(dirResults, "GSE32863_metadata_matched.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ===============================
# 3. FILTRAGGIO - DIFFERENTIAL EXPRESSION ANALYSIS
# ===============================

# Ordinamento delle coppie matched
metadata_normal <- metadata[metadata$Group == "Normal", ]
metadata_tumor  <- metadata[metadata$Group == "Tumor", ]

metadata_normal <- metadata_normal[order(metadata_normal$Patient), ]
metadata_tumor  <- metadata_tumor[order(metadata_tumor$Patient), ]

all(metadata_normal$Patient == metadata_tumor$Patient)

dataN <- expr_filtered[, metadata_normal$Sample]
dataT <- expr_filtered[, metadata_tumor$Sample]

dim(dataN)
dim(dataT)

# Calcolo log Fold Change
logFC <- rowMeans(dataT) - rowMeans(dataN)

summary(logFC)

# Calcolo p-value con t-test paired
pval <- numeric(nrow(expr_filtered))

for (i in 1:nrow(expr_filtered)) {
  
  pval[i] <- tryCatch(
    t.test(
      as.numeric(dataT[i, ]),
      as.numeric(dataN[i, ]),
      paired = TRUE
    )$p.value,
    error = function(e) NA
  )
}

summary(pval)

# Correzione FDR
pval_adj <- p.adjust(pval, method = "fdr")

summary(pval_adj)

# Tabella risultati
results <- data.frame(
  Gene = rownames(expr_filtered),
  pvalue = pval,
  pval_adj = pval_adj,
  logFC = logFC
)

results$direction <- ifelse(
  results$logFC > 0,
  "UP",
  "DOWN"
)

results <- results[order(results$pval_adj), ]

head(results)
dim(results)

# Soglie
fc_threshold <- 2
logFC_threshold <- log2(fc_threshold)

alpha <- 0.05

# Filtraggio DEG
deg <- results[
  abs(results$logFC) >= logFC_threshold &
    results$pval_adj <= alpha,
]

dim(deg)

table(deg$direction)

# Volcano plot
volcano_color <- rep("gray", nrow(results))

volcano_color[
  results$logFC >= logFC_threshold &
    results$pval_adj <= alpha
] <- "red"

volcano_color[
  results$logFC <= -logFC_threshold &
    results$pval_adj <= alpha
] <- "blue"

negLogFDR <- -log10(results$pval_adj)

plot(
  results$logFC,
  negLogFDR,
  pch = 20,
  col = volcano_color,
  main = "Volcano Plot",
  xlab = "log2 Fold Change",
  ylab = "-log10 adjusted p-value"
)

abline(
  v = c(-logFC_threshold, logFC_threshold),
  col = "blue",
  lty = 2,
  lwd = 2
)

abline(
  h = -log10(alpha),
  col = "red",
  lty = 2,
  lwd = 2
)

legend(
  "topright",
  legend = c("UP", "DOWN", "Not significant"),
  col = c("red", "blue", "gray"),
  pch = 20
)

# Export risultati
write.table(
  results,
  file = file.path(dirResults, "GSE32863_all_genes_results.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  deg,
  file = file.path(dirResults, "GSE32863_DEG_filtered.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ===============================
# 4. GRAFICI E ULTERIORI ANALISI
# ===============================

dim(deg)
head(deg)
table(deg$direction)

# 4.1 SELEZIONE DEI GENI TOP
top_up <- deg$Gene[which.max(deg$logFC)]
top_down <- deg$Gene[which.min(deg$logFC)]

top_up
top_down


# 4.2 BOXPLOT GENE PIU' UP-REGULATED
up_values <- data.frame(
  Expression = as.numeric(expr_filtered[top_up, ]),
  Group = metadata$Group,
  Patient = metadata$Patient
)

boxplot(
  Expression ~ Group,
  data = up_values,
  col = c("lightblue", "salmon"),
  main = paste("Top Up-regulated Gene:", top_up),
  xlab = "Group",
  ylab = "Expression"
)

stripchart(
  Expression ~ Group,
  data = up_values,
  vertical = TRUE,
  method = "jitter",
  pch = 20,
  col = "black",
  add = TRUE
)


# 4.3 BOXPLOT GENE PIU' DOWN-REGULATED
down_values <- data.frame(
  Expression = as.numeric(expr_filtered[top_down, ]),
  Group = metadata$Group,
  Patient = metadata$Patient
)

boxplot(
  Expression ~ Group,
  data = down_values,
  col = c("lightblue", "salmon"),
  main = paste("Top Down-regulated Gene:", top_down),
  xlab = "Group",
  ylab = "Expression"
)

stripchart(
  Expression ~ Group,
  data = down_values,
  vertical = TRUE,
  method = "jitter",
  pch = 20,
  col = "black",
  add = TRUE
)


# 4.4 HEATMAP
library(pheatmap)

top_n <- 50

top_deg_genes <- deg$Gene[1:top_n]

heatmap_matrix <- expr_filtered[top_deg_genes, ]

annotation_col <- data.frame(
  Group = metadata$Group
)

rownames(annotation_col) <- metadata$Sample

pheatmap(
  heatmap_matrix,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  main = paste("Heatmap of Top", top_n, "Differentially Expressed Genes")
)


# 4.5 PCA
pca_matrix <- t(expr_filtered[deg$Gene, ])

pca_res <- prcomp(pca_matrix, scale. = TRUE)

summary(pca_res)

pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Group = metadata$Group,
  Sample = metadata$Sample
)

plot(
  pca_df$PC1,
  pca_df$PC2,
  col = ifelse(pca_df$Group == "Tumor", "red", "blue"),
  pch = 19,
  xlab = "PC1",
  ylab = "PC2",
  main = "PCA of Differentially Expressed Genes"
)

legend(
  "topright",
  legend = c("Normal", "Tumor"),
  col = c("blue", "red"),
  pch = 19
)


# ===============================
# 5. ESPORTAZIONE DEI RISULTATI
# ===============================

final_results <- data.frame(
  Gene = deg$Gene,
  adj.P.Val = deg$pval_adj,
  logFC = deg$logFC,
  Regulation = ifelse(
    deg$logFC > 0,
    "Upregulated",
    "Downregulated"
  )
)

final_results <- final_results[order(final_results$adj.P.Val), ]

head(final_results)
dim(final_results)
table(final_results$Regulation)

write.table(
  final_results,
  file = file.path(dirResults, "Differentially_Expressed_Genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

file.exists(
  file.path(dirResults, "Differentially_Expressed_Genes.txt")
)

# ===============================
# 6. ANALISI DI ARRICCHIMENTO FUNZIONALE
# ===============================

library(AnnotationDbi)
library(illuminaHumanv3.db)
library(enrichR)

# 6.1 CONVERSIONE PROBE ID -> GENE SYMBOL
probe_ids <- final_results$Gene

gene_symbols <- mapIds(
  illuminaHumanv3.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

gene_symbols <- na.omit(gene_symbols)
gene_symbols <- unique(gene_symbols)

length(gene_symbols)
head(gene_symbols)

# 6.2 DATABASE DISPONIBILI IN ENRICHR
dbs <- listEnrichrDbs()

head(dbs)

databases <- c(
  "GO_Biological_Process_2023",
  "KEGG_2021_Human",
  "Reactome_2022"
)

# 6.3 ENRICHMENT ANALYSIS
enrich_results <- enrichr(
  gene_symbols,
  databases
)

# 6.4 VISUALIZZAZIONE RISULTATI
GO_results <- enrich_results[["GO_Biological_Process_2023"]]
KEGG_results <- enrich_results[["KEGG_2021_Human"]]
Reactome_results <- enrich_results[["Reactome_2022"]]

head(GO_results)
head(KEGG_results)
head(Reactome_results)

# 6.5 EXPORT RISULTATI ENRICHR
write.table(
  GO_results,
  file = file.path(dirResults, "GO_results_GSE32863.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  KEGG_results,
  file = file.path(dirResults, "KEGG_results_GSE32863.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  Reactome_results,
  file = file.path(dirResults, "Reactome_results_GSE32863.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
