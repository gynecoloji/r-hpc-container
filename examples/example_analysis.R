#!/usr/bin/env Rscript
# Example Analysis Script for Bioinformatics R 4.5.1 Container
# This script demonstrates basic usage of key packages included in the container

# =============================================================================
# SETUP
# =============================================================================

cat("Starting example analysis...\n")
cat("R version:", R.version.string, "\n\n")

# Set random seed for reproducibility
set.seed(123)

# =============================================================================
# EXAMPLE 1: Check Package Availability
# =============================================================================

cat("=" %R% 70, "\n")
cat("EXAMPLE 1: Checking installed packages\n")
cat("=" %R% 70, "\n\n")

# List of key packages to check
key_packages <- c(
  "Seurat", "DESeq2", "DiffBind", "ChIPseeker",
  "ggplot2", "dplyr", "GenomicRanges", "monocle3"
)

for (pkg in key_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    version <- packageVersion(pkg)
    cat(sprintf("✓ %s (v%s) - Available\n", pkg, version))
  } else {
    cat(sprintf("✗ %s - Not available\n", pkg))
  }
}

# =============================================================================
# EXAMPLE 2: Simple DESeq2 Analysis
# =============================================================================

cat("\n", "=" %R% 70, "\n")
cat("EXAMPLE 2: DESeq2 differential expression analysis\n")
cat("=" %R% 70, "\n\n")

library(DESeq2)

# Create simulated count data
n_genes <- 1000
n_samples <- 6

# Simulate counts
counts <- matrix(rnbinom(n = n_genes * n_samples, mu = 100, size = 1/0.5),
                 ncol = n_samples)
rownames(counts) <- paste0("Gene", 1:n_genes)
colnames(counts) <- paste0("Sample", 1:n_samples)

# Create sample metadata
colData <- data.frame(
  condition = factor(rep(c("Control", "Treatment"), each = 3)),
  row.names = colnames(counts)
)

cat("Count matrix dimensions:", dim(counts), "\n")
cat("Sample metadata:\n")
print(colData)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData,
  design = ~ condition
)

# Run DESeq2 analysis
cat("\nRunning DESeq2 analysis...\n")
dds <- DESeq(dds, quiet = TRUE)

# Get results
res <- results(dds)
cat("\nDifferential expression results:\n")
cat("Total genes:", nrow(res), "\n")
cat("Significant genes (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")

# =============================================================================
# EXAMPLE 3: GenomicRanges Operations
# =============================================================================

cat("\n", "=" %R% 70, "\n")
cat("EXAMPLE 3: GenomicRanges operations\n")
cat("=" %R% 70, "\n\n")

library(GenomicRanges)

# Create sample genomic ranges
gr1 <- GRanges(
  seqnames = c("chr1", "chr1", "chr2"),
  ranges = IRanges(start = c(100, 500, 200), width = 100),
  strand = c("+", "-", "+"),
  score = c(5, 10, 15)
)

gr2 <- GRanges(
  seqnames = c("chr1", "chr1", "chr2"),
  ranges = IRanges(start = c(150, 450, 250), width = 50),
  strand = c("+", "-", "+"),
  score = c(7, 12, 18)
)

cat("Genomic ranges 1:\n")
print(gr1)

cat("\nGenomic ranges 2:\n")
print(gr2)

# Find overlaps
overlaps <- findOverlaps(gr1, gr2)
cat("\nNumber of overlaps:", length(overlaps), "\n")

# =============================================================================
# EXAMPLE 4: Data Visualization with ggplot2
# =============================================================================

cat("\n", "=" %R% 70, "\n")
cat("EXAMPLE 4: Creating visualizations\n")
cat("=" %R% 70, "\n\n")

library(ggplot2)
library(dplyr)

# Create sample data
set.seed(42)
plot_data <- data.frame(
  gene_id = paste0("Gene", 1:100),
  log2FC = rnorm(100, mean = 0, sd = 2),
  pvalue = runif(100, 0, 1)
) %>%
  mutate(
    padj = p.adjust(pvalue, method = "BH"),
    significant = padj < 0.05
  )

cat("Creating volcano plot...\n")
cat("Sample data (first 5 rows):\n")
print(head(plot_data, 5))

# Note: In container usage, save plots to file instead of displaying
# This example just demonstrates the plotting code
p <- ggplot(plot_data, aes(x = log2FC, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10(p-value)"
  )

# In real usage, save to file:
# ggsave("volcano_plot.pdf", p, width = 8, height = 6)

cat("Plot created successfully (not displayed in container mode)\n")

# =============================================================================
# EXAMPLE 5: Single-cell Analysis with Seurat
# =============================================================================

cat("\n", "=" %R% 70, "\n")
cat("EXAMPLE 5: Seurat single-cell analysis\n")
cat("=" %R% 70, "\n\n")

library(Seurat)

# Create simulated single-cell data
set.seed(123)
n_cells <- 100
n_features <- 200

# Simulate UMI counts
counts <- matrix(
  rpois(n_cells * n_features, lambda = 5),
  nrow = n_features,
  ncol = n_cells
)
rownames(counts) <- paste0("Gene", 1:n_features)
colnames(counts) <- paste0("Cell", 1:n_cells)

# Create Seurat object
cat("Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "ExampleProject",
  min.cells = 3,
  min.features = 10
)

cat("Seurat object summary:\n")
print(seurat_obj)

# Basic QC metrics
cat("\nNumber of cells:", ncol(seurat_obj), "\n")
cat("Number of features:", nrow(seurat_obj), "\n")

# =============================================================================
# EXAMPLE 6: Working with mlr3
# =============================================================================

cat("\n", "=" %R% 70, "\n")
cat("EXAMPLE 6: Machine learning with mlr3\n")
cat("=" %R% 70, "\n\n")

library(mlr3)

# Create a simple classification task
data <- data.frame(
  feature1 = rnorm(100),
  feature2 = rnorm(100),
  target = factor(sample(c("A", "B"), 100, replace = TRUE))
)

cat("Creating mlr3 classification task...\n")
task <- TaskClassif$new(
  id = "example",
  backend = data,
  target = "target"
)

cat("Task summary:\n")
print(task)
cat("Number of observations:", task$nrow, "\n")
cat("Number of features:", task$ncol - 1, "\n")
cat("Target classes:", paste(task$class_names, collapse = ", "), "\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n", "=" %R% 70, "\n")
cat("ANALYSIS COMPLETE\n")
cat("=" %R% 70, "\n\n")

cat("All examples completed successfully!\n")
cat("This demonstrates that key bioinformatics packages are working.\n\n")

cat("Next steps:\n")
cat("1. Replace simulated data with your real data\n")
cat("2. Modify analysis parameters for your specific use case\n")
cat("3. Save outputs to files (plots, tables, RDS objects)\n")
cat("4. Use --bind flag to mount your data directories\n\n")

cat("Example command:\n")
cat("singularity exec --bind /data:/data bioinformatic_r_4_5_1_v1.sif \\\n")
cat("  Rscript /data/your_analysis.R\n\n")

sessionInfo()