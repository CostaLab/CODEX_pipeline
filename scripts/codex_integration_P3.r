library(Seurat)
library(tidyverse)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(Matrix)

## ---- load shared config ----
this_file <- sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
script_dir <- dirname(normalizePath(this_file))
source(file.path(script_dir, "r_config.R"))

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else file.path(script_dir, "..", "config.yaml")
cfg <- load_config(config_path)
pth <- get_pipeline_paths(cfg)

fig_dir <- pth$fig_dir
out_data_dir <- pth$out_data_dir

chosen_res <- cfg$integration$chosen_resolution
cluster_palette <- cfg$integration$cluster_palette

## ---- configure python environment ----
# 1. Extract the R-specific conda environment name from your yaml structure
#    Prefer an explicit path (some HPC conda installs live in locations
#    reticulate's name-based lookup can't find) but fall back to the name,
#    and if neither is set, skip the explicit bind entirely - reticulate
#    will then auto-detect whatever python is already on PATH, which is
#    normally fine since the env was already `conda activate`-d before this
#    script started (see slurm/*.sbatch / run_pipeline.sh).
conda_env_name <- cfg$slurm$r_conda_env
conda_env_path <- cfg$slurm$r_conda_env_path

# 2. Bind reticulate to this environment before any python commands run
if (!is.null(conda_env_path) && nzchar(conda_env_path)) {
  reticulate::use_condaenv(conda_env_path, required = TRUE)
} else if (!is.null(conda_env_name) && nzchar(conda_env_name)) {
  reticulate::use_condaenv(conda_env_name, required = TRUE)
}


library(scCustomize)

###########

cur_seurat <- readRDS(paste0(out_data_dir, 'seurat_object_harmony_clusts_filt_bg.rds'))

Idents(cur_seurat) <- chosen_res

n_clusters <- length(unique(Idents(cur_seurat)))
cluster_colors <- colorRampPalette(brewer.pal(min(n_clusters, 8), cluster_palette))(n_clusters)

pdf(paste0(fig_dir, '/umap_clusters_final.pdf'), width = 7, height = 7)
p <- DimPlot(cur_seurat, reduction = 'umap', cols = cluster_colors, raster = FALSE, label = TRUE) +
  ggtitle(paste('Clusters @', chosen_res))
print(p)
dev.off()

pdf(paste0(fig_dir, '/umap_clusters_by_condition.pdf'), width = 9, height = 7)
p2 <- DimPlot(cur_seurat, reduction = 'umap', group.by = 'condition', raster = FALSE)
print(p2)
dev.off()

## marker expression per cluster
markers <- rownames(cur_seurat)
pdf(paste0(fig_dir, '/dotplot_markers_by_cluster.pdf'), width = 12, height = 6)
DotPlot(cur_seurat, features = markers) + RotatedAxis()
dev.off()

## cluster markers
clustMarkers <- FindAllMarkers(cur_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(clustMarkers, paste0(out_data_dir, 'clusterMarkers.csv'))

pdf(paste0(fig_dir, '/dotplot_cluster_markers.pdf'), height=4, width=9)
DotPlot(cur_seurat, features = unique(as.character(unlist(tapply(clustMarkers$gene, clustMarkers$cluster, function(x){head(x, 2)}))))) + RotatedAxis()
dev.off()
##

pdf(paste0(fig_dir, '/featureplots_markers.pdf'), width = 16, height = 24)
FeaturePlot(cur_seurat, features = markers, raster = FALSE, ncol = 6)
dev.off()

## cluster composition per condition/sample
comp_table <- table(cur_seurat$condition, Idents(cur_seurat))
write.csv(comp_table, paste0(out_data_dir, 'cluster_composition_by_condition.csv'))

comp_table_sample <- table(cur_seurat$sample, Idents(cur_seurat))
write.csv(comp_table_sample, paste0(out_data_dir, 'cluster_composition_by_sample.csv'))

saveRDS(cur_seurat, file = paste0(out_data_dir, 'seurat_object_final.rds'))
write.csv(cur_seurat@meta.data, paste0(out_data_dir,'seurat_object_final_metadata.csv'))
counts_matrix <- GetAssayData(object = cur_seurat, assay='RNA')
writeMM(counts_matrix, paste0(out_data_dir,'seurat_object_final_counts.mtx'))
write.csv(Features(cur_seurat), paste0(out_data_dir,'seurat_object_final_features.csv'))
as.anndata(x = cur_seurat, file_path = out_data_dir, file_name = "seurat_object_final.h5ad")

print('END P3!!!')
