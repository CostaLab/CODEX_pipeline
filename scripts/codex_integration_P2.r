library(Seurat)
library(sctransform)
library(harmony)
library(tidyverse)
library(hdf5r)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(clustree)

options(future.globals.maxSize = 5000 * 1024^2)

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

vars_to_regress <- unlist(cfg$integration$vars_to_regress)
npcs <- cfg$integration$npcs
res_init <- cfg$integration$clustering$res_init
res_final <- cfg$integration$clustering$res_final
res_step <- cfg$integration$clustering$res_step
###########

cur_seurat <- readRDS(paste0(out_data_dir, 'seurat_object_sct_filt_bg.rds'))

# compute PCA:
cur_seurat <- RunPCA(cur_seurat)

png(paste0(fig_dir, "/pca_elbow_plot_filt_bg.png"), width = 5, height = 3, res = 300, units = 'in')
ElbowPlot(cur_seurat)
dev.off()

## % variance explained
pct <- cur_seurat@reductions$pca@stdev / sum(cur_seurat@reductions$pca@stdev)
print(paste('variance explained by', npcs, 'PCs:', sum(pct[1:npcs])))

# UMAP and clustering with top PCs (pre-batch-correction, for comparison)
cur_seurat <- RunUMAP(cur_seurat, reduction = 'pca', dims = 1:npcs)
cur_seurat <- FindNeighbors(cur_seurat, reduction = 'pca')

pdf(paste0(fig_dir, "/umap_sample_sct_filt_bg.pdf"), width = 7, height = 7)
p3 <- DimPlot(cur_seurat, group.by = 'sample', reduction = 'umap', raster = FALSE) +
  ggtitle('Sample')
print(p3)
dev.off()

pdf(paste0(fig_dir, "/umap_run_sct_filt_bg.pdf"), width = 7, height = 7)
p4 <- DimPlot(cur_seurat, group.by = 'run', reduction = 'umap', raster = FALSE) +
  ggtitle('Run')
print(p4)
dev.off()

## Harmony
cur_seurat <- cur_seurat %>%
  RunHarmony(vars_to_regress, assay.use = "SCT")

saveRDS(cur_seurat, file = paste0(out_data_dir, 'seurat_object_harmony_filt_bg.rds'))
print('Checkpoint 2!')

# UMAP and clustering on Harmony-corrected space
cur_seurat <- RunUMAP(cur_seurat, reduction = 'harmony', dims = 1:npcs)
cur_seurat <- FindNeighbors(cur_seurat, reduction = 'harmony')

## Multiple resolutions, for the clustree sweep
for (res in seq(res_init, res_final, by = res_step)) {
  cur_seurat <- FindClusters(cur_seurat, resolution = res)
}

pdf(paste0(fig_dir, '/clustree_harmony_filt_bg.pdf'), height = 9, width = 7)
clustree(cur_seurat, prefix = "SCT_snn_res.")
dev.off()

saveRDS(cur_seurat, file = paste0(out_data_dir, 'seurat_object_harmony_clusts_filt_bg.rds'))
print('Checkpoint 3!')
###### Checkpoint 3

print('END!!!')
