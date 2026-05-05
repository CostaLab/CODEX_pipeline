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

## Folders
fig_dir <- "/beegfs/data/SchneiderLab/CODEX_processing/pipeline_run/visual_output/"
out_data_dir <- "/beegfs/data/SchneiderLab/CODEX_processing/pipeline_run/"

## variable to regress for batch correction
vars_to_regress=c('sample', 'run')

## number of PCs to use
npcs=25

## clustering resolution range to try
res_init=0.1
res_final=0.5
chosen_res=0.2
###########

cur_seurat<-readRDS(paste0(out_data_dir, 'seurat_object_sct_filt_bg.rds'))

# compute PCA:
cur_seurat <- RunPCA(cur_seurat)

png(paste0(fig_dir, "pca_elbow_plot_filt_bg.png"), width=5, height=3, res=300, units='in')
ElbowPlot(cur_seurat)
dev.off()

## % variance explained
pct=cur_seurat@reductions$pca@stdev/sum(cur_seurat@reductions$pca@stdev)
print(paste('variance explained by',npcs,'PCs:', sum(pct[1:npcs])))

# UMAP and clustering with top PCs
cur_seurat <- RunUMAP(cur_seurat, reduction='pca', dims = 1:npcs)
cur_seurat <- FindNeighbors(cur_seurat, reduction='pca')
cur_seurat <- FindClusters(cur_seurat, resolution = chosen_res)

pdf(paste0(fig_dir, "umap_clusters_sct_filt_bg.pdf"), width=7, height=7)
p1 <- DimPlot(cur_seurat, group.by='seurat_clusters', reduction='umap', raster = FALSE) +
  ggtitle('seurat clusters') 
print(p1)
dev.off()

pdf(paste0(fig_dir, "umap_sample_sct_filt_bg.pdf"), width=7, height=7)
p3 <- DimPlot(cur_seurat, group.by='sample', reduction='umap', raster = FALSE) +
  ggtitle('Sample') 
print(p3)
dev.off()

pdf(paste0(fig_dir, "umap_run_sct_filt_bg.pdf"), width=7, height=7)
p4 <- DimPlot(cur_seurat, group.by='run', reduction='umap', raster = FALSE) +
  ggtitle('Run') 
print(p4)
dev.off()

## Harmony
cur_seurat <- cur_seurat %>%
  RunHarmony(vars_to_regress, assay.use="SCT")

saveRDS(cur_seurat, file=paste0(out_data_dir, 'seurat_object_harmony_filt_bg.rds'))
print('Checkpoint 2!')

#cur_seurat<-readRDS(paste0(out_data_dir, 'seurat_object_harmony_filt_bg_rev.rds'))

# UMAP and clustering
cur_seurat <- RunUMAP(cur_seurat, reduction='harmony', dims = 1:npcs)
cur_seurat <- FindNeighbors(cur_seurat, reduction='harmony')
#cur_seurat <- FindClusters(cur_seurat, resolution = 0.2)

## Multiple resolutions
for(res in seq(res_init,res_final,by=0.1)){
    cur_seurat <- FindClusters(cur_seurat, resolution = res)
}

pdf(paste0(fig_dir, 'clustree_harmony_filt_bg.pdf'), height=9, width=7)
clustree(cur_seurat, prefix = "SCT_snn_res.")
dev.off()

saveRDS(cur_seurat, file=paste0(out_data_dir, 'seurat_object_harmony_clusts_filt_bg.rds'))
print('Checkpoint 3!')
###### Checkpoint 3

print('END!!!')
