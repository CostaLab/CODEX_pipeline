library(Seurat)
library(sctransform)
library(harmony)
library(tidyverse)
library(hdf5r)
library(patchwork)

library(cowplot)
library(RColorBrewer)
library(viridis)
library(Matrix)

options(future.globals.maxSize = 5000 * 1024^2)

## Folders
fig_dir <- "/beegfs/data/SchneiderLab/CODEX_processing/pipeline_run/visual_output/"
out_data_dir <- "/beegfs/data/SchneiderLab/CODEX_processing/pipeline_run/"

## resolution to use
chosen_res='SCT_snn_res.0.3'

## cluster color palette
cluster_palette="Paired"

#######################

cur_seurat<-readRDS(paste0(out_data_dir, 'seurat_object_harmony_clusts_filt_bg.rds'))
cur_seurat@meta.data$seurat_clusters<-cur_seurat@meta.data[,chosen_res]	
Idents(object = cur_seurat) <- "seurat_clusters"
# Classic palette BuPu, with 4 colors
newpal <- brewer.pal(12, cluster_palette) 
# Add more colors to this palette :
newpal <- colorRampPalette(newpal)(20)

#cluster markers with Wilcoxon test
clustMarkers <- FindAllMarkers(cur_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(clustMarkers, paste0(out_data_dir, 'clusterMarkers_filt_bg_r0.3.csv'))

pdf(paste0(fig_dir, 'posClustMarkers_filt_bg_r0.3.pdf'), height=4, width=9)
#png(paste0(fig_dir, 'pngs/posClustMarkers.png'), height=4, width=9, res=250, units='in')
DotPlot(cur_seurat, features = unique(as.character(unlist(tapply(clustMarkers$gene, clustMarkers$cluster, function(x){head(x, 2)}))))) + RotatedAxis()
dev.off()

pdf(paste0(fig_dir, 'allMarkers_bg_r0.3.pdf'), height=4, width=9)
#png(paste0(fig_dir, 'pngs/posClustMarkers.png'), height=4, width=9, res=250, units='in')
DotPlot(cur_seurat, features = rownames(cur_seurat)) + RotatedAxis()
dev.off()

print('END!!!')