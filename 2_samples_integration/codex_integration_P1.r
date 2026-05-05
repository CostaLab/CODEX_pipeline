library(Seurat)
library(sctransform)
library(harmony)
library(tidyverse)
library(hdf5r)
library(patchwork)

library(cowplot)
library(RColorBrewer)
library(viridis)

options(future.globals.maxSize = 5000 * 1024^2)

## Project name
project_name='CODEX_fusion_run21-24'

## Folders
WDIR_PATH<-"/beegfs/data/SchneiderLab/CODEX_processing"
run_folders_prefix=""
run_folders_suffix=""
runs_folder='/pipeline_run'
in_data_dir_pattern <- "/quantifications/"
intensities_file_prefix="test_"
intensities_file_suffix="_asinh.csv"


fig_dir <- "/beegfs/data/SchneiderLab/CODEX_processing/pipeline_run/visual_output/"
out_data_dir <- "/beegfs/data/SchneiderLab/CODEX_processing/pipeline_run/"
## markers that won't be used
not_working=c('X050.CD138', 'X121.c.kit', 'DAPI', 'X128.Na..K.ATPase', 'X007.CD4', 'X150.VCAM.1', 'X023.CD11c', 'X008.CD3e', 'X062.CXCL13', 'X125.MCT', 'X071.Podoplanin', 'X124.TGFB1', 'X129.SPP1', 'X151.ISG15', 'X152.MGP', 'X153.PDE5a')
## runs to integrate
runs=c('run21', 'run22', 'run23', 'run24')
## one condition_vec per run, named as the runs
condition_vecs=list(run21=c('d180', 'd0', 'd0', 'd60', 'd180'),
                    run22=c('d365', 'd0', 'd60'),
                    run23=c('d0', 'd365', 'd60', 'd100'),
                    run24=c('d0', 'd270', 'd0', 'd365', 'd180'))
## working samples per run , named as the runs
working_samples=list(run21=c(1:5),
                     run22=c(1:3),
                     run23=c(1:4),
                     run24=c(1:5))

## Aberrant protein combinations
combs=c('c1')
ab_prot_combs=list(c1=c('X125.MCT', 'X125.MCT'))
#ab_prot_combs_thr=5
qthr=0.9

## filtering thresholds
min_cells_thr=3
min_features_thr=3

## variable to regress for batch correction
vars_to_regress=c('sample', 'run')

###############
cur_data<-data.frame()
for(run in runs){
    cur_t<-data.frame()
    ##
    file_list <- list.files(path = paste0(WDIR_PATH, runs_folder, in_data_dir_pattern), pattern = paste0(intensities_file_prefix, run), full.names = TRUE)
    ##
    for(i in working_samples[[run]]){
        #ints<-read.csv(paste0(WDIR_PATH, run_folders_prefix,run, run_folders_suffix, in_data_dir_pattern, intensities_file_prefix, i, intensities_file_suffix), row.names=1)
        ints<-read.csv(file_list[i], row.names=1)
        ints[,'condition']=condition_vecs[[run]][i]
        rownames(ints)=paste0(rownames(ints), '_', condition_vecs[[run]][i], '_', i)
        ints[,'sample']<-paste0(run,'_',condition_vecs[[run]][i], '_', i)
        ints[,'run']<-run
        #ints[,'sample']<-i
        rownames(ints)=paste0(0:(dim(ints)[1]-1), '_', ints$sample)
        cur_t<-rbind(cur_t, ints)
    }
    if((dim(cur_data)[1]>0)){
        cur_t[setdiff(names(cur_data), names(cur_t))] <- 0
        cur_data[setdiff(names(cur_t), names(cur_data))] <- 0
    }
    cur_data=rbind(cur_data, cur_t)
}


rm(cur_t)
gc()

## cells per sample
pdf(paste0(fig_dir,'ncells_sample.pdf'),8,5)
barplot(table(cur_data$sample))
dev.off()

table(cur_data$sample)

## Aberrant protein combinations filtering
for(comb in combs){
    pdf(paste0(fig_dir,comb,'_cells_filtering.pdf'),8,5)
    plot(cur_data[,ab_prot_combs[[comb]][1]], cur_data[,ab_prot_combs[[comb]][2]])
    ab_prot_combs_thr1=quantile(cur_data[,ab_prot_combs[[comb]][1]],na.rm = T,probs = qthr)
    ab_prot_combs_thr2=quantile(cur_data[,ab_prot_combs[[comb]][2]],na.rm = T,probs = qthr)
    abline(v=ab_prot_combs_thr1, col='red')
    abline(h=ab_prot_combs_thr2, col='red')
    dev.off()
    
    ## Filter
    cur_data<-cur_data[-which((cur_data[,ab_prot_combs[[comb]][1]]>ab_prot_combs_thr1)
                              &(cur_data[,ab_prot_combs[[comb]][2]]>ab_prot_combs_thr2)),]    
}

## cells per sample after filtering
pdf(paste0(fig_dir,'ncells_sample_filt.pdf'),8,5)
barplot(table(cur_data$sample))
dev.off()

table(cur_data$sample)

## Protein filtering
cur_data<-cur_data[,-which(colnames(cur_data) %in% not_working)]

write.csv(cur_data, paste0(out_data_dir,'cur_data.csv'))
print('Checkpoint 0!')

### BG filtering
for(smpl in unique(cur_data$sample)){
    thr=apply(cur_data[cur_data$sample==smpl,setdiff(colnames(cur_data), c('condition', 'run', 'sample'))], 2, median)
    for(x in setdiff(colnames(cur_data), c('condition', 'run', 'sample'))){
        cur_data[cur_data$sample==smpl,x]=cur_data[cur_data$sample==smpl,x]-as.numeric(thr[x])
        }
    }
cur_data[cur_data<0]=0

### Seurat object
cur_seurat <- CreateSeuratObject(
    counts = t(cur_data[,setdiff(colnames(cur_data), c('condition', 'run', 'sample'))]),
    min.cells=min_cells_thr,
    min.features=min_features_thr,
    project=project_name,
    meta.data=cur_data[,c('condition', 'sample', 'run')]
)

pdf(paste0(fig_dir,'qc.pdf'),8,5)
VlnPlot(cur_seurat, group.by="sample", features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size=0, )
dev.off()

####New
saveRDS(cur_seurat, file=paste0(out_data_dir, 'seurat_object_filt_bg_raw.rds'))
#####

## SCTransform
cur_seurat <- SCTransform(
  cur_seurat,
  vars.to.regress = vars_to_regress,
  verbose=TRUE,
  vst.flavor = 'v1'
)
gc()

saveRDS(cur_seurat, file=paste0(out_data_dir, 'seurat_object_sct_filt_bg.rds'))
print('Checkpoint 1!')

