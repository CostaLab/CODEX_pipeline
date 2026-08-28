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

## ---- load shared config ----
this_file <- sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
script_dir <- dirname(normalizePath(this_file))
source(file.path(script_dir, "r_config.R"))

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else file.path(script_dir, "..", "config.yaml")
cfg <- load_config(config_path)
pth <- get_pipeline_paths(cfg)
rl <- get_run_lists(cfg)

project_name <- cfg$project$name
fig_dir <- pth$fig_dir
out_data_dir <- pth$out_data_dir
quant_dir <- pth$quant_dir

intensities_file_prefix <- cfg$integration$intensities_file_prefix
intensities_file_suffix <- cfg$integration$intensities_file_suffix
not_working <- unlist(cfg$integration$not_working_markers)
combs <- names(cfg$integration$aberrant_combinations)
ab_prot_combs <- lapply(cfg$integration$aberrant_combinations, unlist)
qthr <- cfg$integration$aberrant_quantile_threshold
min_cells_thr <- cfg$integration$min_cells_threshold
min_features_thr <- cfg$integration$min_features_threshold
vars_to_regress <- unlist(cfg$integration$vars_to_regress)

runs <- paste0("run", sapply(cfg$runs, function(r) r$id))
condition_vecs <- rl$condition_vecs
working_samples <- rl$working_samples

###############
cur_data <- data.frame()
for (run in runs) {
  cur_t <- data.frame()
  file_list <- list.files(path = quant_dir, pattern = paste0(intensities_file_prefix, run), full.names = TRUE)

  for (i in working_samples[[run]]) {
    ints <- read.csv(file_list[i], row.names = 1)
    ints[, 'condition'] <- condition_vecs[[run]][i]
    rownames(ints) <- paste0(rownames(ints), '_', condition_vecs[[run]][i], '_', i)
    ints[, 'sample'] <- paste0(run, '_', condition_vecs[[run]][i], '_', i)
    ints[, 'run'] <- run
    rownames(ints) <- paste0(0:(dim(ints)[1] - 1), '_', ints$sample)
    cur_t <- rbind(cur_t, ints)
  }
  if ((dim(cur_data)[1] > 0)) {
    cur_t[setdiff(names(cur_data), names(cur_t))] <- 0
    cur_data[setdiff(names(cur_t), names(cur_data))] <- 0
  }
  cur_data <- rbind(cur_data, cur_t)
}

rm(cur_t)
gc()

## cells per sample
pdf(paste0(fig_dir, '/ncells_sample.pdf'), 8, 5)
barplot(table(cur_data$sample))
dev.off()

table(cur_data$sample)

## Aberrant protein combinations filtering
for (comb in combs) {
  pdf(paste0(fig_dir, '/', comb, '_cells_filtering.pdf'), 8, 5)
  plot(cur_data[, ab_prot_combs[[comb]][1]], cur_data[, ab_prot_combs[[comb]][2]])
  ab_prot_combs_thr1 <- quantile(cur_data[, ab_prot_combs[[comb]][1]], na.rm = TRUE, probs = qthr)
  ab_prot_combs_thr2 <- quantile(cur_data[, ab_prot_combs[[comb]][2]], na.rm = TRUE, probs = qthr)
  abline(v = ab_prot_combs_thr1, col = 'red')
  abline(h = ab_prot_combs_thr2, col = 'red')
  dev.off()

  ## Filter
  cur_data <- cur_data[-which((cur_data[, ab_prot_combs[[comb]][1]] > ab_prot_combs_thr1)
                               & (cur_data[, ab_prot_combs[[comb]][2]] > ab_prot_combs_thr2)), ]
}

## cells per sample after filtering
pdf(paste0(fig_dir, '/ncells_sample_filt.pdf'), 8, 5)
barplot(table(cur_data$sample))
dev.off()

table(cur_data$sample)

## Protein filtering
cur_data <- cur_data[, -which(colnames(cur_data) %in% not_working)]

write.csv(cur_data, paste0(out_data_dir, 'cur_data.csv'))
print('Checkpoint 0!')

### BG filtering
for (smpl in unique(cur_data$sample)) {
  thr <- apply(cur_data[cur_data$sample == smpl, setdiff(colnames(cur_data), c('condition', 'run', 'sample'))], 2, median)
  for (x in setdiff(colnames(cur_data), c('condition', 'run', 'sample'))) {
    cur_data[cur_data$sample == smpl, x] <- cur_data[cur_data$sample == smpl, x] - as.numeric(thr[x])
  }
}
cur_data[cur_data < 0] <- 0

### Seurat object
cur_seurat <- CreateSeuratObject(
  counts = t(cur_data[, setdiff(colnames(cur_data), c('condition', 'run', 'sample'))]),
  min.cells = min_cells_thr,
  min.features = min_features_thr,
  project = project_name,
  meta.data = cur_data[, c('condition', 'sample', 'run')]
)

pdf(paste0(fig_dir, '/qc.pdf'), 8, 5)
VlnPlot(cur_seurat, group.by = "sample", features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0)
dev.off()

saveRDS(cur_seurat, file = paste0(out_data_dir, 'seurat_object_filt_bg_raw.rds'))

## SCTransform
cur_seurat <- SCTransform(
  cur_seurat,
  vars.to.regress = vars_to_regress,
  verbose = TRUE,
  vst.flavor = 'v1'
)
gc()

saveRDS(cur_seurat, file = paste0(out_data_dir, 'seurat_object_sct_filt_bg.rds'))
print('Checkpoint 1!')
