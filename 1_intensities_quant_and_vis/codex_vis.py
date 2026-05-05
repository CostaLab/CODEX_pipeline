import sys
import os
from glob import glob
import pandas as pd
import spatialproteomics
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import tifffile
import scanpy as sc
import anndata as ad

sc.set_figure_params(dpi=100, dpi_save=900)

### variables needed
VISUAL_OUTPUT_PATH="../visual_output/"
ZARR_OUTPUT_PATH="../resulting_data/zarr/test_"
QUANT_OUTPUT_PATH="../resulting_data/quantifications/"
## scanpy plots are stored in a 'figures' folder that will be created where the script is stored 
#marker_list_patt='/data/Schneider_lab/IRF8_project_AG/Segmentation_run62_65_66_Anna/MarkerList'
marker_list_file='/beegfs/data/SchneiderLab/CODEX_processing/data/MarkerList.txt'
runs=['21', '22', '23', '24']
intensities_file_suffix="_asinh.csv"

#conditions_list=['MPL', 'Ctrl', 'ThPO', 'ThPO', 'Ctrl']
condition_vecs={'21':['d180', 'd0', 'd0', 'd60', 'd180'],
                   '22':['d365', 'd0', 'd60'],
                   '23':['d0', 'd365', 'd60', 'd100'],
                   '24':['d0', 'd270', 'd0', 'd365', 'd180']}

scatter_plots_mkrs=['025 CD34 human', '010 CD14', '084 CD56 (NCam1)', '004 Gycophorin A', '072 Perilipin']

######

for run in runs:

    dat = { f.split('/')[-1].split('.')[0]: xr.open_zarr(f) for f in glob(ZARR_OUTPUT_PATH+"run"+run+"*.zarr") }
    
    #channels=list(dat[list(dat.keys())[0]].channels.to_dataframe().index)
    #channels=pd.read_csv(marker_list_patt+run+'.txt', header=None)
    channels=pd.read_csv(marker_list_file, header=None)
    channels=list(channels[0])
    
    ## store intensity values per cell
    for i in list(dat.keys()):
        pd.DataFrame(dat[i]._intensity.data, columns=channels).to_csv(QUANT_OUTPUT_PATH+i+intensities_file_suffix)
    
    #condition_vec=conditions_list
    condition_vec=condition_vecs[run]
    
    j=0
    for i in list(dat.keys()):
        print(i)
    
        ### Create anndata from sample in  zarr 
        anndata=ad.AnnData(pd.DataFrame(dat[i]._intensity, columns=channels).iloc[:,1:len(channels)])
        #sp=pd.DataFrame(dat[i]._obs).loc[:,1:2]
        sp=pd.DataFrame(dat[i]._obs).iloc[:,-2:] ## last 2 columns
        sp.columns=['0', '1']
        sp.index=anndata.obs.index
        anndata.obsm['spatial']=sp.to_numpy()
        
        print(anndata.obs.shape)
        
        ### add condition to cell ids
        anndata.obs_names=[n+'_'+condition_vec[j]+'_'+str(j) for n in anndata.obs_names]
        
    
        ### All markers in whole slice
        #sc.pl.spatial(anndata, color=anndata.var_names, spot_size = 20, use_raw=False, save="run"+run+"_"+condition_vec[j]+'_'+str(j)+'_mks.png')
    
        ### Histogram per protein
        df = pd.DataFrame(dat[i]._intensity.data, columns=channels)
        fig, axes = plt.subplots(8, 6, figsize=(16,24), sharex=True, sharey=True)
        ax=axes.flatten()
        
        for x in range(len(channels)):
            _ = df.iloc[:,x].hist(bins=100, ax=ax[x], legend=True)
        #plt.savefig('./visual_output_rev/allgenes_hist_'+i+'.png')
        plt.savefig(VISUAL_OUTPUT_PATH+'allmks_hist'+"run"+run+"_"+condition_vec[j]+'_'+str(j)+'.png')
    
        ### Scatter plots for protein that should not be with all others in the same cell (CD3 here)
        ## make grid for plots
        for mkr in scatter_plots_mkrs:
            gridCoords=[[list(range(8))[i], list(range(6))[j]] for i in range(8) for j in range(6)]
            fig, axes = plt.subplots(8, 6, figsize=(16, 24))
            for i in range(len(anndata.var_names)): 
                axes[gridCoords[i][0]][gridCoords[i][1]].scatter(anndata.X[:,anndata.var_names.values==mkr], anndata.X[:,i])
                axes[gridCoords[i][0]][gridCoords[i][1]].set_title(anndata.var_names[i])
            fig.savefig(VISUAL_OUTPUT_PATH+"run"+run+"_"+condition_vec[j]+'_'+str(j)+mkr+'_scatter.png')
        
        j=j+1