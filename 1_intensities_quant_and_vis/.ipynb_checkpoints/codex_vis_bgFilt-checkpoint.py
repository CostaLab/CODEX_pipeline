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

#conditions_list=['MPL', 'Ctrl', 'ThPO', 'ThPO', 'Ctrl']
condition_vecs={'21':['d180', 'd0', 'd0', 'd60', 'd180'],
                   '22':['d365', 'd0', 'd60'],
                   '23':['d0', 'd365', 'd60', 'd100'],
                   '24':['d0', 'd270', 'd0', 'd365', 'd180']}

for run in runs:

    dat = { f.split('/')[-1].split('.')[0]: xr.open_zarr(f) for f in glob(ZARR_OUTPUT_PATH+"run"+run+"*.zarr") }
    
    #channels=list(dat[list(dat.keys())[0]].channels.to_dataframe().index)
    #channels=pd.read_csv(marker_list_patt+run+'.txt', header=None)
    channels=pd.read_csv(marker_list_file, header=None)
    channels=list(channels[0])
    
    ## sorted as dict keys
    #condition_vec=conditions_list
    condition_vec=condition_vecs[run]
    
    j=0
    for i in list(dat.keys()):
    #for i in list(range(2, 5)):
        print(i)
    
        anndata=ad.AnnData(pd.DataFrame(dat[i]._intensity, columns=channels).iloc[:,1:len(channels)])
        #sp=pd.DataFrame(dat[i]._obs).loc[:,1:2]
        sp=pd.DataFrame(dat[i]._obs).iloc[:,-2:] ## last 2 columns
        #sp.columns=['x', 'y']
        sp.columns=['0', '1']
        sp.index=anndata.obs.index
        anndata.obsm['spatial']=sp.to_numpy()
        
        print(anndata.obs.shape)
        anndata.obs_names=[n+'_'+condition_vec[j]+'_'+str(j) for n in anndata.obs_names]
    
        ### Median filtering
        thr=anndata.to_df().median()
        dat1=anndata.to_df().copy()
        
        for c in dat1.columns:
            dat1[c]=dat1[c]-thr[c]
            dat1[c][dat1[c]<0]=0
    
        anndata.X=dat1.to_numpy()
        ### the following figures will be saved in a "figures" folder
        sc.pl.spatial(anndata, color=anndata.var_names, spot_size = 20, use_raw=False, save="run"+run+"_"+condition_vec[j]+'_'+str(j)+'_mks_bgFilt.png')
        j=j+1