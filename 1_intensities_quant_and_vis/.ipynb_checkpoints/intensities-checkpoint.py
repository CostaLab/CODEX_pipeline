import sys
import os
from glob import glob
import pandas as pd
import spatialproteomics
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import tifffile
import zarr
import dask.array as da

def fix_for_zarr(obj):
    # 1. Fix Global Attributes
    obj.attrs = {str(k): (v.item() if hasattr(v, 'item') else v) 
                 for k, v in obj.attrs.items()}
    
    # 2. Fix Variable Attributes AND Encoding (The "Dask/Tiff" fix)
    for var in obj.variables:
        # Clear encoding - this is often where the hidden int32 keys live
        obj[var].encoding = {} 
        
        # Clean variable-specific attrs
        obj[var].attrs = {str(k): (v.item() if hasattr(v, 'item') else v) 
                          for k, v in obj[var].attrs.items()}
    return obj

### different quantification methods of marker intensity per cell
def intensity_sum(regionmask, intensity):
    return np.sum(intensity[regionmask])

def intensity_median(regionmask, intensity):
    return np.median(intensity[regionmask])

def intensity_sum_arcsinh_scaled(regionmask, intensity):
    x=np.sum(intensity[regionmask])
    return np.arcsinh(1+x/5)

### variables needed
INPUT_PATH="/beegfs/data/SchneiderLab/CODEX_processing/data/cropped_images/"
SEG_PATH="/beegfs/data/SchneiderLab/CODEX_processing/data/Cell_masks_cropped/"
## second segmentation mask
SEG2_PATH='/beegfs/data/SchneiderLab/CODEX_processing/data/MK_masks_cropped/'

VISUAL_OUTPUT_PATH="../visual_output/"
ZARR_OUTPUT_PATH="..resulting_data/zarr/"

## keys of data, seg and seg2 dicts must be the same. Adjust the next variables accordingly
seg_files_suffix='Cell_mask_'
seg2_files_suffix='MK_mask_'
img_files_suffix='Fusion_run'

seg_files_suffix_replace='run'
seg2_files_suffix_replace='run'
img_files_suffix_replace='run'
##

marker_list_file='/beegfs/data/SchneiderLab/CODEX_processing/data/MarkerList.txt'

intensity_func='intensity_mean'
## If you prefer to use one of the functions above (median, sum...) then the name of the function goes without quotes
#intensity_func=intensity_sum

segm_mks=["DAPI", "127 CD61"]


###

data={ f.split('/')[-1].replace('.tif', ''): f for f in glob(INPUT_PATH+'*.tif') }
for key in list(data.keys()):
    data[key.replace(img_files_suffix, img_files_suffix_replace)] = data.pop(key)


#seg = { f.split('/')[-1].replace(seg_files_suffix, ''): f for f in glob(SEG_PATH+'*tif') }
seg={ f.split('/')[-1].replace('.tif', ''): f for f in glob(SEG_PATH+'*.tif') }
for key in list(seg.keys()):
    seg[key.replace(seg_files_suffix, seg_files_suffix_replace)] = seg.pop(key)

## NEW
seg2={ f.split('/')[-1].replace('.tif', ''): f for f in glob(SEG2_PATH+'*.tif') }
for key in list(seg2.keys()):
    seg2[key.replace(seg2_files_suffix, seg2_files_suffix_replace)] = seg2.pop(key)

channels=pd.read_csv(marker_list_file, header=None)
channels=list(channels[0])

dat = {}

for reg in list(data.keys()):
    #img = tifffile.imread(data[reg])

    segm = tifffile.imread(seg[reg])    
    ## New
    segm2 = tifffile.imread(seg2[reg])
    # Open specifically in read-only mode ('r')
    with tifffile.TiffFile(data[reg], mode='r') as tif:
        # Use asarray with out='memmap' to avoid explicit permission requests
        mmap = tif.asarray(out='memmap') 
        d_array = da.from_array(mmap, chunks=(1, mmap.shape[1], mmap.shape[2]))       
    
    # 4. Final check: this should NOT be NaN anymore
    print(f"Dask array chunks: {d_array.chunks}")
                
    sdata = spatialproteomics.load_image_data(d_array, channel_coords=channels, segmentation=segm)
    ## For Two segmentation masks
    sdata['_seg1']=sdata._segmentation
    sdata=sdata.pp.drop_layers('_segmentation')
    ##second mask
    sdata=sdata.pp.add_segmentation(segm2)
    sdata['_seg2']=sdata._segmentation
    sdata=sdata.pp.drop_layers('_segmentation')
    
    #dat[reg] = sdata.copy()
    ### Plot DAPI
    #plt.figure(figsize=(12,12))
    #_ = dat[reg].pp["DAPI"].pl.show()
    #plt.savefig(VISUAL_OUTPUT_PATH+'DAPI_'+reg+'.png', dpi=2400, bbox_inches='tight')
    #plt.close()
    
    ### Plot DAPI & CD42a segmentation masks (segm_mks)
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    _ = (
        sdata.pp[segm_mks]
        .pl.colorize(["blue", "red"])
        .pl.show(render_segmentation=True, segmentation_kwargs={"layer_key": "_seg1"}, ax=ax[0])
    )
    
    _ = (
        sdata.pp[segm_mks]
        .pl.colorize(["blue", "red"])
        .pl.show(render_segmentation=True, segmentation_kwargs={"layer_key": "_seg2"}, ax=ax[1])
    )
    plot_path = os.path.join(VISUAL_OUTPUT_PATH, 'segmentation_plot_pro_ind_'+reg+'.png')
    plt.savefig(plot_path, dpi=2400, bbox_inches='tight')
    plt.close()
    ### merge segmentation masks
    sdata = sdata.pp.merge_segmentation(
        layer_key=["_seg2", "_seg1"], ## big cells mask goes first
        labels=["127 CD61", "DAPI"],
        threshold=1,
    ).pp.add_segmentation("_merged_segmentation")
    ##
    sdata = sdata.pp.add_quantification(func=intensity_func).pp.transform_expression_matrix(method="arcsinh")
    sdata = sdata.pp.add_observations('area')
    ##
    ### try zarr data correction
    sdata = fix_for_zarr(sdata)
    dat[reg] = sdata.copy()
    ### Plot segmentation (merged)
    plt.figure(figsize=(12, 12))
    colors = ['blue', 'red']
    sdata = sdata.la.set_label_colors(segm_mks, colors)
    # merged masks
    _ = sdata.pp[segm_mks].pl.colorize(colors=colors).pl.show(render_segmentation=True)
    plt.axis('off')  # Hide axes for a cleaner output
    plot_path = os.path.join(VISUAL_OUTPUT_PATH, 'segmentation_plot_pro_'+reg+'.png')
    plt.savefig(plot_path, dpi=2400, bbox_inches='tight')
    plt.close()  # Close the plot to free up memory
    
    ### Plot segmentation 
    #plt.figure(figsize=(12, 12))
    #_ = dat[reg].pp['DAPI'].pl.colorize('gold').pl.show(render_segmentation=True)
    #plt.axis('off')  # Hide axes for a cleaner output
    #plot_path = os.path.join(VISUAL_OUTPUT_PATH, 'segmentation_plot_pro_'+reg+'.png')
    #plt.savefig(plot_path, dpi=2400, bbox_inches='tight')
    #plt.close()  # Close the plot to free up memory
    
    del d_array, segm, sdata

for k, v in dat.items():
    v.to_zarr(ZARR_OUTPUT_PATH+'test_'+k+'.zarr')