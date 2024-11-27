from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import glob
import os
import openslide
import json
import cv2
import numpy as np
import pandas as pd
import anndata as ad
import openslide
from scipy.spatial import KDTree
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import glob
import os

#### upload data anf infrastructure #####
# from json file to anndata

def process_slide(slide_path, hovernet_output_path, wsi_ext=".svs"):
    # Extract slide basename
    wsi_basename = os.path.basename(slide_path)
    wsi_basename = wsi_basename[:-(len(wsi_ext))]
    
    # Load slide
    slide = openslide.OpenSlide(slide_path)
    
    # Load the thumbnail image
    thumb_path_wsi = os.path.join(hovernet_output_path, 'thumb', wsi_basename + '.png')
    
    thumb = cv2.cvtColor(cv2.imread(thumb_path_wsi), cv2.COLOR_BGR2RGB)
    
    # Load JSON data
    json_path_wsi = os.path.join(hovernet_output_path, 'json', wsi_basename + '.json')

    
    # Prepare lists for storing information
    centroid_list_wsi = []
    bbox_list_wsi = []
    type_list_wsi = []

    # Add results to individual lists
    with open(json_path_wsi) as json_file:
        data = json.load(json_file)
        nuc_info = data['nuc']
        for inst in nuc_info.values():
            centroid_list_wsi.append(inst['centroid'])
            bbox_list_wsi.append(inst['bbox'])
            type_list_wsi.append(inst['type'])
            
    # Convert nucleus types to their respective names
    nucleus_types = {
        0: "no label",
        1: "neoplastic",
        2: "inflammatory",
        3: "stromal",
        4: "necrotic",
        5: "benign epithelial"
    }
    type_list_wsi = [nucleus_types[t] for t in type_list_wsi]

    # Create a DataFrame for each list
    df_centroid = pd.DataFrame(centroid_list_wsi, columns=['centroid_x', 'centroid_y'])
    df_type = pd.DataFrame(type_list_wsi, columns=['type'])

    # Create a data matrix with zeros 
    X = np.zeros((len(centroid_list_wsi), 1))

    # Create the anndata object
    adata = ad.AnnData(X=X, obs=df_type)

    # Store the spatial coordinates in the obsm attribute
    adata.obsm["spatial"] = df_centroid.values

    # For contours and bounding boxes, we store them in adata.uns as a dictionary for easier retrieval.
    contour_bbox_dict = {
        "bbox": bbox_list_wsi}
    # Save the bounding box info
    adata.uns['spatial'] = {}
    adata.uns['spatial']["HoVer_Net"] = contour_bbox_dict

    # Store the thumbnail image in adata.uns
    # adata.uns['spatial']["HoVer_Net"]["images"] = {'hires': thumb}
    # nov 3
    adata.uns['spatial']["HoVer_Net"]["image_paths"] = {'hires': thumb_path_wsi}

    # Convert the type to categorical instead of numeric
    adata.obs['type'] = adata.obs['type'].astype('category')

    # Close the slide after processing
    slide.close()
    
    return wsi_basename, adata

def add_intensity_features(adata, slide_path):
    slide = openslide.OpenSlide(slide_path)

    mean_intensities = []
    variance_intensities = []
    max_intensities = []
    min_intensities = []

    for bbox in adata.uns['spatial']["HoVer_Net"]["bbox"]:
        # Unpack the top-left and bottom-right points
        (x1, y1), (x2, y2) = bbox

        # Calculate width and height
        w, h = x2 - x1, y2 - y1

        # Read the region from the slide
        region = slide.read_region((x1, y1), 0, (w, h))
        region = np.array(region)[:, :, 0]  # Convert to grayscale if needed

        # Calculate intensity features for the region
        mean_intensities.append(np.mean(region))
        variance_intensities.append(np.var(region))
        max_intensities.append(np.max(region))
        min_intensities.append(np.min(region))

    # Add these features to the AnnData object
    adata.obs['mean_intensity'] = mean_intensities
    adata.obs['variance_intensity'] = variance_intensities
    adata.obs['max_intensity'] = max_intensities
    adata.obs['min_intensity'] = min_intensities
    return adata


def add_spatial_distribution_features(adata):
    # Extract centroid coordinates
    centroids = adata.obsm['spatial']
    
    # Build a KDTree for efficient nearest neighbor search
    tree = KDTree(centroids)
    
    # Compute nearest neighbor distance for each nucleus
    distances, _ = tree.query(centroids, k=2)  # k=2 because the first result is the point itself
    nearest_neighbor_distance = distances[:, 1]  # Take the second closest point
    
    # Add the feature to the AnnData object
    adata.obs['nearest_neighbor_distance'] = nearest_neighbor_distance
    return adata

def slide_has_outputs(slide_path, hovernet_output_path, wsi_ext=".svs"):
    wsi_basename = os.path.basename(slide_path)
    wsi_basename = wsi_basename[:-(len(wsi_ext))]
    thumb_path_wsi = os.path.join(hovernet_output_path, 'thumb', wsi_basename + '.png')
    json_path_wsi = os.path.join(hovernet_output_path, 'json', wsi_basename + '.json')
    return os.path.exists(thumb_path_wsi) and os.path.exists(json_path_wsi)
