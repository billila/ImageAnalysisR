import os
os.environ['NUMEXPR_MAX_THREADS'] = '32'
from argparse import ArgumentParser
import openslide
import numpy as np
from PIL import Image
from tiatoolbox.wsicore.wsireader import WSIReader
from tqdm import tqdm 
import tifffile
from skimage.color import rgb2gray
from skimage.filters import threshold_otsu
import cv2


# For tissue segmentation
def tissue_segmentation(image, sthresh=20, sthresh_up=255, mthresh=7, close=0, use_otsu=False):
    """
    Enhanced tissue segmentation to be more robust and to include steps like HSV conversion,
    median blurring, and optional morphological operations.
    """
    # Convert to HSV and use the saturation channel for thresholding
    img_hsv = cv2.cvtColor(image, cv2.COLOR_RGB2HSV)
    img_saturation = cv2.medianBlur(img_hsv[:, :, 1], mthresh)

    # Apply Otsu's thresholding or fixed threshold
    if use_otsu:
        _, thresholded = cv2.threshold(img_saturation, 0, sthresh_up, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    else:
        _, thresholded = cv2.threshold(img_saturation, sthresh, sthresh_up, cv2.THRESH_BINARY)

    # Morphological closing to fill small holes in the detected tissue
    if close > 0:
        kernel = np.ones((close, close), np.uint8)
        thresholded = cv2.morphologyEx(thresholded, cv2.MORPH_CLOSE, kernel)
    
    return thresholded

def save_tissue_mask(mask, output_path):
    kernel = np.ones((5,5),np.uint8)
    dilated_mask = cv2.dilate(mask,kernel,iterations = 1)
    cv2.imwrite(output_path, dilated_mask)

    
########################################################
# for tiling

def get_tiles(img, tile_size):
    h, w, _ = img.shape
    pad_h = (tile_size - h % tile_size) % tile_size
    pad_w = (tile_size - w % tile_size) % tile_size
    padded_img = np.pad(img, [(0, pad_h), (0, pad_w), (0, 0)], mode='constant', constant_values=255)
    result = []
    padded_h, padded_w, _ = padded_img.shape
    for i in range(0, padded_h, tile_size):
        for j in range(0, padded_w, tile_size):
            tile = padded_img[i:i+tile_size, j:j+tile_size, :]
            result.append({'img': tile, 'idx': (i//tile_size, j//tile_size)})
    return result

def get_tissue_tiles(img, tile_size, tissue_mask):
    """Generate tiles from tissue areas based on a given mask."""
    h, w = tissue_mask.shape
    result = []
    for i in range(0, h, tile_size):
        for j in range(0, w, tile_size):
            if tissue_mask[i:i+tile_size, j:j+tile_size].any():  # Check if the tile contains tissue
                tile = img[i:i+tile_size, j:j+tile_size, :]
                result.append({'img': tile, 'idx': (i//tile_size, j//tile_size)})
    return result

def save_tiles(tiles, tiles_dir, prefix):
    # os.makedirs(tiles_dir, exist_ok=True) is no longer needed here as directories are created in process_svs_slide and process_tiff_slide
    for tile in tiles:
        img, idx = tile['img'], tile['idx']
        tile_filename = f"{prefix}_tile_{idx[0]}_{idx[1]}.png"
        tile_path = os.path.join(tiles_dir, tile_filename)  # Corrected to use tiles_dir
        Image.fromarray(img).save(tile_path)


def process_svs_slide(slide_path, output_dir, tile_size, level, remove_background=False):
    slide_name = os.path.splitext(os.path.basename(slide_path))[0]
    slide_output_dir = os.path.join(output_dir, slide_name)
    tiles_dir = os.path.join(slide_output_dir, "tiles")
    masks_dir = os.path.join(slide_output_dir, "masks")
    os.makedirs(tiles_dir, exist_ok=True)
    os.makedirs(masks_dir, exist_ok=True)
    mask_path = os.path.join(masks_dir, f"{slide_name}_mask.png") 

    if not os.path.exists(slide_output_dir):
        os.makedirs(slide_output_dir, exist_ok=True)

    slide = openslide.OpenSlide(slide_path)
    mpp_x = float(slide.properties.get('openslide.mpp-x', 0))
    mpp_y = float(slide.properties.get('openslide.mpp-y', 0))
    downsample = slide.level_downsamples[level]
    mpp_x_level = mpp_x * downsample
    mpp_y_level = mpp_y * downsample
    print(f"Using level {level} with MPP X: {mpp_x_level}, MPP Y: {mpp_y_level}")

    # Here, we extract tiles from the desired level
    wsi = WSIReader.open(slide_path)
    wsi_overview = wsi.slide_thumbnail(resolution=mpp_x_level, units='mpp',)

    if remove_background:
        tissue_mask = tissue_segmentation(wsi_overview, use_otsu=True)
        mask_path = os.path.join(masks_dir, f"{slide_name}_mask.png")
        save_tissue_mask(tissue_mask, mask_path) 
        tiles = get_tissue_tiles(wsi_overview, tile_size, tissue_mask)
    else:
        tiles = get_tiles(wsi_overview, tile_size)

    save_tiles(tiles, tiles_dir, slide_name)
    print(f"Processed SVS slide: {slide_name}")

def read_tiff(path):
    try:
        with tifffile.TiffFile(path) as tif:
            return tif.asarray()
    except Exception as e:
        print(f"Error reading {path}: {e}")
        return None
        
def process_tiff_slide(slide_path, output_dir, tile_size, remove_background=False):
    slide_name = os.path.splitext(os.path.basename(slide_path))[0]
    slide_output_dir = os.path.join(output_dir, slide_name)
    tiles_dir = os.path.join(slide_output_dir, "tiles")
    masks_dir = os.path.join(slide_output_dir, "masks")
    # Ensure the directories exist
    os.makedirs(tiles_dir, exist_ok=True)
    os.makedirs(masks_dir, exist_ok=True)

    # Define the path for the tissue mask file
    mask_path = os.path.join(masks_dir, f"{slide_name}_mask.png")

    # Read the TIFF image
    image = read_tiff(slide_path)
    
    # Tissue segmentation and tile generation
    if remove_background:
        print('Removing background')
        tissue_mask = tissue_segmentation(image, use_otsu=True)
        print(f"Tissue mask unique values: {np.unique(tissue_mask)}")  
        # Save the tissue mask using the correct mask_path
        save_tissue_mask(tissue_mask, mask_path)
        tiles = get_tissue_tiles(image, tile_size, tissue_mask)
    else:
        tiles = get_tiles(image, tile_size)

    # Save generated tiles
    save_tiles(tiles, tiles_dir, slide_name)
    print(f"Processed TIFF slide: {slide_name}")


def main(slide_folder, output_dir, tile_size, level, slide_ext, remove_background):
    slides = [f for f in os.listdir(slide_folder) if f.endswith(slide_ext)]
    total_slides = len(slides)
    
    for i, filename in enumerate(tqdm(slides, desc="Processing slides"), start=1):
        slide_path = os.path.join(slide_folder, filename)
        print(f"Processing slide {i}/{total_slides}: {filename}")
        if slide_ext == ".svs":
            process_svs_slide(slide_path, output_dir, tile_size, level, remove_background)
        elif slide_ext in [".tiff", ".tif"]:
            process_tiff_slide(slide_path, output_dir, tile_size, remove_background)

if __name__ == "__main__":
    parser = ArgumentParser(description="Generate and save tiles from whole-slide images.")
    parser.add_argument("--slide_folder", type=str, required=True, help="Folder containing whole-slide images.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save the tiles.")
    parser.add_argument("--remove_background", action='store_true', help="Enable background removal in tissue segmentation. Tiles will be generated from tissue areas only.")
    parser.add_argument("--tile_size", type=int, default=1024, help="Size of the square tiles.")
    parser.add_argument("--level", type=int, default=0, help="Slide level to use for tiling.")
    parser.add_argument("--slide_ext", type=str, choices=[".svs", ".tiff", ".tif"], required=True, help="Slide file extension.")
    
    args = parser.parse_args()
    
    main(args.slide_folder, args.output_dir, args.tile_size, args.level, args.slide_ext, args.remove_background)
