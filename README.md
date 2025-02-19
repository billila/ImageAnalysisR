# ImageAnalysisR

Code repo for **ImageAnalysisR** project.

## Progression update
![Progress](https://progress-bar.xyz/41/?title=json)
![h5ad](https://progress-bar.xyz/1/?title=h5ad)


## Case Study: Ovarian Cancer 
file html: [Ovarian Cancer Purity Analysis](https://html-preview.github.io/?url=https://github.com/billila/ImageAnalysisR/blob/main/case_study_OV/ov_purity_analysis.html)

## How to download TCGA image in R
To populate the `wsi_image` folder with TCGA images, follow these steps:
```r
## Make sure BiocManager is installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

## Make sure GenomicDataCommons is installed
# BiocManager::install("GenomicDataCommons")
library("GenomicDataCommons")

# File IDs to download
file_ids <- c(
# Insert your file IDs here
  )

# Download files
lapply(file_ids, gdcdata)
```

## HoVerNet 
To use HoVerNet for segmentation:

1. Set up the Conda environment as described in the [HoverNet GitHub repository](https://github.com/vqdang/hover_net).
2. Run the script [run_infer.sh](code/hovernet/slurm_run_hovernet.sh) in a GPU-enabled environment.
The output will be a JSON file containing the nuclei segmentation and classification results.

## Feature Extraction on HoVerNet segmentation output (JSON file)

To extract features from the HoVerNet segmentation output (JSON file):

* Use the scripts [hovernet_py.qmd](code/hovernet_py.qmd) to process the JSON file and generate an .h5ad file.

## Example file 
In this [folder](example) you can find example for:

* HoVerNet output in `.json`, mask in `.png`, thumb in `.png`
* Features Extraction output in `.h5ad`.

## Literature and useful links:

### General Understanding and open problems in Digital Pathology
- [Annual Review on Digital Pathology](https://www.annualreviews.org/content/journals/10.1146/annurev-cancerbio-062822-010523?crawler=true)  
- [Scientific Review on Digital Pathology](https://www.sciencedirect.com/science/article/pii/S2153353923001712)  

### Systematic evaluation for cell segmentation algorithms
- [arXiv: Systematic Evaluation of Cell Segmentation Algorithms](https://arxiv.org/abs/2310.18689)

### HoVerNet and Related Papers
- Paper using HoVerNet outputs: [PMCID: PMC10985477](https://pmc.ncbi.nlm.nih.gov/articles/PMC10985477/pdf/347.pdf)  
- Original HoVerNet paper: [HoVerNet](https://www.sciencedirect.com/science/article/pii/S1361841519301045)  
- HoVerNet GitHub repository: [GitHub - HoverNet](https://github.com/vqdang/hover_net)
- HoVerNext paper: [HoVerNext](https://openreview.net/pdf?id=3vmB43oqIO)  
- HoVerNext GitHub repository: [GitHub - HoverNext](https://github.com/digitalpathologybern/hover_next_inference)

### Benchmarks for foundation models in Digital Pathology
- [UNI Benchmarks](https://github.com/mahmoodlab/UNI?tab=readme-ov-file#slide-benchmarks)

### Foundation model papers
- **UNI:**  
  - Paper: [Nature Medicine - UNI](https://www.nature.com/articles/s41591-024-02857-3)  
  - GitHub: [GitHub - UNI](https://github.com/mahmoodlab/UNI)
- **Prov-GigaPath:**  
  - Paper: [Nature - Prov-GigaPath](https://www.nature.com/articles/s41586-024-07441-w)  
  - GitHub: [GitHub - Prov-GigaPath](https://github.com/prov-gigapath/prov-gigapath)

### Public H&E image databases
- [TCGA Repository](https://portal.gdc.cancer.gov/analysis_page?app=Downloads)  
- [TCIA Repository](https://pathdb.cancerimagingarchive.net/eaglescope/dist/?configurl=%2Fsystem%2Ffiles%2Fcollectionmetadata%2F202401%2Fcohort_builder_01-27-2024.json)

### GitHub Repositories
- [ImageAnalysisR](https://github.com/billila/ImageAnalysisR) (Shared Code/Material)  
- [imageTCGA](https://github.com/billila/imageTCGA) (Shiny Package)
- [TCIAAPI](https://github.com/billila/TCIAAPI) (R Package)

### R Packages
- **EBImage:**  
  - [EBImage](https://www.bioconductor.org/packages/release/bioc/html/EBImage.html)
  - GitHub: [EBImage](https://github.com/aoles/EBImage)
- **MSIreg:**  
  - Paper: [Oxford Bioinformatics - MSIreg](https://academic.oup.com/bioinformatics/article/40/11/btae624/7825357)  
  - GitHub: [MSIreg](https://github.com/sslakkimsetty/msireg)

