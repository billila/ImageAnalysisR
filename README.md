# ImageAnalysisR

Code repo for ImageAnalysisR project.

## How to download TCGA image in R
How to populate wsi_image folder.
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
## Feature Extraction on HoVerNet segmentation output (JSON file) 


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

### R Packages
- **EBImage:**  
  - [EBImage](https://www.bioconductor.org/packages/release/bioc/html/EBImage.html)
  - GitHub: [EBImage](https://github.com/aoles/EBImage)
- **MSIreg:**  
  - Paper: [Oxford Bioinformatics - MSIreg](https://academic.oup.com/bioinformatics/article/40/11/btae624/7825357)  
  - GitHub: [MSIreg](https://github.com/sslakkimsetty/msireg)

