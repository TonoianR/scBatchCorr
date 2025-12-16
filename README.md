# scBatchCorr

**scBatchCorr** is an R package providing modular utilities for batch correction,
integration, clustering, and inspection of single-cell RNA-seq data using Seurat.

The package is designed to separate:
- reusable computational methods (package)
- project-specific biological decisions (analysis scripts)

## Installation

```r
# install.packages("remotes")
remotes::install_github("TonoianR/scBatchCorr")
```

**Typical workflow**
```{r }
library(Seurat)
library(scBatchCorr)

# Step 1: Integration and PCA
obj <- reprocessing1_integration_pca(obj, "sample")

# Step 2: UMAP and clustering
obj <- reprocessing2_umap_clustering(obj, "sample")

# Step 3: Choose final resolution
obj <- reprocessing3_export_resolution(obj, "sample", resolution = 1)

# Step 4: Inspect and summarize
summarize_umap(obj, "sample")
```
**Package structure**
```r
scBatchCorr/
│
├── R/
│   ├── reprocessing1_integration_pca.R        # Function for integration and PCA
│   ├── reprocessing2_umap_clustering.R        # Function for UMAP and clustering
│   ├── reprocessing3_export_resolution.R     # Function for exporting final resolution
│   ├── reprocessing4_summarize_umap.R        # Function for summarizing UMAP and cell clustering
│
└── DESCRIPTION
└── NAMESPACE
```

**Functions' descriptions
	•	R/reprocessing1_integration_pca.R: This function performs full single-cell data integration and PCA. It normalizes the data, finds integration anchors, and runs PCA on the integrated object.
	•	R/reprocessing2_umap_clustering.R: This function takes the integrated object and performs UMAP dimensionality reduction and multi-resolution clustering. It also generates UMAP plots at different resolutions and checks for batch mixing.
	•	R/reprocessing3_export_resolution.R: This function allows you to select a resolution from clustering results and export the final UMAP for a specified resolution.
	•	R/reprocessing4_summarize_umap.R: This function inspects the integrated data, visualizes UMAP, lists unique clusters, and saves a summary table of the clustering results for downstream analysis.

**Author**
Robert Tonoian
