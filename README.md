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

**Functions' descriptions**
```r
I. reprocessing1_integration_pca

This function integrates multiple datasets and performs PCA. It normalizes the data using SCTransform and finds integration anchors using RPCA. It then runs PCA on the integrated object to reduce dimensionality.

Parameters:
	•	obj: Seurat object containing single-cell RNA-seq data.
	•	name: A string representing the name of the sample for saving output files.
	•	nfeatures: Number of features to select for integration (default: 3000).
	•	npcs: Number of principal components to use for downstream analysis (default: 50).
	•	base_path: The working directory for saving output files (default: current working directory).

Output:
	•	Integrated Seurat object with PCA results.
	•	Elbow plot saved as .png.

II. reprocessing2_umap_clustering

This function performs UMAP embedding and clustering at multiple resolutions using the Louvain algorithm. It visualizes the results and checks for batch mixing by plotting the UMAP colored by sample identity.

Parameters:
	•	object: Seurat object.
	•	name: Sample name for saving output.
	•	dims: Principal components to use for UMAP (default: 1:20).
	•	resolutions: List of resolutions to cluster (default: c(0.2, 0.4, 0.5, 0.6, 0.8, 1.0)).
	•	base_path: Working directory for output files.

Output:
	•	UMAP plot saved for each resolution.
	•	Final clustered Seurat object saved for downstream analyses.

III. reprocessing3_export_resolution

This function exports a final resolution and generates a UMAP plot labeled by clusters.

Parameters:
	•	object: Seurat object.
	•	name: Sample name for saving output.
	•	resolution: Resolution to choose from clustering results (default: 1).

Output:
	•	UMAP plot saved for the selected resolution.

IV. reprocessing4_summarize_umap

This function summarizes clustering results, visualizes UMAP, and generates a table with cluster counts and metadata for cells.

Parameters:
	•	obj: Seurat object.
	•	label: Sample name for saving output.
	•	reduction: The reduction method to use for plotting (default: “umap”).
	•	base_path: Directory to save summary and output files.

Output:
	•	Summary table saved as a .csv file.
	•	Cluster summary and per-sample data printed to the console.

```
**Author**
Robert Tonoian
