#' SCT-based RPCA integration followed by PCA (Seurat v5 Robust Version)
#'
#' Performs SCT normalization, RPCA-based integration, PCA, and ElbowPlot
#' with aggressive cleaning to prevent Seurat v5 layer-mismatch errors.
#'
#' @param obj A Seurat object.
#' @param name Character scalar used as filename prefix.
#' @param nfeatures Number of integration features.
#' @param npcs Number of principal components.
#' @param out_dir Base output directory.
#' @param seed Random seed for reproducibility.
#' @param min_cells Minimum cells required per sample (default: 10).
#'
#' @return A Seurat object after integration and PCA.
#'
#' @export
reprocessing1_integration_pca <- function(
    obj,
    name,
    nfeatures = 3000,
    npcs = 50,
    out_dir,
    seed = NULL,
    min_cells = 5
) {
  if (!is.null(seed)) set.seed(seed)
  out_outputs <- file.path(out_dir, "Outputs"); out_plots <- file.path(out_dir, "Plots")
  dir.create(out_outputs, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
  
  # 2. Reset and STRIP old assays (Fixes the "Different cells" warning)
  DefaultAssay(obj) <- "RNA"
  # Keep only RNA, remove everything else
  for (assay_name in names(obj@assays)) {
    if (assay_name != "RNA") obj[[assay_name]] <- NULL
  }
  obj <- JoinLayers(obj)
  
  # 3. Split and Filter
  obj_list <- SplitObject(obj, split.by = "orig.ident")
  cell_counts <- sapply(obj_list, ncol)
  obj_list <- obj_list[names(cell_counts)[cell_counts >= min_cells]]
  if (length(obj_list) < 2) stop("Not enough samples.")
  
  # 4. Uniform Dimension Logic
  min_size <- min(sapply(obj_list, ncol))
  safe_npcs <- max(2, min(npcs, (min_size - 1)))
  
  # 5. SCT and PCA
  obj_list <- lapply(obj_list, function(x) {
    x <- SCTransform(x, vst.flavor = "v2", method = "glmGamPoi", verbose = FALSE)
    x <- RunPCA(x, npcs = safe_npcs, verbose = FALSE)
    return(x)
  })
  
  # 6. Integration (The "Small Batch" Neighborhood Fix)
  features <- SelectIntegrationFeatures(obj_list, nfeatures = nfeatures)
  obj_list <- PrepSCTIntegration(obj_list, anchor.features = features)
  
  # Neighborhood parameters must be < smallest batch size
  safe_k <- min(5, (min_size - 1))
  
  anchors <- FindIntegrationAnchors(
    object.list = obj_list, 
    normalization.method = "SCT",
    anchor.features = features, 
    dims = 1:safe_npcs, 
    reduction = "rpca", 
    k.anchor = safe_k,
    k.filter = NA # This is the critical fix for the index error
  )
  
  integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  
  # 7. Consolidate RNA
  integrated <- JoinLayers(integrated, assay = "RNA")
  integrated <- NormalizeData(integrated, assay = "RNA", verbose = FALSE)
  integrated <- ScaleData(integrated, assay = "RNA", verbose = FALSE)
  
  # 8. Final PCA
  DefaultAssay(integrated) <- "integrated"
  final_npcs <- min(npcs, (ncol(integrated) - 1))
  integrated <- RunPCA(integrated, npcs = final_npcs, verbose = FALSE)
  
  qs::qsave(integrated, file.path(out_outputs, paste0(name, "_integrated_after_pca.qs")))
  return(integrated)
}