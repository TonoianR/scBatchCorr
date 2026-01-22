#' SCT-based RPCA integration followed by PCA (Adaptive Version)
#'
#' Performs SCT normalization, RPCA-based integration, PCA, and ElbowPlot
#' with automatic parameter scaling for small cell counts.
#'
#' @param obj A Seurat object.
#' @param name Character scalar used as filename prefix.
#' @param nfeatures Number of integration features.
#' @param npcs Number of principal components (max target).
#' @param out_dir Base output directory.
#' @param seed Random seed for reproducibility (default: NULL).
#' @param min_cells Minimum cells required per sample to include in integration (default: 10).
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
    min_cells = 10
) {
  
  # -------------------- validation --------------------
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (!"orig.ident" %in% colnames(obj@meta.data)) stop("`orig.ident` column not found.")
  if (!dir.exists(out_dir)) stop("`out_dir` does not exist: ", out_dir)
  
  if (!is.null(seed)) set.seed(seed)
  
  # -------------------- directories --------------------
  out_outputs <- file.path(out_dir, "Outputs")
  out_plots   <- file.path(out_dir, "Plots")
  dir.create(out_outputs, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
  
  message(">>> Starting Adaptive SCT + RPCA integration for: ", name)
  
  # ---- 1. HARD RESET & CLEANING ----
  DefaultAssay(obj) <- "RNA"
  # Seurat v5 compatibility: extract counts and rebuild
  rna_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  obj[["RNA"]] <- CreateAssayObject(counts = rna_counts)
  if ("SCT" %in% names(obj@assays)) obj[["SCT"]] <- NULL
  VariableFeatures(obj) <- character(0)
  
  # ---- 2. SPLIT & FILTER SAMPLES ----
  obj_list <- SplitObject(obj, split.by = "orig.ident")
  
  # Filter out samples that are too small for meaningful normalization/PCA
  cell_counts <- sapply(obj_list, ncol)
  valid_samples <- names(cell_counts)[cell_counts >= min_cells]
  
  if (length(valid_samples) < length(obj_list)) {
    dropped <- setdiff(names(obj_list), valid_samples)
    message(">>> Warning: Dropping samples with < ", min_cells, " cells: ", paste(dropped, collapse = ", "))
    obj_list <- obj_list[valid_samples]
  }
  
  if (length(obj_list) < 2) {
    stop("Integration requires at least 2 samples with more than ", min_cells, " cells.")
  }
  
  # ---- 3. ADAPTIVE SCT & PCA PER SAMPLE ----
  message(">>> Running SCTransform and PCA on individual samples...")
  obj_list <- lapply(obj_list, function(x) {
    # SCTransform with glmGamPoi for efficiency
    x <- SCTransform(x, vst.flavor = "v2", method = "glmGamPoi", verbose = FALSE)
    
    # DYNAMIC PC CALCULATION:
    # Cannot compute more PCs than (cells - 1). 
    # We take the minimum of your target (50) and the physical limit of the subset.
    current_npcs <- min(npcs, (ncol(x) - 1))
    
    if (current_npcs > 1) {
      x <- RunPCA(x, npcs = current_npcs, verbose = FALSE)
    }
    return(x)
  })
  
  invisible(gc(full = TRUE))
  
  # ---- 4. RPCA INTEGRATION ----
  # Ensure target integration features do not exceed total genes in SCT assay
  actual_features <- min(nfeatures, nrow(obj_list[[1]][["SCT"]]))
  
  features <- SelectIntegrationFeatures(
    object.list = obj_list, 
    nfeatures = actual_features
  )
  
  obj_list <- PrepSCTIntegration(
    object.list = obj_list, 
    anchor.features = features
  )
  
  # DYNAMIC K.ANCHOR: 
  # FindIntegrationAnchors default k.anchor is 5. If a sample has very few cells,
  # we must reduce this to avoid "out of bounds" errors.
  min_cells_found <- min(sapply(obj_list, ncol))
  adaptive_k_anchor <- min(5, min_cells_found - 1)
  
  message(">>> Finding anchors (k.anchor set to ", adaptive_k_anchor, ")")
  anchors <- FindIntegrationAnchors(
    object.list = obj_list,
    normalization.method = "SCT",
    anchor.features = features,
    reduction = "rpca",
    k.anchor = adaptive_k_anchor
  )
  
  integrated <- IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT"
  )
  
  # ---- 5. FINAL REDUCTION & EXPORT ----
  DefaultAssay(integrated) <- "integrated"
  
  # Ensure final PCA doesn't exceed total integrated cell count
  final_npcs <- min(npcs, (ncol(integrated) - 1))
  integrated <- RunPCA(integrated, npcs = final_npcs, verbose = FALSE)
  
  # Save object
  qs::qsave(
    integrated, 
    file.path(out_outputs, paste0(name, "_integrated_after_pca.qs"))
  )
  
  # Generate Elbow Plot
  elbow_plot <- ElbowPlot(integrated, ndims = final_npcs)
  ggplot2::ggsave(
    filename = file.path(out_plots, paste0(name, "_ElbowPlot.png")),
    plot = elbow_plot, width = 7, height = 5, dpi = 300
  )
  
  message(">>> Integration successfully completed for ", name)
  return(integrated)
}