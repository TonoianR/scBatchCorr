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
    min_cells = 10
) {
  
  # -------------------- Validation & Setup --------------------
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (!"orig.ident" %in% colnames(obj@meta.data)) stop("`orig.ident` column not found.")
  if (!dir.exists(out_dir)) stop("`out_dir` does not exist: ", out_dir)
  if (!is.null(seed)) set.seed(seed)
  
  out_outputs <- file.path(out_dir, "Outputs")
  out_plots   <- file.path(out_dir, "Plots")
  dir.create(out_outputs, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
  
  message(">>> Starting Seurat v5 Robust SCT + RPCA integration for: ", name)
  
  # ---- 1. AGGRESSIVE SANITIZATION (The "Subscript Out of Bounds" Fix) ----
  DefaultAssay(obj) <- "RNA"
  
  # Join all layers (counts.1, counts.2, etc) into a single matrix
  # This is critical in Seurat v5 after subsetting.
  obj <- JoinLayers(obj)
  
  # Extract raw counts and metadata to build a clean 'Reset' object
  raw_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  metadata <- obj@meta.data
  
  # Create a fresh object to break all links to stale SCT assays or PCA slots
  new_obj <- CreateSeuratObject(counts = raw_counts, meta.data = metadata)
  
  # ---- 2. SPLIT & FILTER SAMPLES ----
  obj_list <- SplitObject(new_obj, split.by = "orig.ident")
  
  # Filter out samples below the threshold
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
    # Wrapped in try() to ensure one bad sample doesn't kill the whole pipeline
    x <- tryCatch({
      x <- SCTransform(x, vst.flavor = "v2", method = "glmGamPoi", verbose = FALSE)
      
      # DYNAMIC PC CALCULATION
      current_npcs <- min(npcs, (ncol(x) - 1))
      if (current_npcs > 1) {
        x <- RunPCA(x, npcs = current_npcs, verbose = FALSE)
      }
      return(x)
    }, error = function(e) {
      message(">>> Error in sample ", unique(x$orig.ident), ": ", e$message)
      return(NULL)
    })
    return(x)
  })
  
  # Remove any samples where SCT failed
  obj_list <- Filter(Negate(is.null), obj_list)
  
  if (length(obj_list) < 2) stop("Fewer than 2 samples survived SCTransform.")
  
  invisible(gc(full = TRUE))
  
  # ---- 4. RPCA INTEGRATION ----
  # Ensure feature count doesn't exceed total genes available in the smallest SCT assay
  actual_features <- min(nfeatures, nrow(obj_list[[1]][["SCT"]]))
  
  features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = actual_features)
  obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)
  
  # DYNAMIC K.ANCHOR: adjust for small neighborhoods
  min_cells_found <- min(sapply(obj_list, ncol))
  adaptive_k_anchor <- min(5, min_cells_found - 1)
  
  message(">>> Finding anchors (k.anchor: ", adaptive_k_anchor, ")")
  anchors <- FindIntegrationAnchors(
    object.list = obj_list,
    normalization.method = "SCT",
    anchor.features = features,
    reduction = "rpca",
    k.anchor = adaptive_k_anchor
  )
  
  integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  
  # ---- 5. FINAL REDUCTION & EXPORT ----
  DefaultAssay(integrated) <- "integrated"
  
  final_npcs <- min(npcs, (ncol(integrated) - 1))
  integrated <- RunPCA(integrated, npcs = final_npcs, verbose = FALSE)
  
  # Save object
  qs::qsave(integrated, file.path(out_outputs, paste0(name, "_integrated_after_pca.qs")))
  
  # Elbow Plot
  elbow_plot <- ElbowPlot(integrated, ndims = final_npcs)
  ggplot2::ggsave(filename = file.path(out_plots, paste0(name, "_ElbowPlot.png")),
                  plot = elbow_plot, width = 7, height = 5, dpi = 300)
  
  message(">>> Integration successfully completed for ", name)
  return(integrated)
}