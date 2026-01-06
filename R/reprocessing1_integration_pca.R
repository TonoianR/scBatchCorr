#' SCT-based RPCA integration followed by PCA
#'
#' Performs SCT normalization, RPCA-based integration, PCA, and ElbowPlot
#' generation for a Seurat object split by `orig.ident`.
#'
#' @param obj A Seurat object.
#' @param name Character scalar used as filename prefix.
#' @param nfeatures Number of integration features.
#' @param npcs Number of principal components.
#' @param out_dir Base output directory.
#' @param seed Random seed for reproducibility (default: NULL).
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
    seed = NULL
) {
  
  # -------------------- validation --------------------
  if (inherits(obj, "Seurat")) {
    
    if (!"orig.ident" %in% colnames(obj@meta.data)) {
      stop("`orig.ident` column not found in Seurat object metadata.")
    }
    
    DefaultAssay(obj) <- "RNA"
    obj_list <- SplitObject(obj, split.by = "orig.ident")
    
  } else if (is.list(obj) && all(vapply(obj, inherits, logical(1), "Seurat"))) {
    
    obj_list <- obj
    
  } else {
    stop("`obj` must be a Seurat object or a list of Seurat objects.")
  }
  
  if (!is.character(name) || length(name) != 1) {
    stop("`name` must be a single character string.")
  }
  
  if (!is.character(out_dir) || length(out_dir) != 1) {
    stop("`out_dir` must be a single directory path.")
  }
  
  if (!dir.exists(out_dir)) {
    stop("`out_dir` does not exist: ", out_dir)
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # -------------------- directories --------------------
  out_outputs <- file.path(out_dir, "Outputs")
  out_plots   <- file.path(out_dir, "Plots")
  
  dir.create(out_outputs, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
  
  message(">>> Starting SCT + RPCA integration for ", name)
  message(">>> Output directory: ", normalizePath(out_dir))
  
  # -------------------- SCT normalization --------------------
  obj_list <- lapply(
    obj_list,
    function(x) {
      
      DefaultAssay(x) <- "RNA"
      
      # CRITICAL: remove any existing SCT state
      if ("SCT" %in% names(x@assays)) {
        x[["SCT"]] <- NULL
      }
      
      SCTransform(
        x,
        vst.flavor = "v2",
        method = "glmGamPoi",
        verbose = FALSE
      )
    }
  )
  
  rm(obj)
  invisible(gc(full = TRUE))
  
  # -------------------- PCA per subset --------------------
  obj_list <- lapply(obj_list, function(x) {
    if (!"pca" %in% names(x@reductions)) {
      x <- RunPCA(x, verbose = FALSE)
    }
    x
  })
  
  invisible(gc(full = TRUE))
  
  # -------------------- RPCA integration --------------------
  features <- SelectIntegrationFeatures(
    object.list = obj_list,
    nfeatures = nfeatures
  )
  
  obj_list <- PrepSCTIntegration(
    object.list = obj_list,
    anchor.features = features
  )
  
  anchors <- FindIntegrationAnchors(
    object.list = obj_list,
    normalization.method = "SCT",
    anchor.features = features,
    reduction = "rpca"
  )
  
  rm(obj_list)
  invisible(gc(full = TRUE))
  
  integrated <- IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT"
  )
  
  DefaultAssay(integrated) <- "integrated"
  
  # -------------------- PCA + ElbowPlot --------------------
  integrated <- RunPCA(
    integrated,
    npcs = npcs,
    verbose = FALSE
  )
  
  # SAVE integrated object after PCA  ✅
  qs::qsave(
    integrated,
    file.path(out_outputs, paste0(name, "_integrated_after_pca.qs"))
  )
  
  elbow_plot <- ElbowPlot(
    integrated,
    ndims = npcs
  )
  
  ggplot2::ggsave(
    filename = file.path(out_plots, paste0(name, "_ElbowPlot.png")),
    plot = elbow_plot,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  invisible(gc(full = TRUE))
  
  message(">>> Integration and PCA completed for ", name)
  message(">>> Integrated object saved to: ",
          file.path(out_outputs, paste0(name, "_integrated_after_pca.qs")))
  message(">>> Elbow plot saved to: ", file.path(out_plots, paste0(name, "_ElbowPlot.png")))
  
  return(integrated)
}
