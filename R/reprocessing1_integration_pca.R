#' SCT-based RPCA integration followed by PCA
#'
#' Performs SCT normalization, RPCA-based integration, PCA, and ElbowPlot
#' generation for a Seurat object split by `orig.ident`.
#'
#' @param obj A Seurat object.
#' @param name Character scalar used as filename prefix.
#' @param nfeatures Number of integration features.
#' @param npcs Number of principal components.
#' @param out_dir Base output directory (default: current directory).
#' @param seed Random seed for reproducibility (default: NULL = no seed set).
#'
#' @return A Seurat object after integration and PCA.
#'
#' @export
reprocessing1_integration_pca <- function(
    obj,
    name,
    nfeatures = 3000,
    npcs = 50,
    out_dir = ".",
    seed = NULL
) {
  
  # -------------------- validation --------------------
  if (!inherits(obj, "Seurat")) {
    stop("`obj` must be a Seurat object.")
  }
  
  if (!"orig.ident" %in% colnames(obj@meta.data)) {
    stop("`orig.ident` column not found in Seurat object metadata.")
  }
  
  if (!is.character(name) || length(name) != 1) {
    stop("`name` must be a single character string.")
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # -------------------- directories --------------------
  out_outputs <- file.path(out_dir, "Outputs")
  out_plots   <- file.path(out_dir, "Plots")
  
  dir.create(out_outputs, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
  
  message(">>> Starting integration and PCA for ", name)
  
  # -------------------- SCT normalization --------------------
  DefaultAssay(obj) <- "RNA"
  
  obj_list <- SplitObject(obj, split.by = "orig.ident")
  obj_list <- lapply(
    obj_list,
    function(x) {
      SCTransform(
        x,
        vst.flavor = "v2",
        method = "glmGamPoi",
        verbose = FALSE
      )
    }
  )
  
  qs::qsave(
    obj_list,
    file.path(out_outputs, paste0(name, "_list_SCT.qs"))
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
  features <- SelectIntegrationFeatures(obj_list, nfeatures = nfeatures)
  obj_list <- PrepSCTIntegration(obj_list, anchor.features = features)
  
  anchors <- FindIntegrationAnchors(
    object.list = obj_list,
    normalization.method = "SCT",
    anchor.features = features,
    reduction = "rpca"
  )
  
  rm(obj_list)
  invisible(gc(full = TRUE))
  
  integrated <- IntegrateData(
    anchors,
    normalization.method = "SCT"
  )
  
  qs::qsave(
    integrated,
    file.path(out_outputs, paste0(name, "_integrated.qs"))
  )
  
  # -------------------- PCA + ElbowPlot --------------------
  integrated <- RunPCA(integrated, npcs = npcs, verbose = FALSE)
  
  qs::qsave(
    integrated,
    file.path(out_outputs, paste0(name, "_integrated_after_pca.qs"))
  )
  
  elbow_plot <- ElbowPlot(integrated, ndims = npcs)
  
  ggplot2::ggsave(
    filename = file.path(out_plots, paste0(name, "_ElbowPlot.png")),
    plot = elbow_plot,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  invisible(gc(full = TRUE))
  
  message(">>> Integration and PCA completed for ", name)
  
  return(integrated)
}

#' Example of usage
#' integrated_obj <- reprocessing1_integration_pca(
#' obj    = AC_hem,
#' name   = "AC_hem",
#' out_dir = "/path/to/project",
#' seed   = 1234
#' )