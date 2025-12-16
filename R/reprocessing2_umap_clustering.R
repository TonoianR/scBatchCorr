#' UMAP embedding and multi-resolution clustering
#'
#' Computes UMAP embedding and performs Louvain clustering
#' across multiple resolutions on an integrated Seurat object.
#'
#' @param object A Seurat object with PCA computed on the integrated assay.
#' @param name Character scalar used as filename prefix.
#' @param dims Integer vector of PCA dimensions to use.
#' @param resolutions Numeric vector of clustering resolutions.
#' @param out_dir Base output directory (default: current directory).
#' @param seed Random seed for reproducibility (default: NULL).
#'
#' @return A Seurat object with UMAP and clustering metadata added.
#'
#' @export
reprocessing2_umap_clustering <- function(
    object,
    name,
    dims = 1:20,
    resolutions = c(0.2, 0.4, 0.5, 0.6, 0.8, 1.0),
    out_dir = ".",
    seed = NULL
) {
  
  # -------------------- validation --------------------
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }
  
  if (!"pca" %in% names(object@reductions)) {
    stop("PCA reduction not found. Run PCA before UMAP.")
  }
  
  if (!"orig.ident" %in% colnames(object@meta.data)) {
    stop("`orig.ident` column not found in metadata.")
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
  
  message(">>> Running UMAP and clustering for ", name)
  
  # -------------------- UMAP --------------------
  DefaultAssay(object) <- "integrated"
  
  object <- FindNeighbors(
    object,
    dims = dims,
    graph.name = "integrated_snn",
    verbose = FALSE
  )
  
  object <- RunUMAP(
    object,
    dims = dims,
    reduction = "pca",
    verbose = FALSE
  )
  
  qs::qsave(
    object,
    file.path(out_outputs, paste0(name, "_after_umap.qs"))
  )
  
  # -------------------- multi-resolution clustering --------------------
  plots <- vector("list", length(resolutions))
  names(plots) <- as.character(resolutions)
  
  for (res in resolutions) {
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    object <- FindClusters(
      object,
      resolution = res,
      graph.name = "integrated_snn",
      algorithm = 1,
      verbose = FALSE
    )
    
    colname <- paste0("integrated_snn_res.", res)
    
    plots[[as.character(res)]] <-
      DimPlot(
        object,
        group.by = colname,
        label = TRUE,
        repel = TRUE,
        raster = FALSE,
        pt.size = 0.3,
        shuffle = TRUE
      ) +
      ggplot2::ggtitle(paste("Resolution", res))
  }
  
  comparison <- patchwork::wrap_plots(plots, ncol = 3)
  
  ggplot2::ggsave(
    filename = file.path(out_plots, paste0(name, "_UMAP_by_resolution.png")),
    plot = comparison,
    width = 16,
    height = 9,
    dpi = 300
  )
  
  message(">>> Multi-resolution UMAP saved.")
  
  # -------------------- batch mixing check --------------------
  p_mix <- DimPlot(
    object,
    group.by = "orig.ident",
    raster = FALSE,
    pt.size = 0.3,
    shuffle = TRUE
  )
  
  ggplot2::ggsave(
    filename = file.path(out_plots, paste0(name, "_UMAP_by_sample.png")),
    plot = p_mix,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  qs::qsave(
    object,
    file.path(out_outputs, paste0(name, "_final.qs"))
  )
  
  message(">>> UMAP + clustering completed for ", name)
  
  return(object)
}

#' Example of usage
#' AC_hem <- reprocessing2_umap_clustering(
#' object   = AC_hem,
#' name     = "AC_hem",
#' dims     = 1:25,
#' out_dir  = "/srv/projects/E16_liver",
#' seed     = 1234
#' )