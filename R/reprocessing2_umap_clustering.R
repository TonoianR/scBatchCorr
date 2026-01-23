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
    dims = 1:30, # Increased default to catch subtle HSC signals
    resolutions = c(0.2, 0.4, 0.5, 0.6, 0.8, 1.0),
    out_dir = ".",
    seed = NULL
) {
  
  # -------------------- validation --------------------
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.")
  if (!"pca" %in% names(object@reductions)) stop("PCA reduction not found.")
  if (!is.null(seed)) set.seed(seed)
  
  out_outputs <- file.path(out_dir, "Outputs")
  out_plots   <- file.path(out_dir, "Plots")
  dir.create(out_outputs, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
  
  message(">>> Running UMAP and clustering for ", name)
  
  # -------------------- UMAP & Neighbors --------------------
  # USE SCT: It's the best assay for high-resolution blood/stem cell work
  DefaultAssay(object) <- "SCT"
  
  message(">>> Finding Neighbors using integrated PCA...")
  object <- FindNeighbors(
    object,
    dims = dims,
    reduction = "pca", # This is our integrated space from Function 1
    graph.name = "sct_snn", 
    verbose = FALSE
  )
  
  message(">>> Running UMAP...")
  object <- RunUMAP(
    object,
    dims = dims,
    reduction = "pca",
    verbose = FALSE
  )
  
  # -------------------- multi-resolution clustering --------------------
  plots <- vector("list", length(resolutions))
  names(plots) <- as.character(resolutions)
  
  for (res in resolutions) {
    if (!is.null(seed)) set.seed(seed)
    
    message("    Clustering at resolution: ", res)
    object <- FindClusters(
      object,
      resolution = res,
      graph.name = "sct_snn", # Matching the graph name above
      algorithm = 1,
      verbose = FALSE
    )
    
    # In Seurat v5, clustering columns are named: [Assay]_snn_res.[Res]
    colname <- paste0("SCT_snn_res.", res)
    
    plots[[as.character(res)]] <- DimPlot(
      object,
      group.by = colname,
      label = TRUE,
      repel = TRUE,
      pt.size = 0.3,
      shuffle = TRUE
    ) + ggplot2::ggtitle(paste("Res", res)) + NoLegend()
  }
  
  # Save plots
  comparison <- patchwork::wrap_plots(plots, ncol = 3)
  ggplot2::ggsave(
    filename = file.path(out_plots, paste0(name, "_UMAP_by_resolution.png")),
    plot = comparison, width = 16, height = 9, dpi = 300
  )
  
  # Batch mixing check
  p_mix <- DimPlot(object, group.by = "orig.ident", pt.size = 0.3, shuffle = TRUE)
  ggplot2::ggsave(
    filename = file.path(out_plots, paste0(name, "_UMAP_by_sample.png")),
    plot = p_mix, width = 10, height = 7, dpi = 300
  )
  
  message(">>> Saving final object...")
  qs::qsave(object, file.path(out_outputs, paste0(name, "_final.qs")))
  
  return(object)
}