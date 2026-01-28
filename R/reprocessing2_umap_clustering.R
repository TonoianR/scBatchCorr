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
    dims = 1:30, 
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
  DefaultAssay(object) <- "SCT"
  
  # FORCE the graph name to be consistent
  graph_name <- "sct_snn"
  
  message(">>> Finding Neighbors using integrated PCA...")
  object <- FindNeighbors(
    object,
    dims = dims,
    reduction = "pca", 
    graph.name = c("sct_nn", graph_name), # Explicitly naming the graphs
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
      graph.name = graph_name, 
      algorithm = 1,
      verbose = FALSE
    )
    
    # --- CASE-INSENSITIVE SMART SEARCH ---
    all_cols <- colnames(object@meta.data)
    res_pattern <- paste0("res.", res, "$")
    
    # 1. Find all columns ending in this resolution (ignoring case)
    potential_cols <- all_cols[grep(res_pattern, all_cols, ignore.case = TRUE)]
    
    # 2. Filter for columns containing 'SCT' or 'sct' (ignoring case)
    # This prevents grabbing old 'integrated' columns
    sct_cols <- potential_cols[grep("SCT", potential_cols, ignore.case = TRUE)]
    
    if (length(sct_cols) > 0) {
      res_col <- tail(sct_cols, 1) # Grab the newest SCT-based column
    } else if (length(potential_cols) > 0) {
      res_col <- tail(potential_cols, 1) # Fallback to any match
    } else {
      res_col <- "seurat_clusters" # Final emergency fallback
    }
    # -------------------------------------
    
    message("    [Plotting] Successfully identified metadata column: ", res_col)
    
    plots[[as.character(res)]] <- DimPlot(
      object,
      group.by = res_col,
      label = TRUE,
      repel = TRUE,
      pt.size = 0.3,
      shuffle = TRUE
    ) + ggplot2::ggtitle(paste("Res", res)) + NoLegend()
  }
  
  # -------------------- saving and plotting --------------------
  comparison <- patchwork::wrap_plots(plots, ncol = 3)
  ggplot2::ggsave(
    filename = file.path(out_plots, paste0(name, "_UMAP_by_resolution.png")),
    plot = comparison, width = 16, height = 9, dpi = 300
  )
  
  p_mix <- DimPlot(object, group.by = "orig.ident", pt.size = 0.3, shuffle = TRUE)
  ggplot2::ggsave(
    filename = file.path(out_plots, paste0(name, "_UMAP_by_sample.png")),
    plot = p_mix, width = 10, height = 7, dpi = 300
  )
  
  message(">>> Saving final object...")
  qs::qsave(object, file.path(out_outputs, paste0(name, "_final.qs")))
  
  return(object)
}