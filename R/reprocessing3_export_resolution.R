#' Export final UMAP for a selected clustering resolution
#'
#' Sets active identities to a chosen clustering resolution and exports
#' a labeled UMAP plot.
#'
#' @param object A Seurat object with UMAP and clustering metadata.
#' @param name Character scalar used as filename prefix.
#' @param resolution Numeric clustering resolution to export (e.g. 0.8 or 1).
#' @param out_dir Base output directory (default: current directory).
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Plot resolution in dots per inch.
#'
#' @return The Seurat object with updated active identities.
#'
#' @export
reprocessing3_export_resolution <- function(
    object,
    name,
    resolution,
    out_dir = ".",
    width = 8,
    height = 6,
    dpi = 300
) {
  
  # -------------------- validation --------------------
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }
  
  if (!"umap" %in% names(object@reductions)) {
    stop("UMAP reduction not found. Run UMAP before exporting.")
  }
  
  if (!is.character(name) || length(name) != 1) {
    stop("`name` must be a single character string.")
  }
  
  res_col <- paste0("integrated_snn_res.", resolution)
  
  if (!res_col %in% colnames(object@meta.data)) {
    stop(
      "Resolution column not found in metadata: ",
      res_col
    )
  }
  
  # -------------------- directories --------------------
  out_plots <- file.path(out_dir, "Plots")
  dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
  
  message(">>> Exporting UMAP for resolution ", resolution, " (", name, ")")
  
  # -------------------- set active identities --------------------
  Idents(object) <- factor(object[[res_col]][, 1])
  attr(object@active.ident, "name") <- res_col
  
  # -------------------- plot --------------------
  p <- DimPlot(
    object,
    reduction = "umap",
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) +
    ggplot2::ggtitle(
      paste0(name, " — resolution ", resolution)
    )
  
  ggplot2::ggsave(
    filename = file.path(
      out_plots,
      paste0(name, "_UMAP_resolution_", resolution, ".png")
    ),
    plot = p,
    width = width,
    height = height,
    dpi = dpi
  )
  
  message(">>> Final UMAP exported.")
  
  return(object)
}

#' Example of usage
#' # Choose resolution based on multi-resolution UMAPs
#' AC_hem <- reprocessing3_export_resolution(
#' object     = AC_hem,
#' name       = "AC_hem",
#' resolution = 1
#' )