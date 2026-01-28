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
#' Export final UMAP for a selected clustering resolution
#'
#' Sets active identities to a chosen clustering resolution and exports
#' a labeled UMAP plot, specifically handling SCT or integrated prefixes.
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
  
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.")
  
  # -------------------- Priority Search --------------------
  all_cols <- colnames(object@meta.data)
  res_pattern <- paste0("res.", resolution, "$")
  
  # Find all columns matching the resolution
  potential_cols <- all_cols[grep(res_pattern, all_cols, ignore.case = TRUE)]
  
  # STRICT PRIORITY: 1. SCT (Integrated) -> 2. RNA (Fallback)
  sct_match <- potential_cols[grep("sct", potential_cols, ignore.case = TRUE)]
  
  if (length(sct_match) > 0) {
    res_col <- tail(sct_match, 1)
    message(">>> Success: Found Integrated SCT column: ", res_col)
  } else if (length(potential_cols) > 0) {
    res_col <- tail(potential_cols, 1)
    message(">>> Warning: SCT column not found. Using fallback: ", res_col)
  } else {
    stop("Could not find any metadata column matching resolution: ", resolution)
  }
  
  # -------------------- Export Logic --------------------
  out_plots <- file.path(out_dir, "Plots")
  dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
  
  # Lock identities
  object <- SetIdent(object, value = res_col)
  
  p <- DimPlot(
    object,
    reduction = "umap",
    label = TRUE,
    repel = TRUE,
    pt.size = 0.5
  ) +
    ggplot2::ggtitle(paste0(name, " — Resolution ", resolution, " (Integrated)")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  ggplot2::ggsave(
    filename = file.path(out_plots, paste0(name, "_UMAP_final_res_", resolution, ".png")),
    plot = p, width = width, height = height, dpi = dpi
  )
  
  return(object)
}