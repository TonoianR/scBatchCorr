#' Inspect integrated Seurat object and summarize clustering
#'
#' Generates a UMAP for inspection and produces summary tables
#' describing cluster composition and per-sample cell counts.
#'
#' @param object A Seurat object with UMAP and active identities set.
#' @param name Character scalar used as filename prefix.
#' @param reduction Dimensional reduction to visualize (default: "umap").
#' @param out_dir Base output directory (default: current directory).
#' @param export Logical; whether to write summary CSV files (default: TRUE).
#'
#' @return A list containing summary tables (returned invisibly).
#'
#' @export
reprocessing4_summarize_umap <- function(
    object,
    name,
    reduction = "umap",
    out_dir = ".",
    export = TRUE
) {
  
  # -------------------- validation --------------------
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }
  
  if (!reduction %in% names(object@reductions)) {
    stop("Reduction not found: ", reduction)
  }
  
  if (!is.character(name) || length(name) != 1) {
    stop("`name` must be a single character string.")
  }
  
  if (length(Idents(object)) == 0) {
    stop("Active identities are not set. Run Function 3 first.")
  }
  
  # -------------------- directories --------------------
  out_plots   <- file.path(out_dir, "Plots")
  out_outputs <- file.path(out_dir, "Outputs")
  
  dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_outputs, recursive = TRUE, showWarnings = FALSE)
  
  message(">>> Inspecting integrated object: ", name)
  
  # -------------------- UMAP inspection plot --------------------
  p <- DimPlot(
    object,
    reduction = reduction,
    label = TRUE,
    repel = TRUE,
    raster = FALSE,
    pt.size = 0.3
  ) +
    ggplot2::ggtitle(paste0(name, " — ", reduction))
  
  print(p)
  
  # -------------------- cluster summary --------------------
  cluster_counts <- as.data.frame(table(Idents(object)))
  colnames(cluster_counts) <- c("Cluster", "Cells")
  
  cluster_summary <- cluster_counts |>
    dplyr::arrange(desc(Cells))
  
  # -------------------- general summary --------------------
  general_summary <- data.frame(
    Metric = c(
      "Total cells",
      "Total clusters",
      "Number of samples"
    ),
    Value = c(
      ncol(object),
      length(unique(Idents(object))),
      length(unique(object$orig.ident))
    )
  )
  
  # -------------------- per-sample summary --------------------
  sample_summary <- object@meta.data |>
    dplyr::group_by(orig.ident) |>
    dplyr::summarise(
      Cells = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(orig.ident)
  
  # -------------------- optional export --------------------
  if (export) {
    
    readr::write_csv(
      general_summary,
      file.path(out_outputs, paste0(name, "_summary_general.csv"))
    )
    
    readr::write_csv(
      cluster_summary,
      file.path(out_outputs, paste0(name, "_summary_clusters.csv"))
    )
    
    readr::write_csv(
      sample_summary,
      file.path(out_outputs, paste0(name, "_summary_samples.csv"))
    )
    
    message(">>> Summary tables exported to Outputs/")
  }
  
  # -------------------- console preview --------------------
  message("\n>>> General summary:")
  print(knitr::kable(general_summary, align = "l"))
  
  message("\n>>> Cluster composition:")
  print(knitr::kable(cluster_summary, align = "l"))
  
  message("\n>>> Per-sample summary:")
  print(knitr::kable(sample_summary, align = "l"))
  
  return(invisible(list(
    general  = general_summary,
    clusters = cluster_summary,
    samples  = sample_summary
  )))
}

#' Example of usage
#' AC_hem_summary <- reprocessing4_summarize_umap(
#' object = AC_hem,
#' name   = "AC_hem"
#' )