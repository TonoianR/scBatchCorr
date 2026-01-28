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
#' Inspect integrated Seurat object and summarize clustering with Clonal Info
#'
#' @param object A Seurat object with UMAP and active identities set.
#' @param name Character scalar used as filename prefix.
#' @param reduction Dimensional reduction to visualize (default: "umap").
#' @param out_dir Base output directory.
#' @param clone_col The column name in metadata containing clone IDs (default: "cloneid").
#' @param export Logical; whether to write summary CSV files.
#'
#' @export
reprocessing4_summarize_umap <- function(
    object,
    name,
    reduction = "umap",
    out_dir = ".",
    clone_col = "cloneid",
    export = TRUE
) {
  
  # -------------------- validation --------------------
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.")
  if (!reduction %in% names(object@reductions)) stop("Reduction not found: ", reduction)
  if (length(Idents(object)) == 0) stop("Active identities are not set. Run Function 3 first.")
  
  # Ensure clone column exists, if not, create a dummy to avoid crash
  if (!clone_col %in% colnames(object@meta.data)) {
    warning("Clone column '", clone_col, "' not found. Counting as 0.")
    object[[clone_col]] <- NA
  }
  
  out_plots   <- file.path(out_dir, "Plots")
  out_outputs <- file.path(out_dir, "Outputs")
  dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_outputs, recursive = TRUE, showWarnings = FALSE)
  
  message(">>> Summarizing object: ", name)
  
  # Helper: Count non-NA/non-empty clones
  count_clones <- function(x) sum(!is.na(x) & x != "" & x != "None", na.rm = TRUE)
  
  # -------------------- 1. Cluster Summary --------------------
  cluster_summary <- object@meta.data %>%
    dplyr::mutate(CurrentIdent = Idents(object)) %>%
    dplyr::group_by(Cluster = CurrentIdent) %>%
    dplyr::summarise(
      Total_Cells = dplyr::n(),
      Cells_with_CloneID = count_clones(.data[[clone_col]]),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Cluster)
  
  # -------------------- 2. General Summary --------------------
  general_summary <- data.frame(
    Metric = c(
      "Total cells",
      "Total cells with cloneid",
      "Total clusters",
      "Number of merged batches"
    ),
    Value = c(
      ncol(object),
      count_clones(object[[clone_col]]),
      length(unique(Idents(object))),
      length(unique(object$orig.ident))
    )
  )
  
  # -------------------- 3. Per-Sample Summary (The Easy Way) --------------------
  # We use 'initial_sample_ident' if it exists, otherwise fall back to 'orig.ident'
  group_col <- if("initial_sample_ident" %in% colnames(object@meta.data)) "initial_sample_ident" else "orig.ident"
  
  sample_summary <- object@meta.data %>%
    dplyr::group_by(Sample = .data[[group_col]]) %>%
    dplyr::summarise(
      Current_Batch = dplyr::first(orig.ident), # Shows what it was merged into
      Initial_Cells = dplyr::n(),
      Cells_with_CloneID = count_clones(.data[[clone_col]]),
      .groups = "drop"
    )
  
  # -------------------- Export --------------------
  if (export) {
    readr::write_csv(general_summary, file.path(out_outputs, paste0(name, "_summary_general.csv")))
    readr::write_csv(cluster_summary, file.path(out_outputs, paste0(name, "_summary_clusters.csv")))
    readr::write_csv(sample_summary, file.path(out_outputs, paste0(name, "_summary_samples.csv")))
    message(">>> Exported: General, Clusters, and Samples CSVs.")
  }
  
  return(invisible(list(general = general_summary, clusters = cluster_summary, samples = sample_summary)))
}