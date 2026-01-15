#' Rebuild Seurat object from raw RNA counts (clean state)
#'
#' Reconstructs a clean Seurat object for a selected cell set using
#' raw RNA counts from a parent object. All derived data (SCT, data,
#' scale.data, reductions, graphs) are removed.
#'
#' This function is intended as Step 0 prior to SCT-based reprocessing.
#'
#' @param subset_obj A Seurat object representing a subset of cells.
#' @param parent_obj The original Seurat object containing raw RNA counts.
#'
#' @return A clean Seurat object with a single RNA assay containing only counts.
#'
#' @export
reprocessing0_rebuild_from_raw <- function(
    subset_obj,
    parent_obj
) {
  
  # -------------------- validation --------------------
  if (!inherits(subset_obj, "Seurat")) {
    stop("`subset_obj` must be a Seurat object.")
  }
  
  if (!inherits(parent_obj, "Seurat")) {
    stop("`parent_obj` must be a Seurat object.")
  }
  
  if (!"RNA" %in% names(parent_obj@assays)) {
    stop("`parent_obj` must contain an RNA assay.")
  }
  
  # -------------------- cell intersection --------------------
  cells <- intersect(Cells(subset_obj), Cells(parent_obj))
  
  if (length(cells) == 0) {
    stop("No overlapping cells between subset and parent object.")
  }
  
  # -------------------- rebuild object --------------------
  counts <- GetAssayData(
    parent_obj,
    assay = "RNA",
    layer = "counts"
  )[, cells, drop = FALSE]
  
  meta <- parent_obj@meta.data[cells, , drop = FALSE]
  
  clean_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = meta
  )
  
  # -------------------- sanity checks --------------------
  DefaultAssay(clean_obj) <- "RNA"
  
  if (!identical(Layers(clean_obj[["RNA"]]), "counts")) {
    stop("RNA assay contains layers other than counts.")
  }
  
  if (length(clean_obj@reductions) != 0) {
    stop("Reductions detected in clean object (should be empty).")
  }
  
  if (length(clean_obj@graphs) != 0) {
    stop("Graphs detected in clean object (should be empty).")
  }
  
  VariableFeatures(clean_obj) <- character(0)
  
  return(clean_obj)
}
