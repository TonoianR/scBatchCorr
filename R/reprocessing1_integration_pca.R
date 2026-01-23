#' SCT-based RPCA integration followed by PCA (Seurat v5 Robust Version)
#'
#' Performs SCT normalization, RPCA-based integration, PCA, and ElbowPlot
#' with aggressive cleaning to prevent Seurat v5 layer-mismatch errors.
#'
#' @param obj A Seurat object.
#' @param name Character scalar used as filename prefix.
#' @param nfeatures Number of integration features.
#' @param npcs Number of principal components.
#' @param out_dir Base output directory.
#' @param seed Random seed for reproducibility.
#' @param min_cells Minimum cells required per sample (default: 10).
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
    seed = NULL,
    min_cells = 30 # Standard statistical floor
) {
  if (!is.null(seed)) set.seed(seed)
  out_outputs <- file.path(out_dir, "Outputs"); out_plots <- file.path(out_dir, "Plots")
  dir.create(out_outputs, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
  
  # 1. Clean Assays
  DefaultAssay(obj) <- "RNA"
  for (assay_name in names(obj@assays)) {
    if (assay_name != "RNA") obj[[assay_name]] <- NULL
  }
  obj <- JoinLayers(obj)
  
  # 2. SMART MERGING LOGIC
  meta <- obj@meta.data
  counts <- table(meta$orig.ident)
  small_batches <- names(counts)[counts < min_cells]
  large_batches <- names(counts)[counts >= min_cells]
  
  if(length(small_batches) > 0) {
    message(">>> Detected small batches: ", paste(small_batches, collapse=", "))
    
    for(sb in small_batches) {
      # Extract prefix (e.g., "E16liver" from "E16liverAC21")
      # We assume the prefix is the alphabetic part before the numbers
      prefix <- gsub("[0-9]+$", "", sb)
      
      # Find a large batch with the same prefix
      target <- grep(prefix, large_batches, value = TRUE)[1]
      
      if(!is.na(target)) {
        new_name <- paste0(target, "+merged_", sb)
        # Update metadata for all cells originally in the target OR the small batch
        meta$orig.ident[meta$orig.ident == sb | meta$orig.ident == target] <- new_name
        # Update the large_batches list to reflect the new name for subsequent small batches
        large_batches[large_batches == target] <- new_name
        message(">>> Merged ", sb, " into ", target, " -> New Batch: ", new_name)
      } else {
        message(">>> Warning: No matching prefix found for ", sb, ". It will be dropped.")
      }
    }
    obj@meta.data <- meta
  }
  
  # 3. Split and Filter (using the updated metadata)
  obj_list <- SplitObject(obj, split.by = "orig.ident")
  cell_counts <- sapply(obj_list, ncol)
  obj_list <- obj_list[names(cell_counts)[cell_counts >= min_cells]]
  
  if (length(obj_list) < 2) stop("Not enough samples after merging/filtering.")
  
  # 4. Standard PCA & Integration (Now safe at npcs = 50)
  obj_list <- lapply(obj_list, function(x) {
    x <- SCTransform(x, vst.flavor = "v2", method = "glmGamPoi", verbose = FALSE)
    x <- RunPCA(x, npcs = npcs, verbose = FALSE) # Back to full npcs!
    return(x)
  })
  
  features <- SelectIntegrationFeatures(obj_list, nfeatures = nfeatures)
  obj_list <- PrepSCTIntegration(obj_list, anchor.features = features)
  
  anchors <- FindIntegrationAnchors(
    object.list = obj_list, normalization.method = "SCT",
    anchor.features = features, dims = 1:npcs, reduction = "rpca"
  )
  
  integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  
  # 5. Finalize Structure
  integrated <- JoinLayers(integrated, assay = "RNA")
  integrated <- NormalizeData(integrated, assay = "RNA", verbose = FALSE)
  integrated <- ScaleData(integrated, assay = "RNA", verbose = FALSE)
  
  DefaultAssay(integrated) <- "integrated"
  integrated <- RunPCA(integrated, npcs = npcs, verbose = FALSE)
  
  qs::qsave(integrated, file.path(out_outputs, paste0(name, "_integrated_after_pca.qs")))
  return(integrated)
}