#' Seurat v5 Hybrid Integration (SCT + RPCA)
#'
#' Performs SCTransform normalization (v2), batch-specific PCA, and 
#' RPCA-based integration using the Seurat v5 IntegrateLayers engine. 
#' Includes smart merging of small batches (< min_cells) into larger ones 
#' based on sample prefixes to ensure statistical stability.
#'
#' @param obj A Seurat object.
#' @param name Character scalar used as filename prefix for saving.
#' @param nfeatures Number of integration features.
#' @param npcs Number of principal components to use (default: 50).
#' @param out_dir Base output directory for the .qs file.
#' @param seed Random seed for reproducibility.
#' @param min_cells Minimum cells required per batch to avoid RPCA errors.
#'
#' @return A Seurat object with an integrated 'pca' reduction and 
#'         consolidated RNA layers (counts, data, scale.data).
#'
#' @export
reprocessing1_integration_pca <- function(
    obj,
    name,
    nfeatures = 3000,
    npcs = 50,
    out_dir,
    seed = 1234,
    min_cells = 30
) {
  set.seed(seed)
  out_outputs <- file.path(out_dir, "Outputs")
  dir.create(out_outputs, recursive = TRUE, showWarnings = FALSE)
  
  # 1. Clean & Reset
  message(">>> Step 1: Cleaning assays...")
  DefaultAssay(obj) <- "RNA"
  obj[["SCT"]] <- NULL
  obj[["integrated"]] <- NULL
  obj <- JoinLayers(obj)
  
  # 2. Smart Merging
  message(">>> Step 2: Running Smart Merging logic...")
  meta <- obj@meta.data
  counts <- table(meta$orig.ident)
  small_batches <- names(counts)[counts < min_cells]
  large_batches <- names(counts)[counts >= min_cells]
  
  if(length(small_batches) > 0) {
    for(sb in small_batches) {
      prefix <- gsub("[0-9]+$", "", sb)
      target <- grep(prefix, large_batches, value = TRUE)[1]
      if(!is.na(target)) {
        new_name <- paste0(target, "+merged_", sb)
        meta$orig.ident[meta$orig.ident == sb | meta$orig.ident == target] <- new_name
        large_batches[large_batches == target] <- new_name
        message("    [Merge] ", sb, " -> ", target)
      }
    }
    obj@meta.data <- meta
  }
  
  # 3. SCT Multi-layer Workflow
  message(">>> Step 3: Splitting Layers and SCT Normalization...")
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
  obj <- SCTransform(obj, vst.flavor = "v2", method = "glmGamPoi", verbose = FALSE)
  
  # 4. PCA per batch
  message(">>> Step 4: Running PCA per batch...")
  obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  
  # 5. RPCA Integration (Seurat v5 style)
  message(">>> Step 5: Integrating Layers using RPCA math...")
  obj <- IntegrateLayers(
    object = obj, 
    method = RPCAIntegration, 
    orig.reduction = "pca", 
    new.reduction = "integrated_pca",
    normalization.method = "SCT",
    verbose = TRUE
  )
  
  # 6. Finalize Reductions and RNA Layers
  message(">>> Step 6: Consolidating RNA for visualization...")
  
  # Set the integrated reduction
  obj[["pca"]] <- obj[["integrated_pca"]]
  obj[["integrated_pca"]] <- NULL
  
  # CRITICAL FIX: Switch to RNA before JoinLayers
  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj)
  
  # Standardize RNA
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
  # Scale only variable features for speed/memory
  obj <- ScaleData(obj, assay = "RNA", features = VariableFeatures(obj), verbose = FALSE)
  
  # Re-set Default to SCT (best for finding markers later)
  DefaultAssay(obj) <- "SCT"
  
  # ElbowPlot needs this to know the height of the bars
  obj[["pca"]]@stdev <- apply(Embeddings(obj, "pca"), 2, sd)
  
  message(">>> Step 7: Generating ElbowPlot...")
  p_elbow <- ElbowPlot(obj, ndims = npcs) + 
    ggplot2::ggtitle(paste("PCA Variance -", name))
  
  out_plots <- file.path(out_dir, "Plots")
  dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
  
  ggplot2::ggsave(
    filename = file.path(out_plots, paste0(name, "_PCA_ElbowPlot.png")),
    plot = p_elbow, width = 8, height = 6, dpi = 300
  )
  
  message(">>> Step 8: Saving...")
  qs::qsave(obj, file.path(out_outputs, paste0(name, "_integrated.qs")))
  
  # Final Audit Report (Fixed object names)
  message("--- INTEGRATION COMPLETE ---")
  message("Object Name:    ", name)
  message("Final Assays:   ", paste(Assays(obj), collapse = ", "))
  message("RNA Layers:     ", paste(Layers(obj[["RNA"]]), collapse = ", "))
  message("Default Assay:  ", DefaultAssay(obj))
  message("Total Cells:    ", ncol(obj))
  message("-----------------------------")
  
  return(obj)
}