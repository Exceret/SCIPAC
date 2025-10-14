
#' Normalize scRNA-seq data and find highly variable genes without scaling using Seurat
#'
#' @param exprs.data single-cell data matrix. Rows stand for genes and columns stand for cells.
#' @param hvg the number of highly variable genes. The default is 1000.
#'
#' @return a single-cell data matrix. Rows stand for genes and columns stand for cells.
#' @importFrom dplyr %>%
obtain.preprocessed.data <- function(exprs.data, hvg = 1000){
  exprs.data <- Seurat::CreateSeuratObject(exprs.data, project = "CreateSeuratObject", min.cells = 3, min.features = 200)
  exprs.data <- Seurat::NormalizeData(exprs.data, normalization.method = "LogNormalize", scale.factor = 10000)
  exprs.data <- Seurat::FindVariableFeatures(exprs.data, selection.method = "vst", nfeatures = hvg)

  exprs.data.norm <- exprs.data[["RNA"]]$data
  exprs.data.variable.genes <- Seurat::VariableFeatures(exprs.data)

  processed.data <- exprs.data.norm[exprs.data.variable.genes, ] %>%
    as.matrix()

  return(processed.data)
}



#' Pre-process both single-cell and bulk data
#'
#' This function first finds overlap genes between scRNA-seq and bulk data,
#' then normalizes scRNA-seq data and find variable genes,
#' log-transforms bulk RNA-seq data,
#' and uses high variable genes for both scRNA-seq and bulk data
#'
#' @param sc.dat single-cell data matrix. Rows stand for genes and columns stand for cells.
#' @param bulk.dat bulk data matrix. Rows stand for genes and columns stand for cells.
#' @param hvg the number of highly variable genes. The default is \code{1000}.
#'
#' @return A list object contains pre-processed single-cell and bulk data
#' @importFrom dplyr %>%
#' @export

preprocess.sc.bulk.dat <- function(sc.dat, bulk.dat, hvg = 1000){
  overlap.genes <- intersect(rownames(sc.dat), rownames(bulk.dat))

  sc.dat.new <- sc.dat[overlap.genes, ]

  bulk.dat.new <- bulk.dat[overlap.genes, ] %>% as.matrix()
  bulk.dat.new <- log(bulk.dat.new + 1)


  sc.dat.preprocessed <- obtain.preprocessed.data(exprs.data = sc.dat.new, hvg = hvg)
  bulk.dat.preprocessed <- bulk.dat.new[rownames(sc.dat.preprocessed), ]

  return(list("sc.dat.preprocessed" = sc.dat.preprocessed,
              "bulk.dat.preprocessed" = bulk.dat.preprocessed))
}


#' Normalize spatial data and find highly variable genes using Seurat
#'
#' @param spatial.obj Seurat spatial object (e.g., from Load10X_Spatial)
#' @param hvg the number of highly variable genes. The default is 1000.
#' @param assay the assay to use. If NULL, will auto-detect from available assays.
#' Common values: "Spatial" (Visium), "SCT" (spatial SCTransform), "RNA" (regular spatial)
#'
#' @return a spatial expression matrix. Rows stand for genes and columns stand for spots/cells.
#' @importFrom dplyr %>%
obtain.preprocessed.spatial.data <- function(spatial.obj, hvg = 1000, assay = NULL){

  # Auto-detect assay if not specified
  if(is.null(assay)){
    available.assays <- Seurat::Assays(spatial.obj)
    if("Spatial" %in% available.assays){
      assay <- "Spatial"
    } else if("SCT" %in% available.assays){
      assay <- "SCT"
    } else if("RNA" %in% available.assays){
      assay <- "RNA"
    } else {
      stop("Could not auto-detect appropriate assay. Please specify assay parameter.")
    }
    message(paste0("Using assay: ", assay))
  }

  # Set default assay
  Seurat::DefaultAssay(spatial.obj) <- assay

  # Normalize and find variable features
  # For Spatial assay, use standard normalization
  if(assay %in% c("Spatial", "RNA")){
    spatial.obj <- Seurat::NormalizeData(spatial.obj, normalization.method = "LogNormalize", scale.factor = 10000)
    spatial.obj <- Seurat::FindVariableFeatures(spatial.obj, selection.method = "vst", nfeatures = hvg)

    # Extract normalized data for variable genes
    spatial.data.norm <- Seurat::GetAssayData(spatial.obj, assay = assay, layer = "data")
    variable.genes <- Seurat::VariableFeatures(spatial.obj)

  } else if(assay == "SCT"){
    # For SCT assay, normalization is already done
    variable.genes <- Seurat::VariableFeatures(spatial.obj)
    if(length(variable.genes) < hvg){
      warning(paste0("SCT assay only has ", length(variable.genes),
                     " variable features. Using all available."))
    }
    spatial.data.norm <- Seurat::GetAssayData(spatial.obj, assay = "SCT", layer = "data")
  }

  processed.data <- spatial.data.norm[variable.genes, ] %>%
    as.matrix()

  return(processed.data)
}


#' Pre-process both spatial transcriptomics and bulk data
#'
#' This function first finds overlap genes between spatial and bulk data,
#' then normalizes spatial data and finds variable genes,
#' log-transforms bulk RNA-seq data,
#' and uses high variable genes for both spatial and bulk data
#'
#' @param spatial.obj Seurat spatial object (e.g., from Load10X_Spatial)
#' @param bulk.dat bulk data matrix. Rows stand for genes and columns stand for samples.
#' @param hvg the number of highly variable genes. The default is \code{1000}.
#' @param assay the assay to use from spatial object. If NULL, will auto-detect.
#'
#' @return A list object contains pre-processed spatial and bulk data
#' @importFrom dplyr %>%
#' @export
#'
preprocess.spatial.bulk.dat <- function(spatial.obj, bulk.dat, hvg = 1000, assay = NULL){

  # Extract gene names from spatial object
  if(is.null(assay)){
    available.assays <- Seurat::Assays(spatial.obj)
    if("Spatial" %in% available.assays){
      assay <- "Spatial"
    } else if("SCT" %in% available.assays){
      assay <- "SCT"
    } else if("RNA" %in% available.assays){
      assay <- "RNA"
    }
  }

  spatial.genes <- rownames(spatial.obj[[assay]])
  overlap.genes <- intersect(spatial.genes, rownames(bulk.dat))

  message(paste0("Found ", length(overlap.genes), " overlapping genes between spatial and bulk data"))

  # Subset spatial object to overlap genes
  spatial.obj.sub <- spatial.obj[overlap.genes, ]

  # Subset and log-transform bulk data
  bulk.dat.new <- bulk.dat[overlap.genes, ] %>% as.matrix()
  bulk.dat.new <- log(bulk.dat.new + 1)

  # Preprocess spatial data
  spatial.dat.preprocessed <- obtain.preprocessed.spatial.data(
    spatial.obj = spatial.obj.sub,
    hvg = hvg,
    assay = assay
  )

  # Use same genes for bulk data
  bulk.dat.preprocessed <- bulk.dat.new[rownames(spatial.dat.preprocessed), ]

  return(list("spatial.dat.preprocessed" = spatial.dat.preprocessed,
              "bulk.dat.preprocessed" = bulk.dat.preprocessed,
              "spatial.obj" = spatial.obj))  # Return original object for downstream use
}
