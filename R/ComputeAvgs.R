#' Compute Cell-Type-Specific Average Expression of Spatial Clusters
#'
#' @param normdata Numeric matrix of normalized expression data, where rows
#' represent genes and columns represent cells.
#' @param scmeta Data frame containing metadata associated with each cell,
#' including spatial cluster and cell type annotations.
#' @param cluster Character string specifying the column name in 'scmeta'
#' containing spatial cluster annotations.
#' @param Region Character string specifying the column name in metadata data
#' frames containing region annotations (default: NULL).
#' @param scale A boolean specifying whether to do univariance normalization.
#' @param ncores Integer specifying the number of CPU cores to use for parallel processing.
#'
#' @return A matrix of average expression, where rows represent genes and columns
#' represent spatial clusters from the sample.
#'
#' @examples
#' library(data.table)
#' library(Seurat)
#' library(SpatialEcoTyper)
#' library(googledrive)
#' drive_deauth() # no Google sign-in is required
#' drive_download(as_id("1CoQmU3u8MoVC8RbLUvTDQmOuJJ703HHB"),
#'               "HumanMelanomaPatient1_subset_counts.tsv.gz", overwrite = TRUE)
#' drive_download(as_id("1nSPj2zRywFUdbo1fwiz77ds4NuM6bmV2"),
#'               "Melanoma1_subset_SpatialEcoTyper_results.rds", overwrite = TRUE)
#' scdata <- fread("HumanMelanomaPatient1_subset_counts.tsv.gz",
#'                 sep = "\t",header = TRUE, data.table = FALSE)
#' rownames(scdata) <- scdata[, 1]
#' scdata <- as.matrix(scdata[, -1])
#' tmpobj <- CreateSeuratObject(scdata) %>%
#'         SCTransform(clip.range = c(-10, 10), verbose = FALSE)
#' seurat_version = as.integer(gsub("\\..*", "", as.character(packageVersion("SeuratObject"))))
#' if(seurat_version<5){
#'   normdata <- GetAssayData(tmpobj, "data")
#' }else{
#'   normdata <- tmpobj[["SCT"]]$data
#' }
#' metadata = readRDS("Melanoma1_subset_SpatialEcoTyper_results.rds")$metadata
#'
#' # Construct cell-type-specific gene expression signatures of SEs
#' avgexprs <- ComputeAvgs(normdata = normdata, scmeta = metadata)
#' head(avgexprs)
#'
#' @import Seurat
#' @export
#'
ComputeAvgs <- function(normdata, scmeta, cluster = "SE",
                       Region = NULL, scale = TRUE, ncores = 4){
  scmeta$SE = scmeta[, cluster]
  scmeta = scmeta[!is.na(scmeta$SE), ]
  cids = na.omit(intersect(rownames(scmeta), colnames(normdata)))
  scmeta = scmeta[cids, ]
  normdata = normdata[, cids]
  ## Seurat object
  obj = CreateSeuratObject(normdata, meta.data = scmeta)
  seurat_version = as.integer(gsub("\\..*", "", as.character(packageVersion("SeuratObject"))))
  if(seurat_version>=5) obj[["RNA"]]$data = obj[["RNA"]]$counts
  rm(normdata)
  ##### Within sample normalization ######
  lfcs <- mclapply(unique(obj$CellType), function(ct){
    if(sum(obj$CellType==ct)<10) return(NULL)
    tmpobj <- subset(obj, CellType==ct)
    slot = "data"
    if(scale){
      if(!is.null(Region) && (Region%in%colnames(scmeta))){
        if(seurat_version>=5){
          tmpdat <- Znorm(tmpobj[["RNA"]]$data, groups = tmpobj@meta.data[, Region])
        }else{
          tmpdat <- Znorm(GetAssayData(tmpobj, "data"), groups = tmpobj@meta.data[, Region])
        }
        tmpobj <- CreateSeuratObject(tmpdat, meta.data = tmpobj@meta.data)
        if(seurat_version>=5) tmpobj[["RNA"]]$data = tmpobj[["RNA"]]$counts
        slot <- "data"
      }else{
        tmpobj <- ScaleData(tmpobj, verbose = FALSE)
        slot <- "scale.data"
      }
    }
    Idents(tmpobj) <- tmpobj$SE
    if(length(unique(tmpobj$SE))<2) return(NULL)
    lfcs <- AverageExpression(tmpobj, slot = slot)$RNA
    rownames(lfcs) <- paste0(ct, "..", rownames(lfcs))
    lfcs
  }, mc.cores = ncores)
  names(lfcs) <- unique(obj$CellType)
  lfcs <- mclapply(lfcs, function(x){
    x = as.matrix(x)
    x <- x[, match(unique(obj$SE), colnames(x))]
    colnames(x) <- unique(obj$SE)
    x
  }, mc.cores = ncores)
  lfcs <- do.call(rbind, lfcs)
  return(lfcs)
}


#' #' Compute Cell-Type-Specific Fold Changes (FCs) for Spatial Clusters
#' #'
#' #' This function computes fold changes (FCs) for spatial transcriptomic data,
#' #' comparing expression levels between spatial clusters within each cell type.
#' #'
#' #' @param normdata Numeric matrix of normalized expression data, where rows
#' #' represent genes and columns represent cells.
#' #' @param scmeta Data frame containing metadata associated with each cell,
#' #' including spatial cluster and cell type annotations.
#' #' @param cluster Character string specifying the column name in 'scmeta'
#' #' containing spatial cluster annotations.
#' #' @param Region Character string specifying the column name in metadata data
#' #' frames containing region annotations (default: NULL).
#' #' @param scale A boolean specifying whether to do univariance normalization.
#' #' @param ncores Integer specifying the number of CPU cores to use for parallel processing.
#' #'
#' #' @return A matrix of fold changes (FCs), where rows represent genes and columns
#' #' represent spatial clusters from the sample.
#' #'
#' #'
#' #' @import Seurat
#' #' @export
#' #'
#' ComputeFCs <- function(normdata, scmeta, cluster = "SE",
#'                        Region = NULL, scale = FALSE,
#'                        ncores = parallel::detectCores()){
#'   if(!"CellType" %in% colnames(scmeta)){
#'     stop("the meta data have to include a column (CellType) for cell type annotations")
#'   }
#'   scmeta$SE = scmeta[, cluster]
#'   scmeta <- scmeta[!is.na(scmeta$SE), ]
#'   normdata <- normdata[, match(rownames(scmeta), colnames(normdata))]
#'
#'   ## Seurat object
#'   obj <- CreateSeuratObject(normdata, meta.data = scmeta)
#'   seurat_version = as.integer(gsub("\\..*", "", as.character(packageVersion("SeuratObject"))))
#'   if(seurat_version>=5) obj[["RNA"]]$data = obj[["RNA"]]$counts
#'   rm(normdata)
#'   ##### Within sample normalization ######
#'   lfcs <- mclapply(unique(obj$CellType), function(ct){
#'     if(sum(obj$CellType==ct)<10) return(NULL)
#'     tmpobj = subset(obj, CellType==ct)
#'     slot = "data"
#'     if(scale){
#'       if(!is.null(Region) && (Region%in%colnames(scmeta))){
#'         if(seurat_version>=5){
#'           tmpdat <- Znorm(tmpobj[["RNA"]]$data, groups = tmpobj@meta.data[, Region])
#'         }else{
#'           tmpdat <- Znorm(GetAssayData(tmpobj, "data"), groups = tmpobj@meta.data[, Region])
#'         }
#'         tmpobj <- CreateSeuratObject(tmpdat, meta.data = tmpobj@meta.data)
#'       }else{
#'         tmpobj <- ScaleData(tmpobj)
#'         slot <- "scale.data"
#'       }
#'       if(seurat_version>=5) tmpobj[["RNA"]]$data = tmpobj[["RNA"]]$counts
#'     }
#'     Idents(tmpobj) <- tmpobj$SE
#'     if(length(unique(tmpobj$SE))<2) return(NULL)
#'     lfcs <- mclapply(unique(tmpobj$SE), function(cl){
#'       lfcs <- FoldChange(tmpobj, ident.1 = cl, slot = slot)
#'       lfcs[,1]
#'     }, mc.cores = ncores)
#'     lfcs <- do.call(cbind, lfcs)
#'     colnames(lfcs) <- unique(tmpobj$SE)
#'     rownames(lfcs) <- paste0(ct, "..", rownames(tmpobj))
#'     lfcs <- lfcs[, match(unique(obj$SE), colnames(lfcs))]
#'     colnames(lfcs) <- unique(obj$SE)
#'     lfcs[is.na(lfcs)] <- 0
#'     lfcs
#'   }, mc.cores = ncores)
#'   lfcs <- do.call(rbind, lfcs)
#'   return(lfcs)
#' }
