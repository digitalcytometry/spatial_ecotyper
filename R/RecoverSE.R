#' Recovery of SEs Using Pretrained NMF Models
#'
#' This function can recover SEs from Visium, scRNA-seq data or single-cell spatial data.
#'
#' @param dat A numeric matrix of gene expression data, from single-cell spatial
#' transcriptomics, scRNA-seq or spot-resolution spatial transcriptomics data.
#' @param scale Logical indicating whether to scale the input data for predictions.
#' @param Ws A list of cell-type-specific W matrices used to recover SE-specific
#' cell states. Each element in the list should be named after the corresponding
#' cell type and contain a W matrix from an NMF model.
#' @param celltypes Character string specifying the cell type annotations, which
#' is required for scRNA-seq or single-cell spatial data.
#' @param se_results A list including a seurat object and a metadata with spatial
#' cluster annotations (SE column) returned by SpatialEcoTyper. When supplied, the `dat`
#' should be single cell gene expression data used for the SpatialEcoTyper analysis.
#'
#' @return Depending on the input data:
#' - For single-cell data: A vector of predicted SEs for each cell.
#' - For spot-resolution data (e.g., Visium): A matrix of SE abundances across spots.
#'
#' @examples
#' # see https://digitalcytometry.github.io/SpatialEcoTyper_dev/articles/Recovery_scRNA.html
#' # see https://digitalcytometry.github.io/SpatialEcoTyper_dev/articles/Recovery_Spatial.html
#'
#' @export
#'
#' @importFrom parallel mclapply
#'

RecoverSE <- function(dat, scale = TRUE,
                      Ws = NULL,
                      celltypes = NULL,
                      se_results = NULL){
  flag = ifelse(is.null(Ws), "default", "custom")
  ## Load NMF models
  if(is.null(Ws)){
    Wfiles <- list.files(system.file("extdata", package = "SpatialEcoTyper"),
                         "MERSCOPE_W_v1.*.rds", full.names = TRUE)
    Ws <- lapply(Wfiles, readRDS)
    names(Ws) <- gsub(".*v1_|.rds", "", Wfiles)
  }

  if(!is.null(se_results)){ ## SE recovery for single cell spatial data
    scmeta = se_results$metadata
    scmeta = scmeta[match(colnames(dat), rownames(scmeta)), ]
    scmeta$State = paste0(scmeta$SE, "_", scmeta$CellType)

    cts <- unique(intersect(names(Ws), unique(scmeta$CellType)))
    resDF <- mclapply(cts, function(x){
      idx = !(scmeta$CellType!=x | is.na(scmeta$SE))
      if(sum(idx)<2) return(NULL)
      tmpdat = dat[, idx]
      tmpmeta = scmeta[idx, ]
      if(scale) tmpdat = Seurat::ScaleData(tmpdat, verbose = FALSE)
      tmpdat = as.matrix(tmpdat)
      cell2cluster = matrix(0, nrow = ncol(tmpdat), ncol = length(unique(tmpmeta$SE)),
                            dimnames = list(colnames(tmpdat), unique(tmpmeta$SE)))
      idx = cbind(match(colnames(tmpdat), rownames(cell2cluster)),
                  match(tmpmeta$SE, colnames(cell2cluster)))
      idx = as.array(idx)
      cell2cluster[idx] = 1
      cell2cluster = t(t(cell2cluster) / colSums(cell2cluster))
      sedat = tmpdat %*% cell2cluster

      H = NMFpredict(Ws[[x]], sedat, ncell.per.run = 500, scale = FALSE)
      preds = rownames(H)[apply(H, 2, which.max)]
      names(preds) = colnames(H)
      resDF = data.frame(Cluster = names(preds), SE = preds,
                          CellType = x, PredScore = apply(H, 2, max))
      return(resDF)
    })
    resDF <- do.call(rbind, resDF)
    resDF$State <- paste0(resDF$SE, "_", resDF$CellType)
    resDF$SE[resDF$PredScore<0.6] = "nonSE"
    if(flag=="default"){
      states <- readRDS(file.path(system.file("extdata", package = "SpatialEcoTyper"),
                                  "SE_CellStates.rds"))
      states <- setdiff(states, "SE01_B")
      resDF$State <- paste0(resDF$SE, "_", resDF$CellType)
      resDF$SE[!resDF$State %in% states] <- "nonSE"
      semap = c(SE01 = "SE1", SE02 = "SE2", SE03 = "SE3", SE04 = "SE4", SE05 = "SE5",
                SE06 = "nonSE", SE07 = "SE6", SE08 = "SE7", SE09 = "nonSE",
                SE10 = "SE8", SE11 = "SE9", nonSE = "nonSE")
      resDF$SE = semap[resDF$SE]
      resDF$State <- paste0(resDF$SE, "_", resDF$CellType)
    }
    preds <- resDF %>% count(Cluster, SE) %>%
      group_by(Cluster) %>% mutate(Frac = n / sum(n))
    wts <- 1 / table(preds$SE[preds$Cluster %in% scmeta$SE])
    preds$Frac <- preds$Frac * wts[preds$SE]
    preds <- preds %>% arrange(Cluster, -Frac) %>% distinct(Cluster, .keep_all = TRUE)
    idx <- match(scmeta$SE, preds$Cluster)
    preds <- preds$SE[idx]
    preds[is.na(preds)] = "nonSE"
    names(preds) <- rownames(scmeta)
    return(preds)
  }else if(!is.null(celltypes)){ ## SE recovery for scRNA-seq data
    if(ncol(dat)!=length(celltypes)){
      stop("The cell type annotations do not match the columns of expression data")
    }
    if(!is.null(names(celltypes))) celltypes <- celltypes[colnames(dat)]
    cts <- unique(intersect(names(Ws), unique(celltypes)))
    resDF <- mclapply(cts, function(x){
      message(Sys.time(), " Predict ", x, " cell states ")
      if(sum(celltypes==x)>20){
        tmpdat <- dat[, celltypes==x]
        if(scale) tmpdat <- Seurat::ScaleData(tmpdat, verbose = FALSE)
        tmpdat <- as.matrix(tmpdat)
        H <- NMFpredict(Ws[[x]], tmpdat, ncell.per.run = 500, scale = FALSE)
        preds <- rownames(H)[apply(H, 2, which.max)]
        names(preds) <- colnames(H)
        PredScore <- apply(H, 2, max)
        resDF <- data.frame(CID = names(preds), CellType = x, SE = preds, PredScore = PredScore)
        return(resDF)
      }else return(NULL)
    })
    resDF <- do.call(rbind, resDF)
    idx <- !(colnames(dat) %in% resDF$CID)
    if(sum(idx)>1){
      tmpDF <- data.frame(CID = colnames(dat)[idx], CellType = celltypes[idx],
                          SE = "nonSE", PredScore = 1)
      resDF <- rbind(resDF, tmpDF)
    }
    resDF$SE[resDF$PredScore<0.6] <- "nonSE"
    if(flag=="default"){
      states <- readRDS(file.path(system.file("extdata", package = "SpatialEcoTyper"),
                                  "SE_CellStates.rds"))
      states <- setdiff(states, "SE01_B")
      resDF$State <- paste0(resDF$SE, "_", resDF$CellType)
      resDF$SE[!resDF$State %in% states] <- "nonSE"
      semap = c(SE01 = "SE1", SE02 = "SE2", SE03 = "SE3", SE04 = "SE4", SE05 = "SE5",
                SE06 = "nonSE", SE07 = "SE6", SE08 = "SE7", SE09 = "nonSE",
                SE10 = "SE8", SE11 = "SE9", nonSE = "nonSE")
      resDF$SE = semap[resDF$SE]
    }
    resDF <- resDF[match(colnames(dat), resDF$CID), ]
    preds <- resDF$SE
    names(preds) <- colnames(dat)
    return(preds)
  }else{ ## SE recovery for Visium data
    if(scale) dat <- Seurat::ScaleData(dat, verbose = FALSE)
    statepreds <- mclapply(names(Ws), function(x){
      tmp = NMFpredict(W = Ws[[x]], dat, scale = FALSE)
      rownames(tmp) <- paste0(rownames(tmp), "_", x)
      tmp
    })
    statepreds <- do.call(rbind, statepreds)
    ses <- gsub("_.*", "", rownames(statepreds))
    state2se <- matrix(0, nrow = nrow(statepreds), ncol = length(unique(ses)),
                       dimnames = list(rownames(statepreds), unique(ses)))
    idx <- cbind(match(rownames(statepreds), rownames(state2se)),
                 match(ses, colnames(state2se)))
    state2se[as.array(idx)] <- 1
    state2se <- t(t(state2se) / colSums(state2se))
    preds <- t(statepreds) %*% state2se
    preds <- preds / rowSums(preds)

    if(flag=="default"){
      semap = c(SE01 = "SE1", SE02 = "SE2", SE03 = "SE3", SE04 = "SE4", SE05 = "SE5",
                SE06 = "nonSE", SE07 = "SE6", SE08 = "SE7", SE09 = "nonSE",
                SE10 = "SE8", SE11 = "SE9", nonSE = "nonSE")
      nonse = rowSums(preds[, colnames(preds)%in%c("SE06", "SE09", "nonSE")])
      tmp = preds[, !colnames(preds)%in%c("SE06", "SE09", "nonSE")]
      colnames(tmp) = semap[colnames(tmp)]
      preds = cbind(tmp, nonSE = nonse)
    }
    return(preds)
  }
}
