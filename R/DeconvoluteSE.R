#' Infer SE Abundances Using Pretrained NMF Model
#'
#' This function predicts SE abundances in each mixture.
#'
#' @param dat A numeric matrix of bulk gene expression data.
#' @param scale Logical indicating whether to scale the input data for predictions.
#' @param W A W matrix for inferring SE abundances from bulk expression data.
#'
#' @return A matrix of SE abundances in bulk tumors.
#'
#' @examples
#' library(googledrive)
#' drive_deauth() # no Google sign-in is required
#' drive_download(as_id("14QvmgISxaArTzWt_UHvf55aAYN2zm84Q"), "SKCM_RNASeqV2.geneExp.rds",
#'                     overwrite = TRUE)
#' bulkdata <- readRDS("SKCM_RNASeqV2.geneExp.rds")
#'
#' # Predict SE abundances in bulk tumors
#' se_abundances <- DeconvoluteSE(dat = bulkdata)
#' head(se_abundances[, 1:5])
#'
#' @export
#'
#' @importFrom parallel mclapply
#'
DeconvoluteSE <- function(dat, scale = TRUE, W = NULL){
  if(all(dat>=0) & max(dat)>80){
    dat = log2(dat+1)
  }
  if(is.null(W)){
    W <- readRDS(file.path(system.file("extdata", package = "SpatialEcoTyper"),
                           "Bulk_SE_Recovery_W.rds"))
  }
  preds <- NMFpredict(W = W, dat, scale = scale)
  return(t(preds))
}
