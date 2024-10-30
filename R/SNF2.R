#' Enhanced Similarity Network Fusion
#'
#' Similarity Network Fusion (SNF) integrates multiple views (similarity matrices) to
#' construct an overall status matrix. This function is adopted from the SNFtool
#' (https://github.com/maxconway/SNFtool/) and has been enhanced for unsupervised analysis of spatial ecosystem.
#' The new function supports sparse matrix and missing data in the input matrices.
#'
#' @param Wall List of similarity matrices. Each element of the list is a square,
#' symmetric matrix that shows affinities of the data points from a certain view.
#' @param K Number of neighbors in K-nearest neighbors part of the algorithm.
#' @param t Number of iterations for the diffusion process.
#' @param minibatch Integer specifying the number of columns to process in each minibatch.
#' Default is 5000. This option splits the matrix into smaller chunks (minibatch),
#' thus reducing memory usage.
#' @param ncores Integer specifying the number of CPU cores to use for parallel processing.
#' @param verbose Boolean specifying whether to show the progress messages.
#' @return A fused matrix.
#'
#' @import Matrix
#' @export
#'
SNF2 <- function(Wall, K = 50, t = 5,
                 minibatch = 5000, ncores = 1,
                 verbose = FALSE){ #
  require("Matrix")
  # Similarity Network Fusion takes multiple views of a network (Wall) and
  # fuses them together to create a overall affinity matrix.
  #
  # Args:
  #   Wall: List of matrices, each element is a square symmetric affinity
  #       matrix.
  #   K: Number of neighbors used in the K-nearest neighbours step
  #   t: Number of iterations for the diffusion process
  #
  # Returns:
  #   W: Unified similarity graph of all data types in Wall.

  check_wall_names <- function(Wall){
    # Checks if dimnames are consistent across all matrices in Wall
    #   #Move to internal functions?
    # Args:
    #   Wall: List of matrices
    # Returns:
    #   logical: True/False indicator of dimnames equivalence
    name_match <- function(names_A, names_B){
      return(identical(dimnames(names_A), dimnames(names_B)))
    }

    return(all(unlist(lapply(Wall, FUN=name_match, Wall[[1]]))))
  }

  #Check if Wall names are consistant across all matrices in Wall
  wall.name.check <- check_wall_names(Wall)
  wall.names <- dimnames(Wall[[1]])
  if(!wall.name.check){
    warning("Dim names not consistent across all matrices in Wall.
            Returned matrix will have no dim names.")
  }

  LW <- length(Wall)

  #Normalization method for affinity matrices
  normalize <- function(X){
    X = as(X, "sparseMatrix")
    row.sum.mdiag <- Matrix::rowSums(X) - diag(X)
    #If rowSumx(X) == diag(X), set row.sum.mdiag to 1 to avoid div by zero
    row.sum.mdiag[row.sum.mdiag == 0] <- 1
    X <- X/(2*(row.sum.mdiag))
    diag(X) <- 0.5
    X <- (X + Matrix::t(X))/2
    return(X)
  }

  #Normalize different networks to avoid scale problems.
  if(verbose) message(Sys.time(), " Normalize networks ...")
  Wall <- mclapply(Wall, normalize, mc.cores = ncores)

  ### Calculate the local transition matrix.
  if(verbose) message(Sys.time(), " Calculate the local transition matrix ...")
  newW <- lapply(Wall, function(X){ (.dominateset(X, K, ncores)) })
  newW <- mclapply(newW, normalize, mc.cores = ncores)

  #Perform the diffusion for t iterations
  if(verbose) message(Sys.time(), " Perform the diffusion ...")

  sum_matrix <- function(Mats){ ## Sum up a list of matrices
    sumres = Mats[[1]]
    for (i in 2:length(Mats)) {
      sumres <- sumres + Mats[[i]]
    }
    return(sumres)
  }

  for (i in 1:t){
    if(verbose) message("\t", Sys.time(), " Iteration: ", i)
    Wall <- mclapply(1:LW, function(j){
      sumWJ <- sum_matrix(Wall[-j]) / (LW-1)
      x <- matrixMultiply(newW[[j]], sumWJ, minibatch = minibatch)
      x <- matrixMultiply(x, Matrix::t(newW[[j]]), minibatch = minibatch)
      normalize(x)
    }, mc.cores = ncores)
  }

  # Construct the combined affinity matrix by summing diffused matrices
  Wall_sum <- sum_matrix(Wall)
  LW_sum <- Matrix::sparseMatrix(i = integer(0),
                                 j = integer(0),
                                 dims = dim(Wall[[1]]))

  for (i in seq_along(Wall)){
    idx <- rowSums(Wall[[i]]) == 0
    if (sum(idx) > 0){
      Wall[[i]][idx,] <- NA
      Wall[[i]][,idx] <- NA
    }
    LW_sum <- LW_sum + !is.na(Wall[[i]])
  }
  Wall <- Wall_sum / LW_sum
  Wall <- (Wall + Matrix::t(Wall)) / 2
  return(Wall)
}

.dominateset <- function(xx, K = 50, ncores = 1){
  require("Matrix")
  zero <- function(x) {
    if(K>=length(x)) return(x)
    s = sort(x, index.return=TRUE)
    x[s$ix[1:(length(x)-K)]] = 0
    return(x)
  }
  ###This function outputs the top K neighbors.
  xx <- as(xx, "sparseMatrix")
  non_zero_indices <- which(xx != 0, arr.ind = TRUE)
  non_zero_values <- xx@x
  xx@x <- ave(non_zero_values, non_zero_indices[, "col"], FUN = zero)
  # xx = drop0(xx)
  return(xx)
}

