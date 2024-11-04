#' Construct a pairwise similarity matrix
#'
#' `similarity_matrix` is a helper function used within the `rips_similarity_matrix`
#' function to construct a pairwise similarity matrix comparing pairs of
#' persistence diagrams.
#'
#' @param rips.list A list of persistence diagrams for each sample
#' @param n.sample The number of samples in the dataset (the length of `rips.list`)
#' @param PIDs Unique identifiers for each sample. This will be used to name the
#' resulting similarity matrix.
#' @param print.progress Should a progress bar be printed?
#'
#' @return Returns a list of kernel (similarity matrices) for each homology group.
#'
#' @import RcppAlgos
#' @import dplyr
#' @import MiRKAT
#' @importFrom stats complete.cases
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' require(magrittr)
#' require(TDAstats)
#'
#' # Save the PIDs
#' PIDs <- unique(data1.df$PID)
#'
#' # Save the number of samples
#' n.sample <- length(PIDs)
#'
#' # Initialize a list to store the Rips filtrations
#' rips.list <- lapply(1:n.sample, function(i) list()); names(rips.list) <- PIDs
#'
#' # Iterate through the samples to construct Rips filtration
#' for (i in PIDs) {
#'
#'
#'  # Subset the data to just this PID
#'  data.i <- data1.df %>%
#'    dplyr::filter(PID == i) %>%
#'    dplyr::select(x,y)
#'
#'  # Construct a Rips filtration using TDAstats
#'  rips.i <- TDAstats::calculate_homology(data.i, dim = 1, threshold = 10)
#'
#'  # Save
#'  rips.list[[i]] <- rips.i
#' }
#'  K.list <- similarity_matrix(rips.list, n.sample, PIDs)
similarity_matrix <- function(rips.list, n.sample, PIDs, print.progress = FALSE) {

  # Create all pairs of samples
  pairs <- RcppAlgos::comboGeneral(PIDs, 2)

  # Initialize a distance matrix
  dist.deg0 <- dist.deg1 <- matrix(0, nrow = n.sample, ncol = n.sample)
  rownames(dist.deg0) <- colnames(dist.deg0) <- rownames(dist.deg1) <- colnames(dist.deg1) <- PIDs

  # Iterate through pairs and calculate distance
  for (i in 1:nrow(pairs)) {

    # Print progress
    if (print.progress) print(paste0(i, "/", nrow(pairs)))

    # Save the current pair
    ids <- pairs[i,]
    id.1 <- ids[1]
    id.2 <- ids[2]

    rips.1 <- rips.list[[id.1]]
    rips.2 <- rips.list[[id.2]]

    # Calculate distance
    dist.i <- TDAstats::phom.dist(rips.1, rips.2)

    # Save the results
    row.ind <- which(rownames(dist.deg0) == id.1)
    col.ind <- which(rownames(dist.deg0) == id.2)

    dist.deg0[row.ind, col.ind] <- dist.i[1]
    dist.deg1[row.ind, col.ind] <- dist.i[2]
  }

  # Convert the distance matrices to symmetric matrices
  dist.deg0 <- dist.deg0 + t(dist.deg0)
  dist.deg1 <- dist.deg1 + t(dist.deg1)

  # Duplicate the distance matrices to remove NAs
  dist.deg0.v2 <- dist.deg0
  dist.deg1.v2 <- dist.deg1

  # Check for NAs. Remove samples with missing homologies
  if (any(is.na(dist.deg0)) | any(is.na(dist.deg1))) {

    # Save just the rows with complete cases
    dist.deg0.v2 <- dist.deg0[stats::complete.cases(dist.deg0) & stats::complete.cases(dist.deg1),
                              stats::complete.cases(t(dist.deg0)) & stats::complete.cases(t(dist.deg1))]
    dist.deg1.v2 <- dist.deg1[stats::complete.cases(dist.deg0) & stats::complete.cases(dist.deg1),
                              stats::complete.cases(t(dist.deg0)) & stats::complete.cases(t(dist.deg1))]

    # Check
    if (!all(rownames(dist.deg0.v2) == colnames(dist.deg0.v2) &
             rownames(dist.deg1.v2) == colnames(dist.deg1.v2) &
             rownames(dist.deg0.v2) == rownames(dist.deg1.v2) &
             colnames(dist.deg0.v2) == colnames(dist.deg1.v2))) {
      stop("Dimension naming mismatch in similarity matrix construction. Please check.")
    }
  }

  # Construct kernel (similarity) matrix
  K.list <- list(dim0 = MiRKAT::D2K(dist.deg0.v2), dim1 = MiRKAT::D2K(dist.deg1.v2))

  # Add the PIDs back to the rows/columns
  for (i in 1:length(K.list)) {
    rownames(K.list[[i]]) <- rownames(dist.deg0.v2)
    colnames(K.list[[i]]) <- colnames(dist.deg0.v2)
  }

  # Return
  K.list

}

#' Calculate a similarity matrix based on the Rips filtration
#'
#' Given an input of data, computes a Rips filtration and persistence diagram for
#' each image. Then constructs a pairwise similarity matrix for each homology group
#' (connected components and loops) based on the resulting persistence diagrams.
#'
#' @param images.df A data.frame containing the image information. See details.
#' @param max.threshold The maximum radius of the circles in the Rips filtration
#' @param print.progress Boolean indicating whether progress in constructing the
#' Rips filtrations across images should be printed.
#'
#' @details `images.df` should contain a `PID` column indicating which sample
#' each image corresponds to. It should also contain columns `x` and `y` indicating
#' the location of each cell. Each row in `images.df` corresponds to a cell within
#' each image. See the package vignettes for an example of the structure.
#'
#' To choose `max.threshold`, a reasonable first choice is to select the maximum
#' possible radius that would cover all cells in each image. This could be the
#' maximum width of the largest image, for example.
#'
#' This function should be run as a first step of the TopKAT analysis. After that,
#' the resulting kernel matrices may be fed into the `TopKAT` function.
#'
#' @return Returns a list of similarity (kernel) matrices and a list of the
#' resulting persistence diagrams for each image.
#'
#' @export
#'
#' @examples
#' # Construct a Rips filtration for each simulated image
#' simmat <- rips_similarity_matrix(data1.df, max.threshold = 100, print.progress = TRUE)
rips_similarity_matrix <- function(images.df, max.threshold, print.progress = TRUE) {

  # Save the PIDs
  PIDs <- unique(images.df$PID)

  # Save the number of samples
  n.sample <- length(PIDs)

  # Initialize a list to store the Rips filtrations
  rips.list <- lapply(1:n.sample, function(i) list())

  # Name the list
  names(rips.list) <- PIDs

  # Iterate through the samples to construct Rips filtration
  # for (i in 1:n.sample) {
  for (i in PIDs) {

    # Print progress
    if (print.progress) print(paste0("Rips diagram: ", i))

    # Subset the data to just this PID
    data.i <- images.df %>%
      dplyr::filter(PID == i) %>%
      dplyr::select(x,y)

    # Construct a Rips filtration using TDAstats
    rips.i <- TDAstats::calculate_homology(data.i, dim = 1, threshold = max.threshold)

    # Save
    rips.list[[i]] <- rips.i
  }

  # Contruct similarity matrix
  K.list <- similarity_matrix(rips.list, n.sample, PIDs)

  # Return
  return(list(K.list = K.list, rips.list = rips.list))
}


#' Topological Kernel Association Test (TopKAT)
#'
#' Perform global test of association between the geometric (topological)
#' structures of spatially-resolved images of cells with continuous, binary,
#' or survival outcomes.
#'
#' @param y Vector of outcomes. Must be numeric. For continuous and survival data,
#' should be a numeric vector. For a binary outcome, must consist of 0s and 1s.
#' @param X Matrix of covariates to adjust for. May be left as NULL.
#' @param cens Vector of event indicators for a survival outcome. 1 indicates a
#' sample experienced the event, 0 otherwise. If not using a survival outcome,
#' leave as NULL.
#' @param K.list List of 2 kernel matrices corresponding to similarities
#' among connected components and among loops. May provide only 1 kernel matrix if
#' interested in a specific homology.
#' @param omega.list Vector of weights to create different combinations of the kernel matrices.
#' Some suggested options: c(0, 1) then TopKAT will combine the p-values for a test associating the
#' connected components with y and the loops with y separately. c(0, 0.5, 1) will combine p-values
#' across just the connected components, an even split of connected components and loops,
#' and just loops.
#' @param outcome.type What kind of outcome is y? Options include "continuous",
#' "binary", or "survival".
#' @param method Specifies how the kernel p-value should be calculated. Options
#' include "davies" for an analytical p-value calculation or "permutation". Only
#' for continuous or binary outcomes.
#' @param omnibus If providing multiple kernel matrices, specify how the p-values
#' should be combined. Defaults to "cauchy" for the Cauchy combination test.
#' @param perm Required only for survival outcomes. Specifies if permutations should
#' be used. Default is FALSE and Davies method is used.
#' @param nperm If using permutations, how many permutations should be used?
#'
#' @import TDAstats
#' @import MiRKAT
#'
#' @return Returns a list of the following objects:
#' overall.pval: the overall p-value describing the association between similarities
#' in topological structures and clinical outcomes,
#' p.vals: the vector of individual p-values for each weight in omega.list,
#' y: the outcome provided,
#' X: the covariates provided,
#' omega.list: the vector of weights provided,
#' outcome.type: the outcome type specified
#'
#' @export
#'
#' @examples
#' # First, construct the similarity matrix
#' simmat <- rips_similarity_matrix(data1.df, max.threshold = 100, print.progress = TRUE)
#'
#' # Then, run TopKAT
#' res <- TopKAT(y = y,
#'               cens = cens,
#'               K.list = simmat$K.list,
#'               omega.list = c(0, 0.5, 1),
#'               outcome.type = "survival")
#'
#' # Check result
#' res$overall.pval
TopKAT <- function(y, X = NULL, cens = NULL, K.list, omega.list, outcome.type = "continuous",
                   method = "davies", omnibus = "cauchy", perm = FALSE, nperm = 999) {

  # Initialize a list of kernel matrices
  K.aggregate <- lapply(1:length(omega.list), function(i) list())
  names(K.aggregate) <- paste0("omega.", omega.list)

  # Compute the trace of each matrix
  trace.list <- sapply(1:length(K.list), function(i) sum(diag(K.list[[i]])))

  # Iterate through the omegas
  for (omega in omega.list) {

    # Create an aggregate kernel with scaling
    K.aggregate[which(omega.list %in% omega)][[1]] <- (1 - omega) * (1/trace.list[1]) * K.list[[1]] + omega * (1/trace.list[2]) * K.list[[2]]

  }

  # Run the kernel machine regression framework through MiRKAT
  if (outcome.type == "continuous") {
    res <- MiRKAT::MiRKAT(y = y, X = X, Ks = K.aggregate, out_type = "C", method = method, omnibus = omnibus, nperm = nperm)
  }

  if (outcome.type == "binary") {
    res <- MiRKAT::MiRKAT(y = y, X = X, Ks = K.aggregate, out_type = "D", method = method, omnibus = omnibus, nperm = nperm)
  }

  if (outcome.type == "survival") {
    res <- MiRKAT::MiRKATS(obstime = y, delta = cens, X = X, Ks = K.aggregate, perm = perm, omnibus = omnibus, nperm = nperm)
  }

  # Save the results
  p.vals <- res$p_values
  overall.pval <- res$omnibus_p

  # Return the results
  list(overall.pval = overall.pval,
       p.vals = p.vals,
       y = y,
       X = X,
       K.list = K.list,
       omega.list = omega.list,
       outcome.type = outcome.type)

}
