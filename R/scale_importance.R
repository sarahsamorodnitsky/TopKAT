#' Generate a list of persistence diagrams
#'
#' Construct a Rips filtration for each image and return a list of persistence diagrams.
#' This is a helper function for `scale_importance`.
#'
#' @param images.df A data.frame containing the image information. See details.
#' @param max.threshold The maximum radius of the circles in the Rips filtration
#' @param print.progress Boolean indicating whether progress in constructing the
#' Rips filtrations across images should be printed.
#'
#' @return Returns a list of persistence diagrams for each image.
#' @export
#' @examples
#' # Generate a persistence diagram based on a Rips filtration for each image
#' pd.list <- generate_rips(data1.df, 100)
generate_rips <- function(images.df, max.threshold, print.progress = TRUE) {

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

  # Return
  rips.list

}

#' Quantifying the importance of various scale parameters
#'
#' Identify which radii are most important in characterizing an association
#' between the topological structure in each image and outcomes.
#'
#' @param pd.list List of persistence diagrams
#' @param y Outcome vector
#' @param X Covariates to adjust for, if desired. May be left NULL.
#' @param cens Censoring vector for survival outcomes. May be left NULL.
#' @param omega.list Vector of weights to combine kernel matrices
#' @param threshold Maximum radius for Rips filtration
#' @param PIDs Vector of patient IDs
#' @param outcome.type Outcome type, options include "continuous", "binary", or
#' "survival"
#' @param print.progress Boolean, should progress be printed throughout analysis?
#'
#' @details The arguments for `method`, `omnibus`, `perm`, and `nperm` should match
#' what was used in the TopKAT analysis.
#'
#' @return Returns a list with the following elements:
#' \describe{
#' \item{min.thresh}{The radius value at which the lowest TopKAT p-value was obtained}
#' \item{pvals}{The vector of TopKAT p-values for each radius value}
#' \item{threshold.seq}{The vector of maximum radii considered}
#' }
#' @export
#'
#' @examples
#' # Generate a persistence diagram based on a Rips filtration for each image
#' pd.list <- generate_rips(data1.df, 100)
#' # Run the scale importance analysis
#' data1.scale <- scale_importance(pd.list = pd.list,
#'   y = y,
#'   cens = cens,
#'   omega.list = c(0, 0.5, 1),
#'   threshold = 100,
#'   PIDs = 1:100,
#'   outcome.type = "survival")
#'
#' # Plot the results
#' plot(data1.scale$threshold.seq, data1.scale$pvals); abline(v = data1.scale$min.thresh)
scale_importance <- function(pd.list, y, X = NULL, cens = NULL, omega.list, threshold,
                             PIDs, outcome.type = "continuous", n.thresh = 50,
                             print.progress = FALSE) {

  # How many samples are there?
  n.sample <- length(PIDs)

  # Name the list
  if (!all(1:n.sample %in% PIDs)) {
    PIDs.new <- paste0("PID.", PIDs)
  } else {
    PIDs.new <- PIDs
  }

  # Compute a sequence of thresholds
  threshold.seq <- seq(0, threshold, length.out = n.thresh)

  # Initialize a vector of p-values for each threshold
  pvals <- c()

  # Iterate through the thresholds
  for (tt in threshold.seq) {

    if (print.progress) print(paste0(which(threshold.seq %in% tt), "/",
                                     length(threshold.seq)))

    # For each PD, remove all features with birth/deaths after tt
    pd.list.thresh <- lapply(pd.list, function(pd) {

      # Threshold
      pd.temp <- pd[pd[,2] <= tt & pd[,3] <= tt,,drop=FALSE]

      # Check if any CCs or loops are missing
      if (nrow(pd.temp[pd.temp[,1] == 0,,drop=FALSE]) == 0) {
        pd.temp <- rbind(pd.temp, c(0, 0, 0))
      }

      if (nrow(pd.temp[pd.temp[,1] == 1,,drop=FALSE]) == 0) {
        pd.temp <- rbind(pd.temp, c(1, 0, 0))
      }

      # Return
      pd.temp

    })
    names(pd.list.thresh) <- PIDs.new

    # Construct kernel matrices
    K.list <- similarity_matrix(pd.list.thresh, n.sample, PIDs.new)

    # Check if any kernel matrices were entirely 0
    all.0 <- sapply(K.list, function(K) all(K == 0))

    # If both were 0, return an NA
    if (all(all.0)) {
      pvals <- c(pvals, NA)
    }

    # If just one was entirely 0
    if (any(all.0) & !all(all.0)) {

      # Remove it from the list
      K.list <- K.list[!all.0]

      # Run the kernel machine regression framework through MiRKAT
      if (outcome.type == "continuous") {
        res <- MiRKAT::MiRKAT(y = y, X = X, Ks = K.list[[1]], out_type = "C")
      }

      if (outcome.type == "binary") {
        res <- MiRKAT::MiRKAT(y = y, X = X, Ks = K.list[[1]], out_type = "D")
      }

      if (outcome.type == "survival") {
        res <- MiRKAT::MiRKATS(obstime = y, delta = cens, X = X, Ks = K.list[[1]])
      }

      # Save the result
      pvals <- c(pvals, res$p_values)

    }

    # If both were non-zero
    if (!any(all.0)) {

      # Run TopKAT
      res <- TopKAT(y = y,
                    X = X,
                    cens = cens,
                    K.list = K.list,
                    omega.list = omega.list,
                    outcome.type = outcome.type)

      # Save the resulting p-value
      pvals <- c(pvals, res$overall.pval)
    }
  }

  # Remove NAs
  pvals2 <- pvals[!is.na(pvals)]
  threshold.seq2 <- threshold.seq[!is.na(pvals)]

  # Which scale yielded the smallest p-value?
  min.thresh <- threshold.seq2[which.min(pvals2)]

  # Return
  list(min.thresh = min.thresh, pvals = pvals2, threshold.seq = threshold.seq2)
}
