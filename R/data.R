#' Simulated data 1
#'
#' This dataset contains simulated data. The images were split into two groups.
#' The first group was simulated to have large squares. The second group was
#' simulated to have loops. This dataset contains 100 samples.
#'
#' @format ## `data1.df`
#' A data frame with 29951 rows and 4 columns:
#' \describe{
#'   \item{PID}{Patient ID}
#'   \item{id}{Image ID (which matches Patient ID)}
#'   \item{x}{x-coordinate of each cell}
#'   \item{y}{y-coordinate of each cell}
#' }
"data1.df"

#' Simulated data 2
#'
#' This dataset contains simulated data. The images were simulated using the
#' `scSpatialSIM` package to represent different tissue structures. Specifically,
#' we used `scSpatialSIM` to simulate two different tissues and separated these
#' into different images. This dataset contains 100 samples.
#'
#' @format ## `data2.df`
#' A data frame with 100481 rows and 4 columns:
#' \describe{
#'   \item{PID}{Patient ID}
#'   \item{id}{Image ID (which matches Patient ID)}
#'   \item{x}{x-coordinate of each cell}
#'   \item{y}{y-coordinate of each cell}
#' }
"data2.df"

#' Outcome vector y
#'
#' This vector contains outcomes for the simulated data. These are simulated
#' survival outcomes from the exponential distribution. For 100 samples, we
#' simulated event times for two different exponential distributions.
#'
#' @format ## `y`
#' A vector of length 100 of event times.
"y"

#' Event indicator vector cens
#'
#' This vector contains binary indicators reflecting whether or not a sample
#' had the event. If the vector equals 1, this indicates the corresponding
#' patient had the event and is 0 otherwise.
#'
#' @format ## `cens`
#' A vector of length 100 of event indicators.
"cens"
