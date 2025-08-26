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

#' MIBI Data from study in Triple Negative Breast Csancer
#'
#' This dataset contains simulated data. The images were split into two groups.
#' The first group was simulated to have large squares. The second group was
#' simulated to have loops. This dataset contains 100 samples.
#'
#' @format ## `tnbc`
#' A data frame with 175497 rows and 8 columns:
#' \describe{
#'   \item{SampleID}{Patient ID}
#'   \item{Group}{Indicates cell type. 2=immune cell, 6=keratin-positive tumor cell}
#'   \item{immuneGroup}{Indicates more granular immune cell type defintion}
#'   \item{x}{x-coordinate of each cell}
#'   \item{y}{y-coordinate of each cell}
#'   \item{Survival_days_capped*}{Time-to-death in days}
#'   \item{Censored}{Survival event indicator (1=had event; 0=censored)}
#'   \item{Class}{Tumor microenvironment class. 0=mixed; 1=compartmentalized; 2=cold}
#' }
#' @references Keren, L., Bosse, M., Marquez, D., Angoshtari, R., Jain, S., Varma, S., ... & Angelo, M. (2018). A structured tumor-immune microenvironment in triple negative breast cancer revealed by multiplexed ion beam imaging. Cell, 174(6), 1373-1387.
"tnbc"
