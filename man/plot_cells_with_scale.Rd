% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_cells_with_scale}
\alias{plot_cells_with_scale}
\title{Plot cells with radius value}
\usage{
plot_cells_with_scale(image, threshold, title)
}
\arguments{
\item{image}{A data.frame containing the cell information for a single image.}

\item{threshold}{Radius to plot}

\item{title}{Title for the plot}
}
\value{
Returns a ggplot scatterplot of the data with edges of length
\code{threshold} connecting the cells.
}
\description{
This function plots a simplicial complex with the cells as points and the
given radius value connecting them.
}
\examples{
require(magrittr)
plot_cells_with_scale(
  image = data1.df \%>\% dplyr::filter(PID == 51),
  threshold = 10,
  title = "Example Plot of Cells"
)
}
