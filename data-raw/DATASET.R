# -----------------------------------------------------------------------------
# Generate example data for `TopKAT` R package
# -----------------------------------------------------------------------------

# Packages
library(scSpatialSIM)

# Function for generating the data
image_generation <- function(n.sample, dims, group1.shape, group2.shape, Sigma = NULL, sdmax = 15) {

  # Initialize a list of images
  images <- lapply(1:n.sample, function(i) list())

  # Separate the image generation based on real tissue or otherwise --

  # Simulate shapes, CSR, or clusters
  if (group1.shape != "real tissue" & group2.shape != "real tissue") {

    # Initialize a data.frame that stores which shape each image should take
    shape.df <- data.frame(id = 1:n.sample,
                           shape = c(rep(group1.shape, n.sample/2), rep(group2.shape, n.sample/2)))

    # Iterate through the samples
    for (i in 1:n.sample) {

      # Save the current shape needed
      shape.i <- shape.df[shape.df$id == i,]$shape

      # Simulate squares
      if (shape.i == "square") {

        # Simulate the indices for a random number square clusters
        n.clust <- sample(2:6, 1)
        max.width <- round((dims[1]/n.clust))
        left.corners <- sapply(1:n.clust, function(i) c(runif(1, min = round(0.1*dims[1]), max = round(0.9*dims[2])),
                                                        runif(1, min = round(0.1*dims[1]), max = round(0.9*dims[2]))))
        inds <- do.call(rbind, lapply(1:n.clust, function(i) expand.grid(left.corners[1,i] + 0:10, left.corners[2,i] + 0:10)))

      }

      # Simulate loops
      if (shape.i == "loop") {

        # Simulate the number of loops
        n.loop <- sample(2:6, 1)
        max.radius <- round((dims[1]/n.loop))/2
        centers <- sapply(1:n.loop, function(i) c(runif(1, min = round(0.1*dims[1]), max = round(0.9*dims[2])),
                                                  runif(1, min = round(0.1*dims[1]), max = round(0.9*dims[2]))))
        angles <- seq(0, 2 * pi, length.out = 30)
        inds <- do.call(rbind, lapply(1:n.loop, function(i) cbind(centers[1,i] + max.radius * cos(angles),
                                                                  centers[2,i] + max.radius * sin(angles))))

      }

      # Simulate clusters
      if (shape.i == "cluster") {

        # Randomly select means
        mu <- matrix(c(runif(2, 0, dims[1]), runif(2, 0, dims[2])), nrow = 2, byrow = T)

        # Simulate two clusters
        inds <- rbind(MASS::mvrnorm(20^2, mu = mu[,1], Sigma = Sigma),
                      MASS::mvrnorm(20^2, mu = mu[,2], Sigma = Sigma))

      }

      # Simulate CSR
      if (shape.i == "CSR") {

        # Add the cell locations
        inds <- cbind(runif(30^2, 0, dims[1]), runif(30^2, 0, dims[2]))

      }

      # Remove duplicated points
      inds <- unique(inds)

      # Add the cell locations
      image.i <- data.frame(PID = i,
                            id = i,
                            x = inds[,1],
                            y = inds[,2])

      # Add to the data set
      images[[i]] <- image.i

    }

  }

  # Simulate real tissue
  if (group1.shape == "real tissue" & group2.shape == "real tissue") {

    # Iterate through only half the samples to split up Tissues 1 and 2
    for (i in 1:(n.sample/2)) {

      # Create a window
      W <- spatstat.geom::owin(xrange = c(0, dims[1]), yrange = c(0, dims[2]))

      # Initialize image
      sim_object <- scSpatialSIM::CreateSimulationObject(sims = 1, cell_types = 1, window = W) %>%
        scSpatialSIM::GenerateSpatialPattern(lambda = 0.2) %>%
        scSpatialSIM::GenerateTissue(density_heatmap = FALSE, use_window = TRUE, sdmax = sdmax) %>%
        scSpatialSIM::GenerateCellPositivity(density_heatmap = FALSE, probs = c(1, 1), shift = 0)

      # Save the cell indices
      inds1 <- sim_object@`Spatial Files`[[1]] %>%
        dplyr::filter(`Tissue Assignment` == "Tissue 1") %>%
        dplyr::select(x, y)

      inds2 <- sim_object@`Spatial Files`[[1]] %>%
        dplyr::filter(`Tissue Assignment` == "Tissue 2") %>%
        dplyr::select(x, y)

      # Add the cell locations
      image.i1 <- data.frame(PID = i,
                             id = i,
                             x = inds1[,1],
                             y = inds1[,2])

      image.i2 <- data.frame(PID = i + (n.sample/2),
                             id = i + (n.sample/2),
                             x = inds2[,1],
                             y = inds2[,2])

      # Add to the data set
      images[[i]] <- image.i1
      images[[i+(n.sample/2)]] <- image.i2

    }

  }

  # Return the images
  images

}

# -----------------------------------------------------------------------------
# Simulate two versions of data
# (1) connected components vs. loops
# (2) simulated tissue
# -----------------------------------------------------------------------------

# Set a seed
set.seed(1996)

# -------------------------------------
# (1) connected components vs. loops
# -------------------------------------

# Simulate the images
data1 <- image_generation(
  n.sample = 100,
  dims = c(100, 100),
  group1.shape = "square",
  group2.shape = "loop"
)

# Combine into a data.frame
data1.df <- do.call(rbind, data1)

# Add cell types
data1.df$type <- sample(c("cell type 1", "cell type 2", "cell type 3", "cell type 4"),
                        size = nrow(data1.df), replace = TRUE)

# -------------------------------------
# (2) simulated tissue
# -------------------------------------

# This will take a few seconds to run
data2 <- image_generation(
  n.sample = 100,
  dims = c(100, 100),
  group1.shape = "real tissue",
  group2.shape = "real tissue"
)

# Combine into a data.frame
data2.df <- do.call(rbind, data2)

# Add cell types
data2.df$type <- sample(c("cell type 1", "cell type 2", "cell type 3", "cell type 4"),
                        size = nrow(data2.df), replace = TRUE)

# -------------------------------------
# Simulate the outcomes
# -------------------------------------

# Simulate the outcomes
y <- c(rexp(50, log(2)), rexp(50, log(2)/2))
cens <- rbinom(100, 1, 0.9)

# Save the data
usethis::use_data(data1.df, overwrite = TRUE)
usethis::use_data(data2.df, overwrite = TRUE)
usethis::use_data(y, cens, overwrite = TRUE)
