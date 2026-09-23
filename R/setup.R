
# -------------------------------------------------------------------------
# Setup
# -------------------------------------------------------------------------


# Libraries ----


library(terra)
library(mrangr)
library(mirai)
library(pbapply)
library(parallel)
library(qs2)

qopt("nthreads", 4L)

# Helper functions ----
source("R/utils.R")
source("R/functions.R")

# Graphical parameters ----
library(RColorBrewer)
# display.brewer.all(type = "qual")
pal <- brewer.pal(8, "Dark2")
col <- adjustcolor(pal, alpha.f = 0.5)
colt <- adjustcolor(pal, alpha.f = 0.2)


# Simulation environment --------------------------------------------------

# Template raster ----
nrows <- ncols <- 100
xmin <- 250000; xmax <- xmin + nrows * 1000
ymin <- 600000; ymax <- ymin + ncols * 1000
id <- rast(nrows = nrows, ncols = ncols, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
crs(id) <- "epsg:2180"
id
idw <- wrap(id)


# A list of dispersal distances ----
fname <- "data/dlist.rds"
recompute <- TRUE

if (file.exists(fname)) {
  dlist <- readRDS(fname)
  if (length(dlist) == nrows * ncols) recompute <- FALSE
}

if (recompute) {
  dlist <- calculate_dlist(id)
  saveRDS(dlist, file = fname)
}


# Community parameters ----------------------------------------------------

# Number of species in the community
nspec <- 5

# No. of environmental variables
k <- 4

# Simulation time
time <- 40

# Range of sd-s for the interaction matrix tested.
sigmas <- seq(0, 3, 0.25)

