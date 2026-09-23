
# -------------------------------------------------------------------------
# Simulation of metacommunities
# -------------------------------------------------------------------------

source("R/setup.R")


# Simulation parameters ---------------------------------------------------

# The number of replicates of the experiment.
R <- 100

detectCores(logical = FALSE)
cl <-  makeCluster(8, type = "MIRAI")
# cl <- NULL

# Exporting data to workers
clusterExport(cl, c("idw", "niches", "sigmas", "k", "dlist", "nspec", "a_sim", "time"))

clusterEvalQ(cl, {
  library(terra)
  library(mrangr)
  id <- unwrap(idw)
  ls()
})

# Community simulation
fname_sim <- paste0("data/niches_", nspec, ".qs2")

com_rep <- chunked_sim(X = R, fun = niches, param = sigmas, fname = fname_sim, cl = cl)
dim(com_rep)

stopCluster(cl)
