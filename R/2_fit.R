
# -------------------------------------------------------------------------
# Reconstruction of fundamental niches: fitting GLMM do VE data
# -------------------------------------------------------------------------

# Setup -------------------------------------------------------------------

source("R/setup.R")

# Libraries ----
library(mori)
# install.packages("glmmTMB", type = "source")
library(glmmTMB)


# Virtual sampling and model fitting --------------------------------------

# Sample sizes
sim_time <- 30

props <- c(0.001, 0.003, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1)
log(props)
n_props <- length(props)


# Process in chunks ----

# 1. Locate chunked qs2 files
file_list <- list.files(path = "data", pattern = "niches_\\d+_batch_\\d+\\.qs2$", full.names = TRUE)
file_list <- gtools::mixedsort(file_list)

# 2. Initialize workers
daemons(8)

# 3. Setup the workers' environment
everywhere({
  library(mrangr)
  library(glmmTMB)
  library(qs2)
}, fit_func = fit)

prior <- "normal(0, 3)"

# Set the proportion of sampled sites
props
prop <- props[1]

# 4. Process the chunks sequentially, but parallelize the fits within them

# Create a temporary folder to hold the results
if(!dir.exists("data/models")) {
  dir.create("data/models")
}


for (file_path in file_list) {

  # Extract the file name (e.g., "niches_5_batch_1") for unique naming
  file_label <- tools::file_path_sans_ext(basename(file_path))

  # Load the chunk into the main process's RAM
  # file_path <- file_list[1] # For tests
  chunk <- qs2::qs_read(file_path)
  cat("Processing:", file_path, "\n")

  # Place the chunk into OS-level shared memory.
  # This returns an ALTREP wrapper pointing to the shared pages.
  shared_chunk <- share(chunk)

  # Distribute the fits for this chunk across workers
  m_tasks <- mirai_map(
    .x = seq_along(shared_chunk),
    .f = function(i) {

      # Extract the specific object
      obj <- shared_data[[i]]

      # Fit the model
      model_fit <- fit_func(obj = obj, prop = prop, prior = prior)

      # Construct a unique filename for this specific model
      file_name <- sprintf("data/models/%s_model_%03d.qs2", chunk_label, i)

      # 4. Save the model to disk
      qs2::qs_save(model_fit, file_name)

      # 5. Return the filename string
      return(file_name)
    },
    shared_data = shared_chunk,
    prop = prop,
    prior = prior,
    chunk_label = file_label
  )

  # Wait for this chunk's tasks to complete and display the progress bar
  m_tasks[.progress]

  # Memory management
  rm(chunk, shared_chunk)
  gc()
}

# 5. Shut down the daemons
daemons(0)


