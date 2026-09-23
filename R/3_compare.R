
# -------------------------------------------------------------------------
# Reconstruction of fundamental niches: comparing estimates to the truth
# -------------------------------------------------------------------------

# Setup -------------------------------------------------------------------

library(qs2)
library(abind)

# Helper functions ----
source("R/setup.R")
source("R/utils.R")
source("R/functions.R")


# Comparing fundamental and realised niches -------------------------------
# Process in chunks.

# Locate model files
model_path <- "data/models_001"

# Locate simulation chunks
file_list <- list.files(path = "data", pattern = "niches_\\d+_batch_\\d+\\.qs2$", full.names = TRUE)

# Sort files so batch_2 comes before batch_10
file_list <- gtools::mixedsort(file_list)

# Initialise lists to hold the chunked results
a_list <- magnitude_list <- sign_agreement_list <- spat_pred_list <- list()


for (file_path in file_list) {

  # Extract the clean file name (e.g., "niches_5_batch_1")
  file_label <- tools::file_path_sans_ext(basename(file_path))
  cat("Processing:", file_label, "\n")

  # Load the `sim_com_results` chunk
  com_chunk <- qs2::qs_read(file_path)

  # Extract desired slots from simulation objects
  extracted_a <- get_sim_results(obj = com_chunk, what = "a")

  # Load the corresponding model files for this batch
  model_pattern <- paste0("^", file_label, "_model_\\d+\\.qs2$")
  model_files <- list.files(path = model_path, pattern = model_pattern, full.names = TRUE)
  model_files <- sort(model_files) # Ensure models  are in exact order

  # Read the models into a list
  model_chunk <- lapply(model_files, safe_qs_read)

  # Prepare arrays for mapply_array
  # `mapply_array` requires X and Y to have identical dimensions.
  dim(model_chunk) <- dim(com_chunk)

  # model_chunk[[1, 1]]
  # com_chunk[[1, 1]]

  # Apply comparison function pairwise
  magnitude_compare <- mapply_array(
    FUN = compare_niches,
    X = com_chunk,
    Y = model_chunk,
    type = "beta",
    measure = mae
  )

  sign_compare <- mapply_array(
    FUN = compare_niches,
    X = com_chunk,
    Y = model_chunk,
    type = "beta",
    measure =  \(x, y) mean(sign(y) == sign(x))
  )

  spat_compare <- mapply_array(
    FUN = compare_niches,
    X = com_chunk,
    Y = model_chunk,
    type = "spatial"
  )

  # Store the result
  a_list[[file_label]] <- extracted_a
  magnitude_list[[file_label]] <- magnitude_compare
  sign_agreement_list[[file_label]] <- sign_compare
  spat_pred_list[[file_label]] <- spat_compare

  # Memory management
  rm(com_chunk, model_chunk, magnitude_compare, sign_compare, spat_compare, extracted_a)
  gc()
}


# Strength of biotic interactions ----
a <- abind(a_list, along = 4)
dim(a)

# Absolute mean for each simulation.
# This value measures the average strength of biotic interactions within a community.
inter <- apply(a, 3:4, \(x) mean(abs(x), na.rm = TRUE))
dim(inter)

# The magnitude of the error in niche estimation ----
magnitude <- abind(magnitude_list, along = 3)
dim(magnitude)

# The average magnitude of the bias in slope estimation across species and simulations
magnitude <- apply(magnitude, 2:3, \(x) mean(x, na.rm = TRUE))
dim(magnitude)
range(magnitude, na.rm = TRUE)

# Directional accuracy (sign agreement rate) ----
sign_agreement <- abind(sign_agreement_list, along = 3)
sign_agreement <- apply(sign_agreement, 2:3, \(x) mean(x, na.rm = TRUE))
dim(sign_agreement)
range(sign_agreement, na.rm = TRUE)

# Correlation of spatial distribution ----
spat_pred <- abind(spat_pred_list, along = 3)
dim(spat_pred)
spat_pred <- apply(spat_pred, 2:3, \(x) mean(abs(x), na.rm = TRUE))
dim(spat_pred)


# Save -------------------------------------------------------------------

suffix <- strsplit(basename(model_path), "_", fix = TRUE)[[1]][2]
to_save <- c("inter", "magnitude", "sign_agreement", "spat_pred")
fnames <- paste0("results/", to_save, "_", suffix, ".rds")
mapply(\(x, y) saveRDS(get(x), file = y), to_save, fnames)
