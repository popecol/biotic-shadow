
# -------------------------------------------------------------------------
# Utility functions
# -------------------------------------------------------------------------

chunked_sim <- function(X, fun, param, fname, step_size = NULL, rm_batches = FALSE, cl = NULL, ...) {
  # Manages long-running tasks by batching them into temporary files and aggregating upon completion.
  #
  # Arguments:
  #          X: Integer (replicates) or collection (list/array) to iterate over.
  #        fun: Function. The function to execute on parameters or data.
  #      param: Vector or list. Parameters mapped over in map mode or passed to the function in
  #             simulation mode.
  #      fname: Character. The file path to save/resume results (.qs2 format).
  #  step_size: Integer (optional). Number of items to process per batch. If NULL, it is
  #             calculated automatically.
  # rm_batches: Logical. Whether to delete batch files after aggregation.
  #         cl: Cluster object or integer (optional). For parallel execution.
  #        ...: Additional arguments passed to the execution function.
  #
  # Returns:
  #   An array (or list-array) of results. Dimensions of input 'X' are preserved in
  #   map mode where possible.
  #
  # Details:
  #   The function runs in simulation mode if 'X' is a single integer, replicating 'fun' over 'param'.
  #   If 'X' is a collection, it runs in map mode, applying 'fun' to elements of 'X'.
  #   Additional arguments provided via the ellipsis (...) are explicitly forwarded to 'fun'
  #   through the parallel apply functions.

  # 1. Mode detection and setup ----
  is_sim_mode <- is.numeric(X) && length(X) == 1L
  R <- if (is_sim_mode) X else length(X)

  # Recycles or validates 'param' for map mode to ensure safe and parallel subsetting.
  if (!is_sim_mode) {
    if (length(param) == 1L) {
      param <- rep(param, R)
    } else if (length(param) != R) {
      stop("In map mode, 'param' must be of length 1 or equal to the length of 'X'.")
    }
  }

  # 2. Checks if final result already exists ----
  if (file.exists(fname)) {
    message("Found completed file: ", fname)
    return(qs2::qs_read(fname))
  }

  # 3. Step size optimization ----
  if (is.null(step_size)) {
    # Determines worker count.
    if (is.null(cl)) {
      n_workers <- 1L
    } else if (inherits(cl, "cluster")) {
      n_workers <- length(cl)
    } else {
      # Fallback if cl is just a number.
      n_workers <- as.integer(cl)
    }

    # Calculates optimal step size (aims for ~10 batches total).
    raw_step <- ceiling(R / 10)

    # Aligns to cluster size (rounds up to the nearest multiple of n_workers).
    if (n_workers > 1L) {
      step_size <- ceiling(raw_step / n_workers) * n_workers
    } else {
      step_size <- raw_step
    }

    # Enforces minimums and maximums.
    step_size <- max(step_size, n_workers)
    if (step_size > R) {
      step_size <- ceiling(R / n_workers) * n_workers
    }

    step_size <- as.integer(step_size)
    message(sprintf("Calculated optimal step_size: %d (Workers: %d, Total: %d)",
                    step_size, n_workers, R))
  }

  # 4. Chunked loop (split-file) ----
  chunks <- seq(1L, R, by = step_size)
  base_path <- tools::file_path_sans_ext(fname)

  if (!requireNamespace("pbapply", quietly = TRUE)) stop("Package 'pbapply' is required.")
  if (!requireNamespace("qs2", quietly = TRUE)) stop("Package 'qs2' is required.")

  for (i in chunks) {
    batch_fname <- paste0(base_path, "_batch_", i, ".qs2")
    idx <- i:min(i + step_size - 1L, R)

    if (file.exists(batch_fname)) {
      # Uses idx[1L] and idx[length(idx)] to avoid unnecessary array scans.
      message(sprintf("Batch %d-%d found. Skipping.", idx[1L], idx[length(idx)]))
      next
    }

    message(sprintf("Running batch %d-%d...", idx[1L], idx[length(idx)]))

    batch_res <- tryCatch({
      if (is_sim_mode) {
        # Use pblapply to properly serialize the ellipsis to workers.
        pbapply::pblapply(seq_along(idx), function(dummy, ...) {
          lapply(param, fun, ...)
        }, ..., cl = cl)
      } else {
        # Pairs elements locally into chunks before parallel execution.
        chunk_data <- lapply(idx, function(j) list(X[[j]], param[[j]]))

        # Passes the ellipsis into the worker environment.
        pbapply::pblapply(chunk_data, function(args, ...) {
          fun(args[[1L]], args[[2L]], ...)
        }, ..., cl = cl)
      }
    }, error = function(e) {
      message("  Error in this batch: ", e$message)
      replicate(length(idx), e, simplify = FALSE)
    })

    qs2::qs_save(batch_res, batch_fname)
  }

  # 5. Finalize (aggregate) ----
  message("All batches complete. Aggregating into final file...")

  # Reconstructs the full list from batch files.
  com_rep_list <- vector("list", length(chunks))
  for (k in seq_along(chunks)) {
    i <- chunks[k]
    batch_fname <- paste0(base_path, "_batch_", i, ".qs2")
    com_rep_list[[k]] <- qs2::qs_read(batch_fname)
  }

  # Flattens (list of batches to list of results).
  com_rep_list <- unlist(com_rep_list, recursive = FALSE)

  # Checks for errors safely using vapply.
  errors_found <- vapply(com_rep_list, function(x) inherits(x, "error"), logical(1L))
  if (any(errors_found)) {
    warning(sprintf("Process finished but %d items contain errors. Returning list.", sum(errors_found)))
    return(com_rep_list)
  }

  message("Simplifying...")

  if (is_sim_mode) {
    # Simulation mode (standard simplification from list to matrix/array).
    res <- simplify2array(com_rep_list)

  } else {
    # Map mode (structure preservation).
    if (!is.null(dim(X))) {
      # Tries simplifying.
      res <- simplify2array(com_rep_list)

      # Discards simplification and creates a list-array if dimensions mismatch.
      if (!identical(dim(res), dim(X))) {
        message("Simplified result differs from input dimensions. Creating list-array.")
        res <- com_rep_list
        dim(res) <- dim(X)
        if (!is.null(dimnames(X))) dimnames(res) <- dimnames(X)
      }
    } else {
      # Flat input (vector or list).
      res <- simplify2array(com_rep_list)
    }
  }

  # 6. Saves and cleans up ----
  qs2::qs_save(res, fname)
  message("Final file saved: ", fname)

  if (rm_batches) {
    message("Removing temporary batch files...")
    batch_files <- paste0(base_path, "_batch_", chunks, ".qs2")
    file.remove(batch_files)
  }

  return(invisible(res))
}


get_sim_results <- function(obj, what) {
  # Extracts a specific slot from a collection of 'sim_com_results' objects.

  # Arguments:
  #    obj: multiple 'sim_com_results' objects (as a list, matrix, or array)
  #   what: slot name; one of: "extinction", "sim_time", "id", "N_map", but
  #         can be any slot present in 'sim_com_results'.

  # Returns the contents of the slot, simplified into a vector, matrix, or array.
  #   If 'obj' has dimensions (e.g., a matrix), these are preserved in the output.
  #   Note: if 'what' is "id", only the slot from the first object is returned.

  # Check if all elements inherit from "sim_com_results"
  if (!all(vapply(obj, inherits, logical(1), "sim_com_results"))) {
    stop("All elements in 'obj' must inherit from the 'sim_com_results' class")
  }

  # Validate the 'what' argument
  what <- match.arg(what, names(obj[[1]]))

  # If the template map is required, it is always returned as a single wrapped raster.
  if (what == "id") {
    return(obj[[1]][[what]])
  }

  extracted_list <- lapply(obj, `[[`, what)
  result <- simplify2array(extracted_list)

  # Restore dimensions if 'obj' was a matrix/array
  if (!is.null(dim(obj))) {
    # Check if result is just a vector (scalars extracted from a matrix obj)
    if (is.null(dim(result))) {
      dim(result) <- dim(obj)
    } else {
      internal_dims <- dim(result)[-length(dim(result))]
      dim(result) <- c(internal_dims, dim(obj))
    }
  }
  return(result)
}


calculate_dlist <- function(x) {
  # Calculates dlist for later use
  # Argument
  #   x: a template raster

  values(x) <- 1
  dat <- rangr::initialise(x, x, r = 1, max_dist = diagonal(x))
  dat$dlist
}


safe_qs_read <- function(file) {
  # Safely reads a qs2 file and handles potential errors.
  #
  # Arguments:
  #   file: Path to the file to be read (character).
  #
  # Returns:
  #   The parsed object from the file, or NA if an error occurs.
  #
  # Details:
  #   This function wraps qs2::qs_read within a tryCatch block. It intercepts
  #   any reading errors (which typically occur with corrupted or empty files)
  #   and outputs a descriptive warning rather than halting execution.

  tryCatch({
    qs2::qs_read(file)
  }, error = function(e) {
    # If the file is corrupted or empty, print a warning and return NA.
    warning("Skipping corrupted file: ", basename(file), call. = FALSE)
    NA
  })
}


mapply_array <- function(FUN, X, Y, pb = FALSE, MoreArgs = NULL, ...) {

  # Applies a function pairwise over two matrices with a progress bar and
  # reshapes the output into a higher-dimensional array [Species, Parameters, Replicate].
  #
  # Arguments:
  #        FUN: Function to apply. Must return a vector of fixed length (S).
  #       X, Y: Input matrices. Must have identical dimensions.
  #         pb: Logical. If TRUE, uses 'pbapply::pbmapply' to show a progress bar.
  #   MoreArgs: A list of static arguments passed to 'FUN'.
  #        ...: Additional *static* arguments passed to 'FUN' (e.g. flags, constants).
  #             These are NOT iterated over; they are added to MoreArgs.
  #
  # Returns:
  #   An array of dimensions c(S, P, R).
  #   The first dimension corresponds to the vector returned by 'FUN'.
  #   Preserves dimension names from 'X'.
  #
  # Details:
  #   This function acts as a wrapper around 'mapply' (or 'pbmapply'). It ensures
  #   that scalar results from 'FUN' are treated as 1-row matrices, preventing
  #   dimension collapse. Unlike standard 'mapply', '...' arguments are treated
  #   as static (constant) inputs to 'FUN'.

  # 1. Validation ----
  if (!identical(dim(X), dim(Y))) {
    stop("Dimensions of 'X' and 'Y' must be identical.")
  }

  if (pb && !requireNamespace("pbapply", quietly = TRUE)) {
    warning("Package 'pbapply' is not installed. Defaulting to standard mapply.")
    pb <- FALSE
  }

  # 2. Argument Handling ----
  # Treat ... as static arguments (like lapply), not iterated vectors.
  extra_args <- c(as.list(MoreArgs), list(...))

  # 3. Run mapply ----
  if (pb) {
    res <- pbapply::pbmapply(FUN, X, Y, MoreArgs = extra_args, SIMPLIFY = TRUE)
  } else {
    res <- mapply(FUN, X, Y, MoreArgs = extra_args, SIMPLIFY = TRUE)
  }

  # 4. Normalize output ----
  # If FUN returns a scalar, mapply returns a vector of length N.
  # Convert this to a 1 x N matrix to maintain the species dimension.
  if (is.vector(res)) {
    res <- matrix(res, nrow = 1)
  }

  if (!is.matrix(res)) {
    stop("Output is ragged (FUN returned variable lengths) and cannot be reshaped.")
  }

  # 5. Reconstruct dimensions ----
  # Target: c(Species, Original_Rows, Original_Cols)
  final_dims <- c(nrow(res), dim(X))

  # 6. Reconstruct names ----
  input_names <- dimnames(X)
  if (is.null(input_names)) {
    input_names <- rep(list(NULL), length(dim(X)))
  }

  # Combine FUN output names (if any) with Input names
  final_names <- c(list(rownames(res)), input_names)

  # 7. Reshape ----
  array(res, dim = final_dims, dimnames = final_names)
}


read_rds_array <- function(vname, path = ".", simplify = TRUE) {
  # Reads multiple .rds files based on variable names and optionally simplifies them into an array.
  #
  # Arguments:
  #      vname: A character vector of variable names to match (character).
  #       path: The directory path where the .rds files are located (character).
  #   simplify: A logical value indicating whether to simplify the list to an array (logical).
  #
  # Returns:
  #   An array if simplification is requested and dimensions are compatible, or a named list otherwise.
  #
  # Details:
  #   The function constructs a regular expression from the provided variable names
  #   to find matching files. It reads the files into a list and assigns file names as
  #   list names. If requested, it attempts to combine them into an array using 'simplify2array'.

  # Create a regular expression pattern to match filenames.
  pattern <- paste0("^(", paste(vname, collapse = "|"), ").*\\.rds$")

  # List matching files in the specified directory.
  file_list <- list.files(path = path, pattern = pattern, full.names = TRUE)

  # Return early if no files are found to avoid errors.
  if (length(file_list) == 0L) {
    warning("No matching files found.")
    return(NULL)
  }

  # Read all files into a list.
  rds_list <- lapply(file_list, readRDS)

  # Assign file names to the list elements for better identification.
  names(rds_list) <- sub("\\.rds$", "", basename(file_list), ignore.case = TRUE)

  # Simplify the list to an array or return as a list.
  if (simplify) {
    simplify2array(rds_list)
  } else {
    rds_list
  }
}


reshape_2_long <- function(parm_list, col_name = "value") {
  # Converts a nested list of simulation matrices into a long-format data frame.
  #
  # Arguments:
  #   parm_list: A named list of numeric matrices containing simulation results.
  #    col_name: Name to assign to the dependent variable column (character).
  #
  # Returns:
  #   A long-format data frame containing parameter values, trial labels, and matrix contents.
  #
  # Details:
  #   Iterates over the named matrices in 'parm_list', extracts proportion metadata from string
  #   names, constructs a matching long-format grid per matrix using 'sigmas' from the parent
  #   environment, and row-binds them into a single data frame.

  df_list <- lapply(names(parm_list), function(p) {
    mat <- parm_list[[p]]

    # Extract the proportion of sites sampled and generate unique trial labels
    prop_str <- strsplit(p, "_")[[1]]
    prop_str <- prop_str[length(prop_str)]
    # trial_vals <- paste(prop_str, seq_len(ncol(mat)), sep = "_")
    trial_vals <- seq_len(ncol(mat))

    # Build long grid matching R column-major matrix ordering
    df <- expand.grid(sigma = sigmas, trial = trial_vals)

    # Assign dependent variable dynamically using standard subset indexing
    df[[col_name]] <- as.vector(mat)
    df$prop <- as.numeric(paste0("0.", prop_str))

    # Reorder columns
    df[, c("sigma", col_name, "trial", "prop")]
  })

  # Combine elements and format variables
  data_long <- do.call(rbind, df_list)
  # data_long$trial <- as.factor(data_long$trial)
  rownames(data_long) <- NULL

  data_long
}


a_sim <- function(nspec, mu, sigma, digits = 2) {
  # Simulate interaction matrix from a normal distribution
  # Arguments:
  #   nspec: number of species
  #      mu: mean strength of pairwise interactions
  #   sigma: sd of pairwise interactions

  a <- matrix(rnorm(nspec^2, mu, sigma), nrow = nspec, ncol = nspec)
  diag(a) <- NA
  round(a, digits)
}


extreme_betas <- function(x) {
  # Checks for extreme fixed-effect coefficients in a glmmTMB model.

  # This function assesses whether the absolute sum of the fixed-effect
  # coefficients (betas) for any single species within a glmmTMB model
  # exceeds a predefined threshold. Extreme coefficients often indicate
  # issues like quasi-complete separation, particularly when a species is
  # very rare.

  # Argument:
  #   x: a glmmTMB object

  # Returns:
  #   TRUE if the absolute sum of the fixed-effect coefficients for any species
  #   is greater than the arbitrary threshold (currently 30), and FALSE otherwise.

  cf <- glmmTMB::fixef(x)$cond
  cf <- matrix(cf, nrow = nspec, ncol = k + 1)
  any(rowSums(abs(cf)) > 30)
}


is_na_se <- function(x) {
  # Checks for missing standard errors in fixed effects

  # Argument:
  #   x: A fitted object of class "glmmTMB".

  # Returns:
  #   A logical value: TRUE if one or more fixed effect SEs are NA, FALSE otherwise.

  any(is.na(summary(x)[["coefficients"]][["cond"]][, "Std. Error"]))
}


mae <- function(x, y) {
  # Mean absolute deviation (fitted vs. true)
  mean(abs(y - x))
}


quant_gam <- function(x, y, fun = identity, inv_fun = identity,
                      qu = 0.5, k = NULL, newdata = NULL, nthreads = 1) {

  # Estimates conditional quantiles using a Gaussian Location-Scale GAM (GAMLSS).
  #
  # Arguments:
  #         x: Numeric vector of the predictor variable.
  #         y: Numeric vector of the response variable.
  #       fun: Function to transform 'y' to normality (default: identity).
  #   inv_fun: Function to back-transform predictions (default: identity).
  #        qu: Numeric vector of probabilities for quantiles (default: 0.5).
  #         k: Integer (optional). Basis dimension for smooth terms. Defaults to
  #            a heuristic based on sample size.
  #   newdata: Data frame with column 'x' for prediction. If NULL, generates a
  #            regular grid of 100 points.
  #  nthreads: Integer. Number of threads to use for parallel computation
  #            (via OpenMP) during model fitting. Default is 1.
  #
  # Returns:
  #   A data.frame containing the 'newdata' 'x' column and the estimated
  #   quantiles (columns named "q_0.5", etc.).
  #
  # Details:
  #   Fits a distributional model using 'mgcv::gaulss', where both the mean
  #   (mu) and the precision (1/sigma) of the transformed response are modeled
  #   as smooth functions of 'x'. Quantiles are derived analytically from the
  #   predicted normal distribution parameters and back-transformed using 'inv_fun'.
  #   While 'bam' is not supported for this family, 'gam' can use multiple threads
  #   via the 'nthreads' control argument.

  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required for this function.")
  }

  # 1. Prepare data ----
  x <- as.vector(x)
  y <- as.vector(y)

  # Apply transformation
  y_trans <- fun(y)

  data <- data.frame(x = x, y = y_trans)
  data <- data[complete.cases(data), , drop = FALSE]
  n_obs <- nrow(data)

  if (n_obs < 4) stop("Not enough data points to fit a model.")

  # 2. Heuristics for basis dimension (k) ----
  if (is.null(k)) {
    k <- min(10, floor(n_obs / 2))
  }

  # 3. Fit GAMLSS model (location-scale) ----
  # formula: list(mean_formula, precision_formula)
  form <- list(y ~ s(x, bs = "ad", k = k), ~ s(x, bs = "ad", k = k))

  # Use gam.control to set nthreads
  ctrl <- mgcv::gam.control(nthreads = nthreads)

  fit <- mgcv::gam(form, family = mgcv::gaulss(), data = data, control = ctrl)

  # 4. Prediction ----
  if (is.null(newdata)) {
    newdata <- data.frame(x = seq(min(x), max(x), length.out = 100))
  } else {
    if (!"x" %in% names(newdata)) stop("'newdata' must contain a column named 'x'.")
  }

  # type = "response" returns:
  # Col 1: mu (mean)
  # Col 2: 1/sigma (inverse standard deviation)
  pre <- predict(fit, newdata = newdata, type = "response")

  mu_hat  <- pre[, 1]
  sig_hat <- 1 / pre[, 2] # Convert precision to SD

  # 5. Calculate quantiles (analytical) ----
  # We work on the transformed scale (Gaussian assumption applies here)
  q_list <- lapply(qu, function(p) {
    stats::qnorm(p, mean = mu_hat, sd = sig_hat)
  })

  # Bind into a matrix/df
  q_mat <- do.call(cbind, q_list)

  # 6. Back-transformation ----
  # Apply inverse function to the predicted quantiles
  q_mat[] <- inv_fun(q_mat)

  # 7. Format Output ----
  df_res <- as.data.frame(q_mat)
  colnames(df_res) <- paste0("q_", qu)

  return(cbind(newdata, df_res))
}


inv_logit <- function(x, percents = FALSE, adjust = 0) {

  # Calculates the inverse of the 'car::logit' transformation.
  #
  # Arguments:
  #         x: Numeric vector of logit values.
  #  percents: Logical. If TRUE, returns values on [0, 100] scale; otherwise [0, 1].
  #    adjust: Numeric. The adjustment factor used to handle 0/1 bounds in the
  #            original transformation (default 0).
  #
  # Returns:
  #   A numeric vector of proportions (or percentages) corresponding to 'x'.
  #
  # Details:
  #   Reverses the logic of `car::logit` by applying the standard logistic sigmoid
  #   (plogis), reversing the linear interval adjustment, and scaling
  #   for percentages if requested. The output is clamped to valid ranges to
  #   handle floating-point precision errors.

  # 1. Standard inverse logit
  p <- stats::plogis(x)

  # 2. Reverse the adjustment
  if (adjust != 0) {
    p <- (p - adjust) / (1 - 2 * adjust)
  }

  # 3. Handle percents
  if (percents) {
    p <- p * 100
  }

  # 4. Safety clamp
  # Removes floating point artifacts (e.g., -1e-16)
  upper_limit <- if (percents) 100 else 1
  p <- pmax(0, pmin(upper_limit, p))

  return(p)
}


plot_qrt <- function(x, y, q = NULL, col = "black", alpha = 0.05,
                     fun = identity, inv_fun = identity, ...) {
  # Visualizes simulation results with a scatter plot and overlaid quartiles.
  #
  # Arguments:
  #         x: Numeric vector of the independent variable (e.g., distances).
  #         y: Numeric matrix (or vector) of simulation results. If a matrix,
  #            rows correspond to 'x' and columns to replicates.
  #         q: Optional data frame. If provided, must contain columns "x",
  #            "q_0.25", "q_0.5", and "q_0.75". If NULL (default), these are
  #            calculated internally using 'quant_gam'.
  #       col: Color identifier (default "black") for the trend line.
  #     alpha: Numeric. Transparency level for raw data points (0-1). Default 0.05.
  #       fun: Transformation function passed to 'quant_gam' if 'q' is NULL.
  #   inv_fun: Inverse transformation function passed to 'quant_gam' if 'q' is NULL.
  #       ...: Additional graphical parameters passed to 'matplot'.
  #
  # Returns:
  #   NULL. The function produces a base graphics plot.
  #
  # Details:
  #   This function produces a composite plot:
  #   1. Raw data: Plotted as a "cloud" using 'matplot' with high transparency.
  #      Handles missing (NA) or NaN values in 'y'.
  #   2. IQR Ribbon: A semi-transparent polygon spanning the 0.25 to 0.75 quantiles.
  #   3. Median: A solid line representing the 0.50 quantile.

  # 1. Validation ----
  if (is.matrix(y)) {
    if (is.matrix(x)) {
      # Case: Both are matrices (e.g., parameter matrices vary per replicate)
      if (!identical(dim(x), dim(y))) {
        stop("Dimensions of 'x' and 'y' must match when both are matrices.")
      }
    } else {
      # Case: x is a vector (shared across all replicates of y)
      if (length(x) != nrow(y)) {
        stop("Length of vector 'x' must match the number of rows in matrix 'y'.")
      }
    }
  } else {
    # Case: y is a vector
    # We treat x as a vector (flattening if it's a matrix) to ensure length match
    if (length(as.vector(x)) != length(y)) {
      stop("Lengths of 'x' and 'y' must match.")
    }
  }

  # 2. Handle quantiles (q) ----
  if (is.null(q)) {
    # quant_gam internally flattens x and y, so matrix inputs are safe.
    q <- quant_gam(x, y, fun = fun, inv_fun = inv_fun, qu = c(0.25, 0.5, 0.75))
  } else {
    req_cols <- c("x", "q_0.25", "q_0.5", "q_0.75")
    if (!all(req_cols %in% names(q))) {
      stop("Argument 'q' must contain columns: 'x', 'q_0.25', 'q_0.5', 'q_0.75'.")
    }
  }

  # 3. Define colors ----
  point_col  <- grDevices::adjustcolor(col, alpha.f = alpha)
  ribbon_col <- grDevices::adjustcolor(col, alpha.f = 0.3)

  # 4. Plot Raw Data
  matplot(x, y, type = "p", pch = 20, cex = 0.8, col = point_col, xlab = "", ...)

  # 5. Plot IQR ribbon ----
  graphics::polygon(
    x = c(q$x, rev(q$x)),
    y = c(q$q_0.25, rev(q$q_0.75)),
    border = NA,
    col = ribbon_col
  )

  # 6. Plot median line ----
  graphics::lines(q$x, q$q_0.5, lty = 1, lwd = 2, col = col)
}
