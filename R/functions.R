
niches <- function(sigma, progress_bar = FALSE) {
  # This function is the core engine used for running replicated metacommunity
  # simulations. It creates a complete simulation environment by defining
  # habitat conditions, species-habitat relationships, interspecific interactions,
  # and then simulating the community dynamics over time.
  # The primary purpose is to integrate environmental filtering and biotic interactions
  # to model realised metacommunity structure in a spatially explicit context.

  # The simulation proceeds in two main stages:
  #   1. Fundamental niche definition
  #     - Habitat maps: Generates spatially autocorrelated environmental maps (env_map)
  #       using a reference raster (id).
  #     - Niche coefficients: Generates species-specific coefficients (beta) that
  #       define how each species responds to the environmental variables.
  #       These coefficients establish the fundamental niche of each species.
  #     - Carrying capacity: Calculates the carrying capacity maps (K_map) for
  #       all species across the landscape. This is done by applying the coefficients
  #       to the environmental values via an inverse link function. The intercept
  #       is fixed (to 2) to normalise initial abundance.
  #   2. Community simulation
  #     - Interaction matrix: Generates the interspecific interaction matrix (a)
  #       using random numbers from a normal distribution, where the sigma argument
  #       controls the strength (magnitude) of these interactions.
  #     - Initialization: Initializes the simulation data object (data) using the
  #       generated K_map, the interaction matrix (a), and a predefined dispersal
  #       kernel (e.g., exponential distribution).
  #     - Simulation run: Executes the metacommunity simulation for a specified
  #       time period, including a brief burn-in phase.

  # The function returns the complete simulated community dynamics object,
  # enhanced with the generated simulation parameters: the environmental data (env),
  # the fundamental niche coefficients (beta), carrying capacity maps for each species (K),
  # and the interaction matrix (a).

  # 1. Define fundamental niches ----

  # Generate habitat maps.
  env_map <- K_sim(k, id, range = diagonal(id))
  # plot(env_map)

  # Generate the coefficients, i.e. the slopes against the habitat variables.
  # These coefficients define the fundamental niche.
  # Later, they will be reconstructed based on data collected by a virtual ecologist.
  beta <- rnorm(nspec * k, 0, 0.2)
  beta <- matrix(beta, nrow = nspec, ncol = k, byrow = TRUE)

  # The intercept of the linear model used to generate carrying capacity maps
  # is fixed to ensure a similar abundance for all species.
  intercept <- 2
  beta <- cbind(intercept, beta)

  # Calculate the carrying capacity maps.
  env <- values(env_map)
  design_matrix <- cbind(rep(1, nrow(env)), env)
  K <- design_matrix %*% t(beta)
  K <- exp(K) # inverse link
  K <- data.frame(K)
  names(K) <- LETTERS[1:nspec]
  # summary(K)

  # Convert to raster
  K_map <- id
  nlyr(K_map) <- nspec
  values(K_map) <- K
  names(K_map) <- names(K)
  # K_map; plot(K_map)
  # pairs(log(K_map))

  # 2. Simulate the community ----

  # Parameters of the interaction matrix (random numbers from a normal distribution).
  a <- a_sim(nspec, 0, sigma)

  # Intrinsic growth rates
  r <- rnorm(nspec, 1, 0.2)

  # Dispersal kernel
  # curve(dexp(x, rate = 1 / 500), from = 0, to = diagonal(id) / 10, xlab = "Distance")

  # Initialise 'sim_com_data' object
  data <- initialise_com(K_map = K_map, r = r, a = a, rate = 1 / 500, max_dist = diagonal(id), dlist = dlist)
  # summary(data)

  com <- sim_com(data, time, burn = 10, progress_bar = progress_bar)
  # summary(com)
  # plot(com)
  # plot_series(com, col = col, lty = 1, cex.lab = 1.3)

  # Additional slots
  com$env <- env
  com$beta <- beta
  com$K <- K
  com$a <- a

  return(com)
}


fit <- function(obj, prop, prior = NULL) {

  # Fits a GLMM to community data sampled by a virtual ecologist.
  #
  # Arguments:
  #     obj: A `sim_com_results` object, containing an `env` slot (environmental covariates).
  #    prop: A sampling proportion for the virtual ecologist (numeric).
  #   prior: An optional prior specification string (character or NULL).
  #
  # Details:
  #   The function generates synthetic observation data using virtual_ecologist with the
  #   sampling proportion provided by the 'prop' argument and a binomial observation error.
  #   It merges these observations with environmental covariates. Finally, it fits a GLMM
  #   with fixed species-specific environmental effects and a hierarchical random effects
  #   structure (species-specific trends, trends within sites, and random time intercepts).
  #
  #   If a prior is provided, the model performs penalized regression. Specifying a 'normal'
  #   prior corresponds to L2 regularization (ridge penalty), while a 't' prior serves as an
  #   approximation of L1 regularization (lasso penalty). Examples of how to use this include
  #   setting prior = "normal(0, 3)" or prior = "t(0, 1, 3)". Single-threaded execution is
  #   enforced to prevent cluster deadlocks during parallel processing.
  #
  #   Memory management:
  #   Before the GLMM is fitted, the `obj`, `ve`, and `env` variables are deleted from the
  #   local environment. This step is a safeguard when distributing data using the `mori`
  #   package [https://cran.r-project.org/package=mori], which relies on ALTREP pointers.
  #   `glmmTMB` capture a "snapshot" of its local environment (a closure) and store it within
  #   the final model object. If the ALTREP object remains in the environment during fitting,
  #   the memory pointer is saved to the disk alongside the model. When the shared
  #   memory is subsequently purged by the garbage collector, loading the model triggers a
  #   fatal memory error. Explicitly deleting these objects ensures they physically cannot
  #   be captured by `glmmTMB`'s internal closures.

  # Disable OpenMP for this worker
  # This prevents glmmTMB from spawning threads inside the parallel worker.
  if (requireNamespace("TMB", quietly = TRUE)) {
    TMB::openmp(n = 1)
  }

  # Virtual Ecologist
  ve <- virtual_ecologist(obj, prop = prop, obs_error = "rbinom", obs_error_param = 0.5)
  # summary(ve)
  # aggregate(n ~ species, ve, range)

  # Get environmental variables
  env <- obj$env

  # Merges the environmental data and the survey data.
  data <- data.frame(id = seq_len(nrow(env)), env)
  data <- merge(ve, data, by = "id")
  data <- transform(data, ftime = factor(time))
  # summary(data)

  # Constructs the prior penalty data frame if a prior is provided.
  if (!is.null(prior)) {
    penalty_df <- data.frame(prior = prior, class = "fixef", coef = "")
  } else {
    penalty_df <- NULL
  }

  # The model formula with fixed and random effects.
  f <- paste0("n ~ 0 + species + species:(", paste(paste0("X", seq(ncol(env))), collapse = " + "), ")")
  f <- paste(f, "+ (0+time|species) + (0+time|species:id) + (1|ftime:species)")
  f <- formula(f)

  # Delete the ALTREP object (and all intermediates) from the worker's
  # local environment right before running the model.
  # glmmTMB can no longer capture the mori pointer because it is gone.
  rm(obj, ve, env)

  # Fits the generalized linear mixed model.
  # m <-  glm(f, data = data, family = poisson, control = glm.control(maxit = 100))
  # m <- glmmTMB(formula = f, data = data, family = poisson, ziformula = ~ (1|species))
  m <- glmmTMB(formula = f, data = data, family = poisson, priors = penalty_df)

  # summary(m)
  return(m)
}


compare_niches <- function(obj, m, type = c("beta", "spatial"), measure = cor) {

  # Compares fundamental niches (simulation parameters) against realized niches (model estimates).
  #
  # Arguments:
  #   obj: A 'sim_com_results' object containing simulation truth.
  #        - For type="beta": Uses obj$beta (environmental slopes).
  #        - For type="spatial": Uses obj$K (data.frame of carrying capacities).
  #     m: A fitted 'glmmTMB' model object.
  #  type: Character string.
  #        - "beta": Compares true environmental slopes against estimated slopes.
  #        - "spatial": Compares log-transformed K against model predictions (link scale).
  # measure: Function. The metric used to calculate similarity (default: stats::cor).
  #
  # Returns:
  #   A named numeric vector containing the similarity score for each species.
  #
  # Details:
  #   The function delegates to internal helpers based on 'type'.
  #   - 'beta' assumes model coefficients follow a ":X" pattern for slopes.
  #   - 'spatial' assumes obj$K rows correspond to the 'id' factor in the model frame.
  #   Missing or non-finite values in coefficients or predictions are excluded pairwise.

  # Handle case where model fitting failed (m is NA)
  if (length(m) == 1 && is.na(m)) {
    spec_names <- obj$spec_names
    res <- rep(NA_real_, length(spec_names))
    names(res) <- spec_names
    return(res)
  }

  type <- match.arg(type)

  if (type == "beta") {
    return(.compare_beta(obj, m, measure))
  } else {
    return(.compare_spatial(obj, m, measure))
  }
}

.compare_beta <- function(obj, m, measure) {

  # Internal helper for beta comparison.

  spec_names <- obj$spec_names

  # 1. Fundamental niche coefficients (remove intercept)
  true_slopes <- obj$beta[, -1, drop = FALSE]

  # 2. Estimated niche coefficients (extract interaction terms :X)
  cf <- glmmTMB::fixef(m)$cond
  cfx <- cf[grep(":X", names(cf))]

  if (length(cfx) != length(true_slopes)) {
    stop("Dimension mismatch: Model coefficients do not match simulation parameters.")
  }

  estimated_slopes <- matrix(cfx,
                             nrow = nrow(true_slopes),
                             ncol = ncol(true_slopes))

  # 3. Calculate similarity
  similarity <- mapply(measure,
                       asplit(estimated_slopes, 1),
                       asplit(true_slopes, 1))
  names(similarity) <- spec_names

  return(similarity)
}

.compare_spatial <- function(obj, m, measure) {

  # Internal helper for spatial comparison.

  spec_names <- obj$spec_names

  # 1. Fundamental niche (Reference)
  # K is a data.frame: columns = species, rows = spatial ID
  K <- obj$K
  log_K <- log(K)

  # Create an explicit ID column for merging.
  # Assumes rows in K correspond sequentially to IDs (1 to N).
  log_K$id <- factor(seq_len(nrow(K)))

  # 2. Realised niche (Predictions)
  data <- m$frame
  data$pre <- predict(m, type = "link")

  # Merge wide reference data (log_K) into long model data
  data <- merge(data, log_K, by = "id")

  # 3. Calculate similarity
  # Split by species for faster processing than repeated subset()
  data_split <- split(data, data$species)

  similarity <- vapply(spec_names, function(i) {
    dat <- data_split[[i]]

    if (is.null(dat)) return(NA_real_)

    # Compare predictions vs the specific species column from log_K
    measure(dat[["pre"]], dat[[i]])
  }, numeric(1))

  return(similarity)
}
