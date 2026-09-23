
# -------------------------------------------------------------------------
# Analyse: fitting GAM to discrepancy measures
# -------------------------------------------------------------------------

# Setup -------------------------------------------------------------------

# Libraries ----
library(mgcv)
library(car)

ctrl <- mgcv::gam.control(nthreads = 8L)

# Helper functions ----
source("R/utils.R")

var_names <- c("inter", "magnitude", "sign_agreement", "spat_pred")

sigmas <- seq(0, 3, 0.25)


# Data --------------------------------------------------------------------

# Load discrepancy measures
var_list <- lapply(var_names, read_rds_array, path = "results", simplify = FALSE)
names(var_list) <- var_names
list2env(var_list, .GlobalEnv)

# Convert list of matrices to a long data frames
inter_long <- reshape_2_long(inter, "inter")
magnitude_long <- reshape_2_long(magnitude, "magnitude")
sign_agreement_long <- reshape_2_long(sign_agreement, "sign_agreement")
spat_pred_long <- reshape_2_long(spat_pred, "spat_pred")

data <- cbind(inter_long, magnitude_long["magnitude"], sign_agreement_long["sign_agreement"], spat_pred_long["spat_pred"])

rm(inter, inter_long, magnitude, magnitude_long, sign_agreement, sign_agreement_long, spat_pred, spat_pred_long, var_list)

data <- transform(data, trial = factor(trial), log_mag = log(magnitude), log_prop = log(prop), logit_sign = car::logit(sign_agreement))
summary(data)

save(data, file = "results/data.RData")


# GAMM --------------------------------------------------------------------

length(unique(data$prop))

# Magnitude ----

summary(data$magnitude); hist(data$magnitude, breaks = 50)
hist(data$log_mag, breaks = 50)

fm <- log_mag ~ s(inter, bs = "ad", k = 20) + s(log_prop, bs = "cr", k = 8) + s(trial, bs = "re")
ma <- gam(fm, data = data, control = ctrl, gamma = 5)

summary(ma)
plot.gam(ma, pages = 1, all.terms = TRUE, residuals = TRUE, scale = 0, seWithMean = TRUE)
plot.gam(ma, pages = 1, all.terms = TRUE, residuals = FALSE, rug = TRUE, scale = 0, seWithMean = TRUE, trans = exp)


# Sign agreement (probability) ----

summary(data$sign_agreement); hist(data$sign_agreement, breaks = 50)

# 5 species * 4 environmental variables
weights <- rep(20, nrow(data))

fs <- sign_agreement ~ s(inter, bs = "ad", k = 20) + s(log_prop, bs = "cr", k = 8) + s(trial, bs = "re")
ms <- gam(fs, family = binomial, weights = weights, data = data, control = ctrl, gamma = 5)

summary(ms)
plot.gam(ms, pages = 1, all.terms = TRUE, residuals = TRUE, scale = 0, seWithMean = TRUE)
plot.gam(ms, pages = 1, all.terms = TRUE, residuals = FALSE, scale = 0, seWithMean = TRUE, trans = plogis)


# Spatial prediction (correlation) ----

summary(data$spat_pred); hist(data$spat_pred, breaks = 50)
# Strictly positive, no edge cases.

data$z <- atanh(data$spat_pred)

fr <- z ~ s(inter, bs = "ad", k = 20) + s(log_prop, bs = "cr", k = 8) + s(trial, bs = "re")
mr <- gam(fr, data = data, control = ctrl, gamma = 5)

summary(mr)
plot.gam(mr, pages = 1, all.terms = TRUE, residuals = TRUE, scale = 0, seWithMean = TRUE)
plot.gam(mr, pages = 1, all.terms = TRUE, residuals = FALSE, scale = 0, seWithMean = TRUE, trans = tanh)


# Save models -------------------------------------------------------------

save(ma, ms, mr, file = "results/GAMM_models.RData")
