
# -------------------------------------------------------------------------
# Figures
# -------------------------------------------------------------------------


# Setup -------------------------------------------------------------------

library(RColorBrewer)

pal <- brewer.pal(8, "Dark2")
col <- adjustcolor(pal, alpha.f = 0.5)
ribbon_col <- adjustcolor(col, alpha.f = 0.3)
point_col  <- adjustcolor(col, alpha.f = 0.02)

plot_model <- function(obj, trans = identity, ref = c(0, 1), col = "black", ...) {

  model_data <- obj$model
  log_props <- sort(unique(model_data$log_prop))
  props <- exp(log_props)
  n_sites <- props * 1e4

  linkinv <- obj$family$linkinv
  se_trial <- gam.vcomp(obj)$vc[["s(trial)"]]

  op <- par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))

  # (a) Interaction strength ----
  x <- model_data$inter
  x_grid <- seq(min(x), max(x), length.out = 100)
  newdata <- data.frame(inter = x_grid, log_prop = median(model_data$log_prop), trial = model_data$trial[1])
  pre <- predict.gam(obj, newdata = newdata, exclude = "s(trial)", se.fit = TRUE)
  fit <- pre[["fit"]]
  se <- pre$se.fit
  se_rand <- sqrt(se^2 + se_trial^2)
  lci <- fit - 1.96 * se; uci <- fit + 1.96 * se
  lci_rand <- fit - 1.96 * se_rand; uci_rand <- fit + 1.96 * se_rand
  pred_mat <- cbind(fit = fit, lci = lci, uci = uci, lci_rand = lci_rand, uci_rand = uci_rand)
  pred_mat <- trans(linkinv(pred_mat))
  nd <- data.frame(newdata, pred_mat)

  plot(jitter(trans(model_data[[1]])) ~ inter, model_data, pch = 20, cex = 0.8, col = point_col, xlab = "Interaction strength", ...)
  polygon(x = c(nd$inter, rev(nd$inter)), y = c(nd$lci_rand, rev(nd$uci_rand)), border = NA, col = ribbon_col)
  # polygon(x = c(nd$inter, rev(nd$inter)), y = c(nd$lci, rev(nd$uci)), border = NA, col = ribbon_col)
  lines(nd$inter, nd$fit, lwd = 2, col = col)
  if(!is.null(ref))
    abline(h = ref, lty = 3, col = "grey30")

  mtext("(a)", side = 3, line = 1, adj = -0.2, cex = 1.5, font = 1)


  # (b) Sample size ----
  x <- model_data$log_prop
  x_grid <- seq(min(x), max(x), length.out = 100)
  newdata <- data.frame(inter = median(model_data$inter), log_prop = x_grid, trial = model_data$trial[1])
  pre <- predict.gam(obj, newdata = newdata, exclude = "s(trial)", se.fit = TRUE)
  fit <- pre[["fit"]]
  se <- pre$se.fit
  se_rand <- sqrt(se^2 + se_trial^2)
  lci <- fit - 1.96 * se; uci <- fit + 1.96 * se
  lci_rand <- fit - 1.96 * se_rand; uci_rand <- fit + 1.96 * se_rand
  pred_mat <- cbind(fit = fit, lci = lci, uci = uci, lci_rand = lci_rand, uci_rand = uci_rand)
  pred_mat <- trans(linkinv(pred_mat))
  nd <- data.frame(newdata, pred_mat)

  plot(jitter(trans(model_data[[1]])) ~ jitter(log_prop), model_data, pch = 20, cex = 0.8, col = point_col, xlab = "Sample size", xaxt = "n", ...)
  axis(1, at = log(props), labels = n_sites)
  polygon(x = c(nd$log_prop, rev(nd$log_prop)), y = c(nd$lci_rand, rev(nd$uci_rand)), border = NA, col = ribbon_col)
  # polygon(x = c(nd$log_prop, rev(nd$log_prop)), y = c(nd$lci, rev(nd$uci)), border = NA, col = ribbon_col)
  lines(nd$log_prop, nd$fit, lwd = 2, col = col)
  if(!is.null(ref))
    abline(h = ref, lty = 3, col = "grey30")

  mtext("(b)", side = 3, line = 1, adj = -0.2, cex = 1.5, font = 1)

  par(op)
}


# Data & fitted GAMMs -----------------------------------------------------

load("results/data.RData")
load("results/GAMM_models.RData")


# Figure 1 ----------------------------------------------------------------

op <- par(mfrow = c(1, 3), mar = c(5, 3, 4, 1), oma = c(0, 3, 0, 0))

# Magnitude of bias
hist(data$magnitude, col = col[1], main = "", xlab = expression(paste("Absolute bias of ", beta)), ylab = "", cex.lab = 2, cex.axis = 1.3)
mtext("(a)", side = 3, line = 2, adj = -0.2, cex = 1.5)

# Sign agreement
hist(data$sign_agreement, col = col[1], main = "", xlab = "Sign agreement", ylab = "", cex.lab = 2, cex.axis = 1.3)
mtext("(b)", side = 3, line = 2, adj = -0.2, cex = 1.5, font = 1)

# Correlation with K
hist(data$spat_pred, col = col[1], main = "", xlab = "Correlation with K", ylab = "", cex.lab = 2, cex.axis = 1.3)
mtext("(c)", side = 3, line = 2, adj = -0.2, cex = 1.5, font = 1)

mtext("Frequency", side = 2, line = 1, adj = 0.55, outer = TRUE, cex = 1.5)

par(op)



# Fig. 1 ----
plot_model(ma, trans = exp, ylab = expression(paste("Absolute bias of ", beta)), col = pal[1], cex.lab = 1.5, ylim = c(0, 1.1), ref = 0)

# Fig. 2 ----
plot_model(ms, ylab = "Sign agreement", col = pal[1], cex.lab = 1.5, ylim = c(0.45, 1))

# Fig. 3 ----
plot_model(mr, trans = tanh, ylab = "Correlation with K", col = pal[1], cex.lab = 1.5)

