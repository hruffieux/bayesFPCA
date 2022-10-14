plot_fpca_data <- function(N_sample, time_obs, Y, plot_dim, data_col) {

  Y_sample <- Y[N_sample]
  time_obs_sample <- time_obs[N_sample]
  T_sample <- sapply(Y_sample, length)
  length_N <- length(N_sample)

  Y_vec <- Reduce(c, Y_sample)
  time_vec <- Reduce(c, time_obs_sample)

  curve_labels <- vector("list", length=length_N)
  curve_id <- rep(NA, length_N)
  for(i in 1:length_N) {

    N_i <- N_sample[i]
    curve_id[i] <- parse(text=paste("y[", N_i, "] (t)", sep=""))
    curve_val <- eval(bquote(expression(y[.(N_i)] (t))))
    curve_labels[[i]] <- rep(curve_val, T_sample[i])
  }
  curve_labels <- do.call(c, curve_labels)
  curve_labels <- factor(curve_labels, levels=curve_id)

  strip.math <- function(
    which.given, which.panel, var.name, factor.levels, ...
  ) {

    fl <- curve_id

    strip.default(which.given,which.panel,var.name,fl,...)
  }

  raw_data_plots <- xyplot(
    Y_vec ~ time_vec | curve_labels, groups=curve_labels,
    data=data.frame(
      time_vec=time_vec, Y_vec=Y_vec,
      curve_labels=curve_labels
    ),
    layout=plot_dim, main="",
    strip=strip.math,
    par.strip.text=list(cex=0.8),
    par.settings = list(layout.heights = list(strip = 1)),
    xlab="time",
    ylab="response curves",
    as.table=TRUE,
    panel=function(x, y, subscripts, groups) {

      iPan <- panel.number()
      i <- rep(N_sample, each=1)[iPan]
      panel.grid()
      panel.superpose(
        x[order(x)], y[order(x)], subscripts, groups,
        type="p", col=data_col, pch=16, cex=0.4
      )
    }
  )

  print(raw_data_plots)
}

plot_gauss_mfpca_data <- function(N_sample, time_obs, Y, plot_dim, data_col) {

  p <- length(Y[[1]])
  n <- Reduce(rbind, lapply(Y, function(x) sapply(x, length)))
  n_sample <- length(N_sample)

  Y_vec <- unlist(Y[N_sample])
  time_vec <- unlist(time_obs[N_sample])

  curve_labels_1 <- vector("list", length = n_sample)
  for(i in 1:n_sample) {

    N_i <- N_sample[i]
    curve_labels_1[[i]] <- rep.int(1:p, times = n[N_i, ])
  }
  curve_labels_1 <- Reduce(c, curve_labels_1)
  curve_id_1 <- rep(NA, p)
  for(j in 1:p) {

    curve_id_1[j] <- parse(text = as.character(j))
  }
  curve_labels_1 <- factor(curve_labels_1, levels = curve_id_1)

  curve_labels_2 <- vector("list", length = n_sample)
  curve_id_2 <- rep(NA, n_sample)
  for(i in 1:n_sample) {

    N_i <- N_sample[i]
    curve_id_2[i] <- parse(text=paste("Y [", N_i, "] (t)", sep=""))
    curve_val <- eval(bquote(expression(Y[.(N_i)] (t))))
    curve_labels_2[[i]] <- rep(curve_val, sum(n[N_i, ]))
  }
  curve_labels_2 <- do.call(c, curve_labels_2)
  curve_labels_2 <- factor(curve_labels_2, levels = curve_id_2)

  strip.math <- function(
    which.given, which.panel, var.name, factor.levels, ...
  ) {

    if(which.given==1) {

      fl <- curve_id_1

      strip.default(which.given, which.panel, var.name, fl, ...)
    }

    if (which.given==2) {

      fl <- curve_id_2

      strip.default(which.given, which.panel, var.name, fl, ...)
    }
  }

  raw_data_plots <- xyplot(
    Y_vec ~ time_vec | curve_labels_1*curve_labels_2, groups = curve_labels_1,
    data = data.frame(
      time_vec = time_vec, Y_vec = Y_vec,
      curve_labels_1 = curve_labels_1, curve_labels_2 = curve_labels_2
    ), layout = c(p, n_sample), main = "",
    strip=strip.math,
    par.strip.text=list(cex=0.8),
    par.settings = list(layout.heights = list(strip = 1)),
    xlab = "time", ylab = "functional responses", as.table = TRUE,
    panel=function(x,y,subscripts,groups) {

      panel.grid()
      panel.superpose(
        x[order(x)], y[order(x)], subscripts, groups,
        type="p", col = data_col, pch = 16, cex = 0.4
      )
    }
  )

  print(raw_data_plots)
}

plot_fpca_fits <- function(
    N_sample, time_obs, Y,
    time_g, Y_summary,
    plot_dim, data_col, model_col,
    ylab = NULL
) {

  Y_sample <- Y[N_sample]
  time_obs_sample <- time_obs[N_sample]
  T_sample <- sapply(Y_sample, length)
  length_N <- length(N_sample)

  Y_vec <- Reduce(c, Y_sample)
  time_vec <- Reduce(c, time_obs_sample)

  curve_labels <- vector("list", length=length_N)
  curve_id <- rep(NA, length_N)
  for(i in 1:length_N) {

    N_i <- N_sample[i]
    curve_id[i] <- parse(text=paste("y[", N_i, "] (t)", sep=""))
    curve_val <- eval(bquote(expression(y[.(N_i)] (t))))
    curve_labels[[i]] <- rep(curve_val, T_sample[i])
  }
  curve_labels <- do.call(c, curve_labels)
  curve_labels <- factor(curve_labels, levels=curve_id)

  strip.math <- function(
    which.given, which.panel, var.name, factor.levels, ...
  ) {

    fl <- curve_id

    strip.default(which.given,which.panel,var.name,fl,...)
  }

  fitted_data_plots <- xyplot(
    Y_vec ~ time_vec | curve_labels, groups=curve_labels,
    data=data.frame(
      time_vec=time_vec, Y_vec=Y_vec,
      curve_labels=curve_labels
    ),
    layout=plot_dim, main="",
    strip=strip.math,
    par.strip.text=list(cex=0.8),
    par.settings = list(layout.heights = list(strip = 1)),
    xlab="time",
    ylab=ifelse(is.null(ylab), "response curves", ylab),
    as.table=TRUE,
    panel=function(x, y, subscripts, groups) {

      iPan <- panel.number()
      i <- rep(N_sample, each=1)[iPan]
      panel.grid()
      panel.superpose(
        x[order(x)], y[order(x)], subscripts, groups,
        type="p", col=data_col, pch=16, cex=0.4
      )
      panel.xyplot(
        time_g, Y_summary[[i]][,2],
        col=model_col, type="l", lwd=1
      )
      panel.xyplot(
        time_g, Y_summary[[i]][,1],
        col=model_col, type="l", lwd=1, lty=2
      )
      panel.xyplot(
        time_g, Y_summary[[i]][,3],
        col=model_col, type="l", lwd=1, lty=2
      )
    }
  )

  print(fitted_data_plots)
}

plot_fit_comparisons <- function(
    N_sample, time_obs, time_g,
    Y, Y_vmp_summary, Y_mcmc_summary,
    plot_dim, vmp_col, mcmc_col, data_col
) {

  Y_sample <- Y[N_sample]
  time_obs_sample <- time_obs[N_sample]
  T_sample <- sapply(Y_sample, length)
  length_N <- length(N_sample)

  Y_vec <- Reduce(c, Y_sample)
  time_vec <- Reduce(c, time_obs_sample)

  curve_labels <- vector("list", length=length_N)
  curve_id <- rep(NA, length_N)
  for(i in 1:length_N) {

    N_i <- N_sample[i]
    curve_id[i] <- parse(text=paste("y[", N_i, "] (t)", sep=""))
    curve_val <- eval(bquote(expression(y[.(N_i)] (t))))
    curve_labels[[i]] <- rep(curve_val, T_sample[i])
  }
  curve_labels <- do.call(c, curve_labels)
  curve_labels <- factor(curve_labels, levels=curve_id)

  strip.math <- function(
    which.given, which.panel, var.name, factor.levels, ...
  ) {

    fl <- curve_id

    strip.default(which.given,which.panel,var.name,fl,...)
  }

  fit_comparisons <- xyplot(
    Y_vec ~ time_vec | curve_labels, groups=curve_labels,
    data=data.frame(
      time_vec=time_vec, Y_vec=Y_vec,
      curve_labels=curve_labels
    ),
    layout=plot_dim, main="",
    strip=strip.math,
    par.strip.text=list(cex=0.8),
    par.settings = list(layout.heights = list(strip = 1.2)),
    xlab="time",
    ylab="nonlinear curves",
    as.table=TRUE,
    panel=function(x, y, subscripts, groups) {

      iPan <- panel.number()
      i <- rep(N_sample, each=1)[iPan]
      panel.grid()
      panel.superpose(
        x[order(x)], y[order(x)], subscripts, groups,
        type="p", col=data_col, pch=16, cex=0.4
      )
      panel.xyplot(
        time_g, Y_mcmc_summary[[i]][,1],
        col=mcmc_col, type="l", lwd=mcmc_lwd, lty=2
      )
      panel.xyplot(
        time_g, Y_mcmc_summary[[i]][,2],
        col=mcmc_col, type="l", lwd=mcmc_lwd
      )
      panel.xyplot(
        time_g, Y_mcmc_summary[[i]][,3],
        col=mcmc_col, type="l", lwd=mcmc_lwd, lty=2
      )
      panel.xyplot(
        time_g, Y_vmp_summary[[i]][,1],
        col=vmp_col, type="l", lwd=vmp_lwd, lty=2
      )
      panel.xyplot(
        time_g, Y_vmp_summary[[i]][,2],
        col= vmp_col, type="l", lwd=vmp_lwd
      )
      panel.xyplot(
        time_g, Y_vmp_summary[[i]][,3],
        col= vmp_col, type="l", lwd=vmp_lwd, lty=2
      )
    }
  )

  print(fit_comparisons)
}

plot_mfpca_fit_comparisons <- function(
    N_sample, time_obs, time_g,
    Y, Y_vmp_summary, Y_alt_mod_summary,
    plot_dim, vmp_col, alt_mod_col, data_col
) {

  p <- length(Y[[1]])
  n <- Reduce(rbind, lapply(Y, function(x) sapply(x, length)))
  n_sample <- length(N_sample)
  Y_vec <- unlist(Y[N_sample])
  time_vec <- unlist(time_obs[N_sample])

  curve_labels_1 <- vector("list", length = n_sample)
  for(i in 1:n_sample) {

    N_i <- N_sample[i]
    curve_labels_1[[i]] <- rep.int(1:p, times = n[N_i, ])
  }
  curve_labels_1 <- Reduce(c, curve_labels_1)
  curve_id_1 <- rep(NA, p)
  for(j in 1:p) {

    curve_id_1[j] <- parse(text = as.character(j))
  }
  curve_labels_1 <- factor(curve_labels_1, levels = curve_id_1)

  curve_labels_2 <- vector("list", length = n_sample)
  curve_id_2 <- rep(NA, n_sample)
  for(i in 1:n_sample) {

    N_i <- N_sample[i]
    curve_id_2[i] <- parse(text=paste("Y [", N_i, "] (t)", sep=""))
    curve_val <- eval(bquote(expression(Y[.(N_i)] (t))))
    curve_labels_2[[i]] <- rep(curve_val, sum(n[N_i, ]))
  }
  curve_labels_2 <- do.call(c, curve_labels_2)
  curve_labels_2 <- factor(curve_labels_2, levels = curve_id_2)

  strip.math <- function(
    which.given, which.panel, var.name, factor.levels, ...
  ) {

    if(which.given==1) {

      fl <- curve_id_1

      strip.default(which.given, which.panel, var.name, fl, ...)
    }

    if (which.given==2) {

      fl <- curve_id_2

      strip.default(which.given, which.panel, var.name, fl, ...)
    }
  }

  fitted_data_plots <- xyplot(
    Y_vec ~ time_vec | curve_labels_1*curve_labels_2, groups = curve_labels_1,
    data = data.frame(
      time_vec = time_vec, Y_vec = Y_vec,
      curve_labels_1 = curve_labels_1, curve_labels_2 = curve_labels_2
    ), layout = c(p, n_sample), main = "",
    strip=strip.math,
    par.strip.text=list(cex=0.8),
    par.settings = list(layout.heights = list(strip = 1)),
    xlab = "time", ylab = "functional responses", as.table = TRUE,
    panel=function(x,y,subscripts,groups) {

      i_pan <- panel.number()
      i <- rep(1:n_sample, each = p)[i_pan]
      j <- rep(1:p, n_sample)[i_pan]
      panel.grid()

      panel.superpose(
        x[order(x)], y[order(x)], subscripts, groups,
        type="p", col = data_col, pch = 16, cex = 0.4
      )

      panel.xyplot(
        time_g, Y_alt_mod_summary[[N_sample[i]]][[j]][, 1], col = alt_mod_col,
        type = "l", lwd = 2, lty = 2
      )

      panel.xyplot(
        time_g, Y_alt_mod_summary[[N_sample[i]]][[j]][, 2], col = alt_mod_col,
        type = "l", lwd = 2
      )

      panel.xyplot(
        time_g, Y_alt_mod_summary[[N_sample[i]]][[j]][, 3], col = alt_mod_col,
        type = "l", lwd = 2, lty = 2
      )

      panel.xyplot(
        time_g, Y_vmp_summary[[N_sample[i]]][[j]][, 1], col = vmp_col,
        type = "l", lwd = 1, lty = 2
      )

      panel.xyplot(
        time_g, Y_vmp_summary[[N_sample[i]]][[j]][, 2], col = vmp_col,
        type = "l", lwd = 1
      )

      panel.xyplot(
        time_g, Y_vmp_summary[[N_sample[i]]][[j]][, 3], col = vmp_col,
        type = "l", lwd = 1, lty = 2
      )
    }
  )

  print(fitted_data_plots)
}

plot_fpca_scores <- function(N_sample, zeta_summary, zeta, plot_dim, data_col, model_col) {

  length_N <- length(N_sample)

  Zeta_hat <- matrix(NA, length_N, 2)
  zeta_ellipse <- vector("list", length=length_N)
  zeta_labels <- vector("list", length=length_N)
  zeta_id <- rep(NA, length_N)
  for(i in 1:length_N) {

    N_i <- N_sample[i]

    Zeta_hat[i,] <- zeta_summary[[N_i]][[1]]
    zeta_ellipse[[i]] <- zeta_summary[[N_i]][[2]]

    zeta_id[i] <- parse(text=paste("zeta[", N_i, "]", sep=""))
    zeta_val <- eval(bquote(expression(zeta[.(N_i)])))
    zeta_labels[[i]] <- rep(zeta_val, nrow(zeta_ellipse[[i]]))
  }
  zeta_labels <- do.call(c, zeta_labels)
  zeta_labels <- factor(zeta_labels, levels=zeta_id)

  strip.math <- function(
    which.given, which.panel, var.name, factor.levels, ...
  ) {

    fl <- zeta_id

    strip.default(which.given,which.panel,var.name,fl,...)
  }

  zeta_ellipse_mat <- Reduce(rbind, zeta_ellipse)
  zeta_ellipse_x <- zeta_ellipse_mat[,1]
  zeta_ellipse_y <- zeta_ellipse_mat[,2]

  score_plots <- xyplot(
    zeta_ellipse_y ~ zeta_ellipse_x | zeta_labels, groups=zeta_labels,
    data=data.frame(
      zeta_ellipse_x=zeta_ellipse_x, zeta_ellipse_y=zeta_ellipse_y,
      zeta_labels=zeta_labels
    ),
    layout=plot_dim, main="",
    strip=strip.math,
    par.strip.text=list(cex=0.8),
    par.settings = list(layout.heights = list(strip = 1)),
    xlab="first basis function's score",
    ylab="second basis function's score",
    as.table=TRUE,
    panel=function(x, y, subscripts, groups) {

      iPan <- panel.number()
      i <- rep(1:length_N, each=1)[iPan]
      panel.grid()
      panel.xyplot(
        zeta_ellipse[[i]][,1], zeta_ellipse[[i]][,2],
        col=model_col, type="l", lwd=1
      )
      panel.xyplot(
        Zeta_hat[i,1], Zeta_hat[i,2],
        col=model_col, type="p", pch=16, cex=0.4
      )
      panel.xyplot(
        zeta[[N_sample[i]]][1], zeta[[N_sample[i]]][2],
        col=data_col, type="p", pch=16, cex=0.4
      )
    }
  )

  print(score_plots)
}

plot_score_comparisons <- function(
    N_sample, zeta, zeta_vmp_summary, zeta_mcmc_summary,
    data_col, vmp_col, mcmc_col, plot_dim,
    vmp_lwd = 1, mcmc_lwd = 2
) {

  length_N <- length(N_sample)

  zeta_labels <- vector("list", length=length_N)
  zeta_id <- rep(NA, length_N)
  Zeta_vmp_hat <- matrix(NA, length_N, 2)
  Zeta_mcmc_hat <- matrix(NA, length_N, 2)
  zeta_vmp_ellipse <- vector("list", length=length_N)
  zeta_mcmc_ellipse <- vector("list", length=length_N)
  for(i in 1:length_N) {

    N_i <- N_sample[i]

    Zeta_vmp_hat[i,] <- zeta_vmp_summary[[N_i]][[1]]
    zeta_vmp_ellipse[[i]] <- zeta_vmp_summary[[N_i]][[2]]

    Zeta_mcmc_hat[i,] <- zeta_mcmc_summary[[N_i]][[1]]
    zeta_mcmc_ellipse[[i]] <- zeta_mcmc_summary[[N_i]][[2]]

    zeta_id[i] <- parse(text=paste("zeta[", N_i, "]", sep=""))
    zeta_val <- eval(bquote(expression(zeta[.(N_i)])))
    zeta_labels[[i]] <- rep(zeta_val, nrow(zeta_mcmc_ellipse[[i]]))
  }
  zeta_labels <- do.call(c, zeta_labels)
  zeta_labels <- factor(zeta_labels, levels=zeta_id)

  strip.math <- function(
    which.given, which.panel, var.name, factor.levels, ...
  ) {

    fl <- zeta_id

    strip.default(which.given,which.panel,var.name,fl,...)
  }

  zeta_mcmc_ellipse_mat <- Reduce(rbind, zeta_mcmc_ellipse)
  zeta_mcmc_ellipse_x <- zeta_mcmc_ellipse_mat[,1]
  zeta_mcmc_ellipse_y <- zeta_mcmc_ellipse_mat[,2]

  score_comparisons <- xyplot(
    zeta_mcmc_ellipse_y ~ zeta_mcmc_ellipse_x | zeta_labels, groups=zeta_labels,
    data=data.frame(
      zeta_mcmc_ellipse_x=zeta_mcmc_ellipse_x, zeta_mcmc_ellipse_y=zeta_mcmc_ellipse_y,
      zeta_labels=zeta_labels
    ),
    layout=plot_dim, main="",
    strip=strip.math,
    par.strip.text=list(cex=0.8),
    par.settings = list(layout.heights = list(strip = 1)),
    xlab="first basis function's score",
    ylab="second basis function's score",
    as.table=TRUE,
    panel=function(x, y, subscripts, groups) {

      iPan <- panel.number()
      i <- rep(1:length_N, each=1)[iPan]
      panel.grid()
      panel.xyplot(
        zeta_mcmc_ellipse[[i]][,1], zeta_mcmc_ellipse[[i]][,2],
        col=mcmc_col, type="l", lwd=mcmc_lwd
      )
      panel.xyplot(
        Zeta_mcmc_hat[i,1], Zeta_mcmc_hat[i,2],
        col=mcmc_col, type="p", pch=16, cex=0.6
      )
      panel.xyplot(
        zeta_vmp_ellipse[[i]][,1], zeta_vmp_ellipse[[i]][,2],
        col=vmp_col, type="l", lwd=vmp_lwd
      )
      panel.xyplot(
        Zeta_vmp_hat[i,1], Zeta_vmp_hat[i,2],
        col=vmp_col, type="p", pch=16, cex=0.6
      )
      panel.xyplot(
        zeta[[N_sample[i]]][1], zeta[[N_sample[i]]][2],
        col=data_col, type="p", pch=16, cex=0.6
      )
    }
  )

  print(score_comparisons)
}

plot_fpca_global_curves <- function(
    gbl_estimates, L_true, time_g, mu_g, Psi_g,
    model_col, data_col, plot_gbl_dim, ylab_add = ""
) {

  n_g <- length(time_g)

  if (!is.null(mu_g)) {
    nb_col <- L_true+1
  } else { # if mu_g = NULL, mean function not displayed
    nb_col <- L_true
  }

  gbl_estimates <- gbl_estimates[,1:nb_col]

  time_g_gbl <- rep(time_g, nb_col)
  if (!is.null(mu_g)) {
    gbl_g_vec <- c(mu_g, as.vector(Psi_g[,1:L_true]))
  } else {
    gbl_g_vec <- as.vector(Psi_g[,1:L_true])
  }

  gbl_labels <- vector("list", length= nb_col)
  gbl_id <- rep(NA,  nb_col)

  if (!is.null(mu_g)) {
    gbl_id[1] <- expression(mu (t))
    gbl_labels[[1]] <- rep(gbl_id[1], n_g)
  }

  for(l in 1:L_true) {

    if (!is.null(mu_g)) {
      ll <- l+1
    } else {
      ll <- l
    }

    gbl_id[ll] <- parse(text=paste("psi[", l, "] (t)", sep=""))
    bf_val <- eval(bquote(expression(psi[.(l)] (t))))
    gbl_labels[[ll]] <- rep(bf_val, n_g)
  }

  gbl_labels <- do.call(c, gbl_labels)
  gbl_labels <- factor(gbl_labels, levels=gbl_id)

  strip.math <- function(
    which.given, which.panel, var.name, factor.levels, ...
  ) {

    fl <- gbl_id

    strip.default(which.given,which.panel,var.name,fl,...)
  }

  gbl_plots <- xyplot(
    gbl_g_vec ~ time_g_gbl | gbl_labels, groups=gbl_labels,
    data=data.frame(
      time_g_gbl=time_g_gbl, gbl_g_vec=gbl_g_vec,
      gbl_labels=gbl_labels
    ),
    layout=plot_gbl_dim, main="",
    strip=strip.math,
    par.strip.text=list(cex=0.8),
    par.settings = list(layout.heights = list(strip = 1)),
    xlab="time",
    ylab=paste(ifelse(!is.null(mu_g), "mean and basis functions", "basis functions"),
               ylab_add),
    as.table=TRUE,
    panel=function(x, y, subscripts, groups) {

      lPan <- panel.number()
      l <- rep(1:nb_col, each=1)[lPan]
      panel.grid()
      panel.superpose(
        x[order(x)], y[order(x)], subscripts, groups,
        type="l", col=data_col, lwd=2
      )
      panel.xyplot(
        time_g, gbl_estimates[,l],
        col=model_col, type="l", lwd=1
      )
    }
  )

  print(gbl_plots)
}

plot_mfpca_global_curves <- function(
    vmp_gbl_est, alt_mod_gbl_est, time_g, mu_g, Psi_g,
    vmp_col, alt_mod_col, data_col
) {

  n_g <- length(time_g)
  p <- ncol(vmp_gbl_est[[1]])
  L <- length(vmp_gbl_est) - 1

  time_g_gbl <- rep(time_g, p*(L + 1))
  gbl_g_vec <- Reduce(c, lapply(mapply(cbind, mu_g, Psi_g, SIMPLIFY = FALSE), as.vector))

  curve_labels_1 <- rep(rep(1:p, each = n_g), L + 1)
  curve_id_1 <- rep(NA, p)
  for(j in 1:p) {

    curve_id_1[j] <- parse(text = as.character(j))
  }
  curve_labels_1 <- factor(curve_labels_1, levels=curve_id_1)

  curve_labels_2 <- vector("list", length = L + 1)
  curve_id_2 <- rep(NA, L + 1)
  curve_id_2[1] <- expression(mu (t))
  curve_labels_2[[1]] <- rep(curve_id_2[1], n_g*p)
  for(l in 1:L) {

    curve_id_2[l+1] <- parse(text=paste("psi[", l, "] (t)", sep=""))
    curve_val <- eval(bquote(expression(psi[.(l)] (t))))
    curve_labels_2[[l + 1]] <- rep(curve_val, n_g*p)
  }
  curve_labels_2 <- do.call(c, curve_labels_2)
  curve_labels_2 <- factor(curve_labels_2, levels=curve_id_2)

  strip.math <- function(
    which.given, which.panel, var.name, factor.levels, ...
  ) {

    if(which.given==1) {

      fl <- curve_id_1

      strip.default(which.given, which.panel, var.name, fl, ...)
    }

    if (which.given==2) {

      fl <- curve_id_2

      strip.default(which.given, which.panel, var.name, fl, ...)
    }
  }

  fitted_gbl_plots <- xyplot(
    gbl_g_vec ~ time_g_gbl | curve_labels_2*curve_labels_1, groups = curve_labels_2,
    data = data.frame(
      time_g_gbl = time_g_gbl, gbl_g_vec = gbl_g_vec,
      curve_labels_2 = curve_labels_2, curve_labels_1 = curve_labels_1
    ), layout = c(p, L + 1), main = "",
    strip=strip.math,
    par.strip.text=list(cex=0.8),
    par.settings = list(layout.heights = list(strip = 1)),
    xlab = "time", ylab = "mean and eigenfunctions", as.table = TRUE,
    panel=function(x,y,subscripts,groups) {

      l_pan <- panel.number()
      l <- rep(1:(L+1), each = p)[l_pan]
      j <- rep(1:p, L + 1)[l_pan]
      panel.grid()

      panel.superpose(
        x[order(x)], y[order(x)], subscripts, groups,
        type="l", col = data_col, pch = 16, cex = 0.4
      )

      panel.xyplot(
        time_g, alt_mod_gbl_est[[l]][,j], col = alt_mod_col,
        type = "l", lwd = 2
      )

      panel.xyplot(
        time_g, vmp_gbl_est[[l]][, j], col = vmp_col,
        type = "l", lwd = 1
      )
    }
  )

  print(fitted_gbl_plots)
}


panel_plots <- function(vmp_files, mcmc_files, mu_func, Psi_func, plot_dim,
                        plot_height, plot_width, logistic_mod=FALSE) {

  # Set the number of basis functions:

  L <- length(Psi_func)

  # Read the files:

  mu_vmp_res <- read.table(vmp_files[1], header=TRUE)
  mu_vmp_res <- as.matrix(mu_vmp_res)

  mu_mcmc_res <- read.table(mcmc_files[1], header=TRUE)
  mu_mcmc_res <- as.matrix(mu_mcmc_res)

  # Determine necessary parameters:

  n_sims <- nrow(mu_vmp_res)
  n_g <- ncol(mu_vmp_res)
  time_g <- seq(0, 1, length.out = n_g)

  # Summarise the results:

  mu_mcmc_summary <- matrix(NA, n_g, 3)
  mu_mcmc_summary[,1] <- apply(mu_mcmc_res, 2, quantile, 0.025)
  mu_mcmc_summary[,2] <- apply(mu_mcmc_res, 2, mean)
  mu_mcmc_summary[,3] <- apply(mu_mcmc_res, 2, quantile, 0.975)

  psi_vmp_res <- vector("list", length=L)
  psi_mcmc_summary <- vector("list", length=L)
  for(l in 1:L) {

    psi_vmp_res[[l]] <- read.table(vmp_files[l+1], header=TRUE)
    psi_vmp_res[[l]] <- as.matrix(psi_vmp_res[[l]])

    psi_mcmc_res <- read.table(mcmc_files[l+1], header=TRUE)
    psi_mcmc_res <- as.matrix(psi_mcmc_res)

    psi_mcmc_summary[[l]] <- matrix(NA, n_g, 3)
    psi_mcmc_summary[[l]][,1] <- apply(psi_mcmc_res, 2, quantile, 0.025)
    psi_mcmc_summary[[l]][,2] <- apply(psi_mcmc_res, 2, mean)
    psi_mcmc_summary[[l]][,3] <- apply(psi_mcmc_res, 2, quantile, 0.975)
  }

  vmp_res <- vector("list", length=L+1)
  vmp_res[[1]] <- mu_vmp_res
  mcmc_res <- vector("list", length=L+1)
  mcmc_res[[1]] <- mu_mcmc_summary
  for(l in 1:L) {

    vmp_res[[l+1]] <- psi_vmp_res[[l]]
    mcmc_res[[l+1]] <- psi_mcmc_summary[[l]]
  }

  # Establish the grid-based functions:

  mu_g <- mu_func(time_g)
  Psi_g <- matrix(NA, nrow=n_g, ncol=L)
  for(l in 1:L) {

    Psi_g[,l] <- Psi_func[[l]](time_g)
  }

  # Plot the results:

  time_g_gbl <- rep(time_g, L + 1)
  gbl_g_vec <- c(mu_g, as.vector(Psi_g))
  gbl_labels <- vector("list", length=(L+1))
  gbl_id <- rep(NA, L + 1)

  gbl_id[1] <- expression(mu (t))
  gbl_labels[[1]] <- rep(gbl_id[1], n_g)

  for(l in 1:L) {

    gbl_id[l+1] <- parse(text=paste("psi[", l, "] (t)", sep=""))
    bf_val <- eval(bquote(expression(psi[.(l)] (t))))
    gbl_labels[[l+1]] <- rep(bf_val, n_g)
  }

  gbl_labels <- do.call(c, gbl_labels)
  gbl_labels <- factor(gbl_labels, levels=gbl_id)

  strip.math <- function(
    which.given, which.panel, var.name, factor.levels, ...
  ) {

    fl <- gbl_id

    strip.default(which.given,which.panel,var.name,fl,...)
  }

  gbl_plots <- xyplot(
    gbl_g_vec ~ time_g_gbl | gbl_labels, groups=gbl_labels,
    data=data.frame(
      time_g_gbl=time_g_gbl, gbl_g_vec=gbl_g_vec,
      gbl_labels=gbl_labels
    ),
    layout= plot_dim, main="",
    strip=strip.math,
    par.strip.text=list(cex=0.8),
    par.settings = list(layout.heights = list(strip = 1.2)),
    xlab="time",
    ylab="mean and basis functions",
    as.table=TRUE,
    panel=function(x, y, subscripts, groups) {

      lPan <- panel.number()
      l <- rep(1:(L+1), each=1)[lPan]
      panel.grid()

      for(i in 1:n_sims) {

        panel.xyplot(
          time_g, vmp_res[[l]][i,],
          col="red", type="l", lwd=1
        )
      }

      panel.xyplot(
        time_g, mcmc_res[[l]][,2],
        col="deepskyblue2", type="l", lwd=1
      )

      panel.xyplot(
        time_g, mcmc_res[[l]][,1],
        col="deepskyblue2", type="l", lty=2, lwd=1
      )

      panel.xyplot(
        time_g, mcmc_res[[l]][,3],
        col="deepskyblue2", type="l", lty=2, lwd=1
      )

      panel.superpose(
        x[order(x)], y[order(x)], subscripts, groups,
        type="l", col="black", lwd=1
      )
    }
  )

  print(gbl_plots)

}



box_plot_fpca_sims <- function(file_name, plot_width, plot_height, log_acc=FALSE) {

  # Read the file:

  results <- read.table(file_name, header=TRUE)

  # Gather necessary parameters:

  L <- ncol(results) - 3
  N_vals <- unique(results$"N")
  n_sims <- max(results$"sim")

  # Extract accuracy scores:

  acc_vec <- as.vector(as.matrix(results[,-c(1, 2)]))

  if(log_acc) {

    acc_vec <- log(acc_vec)
    y_lab <- "log ISE"
  } else {

    y_lab <- "ISE"
  }

  # Label the curves:

  gbl_labels <- vector("list", length=(L+1))
  gbl_id <- rep(NA, L + 1)

  gbl_id[1] <- expression(mu (t))
  gbl_labels[[1]] <- rep(gbl_id[1], length(N_vals)*n_sims)

  for(l in 1:L) {

    gbl_id[l+1] <- parse(text=paste("psi[", l, "] (t)", sep=""))
    bf_val <- eval(bquote(expression(psi[.(l)] (t))))
    gbl_labels[[l+1]] <- rep(bf_val, length(N_vals)*n_sims)
  }

  gbl_labels <- do.call(c, gbl_labels)
  gbl_labels <- factor(gbl_labels, levels=gbl_id)

  N_labels <- rep(N_vals, each=n_sims)
  N_labels <- rep(N_labels, L + 1)
  N_labels <- factor(N_labels)

  strip.math <- function(
    which.given, which.panel, var.name, factor.levels, ...
  ) {

    fl <- gbl_id

    strip.default(which.given,which.panel,var.name,fl,...)
  }

  box_plots <- bwplot(
    acc_vec ~ N_labels | gbl_labels,
    data = data.frame(acc_vec=acc_vec, N_labels=N_labels, gbl_labels=gbl_labels),
    layout=c(3,1),
    xlab="number of response curves", ylab=y_lab,
    strip=strip.math,
    par.settings = list(layout.heights = list(strip = 1.2)),
    as.table=TRUE,
    panel=function(...) {

      panel.abline(h=0, col="red", lty = 2)
      panel.abline(h=-2, col="red", lty = 2)
      panel.abline(h=-4, col="red", lty = 2)
      panel.abline(h=-6, col="red", lty = 2)
      panel.abline(h=-8, col="red", lty = 2)
      panel.bwplot(...)
    }
  )

  print(box_plots)

}

