#' @export
flip_sign <- function(vec_flip, list_Psi_hat, Zeta_hat, zeta_ellipse) {

  # list_Psi_hat <- lapply(list_Psi_hat, function(Psi_hat_var) {
  #   Psi_hat_var[, seq_along(vec_flip)] <- sweep(Psi_hat_var[, seq_along(vec_flip)], 2, vec_flip, "*") # first two eigenfunctions
  #   Psi_hat_var
  # })

  list_Psi_hat <- lapply(seq_along(vec_flip), function(ll) {
    list_Psi_hat[[ll]] * vec_flip[ll]
  })

  Zeta_hat[,seq_along(vec_flip)] <- sweep(Zeta_hat[,seq_along(vec_flip)], 2, vec_flip, "*")

  zeta_ellipse <- lapply(zeta_ellipse, function(zeta_ellipse_subj) sweep(zeta_ellipse_subj, 2, vec_flip, "*"))

  create_named_list(list_Psi_hat, Zeta_hat, zeta_ellipse)
}

#' @export
display_eigenfunctions <- function(L, time_g, mu_g, Psi_g,
                                   mu_hat, list_Psi_hat, format_univ = FALSE,
                                   mu_hat_add = NULL, list_Psi_hat_add = NULL,
                                   mu_hat_ci = NULL, list_Psi_hat_ci = NULL,
                                   lwd = 2, data_col = "red", p_sample = 1,
                                   vec_col_add = NULL, vec_lwd = NULL) { # perso

  ylim <- c(min(c(unlist(mu_g),
                  as.vector(unlist(mu_hat)),
                  as.vector(unlist(mu_hat_add)),
                  as.vector(unlist(mu_hat_ci)),
                  unlist(Psi_g),
                  as.vector(unlist(list_Psi_hat)),
                  as.vector(unlist(list_Psi_hat_add)),
                  as.vector(unlist(list_Psi_hat_ci)))),
            max(c(unlist(mu_g),
                  as.vector(unlist(mu_hat)),
                  as.vector(unlist(mu_hat_add)),
                  as.vector(unlist(mu_hat_ci)),
                  unlist(Psi_g),
                  as.vector(unlist(list_Psi_hat)),
                  as.vector(unlist(list_Psi_hat_add)),
                  as.vector(unlist(list_Psi_hat_ci)))))

  p_sample <- sort(p_sample)

  par(mfrow = c(1+L, length(p_sample))) #, mar = c(5.1, 4.1, 4.1, 2.1))
  
  if(format_univ){
    mu_g <- list(mu_g)
    Psi_g <- list(Psi_g)
    mu_hat <- matrix(mu_hat)
    list_Psi_hat <- lapply(split(list_Psi_hat, col(list_Psi_hat)), matrix)
  }

  for (j in p_sample) {
    par(mar = c(1,4.5,4,1))
    if (is.list(mu_hat)) {
      mu_hat_j <- lapply(mu_hat, function(ll) ll[,j])
    } else {
      mu_hat_j <- mu_hat[,j]
    }

    if (!is.null(mu_hat_add)) {
      if (is.list(mu_hat_add)) {
        mu_hat_add_j <- lapply(mu_hat_add, function(ll) ll[,j])
      } else {
        mu_hat_add_j <- mu_hat_add[,j]
      }
    } else {
      mu_hat_add_j <- NULL
    }

    if (!is.null(mu_hat_ci)) {
      mu_hat_ci_j <- mu_hat_ci[[j]]
    } else {
      mu_hat_ci_j <- NULL
    }
    display_function(time_g,
                     mu_g[[j]], mu_hat_j, mu_hat_add_j,
                     mu_hat_ci_j,
                     fct_name = ifelse(j == p_sample[1], parse(text=paste0("mu", "(t)")), ""),#expression(mu (t)), ""),
                     main = paste0("Variable ", j, "\n "),
                     ylim= ylim,
                     add_fct_lwd = vec_lwd,
                     add_fct_cols = vec_col_add,
                     col_fct = data_col)
    #               , add_fct_types = c("Simulated", "VB"))

  }

  for (l in 1:L) {
    if (l == L){
      par(mar = c(4,4.5,1,1))
    } else {
      par(mar = c(2.5,4.5,2.5,1))
    }
    for (j in p_sample) {

      if (is.list(mu_hat)) { # keep is.list(mu_hat) as list_Psi_hat will be a list in both cases
        list_Psi_hat_l_j <- lapply(list_Psi_hat, function(ll) ll[[l]][,j])
      } else {
        list_Psi_hat_l_j <- list_Psi_hat[[l]][,j]
      }

      if (!is.null(list_Psi_hat_add)) {
        if (is.list(mu_hat_add)) {
          list_Psi_hat_add_l_j <- lapply(list_Psi_hat_add, function(ll) ll[[l]][,j])
        } else {
          list_Psi_hat_add_l_j <- list_Psi_hat_add[[l]][,j]
        }
      } else {
        list_Psi_hat_add_l_j <- NULL
      }

      if (!is.null(list_Psi_hat_ci)) {
        list_Psi_hat_ci_l_j <- list_Psi_hat_ci[[l]][[j]]
      } else {
        list_Psi_hat_ci_l_j <- NULL
      }

      display_function(time_g, Psi_g[[j]][,l], list_Psi_hat_l_j, list_Psi_hat_add_l_j,
                       list_Psi_hat_ci_l_j,
                       fct_name = ifelse(j == p_sample[1],
                                         parse(text=paste0("psi[", l, "](t)")), ""),
                       ylim =  ylim,
                       add_fct_lwd = vec_lwd,
                       add_fct_cols = vec_col_add,
                       col_fct = data_col)
    }
  }

}

display_function <- function(time_g, fct, add_fct1 = NULL, add_fct2 = NULL,
                             add_fct_ci = NULL,
                             fct_name = NULL,
                             main = NULL,
                             add_fct_types = NULL,
                             add_fct_lwd = NULL,
                             add_fct_cols = NULL,
                             col_fct = "red",
                             ylim = NULL) {

  if (is.null(ylim)) {
    ylim <- c(min(c(fct, unlist(add_fct1), unlist(add_fct2))),
              max(c(fct, unlist(add_fct1), unlist(add_fct2))))
  }

  if (is.null(add_fct_lwd)) {
    add_fct_lwd <- rep(2, 2)
  }

  plot(time_g, fct, type = "l", xlab = "time", ylab = fct_name, main = main,
       col = col_fct, lwd = 2,
       ylim = ylim)

  if (!is.null(add_fct1)) {
    add_fct_cols[1] <- ifelse(is.null(add_fct_cols[1]), "grey40", add_fct_cols[1])
    if (is.list(add_fct1)) {
      for (j in seq_along(add_fct1)) {
        lines(time_g, add_fct1[[j]], col = add_fct_cols[1], lwd = add_fct_lwd[1])
      }
      lines(time_g, fct, col = col_fct, lwd = add_fct_lwd[1])
    } else {
      lines(time_g, add_fct1, col = add_fct_cols[1], lwd = add_fct_lwd[1])
    }
  }

  if (!is.null(add_fct2)) {
    add_fct_cols[2] <- ifelse(is.null(add_fct_cols[2]) | is.na(add_fct_cols[2]), "blue", add_fct_cols[2])

    if (is.list(add_fct2)) {
      for (j in seq_along(add_fct2)) {
        lines(time_g, add_fct2[[j]], col = add_fct_cols[2], lwd = add_fct_lwd[2])
      }
      if (!is.list(add_fct1)) {
        lines(time_g, add_fct1, col = add_fct_cols[1], lwd = add_fct_lwd[1])
      }
      lines(time_g, fct, col = col_fct, lwd = add_fct_lwd[1])
    } else {
      lines(time_g, add_fct2, col = add_fct_cols[2], lwd = add_fct_lwd[2])
    }
  }

  if (!is.null(add_fct_ci)) {
    lines(time_g, add_fct_ci[,1], col = add_fct_cols[1], lwd = add_fct_lwd[1], lty = 2)
    lines(time_g, add_fct_ci[,2], col = add_fct_cols[1], lwd = add_fct_lwd[1], lty = 2)
  }

  if (!is.null(add_fct_types)) {
    legend("topleft", legend = add_fct_types, col = c(col_fct, add_fct_cols),
           bty = "n", lwd = c(add_fct_lwd[1], add_fct_lwd))
  }
}

#' @export
display_fit <- function(N_sample, time_obs, time_g, Y, Y_hat, Y_low, Y_upp, offset = 0.1) {

  N_sample <- sort(N_sample)

  if (is.null(Y_hat)) {
    list_ylim <- c(min(unlist(Y[N_sample]))-offset,
                               max(unlist(Y[N_sample]))+offset)
  } else {
    list_ylim <- c(min(unlist(Y[N_sample]),
                                   unlist(Y_low[N_sample]))-offset,
                        max(unlist(Y[N_sample]),
                            unlist(Y_upp[N_sample]))+offset)
  }

  par(mfrow = c(length(N_sample), 1))
  for (i in N_sample) {
    if (i == N_sample[length(N_sample)]) {
      par(mar = c(4,4.5,1,1))
    } else if (i == N_sample[1]){
      par(mar = c(1,4.5,4,1))
    } else {
      par(mar = c(2.5,4.5,2.5,1))
    }

    plot(time_obs[[i]], Y[[i]],
         main = "",
         xlab = ifelse(i == N_sample[length(N_sample)], "time", ""),
         ylab = parse(text=paste0("Y[", i, "]")),
         pch = 20, col = "grey55",
         xlim = c(0, 1),
         ylim = list_ylim)
    # c(min(Y[[i]], Y_low[[i]]))-offset,
    #          max(Y[[i]], Y_upp[[i]]))+offset))
    if (!is.null(Y_hat)) {
      lwd <- 1.2
      lines(time_g, Y_hat[[i]], col="black", lwd = lwd)
      lines(time_g, Y_low[[i]], col="black",lwd = lwd,lty = 2)
      lines(time_g, Y_upp[[i]], col="black",lwd = lwd,lty = 2)
    }

  }

}


#' @export
display_fit_list <- function(p_sample, N_sample, time_obs, time_g, Y,
                             Y_hat, Y_low, Y_upp,
                             Y_hat_add = NULL, Y_low_add = NULL, Y_upp_add = NULL, offset = 0.1,
                             col_data = "grey55", col = "black", col_add = "blue",
                             lwd = 1.2, lwd_add = 1.2) {

  p_sample <- sort(p_sample)
  N_sample <- sort(N_sample)

  list_ylim <- list()
  for (j in p_sample) {
    if (is.null(Y_hat)) {
      list_ylim <- append(list_ylim,
                          list(c(min(unlist(sapply(Y[N_sample], function(Y_i)Y_i[[j]])))-offset,
                                 max(unlist(sapply(Y[N_sample], function(Y_i)Y_i[[j]])))+offset)))
    } else {
      vec_lim <- c(min(unlist(sapply(Y[N_sample], function(Y_i)Y_i[[j]])),
                       unlist(sapply(Y_low[N_sample], function(Y_i)Y_i[[j]]))),
                   max(unlist(sapply(Y[N_sample], function(Y_i)Y_i[[j]])),
                       unlist(sapply(Y_upp[N_sample], function(Y_i)Y_i[[j]]))))

      if (!is.null(Y_hat_add)) {
        vec_lim <- c(min(vec_lim[1], unlist(sapply(Y_low_add[N_sample], function(Y_i)Y_i[[j]]))),
                     max(vec_lim[2], unlist(sapply(Y_upp_add[N_sample], function(Y_i)Y_i[[j]]))))
      }
      list_ylim <- append(list_ylim,
                          list(c(vec_lim[1]-offset, vec_lim[2]+offset)))
    }
  }

  par(mfrow = c(length(N_sample), length(p_sample)))
  for (i in N_sample) {
    if (i == N_sample[length(N_sample)]) {
      par(mar = c(4,4.5,1.5,1))
    } else if (i == N_sample[1]){
      par(mar = c(1.5,4.5,4,1))
    } else {
      par(mar = c(2.6,4.5,2.6,1))
    }
    jj <- 1
    for (j in p_sample) {
      plot(time_obs[[i]][[j]], Y[[i]][[j]],
           main = ifelse(i == N_sample[1], paste0("Variable ", j, "\n"), ""),
           xlab = ifelse(i == N_sample[length(N_sample)], "time", ""),
           ylab = ifelse(j == p_sample[1], parse(text=paste0("Y[", i, "]")), ""),
           pch = 20, col = col_data,
           xlim = c(0, 1),
           ylim = list_ylim[[jj]])
      # c(min(Y[[i]][[j]], Y_low[[i]][[j]]))-offset,
      #          max(Y[[i]][[j]], Y_upp[[i]][[j]]))+offset))
      if (!is.null(Y_hat)) {
        lines(time_g, Y_hat[[i]][[j]], col=col, lwd = lwd)
        lines(time_g, Y_low[[i]][[j]], col=col,lwd = lwd,lty = 2)
        lines(time_g, Y_upp[[i]][[j]], col=col,lwd = lwd,lty = 2)
      }

      if (!is.null(Y_hat_add)) {
        lines(time_g, Y_hat_add[[i]][[j]], col=col_add, lwd = lwd_add)
        lines(time_g, Y_low_add[[i]][[j]], col=col_add,lwd = lwd_add,lty = 2)
        lines(time_g, Y_upp_add[[i]][[j]], col=col_add,lwd = lwd_add,lty = 2)
      }
      jj <- jj +1
    }

  }

}




#' @export
plot_scores <- function(N_sample, Zeta,
                        Zeta_hat, zeta_ellipse, p = 1,
                        Zeta_hat_add = NULL, zeta_ellipse_add = NULL,
                        vec_col = c("black", "blue"), data_col = "red", mfrow = NULL) {

  n_sample <- length(N_sample)

  par(mfrow = c(1, 1))
  zeta_labels <- vector("list", length = n_sample)
  zeta_id <- rep(NA, n_sample)

  if(is.null(mfrow)) {
    mfrow <- c(floor(sqrt(n_sample)), ceiling(sqrt(n_sample)))
  }

  for(i in 1:n_sample) {

    N_i <- N_sample[i]

    zeta_id[i] <- parse(text=paste("zeta[", N_i, "]", sep=""))
    zeta_val <- eval(bquote(expression(zeta[.(N_i)])))
    zeta_labels[[i]] <- rep(zeta_val, nrow(zeta_ellipse[[N_i]]))
  }
  zeta_labels <- do.call(c, zeta_labels)
  zeta_labels <- factor(zeta_labels, levels=zeta_id)

  strip.math <- function(
    which.given, which.panel, var.name, factor.levels, ...
  ) {

    fl <- zeta_id

    strip.default(which.given,which.panel,var.name,fl,...)
  }

  zeta_ellipse_mat <- Reduce(rbind, zeta_ellipse[N_sample])
  zeta_ellipse_x <- zeta_ellipse_mat[,1]
  zeta_ellipse_y <- zeta_ellipse_mat[,2]


  if (!is.null(Zeta_hat_add)) {
    zeta_ellipse_add_mat <- Reduce(rbind, zeta_ellipse_add[N_sample])
    zeta_ellipse_add_x <- zeta_ellipse_add_mat[,1]
    zeta_ellipse_add_y <- zeta_ellipse_add_mat[,2]
  }


  if (is.list(Zeta)) {
    xtmp_z <- as.vector(sapply(N_sample, function(ns) c(zeta_ellipse[[ns]][,1],
                                                        unlist(lapply(Zeta, function(zz) zz[ns, 1])),
                                                        Zeta_hat[ns,1])))
    ytmp_z <- as.vector(sapply(N_sample, function(ns) c(zeta_ellipse[[ns]][,2],
                                                        unlist(lapply(Zeta, function(zz) zz[ns, 2])),
                                                        Zeta_hat[ns,2])))

  } else {
    xtmp_z <- as.vector(sapply(N_sample, function(ns) c(zeta_ellipse[[ns]][,1],
                                                        Zeta[ns, 1],
                                                        Zeta_hat[ns,1])))
    ytmp_z <- as.vector(sapply(N_sample, function(ns) c(zeta_ellipse[[ns]][,2],
                                                        Zeta[ns, 2],
                                                        Zeta_hat[ns,2])))
  }
  xlim_z <- c(min(xtmp_z), max(xtmp_z))
  ylim_z <- c(min(ytmp_z), max(ytmp_z))

  score_plots <- xyplot(
    zeta_ellipse_y ~ zeta_ellipse_x | zeta_labels, groups = zeta_labels,
    data=data.frame(
      zeta_ellipse_x = zeta_ellipse_x, zeta_ellipse_y = zeta_ellipse_y,
      zeta_labels = zeta_labels
    ),
    layout=mfrow, main="",
    strip=strip.math,
    xlim = 1.1*xlim_z,
    ylim = 1.1*ylim_z,
    par.strip.text=list(cex=0.8),
    par.settings = list(layout.heights = list(strip = 1),
                        strip.background=list(col="white")),
    # key=list(space="top",
    #          lines=list(col=c("grey55","blue"), lty=1, lwd=1),
    #          text=list(c("Simulated scores", "Estimated scores"))),
    xlab="Scores FPC 1",
    ylab="Scores FPC 2",
    as.table=TRUE,
    panel=function(x, y, subscripts, groups) {

      iPan <- panel.number()
      i <- rep(1:n_sample, each=1)[iPan]
      # panel.grid()
      panel.xyplot(
        zeta_ellipse[[N_sample[i]]][,1], zeta_ellipse[[N_sample[i]]][,2],
        col=vec_col[1], type="l", lwd=1.5
      )
      panel.xyplot(
        Zeta_hat[N_sample[i],1], Zeta_hat[N_sample[i],2],
        col= vec_col[1], type="p", pch=16, cex=0.7
      )

      if (is.list(Zeta)) {
        for (j in 1:p) {
          panel.xyplot(
            Zeta[[j]][N_sample[i], 1], Zeta[[j]][N_sample[i], 2],
            # col=grDevices::adjustcolor(data_col, alpha.f = 1/j^0.9),
            col = data_col, type="p", pch=16, cex=0.7
          )
        }

      } else {
        panel.xyplot(
          Zeta[N_sample[i], 1], Zeta[N_sample[i], 2],
          col=data_col, type="p", pch=16, cex=0.7
        )
      }

      if (!is.null(Zeta_hat_add)) {
        panel.xyplot(
          zeta_ellipse_add[[N_sample[i]]][,1], zeta_ellipse_add[[N_sample[i]]][,2],
          col=vec_col[2], type="l", lwd=1.5
        )
        panel.xyplot(
          Zeta_hat_add[N_sample[i],1], Zeta_hat_add[N_sample[i],2],
          col= vec_col[2], type="p", pch=16, cex=0.7
        )
      }
    }
  )

  print(score_plots)
}
