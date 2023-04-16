# temporary script - likely to be removed in the released package which will
# only include variational inference routines.
#
summarise_mcmc_multivariate <- function(stan_obj, C_g, Psi_g) {

  mcmc_samples <- rstan::extract(stan_obj, permuted=FALSE)
  n_mcmc <- dim(mcmc_samples)[1]
  time_g <- C_g[,2]
  n_g <- dim(C_g)[1]
  p <- length(Psi_g)

  fpca_params <- dimnames(mcmc_samples)$parameters
  L <- length(fpca_params[grep("beta_psi", fpca_params, fixed=TRUE)])/2/p
  N <- length(fpca_params[grep("zeta", fpca_params, fixed=TRUE)])/L

  mu_g_mcmc <- Psi_g_mcmc <- vector("list", length=p)
  for (j in 1:p) {

    beta_mu_cols <- fpca_params[grep(paste0("beta_mu[", j), fpca_params, fixed=TRUE)]
    beta_mu_mcmc <- mcmc_samples[, 1, beta_mu_cols]

    u_mu_cols <- fpca_params[grep(paste0("u_mu[", j), fpca_params, fixed=TRUE)]
    u_mu_mcmc <- mcmc_samples[, 1, u_mu_cols]

    nu_mu_mcmc <- t(cbind(beta_mu_mcmc, u_mu_mcmc))
    mu_g_mcmc[[j]] <- C_g%*%nu_mu_mcmc

    Psi_g_mcmc[[j]] <- vector("list", length=L)
    for(l in 1:L) {

      beta_psi_l_cols <- fpca_params[grep(paste0("beta_psi[", j, ",", l), fpca_params, fixed=TRUE)]
      beta_psi_mcmc <- mcmc_samples[, 1, beta_psi_l_cols]

      u_psi_l_cols <- fpca_params[grep(paste0("u_psi[", j, ",", l), fpca_params, fixed=TRUE)]
      u_psi_mcmc <- mcmc_samples[, 1, u_psi_l_cols]

      nu_psi_mcmc <- t(cbind(beta_psi_mcmc, u_psi_mcmc))
      Psi_g_mcmc[[j]][[l]] <- C_g%*%nu_psi_mcmc
    }

  }

  zeta_mcmc <- vector("list", length=N)
  for(i in 1:N) {

    zeta_i <- paste("zeta[", i, ",", sep="")
    zeta_i_cols <- fpca_params[grep(zeta_i, fpca_params, fixed=TRUE)]
    zeta_mcmc[[i]] <- mcmc_samples[, 1, zeta_i_cols]
  }

  one_N <- rep(1, N)
  mu_hat <- vector("list", length=n_mcmc)
  Psi_hat <- vector("list", length=n_mcmc)
  Zeta_hat <- vector("list", length=n_mcmc)
  for(k in 1:n_mcmc) {

    Psi <- vector("list", length=p)
    for (j in 1:p) {

      Psi[[j]] <- matrix(NA, n_g, L)
      for(l in 1:L) {

        Psi[[j]][,l] <- Psi_g_mcmc[[j]][[l]][,k]
      }

    }

    Zeta <- matrix(NA, N, L)
    for(i in 1:N) {

      Zeta[i,] <- zeta_mcmc[[i]][k,]
    }

    # Orthogonalisation:
    #
    mu_g_mcmc_k <- Reduce(c, lapply(mu_g_mcmc, function(ll) ll[,k]))
    Psi_k <- Reduce(rbind, Psi)

    svd_Psi <- svd(Psi_k)
    U_psi <- svd_Psi$u
    D_psi <- diag(svd_Psi$d)
    V_psi <- svd_Psi$v

    zeta_rotn <- t(Zeta %*% V_psi %*% D_psi)
    C_zeta <- cov(t(zeta_rotn))
    eigen_C <- eigen(C_zeta)
    Q <- eigen_C$vectors
    Lambda <- diag(eigen_C$values)

    Psi_tilde <- U_psi %*% Q %*% sqrt(Lambda)
    Zeta_tilde <- crossprod(zeta_rotn, Q %*% solve(sqrt(Lambda)))

    mu_hat[[k]] <- split(mu_g_mcmc_k, rep(1:p, each = n_g))

    Psi_hat[[k]] <- matrix(NA, p*n_g, L)
    Zeta_hat[[k]] <- matrix(NA, N, L)
    norm_vec <- rep(NA, L)
    time_int_vec <- seq(0, p, length = n_g*p)
    for(l in 1:L) {

      norm_vec[l] <- sqrt(trapint(time_int_vec, Psi_tilde[, l]^2))
      Psi_hat[[k]][, l] <- Psi_tilde[, l]/norm_vec[l]
      Zeta_hat[[k]][, l] <- norm_vec[l]*Zeta_tilde[, l]

      Psi_g_comb <- vector("list", length = p)
      for(j in 1:p) {

        Psi_g_comb[[j]] <- Psi_g[[j]][, l]
      }
      Psi_g_comb <- Reduce(c, Psi_g_comb)

      inner_prod_sign <- sign(cprod(Psi_g_comb, Psi_hat[[k]][, l]))
      if(inner_prod_sign == -1) {

        Psi_hat[[k]][, l] <- -Psi_hat[[k]][, l]
        Zeta_hat[[k]][, l] <- -Zeta_hat[[k]][, l]
      }
    }
    Psi_hat[[k]] <- lapply(split(Psi_hat[[k]], rep(1:p, each = n_g)), matrix, nrow = n_g, ncol = L)
  }


  # Summarise the MCMC outputs:
  #
  Y_g_mcmc_summary <- vector("list", length=N)
  for(i in 1:N) {

    Y_g_mcmc_summary[[i]] <- vector("list", length=p)

    for (j in 1:p) {

      Y_g_mcmc <- matrix(NA, n_g, n_mcmc)
      for(k in 1:n_mcmc) {

        # Y_g_mcmc[,k] <- mu_g_mcmc[,k] + Psi_hat[[k]]%*%Zeta_hat[[k]][i,]
        Y_g_mcmc[,k] <- mu_hat[[k]][[j]] + Psi_hat[[k]][[j]]%*%Zeta_hat[[k]][i,]

      }

      Y_g_mcmc_summary[[i]][[j]] <- matrix(NA, nrow=n_g, ncol=3)
      Y_g_mcmc_summary[[i]][[j]][,1] <- apply(Y_g_mcmc, 1, quantile, 0.025)
      Y_g_mcmc_summary[[i]][[j]][,2] <- apply(Y_g_mcmc, 1, mean)
      Y_g_mcmc_summary[[i]][[j]][,3] <- apply(Y_g_mcmc, 1, quantile, 0.975)

    }


  }

  zeta_mcmc_summary <- vector("list", length=N)
  for(i in 1:N) {

    zeta_mcmc_i <- matrix(NA, n_mcmc, 2)
    for(k in 1:n_mcmc) {

      zeta_mcmc_i[k,] <- Zeta_hat[[k]][i,1:2]
    }

    zeta_mcmc_mean <- apply(zeta_mcmc_i, 2, mean)
    zeta_mcmc_cov <- cov(zeta_mcmc_i)

    zeta_mcmc_ellipse <- ellipse(
      zeta_mcmc_cov,
      centre=zeta_mcmc_mean,
      level=0.95
    )

    zeta_mcmc_summary[[i]] <- list(zeta_mcmc_mean, zeta_mcmc_ellipse)
    names(zeta_mcmc_summary[[i]]) <- c("mean", "credible boundary")
  }

  gbl_mcmc_summary <- vector("list", length = L + 1)
  gbl_mcmc_summary[[1]] <- Reduce(cbind, lapply(1:p, function(j) apply(Reduce(cbind, lapply(mu_hat, function(mu_hat_k) mu_hat_k[[j]])), 1, mean)))

  for(l in 1:L) {
    gbl_mcmc_summary[[l+1]] <- matrix(NA, n_g, p)
    for(j in 1:p) {
      gbl_mcmc_summary[[l+1]][, j] <- Reduce("+", lapply(Psi_hat, function(Psi_hat_k) Psi_hat_k[[j]][,l]))/n_mcmc
    }
  }


  gbl_ci_mcmc_summary <- vector("list", length = L + 1)
  gbl_ci_mcmc_summary[[1]] <- vector("list", length = p)
  for(j in 1:p) {
    gbl_ci_mcmc_summary[[1]][[j]] <- t(apply(simplify2array(lapply(mu_hat, function(mu_hat_k) mu_hat_k[[j]])), 1,
                                             quantile, prob = c(0.025, 0.975)))
  }

  for(l in 1:L) {
    gbl_ci_mcmc_summary[[l+1]] <- vector("list", length = p) # matrix(NA, n_g, p)
    for(j in 1:p) {
      gbl_ci_mcmc_summary[[l+1]][[j]] <- t(apply(simplify2array(lapply(Psi_hat, function(Psi_hat_k) Psi_hat_k[[j]][,l])), 1,
                                                 quantile, prob = c(0.025, 0.975)))
    }
  }

  # Summary outputs:

  outputs <- list(Y_g_mcmc_summary, gbl_mcmc_summary, gbl_ci_mcmc_summary, zeta_mcmc_summary)
  names(outputs) <- c("Y_g_mcmc_summary", "gbl_mcmc_summary", "gbl_ci_mcmc_summary", "zeta_mcmc_summary")

  return(outputs)
}



summarise_mcmc <- function(stan_obj, C_g, Psi_g, use_logistic_mod=FALSE) {

  mcmc_samples <- rstan::extract(stan_obj, permuted=FALSE)
  n_mcmc <- dim(mcmc_samples)[1]
  time_g <- C_g[,2]
  n_g <- dim(C_g)[1]

  fpca_params <- dimnames(mcmc_samples)$parameters
  L <- length(fpca_params[grep("beta_psi", fpca_params, fixed=TRUE)])/2
  N <- length(fpca_params[grep("zeta", fpca_params, fixed=TRUE)])/L

  beta_mu_cols <- fpca_params[grep("beta_mu", fpca_params, fixed=TRUE)]
  beta_mu_mcmc <- mcmc_samples[, 1, beta_mu_cols]

  u_mu_cols <- fpca_params[grep("u_mu", fpca_params, fixed=TRUE)]
  u_mu_mcmc <- mcmc_samples[, 1, u_mu_cols]

  nu_mu_mcmc <- t(cbind(beta_mu_mcmc, u_mu_mcmc))
  mu_g_mcmc <- C_g%*%nu_mu_mcmc

  Psi_g_mcmc <- vector("list", length=L)
  for(l in 1:L) {

    beta_psi_l <- paste("beta_psi[", l, ",", sep="")
    beta_psi_l_cols <- fpca_params[grep(beta_psi_l, fpca_params, fixed=TRUE)]
    beta_psi_mcmc <- mcmc_samples[, 1, beta_psi_l_cols]

    u_psi_l <- paste("u_psi[", l, ",", sep="")
    u_psi_l_cols <- fpca_params[grep(u_psi_l, fpca_params, fixed=TRUE)]
    u_psi_mcmc <- mcmc_samples[, 1, u_psi_l_cols]

    nu_psi_mcmc <- t(cbind(beta_psi_mcmc, u_psi_mcmc))
    Psi_g_mcmc[[l]] <- C_g%*%nu_psi_mcmc
  }

  zeta_mcmc <- vector("list", length=N)
  for(i in 1:N) {

    zeta_i <- paste("zeta[", i, ",", sep="")
    zeta_i_cols <- fpca_params[grep(zeta_i, fpca_params, fixed=TRUE)]
    zeta_mcmc[[i]] <- mcmc_samples[, 1, zeta_i_cols]
  }

  one_N <- rep(1, N)
  Psi_hat <- vector("list", length=n_mcmc)
  Zeta_hat <- vector("list", length=n_mcmc)
  for(j in 1:n_mcmc) {

    Psi <- matrix(NA, n_g, L)
    for(l in 1:L) {

      Psi[,l] <- Psi_g_mcmc[[l]][,j]
    }

    Zeta <- matrix(NA, N, L)
    for(i in 1:N) {

      Zeta[i,] <- zeta_mcmc[[i]][j,]
    }

    Psi_svd <- svd(Psi)
    U_orth <- Psi_svd$u
    D_diag <- diag(Psi_svd$d)
    V_orth <- Psi_svd$v

    Zeta_rotn <- Zeta%*%V_orth%*%D_diag
    mu_Zeta_rotn <- apply(Zeta_rotn, 2, mean)

    mu_g_mcmc[,j] <- mu_g_mcmc[,j] + U_orth%*%mu_Zeta_rotn
    Zeta_shift <- Zeta_rotn - tcrossprod(one_N, mu_Zeta_rotn)

    eigen_Zeta_shift <- eigen(cov(Zeta_shift))
    Q <- eigen_Zeta_shift$vectors
    Lambda <- diag(eigen_Zeta_shift$values)
    S <- Q%*%sqrt(Lambda)

    Psi_hat[[j]] <- U_orth%*%S
    Zeta_hat[[j]] <- tcrossprod(Zeta_shift, solve(S))

    norm_const <- rep(NA, L)
    for(l in 1:L) {

      norm_const[l] <- sqrt(trapint(time_g, (Psi_hat[[j]][,l])^2))
      Psi_hat[[j]][,l] <- Psi_hat[[j]][,l]/norm_const[l]
      Zeta_hat[[j]][,l] <- norm_const[l]*Zeta_hat[[j]][,l]

      cprod_sign <- sign(cprod(Psi_hat[[j]][,l], Psi_g[,l]))
      if(cprod_sign==-1) {

        Psi_hat[[j]][,l] <- -Psi_hat[[j]][,l]
        Zeta_hat[[j]][,l] <- -Zeta_hat[[j]][,l]
      }
    }
  }

  # Summarise the MCMC outputs:

  Y_g_mcmc_summary <- vector("list", length=N)
  for(i in 1:N) {

    Y_g_mcmc <- matrix(NA, n_g, n_mcmc)
    for(j in 1:n_mcmc) {

      Y_g_mcmc[,j] <- mu_g_mcmc[,j] + Psi_hat[[j]]%*%Zeta_hat[[j]][i,]
    }

    Y_g_mcmc_summary[[i]] <- matrix(NA, nrow=n_g, ncol=3)
    Y_g_mcmc_summary[[i]][,1] <- apply(Y_g_mcmc, 1, quantile, 0.025)
    Y_g_mcmc_summary[[i]][,2] <- apply(Y_g_mcmc, 1, mean)
    Y_g_mcmc_summary[[i]][,3] <- apply(Y_g_mcmc, 1, quantile, 0.975)
  }

  zeta_mcmc_summary <- vector("list", length=N)
  for(i in 1:N) {

    zeta_mcmc_i <- matrix(NA, n_mcmc, 2)
    for(j in 1:n_mcmc) {

      zeta_mcmc_i[j,] <- Zeta_hat[[j]][i,1:2]
    }

    zeta_mcmc_mean <- apply(zeta_mcmc_i, 2, mean)
    zeta_mcmc_cov <- cov(zeta_mcmc_i)

    zeta_mcmc_ellipse <- ellipse(
      zeta_mcmc_cov,
      centre=zeta_mcmc_mean,
      level=0.95
    )

    zeta_mcmc_summary[[i]] <- list(zeta_mcmc_mean, zeta_mcmc_ellipse)
    names(zeta_mcmc_summary[[i]]) <- c("mean", "credible boundary")
  }

  gbl_mcmc_summary <- matrix(NA, n_g, L+1)
  gbl_mcmc_summary[,1] <- apply(mu_g_mcmc, 1, mean)
  gbl_mcmc_summary[,2:(L+1)] <- Reduce("+", Psi_hat)/n_mcmc

  # Summary outputs:

  outputs <- list(Y_g_mcmc_summary, gbl_mcmc_summary, zeta_mcmc_summary)
  names(outputs) <- c("Y_g_mcmc_summary", "gbl_mcmc_summary", "zeta_mcmc_summary")

  return(outputs)
}
