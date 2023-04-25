#' @export
mu_func <- function(t, j = 1, alpha = 3) (-1)^j*alpha*sin((2*pi+j)*t) # default values = 1 for the case format_univ = T

psi_1 <- function(t, j = 1, p = 1) (-1)^j * sqrt(2/p)*cos(2*pi*t)
psi_2 <- function(t, j = 1, p = 1) (-1)^j * sqrt(2/p)*sin(2*pi*t)

#' @export
Psi_func <- function(time_obs, j = 1, p = 1) {
  ans <- cbind(psi_1(time_obs, j, p), psi_2(time_obs, j, p))
  return(ans)
}

#' @export
generate_fpca_data <- function(N, p, n, L, n_g, vec_sd_eps, mu_func, Psi_func,
                               time_obs = NULL, format_univ = FALSE,
                               generate_from_univ = FALSE,
                               vec_sd_zeta = NULL,
                               vec_rho_Zeta = NULL) {

  if (format_univ & p == 1) { # don't use it, not helpful.
    gauss_fpca_data(N, n, L, n_g, vec_sd_eps, mu_func, Psi_func, time_obs,
                    vec_sd_zeta = vec_sd_zeta)
  } else {
    if (p > 1) {
      stopifnot(!format_univ) # multivariate format required as p > 1
    } else {
      stopifnot(!generate_from_univ | !is.null(vec_rho_Zeta))
    }

    gauss_mfpca_data(N, p, n, L, n_g, vec_sd_eps, mu_func,
                     Psi_func, time_obs,
                     generate_from_univ = generate_from_univ,
                     vec_sd_zeta = vec_sd_zeta,
                     vec_rho_Zeta = vec_rho_Zeta)

  }
}


gauss_fpca_data <- function(N, n, L, n_g, vec_sd_eps, mu_func, Psi_func,
                            vec_sd_zeta = NULL, time_obs = NULL) {

  if (is.null(time_obs)) {
    time_obs <- sapply(n, runif)
    time_obs <- lapply(time_obs, sort)
  }

  # Set the scores:
  if (is.null(vec_sd_zeta)) {
    vec_sd_zeta <- 1/(1:L)
  } else {
    stopifnot(length(vec_sd_zeta) == L)
  }

  Zeta <- vector("list", length=N)
  for(i in 1:N) {

    Zeta[[i]] <- rep(NA, L)

    for(l in 1:L) {

      Zeta[[i]][l] <- rnorm(1, 0, vec_sd_zeta[l])
    }
  }

  Y <- vector("list", length=N)
  for(i in 1:N) {

    mu_t <- mu_func(time_obs[[i]], j = 1)

    Psi_t <- matrix(NA, nrow=n[i], ncol=L)
    for(l in 1:L) {

      Psi_t[,l] <- Psi_func(time_obs[[i]], j = 1, p = 1)[,l]
    }

    epsilon <- rnorm(n[i], 0, vec_sd_eps)

    Y_hat <- mu_t + as.vector(Psi_t%*%Zeta[[i]])
    Y[[i]] <- Y_hat + epsilon
  }

  time_g <- seq(0, 1, length.out=n_g)

  mu_g <- mu_func(time_g)
  Psi_g <- matrix(NA, nrow=n_g, ncol=L)
  for(l in 1:L) {

    Psi_g[,l] <- Psi_func(time_g)[,l]
  }

  Zeta <- Reduce(rbind, Zeta)

  create_named_list(time_obs, Zeta, time_g, mu_g, Psi_g, Y)

}


get_Zeta <- function(N, L, vec_sd_zeta = NULL) {

  if (is.null(vec_sd_zeta)) {
    vec_sd_zeta <- 1/(1:L)
  } else {
    stopifnot(length(vec_sd_zeta) == L)
  }

  Zeta <- matrix(NA, N, L)
  for(i in 1:N) {
    Zeta[i,] <- mvrnorm(1, rep(0, L), diag(vec_sd_zeta^2))
  }

  Zeta

}

get_list_corr_Zeta_univ <- function(N, L, p, vec_sd_zeta = NULL,
                                    vec_rho = 1/(2:(L+1))^0.2) { # vec_rho = correlation of the 1st, 2nd, etc sets of scores across the p variables (expected to decrease as l increases... )

  if (is.null(vec_sd_zeta)) {
    vec_sd_zeta <- 1/(1:L)
  } else {
    stopifnot(length(vec_sd_zeta) == L)
  }

  list_Zeta <- vector("list", p)
  for (l in 1:L) {
    R <- matrix(vec_rho[l], nrow = p,  ncol = p)
    diag(R) <- 1
    L_mat <- t(chol(R))
    tZ <- matrix(sapply(1:p, function(j) rnorm(N, 0, vec_sd_zeta[l])),
                 ncol = N, byrow = TRUE)
    Zeta_l <- as.matrix(t(L_mat %*% tZ))

    list_Zeta <- lapply(1:p, function(j) cbind(list_Zeta[[j]], Zeta_l[,j]))
  }

  list_Zeta
}


get_Y <- function(N, n, p, time_obs, Zeta, vec_sd_eps, mu_func, Psi_func) {

  Y <- vector("list", length = N)
  for(i in 1:N) {

    Y[[i]] <- vector("list", length = p)
    for(j in 1:p) {

      resid_vec <- rnorm(n[i, j], 0, vec_sd_eps[j])
      mean_vec <- mu_func(time_obs[[i]][[j]], j = j) +
        Psi_func(time_obs[[i]][[j]], j = j, p = p) %*% Zeta[i, ]
      Y[[i]][[j]] <- as.vector(mean_vec + resid_vec)
    }
  }

  Y
}


gauss_mfpca_data <- function(N, p, n, L, n_g, vec_sd_eps, mu_func, Psi_func,
                             time_obs = NULL, generate_from_univ = FALSE,
                             vec_sd_zeta = NULL, vec_rho_Zeta = NULL) {

  if (is.null(time_obs)) {
    time_obs <- vector("list", length = N)
    for(i in 1:N) {

      time_obs[[i]] <- vector("list", length = p)
      for(j in 1:p) {

        time_obs[[i]][[j]] <- sort(runif(n[i, j]))
      }
    }
  }

  time_g <- seq(0, 1, length.out = n_g)

  mu_g <- vector("list", length = p)
  Psi_g <- vector("list", length = p)
  for(j in 1:p) {
    mu_g[[j]] <- mu_func(time_g, j = j)
    Psi_g[[j]] <- Psi_func(time_g, j = j, p = p)
  }

  if (generate_from_univ) {

    if (is.null(vec_rho_Zeta) | all(vec_rho_Zeta == 0)) {
      corr_Zeta <- FALSE
      list_Zeta <- NULL
    } else {
      corr_Zeta <- TRUE
      list_Zeta <- get_list_corr_Zeta_univ(N, L, p, vec_sd_zeta = vec_sd_zeta,
                                           vec_rho = vec_rho_Zeta)
    }

    for (j in 1:p) {

      if (corr_Zeta) {
        Zeta_univ <- list_Zeta[[j]]
      } else {
        Zeta_univ <- get_Zeta(N, L, vec_sd_zeta = vec_sd_zeta)
        list_Zeta <- append(list_Zeta, list(Zeta_univ))
      }

      Y_univ <- get_Y(N, n, p, time_obs, Zeta_univ, vec_sd_eps, mu_func, Psi_func)
      if (j == 1) {
        Y <- lapply(Y_univ, function(ll) ll[1])
      } else {
        Y <- lapply(seq_along(Y), function(ii) append(Y[[ii]], Y_univ[[ii]][j]))
      }
    }
    Zeta <- list_Zeta

  } else {

    Zeta <- get_Zeta(N, L, vec_sd_zeta = vec_sd_zeta)

    Y <- get_Y(N, n, p, time_obs, Zeta, vec_sd_eps, mu_func, Psi_func)

  }

  create_named_list(time_obs, Zeta, time_g, mu_g, Psi_g, Y)

}
