#' @export
generate_fpca_data <- function(N, p, n, K, L, n_g, sigma_zeta_vec,
                               sigma_eps, mu_func, Psi_func, time_obs = NULL,
                               format_univ = FALSE) {

  if (format_univ & p == 1) {
    gauss_fpca_data(N, n, K, L, n_g, sigma_zeta_vec, sigma_eps, mu_func, Psi_func,
                    time_obs)
  } else {
    if (p > 1) {
      stopifnot(!format_univ) # multivariate format required as p > 1
    }
    gauss_mfpca_data(N, p, n, K, L, n_g, sigma_zeta_vec, sigma_eps, mu_func,
                     Psi_func, time_obs)
  }
}


gauss_fpca_data <- function(N, n, K, L, n_g, sigma_zeta_vec, sigma_eps,
                            mu_func, Psi_func, time_obs = NULL) {

  # Determine necessary parameters:

  # Set up fixed parameters:

  if (is.null(time_obs)) {
    time_obs <- sapply(n, runif)
    time_obs <- lapply(time_obs, sort)
  }

  unique_time_obs <- sort(unique(Reduce(c, time_obs)))
  int_knots <- quantile(unique_time_obs, seq(0,1,length=K)[-c(1,K)])

  X <- vector("list", length=N)
  Z <- vector("list", length=N)
  C <- vector("list", length=N)
  for(i in 1:N) {

    X[[i]] <- X_design(time_obs[[i]])
    Z[[i]] <- ZOSull(time_obs[[i]], range.x=c(0, 1), intKnots=int_knots)
    C[[i]] <- cbind(X[[i]], Z[[i]])
  }

  # Set the scores:

  Zeta <- vector("list", length=N)
  for(i in 1:N) {

    Zeta[[i]] <- rep(NA, L)

    for(l in 1:L) {

      Zeta[[i]][l] <- rnorm(1, 0, sigma_zeta_vec[l])
    }
  }

  # Set up curve observations:

  mu_t <- vector("list", length=N)
  Psi_t <- vector("list", length=N)
  Y <- vector("list", length=N)
  for(i in 1:N) {

    mu_t[[i]] <- mu_func(time_obs[[i]], j = 1)

    Psi_t[[i]] <- matrix(NA, nrow=n[i], ncol=L)
    for(l in 1:L) {

      Psi_t[[i]][,l] <- Psi_func(time_obs[[i]], j = 1, p = 1)[,l]
    }

    epsilon <- rnorm(n[i], 0, sigma_eps)

    Y_hat <- mu_t[[i]] + as.vector(Psi_t[[i]]%*%Zeta[[i]])
    Y[[i]] <- Y_hat + epsilon
  }

  # Set up plotting grid

  time_g <- seq(0, 1, length.out=n_g)

  X_g <- X_design(time_g)
  Z_g <- ZOSull(time_g, range.x=c(0, 1), intKnots=int_knots)
  C_g <- cbind(X_g, Z_g)

  mu_g <- mu_func(time_g, j = 1)
  Psi_g <- matrix(NA, nrow=n_g, ncol=L)
  for(l in 1:L) {

    Psi_g[,l] <- Psi_func(time_g, j = 1, p = 1)[,l]
  }

  Zeta <- Reduce(rbind, Zeta)

  create_named_list(time_obs, time_g, int_knots, X, Z, C, X_g, Z_g, C_g, Zeta,
                    mu_t, Psi_t, mu_g, Psi_g, Y
  )

}

gauss_mfpca_data <- function(
    N, p, n, K, L, n_g, sigma_zeta_vec,
    sigma_eps, mu_func, Psi_func, time_obs = NULL) {

  # Set up fixed parameters

  if (is.null(time_obs)) {
    time_obs <- vector("list", length = N)
    for(i in 1:N) {

      time_obs[[i]] <- vector("list", length = p)
      for(j in 1:p) {

        time_obs[[i]][[j]] <- sort(runif(n[i, j]))
      }
    }
  }

  time_vec <- unlist(time_obs)
  t_min <- 1.01*min(time_vec) - 0.01*max(time_vec)
  t_max <- 1.01*max(time_vec) - 0.01*min(time_vec)
  int_knots <- quantile(unique(time_vec), seq(0, 1, length = K)[-c(1, K)])

  C <- vector("list", length = N)
  for(i in 1:N) {

    C[[i]] <- vector("list", length = p)
    for(j in 1:p) {

      X <- X_design(time_obs[[i]][[j]])
      Z <- ZOSull(time_obs[[i]][[j]], range.x = c(0, 1), intKnots = int_knots)
      C[[i]][[j]] <- cbind(X, Z)
    }
  }

  # Set up plotting grid

  time_g <- seq(0, 1, length.out = n_g)

  X_g <- X_design(time_g)
  Z_g <- ZOSull(time_g, range.x = c(0, 1), intKnots = int_knots)
  C_g <- cbind(X_g, Z_g)

  mu_g <- vector("list", length = p)
  Psi_g <- vector("list", length = p)
  for(j in 1:p) {

    mu_g[[j]] <- mu_func(time_g, j = j)
    Psi_g[[j]] <- Psi_func(time_g, j = j, p = p)
  }

  # Simulate the scores

  Zeta <- matrix(NA, N, L)
  for(i in 1:N) {

    Zeta[i,] <- mvrnorm(1, rep(0, L), diag(sigma_zeta_vec^2))
  }

  # List of length N, containing n x p matrices (i.e., one matrix per subject)

  Y <- vector("list", length = N)
  for(i in 1:N) {

    Y[[i]] <- vector("list", length = p)
    for(j in 1:p) {

      resid_vec <- rnorm(n[i, j], 0, sigma_eps[j])
      mean_vec <- mu_func(time_obs[[i]][[j]], j = j) +
        Psi_func(time_obs[[i]][[j]], j = j, p = p) %*% Zeta[i, ]
      Y[[i]][[j]] <- as.vector(mean_vec + resid_vec)
    }
  }

  # output the results:

  create_named_list(time_obs, time_g, int_knots, C, C_g, Zeta, mu_g, Psi_g, Y)

}
