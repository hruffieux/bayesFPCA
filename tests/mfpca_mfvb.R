rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

out_dir <- file.path(CORE_DIR, "output/")

require(bayesFPCA)

seed <- 1
set.seed(seed)

n_obs <- list(10:20, 50:70, 30:40, 100:120)
n_int_knots <- c(5, 15, 9, 25)                # number of interior knots

p <- 3                                        # number of responses
N <- 100                                      # number of curves
n <- Reduce(                                  # number of time observations
  cbind, lapply(
    n_obs[1:p], sample,
    size = N, replace = TRUE
  )
)
K <- n_int_knots[1:p] + 2                     # number of spline basis functions
L <- 2                                        # number of FPCA basis functions

n_mfvb <- 25                              # maximum number of vmp iterations (artificially low for test purpose only)

sigma_zeta_vec <- 2/(1:L)                     # sd for first and second scores
sigma_eps <- rep(1, p)                        # sd of the residuals

n_g <- 1000                                   # length of the plotting grid


mu_v2_func <- function(time_obs, j) {

  ans <- ((-1)^j)*3*sin(pi*j*time_obs) - 3/2
  return(ans)
}

psi_1_v2 <- function(time_obs, j) {

  ans <- sqrt(2/p)*sin(2*pi*j*time_obs)
  return(ans)
}

psi_2_v2 <- function(time_obs, j) {

  ans <- sqrt(2/p)*cos(2*pi*j*time_obs)
  return(ans)
}

Psi_v2_func <- function(time_obs, j, p) {

  ans <- cbind(psi_1_v2(time_obs, j), psi_2_v2(time_obs, j))
  return(ans)
}

mfpca_data <- generate_fpca_data(N, p, n, L, n_g, sigma_eps, mu_v2_func,
                                 Psi_v2_func, vec_sd_zeta = sigma_zeta_vec)

time_obs <- mfpca_data$time_obs
Zeta <- mfpca_data$Zeta
mu_g <- mfpca_data$mu_g
Psi_g <- mfpca_data$Psi_g
Y <- mfpca_data$Y


mfpca_res <- run_mfvb_fpca(time_obs, Y, L, K = K, n_mfvb = n_mfvb, n_g = n_g,
                           Psi_g = Psi_g)

time_g <- mfpca_res$time_g
Y_hat <- mfpca_res$Y_hat
Y_low <- mfpca_res$Y_low
Y_upp <- mfpca_res$Y_upp
mu_hat <- mfpca_res$mu_hat
list_Psi_hat <- mfpca_res$list_Psi_hat
Zeta_hat <- mfpca_res$Zeta_hat
list_zeta_ellipse <- mfpca_res$list_zeta_ellipse


n_sample <- 4                                 # number of curves for the plots
N_sample <- sort(sample(1:N, n_sample))       # specific curves for the plots

display_fit_list(1:p, N_sample, time_obs, time_g, Y, Y_hat, Y_low, Y_upp)

display_eigenfunctions(L, time_g, mu_g, Psi_g, mu_hat, list_Psi_hat)

display_scores(N_sample, Zeta, Zeta_hat, list_zeta_ellipse)

