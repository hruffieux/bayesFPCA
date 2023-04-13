rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

out_dir <- file.path(CORE_DIR, "output/")

library(bayesFPCA)

# Establish simulation variables:

N <- 100                            # number of curves
n <- round(runif(N, 20, 30))        # number of time observations for each curve
n_int_knots <- 10                   # number of interior knots
K <- n_int_knots + 2                # number of spline basis functions
L <- 4                              # number of FPCA basis functions
tol <- 1e-5                         # convergence criterion

maxit_vmp <- 25                     # number of VMP iterations - very small for testing purposes only
n_g <- 1000                         # length of the plotting grid

sigma_zeta_vec <- 2/(1:L)           # sd for first and second scores
sigma_eps <- 1                      # sd of the residuals


fpca_data <- generate_fpca_data(N, p = 1, n, L, n_g, sigma_eps, mu_func,
                                Psi_fourier_func, format_univ = TRUE)

time_obs <- fpca_data$time_obs
Zeta <- fpca_data$Zeta
mu_g <- fpca_data$mu_g
Psi_g <- fpca_data$Psi_g
Y <- fpca_data$Y

fpca_res <- run_vmp_fpca(time_obs, Y, K, L, n_g = n_g, tol = tol,
                         maxit = maxit_vmp, Psi_g = Psi_g)

time_g <- fpca_res$time_g
Y_hat <- fpca_res$Y_hat
Y_low <- fpca_res$Y_low
Y_upp <- fpca_res$Y_upp
mu_hat <- fpca_res$mu_hat
Psi_hat <- fpca_res$list_Psi_hat
Zeta_hat <- fpca_res$Zeta_hat
Cov_zeta_hat <- fpca_res$Cov_zeta_hat
list_zeta_ellipse <- fpca_res$list_zeta_ellipse


n_sample <- 6                             # number of curves for the plots
N_sample <- sort(sample(1:N, n_sample))   # specific curves for the plots

display_fit(N_sample, time_obs, time_g, Y, Y_hat, Y_low, Y_upp)

display_eigenfunctions(L, time_g, mu_g, Psi_g, mu_hat, Psi_hat)

display_scores(N_sample, Zeta, Zeta_hat, list_zeta_ellipse)
