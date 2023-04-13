rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

out_dir <- file.path(CORE_DIR, "output/")

library(bayesFPCA)

# Establish simulation variables:

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
n_sample <- 4                                 # number of curves for the plots
N_sample <- sort(sample(1:N, n_sample))       # specific curves for the plots
K <- n_int_knots[1:p] + 2                # number of spline basis functions
L <- 2                              # number of FPCA basis functions

n_vmp <- 200                       # number of vmp iterations
n_g <- 1000                         # length of the plotting grid

sigma_zeta_vec <- 2/(1:L)           # sd for first and second scores
sigma_eps <- rep(1, p)      # sd of the residuals

tol <- 1e-5                       # convergence threshold

# Set the mean function and the FPCA basis functions:

mu <- function(time_obs, j) {

	ans <- ((-1)^j)*3*sin(pi*j*time_obs) - 3/2
	return(ans)
}

psi_1 <- function(time_obs, j) {

	ans <- sqrt(2/p)*sin(2*pi*j*time_obs)
	return(ans)
}

psi_2 <- function(time_obs, j) {

	ans <- sqrt(2/p)*cos(2*pi*j*time_obs)
	return(ans)
}

Psi_func <- function(time_obs, j, p) {

  ans <- cbind(psi_1(time_obs, j), psi_2(time_obs, j))
  return(ans)
}

####################################################
#
#  SIMULATE  THE  DATA
#
####################################################

# Set up the data:

mfpca_data <- generate_fpca_data(N, p, n, K, L, n_g, sigma_eps, mu, Psi_func)

time_obs <- mfpca_data$"time_obs"
Zeta <- mfpca_data$"Zeta"
mu_g <- mfpca_data$"mu_g"
Psi_g <- mfpca_data$"Psi_g"
Y <- mfpca_data$"Y"

####################################################
#
#  VMP  SIMULATIONS
#
####################################################

# list_hyper <- set_hyper()

# sigma_zeta <- list_hyper$"sigma_zeta"
# mu_beta <- list_hyper$"mu_beta"
# Sigma_beta <- list_hyper$"Sigma_beta"
# A <- list_hyper$"A"

# subj_names <- paste0("subj_", 1:N)
# names(Y) <- names(time_obs) <- subj_names

# grid_obj <- get_grid_objects(time_obs, K, n_g = n_g, time_g = NULL, format_univ = FALSE)

# C <- grid_obj$C
# n_g <- grid_obj$n_g
# time_g <- grid_obj$time_g
# C_g <- grid_obj$C_g

# p <- length(Y[[1]])

# tmp_names <- lapply(Y, function(Y_i) names(Y_i))
# tmp_names_time <-  lapply(time_obs, function(time_obs_i) names(time_obs_i))
# if (is.null(unlist(unique(tmp_names)))) {
  # stopifnot(is.null(unlist(unique(tmp_names_time))))
  # var_names <- paste0("variable_", 1:p)
  # Y <- lapply(Y, function(Y_i) { names(Y_i) <- var_names; Y_i})
# } else if (!all_same(tmp_names)) {
  # stop("Variable names in Y must be the same for all individuals.")
# } else {
  # var_names <- names(Y[[1]])
# }

# time_obs <- lapply(time_obs, function(time_obs_i) { names(time_obs_i) <- var_names; time_obs_i})

# eta_vec <- vmp_gauss_mfpca(n_vmp, N, p, L, K, C, Y, mu_beta,
                               # Sigma_beta, A, tol,
                               # plot_elbo = TRUE, verbose = TRUE)

# eta_in <- list(
  # eta_vec$"p(nu|Sigma_nu)->nu", eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
  # eta_vec$"p(zeta)->zeta", eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
# )

# mfpca_res <- mfpc_orthgn(subj_names, L, K, eta_in, time_g, N, p, C_g, Psi_g)

mfpca_res <- run_vmp_fpca(
	time_obs, Y, K, L, tol = tol, maxit = n_vmp,
	plot_elbo = TRUE, Psi_g = Psi_g
)

time_g <- mfpca_res$"time_g"
Y_hat <- mfpca_res$"Y_hat"
Y_low <- mfpca_res$"Y_low"
Y_upp <- mfpca_res$"Y_upp"
mu_hat <- mfpca_res$"mu_hat"
list_Psi_hat <- mfpca_res$"list_Psi_hat"
Zeta_hat <- mfpca_res$"Zeta_hat"
Cov_zeta_hat <- mfpca_res$"Cov_zeta_hat"
list_zeta_ellipse <- mfpca_res$"list_zeta_ellipse"

####################################################
#
#  COMPARISONS
#
####################################################

# Plot the fits:

display_fit_list(1:p, N_sample, time_obs, time_g, Y, Y_hat, Y_low, Y_upp)

wait()

# Plot the eigenfunctions:

display_eigenfunctions(L, time_g, mu_g, Psi_g, mu_hat, list_Psi_hat, p_sample = 1:p)

wait()

# Determine the score accuracies:

plot_scores(N_sample, Zeta, Zeta_hat, list_zeta_ellipse)

############ End of mfpca_vmp.R ############
