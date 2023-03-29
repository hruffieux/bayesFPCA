######### R script: fpca_vmp.R ##########

# For performing Bayesian FPCA via VMP.

# Created: 20 MAR 2023
# Last changed: 21 MAR 2023

# Load libraries:

library(splines)
library(pracma)
library(matrixcalc)
library(lattice)

# Required functions:

setwd("../R")

source("fourier_basis.r")
source("fpca_algs.R")
source("OSullivan_splines.R")
source("plot_functions.R")
source("set_hyper.R")
source("simulation_functions.R")
source("utils.R")
source("vmp_functions.R")
source("wait.r")

setwd("../tests")

# Establish simulation variables:

N <- 100                            # number of curves
n_sample <- 6                       # number of curves for the plots
N_sample <- sort(sample(1:N, n_sample))   # specific curves for the plots
n <- round(runif(N, 20, 30))    # number of time observations for each curve
n_int_knots <- 10                   # number of interior knots
K <- n_int_knots + 2                # number of spline basis functions
L <- 4                              # number of FPCA basis functions
data_col <- "red"                   # colour of the data in the plots
tol <- 1e-5                   # convergence criterion

n_vmp <- 200                        # number of VMP iterations
n_g <- 1000                         # length of the plotting grid

sigma_zeta_vec <- 2/(1:L)           # sd for first and second scores
sigma_eps <- 1                      # sd of the residuals

# Establish hyperparameters:

# sigsq_beta <- 1e10
# Sigma_beta <- sigsq_beta*diag(2)
# mu_beta <- rep(0, 2)
# A <- 1e5
# sigsq_zeta <- 1
# sigma_zeta <- sqrt(sigsq_zeta)
# Sigma_zeta <- sigsq_zeta*diag(L)

# Set the mean function and the FPCA basis functions:

mu <- function(t) return(3*sin(pi*t))
Psi <- fourier_basis(L)

# Generate FPCA data:

fpca_data <- generate_fpca_data(N, 1, n, K, L, n_g, sigma_eps, mu, Psi)

time_obs <- fpca_data$"time_obs"
Zeta <- fpca_data$"Zeta"
mu_g <- fpca_data$"mu_g"
Psi_g <- fpca_data$"Psi_g"
Y <- fpca_data$"Y"

fpca_res <- run_vmp_fpca(time_obs, Y, K, L, list_hyper = NULL, n_g = 1000, time_g = NULL, tol = tol, maxit = n_vmp, plot_elbo = TRUE, Psi_g = Psi_g, verbose = TRUE, seed = NULL)

# subj_names <- paste0("subj_", 1:N)

# grid_obj <- get_grid_objects(time_obs, K, n_g, time_g = NULL, format_univ = TRUE)

# time_g <- grid_obj$"time_g"
# C <- grid_obj$"C"
# C_g <- grid_obj$"C_g"

# # VMP simulations:

# list_hyper <- set_hyper()

# sigma_zeta <- list_hyper$sigma_zeta
# mu_beta <- list_hyper$mu_beta
# Sigma_beta <- list_hyper$Sigma_beta
# A <- list_hyper$A

# eta_vec <- vmp_gauss_fpca(
	# n_vmp, N, L, K, C, Y, sigma_zeta, rep(0, 2),
	# Sigma_beta, A, tol, plot_elbo = TRUE, verbose = TRUE
# )

# # Get the posterior estimates

# eta_in <- list(
	# eta_vec$"p(nu|Sigma_nu)->nu", eta_vec $"p(Y|nu,zeta,sigsq_eps)->nu",
	# eta_vec $"p(zeta)->zeta", eta_vec $"p(Y|nu,zeta,sigsq_eps)->zeta"
# )
# fpc_rotns <- fpc_orthgn(subj_names, L, K, eta_in, time_g, C_g, Psi_g = Psi_g)

# Summarise the VMP results:

time_g <- fpca_res$"time_g"
Y_hat <- fpca_res$"Y_hat"
Y_low <- fpca_res$"Y_low"
Y_upp <- fpca_res$"Y_upp"
mu_hat <- fpca_res$"mu_hat"
Psi_hat <- fpca_res$"list_Psi_hat"
Zeta_hat <- fpca_res$"Zeta_hat"
Cov_zeta_hat <- fpca_res$"Cov_zeta_hat"
list_zeta_ellipse <- fpca_res$"list_zeta_ellipse"

# Y_hat <- fpc_rotns$"Y_hat"
# Y_low <- fpc_rotns$"Y_low"
# Y_upp <- fpc_rotns$"Y_upp"
# mu_hat <- fpc_rotns$"mu_hat"
# Psi_hat <- fpc_rotns$"list_Psi_hat"
# Zeta_hat <- fpc_rotns$"Zeta_hat"
# Cov_zeta_hat <- fpc_rotns$"Cov_zeta_hat"
# list_zeta_ellipse <- fpc_rotns$"list_zeta_ellipse"

# Plot the fitted curves:

display_fit(N_sample, time_obs, time_g, Y, Y_hat, Y_low, Y_upp, offset = 0.1)

wait()

# Plot the mean and basis function estimates:

display_eigenfunctions(L, time_g, mu_g, Psi_g, mu_hat, Psi_hat, format_univ = TRUE)

wait()

# Determine the score accuracies:

plot_scores(N_sample, Zeta, Zeta_hat, list_zeta_ellipse)

############ End of fpca_vmp_12.R ############