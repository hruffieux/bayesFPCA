rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

out_dir <- file.path(CORE_DIR, "output/")

require(bayesFPCA)

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

K <- n_int_knots[1:p] + 2                # number of spline basis functions
L <- 2                              # number of FPCA basis functions

n_mfvb <- 25                       # number of vmp iterations (artificially small for testing purposes)
n_g <- 1000                         # length of the plotting grid

sigma_zeta_vec <- 2*1/(1:L)           # sd for first and second scores
sigma_eps <- rep(1, p)      # sd of the residuals
sigsq_eps <- sigma_eps^2

delta <- 1e-3                       # convergence threshold

sigsq_beta <- 1e10
Sigma_beta <- sigsq_beta*diag(2)
mu_beta <- rep(0, 2)
A <- 1e5
sigsq_zeta <- 1
sigma_zeta <- sqrt(sigsq_zeta)
Sigma_zeta <- sigsq_zeta*diag(L)

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

mfpca_data <- gauss_mfpca_data(N, p, n, L, n_g, sigma_eps, mu, Psi_func)

time_obs <- mfpca_data$"time_obs"
time_g <- mfpca_data$"time_g"
int_knots <- mfpca_data$"int_knots"
C <- mfpca_data$"C"
C_g <- mfpca_data$"C_g"
Zeta <- mfpca_data$"Zeta"
mu_g <- mfpca_data$"mu_g"
Psi_g <- mfpca_data$"Psi_g"
Y <- mfpca_data$"Y"

####################################################
#
#  MFVB  SIMULATIONS
#
####################################################

mfvb_res <- run_mfvb_fpca(
	time_obs, Y, K, L, list_hyper = NULL,
	n_mfvb = n_mfvb, n_g = n_g, time_g = time_g, Psi_g = Psi_g,
	verbose = TRUE, seed = NULL
)

time_g <- mfvb_res$"time_g"
Y_summary <- mfvb_res$"Y_summary"
Y_hat <- mfvb_res$"Y_hat"
Y_low <- mfvb_res$"Y_low"
Y_upp <- mfvb_res$"Y_upp"
gbl_hat <- mfvb_res$"gbl_hat"
mu_hat <- mfvb_res$"mu_hat"
list_Psi_hat <- mfvb_res$"list_Psi_hat"
Zeta_hat <- mfvb_res$"Zeta_hat"
list_zeta_ellipse <- mfvb_res$"list_zeta_ellipse"

wait()

n_sample <- 4                                 # number of curves for the plots
N_sample <- sort(sample(1:N, n_sample))       # specific curves for the plots

display_fit_list(1:p, N_sample, time_obs, time_g, Y,
                             Y_hat, Y_low, Y_upp,
                             Y_hat_add = NULL, Y_low_add = NULL, Y_upp_add = NULL, offset = 0.1,
                             col_data = "grey55", col = "black", col_add = "blue",
                             lwd = 1.2, lwd_add = 1.2)

wait()

# Plot the eigenfunctions:

if(print_pdf) {

	pdf(
		"./images/mfpca_gbl_funcs.pdf",
		width=plot_width, height=plot_height
	)
}

display_eigenfunctions(1:p, L, time_g, mu_g, Psi_g,
                                   mu_hat, list_Psi_hat,
                                   mu_hat_add = NULL, list_Psi_hat_add = NULL,
                                   mu_hat_ci = NULL, list_Psi_hat_ci = NULL,
                                   lwd = 2, data_col = "red",
                                   vec_col_add = NULL, vec_lwd = NULL)

wait()

# Compare first level scores:

display_scores(N_sample, p, Zeta,
                        Zeta_hat, list_zeta_ellipse,
                        Zeta_hat_add = NULL, zeta_ellipse_add = NULL,
                        vec_col = c("black", "blue"), data_col = "red", mfrow = NULL)

############ End of mfpca_vmp.R ############
