rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

out_dir <- file.path(CORE_DIR, "output/")

library(bayesFPCA)

seed <- 1
set.seed(seed)

# Establish simulation variables:

p <- 3                                        # number of responses
N <- 100                                      # number of curves
N_t_min <- 10
N_t_max <- 20
n <- matrix(sample(N_t_min:N_t_max, N*p,
                   replace = TRUE), N, p)     # number of time observations

n_int_knots <- rep(8,p)                              # number of interior knots
K <- n_int_knots + 2                          # number of spline basis functions # if length(K) = 1, then will be set to K <- rep(K, p) within the run_* functions
L <- 2                                        # number of FPCA basis functions

tol  <- 1e-5                                  # convergence tolerance # 1e-3 is not enough! (doesn't match the mfvb run)
maxit_vmp <- 25                              # maximum number of vmp iterations (artificially low for test purpose only)
n_mfvb <- 25                                 # number of mfvb iterations (no convergence criterion implemented yet for mfvb - artificially low for test purpose only)

n_g <- 1000                                   # length of the plotting grid

bool_multiple_repl_eigen <- FALSE             # do multiple replicate for
                                              # estimation for eigenfunctions
if (bool_multiple_repl_eigen) {
  n_repl_eigen <- 4
  n_cpus_eigen <- 4
}

# Parameters for data generation
#
sigma_eps <- rep(1, p)                        # sd of the residuals
sigsq_eps <- sigma_eps^2

# Set the mean function and the FPCA basis functions:

mu_func <- function(t, j) (-1)^j*2*sin((2*pi+j)*t)

psi_1 <- function(t, j, p) (-1)^j * sqrt(2/p)*cos(2*pi*t)
psi_2 <- function(t, j, p) (-1)^j * sqrt(2/p)*sin(2*pi*t)


Psi_func <- function(time_obs, j, p) {
  ans <- cbind(psi_1(time_obs, j, p), psi_2(time_obs, j, p))
  return(ans)
}


format_univ <- FALSE # can be faster if TRUE (although not clear!)
                     # when p = 1 but then the plot functions need to be adapeted
if (format_univ) {
  if (p == 1) {
    n <- as.vector(n)
  } else {
      stop("format_univ must be FALSE if multivariate case.")
  }
}

bool_save <- F
if (bool_save) {
  res_dir <- paste0(out_dir, "/simulations_p_", p, "_N_", N, "_Nt_min_",
                    N_t_min, "_max_", N_t_max, "_K_", K, "_L_", L, "_sigeps_",
                    unique(sigma_eps), "_tol_", tol, "_seed_", seed, "/")
  dir.create(res_dir)
  sink(paste(res_dir, "out.txt", sep = ""), append = F, split = T, type = "output")
  sink(file(paste(res_dir, "err.txt", sep = ""), open = "wt"), type = "message")
} else {
  res_dir <- NULL
}


####################################################
#
#  SIMULATE  THE  DATA
#
####################################################

# Set up the data:
#
data <- generate_fpca_data(N, p, n, K, L, n_g, sigma_eps,
                           mu_func, Psi_func, format_univ = format_univ)

time_obs <- data$time_obs
Zeta <- data$Zeta
Y <- data$Y

time_g <- data$time_g
mu_g <- data$mu_g
Psi_g <- data$Psi_g


####################################################
#
#  VMP  SIMULATIONS
#
####################################################

set.seed(seed)

# doesn't run successfully at the moment as treatment of K and corrections in fpca_rotation and mfpca_rotation functions need to be implemented
#
system.time(vmp_res <- run_vmp_fpca(time_obs, Y, K, L, n_g = NULL, time_g = time_g, # here we use the same time grid as that used to simulate the true mean and eigen- functions
                                    tol = tol, maxit = maxit_vmp,
                                    plot_elbo = TRUE, Psi_g = Psi_g))

Y_hat <- vmp_res$Y_hat
Y_low <- vmp_res$Y_low
Y_upp <- vmp_res$Y_upp
mu_hat <- vmp_res$mu_hat
list_Psi_hat <- vmp_res$list_Psi_hat
Zeta_hat <- vmp_res$Zeta_hat
list_zeta_ellipse <- vmp_res$list_zeta_ellipse

n_sample <- 4                                 # number of curves for the plots
N_sample <- sort(sample(1:N, n_sample))       # specific curves for the plots


if (p == 1 & format_univ) {

  if (bool_save) {
    pdf(paste0(res_dir, "/data_samples_", paste0(N_sample, collapse = "-"), ".pdf"),
        width = 6.3, height = 5, paper='special')
  }
  display_fit(N_sample, time_obs, time_g, Y, Y_hat = NULL,
              Y_low = NULL, Y_upp = NULL)
  if (bool_save) {
    dev.off()
    pdf(paste0(res_dir, "/fits_samples_", paste0(N_sample, collapse = "-"), ".pdf"),
        width = 6.3, height = 5, paper='special')
  }
  display_fit(N_sample, time_obs, time_g, Y, Y_hat, Y_low, Y_upp)
  if (bool_save) {
    dev.off()
  }

} else {

  p_sample <- 1:p                               # variables to display

  if (bool_save) {
    pdf(paste0(res_dir, "/data_samples_", paste0(N_sample, collapse = "-"), ".pdf"),
        width = 6.3, height = 5, paper='special')
  }
  display_fit_list(p_sample, N_sample, time_obs, time_g, Y, Y_hat = NULL,
                   Y_low = NULL, Y_upp = NULL)
  if (bool_save) {
    dev.off()
    pdf(paste0(res_dir, "/fits_vmp_only_samples_", paste0(N_sample, collapse = "-"), ".pdf"),
        width = 6.3, height = 5, paper='special')
  }
  display_fit_list(p_sample, N_sample, time_obs, time_g, Y, Y_hat, Y_low, Y_upp)
  if (bool_save) {
    dev.off()
  }
}



####################################################
#
#  MFVB  SIMULATIONS
#
####################################################

if (!format_univ) { # not implemented for format univ

  set.seed(seed)

  mfvb_res <- run_mfvb_fpca(time_obs, Y, K, L, n_mfvb = n_mfvb, n_g = NULL,
                            time_g = time_g, # here we use the same time grid as that used to simulate the true mean and eigen- functions
                            Psi_g = Psi_g)

  ####################################################
  #
  #  COMPARISONS
  #
  ####################################################

  if (bool_save) {
    pdf(paste0(res_dir, "/fits_vmp_and_mfvb_samples_", paste0(N_sample, collapse = "-"), ".pdf"),
        width = 6.3, height = 5, paper='special')
  }
  display_fit_list(p_sample, N_sample, time_obs, time_g,
                   Y, Y_hat, Y_low, Y_upp, # vmp
                   Y_hat_add = mfvb_res$Y_hat, # mfvb
                   Y_low_add = mfvb_res$Y_low,
                   Y_upp_add = mfvb_res$Y_upp)

  if (bool_save) {
    dev.off()
    pdf(paste0(res_dir, "/eigenfunction_vmp_and_mfvb.pdf"),
          width = 6.3, height = 5, paper='special')
  }

  display_eigenfunctions(p_sample, L, time_g, mu_g, Psi_g,
                         mu_hat, list_Psi_hat, # vmp = grey
                         mu_hat_add = mfvb_res$mu_hat, # mfvb = blue
                         list_Psi_hat_add = mfvb_res$list_Psi_hat,
                         vec_lwd = c(3,1))
  if (bool_save) {
    dev.off()
  }
}

n_plots_scores <- 4
n_sample_ps <- 16

for (ps in 1:n_plots_scores) {

  set.seed(ps)
  N_sample_ps <- sort(sample(1:N, n_sample_ps))

  if (bool_save) {
    pdf(paste0(res_dir, "/scores_samples_", paste0(N_sample_ps, collapse = "-"), ".pdf"),
        width = 6, height = 6, paper='special')
  }
  plot_scores(N_sample_ps, p, Zeta, Zeta_hat, list_zeta_ellipse,
              Zeta_hat_add = mfvb_res$Zeta_hat,
              zeta_ellipse_add = mfvb_res$list_zeta_ellipse,
              mfrow = c(ceiling(sqrt(n_sample_ps)), tail(sqrt(n_sample_ps))))
  if (bool_save) {
    dev.off()
  }
}


# Show multiple replicates for estimation of eigenfunctions
#
if (bool_multiple_repl_eigen & !format_univ) {

  list_res <- parallel::mclapply(1:n_repl_eigen, function(repl) {

    set.seed(repl)

    data_repl <- generate_fpca_data(N, p, n, K, L, n_g, sigma_eps, mu_func,
                                    Psi_func, time_obs = time_obs, # don't re-generate the observations
                                    format_univ = format_univ)

    Y_repl <- data_repl$Y


    vmp_res_repl <- run_vmp_fpca(time_obs, Y_repl, K, L, n_g = NULL, time_g = time_g, # here we use the same time grid as that used to simulate the true mean and eigen- functions
                                 tol = tol, maxit = maxit_vmp, plot_elbo = TRUE,
                                 Psi_g = Psi_g)

    mu_hat <- vmp_res_repl$mu_hat
    list_Psi_hat <- vmp_res_repl$list_Psi_hat

    create_named_list(mu_hat, list_Psi_hat)
  }, mc.cores = n_cpus_eigen)

  list_mu_hat <- lapply(list_res, "[[", "mu_hat")
  list_list_Psi_hat <- lapply(list_res, "[[", "list_Psi_hat")


  if (bool_save) {
    pdf(paste0(res_dir, "/eigenfunction_n_repl_", n_repl_eigen, ".pdf"),
        width = 6.3, height = 5, paper='special')
  }
  display_eigenfunctions(p_sample, L, time_g, mu_g, Psi_g, list_mu_hat,
                         list_list_Psi_hat, lwd = 1)

  if (bool_save) {
    dev.off()
  }

}


# Compare first level scores:

norm_diff <- apply(mfvb_res$Zeta - Zeta, 1, function(x) sqrt(cprod(x)))
mfvb_rmse <- sqrt(mean(norm_diff))

cat("The MFVB-based rmse is:", mfvb_rmse, "\n")

norm_diff <- apply(vmp_res$Zeta - Zeta, 1, function(x) sqrt(cprod(x)))
vmp_rmse <- sqrt(mean(norm_diff))

cat("The VMP-based rmse is:", vmp_rmse, "\n")

