rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

out_dir <- file.path(CORE_DIR, "output/")

require(bayesFPCA)

gen_sin_bf <- function(n, p) {

  sin_bf <- function(x) {

    ans <- sqrt(2/p)*sin(2*n*pi*x)
    return(ans)
  }

  return(sin_bf)
}

gen_cos_bf <- function(n, p) {

  cos_bf <- function(x) {

    ans <- sqrt(2/p)*cos(2*n*pi*x)
    return(ans)
  }

  return(cos_bf)
}

fourier_basis <- function(L, p) {

  L_even <- ((L/2) %% 1 == 0)

  if(L_even) {

    n <- 1:(L/2)
    sin_list <- lapply(n, gen_sin_bf, p)
    cos_list <- lapply(n, gen_cos_bf, p)
    fb <- do.call(c, Map(list, sin_list, cos_list))
  } else {

    n_sin <- 1:ceiling(L/2)
    sin_list <- lapply(n_sin, gen_sin_bf, p)

    if(floor(L/2) > 0) {

      n_cos <- 1:floor(L/2)
      cos_list <- lapply(n_cos, gen_cos_bf, p)

      fb <- vector("list", length = L)
      fb[c(TRUE, FALSE)] <- sin_list
      fb[c(FALSE, TRUE)] <- cos_list
    } else {

      fb <- sin_list
    }
  }

  return(fb)
}

Psi_fourier_func <- function(time_obs, j = 1, p = 1) {
  ans <- sapply(fourier_basis(L, p), function(ff) ff(time_obs))
  return(ans)
}


seed <- 1
set.seed(seed)

N <- 100                            # number of curves
n <- round(runif(N, 20, 30))        # number of time observations for each curve
n_int_knots <- 10                   # number of interior knots
K <- n_int_knots + 2                # number of spline basis functions
L <- 4                              # number of FPCA basis functions

tol <- 1e-5                         # convergence criterion
maxit_vmp <- 250                    # number of VMP iterations

sigma_zeta_vec <- 2/(1:L)           # sd for first and second scores
sigma_eps <- 1                      # sd of the residuals

n_g <- 1000                         # length of the plotting grid

# format univ = T
# ---------------
fpca_data <- generate_fpca_data(N, p = 1, n, L, n_g, sigma_eps, mu_func,
                                Psi_fourier_func, vec_sd_zeta = sigma_zeta_vec,
                                format_univ = T)

time_obs <- fpca_data$time_obs
Zeta <- fpca_data$Zeta
mu_g <- fpca_data$mu_g
Psi_g <- fpca_data$Psi_g
Y <- fpca_data$Y

fpca_res <- run_vmp_fpca(time_obs, Y, L, K = K, n_g = n_g, tol = tol,
                         maxit = maxit_vmp, Psi_g = Psi_g)

time_g <- fpca_res$time_g
Y_hat <- fpca_res$Y_hat
Y_low <- fpca_res$Y_low
Y_upp <- fpca_res$Y_upp
mu_hat <- fpca_res$mu_hat
Psi_hat <- fpca_res$list_Psi_hat
Zeta_hat <- fpca_res$Zeta_hat
list_zeta_ellipse <- fpca_res$list_zeta_ellipse

n_sample <- 6                             # number of curves for the plots
N_sample <- sort(sample(1:N, n_sample))   # specific curves for the plots

display_fit(N_sample, time_obs, time_g, Y, Y_hat, Y_low, Y_upp)

display_eigenfunctions(L, time_g, mu_g, Psi_g, mu_hat, Psi_hat)

display_scores(N_sample, Zeta, Zeta_hat, list_zeta_ellipse)


# format univ = F
# ---------------
n <- matrix(n, nrow = N)        # number of time observations for each curve

fpca_data <- generate_fpca_data(N, p = 1, n, L, n_g, sigma_eps, mu_func,
                                Psi_fourier_func, vec_sd_zeta = sigma_zeta_vec,
                                format_univ = F)

time_obs <- fpca_data$time_obs
Zeta <- fpca_data$Zeta
mu_g <- fpca_data$mu_g
Psi_g <- fpca_data$Psi_g
Y <- fpca_data$Y

fpca_res <- run_vmp_fpca(time_obs, Y, L, K = K, n_g = n_g, tol = tol,
                         maxit = maxit_vmp, Psi_g = Psi_g)

time_g <- fpca_res$time_g
Y_hat <- fpca_res$Y_hat
Y_low <- fpca_res$Y_low
Y_upp <- fpca_res$Y_upp
mu_hat <- fpca_res$mu_hat
Psi_hat <- fpca_res$list_Psi_hat
Zeta_hat <- fpca_res$Zeta_hat
list_zeta_ellipse <- fpca_res$list_zeta_ellipse

n_sample <- 6                             # number of curves for the plots
N_sample <- sort(sample(1:N, n_sample))   # specific curves for the plots

display_fit_list(p_sample = 1, N_sample, time_obs, time_g, Y, Y_hat, Y_low, Y_upp)

display_eigenfunctions(L, time_g, mu_g, Psi_g, mu_hat, Psi_hat)

if (L > 1) { # displays the scores for the first two components
  display_scores(N_sample, Zeta, Zeta_hat, list_zeta_ellipse)
}


run_model_choice_version <- T
if (run_model_choice_version) {

  n_cpus <- 1
  fpca_res_mc_for_K <- run_vmp_fpca_model_choice(time_obs, Y, L, n_g = n_g,
                                                 tol = tol, maxit = maxit_vmp,
                                                 Psi_g = Psi_g, verbose = F,
                                                 n_cpus = n_cpus)

  # selected K
  fpca_res_mc_for_K$K

  time_g_mc <- fpca_res_mc_for_K$time_g
  Y_hat_mc <- fpca_res_mc_for_K$Y_hat
  Y_low_mc <- fpca_res_mc_for_K$Y_low
  Y_upp_mc <- fpca_res_mc_for_K$Y_upp
  mu_hat_mc <- fpca_res_mc_for_K$mu_hat
  Psi_hat_mc <- fpca_res_mc_for_K$list_Psi_hat
  Zeta_hat_mc <- fpca_res_mc_for_K$Zeta_hat
  list_zeta_ellipse_mc <- fpca_res_mc_for_K$list_zeta_ellipse

  display_eigenfunctions(L, time_g, mu_g, Psi_g, mu_hat, Psi_hat,
                         mu_hat_add = mu_hat_mc, list_Psi_hat_add = Psi_hat_mc)

  if (L > 1) { # displays the scores for the first two components
    display_scores(N_sample, Zeta, Zeta_hat, list_zeta_ellipse,
                   Zeta_hat_add = Zeta_hat_mc, zeta_ellipse_add = list_zeta_ellipse_mc)
  }

  # Compare errors with and without model choice for K
  #
  rmse_mfpca <- apply(Zeta - Zeta_hat, 2, function(x) sqrt(mean(x^2)))
  rmse_mfpca_mc <- apply(Zeta - Zeta_hat_mc, 2, function(x) sqrt(mean(x^2)))

  print(rmse_mfpca)
  print(rmse_mfpca_mc)

  ise_mfpca <- ise_mfpca_mc <- rep(NA, L+1)

  for (l in 1:(L+1)) {

    if (l == 1) {
      ise_mfpca[l] <- trapint(time_g, (mu_hat - mu_g[[1]])^2)
      ise_mfpca_mc[l] <- trapint(time_g, (mu_hat_mc - mu_g[[1]])^2)
    } else {
      ise_mfpca[l] <- trapint(time_g, (Psi_hat[[l-1]] - Psi_g[[1]][,l-1])^2)
      ise_mfpca_mc[l] <- trapint(time_g, (Psi_hat_mc[[l-1]] - Psi_g[[1]][,l-1])^2)
    }

  }

  print(ise_mfpca)
  print(ise_mfpca_mc)
}
