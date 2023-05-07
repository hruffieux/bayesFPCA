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

n_mfvb <- 25                    # number of MFVB iterations

sigma_zeta_vec <- 2/(1:L)           # sd for first and second scores
sigma_eps <- 1                      # sd of the residuals

n_g <- 1000                         # length of the plotting grid

# format univ = T # for now, the only possibility for MFVB version, = F not implemented
# ---------------

# mu_func_zero <- function(t, j = 1, alpha = 3) rep(0, length(t))
fpca_data <- generate_fpca_data(N, p = 1, n, L, n_g, sigma_eps, mu_func, # mu_func_zero,
                                Psi_fourier_func, vec_sd_zeta = sigma_zeta_vec,
                                format_univ = T)

time_obs <- fpca_data$time_obs
Zeta <- fpca_data$Zeta
mu_g <- fpca_data$mu_g
Psi_g <- fpca_data$Psi_g
Y <- fpca_data$Y


fpca_res <- run_mfvb_fpca(time_obs, Y, L, K = K, n_mfvb = n_mfvb, n_g = n_g,
                          Psi_g = Psi_g)


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

