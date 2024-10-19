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
n <- matrix(Reduce(                                  # number of time observations
  cbind, lapply(
    n_obs[1:p], sample,
    size = N, replace = TRUE
  )
), nrow = N)
K <- n_int_knots[1:p] + 2                     # number of spline basis functions
L <- 2                                        # number of FPCA basis functions

tol <- 1e-3                                   # convergence threshold
maxit_vmp <- 250                              # maximum number of vmp iterations

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

mfpca_res <- run_vmp_fpca(time_obs, Y, L, K = K, tol = tol, maxit = maxit_vmp,
                          n_g = n_g, Psi_g = Psi_g)

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

if (L > 1) { # displays the scores for the first two components
  display_scores(N_sample, Zeta, Zeta_hat, list_zeta_ellipse)
}


run_model_choice_version <- F
if (run_model_choice_version) {

  model_choice <- "L"
  n_cpus <- 1
  if (model_choice == "K") {

    K_min <- 5
    K_max <- 10

    fpca_res_mc <- run_vmp_fpca_model_choice(time_obs, Y,
                                             model_choice = model_choice,
                                             list_L = list("L" = L),
                                             list_K = list("K_min" = K_min,
                                                           "K_max" = K_max),
                                             n_g = n_g,
                                             tol = tol, maxit = maxit_vmp,
                                             Psi_g = Psi_g, verbose = F,
                                             n_cpus = n_cpus)

    # selected K
    fpca_res_mc$K

  } else {

    lambda_L <- 3
    L_max <- 10

    fpca_res_mc <- run_vmp_fpca_model_choice(time_obs, Y,
                                             model_choice = model_choice,
                                             list_L = list("lambda_L" = lambda_L, "L_max" = L_max),
                                             list_K = list("K" = K),
                                             n_g = n_g,
                                             tol = tol, maxit = maxit_vmp,
                                             Psi_g = Psi_g, verbose = F,
                                             n_cpus = n_cpus)

    # selected L
    fpca_res_mc$L

  }

  time_g_mc <- fpca_res_mc$time_g
  Y_hat_mc <- fpca_res_mc$Y_hat
  Y_low_mc <- fpca_res_mc$Y_low
  Y_upp_mc <- fpca_res_mc$Y_upp
  mu_hat_mc <- fpca_res_mc$mu_hat
  list_Psi_hat_mc <- fpca_res_mc$list_Psi_hat
  Zeta_hat_mc <- fpca_res_mc$Zeta_hat
  list_zeta_ellipse_mc <- fpca_res_mc$list_zeta_ellipse

  display_eigenfunctions(L, time_g, mu_g, Psi_g, mu_hat, list_Psi_hat,
                         mu_hat_add = mu_hat_mc, list_Psi_hat_add = list_Psi_hat_mc)

  if (L > 1) { # displays the scores for the first two components
    display_scores(N_sample, Zeta, Zeta_hat, list_zeta_ellipse,
                   Zeta_hat_add = Zeta_hat_mc, zeta_ellipse_add = list_zeta_ellipse_mc)
  }


  # Compare errors with and without model choice for K
  #
  rmse_mfpca <- apply(Zeta - Zeta_hat, 2, function(x) sqrt(mean(x^2)))
  rmse_mfpca_mc <- apply(Zeta - Zeta_hat_mc[,1:L, drop = F], 2, function(x) sqrt(mean(x^2)))

  print(rmse_mfpca)
  print(rmse_mfpca_mc)

  ise_mfpca <- ise_mfpca_mc <- matrix(NA, L+1, p)

  for (l in 1:(L+1)) {

    for(j in 1:p) {

      if (l == 1) {
        ise_mfpca[l, j] <- trapint(time_g, (mu_hat[,j] - mu_g[[j]])^2)
        ise_mfpca_mc[l, j] <- trapint(time_g, (mu_hat_mc[,j] - mu_g[[j]])^2)
      } else {
        ise_mfpca[l, j] <- trapint(time_g, (list_Psi_hat[[l-1]][,j] - Psi_g[[j]][,l-1])^2)
        ise_mfpca_mc[l, j] <- trapint(time_g, (list_Psi_hat_mc[[l-1]][,j] - Psi_g[[j]][,l-1])^2)
      }

    }

  }
  rownames(ise_mfpca) <- rownames(ise_mfpca_mc) <- c("mu", paste0("Psi_", 1:L))
  colnames(ise_mfpca) <- colnames(ise_mfpca_mc) <- paste0("Variable_", 1:p)

  print(ise_mfpca)
  print(ise_mfpca_mc)
}
