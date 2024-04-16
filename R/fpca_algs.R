#' Variational Message Passing (VMP) inference for univariate or multivariate
#' functional principal component analysis.
#'
#' This function is used to perform FPCA or mFPCA using a VMP algorithm.
#'
#' @param time_obs List of vectors or list of lists of vectors containing time
#'                 of observations for univariate or multivariate curves,
#'                 respectively.
#' @param Y List of vectors (for univariate curves) or list of lists of vectors
#'          (for multivariate curves) with the functional data.
#' @param L Number of eigenfunctions (or latent dimensions).
#' @param K Number of O'Sulivan spline functions to be used in the FPCA or mFPCA
#'          algorithms. If set to \code{NULL} will be set according to the rule
#'          of Ruppert (2002), also enforcing K >=7.
#' @param list_hyper Hyperparameter settings constructed using the function
#'         \code{\link{set_hyper}}. If \code{NULL} default hyperparameters will
#'         be used.
#' @param n_g Desired size for dense grid.
#' @param time_g Dense grid provided as a vector of size \code{n_g}. If provided,
#'               then \code{n_g} must be \code{NULL} as will be taken to be
#'               \code{length(time_g)}.
#' @param tol Tolerance on the relative changes in the ELBO.
#' @param maxit Maximum number of iterations allowed.
#' @param plot_elbo Boolean indicating whether the values of the ELBO should be
#'                  displayed during the run.
#' @param Psi_g Reference eigenfunctions (if available, e.g., in simulations)
#'              used to flip the sign of the resulting scores and eigenfunctions.
#' @param verbose Boolean indicating whether messages should be printed during
#'                the run.
#' @param seed User-specified seed for reproducibilty.
#'
#' @return An object containing the resulting VMP estimates.
#'
#' @seealso  \code{\link{run_mfvb_fpca}}, \code{\link{set_hyper}},
#'           \code{\link{flip_sign}}, \code{\link{display_fit}},
#'           \code{\link{display_fit_list}}, \code{\link{display_scores}},
#'           \code{\link{display_eigenfunctions}}
#'
#' @references
#'
#' Nolan, T. H., Goldsmith, J., & Ruppert, D. (2023). Bayesian Functional
#' Principal Components Analysis via Variational Message Passing with Multilevel
#' Extensions. Bayesian Analysis, 1(1), 1-27.
#'
#' Ruppert, D. (2002). Selecting the number of knots for penalized splines.
#' Journal of computational and graphical statistics, 11(4), 735-757.
#'
#' @export
#'
run_vmp_fpca <- function(time_obs, Y, L, K = NULL,
                         list_hyper = NULL,
                         n_g = 1000, time_g = NULL,
                         tol = 1e-5, maxit = 1e4, plot_elbo = FALSE,
                         Psi_g = NULL, verbose = TRUE, seed = NULL) {

  check_structure(seed, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(seed)) {
    cat(paste0("== Seed set to ", seed, " ==\n\n"))
    set.seed(seed)
  }

  stopifnot(is.list(time_obs))
  stopifnot(is.list(Y))
  stopifnot(isTRUE(all.equal(length(time_obs), length(Y))))

  check_structure(L, "vector", "numeric", 1)
  check_natural(L)

  if (is.null(list_hyper)) {
    list_hyper <- set_hyper()
  } else if (!inherits(list_hyper, "hyper")) {
    stop(paste0("The provided list_hyper must be an object of class ",
                "``hyper''. \n *** you must either use the ",
                "function set_hyper to set your own hyperparameters or ",
                "list_hyper to NULL for automatic choice. ***"))
  }

  sigma_zeta <- list_hyper$sigma_zeta
  # mu_beta <- list_hyper$mu_beta # fixed to zero in mfvb anyway
  mu_beta <- rep(0, 2)
  Sigma_beta <- list_hyper$Sigma_beta
  A <- list_hyper$A

  check_structure(n_g, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(n_g)) check_natural(n_g)

  check_structure(time_g, "vector", "numeric", null_ok = TRUE)

  check_structure(tol, "vector", "numeric", 1)
  check_positive(tol, eps=.Machine$double.eps)

  check_structure(maxit, "vector", "numeric", 1)
  check_natural(maxit)

  check_structure(plot_elbo, "vector", "logical", 1)

  check_structure(verbose, "vector", "logical", 1)

  N <- length(time_obs) # move to vmp_gauss_fpca and vmp_gauss_mfpca

  if (is.null(names(Y))) {
    stopifnot(is.null(names(time_obs)))
    subj_names <- paste0("subj_", 1:N)
    names(Y) <- names(time_obs) <- subj_names
  } else {
    subj_names <- names(Y)
    stopifnot(isTRUE(all.equal(names(time_obs), subj_names)) | is.null(names(time_obs)))
  }
  names(time_obs) <- subj_names

  format_univ <- ifelse(is.list(time_obs[[1]]), FALSE, TRUE) # If TRUE, then p = 1, and vmp_gauss_fpca is used.
  # Note that if FALSE and p = 1, vmp_gauss_mfpca will be used,
  # i.e., special case of multivariate algorithm for p=1 gives the univariate model
  # (produces inference equivalent to vmp_gauss_fpca, the difference is the format of the input)

  if (!format_univ) {
    p <- length(time_obs[[1]])
    if (!is.null(Psi_g)) {
      stopifnot(is.list(Psi_g))
    }
  } else {
    p <- 1
  }

  if (is.null(K)) {  # Ruppert (2002) sets a simple default value for K as min(nobs/4,40), where nobs is the number of observations.
    # here since nobs differs for each i, we take the nobs / 4 = round(median(obs_i)/4), and do this for each variable j = 1, ..., p
    # and we enforce that K>=7

    K <- sapply(1:p, function(j) max(round(min(median(sapply(time_obs, function(time_obs_i) length(time_obs_i[[j]]))/4), 40)), 7))

    # if supplied K is such that length(K) = 1, then will be set to K <- rep(K, p)
  } else if (!format_univ) {
    check_structure(K, "vector", "numeric", c(1,p))
    if (length(K)==1) K <- rep(K, p)
  } else {
    check_structure(K, "vector", "numeric", 1)
  }
  check_natural(K)

  grid_obj <- get_grid_objects(time_obs, K, n_g = n_g, time_g = time_g,
                               format_univ = format_univ)

  C <- grid_obj$C
  n_g <- grid_obj$n_g
  time_g <- grid_obj$time_g
  C_g <- grid_obj$C_g

  if (format_univ) {

    eta_vec <- vmp_gauss_fpca(n_vmp = maxit, N, L, K, C, Y, sigma_zeta, mu_beta,
                              Sigma_beta, A, tol, plot_elbo, verbose = verbose)

    # Orthogonalisation:

    eta_in <- list(
      eta_vec$"p(nu|Sigma_nu)->nu", eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
      eta_vec$"p(zeta)->zeta", eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
    )

    fpc_orthgn(subj_names, L, K, eta_in, time_g, C_g, Psi_g = Psi_g)

  } else {

    p <- length(Y[[1]])

    tmp_names <- lapply(Y, function(Y_i) names(Y_i))
    tmp_names_time <-  lapply(time_obs, function(time_obs_i) names(time_obs_i))
    if (is.null(unlist(unique(tmp_names)))) {
      stopifnot(is.null(unlist(unique(tmp_names_time))))
      var_names <- paste0("variable_", 1:p)
      Y <- lapply(Y, function(Y_i) { names(Y_i) <- var_names; Y_i})
    } else if (!all_same(tmp_names)) {
      stop("Variable names in Y must be the same for all individuals.")
    } else {
      var_names <- names(Y[[1]])
    }

    if (!is.null(unlist(unique(tmp_names_time))) && (!all_same(tmp_names_time) | !isTRUE(all.equal(tmp_names_time[[1]], var_names)))) {
      stop("Variable names for each individual in time_obs must be the same as in Y.")
    }

    time_obs <- lapply(time_obs, function(time_obs_i) { names(time_obs_i) <- var_names; time_obs_i})

    eta_vec <- vmp_gauss_mfpca(n_vmp = maxit, N, p, L, K, C, Y, sigma_zeta, mu_beta,
                               Sigma_beta, A, tol,
                               plot_elbo = plot_elbo, verbose = verbose)

    # Orthogonalisation:

    eta_in <- list(
      eta_vec$"p(nu|Sigma_nu)->nu", eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
      eta_vec$"p(zeta)->zeta", eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
    )

    mfpc_orthgn(subj_names, var_names, Y, L, K, eta_in, time_g, N, p, C_g, Psi_g)

  }

}

vmp_gauss_fpca <- function(n_vmp, N, L, K, C, Y, sigma_zeta, mu_beta,
                           Sigma_beta, A, tol, plot_elbo = FALSE, verbose = TRUE) {

  # Establish necessary parameters:

  Sigma_zeta <- sigma_zeta^2*diag(L)
  T_vec <- sapply(Y, length)
  d <- (K+2)*(L+1)

  # Initialise VMP simulation:

  mu_q_zeta <- vector("list", length=N)
  Sigma_q_zeta <- vector("list", length=N)
  for(i in 1:N) {

    mu_q_zeta[[i]] <- rnorm(L, 0, sigma_zeta)
    Sigma_q_zeta[[i]] <- diag(L)
  }

  eta_vec <- vector("list", length=32)
  names(eta_vec) <- c(
    "nu->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->nu",
    "zeta->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->zeta",
    "sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
    "zeta->p(zeta)", "p(zeta)->zeta",
    "sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
    "a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
    "a_eps->p(a_eps)", "p(a_eps)->a_eps",
    "nu->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->nu",
    "sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
    "sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
    "sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
    "sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
    "a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
    "a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
    "a_m->p(a_m)", "p(a_m)->a_m",
    "a_p->p(a_p)", "p(a_p)->a_p"
  )

  G <- vector("list", length=24)
  names(G) <- c(
    "sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
    "sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
    "a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
    "a_eps->p(a_eps)", "p(a_eps)->a_eps",
    "sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
    "sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
    "sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
    "sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
    "a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
    "a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
    "a_m->p(a_m)", "p(a_m)->a_m",
    "a_p->p(a_p)", "p(a_p)->a_p"
  )

  eta_1_sum <- 0
  eta_2_sum <- 0
  for(i in 1:N) {

    mu_q_zeta_tilde <- c(1, mu_q_zeta[[i]])
    Sigma_q_zeta_tilde <- blkdiag(matrix(0), Sigma_q_zeta[[i]])
    M_q_zeta_zeta_T_tilde <- Sigma_q_zeta_tilde + tcrossprod(mu_q_zeta_tilde)

    sum_val <- cprod(kronecker(t(mu_q_zeta_tilde), C[[i]]), Y[[i]])
    eta_1_sum <- eta_1_sum + sum_val

    sum_val <- as.vector(kronecker(M_q_zeta_zeta_T_tilde, crossprod(C[[i]])))
    eta_2_sum <- eta_2_sum + sum_val
  }
  eta_1 <- 1*eta_1_sum
  eta_2 <- -1/2*eta_2_sum
  eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- c(eta_1, eta_2)

  D_L <- duplication.matrix(L)
  eta_1 <- Reduce(cbind, mu_q_zeta)
  eta_2 <- replicate(N, -0.5*cprod(D_L, as.vector(diag(L) - solve(Sigma_zeta))))
  eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- rbind(eta_1, eta_2)

  eta_1 <- -0.5*sum(T_vec)
  eta_2 <- -0.5*sum(T_vec)
  eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- c(eta_1, eta_2)
  G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- "full"

  eta_vec$"p(zeta)->zeta" <- replicate(
    N, gauss_prior_frag(rep(0, L), Sigma_zeta, use_vech=TRUE)
  )

  eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- c(-3/2, -1/2)
  G$"p(sigsq_eps|a_eps)->sigsq_eps" <- "full"

  eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- c(-1/2, -1/2)
  G$"p(sigsq_eps|a_eps)->a_eps" <- "diag"

  eta_1 <- rep(0, d)
  eta_2 <- -0.5*as.vector(diag(d))
  eta_vec$"p(nu|Sigma_nu)->nu" <- c(eta_1, eta_2)

  eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- c(-K/2, -K/2)
  G$"p(nu|Sigma_nu)->sigsq_m" <- "full"

  eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- replicate(L, c(-K/2, -K/2))
  G$"p(nu|Sigma_nu)->sigsq_p" <- rep("full", L)

  eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- c(-3/2, -1/2)
  G$"p(sigsq_m|a_m)->sigsq_m" <- "full"

  eta_vec$"p(sigsq_p|a_p)->sigsq_p" <- replicate(L, c(-3/2, -1/2))
  G$"p(sigsq_p|a_p)->sigsq_p" <- rep("full", L)

  eta_vec$"p(sigsq_m|a_m)->a_m" <- c(-1/2, -1/2)
  G$"p(sigsq_m|a_m)->a_m" <- "diag"

  eta_vec$"p(sigsq_p|a_p)->a_p" <- replicate(L, c(-1/2, -1/2))
  G$"p(sigsq_p|a_p)->a_p" <- rep("diag", L)

  igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))

  eta_vec$"p(a_eps)->a_eps" <- igw_prior_updates[[2]]
  G$"p(a_eps)->a_eps" <- igw_prior_updates[[1]]

  eta_vec$"p(a_m)->a_m" <- igw_prior_updates[[2]]
  G$"p(a_m)->a_m" <- igw_prior_updates[[1]]

  eta_vec$"p(a_p)->a_p" <- replicate(L, igw_prior_updates[[2]])
  G$"p(a_p)->a_p" <- rep(igw_prior_updates[[1]], L)

  elbo_res <- NULL
  converged <- FALSE
  iter <- 0

  while((!converged) & (iter < n_vmp)) {

    iter <- iter + 1

    if (verbose) cat("Iteration", iter, "\n")

    eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(nu|Sigma_nu)->nu"
    eta_vec$"nu->p(nu|Sigma_nu)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu"

    eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(zeta)->zeta"
    eta_vec$"zeta->p(zeta)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"

    eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
    G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- G$"p(sigsq_eps|a_eps)->sigsq_eps"
    eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
    G$"sigsq_eps->p(sigsq_eps|a_eps)" <- G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"

    eta_vec$"a_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(a_eps)->a_eps"
    G$"a_eps->p(sigsq_eps|a_eps)" <- G$"p(a_eps)->a_eps"
    eta_vec$"a_eps->p(a_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->a_eps"
    G$"a_eps->p(a_eps)" <- G$"p(sigsq_eps|a_eps)->a_eps"

    eta_vec$"sigsq_m->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_m|a_m)->sigsq_m"
    G$"sigsq_m->p(nu|Sigma_nu)" <- G$"p(sigsq_m|a_m)->sigsq_m"
    eta_vec$"sigsq_m->p(sigsq_m|a_m)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_m"
    G$"sigsq_m->p(sigsq_m|a_m)" <- G$"p(nu|Sigma_nu)->sigsq_m"

    eta_vec$"sigsq_p->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_p|a_p)->sigsq_p"
    G$"sigsq_p->p(nu|Sigma_nu)" <- G$"p(sigsq_p|a_p)->sigsq_p"
    eta_vec$"sigsq_p->p(sigsq_p|a_p)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_p"
    G$"sigsq_p->p(sigsq_p|a_p)" <- G$"p(nu|Sigma_nu)->sigsq_p"

    eta_vec$"a_m->p(sigsq_m|a_m)" <- eta_vec$"p(a_m)->a_m"
    G$"a_m->p(sigsq_m|a_m)" <- G$"p(a_m)->a_m"
    eta_vec$"a_m->p(a_m)" <- eta_vec$"p(sigsq_m|a_m)->a_m"
    G$"a_m->p(a_m)" <- G$"p(sigsq_m|a_m)->a_m"

    eta_vec$"a_p->p(sigsq_p|a_p)" <- eta_vec$"p(a_p)->a_p"
    G$"a_p->p(sigsq_p|a_p)" <- G$"p(a_p)->a_p"
    eta_vec$"a_p->p(a_p)" <- eta_vec$"p(sigsq_p|a_p)->a_p"
    G$"a_p->p(a_p)" <- G$"p(sigsq_p|a_p)->a_p"

    # Update p(Y|nu,zeta,sigsq_eps) fragment:

    eta_in <- list(
      eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
      eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
      eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
      eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
      eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
      eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
    )

    G_in <- list(
      G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
      G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
    )

    fpc_lik_fragment <- fpc_lik_frag(
      eta_in, G_in, C, Y, T_vec, L
    )

    eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- fpc_lik_fragment$"eta"[[1]]
    eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- fpc_lik_fragment$"eta"[[2]]
    eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- fpc_lik_fragment$"eta"[[3]]

    G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- fpc_lik_fragment$"G"[[1]]

    # Update p(nu|Sigma_nu) fragment:

    eta_in <- list(
      eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
      eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
      eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
    )

    G_in <- list(
      G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
      G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
    )

    fpc_gauss_pen_fragment <- fpc_gauss_pen_frag(
      eta_in, G_in, L, mu_beta, Sigma_beta
    )

    eta_vec$"p(nu|Sigma_nu)->nu" <- fpc_gauss_pen_fragment$"eta"[[1]]
    eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"eta"[[2]]
    eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"eta"[[3]]

    G$"p(nu|Sigma_nu)->sigsq_m" <- fpc_gauss_pen_fragment$"G"[[1]]
    G$"p(nu|Sigma_nu)->sigsq_p" <- fpc_gauss_pen_fragment$"G"[[2]]

    # Update p(sigsq_eps|a_eps) fragment:

    eta_in <- list(
      eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
      eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
      eta_vec$"a_eps->p(sigsq_eps|a_eps)",
      eta_vec$"p(sigsq_eps|a_eps)->a_eps"
    )

    iter_igw_fragment <- iter_igw_frag(
      eta_in, G$"a_eps->p(sigsq_eps|a_eps)",
      1, G$"sigsq_eps->p(sigsq_eps|a_eps)"
    )

    eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"eta"[[1]]
    eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"eta"[[2]]

    G$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"G"[[1]]
    G$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"G"[[2]]

    # Update p(sigsq_m|a_m) fragment:

    eta_in <- list(
      eta_vec$"sigsq_m->p(sigsq_m|a_m)",
      eta_vec$"p(sigsq_m|a_m)->sigsq_m",
      eta_vec$"a_m->p(sigsq_m|a_m)",
      eta_vec$"p(sigsq_m|a_m)->a_m"
    )

    iter_igw_fragment <- iter_igw_frag(
      eta_in, G$"a_m->p(sigsq_m|a_m)",
      1, G$"sigsq_m->p(sigsq_m|a_m)"
    )

    eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"eta"[[1]]
    eta_vec$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"eta"[[2]]

    G$"p(sigsq_m|a_m)->sigsq_m" <- iter_igw_fragment$"G"[[1]]
    G$"p(sigsq_m|a_m)->a_m" <- iter_igw_fragment$"G"[[2]]

    # Update p(sigsq_p|a_p) fragment:

    for(l in 1:L) {

      eta_in <- list(
        eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
        eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
        eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
        eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
      )

      iter_igw_fragment <- iter_igw_frag(
        eta_in, G$"a_p->p(sigsq_p|a_p)"[l],
        1, G$"sigsq_p->p(sigsq_p|a_p)"[l]
      )

      eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l] <- iter_igw_fragment$"eta"[[1]]
      eta_vec$"p(sigsq_p|a_p)->a_p"[,l] <- iter_igw_fragment$"eta"[[2]]

      G$"p(sigsq_p|a_p)->sigsq_p"[l] <- iter_igw_fragment$"G"[[1]]
      G$"p(sigsq_p|a_p)->a_p"[l] <- iter_igw_fragment$"G"[[2]]
    }

    # Compute the entropy:

    ent <- 0

    eta_nu <- list(
      eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
      eta_vec$"p(nu|Sigma_nu)->nu"
    )
    ent_nu <- entropy_gauss(eta_nu, use_vech=FALSE)

    ent <- ent + ent_nu

    ent_zeta <- 0
    for(i in 1:N) {

      eta_zeta <- list(
        eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"[,i],
        eta_vec$"p(zeta)->zeta"[,i]
      )
      sum_val <- entropy_gauss(eta_zeta, use_vech=TRUE)
      ent_zeta <- ent_zeta + sum_val
    }

    ent <- ent + ent_zeta

    eta_sigsq_eps <- list(
      eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
      eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
    )
    G_sigsq_eps <- c(
      G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
      G$"p(sigsq_eps|a_eps)->sigsq_eps"
    )
    ent_sigsq_eps <- entropy_igw(eta_sigsq_eps, G_sigsq_eps)

    ent <- ent + ent_sigsq_eps

    eta_a_eps <- list(
      eta_vec$"p(sigsq_eps|a_eps)->a_eps",
      eta_vec$"p(a_eps)->a_eps"
    )
    G_a_eps <- c(
      G$"p(sigsq_eps|a_eps)->a_eps",
      G$"p(a_eps)->a_eps"
    )
    ent_a_eps <- entropy_igw(eta_a_eps, G_a_eps)

    ent <- ent + ent_a_eps

    eta_sigsq_m <- list(
      eta_vec$"p(nu|Sigma_nu)->sigsq_m",
      eta_vec$"p(sigsq_m|a_m)->sigsq_m"
    )
    G_sigsq_m <- c(
      G$"p(nu|Sigma_nu)->sigsq_m",
      G$"p(sigsq_m|a_m)->sigsq_m"
    )
    ent_sigsq_m <- entropy_igw(eta_sigsq_m, G_sigsq_m)

    ent <- ent + ent_sigsq_m

    eta_a_m <- list(
      eta_vec$"p(sigsq_m|a_m)->a_m",
      eta_vec$"p(a_m)->a_m"
    )
    G_a_m <- c(
      G$"p(sigsq_m|a_m)->a_m",
      G$"p(a_m)->a_m"
    )
    ent_a_m <- entropy_igw(eta_a_m, G_a_m)

    ent <- ent + ent_a_m

    ent_sigsq_p <- 0
    ent_a_p <- 0
    for(l in 1:L) {

      eta_sigsq_p <- list(
        eta_vec$"p(nu|Sigma_nu)->sigsq_p"[,l],
        eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l]
      )
      G_sigsq_p <- list(
        G$"p(nu|Sigma_nu)->sigsq_p"[l],
        G$"p(sigsq_p|a_p)->sigsq_p"[l]
      )
      sum_val <- entropy_igw(eta_sigsq_p, G_sigsq_p)
      ent_sigsq_p <- ent_sigsq_p + sum_val

      eta_a_p <- list(
        eta_vec$"p(sigsq_p|a_p)->a_p"[,l],
        eta_vec$"p(a_p)->a_p"[,l]
      )
      G_a_p <- c(
        G$"p(sigsq_p|a_p)->a_p"[l],
        G$"p(a_p)->a_p"[l]
      )
      sum_val <- entropy_igw(eta_a_p, G_a_p)
      ent_a_p <- ent_a_p + sum_val
    }

    ent <- ent + ent_sigsq_p + ent_a_p

    # Compute the cross-entropy:

    c_ent <- 0

    eta_in <- list(
      eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
      eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
      eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
      eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
      eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
      eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
    )
    G_in <- list(
      G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
      G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
    )
    c_ent_p_Y <- cross_entropy_fpc_lik_frag(eta_in, G_in, C, Y, T_vec, L)

    c_ent <- c_ent + c_ent_p_Y

    eta_in <- list(
      eta_vec$"nu->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->nu",
      eta_vec$"sigsq_m->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_m",
      eta_vec$"sigsq_p->p(nu|Sigma_nu)", eta_vec$"p(nu|Sigma_nu)->sigsq_p"
    )
    G_in <- list(
      G$"sigsq_m->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_m",
      G$"sigsq_p->p(nu|Sigma_nu)", G$"p(nu|Sigma_nu)->sigsq_p"
    )
    c_ent_p_nu <- cross_entropy_fpc_gauss_pen(eta_in, G_in, L, mu_beta, Sigma_beta)

    c_ent <- c_ent + c_ent_p_nu

    c_ent_p_zeta <- 0
    for(i in 1:N) {

      eta_in <- list(
        eta_vec$"zeta->p(zeta)"[,i],
        eta_vec$"p(zeta)->zeta"[,i]
      )
      sum_val <- cross_entropy_gauss_prior(
        eta_in, rep(0, L),
        Sigma_zeta, use_vech=TRUE
      )
      c_ent_p_zeta <- c_ent_p_zeta + sum_val
    }

    c_ent <- c_ent + c_ent_p_zeta

    eta_in <- list(
      eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
      eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
      eta_vec$"a_eps->p(sigsq_eps|a_eps)",
      eta_vec$"p(sigsq_eps|a_eps)->a_eps"
    )
    G_mess <- G$"sigsq_eps->p(sigsq_eps|a_eps)"
    G_hyper <- G$"a_eps->p(sigsq_eps|a_eps)"
    c_ent_p_sigsq_eps <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)

    c_ent <- c_ent + c_ent_p_sigsq_eps

    eta_in <- list(
      eta_vec$"a_eps->p(a_eps)",
      eta_vec$"p(a_eps)->a_eps"
    )
    G_in <- c(
      G$"a_eps->p(a_eps)",
      G$"p(a_eps)->a_eps"
    )
    c_ent_p_a_eps <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)

    c_ent <- c_ent + c_ent_p_a_eps

    eta_in <- list(
      eta_vec$"sigsq_m->p(sigsq_m|a_m)",
      eta_vec$"p(sigsq_m|a_m)->sigsq_m",
      eta_vec$"a_m->p(sigsq_m|a_m)",
      eta_vec$"p(sigsq_m|a_m)->a_m"
    )
    G_mess <- G$"sigsq_m->p(sigsq_m|a_m)"
    G_hyper <- G$"a_m->p(sigsq_m|a_m)"
    c_ent_p_sigsq_m <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)

    c_ent <- c_ent + c_ent_p_sigsq_m

    eta_in <- list(
      eta_vec$"a_m->p(a_m)",
      eta_vec$"p(a_m)->a_m"
    )
    G_in <- c(
      G$"a_m->p(a_m)",
      G$"p(a_m)->a_m"
    )
    c_ent_p_a_m <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)

    c_ent <- c_ent + c_ent_p_a_m

    c_ent_p_sigsq_p <- 0
    c_ent_p_a_p <- 0
    for(l in 1:L) {

      eta_in <- list(
        eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
        eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
        eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
        eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
      )
      G_mess <- G$"sigsq_p->p(sigsq_p|a_p)"[l]
      G_hyper <- G$"a_p->p(sigsq_p|a_p)"[l]
      sum_val <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)
      c_ent_p_sigsq_p <- c_ent_p_sigsq_p + sum_val

      eta_in <- list(
        eta_vec$"a_p->p(a_p)"[,l],
        eta_vec$"p(a_p)->a_p"[,l]
      )
      G_in <- c(
        G$"a_p->p(a_p)"[l],
        G$"p(a_p)->a_p"[l]
      )
      sum_val <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)
      c_ent_p_a_p <- c_ent_p_a_p + sum_val
    }

    c_ent <- c_ent + c_ent_p_sigsq_p + c_ent_p_a_p

    # Compute the ELBO

    elbo_new <- ent - c_ent
    elbo_res <- c(elbo_res, elbo_new)

    if(plot_elbo) {

      plot(1:iter, elbo_res, pch=16, cex=0.4, xlab="iterations", ylab="ELBO")
    }

    if(iter > 1) {

      elbo_old <- elbo_res[iter - 1]
      tol_1_satisfied <- (abs(elbo_new/elbo_old - 1) < tol)

      if(iter > 2) {

        elbo_old <- elbo_res[iter - 2]
        tol_2_satisfied <- (abs(elbo_new/elbo_old - 1) < tol)
      } else {

        tol_2_satisfied <- FALSE
      }

      tol_satisfied <- (tol_1_satisfied || tol_2_satisfied)

      if(tol_satisfied) {

        converged <- TRUE
      }
    }
  }

  # Get the list of natural parameter vectors:

  return(eta_vec)
}


vmp_gauss_mfpca <- function(n_vmp, N, p, L, K, C, Y, sigma_zeta, mu_beta, Sigma_beta,
                            A, tol, plot_elbo = FALSE, verbose = TRUE) {

  # Establish necessary parameters:

  n <- Reduce(rbind, lapply(Y, function(x) sapply(x, length)))
  d <- (K+2)*(L+1) # this is a vector since K is a vector

  # Initialise VMP simulation:

  E_q_zeta <- vector("list", length = N)
  Cov_q_zeta <- vector("list", length = N)
  for(i in 1:N) {

    E_q_zeta[[i]] <- rnorm(L, 0, sigma_zeta)
    Cov_q_zeta[[i]] <- diag(L)
  }

  eta_vec <- vector("list", length=32)
  names(eta_vec) <- c(
    "nu->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->nu",
    "zeta->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->zeta",
    "sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
    "zeta->p(zeta)", "p(zeta)->zeta",
    "sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
    "a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
    "a_eps->p(a_eps)", "p(a_eps)->a_eps",
    "nu->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->nu",
    "sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
    "sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
    "sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
    "sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
    "a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
    "a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
    "a_m->p(a_m)", "p(a_m)->a_m",
    "a_p->p(a_p)", "p(a_p)->a_p"
  )

  G <- vector("list", length=24)
  names(G) <- c(
    "sigsq_eps->p(Y|nu,zeta,sigsq_eps)", "p(Y|nu,zeta,sigsq_eps)->sigsq_eps",
    "sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
    "a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
    "a_eps->p(a_eps)", "p(a_eps)->a_eps",
    "sigsq_m->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_m",
    "sigsq_p->p(nu|Sigma_nu)", "p(nu|Sigma_nu)->sigsq_p",
    "sigsq_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->sigsq_m",
    "sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
    "a_m->p(sigsq_m|a_m)", "p(sigsq_m|a_m)->a_m",
    "a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
    "a_m->p(a_m)", "p(a_m)->a_m",
    "a_p->p(a_p)", "p(a_p)->a_p"
  )

  eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- vector("list", length = p)
  eta_vec$"p(nu|Sigma_nu)->nu" <- vector("list", length = p)
  for(j in 1:p) {

    eta_1_sum <- 0
    eta_2_sum <- 0
    for(i in 1:N) {

      E_q_zeta_tilde <- c(1, E_q_zeta[[i]])
      Cov_q_zeta_tilde <- adiag(0, Cov_q_zeta[[i]])
      E_q_tcross_zeta_tilde <- Cov_q_zeta_tilde + tcrossprod(E_q_zeta_tilde)

      sum_val <- cprod(kronecker(t(E_q_zeta_tilde), C[[i]][[j]]), Y[[i]][[j]])
      eta_1_sum <- eta_1_sum + sum_val

      sum_val <- as.vector(kronecker(E_q_tcross_zeta_tilde, crossprod(C[[i]][[j]])))
      eta_2_sum <- eta_2_sum + sum_val
    }

    eta_1 <- eta_1_sum
    eta_2 <- -0.5*eta_2_sum
    eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu"[[j]] <- c(eta_1, eta_2)

    eta_1 <- rep(0, d[j])
    eta_2 <- -0.5*as.vector(diag(d[j]))
    eta_vec$"p(nu|Sigma_nu)->nu"[[j]] <- c(eta_1, eta_2)
  }

  D_L <- duplication.matrix(L)
  eta_1 <- Reduce(cbind, E_q_zeta)
  eta_2 <- replicate(N, -0.5*cprod(D_L, as.vector(diag(L))))
  eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- rbind(eta_1, eta_2)

  eta_vec$"p(zeta)->zeta" <- replicate(
    N, gauss_prior_frag(rep(0, L), diag(L), use_vech = TRUE)
  )

  eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- vector("list", length = p)
  G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- vector("list", length = p)
  eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- vector("list", length = p)
  G$"p(sigsq_eps|a_eps)->sigsq_eps" <- vector("list", length = p)
  eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- vector("list", length = p)
  G$"p(sigsq_eps|a_eps)->a_eps" <- vector("list", length = p)
  for(j in 1:p) {

    eta_1 <- -0.5*sum(n[, j])
    eta_2 <- -0.5*sum(n[, j])
    eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"[[j]] <- c(eta_1, eta_2)
    G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"[[j]] <- "full"

    eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"[[j]] <- c(-3/2, -1/2)
    G$"p(sigsq_eps|a_eps)->sigsq_eps"[[j]] <- "full"

    eta_vec$"p(sigsq_eps|a_eps)->a_eps"[[j]] <- c(-1/2, -1/2)
    G$"p(sigsq_eps|a_eps)->a_eps"[[j]] <- "diag"
  }

  eta_vec$"p(nu|Sigma_nu)->sigsq_m" <- vector("list", length = p)
  G$"p(nu|Sigma_nu)->sigsq_m" <- vector("list", length = p)
  eta_vec$"p(nu|Sigma_nu)->sigsq_p" <- vector("list", length = p)
  G$"p(nu|Sigma_nu)->sigsq_p" <- vector("list", length = p)
  eta_vec$"p(sigsq_m|a_m)->sigsq_m" <- vector("list", length = p)
  G$"p(sigsq_m|a_m)->sigsq_m" <- vector("list", length = p)
  eta_vec$"p(sigsq_p|a_p)->sigsq_p" <- vector("list", length = p)
  G$"p(sigsq_p|a_p)->sigsq_p" <- vector("list", length = p)
  eta_vec$"p(sigsq_m|a_m)->a_m" <- vector("list", length = p)
  G$"p(sigsq_m|a_m)->a_m" <- vector("list", length = p)
  eta_vec$"p(sigsq_p|a_p)->a_p" <- vector("list", length = p)
  G$"p(sigsq_p|a_p)->a_p" <- vector("list", length = p)
  for(j in 1:p) {

    eta_vec$"p(nu|Sigma_nu)->sigsq_m"[[j]] <- c(-K[j]/2, -K[j]/2)
    G$"p(nu|Sigma_nu)->sigsq_m"[[j]] <- "full"

    eta_vec$"p(nu|Sigma_nu)->sigsq_p"[[j]] <- replicate(L, c(-K[j]/2, -K[j]/2))
    G$"p(nu|Sigma_nu)->sigsq_p"[[j]] <- rep("full", L)

    eta_vec$"p(sigsq_m|a_m)->sigsq_m"[[j]] <- c(-3/2, -1/2)
    G$"p(sigsq_m|a_m)->sigsq_m"[[j]] <- "full"

    eta_vec$"p(sigsq_p|a_p)->sigsq_p"[[j]] <- replicate(L, c(-3/2, -1/2))
    G$"p(sigsq_p|a_p)->sigsq_p"[[j]] <- rep("full", L)

    eta_vec$"p(sigsq_m|a_m)->a_m"[[j]] <- c(-1/2, -1/2)
    G$"p(sigsq_m|a_m)->a_m"[[j]] <- "diag"

    eta_vec$"p(sigsq_p|a_p)->a_p"[[j]] <- replicate(L, c(-1/2, -1/2))
    G$"p(sigsq_p|a_p)->a_p"[[j]] <- rep("diag", L)
  }

  igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))
  eta_vec$"p(a_eps)->a_eps" <- vector("list", length = p)
  G$"p(a_eps)->a_eps" <- vector("list", length = p)
  eta_vec$"p(a_m)->a_m" <- vector("list", length = p)
  G$"p(a_m)->a_m" <- vector("list", length = p)
  eta_vec$"p(a_p)->a_p" <- vector("list", length = p)
  G$"p(a_p)->a_p" <- vector("list", length = p)
  for(j in 1:p) {

    eta_vec$"p(a_eps)->a_eps"[[j]] <- igw_prior_updates[[2]]
    G$"p(a_eps)->a_eps"[[j]] <- igw_prior_updates[[1]]

    eta_vec$"p(a_m)->a_m"[[j]] <- igw_prior_updates[[2]]
    G$"p(a_m)->a_m"[[j]] <- igw_prior_updates[[1]]

    eta_vec$"p(a_p)->a_p"[[j]] <- replicate(L, igw_prior_updates[[2]])
    G$"p(a_p)->a_p"[[j]] <- rep(igw_prior_updates[[1]], L)
  }

  elbo_res <- NULL
  converged <- FALSE
  iter <- 0
  while((!converged) & (iter < n_vmp)) {

    iter <- iter + 1

    if (verbose) cat("Iteration", iter, "\n")

    eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(nu|Sigma_nu)->nu"
    eta_vec$"nu->p(nu|Sigma_nu)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu"

    eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(zeta)->zeta"
    eta_vec$"zeta->p(zeta)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"

    eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
    G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)" <- G$"p(sigsq_eps|a_eps)->sigsq_eps"
    eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
    G$"sigsq_eps->p(sigsq_eps|a_eps)" <- G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"

    eta_vec$"a_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(a_eps)->a_eps"
    G$"a_eps->p(sigsq_eps|a_eps)" <- G$"p(a_eps)->a_eps"
    eta_vec$"a_eps->p(a_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->a_eps"
    G$"a_eps->p(a_eps)" <- G$"p(sigsq_eps|a_eps)->a_eps"

    eta_vec$"sigsq_m->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_m|a_m)->sigsq_m"
    G$"sigsq_m->p(nu|Sigma_nu)" <- G$"p(sigsq_m|a_m)->sigsq_m"
    eta_vec$"sigsq_m->p(sigsq_m|a_m)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_m"
    G$"sigsq_m->p(sigsq_m|a_m)" <- G$"p(nu|Sigma_nu)->sigsq_m"

    eta_vec$"sigsq_p->p(nu|Sigma_nu)" <- eta_vec$"p(sigsq_p|a_p)->sigsq_p"
    G$"sigsq_p->p(nu|Sigma_nu)" <- G$"p(sigsq_p|a_p)->sigsq_p"
    eta_vec$"sigsq_p->p(sigsq_p|a_p)" <- eta_vec$"p(nu|Sigma_nu)->sigsq_p"
    G$"sigsq_p->p(sigsq_p|a_p)" <- G$"p(nu|Sigma_nu)->sigsq_p"

    eta_vec$"a_m->p(sigsq_m|a_m)" <- eta_vec$"p(a_m)->a_m"
    G$"a_m->p(sigsq_m|a_m)" <- G$"p(a_m)->a_m"
    eta_vec$"a_m->p(a_m)" <- eta_vec$"p(sigsq_m|a_m)->a_m"
    G$"a_m->p(a_m)" <- G$"p(sigsq_m|a_m)->a_m"

    eta_vec$"a_p->p(sigsq_p|a_p)" <- eta_vec$"p(a_p)->a_p"
    G$"a_p->p(sigsq_p|a_p)" <- G$"p(a_p)->a_p"
    eta_vec$"a_p->p(a_p)" <- eta_vec$"p(sigsq_p|a_p)->a_p"
    G$"a_p->p(a_p)" <- G$"p(sigsq_p|a_p)->a_p"

    # Update p(Y|nu,zeta,sigsq_eps) fragment:

    eta_in <- list(
      eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
      eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
      eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
      eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
      eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
      eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
    )

    G_in <- list(
      G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
      G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
    )

    mfpc_lik_fragment <- mfpc_lik_frag(eta_in, G_in, C, Y, L)

    eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu" <- mfpc_lik_fragment$"eta"[[1]]
    eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta" <- mfpc_lik_fragment$"eta"[[2]]
    eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- mfpc_lik_fragment$"eta"[[3]]

    G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps" <- mfpc_lik_fragment$"G"[[1]]

    # For j = 1, ..., p, update p(nu[j]|Sigma_nu[j]) fragment:

    for(j in 1:p) {

      eta_in <- list(
        eta_vec$"nu->p(nu|Sigma_nu)"[[j]], eta_vec$"p(nu|Sigma_nu)->nu"[[j]],
        eta_vec$"sigsq_m->p(nu|Sigma_nu)"[[j]], eta_vec$"p(nu|Sigma_nu)->sigsq_m"[[j]],
        eta_vec$"sigsq_p->p(nu|Sigma_nu)"[[j]], eta_vec$"p(nu|Sigma_nu)->sigsq_p"[[j]]
      )

      G_in <- list(
        G$"sigsq_m->p(nu|Sigma_nu)"[[j]], G$"p(nu|Sigma_nu)->sigsq_m"[[j]],
        G$"sigsq_p->p(nu|Sigma_nu)"[[j]], G$"p(nu|Sigma_nu)->sigsq_p"[[j]]
      )

      fpc_gauss_pen_fragment <- fpc_gauss_pen_frag(
        eta_in, G_in, L, mu_beta, Sigma_beta
      )

      eta_vec$"p(nu|Sigma_nu)->nu"[[j]] <- fpc_gauss_pen_fragment$"eta"[[1]]
      eta_vec$"p(nu|Sigma_nu)->sigsq_m"[[j]] <- fpc_gauss_pen_fragment$"eta"[[2]]
      eta_vec$"p(nu|Sigma_nu)->sigsq_p"[[j]] <- fpc_gauss_pen_fragment$"eta"[[3]]

      G$"p(nu|Sigma_nu)->sigsq_m"[[j]] <- fpc_gauss_pen_fragment$"G"[[1]]
      G$"p(nu|Sigma_nu)->sigsq_p"[[j]] <- fpc_gauss_pen_fragment$"G"[[2]]
    }

    # For j = 1, ..., p, update p(sigsq_m[j]|a_m[j]) fragment:

    for(j in 1:p) {

      eta_in <- list(
        eta_vec$"sigsq_m->p(sigsq_m|a_m)"[[j]],
        eta_vec$"p(sigsq_m|a_m)->sigsq_m"[[j]],
        eta_vec$"a_m->p(sigsq_m|a_m)"[[j]],
        eta_vec$"p(sigsq_m|a_m)->a_m"[[j]]
      )

      iter_igw_fragment <- iter_igw_frag(
        eta_in, G$"a_m->p(sigsq_m|a_m)"[[j]],
        1, G$"sigsq_m->p(sigsq_m|a_m)"[[j]]
      )

      eta_vec$"p(sigsq_m|a_m)->sigsq_m"[[j]] <- iter_igw_fragment$"eta"[[1]]
      eta_vec$"p(sigsq_m|a_m)->a_m"[[j]] <- iter_igw_fragment$"eta"[[2]]

      G$"p(sigsq_m|a_m)->sigsq_m"[[j]] <- iter_igw_fragment$"G"[[1]]
      G$"p(sigsq_m|a_m)->a_m"[[j]] <- iter_igw_fragment$"G"[[2]]
    }

    # For j = 1, ..., p, for l = 1, ..., L, update p(sigsq_p[j][l]|a_p[j][l]) fragment:

    for(j in 1:p) {

      for(l in 1:L) {

        eta_in <- list(
          eta_vec$"sigsq_p->p(sigsq_p|a_p)"[[j]][, l],
          eta_vec$"p(sigsq_p|a_p)->sigsq_p"[[j]][, l],
          eta_vec$"a_p->p(sigsq_p|a_p)"[[j]][, l],
          eta_vec$"p(sigsq_p|a_p)->a_p"[[j]][, l]
        )

        iter_igw_fragment <- iter_igw_frag(
          eta_in, G$"a_p->p(sigsq_p|a_p)"[[j]][l],
          1, G$"sigsq_p->p(sigsq_p|a_p)"[[j]][l]
        )

        eta_vec$"p(sigsq_p|a_p)->sigsq_p"[[j]][, l] <- iter_igw_fragment$"eta"[[1]]
        eta_vec$"p(sigsq_p|a_p)->a_p"[[j]][, l] <- iter_igw_fragment$"eta"[[2]]

        G$"p(sigsq_p|a_p)->sigsq_p"[[j]][l] <- iter_igw_fragment$"G"[[1]]
        G$"p(sigsq_p|a_p)->a_p"[[j]][l] <- iter_igw_fragment$"G"[[2]]
      }
    }

    # Compute the entropy:

    ent <- 0

    for(j in 1:p) {

      eta_in <- list(
        eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu"[[j]],
        eta_vec$"p(nu|Sigma_nu)->nu"[[j]]
      )
      ent_nu <- entropy_gauss(eta_in, use_vech = FALSE)

      ent <- ent + ent_nu
    }

    for(i in 1:N) {

      eta_in <- list(
        eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"[, i],
        eta_vec$"p(zeta)->zeta"[, i]
      )
      ent_zeta <- entropy_gauss(eta_in, use_vech = TRUE)

      ent <- ent + ent_zeta
    }

    for(j in 1:p) {

      eta_sigsq_eps <- list(
        eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"[[j]],
        eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"[[j]]
      )
      G_sigsq_eps <- c(
        G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"[j],
        G$"p(sigsq_eps|a_eps)->sigsq_eps"[j]
      )
      ent_sigsq_eps <- entropy_igw(eta_sigsq_eps, G_sigsq_eps)

      ent <- ent + ent_sigsq_eps

      eta_a_eps <- list(
        eta_vec$"p(sigsq_eps|a_eps)->a_eps"[[j]],
        eta_vec$"p(a_eps)->a_eps"[[j]]
      )
      G_a_eps <- c(
        G$"p(sigsq_eps|a_eps)->a_eps"[[j]],
        G$"p(a_eps)->a_eps"[[j]]
      )
      ent_a_eps <- entropy_igw(eta_a_eps, G_a_eps)

      ent <- ent + ent_a_eps

      eta_sigsq_m <- list(
        eta_vec$"p(nu|Sigma_nu)->sigsq_m"[[j]],
        eta_vec$"p(sigsq_m|a_m)->sigsq_m"[[j]]
      )
      G_sigsq_m <- c(
        G$"p(nu|Sigma_nu)->sigsq_m"[[j]],
        G$"p(sigsq_m|a_m)->sigsq_m"[[j]]
      )
      ent_sigsq_m <- entropy_igw(eta_sigsq_m, G_sigsq_m)

      ent <- ent + ent_sigsq_m

      eta_a_m <- list(
        eta_vec$"p(sigsq_m|a_m)->a_m"[[j]],
        eta_vec$"p(a_m)->a_m"[[j]]
      )
      G_a_m <- c(
        G$"p(sigsq_m|a_m)->a_m"[[j]],
        G$"p(a_m)->a_m"[[j]]
      )
      ent_a_m <- entropy_igw(eta_a_m, G_a_m)

      ent <- ent + ent_a_m

      eta_a_m <- list(
        eta_vec$"p(sigsq_m|a_m)->a_m"[[j]],
        eta_vec$"p(a_m)->a_m"[[j]]
      )
      G_a_m <- c(
        G$"p(sigsq_m|a_m)->a_m"[[j]],
        G$"p(a_m)->a_m"[[j]]
      )
      ent_a_m <- entropy_igw(eta_a_m, G_a_m)

      ent <- ent + ent_a_m

      for(l in 1:L) {

        eta_sigsq_p <- list(
          eta_vec$"p(nu|Sigma_nu)->sigsq_p"[[j]][, l],
          eta_vec$"p(sigsq_p|a_p)->sigsq_p"[[j]][, l]
        )
        G_sigsq_p <- list(
          G$"p(nu|Sigma_nu)->sigsq_p"[[j]][l],
          G$"p(sigsq_p|a_p)->sigsq_p"[[j]][l]
        )
        ent_sigsq_p <- entropy_igw(eta_sigsq_p, G_sigsq_p)

        ent <- ent + ent_sigsq_p

        eta_a_p <- list(
          eta_vec$"p(sigsq_p|a_p)->a_p"[[j]][, l],
          eta_vec$"p(a_p)->a_p"[[j]][,l]
        )
        G_a_p <- c(
          G$"p(sigsq_p|a_p)->a_p"[[j]][l],
          G$"p(a_p)->a_p"[[j]][l]
        )
        ent_a_p <- entropy_igw(eta_a_p, G_a_p)
        ent <- ent + ent_a_p
      }
    }

    # Compute the cross-entropy:

    c_ent <- 0

    eta_in <- list(
      eta_vec$"nu->p(Y|nu,zeta,sigsq_eps)",
      eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
      eta_vec$"zeta->p(Y|nu,zeta,sigsq_eps)",
      eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta",
      eta_vec$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
      eta_vec$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
    )
    G_in <- list(
      G$"sigsq_eps->p(Y|nu,zeta,sigsq_eps)",
      G$"p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
    )

    c_ent_p_Y <- cross_entropy_mfpc_lik_frag(eta_in, G_in, n, C, Y, L)

    c_ent <- c_ent + c_ent_p_Y

    for(j in 1:p) {

      eta_in <- list(
        eta_vec$"nu->p(nu|Sigma_nu)"[[j]], eta_vec$"p(nu|Sigma_nu)->nu"[[j]],
        eta_vec$"sigsq_m->p(nu|Sigma_nu)"[[j]], eta_vec$"p(nu|Sigma_nu)->sigsq_m"[[j]],
        eta_vec$"sigsq_p->p(nu|Sigma_nu)"[[j]], eta_vec$"p(nu|Sigma_nu)->sigsq_p"[[j]]
      )
      G_in <- list(
        G$"sigsq_m->p(nu|Sigma_nu)"[[j]], G$"p(nu|Sigma_nu)->sigsq_m"[[j]],
        G$"sigsq_p->p(nu|Sigma_nu)"[[j]], G$"p(nu|Sigma_nu)->sigsq_p"[[j]]
      )
      c_ent_p_nu <- cross_entropy_fpc_gauss_pen(eta_in, G_in, L, mu_beta, Sigma_beta)

      c_ent <- c_ent + c_ent_p_nu
    }

    for(i in 1:N) {

      eta_in <- list(
        eta_vec$"zeta->p(zeta)"[, i],
        eta_vec$"p(zeta)->zeta"[, i]
      )
      c_ent_p_zeta <- cross_entropy_gauss_prior(eta_in, rep(0, L), diag(L), use_vech = TRUE)

      c_ent <- c_ent + c_ent_p_zeta
    }

    for(j in 1:p) {

      eta_in <- list(
        eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)"[[j]],
        eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"[[j]],
        eta_vec$"a_eps->p(sigsq_eps|a_eps)"[[j]],
        eta_vec$"p(sigsq_eps|a_eps)->a_eps"[[j]]
      )
      G_mess <- G$"sigsq_eps->p(sigsq_eps|a_eps)"[[j]]
      G_hyper <- G$"a_eps->p(sigsq_eps|a_eps)"[[j]]
      c_ent_p_sigsq_eps <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)

      c_ent <- c_ent + c_ent_p_sigsq_eps

      eta_in <- list(
        eta_vec$"a_eps->p(a_eps)"[[j]],
        eta_vec$"p(a_eps)->a_eps"[[j]]
      )
      G_in <- c(
        G$"a_eps->p(a_eps)"[[j]],
        G$"p(a_eps)->a_eps"[[j]]
      )
      c_ent_p_a_eps <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)

      c_ent <- c_ent + c_ent_p_a_eps

      eta_in <- list(
        eta_vec$"sigsq_m->p(sigsq_m|a_m)"[[j]],
        eta_vec$"p(sigsq_m|a_m)->sigsq_m"[[j]],
        eta_vec$"a_m->p(sigsq_m|a_m)"[[j]],
        eta_vec$"p(sigsq_m|a_m)->a_m"[[j]]
      )
      G_mess <- G$"sigsq_m->p(sigsq_m|a_m)"[[j]]
      G_hyper <- G$"a_m->p(sigsq_m|a_m)"[[j]]
      c_ent_p_sigsq_m <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)

      c_ent <- c_ent + c_ent_p_sigsq_m

      eta_in <- list(
        eta_vec$"a_m->p(a_m)"[[j]],
        eta_vec$"p(a_m)->a_m"[[j]]
      )
      G_in <- c(
        G$"a_m->p(a_m)"[[j]],
        G$"p(a_m)->a_m"[[j]]
      )
      c_ent_p_a_m <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)

      c_ent <- c_ent + c_ent_p_a_m

      for(l in 1:L) {

        eta_in <- list(
          eta_vec$"sigsq_p->p(sigsq_p|a_p)"[[j]][, l],
          eta_vec$"p(sigsq_p|a_p)->sigsq_p"[[j]][, l],
          eta_vec$"a_p->p(sigsq_p|a_p)"[[j]][, l],
          eta_vec$"p(sigsq_p|a_p)->a_p"[[j]][, l]
        )
        G_mess <- G$"sigsq_p->p(sigsq_p|a_p)"[[j]][l]
        G_hyper <- G$"a_p->p(sigsq_p|a_p)"[[j]][l]
        c_ent_p_sigsq_p <- cross_entropy_iter_igw(eta_in, G_mess, 1, G_hyper)

        c_ent <- c_ent + c_ent_p_sigsq_p

        eta_in <- list(
          eta_vec$"a_p->p(a_p)"[[j]][, l],
          eta_vec$"p(a_p)->a_p"[[j]][, l]
        )
        G_in <- c(
          G$"a_p->p(a_p)"[[j]][l],
          G$"p(a_p)->a_p"[[j]][l]
        )
        c_ent_p_a_p <- cross_entropy_igw_prior(eta_in, G_in, 1, 1/A^2)

        c_ent <- c_ent + c_ent_p_a_p
      }
    }

    # Compute the ELBO

    elbo_new <- ent - c_ent
    elbo_res <- c(elbo_res, elbo_new)

    if(plot_elbo) {

      plot(1:iter, elbo_res, pch=16, cex=0.4, xlab = "iterations", ylab = "ELBO")
    }

    if(iter > 1) {

      elbo_old <- elbo_res[iter - 1]
      tol_1_satisfied <- (abs(elbo_new/elbo_old - 1) < tol)

      if(iter > 2) {

        elbo_old <- elbo_res[iter - 2]
        tol_2_satisfied <- (abs(elbo_new/elbo_old - 1) < tol)
      } else {

        tol_2_satisfied <- FALSE
      }

      tol_satisfied <- (tol_1_satisfied || tol_2_satisfied)

      if(tol_satisfied) {

        converged <- TRUE
      }
    }
  }

  # Get the list of natural parameter vectors:

  return(eta_vec)
}



#' Mean-field variational Bayes (MFVB) inference for univariate or multivariate
#' functional principal component analysis.
#'
#' This function is used to perform FPCA or mFPCA using a MFVB algorithm.
#'
#' @param time_obs  List of vectors or list of lists of vectors containing time
#'                  of observations for univariate or multivariate curves,
#'                  respectively.
#' @param Y List of vectors (for univariate curves) or list of lists of vectors
#'          (for multivariate curves) with the functional data.
#' @param L Number of eigenfunctions (or latent dimensions).
#' @param K Number of O'Sulivan spline functions to be used in the FPCA or mFPCA
#'          algorithms. If set to \code{NULL} will be set according to the rule
#'          of Ruppert (2002), also enforcing K >=7.
#' @param list_hyper Hyperparameter settings constructed using the function
#'         \code{\link{set_hyper}}. If \code{NULL} default hyperparameters will
#'         be used.
#' @param n_g Desired size for dense grid.
#' @param time_g Dense grid provided as a vector of size \code{n_g}. If provided,
#'               then \code{n_g} must be \code{NULL} as will be taken to be
#'               \code{length(time_g)}.
#' @param tol	Tolerance on the relative changes in the ELBO.
#' @param maxit Maximum number of iterations allowed.
#' @param plot_elbo	Boolean indicating whether the values of the ELBO should be displayed during the run.
#' @param Psi_g Reference eigenfunctions (if available, e.g., in simulations)
#'              used to flip the sign of the resulting scores and eigenfunctions.
#' @param verbose Boolean indicating whether messages should be printed during
#'                the run.
#' @param seed User-specified seed for reproducibilty.
#'
#' @return An object containing the resulting MFVB estimates.
#'
#' @seealso  \code{\link{run_vmp_fpca}}, \code{\link{set_hyper}},
#'           \code{\link{flip_sign}}, \code{\link{display_fit}},
#'           \code{\link{display_fit_list}}, \code{\link{display_scores}},
#'           \code{\link{display_eigenfunctions}}
#'
#' @references
#' Nolan, T. H., Goldsmith, J., & Ruppert, D. (2023). Bayesian Functional
#' Principal Components Analysis via Variational Message Passing with Multilevel
#' Extensions. Bayesian Analysis, 1(1), 1-27.
#'
#' Ruppert, D. (2002). Selecting the number of knots for penalized splines.
#' Journal of computational and graphical statistics, 11(4), 735-757.
#'
#' @export
#'
run_mfvb_fpca <- function(time_obs, Y, L, K = NULL,
                          list_hyper = NULL, tol = 1e-5, maxit = 1e4,
                          plot_elbo = FALSE,
                          n_g = 1000, time_g = NULL,
                          Psi_g = NULL, verbose = TRUE, seed = NULL) {

  check_structure(seed, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(seed)) {
    cat(paste0("== Seed set to ", seed, " ==\n\n"))
    set.seed(seed)
  }

  stopifnot(is.list(time_obs))
  stopifnot(is.list(Y))
  stopifnot(isTRUE(all.equal(length(time_obs), length(Y))))

  check_structure(L, "vector", "numeric", 1)
  check_natural(L)

  if (is.null(list_hyper)) {
    list_hyper <- set_hyper()
  } else if (!inherits(list_hyper, "hyper")) {
    stop(paste0("The provided list_hyper must be an object of class ",
                "``hyper''. \n *** you must either use the ",
                "function set_hyper to set your own hyperparameters or ",
                "list_hyper to NULL for automatic choice. ***"))
  }

  sigma_zeta <- list_hyper$sigma_zeta
  # mu_beta <- list_hyper$mu_beta # fixed to zero in mfpca MFVB anyway
  mu_beta <- rep(0, 2)
  Sigma_beta <- list_hyper$Sigma_beta
  A <- list_hyper$A

  check_structure(n_g, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(n_g)) check_natural(n_g)

  check_structure(time_g, "vector", "numeric", null_ok = TRUE)

  check_structure(tol, "vector", "numeric", 1)
  check_positive(tol, eps=.Machine$double.eps)

  check_structure(maxit, "vector", "numeric", 1)
  check_natural(maxit)

  check_structure(plot_elbo, "vector", "logical", 1)

  check_structure(verbose, "vector", "logical", 1)

  N <- length(time_obs) # move to vmp_gauss_fpca and vmp_gauss_mfpca

  if (is.null(names(Y))) {
    stopifnot(is.null(names(time_obs)))
    subj_names <- paste0("subj_", 1:N)
    names(Y) <- names(time_obs) <- subj_names
  } else {
    subj_names <- names(Y)
    stopifnot(isTRUE(all.equal(names(time_obs), subj_names)) | is.null(names(time_obs)))
  }
  names(time_obs) <- subj_names

  format_univ <- ifelse(is.list(time_obs[[1]]), FALSE, TRUE)  # If TRUE, then p = 1, and vmp_gauss_fpca is used.
                                                              # Note that if FALSE and p = 1, vmp_gauss_mfpca will be used,
                                                              # i.e., special case of multivariate algorithm for p=1 gives the univariate model
                                                              # (produces inference equivalent to vmp_gauss_fpca, the difference is the format of the input)

  if (!format_univ) {
    p <- length(time_obs[[1]])
    if (!is.null(Psi_g)) {
      stopifnot(is.list(Psi_g))
    }
  } else {
    p <- 1
  }


  if (is.null(K)) {  # Ruppert (2002) sets a simple default value for K as min(nobs/4,40), where nobs is the number of observations.
    # here since nobs differs for each i, we take the nobs / 4 = round(median(obs_i)/4), and do this for each variable j = 1, ..., p
    # and we enforce that K>=7

    K <- sapply(1:p, function(j) max(round(min(median(sapply(time_obs, function(time_obs_i) length(time_obs_i[[j]]))/4), 40)), 7))

    # if supplied K is such that length(K) = 1, then will be set to K <- rep(K, p)
  } else if (!format_univ) {
    check_structure(K, "vector", "numeric", c(1,p))
    if (length(K)==1) K <- rep(K, p)
  } else {
    check_structure(K, "vector", "numeric", 1)
  }
  check_natural(K)

  grid_obj <- get_grid_objects(time_obs, K, n_g = n_g, time_g = time_g,
                               format_univ = format_univ)

  C <- grid_obj$C
  n_g <- grid_obj$n_g
  time_g <- grid_obj$time_g
  C_g <- grid_obj$C_g

  if (format_univ) {

    # directly includes the orthnogonalisation step, unlike the vmp_gauss_mfpca function
    mfvb_gauss_fpca(maxit, N, L, K, C, Y, sigma_zeta, mu_beta, Sigma_beta, A,
                    tol, plot_elbo, time_g, C_g, Psi_g, verbose)

  } else {

    tmp_names <- lapply(Y, function(Y_i) names(Y_i))
    tmp_names_time <-  lapply(time_obs, function(time_obs_i) names(time_obs_i))
    if (is.null(unlist(unique(tmp_names)))) {
      stopifnot(is.null(unlist(unique(tmp_names_time))))
      var_names <- paste0("variable_", 1:p)
      Y <- lapply(Y, function(Y_i) { names(Y_i) <- var_names; Y_i})
    } else if (!all_same(tmp_names)) {
      stop("Variable names in Y must be the same for all individuals.")
    } else {
      var_names <- names(Y[[1]])
    }

    if (!is.null(unlist(unique(tmp_names_time))) && (!all_same(tmp_names_time) | !isTRUE(all.equal(tmp_names_time[[1]], var_names)))) {
      stop("Variable names for each individual in time_obs must be the same as in Y.")
    }

    time_obs <- lapply(time_obs, function(time_obs_i) { names(time_obs_i) <- var_names; time_obs_i}) # we should probably return this as well.

    # directly includes the orthogonalisation step, unlike the vmp_gauss_mfpca function
    mfvb_gauss_mfpca(maxit, N, p, L, K, C, Y, sigma_zeta, Sigma_beta, A,
                     tol, plot_elbo, var_names, time_g, C_g, Psi_g, verbose)

  }

}


mfvb_gauss_mfpca <- function(maxit, N, p, L, K, C, Y, sigma_zeta,
                             Sigma_beta, A, tol, plot_elbo, var_names, time_g, C_g,
                             Psi_g, verbose, eps = .Machine$double.eps^0.5,
                             debug = TRUE) {


  n <- Reduce(rbind, lapply(Y, function(x) sapply(x, length))) # rows = samples, columns = variables
  sum_n <-  colSums(n)
  n_g <- length(time_g)

  CTC <- vector("list", length = N)
  for(i in 1:N) {

    CTC[[i]] <- vector("list", length = p)
    for(j in 1:p) {

      CTC[[i]][[j]] <- crossprod(C[[i]][[j]])
    }
  }

  # inv_fixed_var <- 1/sigsq_beta*diag(2) # or :
  # inv_fixed_var <- Sigma_beta
  # diag(inv_fixed_var) <- 1/diag(Sigma_beta)
  # more generally (since Sigma_beta is a user-specified parameter):
  inv_Sigma_beta <- solve(Sigma_beta)

  # Set fixed parameters:

  mu_inds <- vector("list", length = p)
  psi_inds <- vector("list", length = p)
  for(j in 1:p) {

    mu_inds[[j]] <- 1:(K[j] + 2)
    psi_inds[[j]] <- split((K[j] + 3):((K[j] + 2)*(L + 1)), rep(1:L, each = K[j] + 2))
  }

  kappa_q_sigsq_eps <- sum_n + 1 # vector of length p
  kappa_q_a_eps <- 2

  kappa_q_sigsq_mu <- K + 1 # vector of length p
  kappa_q_a_mu <- 2

  kappa_q_sigsq_psi <- matrix(rep(K + 1, each = L), nrow = p, byrow = T) # matrix p x L
  kappa_q_a_psi <- 2

  # Initialisation

  E_q_recip_sigsq_eps <- rep(1, p)
  E_q_recip_a_eps <- rep(1, p)

  E_q_recip_sigsq_mu <- rep(1, p)
  E_q_recip_a_mu <- rep(1, p)
  E_q_recip_sigsq_psi <- matrix(1, nrow = p, ncol = L)
  E_q_recip_a_psi <- matrix(1, nrow = p, ncol = L)

  E_q_zeta <- vector("list", length = N)
  Cov_q_zeta <- vector("list", length = N)
  for(i in 1:N) {

    E_q_zeta[[i]] <- rnorm(L, 0, sigma_zeta)
    Cov_q_zeta[[i]] <- diag(L)
  }

  # Iterations:

  elbo_res <- NULL
  converged <- FALSE
  iter <- 0

  while((!converged) & (iter < maxit)) {

    iter <- iter + 1

    if (verbose) cat("Iteration", iter, "\n")

    # Update q(nu):

    Cov_q_nu <- vector("list", length = p)
    inv_Cov_q_nu <- vector("list", length = p)
    E_q_nu <- vector("list", length = p)
    for(j in 1:p) {

      Cov_sum <- 0
      mu_sum <- 0
      for(i in 1:N) {

        E_q_zeta_tilde <- c(1, E_q_zeta[[i]])
        Cov_q_zeta_tilde <- adiag(0, Cov_q_zeta[[i]])
        E_q_zetazetaT_tilde <- Cov_q_zeta_tilde + tcrossprod(E_q_zeta_tilde)

        Cov_sum <- Cov_sum + kronecker(E_q_zetazetaT_tilde, CTC[[i]][[j]])
        mu_sum <- mu_sum + cprod(kronecker(t(E_q_zeta_tilde), C[[i]][[j]]), Y[[i]][[j]])
      }

      E_q_Sigma_mu <- adiag(inv_Sigma_beta, E_q_recip_sigsq_mu[j]*diag(K[j]))
      E_q_Sigma_psi <- vector("list", length = L)
      for(l in 1:L) {

        E_q_Sigma_psi[[l]] <- adiag(inv_Sigma_beta, E_q_recip_sigsq_psi[j, l]*diag(K[j]))
      }
      E_q_inv_Sigma_nu <- adiag(E_q_Sigma_mu, Reduce(adiag, E_q_Sigma_psi))

      inv_Cov_q_nu[[j]] <- E_q_recip_sigsq_eps[j]*Cov_sum + E_q_inv_Sigma_nu
      Cov_q_nu[[j]] <- solve(inv_Cov_q_nu[[j]])
      E_q_nu[[j]] <- E_q_recip_sigsq_eps[j]*as.vector(Cov_q_nu[[j]] %*% mu_sum)
    }

    E_q_V <- vector("list", length = p)
    E_q_V_psi <- vector("list", length = p)
    for(j in 1:p) {

      E_q_V_psi[[j]] <- matrix(E_q_nu[[j]][-mu_inds[[j]]], K[j] + 2, L)
      E_q_V[[j]] <- matrix(E_q_nu[[j]], K[j] + 2, L + 1)
    }

    E_q_H_psi <- vector("list", length = N)
    E_q_h_nu <- vector("list", length = N)
    E_q_h_mu <- vector("list", length = N)
    E_q_H_nu <- vector("list", length = N)
    for(i in 1:N) {

      E_q_H_psi[[i]] <- vector("list", length = p)
      E_q_h_nu[[i]] <- vector("list", length = p)
      E_q_h_mu[[i]] <- vector("list", length = p)
      E_q_H_nu[[i]] <- vector("list", length = p)
      for(j in 1:p) {

        Cov_j <- Cov_q_nu[[j]]
        E_j <- E_q_nu[[j]]

        E_q_H_psi[[i]][[j]] <- matrix(NA, L, L)
        E_q_h_nu[[i]][[j]] <- rep(NA, L)
        for(l_1 in 1:L) {

          inds_1 <- psi_inds[[j]][[l_1]]
          for(l_2 in 1:L) {

            inds_2 <- psi_inds[[j]][[l_2]]

            tr_term <- tr(Cov_j[inds_2, inds_1] %*% CTC[[i]][[j]])
            cprod_term <- cprod(E_j[inds_1], CTC[[i]][[j]] %*% E_j[inds_2])
            E_q_H_psi[[i]][[j]][l_1, l_2] <- tr_term + cprod_term
          }

          tr_term <- tr(Cov_j[mu_inds[[j]], inds_1] %*% CTC[[i]][[j]])
          cprod_term <- cprod(E_j[inds_1], CTC[[i]][[j]] %*% E_j[mu_inds[[j]]])
          E_q_h_nu[[i]][[j]][l_1] <- tr_term + cprod_term
        }

        tr_term <- tr(Cov_j[mu_inds[[j]], mu_inds[[j]]] %*% CTC[[i]][[j]])
        cprod_term <- cprod(E_j[mu_inds[[j]]], CTC[[i]][[j]] %*% E_j[mu_inds[[j]]])
        E_q_h_mu[[i]][[j]] <- tr_term + cprod_term

        top_block <- cbind(E_q_h_mu[[i]][[j]], t(E_q_h_nu[[i]][[j]]))
        bottom_block <- cbind(E_q_h_nu[[i]][[j]], E_q_H_psi[[i]][[j]])
        E_q_H_nu[[i]][[j]] <- rbind(top_block, bottom_block)
      }
    }

    # For i = 1, ..., N, update q(zeta[i]):

    Cov_q_zeta <- vector("list", length = N)
    inv_Cov_q_zeta <- vector("list", length = N)
    E_q_zeta <- vector("list", length = N)
    for(i in 1:N) {

      Cov_sum <- 0
      mu_sum <- 0
      for(j in 1:p) {

        Cov_sum <- Cov_sum + E_q_recip_sigsq_eps[j]*E_q_H_psi[[i]][[j]]

        E_q_Psi_j <- C[[i]][[j]] %*% E_q_V_psi[[j]]
        freq_scores <- cprod(E_q_Psi_j, Y[[i]][[j]]) - E_q_h_nu[[i]][[j]]
        mu_sum <- mu_sum + E_q_recip_sigsq_eps[j]*freq_scores
      }

      inv_Cov_q_zeta[[i]] <- Cov_sum + 1/sigma_zeta^2*diag(L)
      Cov_q_zeta[[i]] <- solve(inv_Cov_q_zeta[[i]])
      E_q_zeta[[i]] <- as.vector(Cov_q_zeta[[i]] %*% mu_sum)
    }

    # For j = 1, ..., p, update q(sigsq_eps[j]):

    lambda_q_sigsq_eps <- rep(0, p)
    E_q_recip_sigsq_eps <- rep(NA, p)
    for(j in 1:p) {

      for(i in 1:N) {

        E_q_zeta_tilde <- c(1, E_q_zeta[[i]])
        Cov_q_zeta_tilde <- adiag(0, Cov_q_zeta[[i]])
        E_q_zetazetaT_tilde <- Cov_q_zeta_tilde + tcrossprod(E_q_zeta_tilde)

        term_1 <- cprod(Y[[i]][[j]])
        term_2 <- -2*cprod(Y[[i]][[j]], C[[i]][[j]] %*% E_q_V[[j]] %*% E_q_zeta_tilde)
        term_3 <- tr(E_q_zetazetaT_tilde %*% E_q_H_nu[[i]][[j]])

        lambda_q_sigsq_eps[j] <- lambda_q_sigsq_eps[j] + term_1 + term_2 + term_3
      }

      lambda_q_sigsq_eps[j] <- lambda_q_sigsq_eps[j] + E_q_recip_a_eps[j]
      E_q_recip_sigsq_eps[j] <- kappa_q_sigsq_eps[j] / lambda_q_sigsq_eps[j]
    }

    # For j = 1, ..., p, update q(a_eps[j]):

    lambda_q_a_eps <- E_q_recip_sigsq_eps + 1/A^2
    E_q_recip_a_eps <- kappa_q_a_eps / lambda_q_a_eps

    # For j = 1, ..., p, update q(sigsq_mu[j]):

    lambda_q_sigsq_mu <- rep(NA, p)
    for(j in 1:p) {

      u_inds <- mu_inds[[j]][3:(K[j] + 2)]
      tr_term <- tr(Cov_q_nu[[j]][u_inds, u_inds])
      cprod_term <- cprod(E_q_nu[[j]][u_inds])
      lambda_q_sigsq_mu[j] <- tr_term + cprod_term + E_q_recip_a_mu[j]
    }
    E_q_recip_sigsq_mu <- kappa_q_sigsq_mu / lambda_q_sigsq_mu

    # Update q(a_mu):

    lambda_q_a_mu <- E_q_recip_sigsq_mu + 1/A^2
    E_q_recip_a_mu <- kappa_q_a_mu / lambda_q_a_mu

    # For j = 1, ..., p and l = 1, ..., L, update q(sigsq_psi[j, l]):

    lambda_q_sigsq_psi <- matrix(NA, p, L)
    for(l in 1:L) {

      for(j in 1:p) {

        u_inds <- psi_inds[[j]][[l]]
        tr_term <- tr(Cov_q_nu[[j]][u_inds, u_inds])
        cprod_term <- cprod(E_q_nu[[j]][u_inds])
        lambda_q_sigsq_psi[j, l] <- tr_term + cprod_term + E_q_recip_a_psi[j, l]
      }
    }
    E_q_recip_sigsq_psi <- kappa_q_sigsq_psi / lambda_q_sigsq_psi

    # For j = 1, ..., p and l = 1, ..., L, update q(a_psi[j, l]):

    lambda_q_a_psi <- E_q_recip_sigsq_psi + 1/A^2
    E_q_recip_a_psi <- kappa_q_a_psi/lambda_q_a_psi


    ###
    # E_q_log_sigsq objects (for ELBO)
    #

    E_q_log_sigsq_eps <- compute_mu_q_log(kappa_q_sigsq_eps, lambda_q_sigsq_eps) # length p
    E_q_log_sigsq_mu <- compute_mu_q_log(kappa_q_sigsq_mu, lambda_q_sigsq_mu) # length p
    E_q_log_sigsq_psi <- compute_mu_q_log(kappa_q_sigsq_psi, lambda_q_sigsq_psi) # matrix nrow = p, ncol = L

    E_q_log_a_eps <- compute_mu_q_log(kappa_q_a_eps, lambda_q_a_eps) # length p
    E_q_log_a_mu <- compute_mu_q_log(kappa_q_a_mu, lambda_q_a_mu) # length p
    E_q_log_a_psi <- compute_mu_q_log(kappa_q_a_psi, lambda_q_a_psi) # matrix nrow = p, ncol = L


    # ELBO:
    #
    elbo_y <- sum(sapply(1:p, function(j) { e_y(N, sum_n[j],
                                                lapply(Y, function(Y_i) Y_i[[j]]),
                                                lapply(C, function(C_i) C_i[[j]]),
                                                E_q_V[[j]],
                                                lapply(E_q_H_nu, function(E_q_H_nu_i) E_q_H_nu_i[[j]]),
                                                E_q_zeta, Cov_q_zeta,
                                                E_q_log_sigsq_eps[j], E_q_recip_sigsq_eps[j])}))

    elbo_nu <- sum(sapply(1:p, function(j) {e_nu(K[j], L, mu_beta = rep(0, 2), # assumed to be 0 in the mFPCA version of mfvb
                                                 inv_Sigma_beta, # hyperparameter
                                                 E_q_nu[[j]], Cov_q_nu[[j]], inv_Cov_q_nu[[j]],
                                                 E_q_recip_sigsq_mu[j], E_q_recip_sigsq_psi[j,],
                                                 E_q_log_sigsq_mu[j], E_q_log_sigsq_psi[j,])}))


    elbo_zeta <- e_zeta(N, inv_Sigma_zeta =  1/sigma_zeta^2*diag(L), # hyperparameter
                        E_q_zeta, Cov_q_zeta, inv_Cov_q_zeta)

    elbo_sigsq_eps <- sum(sapply(1:p, function(j) {e_sigsq(E_q_recip_sigsq_eps[j], E_q_log_sigsq_eps[j],
                                                           E_q_recip_a_eps[j], E_q_log_a_eps[j],
                                                           kappa_q_sigsq_eps[j], lambda_q_sigsq_eps[j])}))

    elbo_a_eps <- sum(sapply(1:p, function(j) {e_a(E_q_recip_a_eps[j], E_q_log_a_eps[j],
                                                   A, kappa_q_a_eps, lambda_q_a_eps[j])}))

    elbo_sigsq_mu <- sum(sapply(1:p, function(j) {e_sigsq(E_q_recip_sigsq_mu[j], E_q_log_sigsq_mu[j],
                                                          E_q_recip_a_mu[j], E_q_log_a_mu[j],
                                                          kappa_q_sigsq_mu[j], lambda_q_sigsq_mu[j])}))

    elbo_a_mu <- sum(sapply(1:p, function(j) {e_a(E_q_recip_a_mu[j], E_q_log_a_mu[j],
                                                  A, kappa_q_a_mu, lambda_q_a_mu[j])}))

    elbo_sigsq_psi <- sum(sapply(1:p, function(j) {sum(e_sigsq(E_q_recip_sigsq_psi[j,], E_q_log_sigsq_psi[j,],
                                                               E_q_recip_a_psi[j,], E_q_log_a_psi[j,],
                                                               kappa_q_sigsq_psi[j,], lambda_q_sigsq_psi[j,]))})) # sum because vector of length L (L eigenfunctions)

    elbo_a_psi <- sum(sapply(1:p, function(j) {sum(e_a(E_q_recip_a_psi[j,], E_q_log_a_psi[j,],
                                                       A, kappa_q_a_psi, lambda_q_a_psi[j,]))})) # sum because vector of length L (L eigenfunctions)

    elbo_new <- elbo_y + elbo_nu + elbo_zeta + elbo_sigsq_eps + elbo_a_eps +
      elbo_sigsq_mu + elbo_a_mu + elbo_sigsq_psi + elbo_a_psi

    elbo_res <- c(elbo_res, elbo_new)

    if(plot_elbo) {

      plot(1:iter, elbo_res, pch=16, cex=0.4, xlab="iterations", ylab="ELBO")
    }

    if(iter > 1) {

      elbo_old <- elbo_res[iter - 1]

      if (debug && elbo_new + eps < elbo_old)
        stop("ELBO not increasing monotonically. Exit. ")

      converged <- (abs(elbo_new/elbo_old - 1) < tol)

    }

  }

  # Orthogonalisation:

  E_q_Zeta <- Reduce(rbind, E_q_zeta)
  E_q_mu <- vector("list", length = p)
  E_q_Psi <- vector("list", length = p)
  for(j in 1:p) {

    E_q_V <- matrix(E_q_nu[[j]], K[j] + 2, L + 1)
    gbl_post <- C_g[[j]] %*% E_q_V
    E_q_mu[[j]] <- gbl_post[, 1]
    E_q_Psi[[j]] <- gbl_post[, 1:L + 1]
  }
  E_q_mu <- Reduce(c, E_q_mu)
  E_q_Psi <- Reduce(rbind, E_q_Psi)

  svd_Psi <- svd(E_q_Psi)
  U_psi <- svd_Psi$u
  D_psi <- diag(svd_Psi$d)
  V_psi <- svd_Psi$v

  zeta_rotn <- t(E_q_Zeta %*% V_psi %*% D_psi)
  C_zeta <- cov(t(zeta_rotn))
  eigen_C <- eigen(C_zeta)
  Q <- eigen_C$vectors
  Lambda <- diag(eigen_C$values + 1e-10)

  Psi_tilde <- U_psi %*% Q %*% sqrt(Lambda)
  Zeta_tilde <- crossprod(zeta_rotn, Q %*% solve(sqrt(Lambda)))

  mu_hat <- split(E_q_mu, rep(1:p, each = n_g))

  Psi_hat <- matrix(NA, p*n_g, L)
  Zeta_hat <- matrix(NA, N, L)
  norm_vec <- rep(NA, L)
  for(l in 1:L) {

    psi_l <- split(Psi_tilde[, l], rep(1:p, each = n_g))
    norm_vec[l] <- sqrt(sum(sapply(psi_l, function(x) trapint(time_g, x^2))))
    Psi_hat[, l] <- Psi_tilde[, l]/norm_vec[l]
    Zeta_hat[, l] <- norm_vec[l]*Zeta_tilde[, l]

    if(!is.null(Psi_g)) {

      Psi_g_comb <- vector("list", length = p)
      for(j in 1:p) {

        Psi_g_comb[[j]] <- Psi_g[[j]][, l]
      }
      Psi_g_comb <- Reduce(c, Psi_g_comb)

      inner_prod_sign <- sign(cprod(Psi_g_comb, Psi_hat[, l]))
      Psi_hat[, l] <- inner_prod_sign*Psi_hat[, l]
      Zeta_hat[, l] <- inner_prod_sign*Zeta_hat[, l]

    }
  }
  colnames(Zeta_hat) <- paste0("FPC_", 1:L)
  rownames(Zeta_hat) <- paste0("subj_", 1:N)

  Psi_hat <- lapply(split(Psi_hat, rep(1:p, each = n_g)), matrix, nrow = n_g, ncol = L)

  Cov_zeta_hat <- vector("list", length = N)
  rotn_mat <- V_psi %*% D_psi %*% Q %*% solve(sqrt(Lambda)) %*% diag(norm_vec)
  for(i in 1:N) {

    Cov_zeta_hat[[i]] <- crossprod(rotn_mat, Cov_q_zeta[[i]] %*% rotn_mat)
  }

  # Store the results:

  zeta_summary <- list_zeta_ellipse <- vector("list", length=N)

  for(i in 1:N) {

    zeta_mean <- Zeta_hat[i,][1:2]

    zeta_ellipse <- ellipse(
      Cov_zeta_hat[[i]][1:2, 1:2],
      centre=zeta_mean,
      level=0.95
    )

    list_zeta_ellipse[[i]] <- zeta_ellipse

    zeta_summary[[i]] <- list(zeta_mean, zeta_ellipse)
    names(zeta_summary[[i]]) <- c("mean", "credible boundary")
  }

  gbl_hat <- vector("list", length = L + 1)
  list_Psi_hat <- vector("list", length = L)
  gbl_hat[[1]] <- mu_hat <- as.matrix(Reduce(cbind, mu_hat)) # now deals with the case p = 1
  colnames(gbl_hat[[1]]) <- colnames(mu_hat) <- var_names
  for(l in 1:L) {

    gbl_hat[[l+1]] <- list_Psi_hat[[l]] <- matrix(NA, n_g, p)
    for(j in 1:p) {

      gbl_hat[[l+1]][, j] <- list_Psi_hat[[l]][,j] <- Psi_hat[[j]][, l]
    }
    colnames(gbl_hat[[l+1]]) <- colnames(list_Psi_hat[[l]]) <- var_names
  }
  names(list_Psi_hat) <- paste0("FPC_", 1:L)

  E_q_sigsq_eps <- lambda_q_sigsq_eps / (sum(n[, j]) - 1)

  Y_summary <- Y_hat <- Y_low <- Y_upp <- vector("list", length = N)
  for(i in 1:N) {

    Y_summary[[i]] <- Y_hat[[i]] <- Y_low[[i]] <- Y_upp[[i]] <- vector("list", length = p)
    for(j in 1:p) {

      Y_hat_ij <- mu_hat[,j] + Psi_hat[[j]] %*% Zeta_hat[i, ]
      sd_vec_ij <- sqrt(diag(tcrossprod(Psi_hat[[j]] %*% Cov_zeta_hat[[i]], Psi_hat[[j]])) + E_q_sigsq_eps[j])

      Y_summary[[i]][[j]] <- matrix(NA, n_g, 3)

      Y_summary[[i]][[j]][, 1] <- Y_hat_ij + qnorm(0.025)*sd_vec_ij
      Y_summary[[i]][[j]][, 2] <- Y_hat_ij
      Y_summary[[i]][[j]][, 3] <- Y_hat_ij + qnorm(0.975)*sd_vec_ij
      Y_hat[[i]][[j]] <- Y_hat_ij
      Y_low[[i]][[j]] <- Y_hat_ij + qnorm(0.025)*sd_vec_ij
      Y_upp[[i]][[j]] <- Y_hat_ij + qnorm(0.975)*sd_vec_ij
    }
  }

  create_named_list(time_g, K,
                    Y_summary, Y_hat, Y_low, Y_upp,
                    gbl_hat, mu_hat, list_Psi_hat,
                    Zeta_hat, Cov_zeta_hat, list_zeta_ellipse)


}



mfvb_gauss_fpca <- function(maxit, N, L, K, C, Y, sigma_zeta, mu_beta,
                            Sigma_beta, A, tol, plot_elbo, time_g, C_g,
                            Psi_g, verbose, eps = .Machine$double.eps^0.5,
                            debug = TRUE) {

  n_g <- length(time_g)
  T_vec <- sapply(Y, length)
  sum_T <- sum(T_vec)

  inv_Sigma_zeta <- solve(sigma_zeta^2*diag(L))
  inv_Sigma_beta <- solve(Sigma_beta)

  kappa_q_sigsq_mu <- K + 1
  mu_q_recip_sigsq_mu <- 1

  kappa_q_a_mu <- 2
  mu_q_recip_a_mu <- 1

  kappa_q_sigsq_psi <- K + 1
  mu_q_recip_sigsq_psi <- rep(1, L)

  kappa_q_a_psi <- 2
  mu_q_recip_a_psi <- rep(1, L)

  mu_q_zeta <- vector("list", length=N)
  Sigma_q_zeta <- vector("list", length=N)
  for(i in 1:N) {

    mu_q_zeta[[i]] <- rnorm(L, 0, sigma_zeta)
    Sigma_q_zeta[[i]] <- diag(L)
  }

  kappa_q_sigsq_eps <- sum_T + 1
  mu_q_recip_sigsq_eps <- 1

  kappa_q_a_eps <- 2
  mu_q_recip_a_eps <- 1

  elbo_res <- NULL
  converged <- FALSE
  iter <- 0

  while((!converged) & (iter < maxit)) {

    iter <- iter + 1

    if (verbose) cat("Iteration", iter, "\n")

    # Update q(nu):

    E_q_inv_Sigma_nu <- vector("list", length=L+1)
    E_q_inv_Sigma_nu[[1]] <- blkdiag(
      inv_Sigma_beta,
      mu_q_recip_sigsq_mu*diag(K)
    )
    mu_nu <- vector("list", length=L+1)
    mu_nu[[1]] <- c(mu_beta, rep(0, K))
    for(l in 1:L) {

      E_q_inv_Sigma_nu[[l+1]] <- blkdiag(
        inv_Sigma_beta,
        mu_q_recip_sigsq_psi[l]*diag(K)
      )

      mu_nu[[l+1]] <- c(mu_beta, rep(0, K))
    }
    E_q_inv_Sigma_nu <- Reduce(blkdiag, E_q_inv_Sigma_nu)
    mu_nu <- Reduce(c, mu_nu)

    sum_term_Sigma <- 0
    sum_term_mu <- 0
    for(i in 1:N) {

      mu_q_zeta_tilde <- c(1, mu_q_zeta[[i]])
      Sigma_q_zeta_tilde <- blkdiag(matrix(0), Sigma_q_zeta[[i]])

      M_q_zzT_tilde <- Sigma_q_zeta_tilde + tcrossprod(mu_q_zeta_tilde)
      sum_val <- kronecker(M_q_zzT_tilde, crossprod(C[[i]]))
      sum_term_Sigma <- sum_term_Sigma + sum_val

      sum_val <- cprod(kronecker(t(mu_q_zeta_tilde), C[[i]]), Y[[i]])
      sum_term_mu <- sum_term_mu + sum_val
    }
    inv_Sigma_q_nu <- mu_q_recip_sigsq_eps*sum_term_Sigma + E_q_inv_Sigma_nu
    Sigma_q_nu <- solve(inv_Sigma_q_nu)
    mu_q_nu <- as.vector(
      Sigma_q_nu%*%(
        mu_q_recip_sigsq_eps*sum_term_mu
        + as.vector(E_q_inv_Sigma_nu%*%mu_nu)
      )
    )

    mu_inds <- 1:(K+2)
    mu_q_nu_mu <- mu_q_nu[mu_inds]
    Sigma_q_nu_mu <- Sigma_q_nu[mu_inds, mu_inds]

    psi_inds <- vector("list", length=L)
    mu_q_nu_psi <- vector("list", length=L)
    Sigma_q_nu_psi <- vector("list", length=L)
    ind_start <- K+3
    for(l in 1:L) {

      ind_end <- ind_start + K + 1
      psi_inds[[l]] <- ind_start:ind_end

      mu_q_nu_psi[[l]] <- mu_q_nu[psi_inds[[l]]]
      Sigma_q_nu_psi[[l]] <- Sigma_q_nu[psi_inds[[l]], psi_inds[[l]]]

      ind_start <- ind_end + 1
    }

    # For 1 <= i <= N, update q(zeta[i]):

    M_q_H <- vector("list", length=N)
    mu_q_h <- vector("list", length=N)
    for(i in 1:N) {

      M_q_H[[i]] <- matrix(NA, L, L)
      mu_q_h[[i]] <- rep(NA, L)
      for(l_i in 1:L) {

        i_inds <- psi_inds[[l_i]]
        mu_q_nu_psi_i <- mu_q_nu[i_inds]

        Cov_q_nu_mu_psi <- Sigma_q_nu[mu_inds, i_inds]

        tr_term_h <- tr(Cov_q_nu_mu_psi%*%crossprod(C[[i]]))
        cprod_term_h <- cprod(mu_q_nu_psi_i, crossprod(C[[i]])%*%mu_q_nu_mu)
        mu_q_h[[i]][l_i] <- tr_term_h + cprod_term_h

        for(l_j in 1:L) {

          j_inds <- psi_inds[[l_j]]
          mu_q_nu_psi_j <- mu_q_nu[j_inds]

          Cov_q_nu_psi_ji <- Sigma_q_nu[j_inds, i_inds]

          tr_term <- tr(Cov_q_nu_psi_ji%*%crossprod(C[[i]]))
          cprod_term <- cprod(mu_q_nu_psi_i, crossprod(C[[i]])%*%mu_q_nu_psi_j)
          M_q_H[[i]][l_i, l_j] <- tr_term + cprod_term
        }
      }
    }

    M_q_V_psi <- Reduce(cbind, mu_q_nu_psi)

    mu_q_zeta <- vector("list", length=N)
    inv_Sigma_q_zeta <- vector("list", length=N)
    Sigma_q_zeta <- vector("list", length=N)
    for(i in 1:N) {

      M_q_Psi <- C[[i]]%*%M_q_V_psi
      centre_vec <- cprod(M_q_Psi, Y[[i]]) - mu_q_h[[i]]

      inv_Sigma_q_zeta[[i]] <- mu_q_recip_sigsq_eps*M_q_H[[i]] + inv_Sigma_zeta
      Sigma_q_zeta[[i]] <- solve(inv_Sigma_q_zeta[[i]])
      mu_q_zeta[[i]] <- mu_q_recip_sigsq_eps*as.vector(Sigma_q_zeta[[i]]%*%centre_vec)
    }

    # Update q(sigsq_eps):

    M_q_V <- cbind(mu_q_nu_mu, M_q_V_psi)

    lambda_q_sigsq_eps <- mu_q_recip_a_eps
    M_q_H_tilde <- vector("list", length=N)
    for(i in 1:N) {

      tr_term <- tr(Sigma_q_nu_mu%*%crossprod(C[[i]]))
      cprod_term <- cprod(mu_q_nu_mu, crossprod(C[[i]])%*%mu_q_nu_mu)
      mu_q_h_mu <- tr_term + cprod_term

      top_mat <- c(mu_q_h_mu, mu_q_h[[i]])
      bottom_mat <- cbind(mu_q_h[[i]], M_q_H[[i]])
      M_q_H_tilde[[i]] <- rbind(top_mat, bottom_mat)

      mu_q_zeta_tilde <- c(1, mu_q_zeta[[i]])
      Sigma_q_zeta_tilde <- blkdiag(matrix(0), Sigma_q_zeta[[i]])
      M_q_zzT_tilde <- Sigma_q_zeta_tilde + tcrossprod(mu_q_zeta_tilde)

      sum_val <- cprod(Y[[i]])
      sum_val <- sum_val - 2*cprod(C[[i]]%*%M_q_V%*%mu_q_zeta_tilde, Y[[i]])
      sum_val <- sum_val + tr(M_q_zzT_tilde%*%M_q_H_tilde[[i]])
      lambda_q_sigsq_eps <- lambda_q_sigsq_eps + sum_val
    }

    mu_q_recip_sigsq_eps <- kappa_q_sigsq_eps/lambda_q_sigsq_eps

    # Update q(a_eps):

    lambda_q_a_eps <- mu_q_recip_sigsq_eps + 1/A^2
    mu_q_recip_a_eps <- kappa_q_a_eps/lambda_q_a_eps

    # Update q(sigsq_mu):

    mu_q_u_mu <- mu_q_nu_mu[-c(1:2)]
    Sigma_q_u_mu <- Sigma_q_nu_mu[-c(1:2), -c(1:2)]

    lambda_q_sigsq_mu <- cprod(mu_q_u_mu)
    lambda_q_sigsq_mu <- lambda_q_sigsq_mu + tr(Sigma_q_u_mu)
    lambda_q_sigsq_mu <- lambda_q_sigsq_mu + mu_q_recip_a_mu
    mu_q_recip_sigsq_mu <- kappa_q_sigsq_mu/lambda_q_sigsq_mu

    # Update q(a_mu):

    lambda_q_a_mu <- mu_q_recip_sigsq_mu + 1/A^2
    mu_q_recip_a_mu <- kappa_q_a_mu/lambda_q_a_mu

    # Update q(sigsq_psi):

    lambda_q_sigsq_psi <- rep(NA, L)
    for(l in 1:L) {

      mu_q_u_psi <- mu_q_nu_psi[[l]][-c(1:2)]
      Sigma_q_u_psi <- Sigma_q_nu_psi[[l]][-c(1:2), -c(1:2)]

      lambda_q_sigsq_psi[l] <- cprod(mu_q_u_psi)
      lambda_q_sigsq_psi[l] <- lambda_q_sigsq_psi[l] + tr(Sigma_q_u_psi)
      lambda_q_sigsq_psi[l] <- lambda_q_sigsq_psi[l] + mu_q_recip_a_psi[l]
    }

    mu_q_recip_sigsq_psi <- kappa_q_sigsq_psi/lambda_q_sigsq_psi

    # Update q(a_psi):

    lambda_q_a_psi <- mu_q_recip_sigsq_psi + 1/A^2
    mu_q_recip_a_psi <- kappa_q_a_psi/lambda_q_a_psi


    # mu_q_log_sigsq objects (for ELBO)
    #
    mu_q_log_sigsq_eps <- compute_mu_q_log(kappa_q_sigsq_eps, lambda_q_sigsq_eps)
    mu_q_log_sigsq_mu <- compute_mu_q_log(kappa_q_sigsq_mu, lambda_q_sigsq_mu)
    mu_q_log_sigsq_psi <- compute_mu_q_log(kappa_q_sigsq_psi, lambda_q_sigsq_psi)

    mu_q_log_a_eps <- compute_mu_q_log(kappa_q_a_eps, lambda_q_a_eps)
    mu_q_log_a_mu <- compute_mu_q_log(kappa_q_a_mu, lambda_q_a_mu)
    mu_q_log_a_psi <- compute_mu_q_log(kappa_q_a_psi, lambda_q_a_psi)


    # ELBO:
    #
    elbo_y <- e_y(N, sum_T, Y, C, M_q_V, M_q_H_tilde, mu_q_zeta, Sigma_q_zeta,
                  mu_q_log_sigsq_eps, mu_q_recip_sigsq_eps)

    elbo_nu <- e_nu(K, L, mu_beta, # hyperparamter
                    inv_Sigma_beta, # hyperparameter
                    mu_q_nu, Sigma_q_nu, inv_Sigma_q_nu,
                    mu_q_recip_sigsq_mu, mu_q_recip_sigsq_psi,
                    mu_q_log_sigsq_mu, mu_q_log_sigsq_psi)


    elbo_zeta <- e_zeta(N, inv_Sigma_zeta, # hyperparameter
                        mu_q_zeta, Sigma_q_zeta, inv_Sigma_q_zeta)

    elbo_sigsq_eps <- e_sigsq(mu_q_recip_sigsq_eps, mu_q_log_sigsq_eps,
                              mu_q_recip_a_eps, mu_q_log_a_eps,
                              kappa_q_sigsq_eps, lambda_q_sigsq_eps)

    elbo_a_eps <- e_a(mu_q_recip_a_eps, mu_q_log_a_eps,
                      A, kappa_q_a_eps, lambda_q_a_eps)

    elbo_sigsq_mu <- e_sigsq(mu_q_recip_sigsq_mu, mu_q_log_sigsq_mu,
                             mu_q_recip_a_mu, mu_q_log_a_mu,
                             kappa_q_sigsq_mu, lambda_q_sigsq_mu)

    elbo_a_mu <- e_a(mu_q_recip_a_mu, mu_q_log_a_mu,
                     A, kappa_q_a_mu, lambda_q_a_mu)

    elbo_sigsq_psi <- sum(e_sigsq(mu_q_recip_sigsq_psi, mu_q_log_sigsq_psi,
                                  mu_q_recip_a_psi, mu_q_log_a_psi,
                                  kappa_q_sigsq_psi, lambda_q_sigsq_psi)) # sum because vector of length L (L eigenfunctions)

    elbo_a_psi <- sum(e_a(mu_q_recip_a_psi, mu_q_log_a_psi,
                          A, kappa_q_a_psi, lambda_q_a_psi)) # sum because vector of length L (L eigenfunctions)

    elbo_new <- elbo_y + elbo_nu + elbo_zeta + elbo_sigsq_eps + elbo_a_eps +
      elbo_sigsq_mu + elbo_a_mu + elbo_sigsq_psi + elbo_a_psi

    elbo_res <- c(elbo_res, elbo_new)

    if(plot_elbo) {

      plot(1:iter, elbo_res, pch=16, cex=0.4, xlab="iterations", ylab="ELBO")
    }

    if(iter > 1) {

      elbo_old <- elbo_res[iter - 1]

      if (debug && elbo_new + eps < elbo_old)
        stop("ELBO not increasing monotonically. Exit. ")

      converged <- (abs(elbo_new/elbo_old - 1) < tol)

    }

  }

  # Orthogonal decomposition:

  mu_q_mu <- as.vector(C_g%*%mu_q_nu_mu)

  M_q_V_psi <- Reduce(cbind, mu_q_nu_psi)
  M_q_Psi <- C_g %*% M_q_V_psi

  M_q_Zeta <- Reduce(rbind, mu_q_zeta)

  one_N <- rep(1, N)
  mu_mat <- tcrossprod(mu_q_mu, one_N)
  Y_mat <- mu_mat + tcrossprod(M_q_Psi, M_q_Zeta)

  M_q_Psi_svd <- svd(M_q_Psi)
  U_orth <- M_q_Psi_svd$u
  D_diag <- diag(M_q_Psi_svd$d)
  V_orth <- M_q_Psi_svd$v

  M_q_Zeta_rotn <- M_q_Zeta%*%V_orth%*%D_diag

  eigen_M_q_Zeta_shift <- eigen(cov(M_q_Zeta_rotn))
  Q <- eigen_M_q_Zeta_shift$vectors
  Lambda <- diag(eigen_M_q_Zeta_shift$values + 1e-10)
  Lambda_inv <- diag(1/(eigen_M_q_Zeta_shift$values + 1e-10))
  S <- Q%*%sqrt(Lambda)
  S_inv <- tcrossprod(sqrt(Lambda_inv), Q)

  Psi_hat <- U_orth%*%S
  Zeta_hat <- tcrossprod(M_q_Zeta_rotn, S_inv)

  # print(all.equal(tcrossprod(Psi_hat, Zeta_hat),
  #                 tcrossprod(M_q_Psi, M_q_Zeta)))

  norm_const <- rep(NA, L)
  for(l in 1:L) {

    norm_const[l] <- sqrt(trapint(time_g, (Psi_hat[,l])^2))
    if(norm_const[l]!=0) {

      Psi_hat[,l] <- Psi_hat[,l]/norm_const[l]
      Zeta_hat[,l] <- norm_const[l]*Zeta_hat[,l]

      if(!is.null(Psi_g)) {

        cprod_sign <- sign(cprod(Psi_hat[,l], Psi_g[,l]))
        Psi_hat[,l] <- cprod_sign*Psi_hat[,l]
        Zeta_hat[,l] <- cprod_sign*Zeta_hat[,l]

      }
    }
  }
  colnames(Zeta_hat) <- paste0("FPC_", 1:L)
  rownames(Zeta_hat) <- paste0("subj_", 1:N)

  scale_mat <- diag(norm_const)
  Cov_zeta_hat <- vector("list", length = N)
  for(i in 1:N) {

    mat_transform <- S_inv%*%tcrossprod(D_diag, V_orth)
    Cov_zeta_dot_hat <- tcrossprod(mat_transform%*%Sigma_q_zeta[[i]], mat_transform)
    Cov_zeta_hat[[i]] <- tcrossprod(scale_mat %*% Cov_zeta_dot_hat, scale_mat)
  }

  zeta_summary <- list_zeta_ellipse <-  vector("list", length=N)
  for(i in 1:N) {

    zeta_mean <- Zeta_hat[i,][1:2]

    zeta_ellipse <- ellipse(
      Cov_zeta_hat[[i]][1:2, 1:2],
      centre = zeta_mean,
      level = 0.95
    )

    list_zeta_ellipse[[i]] <- zeta_ellipse

    zeta_summary[[i]] <- list(zeta_mean, zeta_ellipse)
    names(zeta_summary[[i]]) <- c("mean", "credible boundary")
  }

  gbl_hat <- vector("list", length = L + 1)
  gbl_hat[[1]] <- mu_hat <- mu_q_mu
  list_Psi_hat <- Psi_hat
  for(l in 1:L) {
    gbl_hat[[l+1]] <- list_Psi_hat[[l]]
  }

  E_q_sigsq_eps <- lambda_q_sigsq_eps / (kappa_q_sigsq_eps - 2)
  Y_summary <- Y_hat <- Y_low <- Y_upp <- vector("list", length = N)
  for(i in 1:N) {

    sd_vec_i <- sqrt(diag(tcrossprod(Psi_hat%*%Cov_zeta_hat[[i]], Psi_hat)) + E_q_sigsq_eps)

    Y_summary[[i]] <- matrix(NA, nrow=n_g, ncol=3)
    Y_summary[[i]][,1] <- Y_mat[,i] + qnorm(0.025)*sd_vec_i
    Y_summary[[i]][,2] <- Y_mat[,i]
    Y_summary[[i]][,3] <- Y_mat[,i] + qnorm(0.975)*sd_vec_i
    Y_hat[[i]] <- Y_mat[,i]
    Y_low[[i]] <- Y_mat[,i] + qnorm(0.025)*sd_vec_i
    Y_upp[[i]] <- Y_mat[,i] + qnorm(0.975)*sd_vec_i
  }


  create_named_list(time_g, K,
                    Y_summary, Y_hat, Y_low, Y_upp,
                    gbl_hat, mu_hat, list_Psi_hat,
                    Zeta_hat, Cov_zeta_hat, list_zeta_ellipse)

}

