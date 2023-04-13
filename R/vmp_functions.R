############### R library: vmp_functions.r ###############

# A library that stores VMP fragment and sufficient
# statistic expectation updates:

# Created: 08 JUN 2020
# Last Updated: 09 DEC 2022

# require(matrixcalc)
# require(pracma)
# require(MASS)
# require(ellipse)
# source("trapint.r")
#source("logistic.r")

# LIST  OF  FUNCTIONS:

# Fragment Updates:
# igw_prior_frag
# iter_igw_frag
# gauss_prior_frag
# gauss_lik_frag
# logistic_lik_frag
# gauss_pen_frag
# mult_gauss_pen_frag
# fpc_gauss_pen_frag
# fpc_lik_frag
# mlfpc_lik_frag
# mfpc_lik_frag
# logistic_fpc_lik_frag

# Natural Parameter to Common Parameter Updates and
# Expectatoin of Sufficient Statistic Computations:
# igw_q
# gauss_q

# Entropy Functions:
# entropy_igw
# entropy_gauss
# entropy_two_lev_gauss
# entropy_two_lev_fpc_prod

# Cross-Entropy Functions:
# cross_entropy_igw_prior
# cross_entropy_iter_igw
# cross_entropy_gauss_prior
# cross_entropy_logistic_lik
# cross_entropy_two_lev_gauss_prior
# cross_entropy_fpc_gauss_pen
# cross_entropy_fpc_lik_frag
# cross_entropy_mfpc_lik_frag
# cross_entropy_mlfpc_lik_frag
# cross_entropy_logistic_fpc_lik_frag

# FPCA functions:
# fpc_orthgn
# mfpc_orthgn
# mlfpc_orthgn
# fpc_orthgn_two_lev
# fpca_mc
# logistic_fpca_mc

# Multilevel functions:
# solve_two_lev_sparse_mat
# two_level_nat_to_comm_parms
# two_level_fpc_orthogonalization
# two_level_fpc_rotation

# Miscellaneous Functions:
# is_int
# tr
# cprod

##########################################
#
#  FRAGMENT  UPDATES
#
##########################################

igw_prior_frag <- function(h_params=list(G_Theta, xi_Theta, Lambda_Theta)) {

  if(!is.list(h_params)) {

    stop("h_params must be a list")
  }

  if(length(h_params)!=3) {

    stop(
      "h_params must have length 3"
     )
  }

  G_Theta <- h_params[[1]]
  xi_Theta <- h_params[[2]]
  Lambda_Theta <- h_params[[3]]

  if(xi_Theta<0) {

    stop("xi_Theta must be non-negative")
  }

  if(is.vector(Lambda_Theta)) {

    if(length(Lambda_Theta)!=1) {

      stop("Lambda_Theta must be square")
    }

    if(Lambda_Theta<=0) {

      stop("Lambda_Theta must be positive-definite")
    }

    d <- 1

    D_d <- 1
  }

  if(is.matrix(Lambda_Theta)) {

    if(nrow(Lambda_Theta)!=ncol(Lambda_Theta)) {

      stop("Lambda_Theta must be square")
    }

    if(any(Lambda_Theta!=t(Lambda_Theta))) {

      stop("Lambda_Theta must be symmetric")
    }

    if(any(eigen(Lambda_Theta)$values<0)) {

      stop("Lambda_Theta must be positive-definite")
    }

    d <- nrow(Lambda_Theta)

    D_d <- duplication.matrix(d)
  }

  eta_1 <- -0.5*(xi_Theta + 2)
  eta_2 <- -0.5*as.vector(crossprod(D_d, as.vector(Lambda_Theta)))
  eta_vec <- c(eta_1, eta_2)

  ans <- list(G_Theta, eta_vec)
  names(ans) <- c("G_pTheta_Theta", "eta_pTheta_Theta")
  return(ans)
}

iter_igw_frag <- function(eta_in, G_mess, xi, G_hyper) {

  # order of eta_in:
  # 1. Sigma -> p(Sigma|A)
  # 2. p(Sigma|A) -> Sigma
  # 3. A -> p(Sigma|A)
  # 4. p(Sigma|A) -> A

  # order of eta_out:
  # 1. p(Sigma|A) -> Sigma
  # 2. p(Sigma|A) -> A

  # order of G_out:
  # 1. p(Sigma|A) -> Sigma
  # 2. p(Sigma|A) -> A

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in)!=4) {

    stop("eta_in must have length 4")
  }

  if(!(is.character(G_mess) & is.character(G_hyper))) {

    stop("G_mess and G_hyper must be a characters")
  } else {

    if((G_mess != "full") & (G_mess != "diag")) {

      stop("G_mess can only be 'full' or 'diag'")
    }

    if((G_hyper != "full") & (G_hyper != "diag")) {

      stop("G_hyper can only be 'full' or 'diag'")
    }
  }

  if(xi<0) {

    stop("xi must be non-negative")
  }

  l <- length(eta_in[[1]][-1])
  d <- (sqrt(1 + 8*l) - 1)/2

  if(d==1) {

    D_d <- 1
  } else {

    D_d <- duplication.matrix(d)
  }

  D_d_plus <- solve(crossprod(D_d), t(D_d))

  G_Sigma <- G_hyper
  G_A <- G_mess

  eta_Sigma <- eta_in[[1]] + eta_in[[2]]
  eta_A <- eta_in[[3]] + eta_in[[4]]

  if(G_A=="full") {

    omega_1 <- (d + 1)/2
  } else {

    omega_1 <- 1
  }

  eta_2_mat <- matrix(as.vector(crossprod(D_d_plus, eta_A[-1])), d, d)
  M_q_A_inv <- (eta_A[1] + omega_1)*solve(eta_2_mat)

  if(G_Sigma=="diag") {

    M_q_A_inv <- diag(diag(M_q_A_inv), nrow=d)
    omega_2 <- 1
  } else {

    omega_2 <- (d + 1)/2
  }

  eta_2_mat <- matrix(as.vector(crossprod(D_d_plus, eta_Sigma[-1])), d, d)
  M_q_Sigma_inv <- (eta_Sigma[1] + omega_2)*solve(eta_2_mat)

  if(G_A=="diag") {

    M_q_Sigma_inv <- diag(diag(M_q_Sigma_inv), nrow=d)
  }

  eta_out <- vector("list", length=2)
  names(eta_out) <- c("p(Sigma|A)->Sigma", "p(Sigma|A)->A")

  eta_1 <- -0.5*(xi + 2)
  eta_2 <- -0.5*as.vector(crossprod(D_d, as.vector(M_q_A_inv)))
  eta_out$"p(Sigma|A)->Sigma" <- c(eta_1, eta_2)

  eta_1 <- -0.5*(xi + 2 - 2*omega_2)
  eta_2 <- -0.5*as.vector(crossprod(D_d, as.vector(M_q_Sigma_inv)))
  eta_out$"p(Sigma|A)->A" <- c(eta_1, eta_2)

  G_out <- vector("list", length=2)
  names(G_out) <- c("p(Sigma|A)->Sigma", "p(Sigma|A)->A")
  G_out[[1]] <- G_Sigma
  G_out[[2]] <- G_A

  ans <- list(G_out, eta_out)
  names(ans) <- c("G", "eta")
  return(ans)
}

gauss_prior_frag <- function(mu, Sigma, use_vech=TRUE) {

  if(nrow(Sigma)!=ncol(Sigma)) {

    stop("Sigma must be square")
  }

  if(any(Sigma!=t(Sigma))) {

    stop("Sigma must be symmetric")
  }

  if(any(eigen(Sigma)$values<0)) {

    stop("Sigma must be positive-definite")
  }

  if(length(mu)!=nrow(Sigma)) {

    stop("length of mu is not equal to the number of rows of Sigma")
  }

  if(!is.logical(use_vech)) {

    stop("use_vech must be of type logical")
  }

  d <- length(mu)
  if(d==1) {

    use_vech <- FALSE
  }

  eta_1 <- solve(Sigma, mu)
  vec_Sigma_inv <- as.vector(solve(Sigma))

  if(use_vech) {

    D_d <- duplication.matrix(d)
    eta_2 <- -0.5*cprod(D_d, vec_Sigma_inv)
  } else{

    eta_2 <- -0.5*vec_Sigma_inv
  }

  eta_out <- c(eta_1, eta_2)
  return(eta_out)
}

gauss_lik_frag <- function(eta_in, G_in, y, A, use_vech = TRUE) {

  # order of eta_in:
  # 1. theta -> p(y|theta,sigsq)
  # 2. p(y|theta,sigsq) -> theta
  # 3. sigsq -> p(y|theta,sigsq)
  # 4. p(y|theta,sigsq) -> sigsq

  # order of G_in:
  # 1. sigsq -> p(y|theta,sigsq)
  # 2. p(y|theta,sigsq) -> sigsq

  # order of eta_out:
  # 1. p(y|theta,sigsq) -> theta
  # 2. p(y|theta,sigsq) -> sigsq

  # order of G_out:
  # 1. p(y|theta,sigsq) -> sigsq

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in) != 4) {

    stop("eta_in must have length 4")
  }

  if(!is.list(G_in)) {

    stop("G_in must be a list")
  }

  if(length(G_in) != 2) {

    stop("G_in must have length 2")
  }

  if(!all(Reduce("c", lapply(G_in, is.character)))) {

    stop("all of the elements of G_in must be characters")
  } else {

    if(!all(Reduce("c", lapply(G_in, function(x) x=="diag" || x=="full")))) {

      stop("The input graphs must be 'diag' or 'full'.")
    }
  }

  if(G_in[[1]]==G_in[[2]]) {

    G_sigsq <- G_in[[1]]
  } else {

    stop("The graph messages involving sigsq are not identical")
  }

  n <- length(y)

  eta_in_theta <- list(eta_in[[1]], eta_in[[2]])
  q_theta <- gauss_q(eta_in_theta, use_vech = use_vech)
  E_q_theta <- q_theta[[1]]
  Cov_q_theta <- q_theta[[2]]

  d <- length(E_q_theta)

  if(use_vech) {

    D_theta <- duplication.matrix(d)
  }

  eta_in_sigsq <- list(eta_in[[3]], eta_in[[4]])
  G_in_sigsq <- c(G_in[[1]], G_in[[2]])
  q_sigsq <- igw_q(eta_in_sigsq, G_in_sigsq)
  E_q_recip_sigsq <- q_sigsq[[3]]

  eta_theta_1 <- E_q_recip_sigsq*cprod(A, y)

  if(use_vech) {

    eta_theta_2 <- -1/2*E_q_recip_sigsq*cprod(D_theta, as.vector(crossprod(A)))
  } else {

    eta_theta_2 <- -1/2*E_q_recip_sigsq*as.vector(crossprod(A))
  }

  eta_theta <- c(eta_theta_1, eta_theta_2)

  mean_vec <- y - as.vector(A %*% E_q_theta)
  Cov_mat <- tcrossprod(A %*% Cov_q_theta, A)

  eta_sigsq_1 <- -n/2
  eta_sigsq_2 <- -1/2*E_cprod(mean_vec, Cov_mat, mean_vec)
  eta_sigsq <- c(eta_sigsq_1, eta_sigsq_2)

  G_sigsq <- "full"

  eta_out <- list(eta_theta, eta_sigsq)
  names(eta_out) <- c("p(y|theta,sigsq)->theta", "p(y|theta,sigsq)->sigsq")

  G_out <- list(G_sigsq)
  names(G_out) <- c("p(y|theta,sigsq)->sigsq")

  ans <- list(eta_out, G_out)
  names(ans) <- c("eta", "G")

  return(ans)
}

logistic_lik_frag <- function(eta_in, y, A) {

  # order of eta_in:
  # 1. theta -> p(y|theta)
  # 2. p(y|theta) -> theta

  # order of eta_out:
  # 1. p(y|theta) -> theta

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in) != 2) {

    stop("eta_in must have length 2")
  }

  # Set up p and s Monahan-Stefanski vectors:

  p <- c(
    0.003246343272134, 0.051517477033972, 0.195077912673858, 0.315569823632818,
    0.274149576158423, 0.131076880695470, 0.027912418727972, 0.001449567805354
  )

  s <- c(
    1.365340806296348, 1.059523971016916, 0.830791313765644, 0.650732166639391,
    0.508135425366489, 0.396313345166341, 0.308904252267995, 0.238212616409306
  )

  n <- length(y)

  eta_in_theta <- list(eta_in[[1]], eta_in[[2]])
  q_theta <- gauss_q(eta_in_theta, use_vech = TRUE)
  E_q_theta <- q_theta[[1]]
  Cov_q_theta <- q_theta[[2]]

  d <- length(E_q_theta)
  D_d <- duplication.matrix(d)

  mu <- as.vector(A %*% E_q_theta)
  sigsq <- diag(tcrossprod(A %*% Cov_q_theta, A))

  Omega <- sqrt(matrix(1, n, 8) + tcrossprod(sigsq, s^2))
  omega_1 <- as.vector(pnorm(tcrossprod(mu, s)/Omega) %*% p)
  omega_2 <- as.vector((dnorm(tcrossprod(mu, s)/Omega)/Omega) %*% (p*s))

  eta_1 <- cprod(A, y - omega_1 + omega_2*mu)
  eta_2 <- -0.5*cprod(D_d, as.vector(crossprod(A, diag(omega_2) %*% A)))
  eta_out <- list(c(eta_1, eta_2))
  names(eta_out) <- "p(y|theta)->theta"

  ans <- list(eta_out)
  names(ans) <- c("eta")

  return(ans)
}

gauss_pen_frag <- function(eta_in, G_in, sigsq_beta, use_vech = TRUE) {

  # order of eta_in:
  # 1. nu -> p(nu|sigsq_nu)
  # 2. p(nu|sigsq_nu) -> nu
  # 3. sigsq_nu -> p(nu|sigsq_nu)
  # 4. p(nu|sigsq_nu) -> sigsq_nu

  # order of G_in:
  # 1. sigsq_nu -> p(nu|sigsq_nu)
  # 2. p(nu|sigsq_nu) -> sigsq_nu

  # order of eta_out:
  # 1. p(nu|sigsq_nu) -> nu
  # 2. p(nu|sigsq_nu) -> sigsq_nu

  # order of G_out:
  # 1. p(nu|sigsq_nu) -> sigsq_nu

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in)!=4) {

    stop("eta_in must have length 4")
  }

  if(!is.list(G_in)) {

    stop("G_in must be a list")
  }

  if(length(G_in)!=2) {

    stop("G_in must have length 2")
  }

  if(!all(Reduce("c", lapply(G_in, is.character)))) {

    stop("all of the elements of G_in must be characters")
  } else {

    if(!all(Reduce("c", lapply(G_in, function(x) x=="diag" || x=="full")))) {

      stop("The input graphs must be 'diag' or 'full'.")
    }
  }

  if(G_in[[1]]==G_in[[2]]) {

    G_sigsq_nu <- G_in[[1]]
  } else {

    stop("The graph messages involving sigsq_nu are not identical")
  }

  if(sigsq_beta <= 0) {

    stop("sigsq_beta must be positive")
  }

  eta_in_nu <- list(eta_in[[1]], eta_in[[2]])
  q_nu <- gauss_q(eta_in_nu, use_vech = use_vech)
  mu_q_nu <- q_nu[[1]]
  Sigma_q_nu <- q_nu[[2]]

  d <- length(mu_q_nu)
  K <- d - 2

  if(use_vech) {

    D_nu <- duplication.matrix(K+2)
  }

  eta_in_sigsq_nu <- list(eta_in[[3]], eta_in[[4]])
  G_in_sigsq_nu <- c(G_in[[1]], G_in[[2]])
  q_sigsq_nu <- igw_q(eta_in_sigsq_nu, G_in_sigsq_nu)
  mu_q_recip_sigsq_nu <- q_sigsq_nu[[3]]
  E_q_inv_Sigma_nu <- adiag(1/sigsq_beta*diag(2), mu_q_recip_sigsq_nu*diag(K))

  eta_1 <- rep(0, K+2)

  if(use_vech) {

    eta_2 <- -0.5*cprod(D_nu, as.vector(E_q_inv_Sigma_nu))
  } else {

    eta_2 <- -0.5*as.vector(E_q_inv_Sigma_nu)
  }

  eta_nu <- c(eta_1, eta_2)

  eta_1 <- -K/2
  eta_2 <- -0.5*(tr(Sigma_q_nu[-c(1, 2), -c(1, 2)]) + cprod(mu_q_nu[-c(1, 2)]))
  eta_sigsq_nu <- c(eta_1, eta_2)

  G_sigsq_nu <- "full"

  eta_out <- list(eta_nu, eta_sigsq_nu)
  names(eta_out) <- c("p(nu|sigsq_nu) -> nu", "p(nu|sigsq_nu) -> sigsq_nu")

  G_out <- list(G_sigsq_nu)
  names(G_out) <- c("p(nu|sigsq_nu) -> sigsq_nu")

  ans <- list(eta_out, G_out)
  names(ans) <- c("eta", "G")

  return(ans)
}

mult_gauss_pen_frag <- function(eta_in, G_in, L, sigsq_beta) {

  # order of eta_in:
  # 1. nu -> p(nu|Sigma_nu)
  # 2. p(nu|Sigma_nu) -> nu
  # 3. sigsq_p -> p(nu|Sigma_nu)
  # 4. p(nu|Sigma_nu) -> sigsq_p

  # order of G_in:
  # 1. sigsq_p -> p(nu|Sigma_nu)
  # 2. p(nu|Sigma_nu) -> sigsq_p

  # order of eta_out:
  # 1. p(nu|Sigma_nu) -> nu
  # 2. p(nu|Sigma_nu) -> sigsq_p

  # order of G_out:
  # 1. p(nu|Sigma_nu) -> sigsq_p

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in)!=4) {

    stop("eta_in must have length 4")
  }

  if(!is.list(G_in)) {

    stop("G_in must be a list")
  }

  if(length(G_in)!=2) {

    stop("G_in must have length 2")
  }

  if(!all(Reduce("c", lapply(G_in, is.character)))) {

    stop("all of the elements of G_in must be characters")
  } else {

    if(!all(Reduce("c", lapply(G_in, function(x) x=="diag" || x=="full")))) {

      stop("The input graphs must be 'diag' or 'full'.")
    }
  }

  if(all(G_in[[1]]==G_in[[2]])) {

    G_sigsq_p <- G_in[[1]]
  } else {

    stop("The graph messages involving sigsq_m are not identical")
  }

  if(sigsq_beta <= 0) {

    stop("sigsq_beta must be positive")
  }

  if(!is_int(L)) {

    stop("L must be a whole number.")
  } else {

    if(L<1) {

      stop("There must be at least one basis function.")
    }
  }

  l_eta_nu <- length(eta_in[[1]])
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/L - 2

  eta_in_nu <- list(eta_in[[1]], eta_in[[2]])
  q_nu <- gauss_q(eta_in_nu, use_vech=FALSE)
  mu_q_nu <- q_nu[[1]]
  Sigma_q_nu <- q_nu[[2]]

  inds <- seq(1, ((K+2)*L))
  groups <- ceiling(inds/(K+2))
  ind_groups <- split(inds, groups)
  ind_groups <- lapply(ind_groups, function(x) x[-c(1,2)])
  Sigma_q_u_p <- vector("list", length=L)
  mu_q_u_p <- vector("list", length=L)
  for(l in 1:L) {

    inds_l <- ind_groups[[l]]
    mu_q_u_p[[l]] <- mu_q_nu[inds_l]
    Sigma_q_u_p[[l]] <- Sigma_q_nu[inds_l, inds_l]
  }

  M_q_inv_Sigma_p <- vector("list", length = L)
  for(l in 1:L) {

    eta_in_sigsq_p <- list(eta_in[[3]][,l], eta_in[[4]][,l])
    G_in_sigsq_p <- list(G_in[[1]][l], G_in[[2]][l])
    q_sigsq_p <- igw_q(eta_in_sigsq_p, G_in_sigsq_p)
    mu_q_recip_sigsq_p <- q_sigsq_p$E_X_inv
    M_q_inv_Sigma_p[[l]] <- adiag(1/sigsq_beta*diag(2), mu_q_recip_sigsq_p*diag(K))
  }
  M_q_inv_Sigma_p <- Reduce(adiag, M_q_inv_Sigma_p)

  eta_1 <- rep(0, L*(K + 2))
  eta_2 <- -0.5*as.vector(M_q_inv_Sigma_p)
  eta_nu <- c(eta_1, eta_2)

  eta_1 <- -K/2
  eta_sigsq_p <- matrix(NA, 2, L)
  for(l in 1:L) {

    tr_term <- tr(Sigma_q_u_p[[l]])
    cprod_term <- cprod(mu_q_u_p[[l]])
    eta_2 <- -0.5*(tr_term + cprod_term)
    eta_sigsq_p[,l] <- c(eta_1, eta_2)
  }

  eta_out <- list(eta_nu, eta_sigsq_p)
  names(eta_out) <- c("p(nu_p|sigsq_p)->nu_p", "p(nu_p|sigsq_p)->sigsq_p")

  G_out <- list(G_sigsq_p)
  names(G_out) <- c("p(nu_p|sigsq_p)->sigsq_p")

  ans <- list(eta_out, G_out)
  names(ans) <- c("eta", "G")

  return(ans)
}

fpc_gauss_pen_frag <- function(eta_in, G_in, L, mu_beta, Sigma_beta) {

  # order of eta_in:
  # 1. nu -> p(nu|Sigma_nu)
  # 2. p(nu|Sigma_nu) -> nu
  # 3. sigsq_m -> p(nu|Sigma_nu)
  # 4. p(nu|Sigma_nu) -> sigsq_m
  # 5. sigsq_p -> p(nu|Sigma_nu)
  # 6. p(nu|Sigma_nu) -> sigsq_p

  # order of G_in:
  # 1. sigsq_m -> p(nu|Sigma_nu)
  # 2. p(nu|Sigma_nu) -> sigsq_m
  # 3. sigsq_p -> p(nu|Sigma_nu)
  # 4. p(nu|Sigma_nu) -> sigsq_p

  # order of eta_out:
  # 1. p(nu|Sigma_nu) -> nu
  # 2. p(nu|Sigma_nu) -> sigsq_m
  # 3. p(nu|Sigma_nu) -> sigsq_p

  # order of G_out:
  # 1. p(nu|Sigma_nu) -> sigsq_m
  # 2. p(nu|Sigma_nu) -> sigsq_p

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in)!=6) {

    stop("eta_in must have length 6")
  }

  if(!is.list(G_in)) {

    stop("G_in must be a list")
  }

  if(length(G_in)!=4) {

    stop("G_in must have length 4")
  }

  if(!all(Reduce("c", lapply(G_in, is.character)))) {

    stop("all of the elements of G_in must be characters")
  } else {

    if(!all(Reduce("c", lapply(G_in, function(x) x=="diag" | x=="full")))) {

      stop("The input graphs must be 'diag' or 'full'.")
    }
  }

  if(G_in[[1]]==G_in[[2]]) {

    G_sigsq_m <- G_in[[1]]
  } else {

    stop("The graph messages involving sigsq_m are not identical")
  }

  if(all(G_in[[3]]==G_in[[4]])) {

    G_sigsq_p <- G_in[[3]]
  } else {

    stop("The graph messages involving sigsq_p are not identical")
  }

  if(nrow(Sigma_beta)!=ncol(Sigma_beta)) {

    stop("Sigma_beta must be square")
  }

  if(any(Sigma_beta!=t(Sigma_beta))) {

    stop("Sigma_beta must be symmetric")
  }

  if(any(eigen(Sigma_beta)$values<0)) {

    stop("Sigma_beta must be positive-definite")
  }

  if(length(mu_beta)!=nrow(Sigma_beta)) {

    stop("length of mu_beta is not equal to the number of rows of Sigma_beta")
  }

  if(!is_int(L)) {

    stop("L must be a whole number.")
  } else {

    if(L<1) {

      stop("There must be at least one basis function.")
    }
  }

  eta_q_nu <- eta_in[[1]] + eta_in[[2]]
  eta_q_sigsq_m <- eta_in[[3]] + eta_in[[4]]
  eta_q_sigsq_p <- eta_in[[5]] + eta_in[[6]]

  l_eta_nu <- length(eta_q_nu)
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/(L+1) - 2

  eta_q_nu_1 <- eta_q_nu[1:d]
  eta_q_nu_2 <- eta_q_nu[-c(1:d)]

  Sigma_q_nu <- -0.5*solve(matrix(eta_q_nu_2, d, d))
  mu_q_nu <- as.vector(Sigma_q_nu%*%eta_q_nu_1)

  mu_inds <- 1:(K+2)
  mu_q_nu_mu <- mu_q_nu[mu_inds]
  Sigma_q_nu_mu <- Sigma_q_nu[mu_inds, mu_inds]

  mu_q_u_mu <- mu_q_nu_mu[-c(1,2)]
  Sigma_q_u_mu <- Sigma_q_nu_mu[-c(1,2), -c(1,2)]

  mu_q_nu_psi <- mu_q_nu[-mu_inds]
  Sigma_q_nu_psi <- Sigma_q_nu[-mu_inds, -mu_inds]

  psi_inds <- seq(1, d-(K+2))
  groups <- ceiling(psi_inds/(K+2))
  psi_groups <- split(psi_inds, groups)
  u_psi_groups <- lapply(psi_groups, function(x) x[-c(1,2)])
  Sigma_q_u_psi <- vector("list", length=L)
  mu_q_u_psi <- vector("list", length=L)
  for(l in 1:L) {

    inds_l <- u_psi_groups[[l]]
    mu_q_u_psi[[l]] <- mu_q_nu_psi[inds_l]
    Sigma_q_u_psi[[l]] <- Sigma_q_nu_psi[inds_l, inds_l]
  }

  mu_q_recip_sigsq_m <- (eta_q_sigsq_m[1] + 1)/eta_q_sigsq_m[2]
  M_q_inv_Sigma_m <- blkdiag(solve(Sigma_beta), mu_q_recip_sigsq_m*diag(K))
  mu_m <- c(mu_beta, rep(0, K))

  mu_q_recip_sigsq_p <- (eta_q_sigsq_p[1,] + 1)/eta_q_sigsq_p[2,]
  M_q_inv_Sigma_p <- vector("list", length=L)
  mu_p <- vector("list", length=L)
  for(l in 1:L) {

    M_q_inv_Sigma_p[[l]] <- blkdiag(
      solve(Sigma_beta),
      mu_q_recip_sigsq_p[l]*diag(K)
    )

    mu_p[[l]] <- c(mu_beta, rep(0, K))
  }

  M_q_inv_Sigma_nu <- blkdiag(M_q_inv_Sigma_m, Reduce(blkdiag, M_q_inv_Sigma_p))
  mu_nu <- c(mu_m, Reduce("c", mu_p))

  eta_1 <- as.vector(M_q_inv_Sigma_nu%*%mu_nu)
  eta_2 <- -0.5*as.vector(M_q_inv_Sigma_nu)
  eta_nu <- c(eta_1, eta_2)

  eta_1 <- -K/2

  eta_2 <- -0.5*(sum(diag(Sigma_q_u_mu)) + as.vector(crossprod(mu_q_u_mu)))
  eta_sigsq_m <- c(eta_1, eta_2)

  eta_sigsq_p <- matrix(NA, 2, L)
  for(l in 1:L) {

    tr_term <- sum(diag(Sigma_q_u_psi[[l]]))
    cprod_term <- as.vector(crossprod(mu_q_u_psi[[l]]))
    eta_2 <- -0.5*(tr_term + cprod_term)
    eta_sigsq_p[,l] <- c(eta_1, eta_2)
  }

  eta_out <- list(eta_nu, eta_sigsq_m, eta_sigsq_p)
  names(eta_out) <- c(
    "p(nu|Sigma_nu)->nu",
    "p(nu|Sigma_nu)->sigsq_mu",
    "p(nu|Sigma_nu)->sigsq_psi"
  )

  G_out <- list(G_sigsq_m, G_sigsq_p)
  names(G_out) <- c("p(nu|Sigma_nu)->sigsq_mu", "p(nu|Sigma_nu)->sigsq_psi")

  ans <- list(eta_out, G_out)
  names(ans) <- c("eta", "G")

  return(ans)
}

fpc_lik_frag <- function(eta_in, G_in, C, Y, T_vec, L) {

  # order of eta_in:
  # 1. nu -> p(Y|nu,zeta,sigsq_eps)
  # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # 3. zeta -> p(Y|nu,zeta,sigsq_eps)
  # 4. p(Y|nu,zeta,sigsq_eps) -> zeta
  # 5. sigsq_eps -> p(Y|nu,zeta,sigsq_eps)
  # 6. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of G_in:
  # 1. sigsq_eps -> p(Y|nu,zeta,sigsq_eps)
  # 2. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of eta_out:
  # 1. p(Y|nu,zeta,sigsq_eps) -> nu
  # 2. p(Y|nu,zeta,sigsq_eps) -> zeta
  # 3. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of G_out:
  # 1. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in)!=6) {

    stop("eta_in must have length 6")
  }

  if(!is.list(G_in)) {

    stop("G_in must be a list")
  }

  if(length(G_in)!=2) {

    stop("G_in must have length 2")
  }

  if(!all(Reduce("c", lapply(G_in, is.character)))) {

    stop("all of the elements of G_in must be characters")
  } else {

    if(!all(Reduce("c", lapply(G_in, function(x) x=="diag" || x=="full")))) {

      stop("The input graphs must be 'diag' or 'full'.")
    }
  }

  if(G_in[[1]]==G_in[[2]]) {

    G_sigsq_eps <- G_in[[1]]
  } else {

    stop("The graph messages involving sigsq_eps are not identical")
  }

  if(!is_int(L)) {

    stop("L must be a whole number.")
  } else {

    if(L<1) {

      stop("There must be at least one basis function.")
    }
  }

  N <- ncol(eta_in[[3]])
  l_eta_nu <- length(eta_in[[1]])
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/(L+1) - 2

  eta_in_nu <- list(eta_in[[1]], eta_in[[2]])
  q_nu <- gauss_q(eta_in_nu, use_vech=FALSE)
  mu_q_nu <- q_nu[[1]]
  Sigma_q_nu <- q_nu[[2]]

  eta_in_sigsq_eps <- list(eta_in[[5]], eta_in[[6]])
  G_in_sigsq_eps <- c(G_in[[1]], G_in[[2]])
  q_sigsq_eps <- igw_q(eta_in_sigsq_eps, G_in_sigsq_eps)

  mu_q_recip_sigsq_eps <- q_sigsq_eps[[3]]

  eta_1_sum <- 0
  eta_2_sum <- 0
  mu_q_zeta <- vector("list", length=N)
  Sigma_q_zeta <- vector("list", length=N)
  for(i in 1:N) {

    eta_in_zeta <- list(
      eta_in[[3]][,i],
      eta_in[[4]][,i]
    )
    q_zeta <- gauss_q(eta_in_zeta, use_vech=TRUE)
    mu_q_zeta[[i]] <- q_zeta[[1]]
    Sigma_q_zeta[[i]] <- q_zeta[[2]]

    mu_q_zeta_tilde <- c(1, mu_q_zeta[[i]])
    Sigma_q_zeta_tilde <- blkdiag(matrix(0), Sigma_q_zeta[[i]])
    M_q_zeta_zeta_T_tilde <- Sigma_q_zeta_tilde + tcrossprod(mu_q_zeta_tilde)

    sum_val <- cprod(kronecker(t(mu_q_zeta_tilde), C[[i]]), Y[[i]])
    eta_1_sum <- eta_1_sum + sum_val

    sum_val <- as.vector(kronecker(M_q_zeta_zeta_T_tilde, crossprod(C[[i]])))
    eta_2_sum <- eta_2_sum + sum_val
  }
  eta_1 <- mu_q_recip_sigsq_eps*eta_1_sum
  eta_2 <- -0.5* mu_q_recip_sigsq_eps*eta_2_sum
  eta_nu <- c(eta_1, eta_2)

  mu_inds <- 1:(K+2)
  mu_q_nu_mu <- mu_q_nu[mu_inds]
  Sigma_q_nu_mu <- Sigma_q_nu[mu_inds, mu_inds]

  psi_inds <- (1:d)[-mu_inds]
  psi_groups <- ceiling(psi_inds/(K+2))-1
  psi_inds <- split(psi_inds, psi_groups)

  M_q_H <- vector("list", length=N)
  mu_q_h <- vector("list", length=N)
  M_q_H_tilde <- vector("list", length=N)
  mu_q_h_mu <- rep(NA, N)
  for(i in 1:N) {

    tr_term <- tr(Sigma_q_nu_mu%*%crossprod(C[[i]]))
    cprod_term <- cprod(mu_q_nu_mu, crossprod(C[[i]])%*%mu_q_nu_mu)
    mu_q_h_mu[i] <- tr_term + cprod_term

    M_q_H[[i]] <- matrix(NA, L, L)
    mu_q_h[[i]] <- rep(NA, L)

    for(l_i in 1:L) {

      i_inds <- psi_inds[[l_i]]
      mu_q_nu_psi_i <- mu_q_nu[i_inds]

      Cov_q_nu_psi_i <- Sigma_q_nu[mu_inds, i_inds]

      tr_term <- tr(Cov_q_nu_psi_i%*%crossprod(C[[i]]))
      cprod_term <- cprod(mu_q_nu_psi_i, crossprod(C[[i]])%*%mu_q_nu_mu)
      mu_q_h[[i]][l_i] <- tr_term + cprod_term

      for(l_j in 1:L) {

        j_inds <- psi_inds[[l_j]]
        mu_q_nu_psi_j <- mu_q_nu[j_inds]

        Cov_q_nu_psi_ji <- Sigma_q_nu[j_inds, i_inds]

        tr_term <- tr(Cov_q_nu_psi_ji%*%crossprod(C[[i]]))
        cprod_term <- cprod(mu_q_nu_psi_i, crossprod(C[[i]])%*%mu_q_nu_psi_j)
        M_q_H[[i]][l_i, l_j] <- tr_term + cprod_term
      }
    }

    top_block <- t(c(mu_q_h_mu[i], mu_q_h[[i]]))
    bottom_block <- cbind(mu_q_h[[i]], M_q_H[[i]])
    M_q_H_tilde[[i]] <- rbind(top_block, bottom_block)
  }

  M_q_V_psi <- matrix(NA, K+2, L)
  for(l in 1:L) {

    psi_l_inds <- psi_inds[[l]]
    M_q_V_psi[,l] <- mu_q_nu[psi_l_inds]
  }

  M_q_V <- cbind(mu_q_nu_mu, M_q_V_psi)

  D_L <- duplication.matrix(L)

  eta_1 <- matrix(NA, L, N)
  eta_2 <- matrix(NA, 0.5*L*(L+1), N)
  for(i in 1:N) {

    M_q_Psi <- C[[i]]%*%M_q_V_psi
    eta_1[,i] <- mu_q_recip_sigsq_eps*(cprod(M_q_Psi, Y[[i]]) - mu_q_h[[i]])
    eta_2[,i] <- -0.5*mu_q_recip_sigsq_eps*cprod(D_L, as.vector(M_q_H[[i]]))
  }
  eta_zeta <- rbind(eta_1, eta_2)

  E_q_resid <- 0
  for(i in 1:N) {

    mu_q_zeta_tilde <- c(1, mu_q_zeta[[i]])
    Sigma_q_zeta_tilde <- blkdiag(matrix(0), Sigma_q_zeta[[i]])
    M_q_zeta_zeta_T_tilde <- Sigma_q_zeta_tilde + tcrossprod(mu_q_zeta_tilde)

    E_q_Y_hat <- C[[i]]%*%M_q_V%*%mu_q_zeta_tilde

    term_1 <- cprod(Y[[i]])
    term_2 <- -2*cprod(E_q_Y_hat, Y[[i]])
    term_3 <- tr(M_q_zeta_zeta_T_tilde%*%M_q_H_tilde[[i]])
    sum_val <- term_1 + term_2 + term_3

    E_q_resid <- E_q_resid + sum_val
  }

  eta_1 <- -0.5*sum(T_vec)
  eta_2 <- -0.5*E_q_resid
  eta_sigsq_eps <- c(eta_1, eta_2)

  eta_out <- list(eta_nu, eta_zeta, eta_sigsq_eps)
  names(eta_out) <- c(
    "p(Y|nu,zeta,sigsq_eps)->nu",
    "p(Y|nu,zeta,sigsq_eps)->zeta",
    "p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
  )

  G_out <- G_sigsq_eps
  names(G_out) <- "p(Y|nu,zeta,sigsq_eps)->sigsq_eps"

  ans <- list(eta_out, G_out)
  names(ans) <- c("eta", "G")

  return(ans)
}

mlfpc_lik_frag <- function(eta_in, G_in, C, Y, T_vec, L_1, L_2) {

  # order of eta_in:
  # 1. nu -> p(Y|nu,zeta,sigsq_eps)
  # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # 3. zeta -> p(Y|nu,zeta,sigsq_eps)
  # 4. p(Y|nu,zeta,sigsq_eps) -> zeta
  # 5. sigsq_eps -> p(Y|nu,zeta,sigsq_eps)
  # 6. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of G_in:
  # 1. sigsq_eps -> p(Y|nu,zeta,sigsq_eps)
  # 2. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of eta_out:
  # 1. p(Y|nu,zeta,sigsq_eps) -> nu
  # 2. p(Y|nu,zeta,sigsq_eps) -> zeta
  # 3. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of G_out:
  # 1. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in)!=6) {

    stop("eta_in must have length 6")
  }

  if(!is.list(G_in)) {

    stop("G_in must be a list")
  }

  if(length(G_in)!=2) {

    stop("G_in must have length 2")
  }

  if(!all(Reduce("c", lapply(G_in, is.character)))) {

    stop("all of the elements of G_in must be characters")
  } else {

    if(!all(Reduce("c", lapply(G_in, function(x) x=="diag" || x=="full")))) {

      stop("The input graphs must be 'diag' or 'full'.")
    }
  }

  if(G_in[[1]]==G_in[[2]]) {

    G_sigsq_eps <- G_in[[1]]
  } else {

    stop("The graph messages involving sigsq_eps are not identical")
  }

  if(!is_int(L_1)) {

    stop("L_1 must be a whole number.")
  } else {

    if(L_1 < 1) {

      stop("There must be at least one basis function for the first level.")
    }
  }

  if(!is_int(L_2)) {

    stop("L_2 must be a whole number.")
  } else {

    if(L_2 < 1) {

      stop("There must be at least one basis function for the second level.")
    }
  }

  L <- L_1 + L_2
  N <- length(eta_in[[3]])
  l_eta_nu <- length(eta_in[[1]])
  l_eta_zeta <- sapply(eta_in[[3]], length)
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/(L+1) - 2
  l_zeta_1 <- L_1 + 0.5*L_1*(L_1 + 1)
  l_zeta_2 <- L_2 + 0.5*L_2*(L_2 + 1) + L_1*L_2
  M <- (l_eta_zeta - l_zeta_1)/l_zeta_2

  D_zeta_1 <- duplication.matrix(L_1)
  D_zeta_2 <- duplication.matrix(L_2)

  YTY <- vector("list", length = N)
  CTY <- vector("list", length = N)
  CTC <- vector("list", length = N)
  for(i in 1:N) {

    YTY[[i]] <- rep(NA, M[i])
    CTY[[i]] <- matrix(NA, K + 2, M[i])
    CTC[[i]] <- vector("list", length = M[i])
    for(j in 1:M[i]) {

      YTY[[i]][j] <- cprod(Y[[i]][[j]])
      CTY[[i]][, j] <- cprod(C[[i]][[j]], Y[[i]][[j]])
      CTC[[i]][[j]] <- crossprod(C[[i]][[j]])
    }
  }

  eta_in_nu <- list(eta_in[[1]], eta_in[[2]])
  q_nu <- gauss_q(eta_in_nu, use_vech = FALSE)
  E_q_nu <- q_nu[[1]]
  Cov_q_nu <- q_nu[[2]]

  eta_in_sigsq_eps <- list(eta_in[[5]], eta_in[[6]])
  q_sigsq_eps <- igw_q(eta_in_sigsq_eps, G_in)
  E_q_recip_sigsq_eps <- q_sigsq_eps[[3]]

  eta_1_sum <- 0
  eta_2_sum <- 0
  E_q_zeta_1 <- vector("list", length = N)
  Cov_q_zeta_1 <- vector("list", length = N)
  E_q_zeta_2 <- vector("list", length = N)
  Cov_q_zeta_2 <- vector("list", length = N)
  Cov_q_zeta_12 <- vector("list", length = N)
  for(i in 1:N) {

    eta_in_zeta <- list(
      eta_in[[3]][[i]],
      eta_in[[4]][[i]]
    )
    q_zeta <- two_level_nat_to_comm_parms(L_1, L_2, M[i], eta_in_zeta)
    E_q_zeta_1[[i]] <- q_zeta[[1]]
    Cov_q_zeta_1[[i]] <- q_zeta[[2]]
    E_q_zeta_2[[i]] <- q_zeta[[3]]
    Cov_q_zeta_2[[i]] <- q_zeta[[4]]
    Cov_q_zeta_12[[i]] <- q_zeta[[5]]

    for(j in 1:M[i]) {

      E_q_zeta <- c(E_q_zeta_1[[i]], E_q_zeta_2[[i]][[j]])
      E_q_zeta_tilde <- c(1, E_q_zeta)

      top_block <- cbind(Cov_q_zeta_1[[i]], Cov_q_zeta_12[[i]][[j]])
      bottom_block <- cbind(t(Cov_q_zeta_12[[i]][[j]]), Cov_q_zeta_2[[i]][[j]])
      Cov_q_zeta <- rbind(top_block, bottom_block)
      Cov_q_zeta_tilde <- adiag(0, Cov_q_zeta)

      E_q_tcross_zeta_tilde <- Cov_q_zeta_tilde + tcrossprod(E_q_zeta_tilde)

      sum_val <- cprod(kronecker(t(E_q_zeta_tilde), C[[i]][[j]]), Y[[i]][[j]])
      eta_1_sum <- eta_1_sum + sum_val

      sum_val <- kronecker(E_q_tcross_zeta_tilde, CTC[[i]][[j]])
      eta_2_sum <- eta_2_sum + sum_val
    }
  }
  eta_1 <- E_q_recip_sigsq_eps*eta_1_sum
  eta_2 <- -0.5*E_q_recip_sigsq_eps*eta_2_sum
  eta_nu <- c(eta_1, eta_2)

  inds_mat <- matrix(1:d, K + 2, L + 1)
  mu_inds <- inds_mat[, 1]
  psi_1_inds <- inds_mat[, 1:L_1 + 1]
  psi_2_inds <- inds_mat[, 1:L_2 + L_1 + 1]

  E_q_V <- matrix(E_q_nu, K+2, L + 1)
  E_q_nu_mu <- E_q_V[, 1]
  E_q_V_psi_1 <- E_q_V[, 1:L_1 + 1]
  E_q_nu_psi_1 <- as.vector(E_q_V_psi_1)
  E_q_V_psi_2 <- E_q_V[, 1:L_2 + L_1 + 1]
  E_q_nu_psi_2 <- as.vector(E_q_V_psi_2)

  Cov_q_nu_mu <- Cov_q_nu[mu_inds, mu_inds]
  Cov_q_nu_psi_1 <- Cov_q_nu[as.vector(psi_1_inds), as.vector(psi_1_inds)]
  Cov_q_nu_psi_2 <- Cov_q_nu[as.vector(psi_2_inds), as.vector(psi_2_inds)]
  Cov_q_nu_mu_psi_1 <- Cov_q_nu[mu_inds, as.vector(psi_1_inds)]
  Cov_q_nu_mu_psi_2 <- Cov_q_nu[mu_inds, as.vector(psi_2_inds)]
  Cov_q_nu_psi_21 <- Cov_q_nu[as.vector(psi_2_inds), as.vector(psi_1_inds)]

  E_q_h_mu <- E_cprod(
    E_q_nu_mu, Cov_q_nu_mu, E_q_nu_mu,
    Reduce("+", lapply(CTC, Reduce, f = "+"))
  )

  eta_zeta_1 <- vector("list", length = N)
  eta_zeta_2 <- vector("list", length = N)
  E_q_h_mu_psi_1 <- vector("list", length = N)
  E_q_H_psi_11 <- vector("list", length = N)
  E_q_h_mu_psi_2 <- vector("list", length = N)
  E_q_H_psi_22 <- vector("list", length = N)
  E_q_H_psi_12 <- vector("list", length = N)
  for(i in 1:N) {

    CTY_i <- rowSums(CTY[[i]])
    CTC_i <- Reduce("+", CTC[[i]])

    E_q_h_mu_psi_1[[i]] <- E_h(L_1, E_q_nu_psi_1, Cov_q_nu_mu_psi_1, E_q_nu_mu, CTC_i)
    E_q_H_psi_11[[i]] <- E_H(L_1, L_1, E_q_nu_psi_1, Cov_q_nu_psi_1, E_q_nu_psi_1, CTC_i)
    eta_11 <- E_q_recip_sigsq_eps*(cprod(E_q_V_psi_1, CTY_i) - E_q_h_mu_psi_1[[i]])
    eta_12 <- -0.5*E_q_recip_sigsq_eps*cprod(D_zeta_1, as.vector(E_q_H_psi_11[[i]]))
    eta_zeta_1[[i]] <- c(eta_11, eta_12)

    E_q_h_mu_psi_2[[i]] <- vector("list", length = M[i])
    E_q_H_psi_22[[i]] <- vector("list", length = M[i])
    E_q_H_psi_12[[i]] <- vector("list", length = M[i])
    eta_21 <- matrix(NA, L_2, M[i])
    eta_22 <- matrix(NA, 0.5*L_2*(L_2 + 1), M[i])
    eta_23 <- matrix(NA, L_1*L_2, M[i])
    for(j in 1:M[i]) {

      E_q_h_mu_psi_2[[i]][[j]] <- E_h(
        L_2, E_q_nu_psi_2, Cov_q_nu_mu_psi_2,
        E_q_nu_mu, CTC[[i]][[j]]
      )

      E_q_H_psi_22[[i]][[j]] <- E_H(
        L_2, L_2, E_q_nu_psi_2, Cov_q_nu_psi_2,
        E_q_nu_psi_2, CTC[[i]][[j]]
      )

      E_q_H_psi_12[[i]][[j]] <- E_H(
        L_1, L_2, E_q_nu_psi_1, Cov_q_nu_psi_21,
        E_q_nu_psi_2, CTC[[i]][[j]]
      )

      eta_21[, j] <- E_q_recip_sigsq_eps*(
        cprod(E_q_V_psi_2, CTY[[i]][, j]) - E_q_h_mu_psi_2[[i]][[j]]
      )
      eta_22[, j] <- -0.5*E_q_recip_sigsq_eps*cprod(
        D_zeta_2,
        as.vector(E_q_H_psi_22[[i]][[j]])
      )
      eta_23[, j] <- -E_q_recip_sigsq_eps*as.vector(E_q_H_psi_12[[i]][[j]])
    }
    eta_zeta_2[[i]] <- as.vector(rbind(eta_21, eta_22, eta_23))
  }
  eta_zeta <- mapply(c, eta_zeta_1, eta_zeta_2, SIMPLIFY = FALSE)

  term_1 <- sum(sapply(YTY, sum))
  term_2 <- -2*cprod(E_q_nu_mu, Reduce("+", lapply(CTY, rowSums)))
  term_3 <- 0
  term_4 <- 0
  for(i in 1:N) {

    sum_val <- -2*cprod(E_q_V_psi_1 %*% E_q_zeta_1[[i]], rowSums(CTY[[i]]))
    term_3 <- term_3 + sum_val

    for(j in 1:M[i]) {

      sum_val <- -2*cprod(E_q_V_psi_2 %*% E_q_zeta_2[[i]][[j]], CTY[[i]][, j])
      term_4 <- term_4 + sum_val
    }
  }
  y_line <- term_1 + term_2 + term_3 + term_4

  term_1 <- E_q_h_mu
  term_2 <- 0
  term_3 <- 0
  for(i in 1:N) {

    sum_val <- 2*cprod(E_q_zeta_1[[i]], E_q_h_mu_psi_1[[i]])
    term_2 <- term_2 + sum_val

    for(j in 1:M[i]) {

      sum_val <- 2*cprod(E_q_zeta_2[[i]][[j]], E_q_h_mu_psi_2[[i]][[j]])
      term_3 <- term_3 + sum_val
    }
  }
  mu_line <- term_1 + term_2 + term_3

  term_1 <- 0
  term_2 <- 0
  for(i in 1:N) {

    E_vec <- E_q_zeta_1[[i]]
    Cov_mat <- Cov_q_zeta_1[[i]]
    A_mat <- E_q_H_psi_11[[i]]
    sum_val <- E_cprod(E_vec, Cov_mat, E_vec, A_mat)
    term_1 <- term_1 + sum_val

    for(j in 1:M[i]) {

      E_vec_1 <- E_q_zeta_1[[i]]
      Cov_mat <- t(Cov_q_zeta_12[[i]][[j]])
      E_vec_2 <- E_q_zeta_2[[i]][[j]]
      A_mat <- E_q_H_psi_12[[i]][[j]]
      sum_val <- 2*E_cprod(E_vec_1, Cov_mat, E_vec_2, A_mat)
      term_2 <- term_2 + sum_val
    }
  }
  zeta_1_line <- term_1 + term_2

  term_1 <- 0
  for(i in 1:N) {

    for(j in 1:M[i]) {

      E_vec <- E_q_zeta_2[[i]][[j]]
      Cov_mat <- Cov_q_zeta_2[[i]][[j]]
      A_mat <- E_q_H_psi_22[[i]][[j]]
      sum_val <- E_cprod(E_vec, Cov_mat, E_vec, A_mat)
      term_1 <- term_1 + sum_val
    }
  }
  zeta_2_line <- term_1

  E_q_norm <- y_line + mu_line + zeta_1_line + zeta_2_line

  eta_sigsq_eps_1 <- -0.5*sum(sapply(T_vec, sum))
  eta_sigsq_eps_2 <- -0.5*E_q_norm
  eta_sigsq_eps <- c(eta_sigsq_eps_1, eta_sigsq_eps_2)

  G_sigsq_eps <- "full"

  eta_out <- list(eta_nu, eta_zeta, eta_sigsq_eps)
  names(eta_out) <- c(
    "p(Y|nu,zeta,sigsq_eps)->nu",
    "p(Y|nu,zeta,sigsq_eps)->zeta",
    "p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
  )

  G_out <- G_sigsq_eps
  names(G_out) <- "p(Y|nu,zeta,sigsq_eps)->sigsq_eps"

  ans <- list(eta_out, G_out)
  names(ans) <- c("eta", "G")

  return(ans)
}

# mlfpc_lik_frag <- function(eta_in, G_in, C, Y, T_vec, L_1, L_2) {

  # # order of eta_in:
  # # 1. nu -> p(Y|nu,zeta,sigsq_eps)
  # # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # # 3. zeta -> p(Y|nu,zeta,sigsq_eps)
  # # 4. p(Y|nu,zeta,sigsq_eps) -> zeta
  # # 5. sigsq_eps -> p(Y|nu,zeta,sigsq_eps)
  # # 6. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # # order of G_in:
  # # 1. sigsq_eps -> p(Y|nu,zeta,sigsq_eps)
  # # 2. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # # order of eta_out:
  # # 1. p(Y|nu,zeta,sigsq_eps) -> nu
  # # 2. p(Y|nu,zeta,sigsq_eps) -> zeta
  # # 3. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # # order of G_out:
  # # 1. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # if(!is.list(eta_in)) {

    # stop("eta_in must be a list")
  # }

  # if(length(eta_in)!=6) {

    # stop("eta_in must have length 6")
  # }

  # if(!is.list(G_in)) {

    # stop("G_in must be a list")
  # }

  # if(length(G_in)!=2) {

    # stop("G_in must have length 2")
  # }

  # if(!all(Reduce("c", lapply(G_in, is.character)))) {

    # stop("all of the elements of G_in must be characters")
  # } else {

    # if(!all(Reduce("c", lapply(G_in, function(x) x=="diag" || x=="full")))) {

      # stop("The input graphs must be 'diag' or 'full'.")
    # }
  # }

  # if(G_in[[1]]==G_in[[2]]) {

    # G_sigsq_eps <- G_in[[1]]
  # } else {

    # stop("The graph messages involving sigsq_eps are not identical")
  # }

  # if(!is_int(L_1)) {

    # stop("L_1 must be a whole number.")
  # } else {

    # if(L_1 < 1) {

      # stop("There must be at least one basis function for the first level.")
    # }
  # }

  # if(!is_int(L_2)) {

    # stop("L_2 must be a whole number.")
  # } else {

    # if(L_2 < 1) {

      # stop("There must be at least one basis function for the second level.")
    # }
  # }

  # L <- L_1 + L_2
  # N <- ncol(eta_in[[3]])
  # l_eta_nu <- length(eta_in[[1]])
  # l_eta_zeta <- nrow(eta_in[[3]])
  # d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  # K <- d/(L+1) - 2
  # l_zeta_1 <- L_1 + 0.5*L_1*(L_1 + 1)
  # l_zeta_2 <- L_2 + 0.5*L_2*(L_2 + 1) + L_1*L_2
  # M <- (l_eta_zeta - l_zeta_1)/l_zeta_2

  # D_zeta_1 <- duplication.matrix(L_1)
  # D_zeta_2 <- duplication.matrix(L_2)

  # YTY <- matrix(NA, N, M)
  # CTY <- vector("list", length = N)
  # CTC <- vector("list", length = N)
  # for(i in 1:N) {

    # CTY[[i]] <- matrix(NA, K + 2, M)
    # CTC[[i]] <- vector("list", length = M)
    # for(j in 1:M) {

      # YTY[i, j] <- cprod(Y[[i]][[j]])
      # CTY[[i]][, j] <- cprod(C[[i]][[j]], Y[[i]][[j]])
      # CTC[[i]][[j]] <- crossprod(C[[i]][[j]])
    # }
  # }

  # eta_in_nu <- list(eta_in[[1]], eta_in[[2]])
  # q_nu <- gauss_q(eta_in_nu, use_vech = FALSE)
  # E_q_nu <- q_nu[[1]]
  # Cov_q_nu <- q_nu[[2]]

  # eta_in_sigsq_eps <- list(eta_in[[5]], eta_in[[6]])
  # q_sigsq_eps <- igw_q(eta_in_sigsq_eps, G_in)
  # E_q_recip_sigsq_eps <- q_sigsq_eps[[3]]

  # eta_1_sum <- 0
  # eta_2_sum <- 0
  # E_q_zeta_1 <- vector("list", length = N)
  # Cov_q_zeta_1 <- vector("list", length = N)
  # E_q_zeta_2 <- vector("list", length = N)
  # Cov_q_zeta_2 <- vector("list", length = N)
  # Cov_q_zeta_12 <- vector("list", length = N)
  # for(i in 1:N) {

    # eta_in_zeta <- list(
      # eta_in[[3]][, i],
      # eta_in[[4]][, i]
    # )
    # q_zeta <- two_level_nat_to_comm_parms(L_1, L_2, M, eta_in_zeta)
    # E_q_zeta_1[[i]] <- q_zeta[[1]]
    # Cov_q_zeta_1[[i]] <- q_zeta[[2]]
    # E_q_zeta_2[[i]] <- q_zeta[[3]]
    # Cov_q_zeta_2[[i]] <- q_zeta[[4]]
    # Cov_q_zeta_12[[i]] <- q_zeta[[5]]

    # for(j in 1:M) {

      # E_q_zeta <- c(E_q_zeta_1[[i]], E_q_zeta_2[[i]][[j]])
      # E_q_zeta_tilde <- c(1, E_q_zeta)

      # top_block <- cbind(Cov_q_zeta_1[[i]], Cov_q_zeta_12[[i]][[j]])
      # bottom_block <- cbind(t(Cov_q_zeta_12[[i]][[j]]), Cov_q_zeta_2[[i]][[j]])
      # Cov_q_zeta <- rbind(top_block, bottom_block)
      # Cov_q_zeta_tilde <- adiag(0, Cov_q_zeta)

      # E_q_tcross_zeta_tilde <- Cov_q_zeta_tilde + tcrossprod(E_q_zeta_tilde)

      # sum_val <- cprod(kronecker(t(E_q_zeta_tilde), C[[i]][[j]]), Y[[i]][[j]])
      # eta_1_sum <- eta_1_sum + sum_val

      # sum_val <- kronecker(E_q_tcross_zeta_tilde, CTC[[i]][[j]])
      # eta_2_sum <- eta_2_sum + sum_val
    # }
  # }
  # eta_1 <- E_q_recip_sigsq_eps*eta_1_sum
  # eta_2 <- -0.5*E_q_recip_sigsq_eps*eta_2_sum
  # eta_nu <- c(eta_1, eta_2)

  # inds_mat <- matrix(1:d, K + 2, L + 1)
  # mu_inds <- inds_mat[, 1]
  # psi_1_inds <- inds_mat[, 1:L_1 + 1]
  # psi_2_inds <- inds_mat[, 1:L_2 + L_1 + 1]

  # E_q_V <- matrix(E_q_nu, K+2, L + 1)
  # E_q_nu_mu <- E_q_V[, 1]
  # E_q_V_psi_1 <- E_q_V[, 1:L_1 + 1]
  # E_q_nu_psi_1 <- as.vector(E_q_V_psi_1)
  # E_q_V_psi_2 <- E_q_V[, 1:L_2 + L_1 + 1]
  # E_q_nu_psi_2 <- as.vector(E_q_V_psi_2)

  # Cov_q_nu_mu <- Cov_q_nu[mu_inds, mu_inds]
  # Cov_q_nu_psi_1 <- Cov_q_nu[as.vector(psi_1_inds), as.vector(psi_1_inds)]
  # Cov_q_nu_psi_2 <- Cov_q_nu[as.vector(psi_2_inds), as.vector(psi_2_inds)]
  # Cov_q_nu_mu_psi_1 <- Cov_q_nu[mu_inds, as.vector(psi_1_inds)]
  # Cov_q_nu_mu_psi_2 <- Cov_q_nu[mu_inds, as.vector(psi_2_inds)]
  # Cov_q_nu_psi_21 <- Cov_q_nu[as.vector(psi_2_inds), as.vector(psi_1_inds)]

  # E_q_h_mu <- E_cprod(
    # E_q_nu_mu, Cov_q_nu_mu, E_q_nu_mu,
    # Reduce("+", lapply(CTC, Reduce, f = "+"))
  # )

  # eta_zeta_1 <- matrix(NA, l_zeta_1, N)
  # eta_zeta_2 <- matrix(NA, M*l_zeta_2, N)
  # E_q_h_mu_psi_1 <- vector("list", length = N)
  # E_q_H_psi_11 <- vector("list", length = N)
  # E_q_h_mu_psi_2 <- vector("list", length = N)
  # E_q_H_psi_22 <- vector("list", length = N)
  # E_q_H_psi_12 <- vector("list", length = N)
  # for(i in 1:N) {

    # CTY_i <- rowSums(CTY[[i]])
    # CTC_i <- Reduce("+", CTC[[i]])

    # E_q_h_mu_psi_1[[i]] <- E_h(L_1, E_q_nu_psi_1, Cov_q_nu_mu_psi_1, E_q_nu_mu, CTC_i)
    # E_q_H_psi_11[[i]] <- E_H(L_1, L_1, E_q_nu_psi_1, Cov_q_nu_psi_1, E_q_nu_psi_1, CTC_i)
    # eta_11 <- E_q_recip_sigsq_eps*(cprod(E_q_V_psi_1, CTY_i) - E_q_h_mu_psi_1[[i]])
    # eta_12 <- -0.5*E_q_recip_sigsq_eps*cprod(D_zeta_1, as.vector(E_q_H_psi_11[[i]]))
    # eta_zeta_1[, i] <- c(eta_11, eta_12)

    # E_q_h_mu_psi_2[[i]] <- vector("list", length = M)
    # E_q_H_psi_22[[i]] <- vector("list", length = M)
    # E_q_H_psi_12[[i]] <- vector("list", length = M)
    # eta_21 <- matrix(NA, L_2, M)
    # eta_22 <- matrix(NA, 0.5*L_2*(L_2 + 1), M)
    # eta_23 <- matrix(NA, L_1*L_2, M)
    # for(j in 1:M) {

      # E_q_h_mu_psi_2[[i]][[j]] <- E_h(
        # L_2, E_q_nu_psi_2, Cov_q_nu_mu_psi_2,
        # E_q_nu_mu, CTC[[i]][[j]]
      # )

      # E_q_H_psi_22[[i]][[j]] <- E_H(
        # L_2, L_2, E_q_nu_psi_2, Cov_q_nu_psi_2,
        # E_q_nu_psi_2, CTC[[i]][[j]]
      # )

      # E_q_H_psi_12[[i]][[j]] <- E_H(
        # L_1, L_2, E_q_nu_psi_1, Cov_q_nu_psi_21,
        # E_q_nu_psi_2, CTC[[i]][[j]]
      # )

      # eta_21[, j] <- E_q_recip_sigsq_eps*(
        # cprod(E_q_V_psi_2, CTY[[i]][, j]) - E_q_h_mu_psi_2[[i]][[j]]
      # )
      # eta_22[, j] <- -0.5*E_q_recip_sigsq_eps*cprod(
        # D_zeta_2,
        # as.vector(E_q_H_psi_22[[i]][[j]])
      # )
      # eta_23[, j] <- -E_q_recip_sigsq_eps*as.vector(E_q_H_psi_12[[i]][[j]])
    # }
    # eta_zeta_2[, i] <- as.vector(rbind(eta_21, eta_22, eta_23))
  # }
  # eta_zeta <- rbind(eta_zeta_1, eta_zeta_2)

  # term_1 <- sum(YTY)
  # term_2 <- -2*cprod(E_q_nu_mu, Reduce("+", lapply(CTY, rowSums)))
  # term_3 <- 0
  # term_4 <- 0
  # for(i in 1:N) {

    # sum_val <- -2*cprod(E_q_V_psi_1 %*% E_q_zeta_1[[i]], rowSums(CTY[[i]]))
    # term_3 <- term_3 + sum_val

    # for(j in 1:M) {

      # sum_val <- -2*cprod(E_q_V_psi_2 %*% E_q_zeta_2[[i]][[j]], CTY[[i]][, j])
      # term_4 <- term_4 + sum_val
    # }
  # }
  # y_line <- term_1 + term_2 + term_3 + term_4

  # term_1 <- E_q_h_mu
  # term_2 <- 0
  # term_3 <- 0
  # for(i in 1:N) {

    # sum_val <- 2*cprod(E_q_zeta_1[[i]], E_q_h_mu_psi_1[[i]])
    # term_2 <- term_2 + sum_val

    # for(j in 1:M) {

      # sum_val <- 2*cprod(E_q_zeta_2[[i]][[j]], E_q_h_mu_psi_2[[i]][[j]])
      # term_3 <- term_3 + sum_val
    # }
  # }
  # mu_line <- term_1 + term_2 + term_3

  # term_1 <- 0
  # term_2 <- 0
  # for(i in 1:N) {

    # E_vec <- E_q_zeta_1[[i]]
    # Cov_mat <- Cov_q_zeta_1[[i]]
    # A_mat <- E_q_H_psi_11[[i]]
    # sum_val <- E_cprod(E_vec, Cov_mat, E_vec, A_mat)
    # term_1 <- term_1 + sum_val

    # for(j in 1:M) {

      # E_vec_1 <- E_q_zeta_1[[i]]
      # Cov_mat <- t(Cov_q_zeta_12[[i]][[j]])
      # E_vec_2 <- E_q_zeta_2[[i]][[j]]
      # A_mat <- E_q_H_psi_12[[i]][[j]]
      # sum_val <- 2*E_cprod(E_vec_1, Cov_mat, E_vec_2, A_mat)
      # term_2 <- term_2 + sum_val
    # }
  # }
  # zeta_1_line <- term_1 + term_2

  # term_1 <- 0
  # for(i in 1:N) {

    # for(j in 1:M) {

      # E_vec <- E_q_zeta_2[[i]][[j]]
      # Cov_mat <- Cov_q_zeta_2[[i]][[j]]
      # A_mat <- E_q_H_psi_22[[i]][[j]]
      # sum_val <- E_cprod(E_vec, Cov_mat, E_vec, A_mat)
      # term_1 <- term_1 + sum_val
    # }
  # }
  # zeta_2_line <- term_1

  # E_q_norm <- y_line + mu_line + zeta_1_line + zeta_2_line

  # eta_sigsq_eps_1 <- -0.5*sum(sapply(T_vec, sum))
  # eta_sigsq_eps_2 <- -0.5*E_q_norm
  # eta_sigsq_eps <- c(eta_sigsq_eps_1, eta_sigsq_eps_2)

  # G_sigsq_eps <- "full"

  # eta_out <- list(eta_nu, eta_zeta, eta_sigsq_eps)
  # names(eta_out) <- c(
    # "p(Y|nu,zeta,sigsq_eps)->nu",
    # "p(Y|nu,zeta,sigsq_eps)->zeta",
    # "p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
  # )

  # G_out <- G_sigsq_eps
  # names(G_out) <- "p(Y|nu,zeta,sigsq_eps)->sigsq_eps"

  # ans <- list(eta_out, G_out)
  # names(ans) <- c("eta", "G")

  # return(ans)
# }

mfpc_lik_frag <- function(eta_in, G_in, C, Y, L) {

  # order of eta_in:
  # 1. nu -> p(Y|nu,zeta,sigsq_eps)
  # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # 3. zeta -> p(Y|nu,zeta,sigsq_eps)
  # 4. p(Y|nu,zeta,sigsq_eps) -> zeta
  # 5. sigsq_eps -> p(Y|nu,zeta,sigsq_eps)
  # 6. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of G_in:
  # 1. sigsq_eps -> p(Y|nu,zeta,sigsq_eps)
  # 2. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of eta_out:
  # 1. p(Y|nu,zeta,sigsq_eps) -> nu
  # 2. p(Y|nu,zeta,sigsq_eps) -> zeta
  # 3. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of G_out:
  # 1. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in)!=6) {

    stop("eta_in must have length 6")
  }

  if(!is.list(G_in)) {

    stop("G_in must be a list")
  }

  if(length(G_in)!=2) {

    stop("G_in must have length 2")
  }

  if(!all(is.character(unlist(G_in)))) {

    stop("all of the elements of G_in must be characters")
  } else {

    if(!all(unlist(G_in) == "diag" | unlist(G_in) == "full")) {

      stop("The input graphs must be 'diag' or 'full'.")
    }
  }

  if(all(unlist(G_in[[1]])==unlist(G_in[[2]]))) {

    G_sigsq_eps <- G_in[[1]]
  } else {

    stop("The graph messages involving sigsq_eps are not identical")
  }

  if(!is_int(L)) {

    stop("L must be a whole number.")
  } else {

    if(L<1) {

      stop("There must be at least one basis function.")
    }
  }

  N <- ncol(eta_in[[3]])
  p <- length(eta_in[[1]])
  l_eta_nu <- length(eta_in[[1]][[1]])
  l_eta_zeta <- L + 0.5*L*(L + 1)
  K <- sapply(C[[1]], ncol) - 2
  d <- (K+2)*(L+1)
  D_L <- duplication.matrix(L)
  n <- Reduce(rbind, lapply(Y, function(x) sapply(x, length)))

  YTY <- matrix(NA, N, p)
  CTY <- vector("list", length = N)
  CTC <- vector("list", length = N)
  for(i in 1:N) {

    CTY[[i]] <- vector("list", length = p)
    CTC[[i]] <- vector("list", length = p)
    for(j in 1:p) {

      YTY[i, j] <- cprod(Y[[i]][[j]])
      CTY[[i]][[j]] <- cprod(Y[[i]][[j]], C[[i]][[j]])
      CTC[[i]][[j]] <- crossprod(C[[i]][[j]])
    }
  }

  E_q_zeta <- vector("list", length = N)
  Cov_q_zeta <- vector("list", length = N)
  E_q_tcross_zeta <- vector("list", length = N)
  for(i in 1:N) {

    eta_in_zeta <- list(
      eta_in[[3]][, i],
      eta_in[[4]][, i]
    )
    q_zeta <- gauss_q(eta_in_zeta, use_vech = TRUE)
    E_q_zeta[[i]] <- q_zeta[[1]]
    Cov_q_zeta[[i]] <- q_zeta[[2]]
    E_q_tcross_zeta[[i]] <- Cov_q_zeta[[i]] + tcrossprod(E_q_zeta[[i]])
  }

  E_q_recip_sigsq_eps <- rep(NA, p)
  for(j in 1:p) {

    eta_in_sigsq_eps <- list(eta_in[[5]][[j]], eta_in[[6]][[j]])
    G_in_sigsq_eps <- c(G_in[[1]][[j]], G_in[[2]][[j]])
    q_sigsq_eps <- igw_q(eta_in_sigsq_eps, G_in_sigsq_eps)
    E_q_recip_sigsq_eps[j] <- q_sigsq_eps[[3]]
  }

  eta_nu <- vector("list", length = p)
  E_q_nu <- vector("list", length = p)
  Cov_q_nu <- vector("list", length = p)
  for(j in 1:p) {

    eta_in_nu <- list(eta_in[[1]][[j]], eta_in[[2]][[j]])
    q_nu <- gauss_q(eta_in_nu, use_vech = FALSE)
    E_q_nu[[j]] <- q_nu[[1]]
    Cov_q_nu[[j]] <- q_nu[[2]]

    eta_1_sum <- 0
    eta_2_sum <- 0
    for(i in 1:N) {

      E_q_zeta_tilde <- c(1, E_q_zeta[[i]])
      Cov_q_zeta_tilde <- adiag(0, Cov_q_zeta[[i]])
      E_q_tcross_zeta_tilde <- Cov_q_zeta_tilde + tcrossprod(E_q_zeta_tilde)

      sum_val <- cprod(kronecker(t(E_q_zeta_tilde), C[[i]][[j]]), Y[[i]][[j]])
      eta_1_sum <- eta_1_sum + sum_val

      sum_val <- as.vector(kronecker(E_q_tcross_zeta_tilde, CTC[[i]][[j]]))
      eta_2_sum <- eta_2_sum + sum_val
    }

    eta_1 <- E_q_recip_sigsq_eps[j]*eta_1_sum
    eta_2 <- -0.5*E_q_recip_sigsq_eps[j]*eta_2_sum
    eta_nu[[j]] <- c(eta_1, eta_2)
  }

  E_q_V <- vector("list", length = p)
  E_q_nu_mu <- vector("list", length = p)
  E_q_V_psi <- vector("list", length = p)
  E_q_nu_psi <- vector("list", length = p)
  Cov_q_nu_mu <- vector("list", length = p)
  Cov_q_nu_psi <- vector("list", length = p)
  Cov_q_nu_mu_psi <- vector("list", length = p)
  for(j in 1:p) {

    inds_mat <- matrix(1:d[j], K[j] + 2, L + 1)
    mu_inds <- inds_mat[, 1]
    psi_inds <- as.vector(inds_mat[, 1:L + 1])

    E_q_V[[j]] <- matrix(E_q_nu[[j]], K[j] + 2, L + 1)
    E_q_nu_mu[[j]] <- E_q_V[[j]][, 1]
    E_q_V_psi[[j]] <- E_q_V[[j]][, 1:L + 1]
    E_q_nu_psi[[j]] <- as.vector(E_q_V_psi[[j]])

    Cov_q_nu_mu[[j]] <- Cov_q_nu[[j]][mu_inds, mu_inds]
    Cov_q_nu_psi[[j]] <- Cov_q_nu[[j]][psi_inds, psi_inds]
    Cov_q_nu_mu_psi[[j]] <- Cov_q_nu[[j]][mu_inds, psi_inds]
  }

  E_q_h_mu <- matrix(NA, N, p)
  E_q_h_mu_psi <- vector("list", length = N)
  E_q_H_psi <- vector("list", length = N)
  eta_zeta <- matrix(NA, l_eta_zeta, N)
  for(i in 1:N) {

    E_q_h_mu_psi[[i]] <- matrix(NA, L, p)
    E_q_H_psi[[i]] <- vector("list", length = p)
    eta_1_sum <- 0
    eta_2_sum <- 0
    for(j in 1:p) {

      E_q_h_mu[i, j] <- E_cprod(
        E_q_nu_mu[[j]], Cov_q_nu_mu[[j]],
        E_q_nu_mu[[j]], CTC[[i]][[j]]
      )

      E_q_h_mu_psi[[i]][, j] <- E_h(
        L, E_q_nu_psi[[j]], Cov_q_nu_mu_psi[[j]],
        E_q_nu_mu[[j]], CTC[[i]][[j]]
      )

      E_q_H_psi[[i]][[j]] <- E_H(
        L, L, E_q_nu_psi[[j]], Cov_q_nu_psi[[j]],
        E_q_nu_psi[[j]], CTC[[i]][[j]]
      )

      freq_scores <- cprod(E_q_V_psi[[j]], CTY[[i]][[j]]) - E_q_h_mu_psi[[i]][, j]
      sum_val <- E_q_recip_sigsq_eps[j]*freq_scores
      eta_1_sum <- eta_1_sum + sum_val

      sum_val <- E_q_recip_sigsq_eps[j]*cprod(D_L, as.vector(E_q_H_psi[[i]][[j]]))
      eta_2_sum <- eta_2_sum + sum_val
    }
    eta_1 <- eta_1_sum
    eta_2 <- -0.5*eta_2_sum
    eta_zeta[, i] <- c(eta_1, eta_2)
  }

  eta_sigsq_eps <- vector("list", length = p)
  for(j in 1:p) {

    eta_2_sum <- 0
    for(i in 1:N) {

      summands <- rep(NA, 6)
      summands[1] <- YTY[i, j]
      summands[2] <- - 2*cprod(E_q_nu_mu[[j]], CTY[[i]][[j]])
      summands[3] <- -2*cprod(E_q_V_psi[[j]] %*% E_q_zeta[[i]], CTY[[i]][[j]])
      summands[4] <- E_q_h_mu[i, j]
      summands[5] <- 2*cprod(E_q_zeta[[i]], E_q_h_mu_psi[[i]][, j])
      summands[6] <- tr(E_q_tcross_zeta[[i]] %*% E_q_H_psi[[i]][[j]])
      sum_val <- sum(summands)
      eta_2_sum <- eta_2_sum + sum_val
    }

    eta_1 <- -0.5*sum(n[, j])
    eta_2 <- -0.5*eta_2_sum
    eta_sigsq_eps[[j]] <- c(eta_1, eta_2)
  }

  G_sigsq_eps <- replicate(p, "full", simplify = FALSE)

  eta_out <- list(eta_nu, eta_zeta, eta_sigsq_eps)
  names(eta_out) <- c(
    "p(Y|nu,zeta,sigsq_eps)->nu",
    "p(Y|nu,zeta,sigsq_eps)->zeta",
    "p(Y|nu,zeta,sigsq_eps)->sigsq_eps"
  )

  G_out <- list(G_sigsq_eps)
  names(G_out) <- "p(Y|nu,zeta,sigsq_eps)->sigsq_eps"

  ans <- list(eta_out, G_out)
  names(ans) <- c("eta", "G")

  return(ans)
}

logistic_fpc_lik_frag <- function(eta_in, C, Y, L) {

  # order of eta_in:
  # 1. nu -> p(Y|nu,zeta)
  # 2. p(Y|nu,zeta) -> nu
  # 3. zeta -> p(Y|nu,zeta)
  # 4. p(Y|nu,zeta) -> zeta

  # order of eta_out:
  # 1. p(Y|nu,zeta,sigsq_eps) -> nu
  # 2. p(Y|nu,zeta,sigsq_eps) -> zeta

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in)!=4) {

    stop("eta_in must have length 4")
  }

  if(!is_int(L)) {

    stop("L must be a whole number.")
  } else {

    if(L<1) {

      stop("There must be at least one basis function.")
    }
  }

  N <- ncol(eta_in[[3]])
  l_eta_nu <- length(eta_in[[1]])
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/(L+1) - 2
  D_L <- duplication.matrix(L)

  eta_in_nu <- list(eta_in[[1]], eta_in[[2]])
  q_nu <- gauss_q(eta_in_nu, use_vech=FALSE)
  mu_q_nu <- q_nu[[1]]
  Sigma_q_nu <- q_nu[[2]]

  mu_q_zeta <- vector("list", length=N)
  Sigma_q_zeta <- vector("list", length=N)
  for(i in 1:N) {

    eta_in_zeta <- list(
      eta_in[[3]][,i],
      eta_in[[4]][,i]
    )
    q_zeta <- gauss_q(eta_in_zeta, use_vech=TRUE)
    mu_q_zeta[[i]] <- q_zeta[[1]]
    Sigma_q_zeta[[i]] <- q_zeta[[2]]
  }

  mu_inds <- 1:(K+2)
  mu_q_nu_mu <- mu_q_nu[mu_inds]
  Sigma_q_nu_mu <- Sigma_q_nu[mu_inds, mu_inds]

  psi_inds <- (1:d)[-mu_inds]
  psi_groups <- ceiling(psi_inds/(K+2)) - 1
  psi_inds <- split(psi_inds, psi_groups)

  M_q_V_psi <- matrix(NA, K+2, L)
  for(l in 1:L) {

    psi_l_inds <- psi_inds[[l]]
    M_q_V_psi[,l] <- mu_q_nu[psi_l_inds]
  }

  Sigma_q_nu_blocks <- vector("list", length = L+1)
  Sigma_q_nu_blocks[[1]] <- vector("list", length = L+1)
  Sigma_q_nu_blocks[[1]][[1]] <- Sigma_q_nu[mu_inds, mu_inds]
  for(l in 1:L) {

    l_inds <- psi_inds[[l]]
    Sigma_q_nu_blocks[[1]][[l+1]] <- Sigma_q_nu[mu_inds, l_inds]

    Sigma_q_nu_blocks[[l+1]] <- vector("list", length = L+1)
    Sigma_q_nu_blocks[[l+1]][[1]] <- Sigma_q_nu[l_inds, mu_inds]
    for(l_p in 1:L) {

      l_p_inds <- psi_inds[[l_p]]
      Sigma_q_nu_blocks[[l+1]][[l_p+1]] <- Sigma_q_nu[l_inds, l_p_inds]
    }
  }

  xi <- vector("list", length=N)
  term_1 <- Sigma_q_nu_mu + tcrossprod(mu_q_nu_mu)
  eta_1_sum <- 0
  eta_2_sum <- 0
  eta_zeta_1 <- matrix(NA, L, N)
  eta_zeta_2 <- matrix(NA, 0.5*L*(L+1), N)
  for(i in 1:N) {

    term_2 <- 0
    for(l in 1:L) {

      Cov_mat <- Sigma_q_nu_blocks[[1]][[l+1]]
      mu_outer_mat <- tcrossprod(mu_q_nu_mu, M_q_V_psi[,l])
      prod_term <- Cov_mat + mu_outer_mat

      term_2 <- term_2 + mu_q_zeta[[i]][l]*prod_term
    }

    term_3 <- t(term_2)

    term_4 <- 0
    for(l in 1:L) {

      for(l_p in 1:L) {

        cov_val <- Sigma_q_zeta[[i]][l, l_p]
        mu_prod <- mu_q_zeta[[i]][l]*mu_q_zeta[[i]][l_p]
        E_q_zeta_prod <- cov_val + mu_prod

        Cov_mat <- Sigma_q_nu_blocks[[l+1]][[l_p+1]]
        nu_outer_term <- tcrossprod(M_q_V_psi[,l], M_q_V_psi[,l_p])
        E_q_outer_nu <- Cov_mat + nu_outer_term

        sum_val <- E_q_zeta_prod*E_q_outer_nu
        term_4 <- term_4 + sum_val
      }
    }

    E_q_outer_theta <- term_1 + term_2 + term_3 + term_4
    xi[[i]] <- sqrt(diag(tcrossprod(C[[i]]%*%E_q_outer_theta, C[[i]])))

    A_xi <- -tanh(xi[[i]]/2)/(4*xi[[i]])

    mu_q_zeta_tilde <- c(1, mu_q_zeta[[i]])
    Sigma_q_zeta_tilde <- blkdiag(matrix(0), Sigma_q_zeta[[i]])
    E_q_outer_zeta_tilde <- Sigma_q_zeta_tilde + tcrossprod(mu_q_zeta_tilde)

    sum_val <- cprod(kronecker(t(mu_q_zeta_tilde), C[[i]]), Y[[i]] - 0.5)
    eta_1_sum <- eta_1_sum + sum_val

    M <- crossprod(C[[i]], diag(A_xi)%*%C[[i]])
    sum_val <- as.vector(kronecker(E_q_outer_zeta_tilde, M))
    eta_2_sum <- eta_2_sum + sum_val

    E_q_h_xi <- rep(NA, L)
    E_q_H_xi <- matrix(NA, L, L)
    for(l in 1:L) {

      tr_term <- tr(Sigma_q_nu_blocks[[1]][[l+1]]%*%M)
      cprod_term <- cprod(M_q_V_psi[,l], M%*%mu_q_nu_mu)
      E_q_h_xi[l] <- tr_term + cprod_term

      for(l_p in 1:L) {

        tr_term <- tr(Sigma_q_nu_blocks[[l_p+1]][[l+1]]%*%M)
        cprod_term <- cprod(M_q_V_psi[,l], M%*%M_q_V_psi[,l_p])
        E_q_H_xi[l, l_p] <- tr_term + cprod_term
      }
    }

    eta_zeta_1[,i] <- 2*E_q_h_xi + cprod(C[[i]]%*%M_q_V_psi, Y[[i]] - 0.5)
    eta_zeta_2[,i] <- cprod(D_L, as.vector(E_q_H_xi))
  }
  eta_1 <- eta_1_sum
  eta_2 <- eta_2_sum
  eta_nu <- c(eta_1, eta_2)
  eta_zeta <- rbind(eta_zeta_1, eta_zeta_2)

  eta_out <- list(eta_nu, eta_zeta)
  names(eta_out) <- c(
    "p(Y|nu,zeta,sigsq_eps)->nu",
    "p(Y|nu,zeta,sigsq_eps)->zeta"
  )

  return(eta_out)
}

##########################################
#
#  NATURAL  PARAMETER  TO  COMMON
#  PARAMETER  UPDATES  AND  EXPECTATION
#  OF  SUFFFICIENT  STATITISTIC
#  COMPUTATIONS.
#
##########################################

igw_q <- function(eta_in, G) {

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in)!=2) {

    stop("eta_in must have length 2")
  }

  if(!is.vector(G)) {

    stop("G must be a vector")
  }

  if(length(G)!=2) {

    stop("G must have length 2")
  }

  if(length(unique(G))!=1) {

    stop("graph messages must be identical")

    if((unique(G) != "full") & (unique(G) != "diag")) {

      stop("G can only be 'full' or 'diag'")
    }
  }

  eta_q <- eta_in[[1]] + eta_in[[2]]
  eta_q_1 <- eta_q[1]
  eta_q_2 <- eta_q[-1]

  l <- length(eta_q_2)
  d <- (sqrt(1 + 8*l) - 1)/2

  if(d==1) {

    D_d <- 1
  } else {

    D_d <- duplication.matrix(d)
  }

  D_d_plus <- solve(crossprod(D_d), t(D_d))

  G <- unique(G)

  xi <- -2*eta_q_1 - 2
  eta_2_mat <- matrix(as.vector(crossprod(D_d_plus, eta_q_2)), d, d)
  Lambda <- -2*eta_2_mat

  if(G=="full") {

    E_X_inv <- (eta_q_1 + 0.5*(d+1))*solve(eta_2_mat)

    sum_val <- 0
    for(j in 1:d) {

      digamma_arg <- (xi + 2 - d - j)/2
      sum_val <- sum_val + digamma(digamma_arg)
    }
    E_log_det_X <- determinant(Lambda/2, logarithm=TRUE)$modulus[1] - sum_val

  } else {

    E_X_inv <- (eta_q_1 + 1)*solve(eta_2_mat)

    sum_val <- 0
    for(j in 1:d) {

      sum_val <- sum_val + log(Lambda[j,j]/2)
    }
    E_log_det_X <- sum_val - d*digamma(xi/2)
  }

  if(d==1) {

    E_X_inv <- as.vector(E_X_inv)
    Lambda <- as.vector(Lambda)
  }

  ans <- list(xi, Lambda, E_X_inv, E_log_det_X)
  names(ans) <- c("xi", "Lambda", "E_X_inv", "E_log_det_X")
  return(ans)
}

gauss_q <- function(eta_in, use_vech=TRUE) {

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in)!=2) {

    stop("eta_in must have length 2")
  }

  if(!is.logical(use_vech)) {

    stop("use_vech must be of type logical")
  }

  eta_q <- eta_in[[1]] + eta_in[[2]]

  l_eta_q <- length(eta_q)

  if(use_vech) {

    d <- (sqrt(9 + 8*l_eta_q) - 3)/2
  } else {

    d <- (sqrt(1 + 4*l_eta_q) - 1)/2
  }

  eta_q_1 <- eta_q[1:d]
  eta_q_2 <- eta_q[-c(1:d)]

  if(use_vech) {

    D_d <- duplication.matrix(d)
    D_d_plus <- solve(crossprod(D_d), t(D_d))
    eta_q_2 <- as.vector(crossprod(D_d_plus, eta_q_2))
  }

  eta_q_2_mat <- matrix(eta_q_2, d, d)
  Sigma_q_x <- -0.5*solve(eta_q_2_mat)
  mu_q_x <- as.vector(Sigma_q_x%*%eta_q_1)

  if(d==1) {

    Sigma_q_x <- as.vector(Sigma_q_x)
  }

  ans <- list(mu_q_x, Sigma_q_x)
  names(ans) <- c("mu_q_x", "Sigma_q_x")
  return(ans)
}

##########################################
#
#  ENTROPY  FUNCTIONS
#
##########################################

entropy_igw <- function(eta_in, G_in) {

  q_X <- igw_q(eta_in, G_in)
  xi <- q_X[[1]]
  Lambda <- q_X[[2]]
  E_X_inv <- q_X[[3]]
  E_log_det_X <- q_X[[4]]

  G <- unique(G_in)

  if(is.vector(Lambda)) {

    d <- 1
    Lambda <- matrix(Lambda)
  } else {

    d <- nrow(Lambda)
  }

  if(G=="full") {

    term_1 <- d/2*(xi - d + 1)*log(2)
    term_2 <- d/4*(d - 1)*log(pi)
    term_3 <- sum(lgamma(1/2*(xi - d + c(1:d)) + 1))
    term_4 <- -1/2*(xi - d + 1)*determinant(Lambda, logarithm=TRUE)$modulus[1]
    term_5 <- 1/2*(xi + 2)* E_log_det_X
    term_6 <- 1/2*tr(Lambda%*%E_X_inv)

    ans <- term_1 + term_2 + term_3 + term_4 + term_5 + term_6

  } else {

    term_1 <- d*xi/2*log(2)
    term_2 <- d*lgamma(xi/2)
    term_3 <- -xi/2*sum(log(diag(Lambda)))
    term_4 <- 1/2*(xi + 2)*E_log_det_X
    term_5 <- 1/2*tr(Lambda%*%E_X_inv)

    ans <- term_1 + term_2 + term_3 + term_4 + term_5
  }

  return(ans)
}

entropy_gauss <- function(eta_in, use_vech=TRUE) {

  q_x <- gauss_q(eta_in, use_vech)
  mu <- q_x[[1]]
  Sigma <- q_x[[2]]

  d <- length(mu)

  ans <- d/2*(1 + log(2*pi)) + 1/2*determinant(Sigma, logarithm=TRUE)$modulus[1]
  return(ans)
}

entropy_two_lev_gauss <- function(eta_in, p, q, m) {

  q_x <- two_level_nat_to_comm_parms(p, q, m, eta_in)
  E_q_x_1 <- q_x[[1]]
  Cov_q_x_1 <- q_x[[2]]
  E_q_x_2 <- q_x[[3]]
  Cov_q_x_2 <- q_x[[4]]
  Cov_q_x_12 <- q_x[[5]]

  d <- length(E_q_x_1) + sum(sapply(E_q_x_2, length))

  log_det_Cov_q_x <- determinant(solve(Cov_q_x_1), logarithm = TRUE)$modulus[1]
  for(i in 1:m) {

    A_22_inv <- Cov_q_x_2[[i]] - crossprod(
      Cov_q_x_12[[i]],
      solve(Cov_q_x_1, Cov_q_x_12[[i]])
    )

    A_22 <- solve(A_22_inv)

    sum_val <- determinant(A_22, logarithm = TRUE)$modulus[1]
    log_det_Cov_q_x <- log_det_Cov_q_x + sum_val
  }

  ans <- d/2*(1 + log(2*pi)) + 1/2*log_det_Cov_q_x
  return(ans)
}

entropy_two_lev_fpc_prod <- function(eta_in, C, L_1, L_2) {

  # order of eta_in:
  # 1. p(nu|Sigma_nu) -> nu
  # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # 3. p(zeta)->zeta
  # 4. p(Y|nu,zeta,sigsq_eps)->zeta

  L <- L_1 + L_2
  N <- ncol(eta_in[[3]])
  l_eta_nu <- length(eta_in[[1]])
  l_eta_zeta <- nrow(eta_in[[3]])
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/(L+1) - 2
  l_zeta_1 <- L_1 + 0.5*L_1*(L_1 + 1)
  l_zeta_2 <- L_2 + 0.5*L_2*(L_2 + 1) + L_1*L_2
  M <- (l_eta_zeta - l_zeta_1)/l_zeta_2

  D_1 <- duplication.matrix(L_1)
  D_1_plus <- solve(crossprod(D_1), t(D_1))

  eta_nu <- list(eta_in[[1]], eta_in[[2]])
  q_nu <- gauss_q(eta_nu, use_vech = FALSE)
  E_q_nu <- q_nu[[1]]
  Cov_q_nu <- q_nu[[2]]

  ent_nu <- d/2*(1 + log(2*pi)) + 1/2*determinant(Cov_q_nu, logarithm = TRUE)$modulus[1]

  nu_inds <- matrix(1:d, K + 2, L_1 + L_2 + 1)
  mu_inds <- nu_inds[, 1]
  psi_1_inds <- nu_inds[, 1:L_1 + 1]
  psi_2_inds <- nu_inds[, 1:L_2 + L_1 + 1]

  E_q_V <- matrix(E_q_nu, K + 2, L + 1)
  E_q_nu_mu <- E_q_V[, 1]
  E_q_V_psi_1 <- E_q_V[, 1:L_1 + 1]
  E_q_V_psi_2 <- E_q_V[, 1:L_2 + L_1 + 1]

  ent_zeta <- 0
  for(i in 1:N) {

    eta_q_zeta <- eta_in[[3]][, i] + eta_in[[4]][, i]
    eta_2 <- eta_q_zeta[(L_1 + 1):(L_1 + 0.5*L_1*(L_1 + 1))]
    Cov_q_zeta_1 <- -2*vecInverse(cprod(D_1_plus, eta_2))

    E_q_H_22 <- vector("list", length = M)
    for(j in 1:M) {

      E_q_H_22[[j]] <- matrix(NA, L_2, L_2)
      for(l_1 in 1:L_2) {

        inds_psi_21 <- psi_2_inds[, l_1]
        E_q_nu_psi_21 <- E_q_nu[inds_psi_21]

        for(l_2 in 1:L_2) {

          inds_psi_22 <- psi_2_inds[, l_2]
          E_q_nu_psi_22 <- E_q_nu[inds_psi_22]

          Cov_q_nu_psi_21 <- Cov_q_nu[inds_psi_22, inds_psi_21]

          tr_term <- tr(Cov_q_nu_psi_21 %*% crossprod(C[[i]][[j]]))
          cprod_term <- cprod(E_q_nu_psi_21, crossprod(C[[i]][[j]]) %*% E_q_nu_psi_22)
          E_q_H_22[[j]][l_1, l_2] <- tr_term + cprod_term
        }
      }
    }

    log_det_Cov_q_zeta_1 <- determinant(Cov_q_zeta_1, logarithm = TRUE)$modulus[1]
    log_det_E_q_H_22 <- sum(
      sapply(
        E_q_H_22,
        function(X) determinant(X, logarithm = TRUE)$modulus[1]
      )
    )
    log_det_Cov_q_zeta <- log_det_Cov_q_zeta_1 - log_det_E_q_H_22

    sum_val <- l_eta_zeta/2*(1 + log(2*pi)) + 1/2*log_det_Cov_q_zeta
    ent_zeta <- ent_zeta + sum_val
  }

  ans <- ent_nu + ent_zeta
  return(ans)
}

##########################################
#
#  CROSS-ENTROPY  FUNCTIONS
#
##########################################

cross_entropy_igw_prior <- function(eta_in, G_in, xi, Lambda) {

  q_X <- igw_q(eta_in, G_in)
  E_X_inv <- q_X[[3]]
  E_log_det_X <- q_X[[4]]

  G <- unique(G_in)

  l <- length(eta_in[[1]][-1])
  d <- (sqrt(1 + 8*l) - 1)/2

  if(G=="full") {

    term_1 <- d/2*(xi - d + 1)*log(2)
    term_2 <- d/4*(d-1)*log(pi)
    term_3 <- sum(lgamma((xi - d - c(1:d))/2 + 1))
    term_4 <- -(xi - d + 1)/2*determinant(Lambda, logarithm=TRUE)$modulus[1]
    term_5 <- (xi + 2)/2*E_log_det_X
    term_6 <- 1/2*tr(Lambda%*%E_X_inv)

    ans <- term_1 + term_2 + term_3 + term_4 + term_5 + term_6

  } else {

    term_1 <- d*xi/2*log(2)
    term_2 <- d*lgamma(xi/2)
    term_3 <- -xi/2*sum(log(diag(Lambda)))
    term_4 <- (xi + 2)/2*E_log_det_X
    term_5 <- 1/2*tr(Lambda%*%E_X_inv)

    ans <- term_1 + term_2 + term_3 + term_4 + term_5
  }

  return(ans)
}

cross_entropy_iter_igw <- function(eta_in, G_mess, xi, G_hyper) {

  # order of eta_in:
  # 1. Sigma -> p(Sigma|A)
  # 2. p(Sigma|A) -> Sigma
  # 3. A -> p(Sigma|A)
  # 4. p(Sigma|A) -> A

  # order of eta_out:
  # 1. p(Sigma|A) -> Sigma
  # 2. p(Sigma|A) -> A

  # order of G_out:
  # 1. p(Sigma|A) -> Sigma
  # 2. p(Sigma|A) -> A

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in)!=4) {

    stop("eta_in must have length 4")
  }

  if(!(is.character(G_mess) & is.character(G_hyper))) {

    stop("G_mess and G_hyper must be a characters")
  } else {

    if((G_mess != "full") & (G_mess != "diag")) {

      stop("G_mess can only be 'full' or 'diag'")
    }

    if((G_hyper != "full") & (G_hyper != "diag")) {

      stop("G_hyper can only be 'full' or 'diag'")
    }
  }

  if(xi<0) {

    stop("xi must be non-negative")
  }

  eta_in_Sigma <- list(eta_in[[1]], eta_in[[2]])
  G_Sigma <- rep(G_mess, 2)
  q_Sigma <- igw_q(eta_in_Sigma, G_Sigma)
  E_Sigma_inv <- q_Sigma[[3]]
  E_log_det_Sigma <- q_Sigma[[4]]

  eta_in_A <- list(eta_in[[3]], eta_in[[4]])
  G_A <- rep(G_hyper, 2)
  q_A <- igw_q(eta_in_A, G_A)
  E_A_inv <- q_A[[3]]
  E_log_det_A <- q_A[[4]]

  if(is.vector(E_Sigma_inv)) {

    d <- 1
  } else {

    d <- nrow(E_Sigma_inv)
  }

  term_1 <- d/2*(xi - d + 1)*log(2)
  term_2 <- d/4*(d - 1)*log(pi)
  term_3 <- sum(lgamma(1/2*(xi - d + c(1:d)) + 1))
  term_4 <- (xi - d + 1)/2*E_log_det_A
  term_5 <- (xi + 2)/2*E_log_det_Sigma
  term_6 <- 1/2*tr(E_A_inv%*%E_Sigma_inv)

  ans <- term_1 + term_2 + term_3 + term_4 + term_5 + term_6
  return(ans)
}

cross_entropy_gauss_prior <- function(eta_in, mu, Sigma, use_vech=TRUE) {

  q_x <- gauss_q(eta_in, use_vech)
  mu_q_x <- q_x[[1]]
  Sigma_q_x <- q_x[[2]]

  d <- length(mu)

  if(any(dim(Sigma)!=c(d,d))) {

    stop("Either Sigma or mu do not have the appropriate dimensions")
  }

  if(length(mu_q_x)!=d) {

    stop("mu_q_x and mu must have the same length")
  }

  if(any(dim(Sigma_q_x)!=dim(Sigma))) {

    stop("Sigma_q_x and Sigma must have the same dimension")
  }

  diff_mean <- mu_q_x - mu

  term_1 <- d/2*log(2*pi)
  term_2 <- 1/2*determinant(Sigma, logarithm=TRUE)$modulus[1]
  term_3 <- 1/2*tr(Sigma_q_x%*%solve(Sigma))
  term_4 <- 1/2*cprod(diff_mean, solve(Sigma, diff_mean))

  ans <- term_1 + term_2 + term_3 + term_4
  return(ans)
}

cross_entropy_logistic_lik <- function(eta_in, y, A, use_vech = TRUE) {

  q_x <- gauss_q(eta_in, use_vech)
  E_q_x <- q_x[[1]]
  Cov_q_x <- q_x[[2]]

  d <- length(E_q_x)
  n <- length(y)

  p <- c(
    0.003246343272134, 0.051517477033972, 0.195077912673858, 0.315569823632818,
    0.274149576158423, 0.131076880695470, 0.027912418727972, 0.001449567805354
  )

  s <- c(
    1.365340806296348, 1.059523971016916, 0.830791313765644, 0.650732166639391,
    0.508135425366489, 0.396313345166341, 0.308904252267995, 0.238212616409306
  )

  d <- length(E_q_x)
  D_d <- duplication.matrix(d)

  mu <- as.vector(A %*% E_q_x)
  sigsq <- diag(tcrossprod(A %*% Cov_q_x, A))

  Omega <- sqrt(matrix(1, n, 8) + tcrossprod(sigsq, s^2))
  B_0 <- as.vector(pnorm(tcrossprod(mu, s)/Omega) %*% p)

  ans <- cprod(y, A %*% E_q_x) - sum(B_0)
  return(ans)
}

cross_entropy_two_lev_gauss_prior <- function(eta_in, sigsq, p, q) {

  d <- length(eta_in[[1]])
  d_1 <- p + 0.5*p*(p + 1)
  d_2 <- q + 0.5*q*(q + 1) + p*q
  m <- (d - d_1)/d_2

  q_x <- two_level_nat_to_comm_parms(p, q, m, eta_in)
  E_q_x_1 <- q_x[[1]]
  Cov_q_x_1 <- q_x[[2]]
  E_q_x_2 <- q_x[[3]]
  Cov_q_x_2 <- q_x[[4]]

  tr_term <- 1/sigsq*tr(Cov_q_x_1)
  cprod_term <- 1/sigsq*cprod(E_q_x_1)
  for(j in 1:m) {

    sum_val <- 1/sigsq*tr(Cov_q_x_2[[j]])
    tr_term <- tr_term + sum_val

    sum_val <- 1/sigsq*cprod(E_q_x_2[[j]])
    cprod_term <- cprod_term + sum_val
  }

  term_1 <- d/2*log(2*pi)
  term_2 <- d/2*log(sigsq)
  term_3 <- 1/2*tr_term
  term_4 <- 1/2*cprod_term

  ans <- term_1 + term_2 + term_3 + term_4
  return(ans)
}

cross_entropy_fpc_gauss_pen <- function(eta_in, G_in, L, mu_beta, Sigma_beta) {

  # order of eta_in:
  # 1. nu -> p(nu|Sigma_nu)
  # 2. p(nu|Sigma_nu) -> nu
  # 3. sigsq_m -> p(nu|Sigma_nu)
  # 4. p(nu|Sigma_nu) -> sigsq_m
  # 5. sigsq_p -> p(nu|Sigma_nu)
  # 6. p(nu|Sigma_nu) -> sigsq_p

  # order of G_in:
  # 1. sigsq_m -> p(nu|Sigma_nu)
  # 2. p(nu|Sigma_nu) -> sigsq_m
  # 3. sigsq_p -> p(nu|Sigma_nu)
  # 4. p(nu|Sigma_nu) -> sigsq_p


  # order of eta_out:
  # 1. p(nu|Sigma_nu) -> nu
  # 2. p(nu|Sigma_nu) -> sigsq_m
  # 3. p(nu|Sigma_nu) -> sigsq_p

  # order of G_out:
  # 1. p(nu|Sigma_nu) -> sigsq_m
  # 2. p(nu|Sigma_nu) -> sigsq_p

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in)!=6) {

    stop("eta_in must have length 6")
  }

  if(!is.list(G_in)) {

    stop("G_in must be a list")
  }

  if(length(G_in)!=4) {

    stop("G_in must have length 4")
  }

  if(!all(Reduce("c", lapply(G_in, is.character)))) {

    stop("all of the elements of G_in must be characters")
  } else {

    if(!all(Reduce("c", lapply(G_in, function(x) x=="diag" | x=="full")))) {

      stop("The input graphs must be 'diag' or 'full'.")
    }
  }

  if(G_in[[1]]==G_in[[2]]) {

    G_sigsq_m <- G_in[[1]]
  } else {

    stop("The graph messages involving sigsq_m are not identical")
  }

  if(all(G_in[[3]]==G_in[[4]])) {

    G_sigsq_p <- G_in[[3]]
  } else {

    stop("The graph messages involving sigsq_p are not identical")
  }

  if(nrow(Sigma_beta)!=ncol(Sigma_beta)) {

    stop("Sigma_beta must be square")
  }

  if(any(Sigma_beta!=t(Sigma_beta))) {

    stop("Sigma_beta must be symmetric")
  }

  if(any(eigen(Sigma_beta)$values<0)) {

    stop("Sigma_beta must be positive-definite")
  }

  if(length(mu_beta)!=nrow(Sigma_beta)) {

    stop("length of mu_beta is not equal to the number of rows of Sigma_beta")
  }

  if(!is_int(L)) {

    stop("L must be a whole number.")
  } else {

    if(L<1) {

      stop("There must be at least one basis function.")
    }
  }

  eta_nu <- list(eta_in[[1]], eta_in[[2]])
  q_nu <- gauss_q(eta_nu, use_vech=FALSE)
  mu_q_nu <- q_nu[[1]]
  Sigma_q_nu <- q_nu[[2]]
  M_q_nu_nuT <- Sigma_q_nu + tcrossprod(mu_q_nu)

  eta_sigsq_m <- list(eta_in[[3]], eta_in[[4]])
  G_sigsq_m <- c(G_in[[1]], G_in[[2]])
  q_sigsq_m <- igw_q(eta_sigsq_m, G_sigsq_m)
  mu_q_recip_sigsq_m <- q_sigsq_m[[3]]
  mu_q_log_sigsq_m <- q_sigsq_m[[4]]

  mu_q_recip_sigsq_p <- rep(NA, L)
  mu_q_log_sigsq_p <- rep(NA, L)
  for(l in 1:L) {

    eta_sigsq_p <- list(eta_in[[5]][,l], eta_in[[6]][,l])
    G_sigsq_p <- c(G_in[[3]][l], G_in[[4]][l])
    q_sigsq_p <- igw_q(eta_sigsq_p, G_sigsq_p)
    mu_q_recip_sigsq_p[l] <- q_sigsq_p[[3]]
    mu_q_log_sigsq_p[l] <- q_sigsq_p[[4]]
  }

  l_eta_nu <- length(eta_in[[1]])
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/(L+1) - 2

  M_q_inv_Sigma_m <- blkdiag(solve(Sigma_beta), mu_q_recip_sigsq_m*diag(K))
  mu_m <- c(mu_beta, rep(0, K))

  M_q_inv_Sigma_p <- vector("list", length=L)
  mu_p <- vector("list", length=L)
  for(l in 1:L) {

    M_q_inv_Sigma_p[[l]] <- blkdiag(
      solve(Sigma_beta),
      mu_q_recip_sigsq_p[l]*diag(K)
    )

    mu_p[[l]] <- c(mu_beta, rep(0, K))
  }
  M_q_inv_Sigma_nu <- blkdiag(M_q_inv_Sigma_m, Reduce(blkdiag, M_q_inv_Sigma_p))
  mu_nu <- c(mu_m, Reduce("c", mu_p))

  term_1 <- d/2*log(2*pi)
  term_2 <- (L+1)/2*determinant(Sigma_beta, logarithm=TRUE)$modulus[1]
  term_3 <- K/2*mu_q_log_sigsq_m
  term_4 <- K/2*sum(mu_q_log_sigsq_p)
  term_5 <- tr(M_q_nu_nuT%*%M_q_inv_Sigma_nu)/2
  term_6 <- -cprod(mu_q_nu, M_q_inv_Sigma_nu%*%mu_nu)
  term_7 <- cprod(mu_nu, M_q_inv_Sigma_nu%*%mu_nu)/2

  ans <- term_1 + term_2 + term_3 + term_4 + term_5 + term_6 + term_7
  return(ans)
}

cross_entropy_fpc_lik_frag <- function(eta_in, G_in, C, Y, T_vec, L) {

  # order of eta_in:
  # 1. nu -> p(Y|nu,zeta,sigsq_eps)
  # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # 3. zeta -> p(Y|nu,zeta,sigsq_eps)
  # 4. p(Y|nu,zeta,sigsq_eps) -> zeta
  # 5. sigsq_eps -> p(Y|nu,zeta,sigsq_eps)
  # 6. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of G_in:
  # 1. sigsq_eps -> p(Y|nu,zeta,sigsq_eps)
  # 2. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of eta_out:
  # 1. p(Y|nu,zeta,sigsq_eps) -> nu
  # 2. p(Y|nu,zeta,sigsq_eps) -> zeta
  # 3. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of G_out:
  # 1. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in)!=6) {

    stop("eta_in must have length 6")
  }

  if(!is.list(G_in)) {

    stop("G_in must be a list")
  }

  if(length(G_in)!=2) {

    stop("G_in must have length 2")
  }

  if(!all(Reduce("c", lapply(G_in, is.character)))) {

    stop("all of the elements of G_in must be characters")
  } else {

    if(!all(Reduce("c", lapply(G_in, function(x) x=="diag" || x=="full")))) {

      stop("The input graphs must be 'diag' or 'full'.")
    }
  }

  if(G_in[[1]]==G_in[[2]]) {

    G_sigsq_eps <- G_in[[1]]
  } else {

    stop("The graph messages involving sigsq_eps are not identical")
  }

  if(!is_int(L)) {

    stop("L must be a whole number.")
  } else {

    if(L<1) {

      stop("There must be at least one basis function.")
    }
  }

  N <- ncol(eta_in[[3]])
  l_eta_nu <- length(eta_in[[1]])
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/(L+1) - 2
  sum_T <- sum(T_vec)

  eta_in_nu <- list(eta_in[[1]], eta_in[[2]])
  q_nu <- gauss_q(eta_in_nu, use_vech=FALSE)
  mu_q_nu <- q_nu[[1]]
  Sigma_q_nu <- q_nu[[2]]

  eta_in_sigsq_eps <- list(eta_in[[5]], eta_in[[6]])
  G_in_sigsq_eps <- c(G_in[[1]], G_in[[2]])
  q_sigsq_eps <- igw_q(eta_in_sigsq_eps, G_in_sigsq_eps)
  mu_q_recip_sigsq_eps <- q_sigsq_eps[[3]]
  mu_q_log_sigsq_eps <- q_sigsq_eps[[4]]

  mu_q_zeta <- vector("list", length=N)
  Sigma_q_zeta <- vector("list", length=N)
  for(i in 1:N) {

    eta_in_zeta <- list(
      eta_in[[3]][,i],
      eta_in[[4]][,i]
    )
    q_zeta <- gauss_q(eta_in_zeta, use_vech=TRUE)
    mu_q_zeta[[i]] <- q_zeta[[1]]
    Sigma_q_zeta[[i]] <- q_zeta[[2]]
  }

  mu_inds <- 1:(K+2)
  mu_q_nu_mu <- mu_q_nu[mu_inds]
  Sigma_q_nu_mu <- Sigma_q_nu[mu_inds, mu_inds]

  psi_inds <- (1:d)[-mu_inds]
  psi_groups <- ceiling(psi_inds/(K+2))-1
  psi_inds <- split(psi_inds, psi_groups)

  M_q_H <- vector("list", length=N)
  mu_q_h <- vector("list", length=N)
  M_q_H_tilde <- vector("list", length=N)
  mu_q_h_mu <- rep(NA, N)
  for(i in 1:N) {

    tr_term <- tr(Sigma_q_nu_mu%*%crossprod(C[[i]]))
    cprod_term <- cprod(mu_q_nu_mu, crossprod(C[[i]])%*%mu_q_nu_mu)
    mu_q_h_mu[i] <- tr_term + cprod_term

    M_q_H[[i]] <- matrix(NA, L, L)
    mu_q_h[[i]] <- rep(NA, L)

    for(l_i in 1:L) {

      i_inds <- psi_inds[[l_i]]
      mu_q_nu_psi_i <- mu_q_nu[i_inds]

      Cov_q_nu_psi_i <- Sigma_q_nu[mu_inds, i_inds]

      tr_term <- tr(Cov_q_nu_psi_i%*%crossprod(C[[i]]))
      cprod_term <- cprod(mu_q_nu_psi_i, crossprod(C[[i]])%*%mu_q_nu_mu)
      mu_q_h[[i]][l_i] <- tr_term + cprod_term

      for(l_j in 1:L) {

        j_inds <- psi_inds[[l_j]]
        mu_q_nu_psi_j <- mu_q_nu[j_inds]

        Cov_q_nu_psi_ji <- Sigma_q_nu[j_inds, i_inds]

        tr_term <- tr(Cov_q_nu_psi_ji%*%crossprod(C[[i]]))
        cprod_term <- cprod(mu_q_nu_psi_i, crossprod(C[[i]])%*%mu_q_nu_psi_j)
        M_q_H[[i]][l_i, l_j] <- tr_term + cprod_term
      }
    }

    top_block <- t(c(mu_q_h_mu[i], mu_q_h[[i]]))
    bottom_block <- cbind(mu_q_h[[i]], M_q_H[[i]])
    M_q_H_tilde[[i]] <- rbind(top_block, bottom_block)
  }

  M_q_V_psi <- matrix(NA, K+2, L)
  for(l in 1:L) {

    psi_l_inds <- psi_inds[[l]]
    M_q_V_psi[,l] <- mu_q_nu[psi_l_inds]
  }

  M_q_V <- cbind(mu_q_nu_mu, M_q_V_psi)

  E_q_resid <- 0
  for(i in 1:N) {

    mu_q_zeta_tilde <- c(1, mu_q_zeta[[i]])
    Sigma_q_zeta_tilde <- blkdiag(matrix(0), Sigma_q_zeta[[i]])
    M_q_zeta_zeta_T_tilde <- Sigma_q_zeta_tilde + tcrossprod(mu_q_zeta_tilde)

    E_q_Y_hat <- C[[i]]%*%M_q_V%*%mu_q_zeta_tilde

    term_1 <- cprod(Y[[i]])
    term_2 <- -2*cprod(E_q_Y_hat, Y[[i]])
    term_3 <- tr(M_q_zeta_zeta_T_tilde%*%M_q_H_tilde[[i]])
    sum_val <- term_1 + term_2 + term_3

    E_q_resid <- E_q_resid + sum_val
  }

  term_1 <- sum_T/2*log(2*pi)
  term_2 <- sum_T/2*mu_q_log_sigsq_eps
  term_3 <- 1/2*mu_q_recip_sigsq_eps*E_q_resid

  ans <- term_1 + term_2 + term_3
  return(ans)
}

cross_entropy_mfpc_lik_frag <- function(eta_in, G_in, C, Y, L) {

  # order of eta_in:
  # 1. nu -> p(Y|nu,zeta,sigsq_eps)
  # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # 3. zeta -> p(Y|nu,zeta,sigsq_eps)
  # 4. p(Y|nu,zeta,sigsq_eps) -> zeta
  # 5. sigsq_eps -> p(Y|nu,zeta,sigsq_eps)
  # 6. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of G_in:
  # 1. sigsq_eps -> p(Y|nu,zeta,sigsq_eps)
  # 2. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of eta_out:
  # 1. p(Y|nu,zeta,sigsq_eps) -> nu
  # 2. p(Y|nu,zeta,sigsq_eps) -> zeta
  # 3. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of G_out:
  # 1. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in) != 6) {

    stop("eta_in must have length 6")
  }

  if(!is.list(G_in)) {

    stop("G_in must be a list")
  }

  if(length(G_in) != 2) {

    stop("G_in must have length 2")
  }

  if(!all(sapply(unlist(G_in), is.character))) {

    stop("all of the elements of G_in must be characters")
  } else {

    if(!all(sapply(unlist(G_in), function(x) x=="diag" || x=="full"))) {

      stop("The input graphs must be 'diag' or 'full'.")
    }
  }

  if(all(unlist(G_in[[1]]) == unlist(G_in[[2]]))) {

    G_sigsq_eps <- G_in[[1]]
  } else {

    stop("The graph messages involving sigsq_eps are not identical")
  }

  if(!is_int(L)) {

    stop("L must be a whole number.")
  } else {

    if(L < 1) {

      stop("There must be at least one basis function for the first level.")
    }
  }

  N <- ncol(eta_in[[3]])
  p <- length(eta_in[[1]])
  l_eta_nu <- length(eta_in[[1]][[1]])
  l_eta_zeta <- L + 0.5*L*(L + 1)
  K <- sapply(C[[1]], ncol) - 2
  d <- (K+2)*(L+1)
  D_L <- duplication.matrix(L)

  YTY <- matrix(NA, N, p)
  CTY <- vector("list", length = N)
  CTC <- vector("list", length = N)
  for(i in 1:N) {

    CTY[[i]] <- vector("list", length = p)
    CTC[[i]] <- vector("list", length = p)
    for(j in 1:p) {

      YTY[i, j] <- cprod(Y[[i]][[j]])
      CTY[[i]][[j]] <- cprod(Y[[i]][[j]], C[[i]][[j]])
      CTC[[i]][[j]] <- crossprod(C[[i]][[j]])
    }
  }

  E_q_zeta <- vector("list", length = N)
  Cov_q_zeta <- vector("list", length = N)
  E_q_tcross_zeta <- vector("list", length = N)
  for(i in 1:N) {

    eta_in_zeta <- list(
      eta_in[[3]][, i],
      eta_in[[4]][, i]
    )
    q_zeta <- gauss_q(eta_in_zeta, use_vech = TRUE)
    E_q_zeta[[i]] <- q_zeta[[1]]
    Cov_q_zeta[[i]] <- q_zeta[[2]]
    E_q_tcross_zeta[[i]] <- Cov_q_zeta[[i]] + tcrossprod(E_q_zeta[[i]])
  }

  E_q_recip_sigsq_eps <- rep(NA, p)
  E_q_log_sigsq_eps <- rep(NA, p)
  for(j in 1:p) {

    eta_in_sigsq_eps <- list(eta_in[[5]][[j]], eta_in[[6]][[j]])
    G_in_sigsq_eps <- c(G_in[[1]][[j]], G_in[[2]][[j]])
    q_sigsq_eps <- igw_q(eta_in_sigsq_eps, G_in_sigsq_eps)
    E_q_recip_sigsq_eps[j] <- q_sigsq_eps[[3]]
    E_q_log_sigsq_eps[j] <- q_sigsq_eps[[4]]
  }

  E_q_nu <- vector("list", length = p)
  Cov_q_nu <- vector("list", length = p)
  for(j in 1:p) {

    eta_in_nu <- list(eta_in[[1]][[j]], eta_in[[2]][[j]])
    q_nu <- gauss_q(eta_in_nu, use_vech = FALSE)
    E_q_nu[[j]] <- q_nu[[1]]
    Cov_q_nu[[j]] <- q_nu[[2]]
  }

  E_q_V <- vector("list", length = p)
  E_q_nu_mu <- vector("list", length = p)
  E_q_V_psi <- vector("list", length = p)
  E_q_nu_psi <- vector("list", length = p)
  Cov_q_nu_mu <- vector("list", length = p)
  Cov_q_nu_psi <- vector("list", length = p)
  Cov_q_nu_mu_psi <- vector("list", length = p)
  for(j in 1:p) {

    inds_mat <- matrix(1:d[j], K[j] + 2, L + 1)
    mu_inds <- inds_mat[, 1]
    psi_inds <- as.vector(inds_mat[, 1:L + 1])

    E_q_V[[j]] <- matrix(E_q_nu[[j]], K[j] + 2, L + 1)
    E_q_nu_mu[[j]] <- E_q_V[[j]][, 1]
    E_q_V_psi[[j]] <- E_q_V[[j]][, 1:L + 1]
    E_q_nu_psi[[j]] <- as.vector(E_q_V_psi[[j]])

    Cov_q_nu_mu[[j]] <- Cov_q_nu[[j]][mu_inds, mu_inds]
    Cov_q_nu_psi[[j]] <- Cov_q_nu[[j]][psi_inds, psi_inds]
    Cov_q_nu_mu_psi[[j]] <- Cov_q_nu[[j]][mu_inds, psi_inds]
  }

  E_q_h_mu <- matrix(NA, N, p)
  E_q_h_mu_psi <- vector("list", length = N)
  E_q_H_psi <- vector("list", length = N)
  for(i in 1:N) {

    E_q_h_mu_psi[[i]] <- matrix(NA, L, p)
    E_q_H_psi[[i]] <- vector("list", length = p)
    for(j in 1:p) {

      E_q_h_mu[i, j] <- E_cprod(
        E_q_nu_mu[[j]], Cov_q_nu_mu[[j]],
        E_q_nu_mu[[j]], CTC[[i]][[j]]
      )

      E_q_h_mu_psi[[i]][, j] <- E_h(
        L, E_q_nu_psi[[j]], Cov_q_nu_mu_psi[[j]],
        E_q_nu_mu[[j]], CTC[[i]][[j]]
      )

      E_q_H_psi[[i]][[j]] <- E_H(
        L, L, E_q_nu_psi[[j]], Cov_q_nu_psi[[j]],
        E_q_nu_psi[[j]], CTC[[i]][[j]]
      )
    }
  }

  E_q_cprod_resids <- rep(0, p)
  for(j in 1:p) {

    for(i in 1:N) {

      summands <- rep(NA, 6)
      summands[1] <- YTY[i, j]
      summands[2] <- - 2*cprod(E_q_nu_mu[[j]], CTY[[i]][[j]])
      summands[3] <- -2*cprod(E_q_V_psi[[j]] %*% E_q_zeta[[i]], CTY[[i]][[j]])
      summands[4] <- E_q_h_mu[i, j]
      summands[5] <- 2*cprod(E_q_zeta[[i]], E_q_h_mu_psi[[i]][, j])
      summands[6] <- tr(E_q_tcross_zeta[[i]] %*% E_q_H_psi[[i]][[j]])
      sum_val <- sum(summands)
      E_q_cprod_resids[j] <- E_q_cprod_resids[j] + sum_val
    }
  }

  summands <- rep(NA, 3)
  summands[1] <- -0.5*log(2*pi)*sum(n)
  summands[2] <- -0.5*sum(n %*% E_q_log_sigsq_eps)
  summands[3] <- -0.5*cprod(E_q_recip_sigsq_eps, E_q_cprod_resids)

  ans <- -sum(summands)
  return(ans)
}

cross_entropy_mlfpc_lik_frag <- function(eta_in, G_in, C, Y, T_vec, L_1, L_2) {

  # order of eta_in:
  # 1. nu -> p(Y|nu,zeta,sigsq_eps)
  # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # 3. zeta -> p(Y|nu,zeta,sigsq_eps)
  # 4. p(Y|nu,zeta,sigsq_eps) -> zeta
  # 5. sigsq_eps -> p(Y|nu,zeta,sigsq_eps)
  # 6. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of G_in:
  # 1. sigsq_eps -> p(Y|nu,zeta,sigsq_eps)
  # 2. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of eta_out:
  # 1. p(Y|nu,zeta,sigsq_eps) -> nu
  # 2. p(Y|nu,zeta,sigsq_eps) -> zeta
  # 3. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # order of G_out:
  # 1. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in)!=6) {

    stop("eta_in must have length 6")
  }

  if(!is.list(G_in)) {

    stop("G_in must be a list")
  }

  if(length(G_in)!=2) {

    stop("G_in must have length 2")
  }

  if(!all(Reduce("c", lapply(G_in, is.character)))) {

    stop("all of the elements of G_in must be characters")
  } else {

    if(!all(Reduce("c", lapply(G_in, function(x) x=="diag" || x=="full")))) {

      stop("The input graphs must be 'diag' or 'full'.")
    }
  }

  if(G_in[[1]]==G_in[[2]]) {

    G_sigsq_eps <- G_in[[1]]
  } else {

    stop("The graph messages involving sigsq_eps are not identical")
  }

  if(!is_int(L_1)) {

    stop("L_1 must be a whole number.")
  } else {

    if(L_1 < 1) {

      stop("There must be at least one basis function for the first level.")
    }
  }

  if(!is_int(L_2)) {

    stop("L_2 must be a whole number.")
  } else {

    if(L_2 < 1) {

      stop("There must be at least one basis function for the second level.")
    }
  }

  L <- L_1 + L_2
  N <- length(eta_in[[3]])
  l_eta_nu <- length(eta_in[[1]])
  l_eta_zeta <- sapply(eta_in[[3]], length)
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/(L+1) - 2
  l_zeta_1 <- L_1 + 0.5*L_1*(L_1 + 1)
  l_zeta_2 <- L_2 + 0.5*L_2*(L_2 + 1) + L_1*L_2
  M <- (l_eta_zeta - l_zeta_1)/l_zeta_2
  sum_T <- sum(sapply(T_vec, sum))

  YTY <- vector("list", length = N)
  CTY <- vector("list", length = N)
  CTC <- vector("list", length = N)
  for(i in 1:N) {

    YTY[[i]] <- rep(NA, M[i])
    CTY[[i]] <- matrix(NA, K + 2, M[i])
    CTC[[i]] <- vector("list", length = M[i])
    for(j in 1:M[i]) {

      YTY[[i]][j] <- cprod(Y[[i]][[j]])
      CTY[[i]][, j] <- cprod(C[[i]][[j]], Y[[i]][[j]])
      CTC[[i]][[j]] <- crossprod(C[[i]][[j]])
    }
  }

  eta_in_nu <- list(eta_in[[1]], eta_in[[2]])
  q_nu <- gauss_q(eta_in_nu, use_vech = FALSE)
  E_q_nu <- q_nu[[1]]
  Cov_q_nu <- q_nu[[2]]

  eta_in_sigsq_eps <- list(eta_in[[5]], eta_in[[6]])
  G_in_sigsq_eps <- c(G_in[[1]], G_in[[2]])
  q_sigsq_eps <- igw_q(eta_in_sigsq_eps, G_in)
  E_q_recip_sigsq_eps <- q_sigsq_eps[[3]]
  E_q_log_sigsq_eps <- q_sigsq_eps[[4]]

  E_q_zeta_1 <- vector("list", length = N)
  Cov_q_zeta_1 <- vector("list", length = N)
  E_q_zeta_2 <- vector("list", length = N)
  Cov_q_zeta_2 <- vector("list", length = N)
  Cov_q_zeta_12 <- vector("list", length = N)
  for(i in 1:N) {

    eta_in_zeta <- list(
      eta_in[[3]][[i]],
      eta_in[[4]][[i]]
    )
    q_zeta <- two_level_nat_to_comm_parms(L_1, L_2, M[i], eta_in_zeta)
    E_q_zeta_1[[i]] <- q_zeta[[1]]
    Cov_q_zeta_1[[i]] <- q_zeta[[2]]
    E_q_zeta_2[[i]] <- q_zeta[[3]]
    Cov_q_zeta_2[[i]] <- q_zeta[[4]]
    Cov_q_zeta_12[[i]] <- q_zeta[[5]]
  }

  inds_mat <- matrix(1:d, K + 2, L + 1)
  mu_inds <- inds_mat[, 1]
  psi_1_inds <- inds_mat[, 1:L_1 + 1]
  psi_2_inds <- inds_mat[, 1:L_2 + L_1 + 1]

  E_q_V <- matrix(E_q_nu, K+2, L + 1)
  E_q_nu_mu <- E_q_V[, 1]
  E_q_V_psi_1 <- E_q_V[, 1:L_1 + 1]
  E_q_nu_psi_1 <- as.vector(E_q_V_psi_1)
  E_q_V_psi_2 <- E_q_V[, 1:L_2 + L_1 + 1]
  E_q_nu_psi_2 <- as.vector(E_q_V_psi_2)

  Cov_q_nu_mu <- Cov_q_nu[mu_inds, mu_inds]
  Cov_q_nu_psi_1 <- Cov_q_nu[as.vector(psi_1_inds), as.vector(psi_1_inds)]
  Cov_q_nu_psi_2 <- Cov_q_nu[as.vector(psi_2_inds), as.vector(psi_2_inds)]
  Cov_q_nu_mu_psi_1 <- Cov_q_nu[mu_inds, as.vector(psi_1_inds)]
  Cov_q_nu_mu_psi_2 <- Cov_q_nu[mu_inds, as.vector(psi_2_inds)]
  Cov_q_nu_psi_21 <- Cov_q_nu[as.vector(psi_2_inds), as.vector(psi_1_inds)]

  E_q_h_mu <- E_cprod(
    E_q_nu_mu, Cov_q_nu_mu, E_q_nu_mu,
    Reduce("+", lapply(CTC, Reduce, f = "+"))
  )

  E_q_h_mu_psi_1 <- vector("list", length = N)
  E_q_H_psi_11 <- vector("list", length = N)
  E_q_h_mu_psi_2 <- vector("list", length = N)
  E_q_H_psi_22 <- vector("list", length = N)
  E_q_H_psi_12 <- vector("list", length = N)
  for(i in 1:N) {

    CTY_i <- rowSums(CTY[[i]])
    CTC_i <- Reduce("+", CTC[[i]])

    E_q_h_mu_psi_1[[i]] <- E_h(L_1, E_q_nu_psi_1, Cov_q_nu_mu_psi_1, E_q_nu_mu, CTC_i)
    E_q_H_psi_11[[i]] <- E_H(L_1, L_1, E_q_nu_psi_1, Cov_q_nu_psi_1, E_q_nu_psi_1, CTC_i)

    E_q_h_mu_psi_2[[i]] <- vector("list", length = M[i])
    E_q_H_psi_22[[i]] <- vector("list", length = M[i])
    E_q_H_psi_12[[i]] <- vector("list", length = M[i])
    for(j in 1:M[i]) {

      E_q_h_mu_psi_2[[i]][[j]] <- E_h(
        L_2, E_q_nu_psi_2, Cov_q_nu_mu_psi_2,
        E_q_nu_mu, CTC[[i]][[j]]
      )

      E_q_H_psi_22[[i]][[j]] <- E_H(
        L_2, L_2, E_q_nu_psi_2, Cov_q_nu_psi_2,
        E_q_nu_psi_2, CTC[[i]][[j]]
      )

      E_q_H_psi_12[[i]][[j]] <- E_H(
        L_1, L_2, E_q_nu_psi_1, Cov_q_nu_psi_21,
        E_q_nu_psi_2, CTC[[i]][[j]]
      )
    }
  }

  term_1 <- sum(sapply(YTY, sum))
  term_2 <- -2*cprod(E_q_nu_mu, Reduce("+", lapply(CTY, rowSums)))
  term_3 <- 0
  term_4 <- 0
  for(i in 1:N) {

    sum_val <- -2*cprod(E_q_V_psi_1 %*% E_q_zeta_1[[i]], rowSums(CTY[[i]]))
    term_3 <- term_3 + sum_val

    for(j in 1:M[i]) {

      sum_val <- -2*cprod(E_q_V_psi_2 %*% E_q_zeta_2[[i]][[j]], CTY[[i]][, j])
      term_4 <- term_4 + sum_val
    }
  }
  y_line <- term_1 + term_2 + term_3 + term_4

  term_1 <- E_q_h_mu
  term_2 <- 0
  term_3 <- 0
  for(i in 1:N) {

    sum_val <- 2*cprod(E_q_zeta_1[[i]], E_q_h_mu_psi_1[[i]])
    term_2 <- term_2 + sum_val

    for(j in 1:M[i]) {

      sum_val <- 2*cprod(E_q_zeta_2[[i]][[j]], E_q_h_mu_psi_2[[i]][[j]])
      term_3 <- term_3 + sum_val
    }
  }
  mu_line <- term_1 + term_2 + term_3

  term_1 <- 0
  term_2 <- 0
  for(i in 1:N) {

    E_vec <- E_q_zeta_1[[i]]
    Cov_mat <- Cov_q_zeta_1[[i]]
    A_mat <- E_q_H_psi_11[[i]]
    sum_val <- E_cprod(E_vec, Cov_mat, E_vec, A_mat)
    term_1 <- term_1 + sum_val

    for(j in 1:M[i]) {

      E_vec_1 <- E_q_zeta_1[[i]]
      Cov_mat <- t(Cov_q_zeta_12[[i]][[j]])
      E_vec_2 <- E_q_zeta_2[[i]][[j]]
      A_mat <- E_q_H_psi_12[[i]][[j]]
      sum_val <- 2*E_cprod(E_vec_1, Cov_mat, E_vec_2, A_mat)
      term_2 <- term_2 + sum_val
    }
  }
  zeta_1_line <- term_1 + term_2

  term_1 <- 0
  for(i in 1:N) {

    for(j in 1:M[i]) {

      E_vec <- E_q_zeta_2[[i]][[j]]
      Cov_mat <- Cov_q_zeta_2[[i]][[j]]
      A_mat <- E_q_H_psi_22[[i]][[j]]
      sum_val <- E_cprod(E_vec, Cov_mat, E_vec, A_mat)
      term_1 <- term_1 + sum_val
    }
  }
  zeta_2_line <- term_1

  E_q_norm <- y_line + mu_line + zeta_1_line + zeta_2_line

  term_1 <- sum_T/2*log(2*pi)
  term_2 <- sum_T/2*E_q_log_sigsq_eps
  term_3 <- 1/2*E_q_recip_sigsq_eps*E_q_norm

  ans <- term_1 + term_2 + term_3
  return(ans)
}

# cross_entropy_mlfpc_lik_frag <- function(eta_in, G_in, C, Y, T_vec, L_1, L_2) {

  # # order of eta_in:
  # # 1. nu -> p(Y|nu,zeta,sigsq_eps)
  # # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # # 3. zeta -> p(Y|nu,zeta,sigsq_eps)
  # # 4. p(Y|nu,zeta,sigsq_eps) -> zeta
  # # 5. sigsq_eps -> p(Y|nu,zeta,sigsq_eps)
  # # 6. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # # order of G_in:
  # # 1. sigsq_eps -> p(Y|nu,zeta,sigsq_eps)
  # # 2. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # # order of eta_out:
  # # 1. p(Y|nu,zeta,sigsq_eps) -> nu
  # # 2. p(Y|nu,zeta,sigsq_eps) -> zeta
  # # 3. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # # order of G_out:
  # # 1. p(Y|nu,zeta,sigsq_eps) -> sigsq_eps

  # if(!is.list(eta_in)) {

    # stop("eta_in must be a list")
  # }

  # if(length(eta_in)!=6) {

    # stop("eta_in must have length 6")
  # }

  # if(!is.list(G_in)) {

    # stop("G_in must be a list")
  # }

  # if(length(G_in)!=2) {

    # stop("G_in must have length 2")
  # }

  # if(!all(Reduce("c", lapply(G_in, is.character)))) {

    # stop("all of the elements of G_in must be characters")
  # } else {

    # if(!all(Reduce("c", lapply(G_in, function(x) x=="diag" || x=="full")))) {

      # stop("The input graphs must be 'diag' or 'full'.")
    # }
  # }

  # if(G_in[[1]]==G_in[[2]]) {

    # G_sigsq_eps <- G_in[[1]]
  # } else {

    # stop("The graph messages involving sigsq_eps are not identical")
  # }

  # if(!is_int(L_1)) {

    # stop("L_1 must be a whole number.")
  # } else {

    # if(L_1 < 1) {

      # stop("There must be at least one basis function for the first level.")
    # }
  # }

  # if(!is_int(L_2)) {

    # stop("L_2 must be a whole number.")
  # } else {

    # if(L_2 < 1) {

      # stop("There must be at least one basis function for the second level.")
    # }
  # }

  # L <- L_1 + L_2
  # N <- ncol(eta_in[[3]])
  # l_eta_nu <- length(eta_in[[1]])
  # l_eta_zeta <- nrow(eta_in[[3]])
  # d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  # K <- d/(L+1) - 2
  # l_zeta_1 <- L_1 + 0.5*L_1*(L_1 + 1)
  # l_zeta_2 <- L_2 + 0.5*L_2*(L_2 + 1) + L_1*L_2
  # M <- (l_eta_zeta - l_zeta_1)/l_zeta_2
  # sum_T <- sum(sapply(T_vec, sum))

  # YTY <- matrix(NA, N, M)
  # CTY <- vector("list", length = N)
  # CTC <- vector("list", length = N)
  # for(i in 1:N) {

    # CTY[[i]] <- matrix(NA, K + 2, M)
    # CTC[[i]] <- vector("list", length = M)
    # for(j in 1:M) {

      # YTY[i, j] <- cprod(Y[[i]][[j]])
      # CTY[[i]][, j] <- cprod(C[[i]][[j]], Y[[i]][[j]])
      # CTC[[i]][[j]] <- crossprod(C[[i]][[j]])
    # }
  # }

  # eta_in_nu <- list(eta_in[[1]], eta_in[[2]])
  # q_nu <- gauss_q(eta_in_nu, use_vech = FALSE)
  # E_q_nu <- q_nu[[1]]
  # Cov_q_nu <- q_nu[[2]]

  # eta_in_sigsq_eps <- list(eta_in[[5]], eta_in[[6]])
  # G_in_sigsq_eps <- c(G_in[[1]], G_in[[2]])
  # q_sigsq_eps <- igw_q(eta_in_sigsq_eps, G_in)
  # E_q_recip_sigsq_eps <- q_sigsq_eps[[3]]
  # E_q_log_sigsq_eps <- q_sigsq_eps[[4]]

  # E_q_zeta_1 <- vector("list", length = N)
  # Cov_q_zeta_1 <- vector("list", length = N)
  # E_q_zeta_2 <- vector("list", length = N)
  # Cov_q_zeta_2 <- vector("list", length = N)
  # Cov_q_zeta_12 <- vector("list", length = N)
  # for(i in 1:N) {

    # eta_in_zeta <- list(
      # eta_in[[3]][, i],
      # eta_in[[4]][, i]
    # )
    # q_zeta <- two_level_nat_to_comm_parms(L_1, L_2, M, eta_in_zeta)
    # E_q_zeta_1[[i]] <- q_zeta[[1]]
    # Cov_q_zeta_1[[i]] <- q_zeta[[2]]
    # E_q_zeta_2[[i]] <- q_zeta[[3]]
    # Cov_q_zeta_2[[i]] <- q_zeta[[4]]
    # Cov_q_zeta_12[[i]] <- q_zeta[[5]]
  # }

  # inds_mat <- matrix(1:d, K + 2, L + 1)
  # mu_inds <- inds_mat[, 1]
  # psi_1_inds <- inds_mat[, 1:L_1 + 1]
  # psi_2_inds <- inds_mat[, 1:L_2 + L_1 + 1]

  # E_q_V <- matrix(E_q_nu, K+2, L + 1)
  # E_q_nu_mu <- E_q_V[, 1]
  # E_q_V_psi_1 <- E_q_V[, 1:L_1 + 1]
  # E_q_nu_psi_1 <- as.vector(E_q_V_psi_1)
  # E_q_V_psi_2 <- E_q_V[, 1:L_2 + L_1 + 1]
  # E_q_nu_psi_2 <- as.vector(E_q_V_psi_2)

  # Cov_q_nu_mu <- Cov_q_nu[mu_inds, mu_inds]
  # Cov_q_nu_psi_1 <- Cov_q_nu[as.vector(psi_1_inds), as.vector(psi_1_inds)]
  # Cov_q_nu_psi_2 <- Cov_q_nu[as.vector(psi_2_inds), as.vector(psi_2_inds)]
  # Cov_q_nu_mu_psi_1 <- Cov_q_nu[mu_inds, as.vector(psi_1_inds)]
  # Cov_q_nu_mu_psi_2 <- Cov_q_nu[mu_inds, as.vector(psi_2_inds)]
  # Cov_q_nu_psi_21 <- Cov_q_nu[as.vector(psi_2_inds), as.vector(psi_1_inds)]

  # E_q_h_mu <- E_cprod(
    # E_q_nu_mu, Cov_q_nu_mu, E_q_nu_mu,
    # Reduce("+", lapply(CTC, Reduce, f = "+"))
  # )

  # E_q_h_mu_psi_1 <- vector("list", length = N)
  # E_q_H_psi_11 <- vector("list", length = N)
  # E_q_h_mu_psi_2 <- vector("list", length = N)
  # E_q_H_psi_22 <- vector("list", length = N)
  # E_q_H_psi_12 <- vector("list", length = N)
  # for(i in 1:N) {

    # CTY_i <- rowSums(CTY[[i]])
    # CTC_i <- Reduce("+", CTC[[i]])

    # E_q_h_mu_psi_1[[i]] <- E_h(L_1, E_q_nu_psi_1, Cov_q_nu_mu_psi_1, E_q_nu_mu, CTC_i)
    # E_q_H_psi_11[[i]] <- E_H(L_1, L_1, E_q_nu_psi_1, Cov_q_nu_psi_1, E_q_nu_psi_1, CTC_i)

    # E_q_h_mu_psi_2[[i]] <- vector("list", length = M)
    # E_q_H_psi_22[[i]] <- vector("list", length = M)
    # E_q_H_psi_12[[i]] <- vector("list", length = M)
    # for(j in 1:M) {

      # E_q_h_mu_psi_2[[i]][[j]] <- E_h(
        # L_2, E_q_nu_psi_2, Cov_q_nu_mu_psi_2,
        # E_q_nu_mu, CTC[[i]][[j]]
      # )

      # E_q_H_psi_22[[i]][[j]] <- E_H(
        # L_2, L_2, E_q_nu_psi_2, Cov_q_nu_psi_2,
        # E_q_nu_psi_2, CTC[[i]][[j]]
      # )

      # E_q_H_psi_12[[i]][[j]] <- E_H(
        # L_1, L_2, E_q_nu_psi_1, Cov_q_nu_psi_21,
        # E_q_nu_psi_2, CTC[[i]][[j]]
      # )
    # }
  # }

  # term_1 <- sum(YTY)
  # term_2 <- -2*cprod(E_q_nu_mu, Reduce("+", lapply(CTY, rowSums)))
  # term_3 <- 0
  # term_4 <- 0
  # for(i in 1:N) {

    # sum_val <- -2*cprod(E_q_V_psi_1 %*% E_q_zeta_1[[i]], rowSums(CTY[[i]]))
    # term_3 <- term_3 + sum_val

    # for(j in 1:M) {

      # sum_val <- -2*cprod(E_q_V_psi_2 %*% E_q_zeta_2[[i]][[j]], CTY[[i]][, j])
      # term_4 <- term_4 + sum_val
    # }
  # }
  # y_line <- term_1 + term_2 + term_3 + term_4

  # term_1 <- E_q_h_mu
  # term_2 <- 0
  # term_3 <- 0
  # for(i in 1:N) {

    # sum_val <- 2*cprod(E_q_zeta_1[[i]], E_q_h_mu_psi_1[[i]])
    # term_2 <- term_2 + sum_val

    # for(j in 1:M) {

      # sum_val <- 2*cprod(E_q_zeta_2[[i]][[j]], E_q_h_mu_psi_2[[i]][[j]])
      # term_3 <- term_3 + sum_val
    # }
  # }
  # mu_line <- term_1 + term_2 + term_3

  # term_1 <- 0
  # term_2 <- 0
  # for(i in 1:N) {

    # E_vec <- E_q_zeta_1[[i]]
    # Cov_mat <- Cov_q_zeta_1[[i]]
    # A_mat <- E_q_H_psi_11[[i]]
    # sum_val <- E_cprod(E_vec, Cov_mat, E_vec, A_mat)
    # term_1 <- term_1 + sum_val

    # for(j in 1:M) {

      # E_vec_1 <- E_q_zeta_1[[i]]
      # Cov_mat <- t(Cov_q_zeta_12[[i]][[j]])
      # E_vec_2 <- E_q_zeta_2[[i]][[j]]
      # A_mat <- E_q_H_psi_12[[i]][[j]]
      # sum_val <- 2*E_cprod(E_vec_1, Cov_mat, E_vec_2, A_mat)
      # term_2 <- term_2 + sum_val
    # }
  # }
  # zeta_1_line <- term_1 + term_2

  # term_1 <- 0
  # for(i in 1:N) {

    # for(j in 1:M) {

      # E_vec <- E_q_zeta_2[[i]][[j]]
      # Cov_mat <- Cov_q_zeta_2[[i]][[j]]
      # A_mat <- E_q_H_psi_22[[i]][[j]]
      # sum_val <- E_cprod(E_vec, Cov_mat, E_vec, A_mat)
      # term_1 <- term_1 + sum_val
    # }
  # }
  # zeta_2_line <- term_1

  # E_q_norm <- y_line + mu_line + zeta_1_line + zeta_2_line

  # term_1 <- sum_T/2*log(2*pi)
  # term_2 <- sum_T/2*E_q_log_sigsq_eps
  # term_3 <- 1/2*E_q_recip_sigsq_eps*E_q_norm

  # ans <- term_1 + term_2 + term_3
  # return(ans)
# }

cross_entropy_logistic_fpc_lik_frag <- function(eta_in, C, Y, L) {

  # order of eta_in:
  # 1. nu -> p(Y|nu,zeta)
  # 2. p(Y|nu,zeta) -> nu
  # 3. zeta -> p(Y|nu,zeta)
  # 4. p(Y|nu,zeta) -> zeta

  # order of eta_out:
  # 1. p(Y|nu,zeta,sigsq_eps) -> nu
  # 2. p(Y|nu,zeta,sigsq_eps) -> zeta

  if(!is.list(eta_in)) {

    stop("eta_in must be a list")
  }

  if(length(eta_in)!=4) {

    stop("eta_in must have length 4")
  }

  if(!is_int(L)) {

    stop("L must be a whole number.")
  } else {

    if(L<1) {

      stop("There must be at least one basis function.")
    }
  }

  N <- ncol(eta_in[[3]])
  l_eta_nu <- length(eta_in[[1]])
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/(L+1) - 2

  eta_in_nu <- list(eta_in[[1]], eta_in[[2]])
  q_nu <- gauss_q(eta_in_nu, use_vech=FALSE)
  mu_q_nu <- q_nu[[1]]
  Sigma_q_nu <- q_nu[[2]]

  mu_q_zeta <- vector("list", length=N)
  Sigma_q_zeta <- vector("list", length=N)
  for(i in 1:N) {

    eta_in_zeta <- list(
      eta_in[[3]][,i],
      eta_in[[4]][,i]
    )
    q_zeta <- gauss_q(eta_in_zeta, use_vech=TRUE)
    mu_q_zeta[[i]] <- q_zeta[[1]]
    Sigma_q_zeta[[i]] <- q_zeta[[2]]
  }

  mu_inds <- 1:(K+2)
  mu_q_nu_mu <- mu_q_nu[mu_inds]
  Sigma_q_nu_mu <- Sigma_q_nu[mu_inds, mu_inds]

  psi_inds <- (1:d)[-mu_inds]
  psi_groups <- ceiling(psi_inds/(K+2)) - 1
  psi_inds <- split(psi_inds, psi_groups)

  M_q_V_psi <- matrix(NA, K+2, L)
  for(l in 1:L) {

    psi_l_inds <- psi_inds[[l]]
    M_q_V_psi[,l] <- mu_q_nu[psi_l_inds]
  }

  Sigma_q_nu_blocks <- vector("list", length = L+1)
  Sigma_q_nu_blocks[[1]] <- vector("list", length = L+1)
  Sigma_q_nu_blocks[[1]][[1]] <- Sigma_q_nu[mu_inds, mu_inds]
  for(l in 1:L) {

    l_inds <- psi_inds[[l]]
    Sigma_q_nu_blocks[[1]][[l+1]] <- Sigma_q_nu[mu_inds, l_inds]

    Sigma_q_nu_blocks[[l+1]] <- vector("list", length = L+1)
    Sigma_q_nu_blocks[[l+1]][[1]] <- Sigma_q_nu[l_inds, mu_inds]
    for(l_p in 1:L) {

      l_p_inds <- psi_inds[[l_p]]
      Sigma_q_nu_blocks[[l+1]][[l_p+1]] <- Sigma_q_nu[l_inds, l_p_inds]
    }
  }

  ans <- 0
  term_1 <- Sigma_q_nu_mu + tcrossprod(mu_q_nu_mu)
  for(i in 1:N) {

    term_2 <- 0
    for(l in 1:L) {

      Cov_mat <- Sigma_q_nu_blocks[[1]][[l+1]]
      mu_outer_mat <- tcrossprod(mu_q_nu_mu, M_q_V_psi[,l])
      prod_term <- Cov_mat + mu_outer_mat

      term_2 <- term_2 + mu_q_zeta[[i]][l]*prod_term
    }

    term_3 <- t(term_2)

    term_4 <- 0
    for(l in 1:L) {

      for(l_p in 1:L) {

        cov_val <- Sigma_q_zeta[[i]][l, l_p]
        mu_prod <- mu_q_zeta[[i]][l]*mu_q_zeta[[i]][l_p]
        E_q_zeta_prod <- cov_val + mu_prod

        Cov_mat <- Sigma_q_nu_blocks[[l+1]][[l_p+1]]
        nu_outer_term <- tcrossprod(M_q_V_psi[,l], M_q_V_psi[,l_p])
        E_q_outer_nu <- Cov_mat + nu_outer_term

        sum_val <- E_q_zeta_prod*E_q_outer_nu
        term_4 <- term_4 + sum_val
      }
    }

    E_q_outer_theta <- term_1 + term_2 + term_3 + term_4
    xi <- sqrt(diag(tcrossprod(C[[i]]%*%E_q_outer_theta, C[[i]])))

    A_xi <- -tanh(xi/2)/(4*xi)
    C_xi <- xi/2 - log(1 + exp(xi)) + xi*tanh(xi/2)/4

    mu_q_zeta_tilde <- c(1, mu_q_zeta[[i]])
    Sigma_q_zeta_tilde <- blkdiag(matrix(0), Sigma_q_zeta[[i]])
    E_q_outer_zeta_tilde <- Sigma_q_zeta_tilde + tcrossprod(mu_q_zeta_tilde)

    sum_val <- cprod(kronecker(t(mu_q_zeta_tilde), C[[i]]), Y[[i]] - 0.5)

    M <- crossprod(C[[i]], diag(A_xi)%*%C[[i]])

    E_q_theta <- mu_q_nu_mu + M_q_V_psi%*%mu_q_zeta[[i]]

    tr_term <- tr(E_q_outer_theta%*%M)
    cprod_term <- cprod(C[[i]]%*%E_q_theta, Y[[i]] - 0.5)
    sum_term <- sum(C_xi)

    ans <- ans + tr_term + cprod_term + sum_term
  }

  return(ans)
}

##########################################
#
#  FPCA  FUNCTIONS
#
##########################################

fpc_orthgn <- function(subj_names, L, K, eta_in, time_g, C_g, Psi_g = NULL) {

  # order of eta_in:
  # 1. p(nu|Sigma_nu) -> nu
  # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # 3. p(zeta) -> zeta
  # 4. p(Y|nu,zeta,sigsq_eps) -> zeta

  N <- ncol(eta_in[[3]])
  l_eta_nu <- length(eta_in[[1]])
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2

  # Determine the original q_nu:

  eta_nu <- list(eta_in[[1]], eta_in[[2]])
  q_nu <- gauss_q(eta_nu, use_vech=FALSE)
  mu_q_nu <- q_nu[[1]]
  Sigma_q_nu <- q_nu[[2]]

  mu_inds <- 1:(K+2)
  psi_inds <- (1:d)[-mu_inds]
  psi_groups <- ceiling(psi_inds/(K+2))-1
  psi_inds <- split(psi_inds, psi_groups)

  mu_q_nu_mu <- mu_q_nu[mu_inds]
  Sigma_q_nu_mu <- Sigma_q_nu[mu_inds, mu_inds]
  mu_q_mu <- as.vector(C_g%*%mu_q_nu_mu)

  mu_q_nu_psi <- vector("list", length=L)
  Sigma_q_nu_psi <- vector("list", length=L)
  for(l in 1:L) {

    psi_l_inds <- psi_inds[[l]]
    mu_q_nu_psi[[l]] <- mu_q_nu[psi_l_inds]
    Sigma_q_nu_psi[[l]] <- Sigma_q_nu[psi_l_inds, psi_l_inds]
  }

  M_q_V_psi <- Reduce(cbind, mu_q_nu_psi)
  M_q_Psi <- C_g%*%M_q_V_psi

  # Determine the original q_zeta:

  mu_q_zeta <- vector("list", length=N)
  Sigma_q_zeta <- vector("list", length=N)
  for(i in 1:N) {

    eta_zeta <- list(eta_in[[3]][,i], eta_in[[4]][,i])
    q_zeta <- gauss_q(eta_zeta, use_vech=TRUE)
    mu_q_zeta[[i]] <- q_zeta[[1]]
    Sigma_q_zeta[[i]] <- q_zeta[[2]]
  }

  M_q_Zeta <- Reduce(rbind, mu_q_zeta)

  # Generate the response curves:

  one_N <- rep(1, N)
  mu_mat <- tcrossprod(mu_q_mu, one_N)
  Y_mat <- mu_mat + tcrossprod(M_q_Psi, M_q_Zeta)

  # Rotate, shift and scale the global curves and scores:

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

  norm_const <- rep(NA, L)
  for(l in 1:L) {

    norm_const[l] <- sqrt(trapint(time_g, (Psi_hat[,l])^2))
    if(norm_const[l]!=0) {

      Psi_hat[,l] <- Psi_hat[,l]/norm_const[l]
      Zeta_hat[,l] <- norm_const[l]*Zeta_hat[,l]

      if(!is.null(Psi_g)) {

        cprod_sign <- sign(cprod(Psi_hat[,l], Psi_g[,l]))
        if(cprod_sign==-1) {

          Psi_hat[,l] <- -Psi_hat[,l]
          Zeta_hat[,l] <- -Zeta_hat[,l]
        }
      }
    }
  }

  mu_q_zeta <- split(Zeta_hat, row(Zeta_hat))

  # Summarise the VMP results:

  scale_mat <- diag(norm_const)
  Cov_zeta_hat <- vector("list", length = N)
  for(i in 1:N) {

    mat_transform <- S_inv%*%tcrossprod(D_diag, V_orth)
    Cov_zeta_dot_hat <- tcrossprod(mat_transform%*%Sigma_q_zeta[[i]], mat_transform)
    Cov_zeta_hat[[i]] <- tcrossprod(scale_mat%*%Cov_zeta_dot_hat, scale_mat)
  }

  Y_hat <- Y_low <- Y_upp <- vector("list", length = N)
  for(i in 1:N) {

      sd_vec_i <- sqrt(diag(tcrossprod(Psi_hat%*%Cov_zeta_hat[[i]], Psi_hat)))

      Y_hat[[i]] <- Y_mat[,i]
      Y_low[[i]] <- Y_mat[,i] + qnorm(0.025)*sd_vec_i
      Y_upp[[i]] <- Y_mat[,i] + qnorm(0.975)*sd_vec_i

  }

  names(Y_hat) <- names(Y_low) <- names(Y_upp) <- names(Cov_zeta_hat) <- subj_names

  list_zeta_ellipse <- vector("list", length=N)
  for(i in 1:N) {

    zeta_mean <- Zeta_hat[i,][1:2]

    zeta_ellipse <- ellipse(
      Cov_zeta_hat[[i]][1:2, 1:2],
      centre=zeta_mean,
      level=0.95
    )

    list_zeta_ellipse[[i]] <- zeta_ellipse
  }
  rownames(Zeta_hat) <- names(list_zeta_ellipse) <- subj_names
  colnames(Zeta_hat) <- paste0("FPC_", 1:L)

  mu_hat <- mu_q_mu
  list_Psi_hat <- Psi_hat

  # Establish the outputs:

  outputs <- create_named_list(time_g, Y_hat, Y_low, Y_upp,
                               mu_hat, list_Psi_hat,
                               Zeta_hat, Cov_zeta_hat,
                               list_zeta_ellipse)
}

mfpc_orthgn <- function(subj_names, var_names, L, K, eta_in,
                        time_g, N, p, C_g, Psi_g = NULL) {

  # order of eta_in:
  # 1. p(nu|Sigma_nu) -> nu
  # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # 3. p(zeta) -> zeta
  # 4. p(Y|nu,zeta,sigsq_eps) -> zeta

  N <- ncol(eta_in[[3]])
  p <- length(eta_in[[1]])
  l_eta_nu <- length(eta_in[[1]][[1]])
  l_eta_zeta <- L + 0.5*L*(L + 1)
  d <- (K+2)*(L+1)
  D_L <- duplication.matrix(L)
  n <- Reduce(rbind, lapply(Y, function(x) sapply(x, length)))
  n_g <- length(time_g)

  E_q_zeta <- vector("list", length = N)
  Cov_q_zeta <- vector("list", length = N)
  for(i in 1:N) {

    eta_in_zeta <- list(
      eta_in[[3]][, i],
      eta_in[[4]][, i]
    )
    q_zeta <- gauss_q(eta_in_zeta, use_vech = TRUE)
    E_q_zeta[[i]] <- q_zeta[[1]]
    Cov_q_zeta[[i]] <- q_zeta[[2]]
  }
  E_q_Zeta <- Reduce(rbind, E_q_zeta)

  E_q_mu <- vector("list", length = p)
  E_q_Psi <- vector("list", length = p)
  for(j in 1:p) {

    eta_in_nu <- list(eta_in[[1]][[j]], eta_in[[2]][[j]])
    q_nu <- gauss_q(eta_in_nu, use_vech = FALSE)
    E_q_nu <- q_nu[[1]]
    E_q_V <- matrix(E_q_nu, K[j] + 2, L + 1)

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
      if(inner_prod_sign == -1) {

        Psi_hat[, l] <- -Psi_hat[, l]
        Zeta_hat[, l] <- -Zeta_hat[, l]
      }
    }
  }
  Psi_hat <- lapply(split(Psi_hat, rep(1:p, each = n_g)), matrix, nrow = n_g, ncol = L)

  Cov_zeta_hat <- vector("list", length = N)
  rotn_mat <- V_psi %*% D_psi %*% Q %*% solve(sqrt(Lambda)) %*% diag(norm_vec)
  for(i in 1:N) {

    Cov_zeta_hat[[i]] <- crossprod(rotn_mat, Cov_q_zeta[[i]] %*% rotn_mat)
  }

  # Store the results:


  list_zeta_ellipse <- vector("list", length = N)
  for(i in 1:N) {

    zeta_mean <- Zeta_hat[i,][1:2]

    zeta_ellipse <- ellipse(
      Cov_zeta_hat[[i]][1:2, 1:2],
      centre=zeta_mean,
      level=0.95
    )

    list_zeta_ellipse[[i]] <- zeta_ellipse

  }

  mu_hat <- as.matrix(Reduce(cbind, mu_hat)) # now deals with the case p = 1
  colnames(mu_hat) <- var_names

  list_Psi_hat <- vector("list", length = L)

  for(l in 1:L) {

    list_Psi_hat[[l]] <- matrix(NA, n_g, p)
    for(j in 1:p) {
      list_Psi_hat[[l]][,j] <- Psi_hat[[j]][, l]
    }
    colnames(list_Psi_hat[[l]]) <- var_names
  }
  names(list_Psi_hat) <- paste0("FPC_", 1:L)

  Y_hat <- Y_low <- Y_upp <- vector("list", length = N)
  for(i in 1:N) {

    Y_hat[[i]] <- Y_low[[i]] <- Y_upp[[i]] <- vector("list", length = p)
    for(j in 1:p) {

      Y_hat_ij <- mu_hat[, j] + Psi_hat[[j]] %*% Zeta_hat[i, ]
      sd_vec_ij <- sqrt(diag(tcrossprod(Psi_hat[[j]] %*% Cov_zeta_hat[[i]], Psi_hat[[j]])))

      Y_hat[[i]][[j]] <- Y_hat_ij
      Y_low[[i]][[j]] <- Y_hat_ij + qnorm(0.025)*sd_vec_ij
      Y_upp[[i]][[j]] <- Y_hat_ij + qnorm(0.975)*sd_vec_ij
    }
    names(Y_hat[[i]]) <- names(Y_low[[i]]) <- names(Y_upp[[i]]) <- var_names

  }

  rownames(Zeta_hat) <- names(list_zeta_ellipse) <- names(Y_hat) <- names(Y_low) <- names(Y_upp) <- subj_names
  colnames(Zeta_hat) <- paste0("FPC_", 1:L)

  # Establish the outputs:

  outputs <- create_named_list(time_g, Y_hat, Y_low, Y_upp,
                               mu_hat, list_Psi_hat,
                               Zeta_hat, Cov_zeta_hat,
                               list_zeta_ellipse)

  return(outputs)
}

mlfpc_orthgn <- function(
  eta_in, time_g, C_g, L_1, L_2,
  N, p, Psi_g = NULL
) {

  # order of eta_in:
  # 1. p(nu|Sigma_nu) -> nu
  # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # 3. p(zeta) -> zeta
  # 4. p(Y|nu,zeta,sigsq_eps) -> zeta

  if(length(eta_in) != 4) {

    stop("eta_in must have length 4")
  }

  if(!is.null(Psi_g)) {

    if(length(Psi_g) != 2) {

      stop("Psi_g must be a list of two matrices.")
    }

    if(ncol(Psi_g[[1]]) != L_1) {

      stop("Psi_g[[1]] must have L_1 columns")
    }

    if(ncol(Psi_g[[2]]) != L_2) {

      stop("Psi_g[[2]] must have L_2 columns")
    }
  }

  N <- length(eta_in[[3]])
  L <- L_1 + L_2
  l_eta_nu <- length(eta_in[[1]])
  l_eta_zeta <- sapply(eta_in[[3]], length)
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/(L+1) - 2
  l_zeta_1 <- L_1 + 0.5*L_1*(L_1 + 1)
  l_zeta_2 <- L_2 + 0.5*L_2*(L_2 + 1) + L_1*L_2
  M <- (l_eta_zeta - l_zeta_1)/l_zeta_2

  eta_nu <- list(eta_in[[1]], eta_in[[2]])
  q_nu <- gauss_q(eta_nu, use_vech = FALSE)
  E_q_nu <- q_nu[[1]]
  Cov_q_nu <- q_nu[[2]]

  E_q_zeta_1 <- vector("list", length = N)
  Cov_q_zeta_1 <- vector("list", length = N)
  E_q_zeta_2 <- vector("list", length = N)
  Cov_q_zeta_2 <- vector("list", length = N)
  Cov_q_zeta_12 <- vector("list", length = N)
  for(i in 1:N) {

    eta_in_zeta <- list(
      eta_in[[3]][[i]],
      eta_in[[4]][[i]]
    )
    q_zeta <- two_level_nat_to_comm_parms(L_1, L_2, M[i], eta_in_zeta)
    E_q_zeta_1[[i]] <- q_zeta[[1]]
    Cov_q_zeta_1[[i]] <- q_zeta[[2]]
    E_q_zeta_2[[i]] <- q_zeta[[3]]
    Cov_q_zeta_2[[i]] <- q_zeta[[4]]
    Cov_q_zeta_12[[i]] <- q_zeta[[5]]
  }

  inds_mat <- matrix(1:d, K + 2, L + 1)
  mu_inds <- inds_mat[, 1]
  psi_inds <- as.vector(inds_mat[, -1])

  E_q_V <- matrix(E_q_nu, nrow = K + 2, ncol = L + 1)
  gbl_mat <- C_g %*% E_q_V
  mu_hat <- gbl_mat[, 1]
  E_q_Psi_1 <- gbl_mat[, 1:L_1 + 1]
  E_q_Psi_2 <- gbl_mat[, 1:L_2 + L_1 + 1]

  Cov_mu_hat <- tcrossprod(C_g %*% Cov_q_nu[mu_inds, mu_inds], C_g)

  E_q_Zeta_1 <- Reduce(rbind, E_q_zeta_1)
  E_q_Zeta_2 <- Reduce(rbind, lapply(E_q_zeta_2, Reduce, f = rbind))

  Psi_1_svd <- svd(E_q_Psi_1)
  U_1 <- Psi_1_svd$u
  D_1 <- diag(Psi_1_svd$d)
  V_1 <- Psi_1_svd$v

  Psi_2_svd <- svd(E_q_Psi_2)
  U_2 <- Psi_2_svd$u
  D_2 <- diag(Psi_2_svd$d)
  V_2 <- Psi_2_svd$v

  Zeta_1_rotn <- E_q_Zeta_1 %*% V_1 %*% D_1
  eig_C_zeta_1 <- eigen(cov(Zeta_1_rotn))
  Q_1 <- eig_C_zeta_1$vectors
  Lambda_1 <- diag(eig_C_zeta_1$values + 1e-10)

  Zeta_2_rotn <- E_q_Zeta_2 %*% V_2 %*% D_2
  eig_C_zeta_2 <- eigen(cov(Zeta_2_rotn))
  Q_2 <- eig_C_zeta_2$vectors
  Lambda_2 <- diag(eig_C_zeta_2$values + 1e-10)

  Psi_1_hat <- U_1 %*% Q_1 %*% sqrt(Lambda_1)
  Zeta_1_hat <- Zeta_1_rotn %*% Q_1 %*% solve(sqrt(Lambda_1))
  norm_psi_1 <- rep(NA, L_1)
  for(l in 1:L_1) {

    norm_val <- sqrt(trapint(time_g, Psi_1_hat[, l]^2))
    Psi_1_hat[, l] <- Psi_1_hat[, l]/norm_val
    Zeta_1_hat[, l] <- norm_val*Zeta_1_hat[, l]

    norm_psi_1[l] <- norm_val

    if(!is.null(Psi_g)) {

      cprod_sign <- sign(cprod(Psi_1_hat[, l], Psi_g[[1]][, l]))
      if(cprod_sign == -1) {

        Psi_1_hat[, l] <- -Psi_1_hat[, l]
        Zeta_1_hat[, l] <- -Zeta_1_hat[, l]
      }
    }
  }
  A_1 <- diag(norm_psi_1)

  Psi_2_hat <- U_2 %*% Q_2 %*% sqrt(Lambda_2)
  Zeta_2_hat <- Zeta_2_rotn %*% Q_2 %*% solve(sqrt(Lambda_2))
  norm_psi_2 <- rep(NA, L_2)
  for(l in 1:L_2) {

    norm_val <- sqrt(trapint(time_g, Psi_2_hat[, l]^2))
    Psi_2_hat[, l] <- Psi_2_hat[, l]/norm_val
    Zeta_2_hat[, l] <- norm_val*Zeta_2_hat[, l]

    norm_psi_2[l] <- norm_val

    if(!is.null(Psi_g)) {

      cprod_sign <- sign(cprod(Psi_2_hat[, l], Psi_g[[2]][, l]))
      if(cprod_sign == -1) {

        Psi_2_hat[, l] <- -Psi_2_hat[, l]
        Zeta_2_hat[, l] <- -Zeta_2_hat[, l]
      }
    }
  }
  Zeta_2_hat <- lapply(
    split(Zeta_2_hat, rep.int(1:N, times = M)),
    matrix, ncol = L_2
  )
  A_2 <- diag(norm_psi_2)

  inds_1 <- 1:L_1
  Cov_zeta_1_hat <- vector("list", length = N)
  Cov_zeta_2_hat <- vector("list", length = N)
  Cov_zeta_12_hat <- vector("list", length = N)
  S_1 <- A_1 %*% solve(sqrt(Lambda_1)) %*% t(Q_1) %*% D_1 %*% V_1
  S_2 <- A_2 %*% solve(sqrt(Lambda_2)) %*% t(Q_2) %*% D_2 %*% V_2
  for(i in 1:N) {

    inds_2 <- matrix(1:(M[i]*L_2) + L_1, L_2, M[i])

    Cov_zeta_1_hat[[i]] <- tcrossprod(S_1 %*% Cov_q_zeta_1[[i]], S_1)
    Cov_zeta_2_hat[[i]] <- vector("list", length = M[i])
    Cov_zeta_12_hat[[i]] <- vector("list", length = M[i])
    for(j in 1:M[i]) {

      Cov_zeta_2_hat[[i]][[j]] <- tcrossprod(S_2 %*% Cov_q_zeta_2[[i]][[j]], S_2)
      Cov_zeta_12_hat[[i]][[j]] <- tcrossprod(S_1 %*% Cov_q_zeta_12[[i]][[j]], S_2)
    }
  }

  # Store the results:

  gbl_mat <- cbind(mu_hat, Psi_1_hat, Psi_2_hat)

  Y_summary <- vector("list", length = N)
  for(i in 1:N) {

    Psi_zeta_1 <- as.vector(Psi_1_hat %*% Zeta_1_hat[i, ])

    Cov_11 <- tcrossprod(Psi_1_hat %*% Cov_zeta_1_hat[[i]], Psi_1_hat)

    Y_summary[[i]] <- vector("list", length = p)
    for(j in 1:p) {

      Psi_zeta_2 <- as.vector(Psi_2_hat %*% Zeta_2_hat[[i]][j, ])

      Cov_22 <- tcrossprod(Psi_2_hat %*% Cov_zeta_2_hat[[i]][[j]], Psi_2_hat)

      Cov_12 <- tcrossprod(Psi_1_hat %*% Cov_zeta_12_hat[[i]][[j]], Psi_2_hat)
      Cov_21 <- t(Cov_12)

      sd_vec <- sqrt(diag(Cov_mu_hat + Cov_11 + Cov_12 + Cov_21) + rep(median(diag(Cov_22)), n_g))
      #sd_vec <- sqrt(diag(Cov_mu_hat + Cov_11 + Cov_22 + Cov_12 + Cov_21))

      Y_hat <- mu_hat + Psi_zeta_1 + Psi_zeta_2

      Y_summary[[i]][[j]] <- matrix(NA, n_g, 3)
      Y_summary[[i]][[j]][, 1] <- Y_hat + qnorm(0.025)*sd_vec
      Y_summary[[i]][[j]][, 2] <- Y_hat
      Y_summary[[i]][[j]][, 3] <- Y_hat + qnorm(0.975)*sd_vec
    }
  }

  ans <- list(Y_summary, gbl_mat, Zeta_1_hat, Zeta_2_hat, Cov_zeta_1_hat, Cov_zeta_2_hat)
  names(ans) <- c("fits", "gbl_curves", "zeta_1", "zeta_2", "Cov_zeta_1", "Cov_zeta_2")
  return(ans)
}

# mlfpca_rotation <- function(
  # eta_in, time_g, C_g, L_1, L_2,
  # N_sample, M_sample, Psi_g = NULL
# ) {

  # # order of eta_in:
  # # 1. p(nu|Sigma_nu) -> nu
  # # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # # 3. p(zeta) -> zeta
  # # 4. p(Y|nu,zeta,sigsq_eps) -> zeta

  # if(length(eta_in) != 4) {

    # stop("eta_in must have length 4")
  # }

  # if(!is.null(Psi_g)) {

    # if(length(Psi_g) != 2) {

      # stop("Psi_g must be a list of two matrices.")
    # }

    # if(ncol(Psi_g[[1]]) != L_1) {

      # stop("Psi_g[[1]] must have L_1 columns")
    # }

    # if(ncol(Psi_g[[2]]) != L_2) {

      # stop("Psi_g[[2]] must have L_2 columns")
    # }
  # }

  # N <- ncol(eta_in[[3]])
  # L <- L_1 + L_2
  # l_eta_nu <- length(eta_in[[1]])
  # l_eta_zeta <- nrow(eta_in[[3]])
  # d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  # K <- d/(L+1) - 2
  # l_zeta_1 <- L_1 + 0.5*L_1*(L_1 + 1)
  # l_zeta_2 <- L_2 + 0.5*L_2*(L_2 + 1) + L_1*L_2
  # M <- (l_eta_zeta - l_zeta_1)/l_zeta_2
  # n_sample <- length(N_sample)
  # m_sample <- length(M_sample)

  # eta_nu <- list(eta_in[[1]], eta_in[[2]])
  # q_nu <- gauss_q(eta_nu, use_vech = FALSE)
  # E_q_nu <- q_nu[[1]]
  # Cov_q_nu <- q_nu[[2]]

  # E_q_zeta_1 <- vector("list", length = N)
  # Cov_q_zeta_1 <- vector("list", length = N)
  # E_q_zeta_2 <- vector("list", length = N)
  # Cov_q_zeta_2 <- vector("list", length = N)
  # Cov_q_zeta_12 <- vector("list", length = N)
  # for(i in 1:N) {

    # eta_in_zeta <- list(
      # eta_in[[3]][, i],
      # eta_in[[4]][, i]
    # )
    # q_zeta <- two_level_nat_to_comm_parms(L_1, L_2, M, eta_in_zeta)
    # E_q_zeta_1[[i]] <- q_zeta[[1]]
    # Cov_q_zeta_1[[i]] <- q_zeta[[2]]
    # E_q_zeta_2[[i]] <- q_zeta[[3]]
    # Cov_q_zeta_2[[i]] <- q_zeta[[4]]
    # Cov_q_zeta_12[[i]] <- q_zeta[[5]]
  # }

  # inds_mat <- matrix(1:d, K + 2, L + 1)
  # mu_inds <- inds_mat[, 1]
  # psi_inds <- as.vector(inds_mat[, -1])

  # E_q_V <- matrix(E_q_nu, nrow = K + 2, ncol = L + 1)
  # gbl_mat <- C_g %*% E_q_V
  # mu_hat <- gbl_mat[, 1]
  # E_q_Psi_1 <- gbl_mat[, 1:L_1 + 1]
  # E_q_Psi_2 <- gbl_mat[, 1:L_2 + L_1 + 1]

  # Cov_mu_hat <- tcrossprod(C_g %*% Cov_q_nu[mu_inds, mu_inds], C_g)

  # E_q_Zeta_1 <- Reduce(rbind, E_q_zeta_1)
  # E_q_Zeta_2 <- Reduce(rbind, lapply(E_q_zeta_2, Reduce, f = rbind))

  # Psi_1_svd <- svd(E_q_Psi_1)
  # U_1 <- Psi_1_svd$u
  # D_1 <- diag(Psi_1_svd$d)
  # V_1 <- Psi_1_svd$v

  # Psi_2_svd <- svd(E_q_Psi_2)
  # U_2 <- Psi_2_svd$u
  # D_2 <- diag(Psi_2_svd$d)
  # V_2 <- Psi_2_svd$v

  # Zeta_1_rotn <- E_q_Zeta_1 %*% V_1 %*% D_1
  # eig_C_zeta_1 <- eigen(cov(Zeta_1_rotn))
  # Q_1 <- eig_C_zeta_1$vectors
  # Lambda_1 <- diag(eig_C_zeta_1$values)

  # Zeta_2_rotn <- E_q_Zeta_2 %*% V_2 %*% D_2
  # eig_C_zeta_2 <- eigen(cov(Zeta_2_rotn))
  # Q_2 <- eig_C_zeta_2$vectors
  # Lambda_2 <- diag(eig_C_zeta_2$values)

  # Psi_1_hat <- U_1 %*% Q_1 %*% sqrt(Lambda_1)
  # Zeta_1_hat <- Zeta_1_rotn %*% Q_1 %*% solve(sqrt(Lambda_1))
  # norm_psi_1 <- rep(NA, L_1)
  # for(l in 1:L_1) {

    # norm_val <- sqrt(trapint(time_g, Psi_1_hat[, l]^2))
    # Psi_1_hat[, l] <- Psi_1_hat[, l]/norm_val
    # Zeta_1_hat[, l] <- norm_val*Zeta_1_hat[, l]

    # norm_psi_1[l] <- norm_val

    # if(!is.null(Psi_g)) {

      # cprod_sign <- sign(cprod(Psi_1_hat[, l], Psi_g[[1]][, l]))
      # if(cprod_sign == -1) {

        # Psi_1_hat[, l] <- -Psi_1_hat[, l]
        # Zeta_1_hat[, l] <- -Zeta_1_hat[, l]
      # }
    # }
  # }
  # A_1 <- diag(norm_psi_1)

  # Psi_2_hat <- U_2 %*% Q_2 %*% sqrt(Lambda_2)
  # Zeta_2_hat <- Zeta_2_rotn %*% Q_2 %*% solve(sqrt(Lambda_2))
  # norm_psi_2 <- rep(NA, L_2)
  # for(l in 1:L_2) {

    # norm_val <- sqrt(trapint(time_g, Psi_2_hat[, l]^2))
    # Psi_2_hat[, l] <- Psi_2_hat[, l]/norm_val
    # Zeta_2_hat[, l] <- norm_val*Zeta_2_hat[, l]

    # norm_psi_2[l] <- norm_val

    # if(!is.null(Psi_g)) {

      # cprod_sign <- sign(cprod(Psi_2_hat[, l], Psi_g[[2]][, l]))
      # if(cprod_sign == -1) {

        # Psi_2_hat[, l] <- -Psi_2_hat[, l]
        # Zeta_2_hat[, l] <- -Zeta_2_hat[, l]
      # }
    # }
  # }
  # Zeta_2_hat <- lapply(
    # split(Zeta_2_hat, ceiling(row(Zeta_2_hat)/M)),
    # matrix, nrow = M, ncol = L_2
  # )
  # A_2 <- diag(norm_psi_2)

  # inds_1 <- 1:L_1
  # inds_2 <- matrix(1:(M*L_2) + L_1, L_2, M)
  # Cov_zeta_1_hat <- vector("list", length = N)
  # Cov_zeta_2_hat <- vector("list", length = N)
  # Cov_zeta_12_hat <- vector("list", length = N)
  # S_1 <- A_1 %*% solve(sqrt(Lambda_1)) %*% t(Q_1) %*% D_1 %*% V_1
  # S_2 <- A_2 %*% solve(sqrt(Lambda_2)) %*% t(Q_2) %*% D_2 %*% V_2
  # for(i in 1:N) {

    # Cov_zeta_1_hat[[i]] <- tcrossprod(S_1 %*% Cov_q_zeta_1[[i]], S_1)
    # Cov_zeta_2_hat[[i]] <- vector("list", length = M)
    # Cov_zeta_12_hat[[i]] <- vector("list", length = M)
    # for(j in 1:M) {

      # Cov_zeta_2_hat[[i]][[j]] <- tcrossprod(S_2 %*% Cov_q_zeta_2[[i]][[j]], S_2)
      # Cov_zeta_12_hat[[i]][[j]] <- tcrossprod(S_1 %*% Cov_q_zeta_12[[i]][[j]], S_2)
    # }
    # }

  # # Store the results:

  # gbl_mat <- cbind(mu_hat, Psi_1_hat, Psi_2_hat)

  # Y_summary <- vector("list", length = n_sample)
  # for(i in 1:n_sample) {

    # N_i <- N_sample[i]

    # Psi_zeta_1 <- as.vector(Psi_1_hat %*% Zeta_1_hat[N_i, ])

    # Cov_11 <- tcrossprod(Psi_1_hat %*% Cov_zeta_1_hat[[N_i]], Psi_1_hat)

    # Y_summary[[i]] <- vector("list", length = m_sample)
    # for(j in 1:m_sample) {

      # M_j <- M_sample[j]

      # Psi_zeta_2 <- as.vector(Psi_2_hat %*% Zeta_2_hat[[N_i]][M_j, ])

      # Cov_22 <- tcrossprod(Psi_2_hat %*% Cov_zeta_2_hat[[N_i]][[M_j]], Psi_2_hat)

      # Cov_12 <- tcrossprod(Psi_1_hat %*% Cov_zeta_12_hat[[N_i]][[M_j]], Psi_2_hat)
      # Cov_21 <- t(Cov_12)

      # sd_vec <- sqrt(diag(Cov_mu_hat + Cov_11 + Cov_22 + Cov_12 + Cov_21))

      # Y_hat <- mu_hat + Psi_zeta_1 + Psi_zeta_2

      # Y_summary[[i]][[j]] <- matrix(NA, n_g, 3)
      # Y_summary[[i]][[j]][, 1] <- Y_hat + qnorm(0.025)*sd_vec
      # Y_summary[[i]][[j]][, 2] <- Y_hat
      # Y_summary[[i]][[j]][, 3] <- Y_hat + qnorm(0.975)*sd_vec
    # }
  # }

  # ans <- list(Y_summary, gbl_mat, Zeta_1_hat, Zeta_2_hat)
  # names(ans) <- c("fits", "gbl_curves", "zeta_1", "zeta_2")
  # return(ans)
# }

fpc_orthgn_two_lev <- function(eta_in, time_g, C_g, Psi_1_g, Psi_2_g) {

  # order of eta_in:
  # 1. p(nu|Sigma_nu) -> nu
  # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # 3. p(zeta) -> zeta
  # 4. p(Y|nu,zeta,sigsq_eps) -> zeta

  N <- ncol(eta_in[[3]])
  L_1 <- ncol(Psi_1_g)
  L_2 <- ncol(Psi_2_g)
  L <- L_1 + L_2
  l_eta_nu <- length(eta_in[[1]])
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/(L+1) - 2

  # Determine the original q_nu:

  eta_nu <- list(eta_in[[1]], eta_in[[2]])
  q_nu <- gauss_q(eta_nu, use_vech = FALSE)
  E_q_nu <- q_nu[[1]]
  Cov_q_nu <- q_nu[[2]]

  # Store mean function and eigenfunctions:

  E_q_V <- matrix(E_q_nu, K + 2, L + 1)
  gbl_mat <- C_g %*% E_q_V
  E_q_mu <- gbl_mat[, 1]
  E_q_Psi_1 <- gbl_mat[, 1:L_1 + 1]
  E_q_Psi_2 <- gbl_mat[, 1:L_2 + L_1 + 1]

  # Store scores:

  one_MN <- rep(1, M*N)
  one_M <- rep(1, M)

  E_q_Zeta_1 <- vector("list", length = N)
  E_q_Zeta_2 <- vector("list", length = N)
  inds_1 <- 1:L_1
  inds_2 <- 1:(M*L_2) + L_1
  for(i in 1:N) {

    E_q_Zeta_1[[i]] <- tcrossprod(one_M, E_q_zeta[[i]][inds_1])
    E_q_Zeta_2[[i]] <- t(matrix(E_q_zeta[[i]][inds_2], L_2, M))
  }
  E_q_Zeta_1 <- Reduce(rbind, E_q_Zeta_1)
  E_q_Zeta_2 <- Reduce(rbind, E_q_Zeta_2)

  # Orthogonalization and Rotation:

  E_q_Psi_1_svd <- svd(E_q_Psi_1)
  U_Psi_1 <- E_q_Psi_1_svd$u
  D_Psi_1 <- diag(E_q_Psi_1_svd$d)
  V_Psi_1 <- E_q_Psi_1_svd$v

  E_q_Psi_2_svd <- svd(E_q_Psi_2)
  U_Psi_2 <- E_q_Psi_2_svd$u
  D_Psi_2 <- diag(E_q_Psi_2_svd$d)
  V_Psi_2 <- E_q_Psi_2_svd$v

  Zeta_1_rotn <- E_q_Zeta_1 %*% V_Psi_1 %*% D_Psi_1
  m_zeta_1 <- apply(Zeta_1_rotn, 2, mean)
  Zeta_1_shift <- Zeta_1_rotn - tcrossprod(one_MN, m_zeta_1)

  Zeta_2_rotn <- E_q_Zeta_2 %*% V_Psi_2 %*% D_Psi_2
  m_zeta_2 <- apply(Zeta_2_rotn, 2, mean)
  Zeta_2_shift <- Zeta_2_rotn - tcrossprod(one_MN, m_zeta_2)

  eigen_Zeta_1_shift <- eigen(cov(Zeta_1_shift))
  Q_1 <- eigen_Zeta_1_shift$vectors
  Lambda_1 <- diag(eigen_Zeta_1_shift$values + 1e-10)
  Lambda_1_inv <- diag(1/(eigen_Zeta_1_shift$values + 1e-10))
  S_1 <- Q_1 %*% sqrt(Lambda_1)
  S_1_inv <- tcrossprod(sqrt(Lambda_1_inv), Q_1)

  eigen_Zeta_2_shift <- eigen(cov(Zeta_2_shift))
  Q_2 <- eigen_Zeta_2_shift$vectors
  Lambda_2 <- diag(eigen_Zeta_2_shift$values + 1e-10)
  Lambda_2_inv <- diag(1/(eigen_Zeta_2_shift$values + 1e-10))
  S_2 <- Q_2 %*% sqrt(Lambda_2)
  S_2_inv <- tcrossprod(sqrt(Lambda_2_inv), Q_2)

  mu_hat <- E_q_mu + as.vector(U_Psi_1 %*% m_zeta_1 + U_Psi_2 %*% m_zeta_2)

  Psi_1_hat <- U_Psi_1 %*% S_1
  Zeta_1_hat <- tcrossprod(Zeta_1_shift, S_1_inv)

  Psi_2_hat <- U_Psi_2 %*% S_2
  Zeta_2_hat <- tcrossprod(Zeta_2_shift, S_2_inv)

  norm_const_1 <- rep(NA, L_1)
  for(l in 1:L_1) {

    norm_const_1[l] <- sqrt(trapint(time_g, (Psi_1_hat[,l])^2))
    if(norm_const_1[l] != 0) {

      Psi_1_hat[,l] <- Psi_1_hat[,l]/norm_const_1[l]
      Zeta_1_hat[,l] <- norm_const_1[l]*Zeta_1_hat[,l]

      cprod_sign <- sign(cprod(Psi_1_hat[,l], Psi_1_g[,l]))
      if(cprod_sign == -1) {

        Psi_1_hat[,l] <- -Psi_1_hat[,l]
        Zeta_1_hat[,l] <- -Zeta_1_hat[,l]
      }
    }
  }

  norm_const_2 <- rep(NA, L_2)
  for(l in 1:L_2) {

    norm_const_2[l] <- sqrt(trapint(time_g, (Psi_2_hat[,l])^2))
    if(norm_const_2[l] != 0) {

      Psi_2_hat[,l] <- Psi_2_hat[,l]/norm_const_2[l]
      Zeta_2_hat[,l] <- norm_const_2[l]*Zeta_2_hat[,l]

      cprod_sign <- sign(cprod(Psi_2_hat[,l], Psi_2_g[,l]))
      if(cprod_sign == -1) {

        Psi_2_hat[,l] <- -Psi_2_hat[,l]
        Zeta_2_hat[,l] <- -Zeta_2_hat[,l]
      }
    }
  }

  # Store the rotated parameters for the scores:

  Zeta_1_hat <- apply(Zeta_1_hat, 2, unique)
  Zeta_2_hat <- lapply(split(Zeta_2_hat, rep(1:N, each = M)), function(x) matrix(x, M, L_2))

  scale_mat_1 <- diag(norm_const_1)
  scale_mat_2 <- diag(norm_const_2)
  transform_mat_1 <- scale_mat_1 %*% S_1_inv %*% tcrossprod(D_Psi_1, V_Psi_1)
  transform_mat_2 <- scale_mat_2 %*% S_2_inv %*% tcrossprod(D_Psi_2, V_Psi_2)
  transform_mat <- adiag(transform_mat_1, Reduce(adiag, rep(list(transform_mat_2), M)))
  Cov_zeta_hat <- vector("list", length = N)
  for(i in 1:N) {

    Cov_zeta_hat[[i]] <- tcrossprod(transform_mat %*% Cov_q_zeta[[i]], transform_mat)
  }

  # Establish the outputs:

  outputs <- list(mu_hat, Psi_1_hat, Psi_2_hat, Zeta_1_hat, Zeta_2_hat, Cov_zeta_hat)
  names(outputs) <- c(
    "mu_hat", "Psi_1_hat", "Psi_2_hat",
    "Zeta_1_hat", "Zeta_2_hat", "Cov_zeta_hat"
  )

  return(outputs)
}

fpca_mc <- function(q_nu, q_zeta, n_mc, time_g, C_g, Psi_g) {

  # components of q_nu:
  # mu_q_nu
  # Sigma_q_nu

  # components of q_zeta:
  # mu_q_zeta
  # Sigma_q_zeta

  mu_q_nu <- q_nu[[1]]
  Sigma_q_nu <- q_nu[[2]]

  mu_q_zeta <- q_zeta[[1]]
  Sigma_q_zeta <- q_zeta[[2]]

  N <- length(mu_q_zeta)
  L <- length(mu_q_nu_psi)
  n_g <- length(time_g)

  Y_mc <- vector("list", length=N)
  zeta_star_mc <- vector("list", length=N)
  for(i in 1:N) {

    Y_mc[[i]] <- matrix(NA, n_mc, n_g)
    zeta_star_mc[[i]] <- matrix(NA, n_mc, L)
  }

  psi_star_mc <- vector("list", length=L)
  for(l in 1:L) {

    psi_star_mc[[l]] <- matrix(NA, n_mc, n_g)
  }

  for(j in 1:n_mc) {

    cat("starting Monte Carlo sample", j, "of", n_mc, "\n")

    nu_mc <- mvrnorm(1, mu_q_nu, Sigma_q_nu)

    mu_inds <- 1:(K+2)
    nu_mu_mc <- nu_mc[mu_inds]
    mu_mc <- C_g%*%nu_mu_mc

    psi_inds <- (1:d)[-mu_inds]
    psi_groups <- ceiling(psi_inds/(K+2))-1
    psi_inds <- split(psi_inds, psi_groups)
    Psi_mc <- matrix(NA, n_g, L)
    for(l in 1:L) {

      psi_inds_l <- psi_inds[[l]]
      nu_psi_mc <- nu_mc[psi_inds_l]
      Psi_mc[,l] <- C_g%*%nu_psi_mc
    }

    Zeta_mc <- matrix(NA, N, L)
    for(i in 1:N) {

      Zeta_mc[i,] <- mvrnorm(1, mu_q_zeta[[i]], Sigma_q_zeta[[i]])
    }

    svd_list <- svd(Psi_mc)
    U_orth <- svd_list$u
    D_diag <- diag(svd_list$d)
    V_orth <- svd_list$v

    Psi_star_mc <- U_orth
    Zeta_star_mc <- Zeta_mc%*%V_orth%*%D_diag
    for(l in 1:L) {

      norm_const <- trapint(time_g, Psi_star_mc[,l]^2)
      Psi_star_mc[,l] <- 1/sqrt(norm_const)*Psi_star_mc[,l]
      Zeta_star_mc[,l] <- sqrt(norm_const)*Zeta_star_mc[,l]

      cprod_test <- cprod(Psi_g[,l], Psi_star_mc[,l])

      if(cprod_test < 0) {

        Psi_star_mc[,l] <- -Psi_star_mc[,l]
        Zeta_star_mc[,l] <- -Zeta_star_mc[,l]
      }

      psi_star_mc[[l]][j,] <- Psi_star_mc[,l]
    }

    for(i in 1:N) {

      zeta_star_mc[[i]][j,] <- Zeta_star_mc[i,]

      zeta_mc <- Zeta_mc[i,]
      Y_mc[[i]][j,] <- mu_mc + as.vector(Psi_mc%*%zeta_mc)
    }
  }

  output <- list(Y_mc, psi_star_mc, zeta_star_mc)
  names(output) <- c("Y", "psi_star", "zeta_star")

  return(output)
}

logistic_fpca_mc <- function(q_nu, q_zeta, n_mc, time_g, C_g, Psi_g) {

  # components of q_nu:
  # mu_q_nu
  # Sigma_q_nu

  # components of q_zeta:
  # mu_q_zeta
  # Sigma_q_zeta

  mu_q_nu <- q_nu[[1]]
  Sigma_q_nu <- q_nu[[2]]

  mu_q_zeta <- q_zeta[[1]]
  Sigma_q_zeta <- q_zeta[[2]]

  N <- length(mu_q_zeta)
  L <- length(mu_q_nu_psi)
  n_g <- length(time_g)

  Y_mc <- vector("list", length=N)
  zeta_star_mc <- vector("list", length=N)
  for(i in 1:N) {

    Y_mc[[i]] <- matrix(NA, n_mc, n_g)
    zeta_star_mc[[i]] <- matrix(NA, n_mc, L)
  }

  psi_star_mc <- vector("list", length=L)
  for(l in 1:L) {

    psi_star_mc[[l]] <- matrix(NA, n_mc, n_g)
  }

  for(j in 1:n_mc) {

    cat("starting Monte Carlo sample", j, "of", n_mc, "\n")

    nu_mc <- mvrnorm(1, mu_q_nu, Sigma_q_nu)

    mu_inds <- 1:(K+2)
    nu_mu_mc <- nu_mc[mu_inds]
    mu_mc <- C_g%*%nu_mu_mc

    psi_inds <- (1:d)[-mu_inds]
    psi_groups <- ceiling(psi_inds/(K+2))-1
    psi_inds <- split(psi_inds, psi_groups)
    Psi_mc <- matrix(NA, n_g, L)
    for(l in 1:L) {

      psi_inds_l <- psi_inds[[l]]
      nu_psi_mc <- nu_mc[psi_inds_l]
      Psi_mc[,l] <- C_g%*%nu_psi_mc
    }

    Zeta_mc <- matrix(NA, N, L)
    for(i in 1:N) {

      Zeta_mc[i,] <- mvrnorm(1, mu_q_zeta[[i]], Sigma_q_zeta[[i]])
    }

    svd_list <- svd(Psi_mc)
    U_orth <- svd_list$u
    D_diag <- diag(svd_list$d)
    V_orth <- svd_list$v

    Psi_star_mc <- U_orth
    Zeta_star_mc <- Zeta_mc%*%V_orth%*%D_diag
    for(l in 1:L) {

      norm_const <- trapint(time_g, Psi_star_mc[,l]^2)
      Psi_star_mc[,l] <- 1/sqrt(norm_const)*Psi_star_mc[,l]
      Zeta_star_mc[,l] <- sqrt(norm_const)*Zeta_star_mc[,l]

      cprod_test <- cprod(Psi_g[,l], Psi_star_mc[,l])

      if(cprod_test < 0) {

        Psi_star_mc[,l] <- -Psi_star_mc[,l]
        Zeta_star_mc[,l] <- -Zeta_star_mc[,l]
      }

      psi_star_mc[[l]][j,] <- Psi_star_mc[,l]
    }

    for(i in 1:N) {

      zeta_star_mc[[i]][j,] <- Zeta_star_mc[i,]

      zeta_mc <- Zeta_mc[i,]
      Y_mc[[i]][j,] <- inv_logit(mu_mc + as.vector(Psi_mc%*%zeta_mc))
    }
  }

  output <- list(Y_mc, psi_star_mc, zeta_star_mc)
  names(output) <- c("Y", "psi_star", "zeta_star")

  return(output)
}

##########################################
#
#  MULTILEVEL  FUNCTIONS
#
##########################################

solve_two_lev_sparse_mat <- function(a_1, A_11, a_2, A_22, A_12) {

  if(!is.list(a_2)) {

    stop("a_2 must be a list")
  }

  if(!is.list(A_22)) {

    stop("A_22 must be a list")
  }

  if(!is.list(A_12)) {

    stop("A_12 must be a list")
  }

  m <- length(a_2)
  q <- length(a_2[[1]])

  omega <- a_1
  Omega <- A_11
  for(i in 1:m) {

    omega <- omega - as.vector(A_12[[i]] %*% solve(A_22[[i]], a_2[[i]]))
    Omega <- Omega - A_12[[i]] %*% solve(A_22[[i]], t(A_12[[i]]))
  }

  A_inv_11 <- solve(Omega)
  x_1 <- as.vector(A_inv_11 %*% omega)

  x_2 <- vector("list", length = m)
  A_inv_12 <- vector("list", length = m)
  A_inv_22 <- vector("list", length = m)
  for(i in 1:m) {

    x_2[[i]] <- solve(A_22[[i]], a_2[[i]] - cprod(A_12[[i]], x_1))
    A_inv_12[[i]] <- -t(solve(A_22[[i]], crossprod(A_12[[i]], A_inv_11)))
    A_inv_22[[i]] <- solve(A_22[[i]], diag(q) - crossprod(A_12[[i]], A_inv_12[[i]]))
  }

  ans <- list(x_1, A_inv_11, x_2, A_inv_22, A_inv_12)
  names(ans) <- c("x_1", "A_inv_11", "x_2", "A_inv_22", "A_inv_12")
  return(ans)
}

two_level_nat_to_comm_parms <- function(p, q, m, eta_nu_in) {

  # p is the length of the first level vector
  # q is the length of the second level vector

  if(!is.list(eta_nu_in)) {

    stop("eta_nu_in must be a list")
  }

  if(length(eta_nu_in) != 2) {

    stop("eta_nu_in must be a list of length 2")
  }

  eta_q_nu <- eta_nu_in[[1]] + eta_nu_in[[2]]

  omega_1 <- eta_q_nu[1:p]
  omega_2 <- eta_q_nu[(p+1):(p + 0.5*p*(p + 1))]

  D_p <- duplication.matrix(p)
  D_p_plus <- solve(crossprod(D_p), t(D_p))

  D_q <- duplication.matrix(q)
  D_q_plus <- solve(crossprod(D_q), t(D_q))

  Omega_1 <- -2*vecInverse(cprod(D_p_plus, omega_2))

  i_stt <- p + 0.5*p*(p + 1) + 1
  i_end <- i_stt + q - 1
  omega_3 <- vector("list", length = m)
  Omega_2 <- vector("list", length = m)
  Omega_3 <- vector("list", length = m)
  for(i in 1:m) {

    omega_3[[i]] <- eta_q_nu[i_stt:i_end]

    i_stt <- i_end + 1
    i_end <- i_stt + 0.5*q*(q + 1) - 1

    omega_4 <- eta_q_nu[i_stt:i_end]

    i_stt <- i_end + 1
    i_end <- i_stt + p*q - 1

    omega_5 <- eta_q_nu[i_stt:i_end]

    i_stt <- i_end + 1
    i_end <- i_stt + q - 1

    Omega_2[[i]] <- -2*vecInverse(cprod(D_q_plus, omega_4))
    Omega_3[[i]] <- -matrix(omega_5, p, q)
  }

  ans <- solve_two_lev_sparse_mat(omega_1, Omega_1, omega_3, Omega_2, Omega_3)
  names(ans) <- c("E_q_nu_1", "Cov_q_nu_1", "E_q_nu_2", "Cov_q_nu_2", "Cov_q_nu_12")
  return(ans)
}

two_level_fpc_orthogonalization <- function(eta_in, p, time_g, C_g, Psi_1_g, Psi_2_g) {

  # order of eta_in:
  # 1. p(nu|Sigma_nu) -> nu
  # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # 3. p(zeta) -> zeta
  # 4. p(Y|nu,zeta,sigsq_eps) -> zeta

  N <- ncol(eta_in[[3]])
  L_1 <- ncol(Psi_1_g)
  L_2 <- ncol(Psi_2_g)
  L <- L_1 + L_2
  l_eta_nu <- length(eta_in[[1]])
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/(L+1) - 2

  # Compute q(nu):

  eta_nu <- list(eta_in[[1]], eta_in[[2]])
  q_nu <- gauss_q(eta_nu, use_vech = FALSE)
  E_q_nu <- q_nu[[1]]
  Cov_q_nu <- q_nu[[2]]

  # Compute q(zeta):

  E_q_zeta_1 <- vector("list", length = N)
  Cov_q_zeta_1 <- vector("list", length = N)
  E_q_zeta_2 <- vector("list", length = N)
  Cov_q_zeta_2 <- vector("list", length = N)
  Cov_q_zeta_12 <- vector("list", length = N)
  for(i in 1:N) {

    eta_in_zeta <- list(eta_in[[3]][, i], eta_in[[4]][, i])
    q_zeta <- two_level_nat_to_comm_parms(L_1, L_2, M, eta_in_zeta)
    E_q_zeta_1[[i]] <- q_zeta[[1]]
    Cov_q_zeta_1[[i]] <- q_zeta[[2]]
    E_q_zeta_2[[i]] <- q_zeta[[3]]
    Cov_q_zeta_2[[i]] <- q_zeta[[4]]
    Cov_q_zeta_12[[i]] <- q_zeta[[5]]
  }

  # Store mean function and eigenfunctions:

  E_q_V <- matrix(E_q_nu, K + 2, L + 1)
  gbl_mat <- C_g %*% E_q_V
  E_q_mu <- gbl_mat[, 1]
  E_q_Psi_1 <- gbl_mat[, 1:L_1 + 1]
  E_q_Psi_2 <- gbl_mat[, 1:L_2 + L_1 + 1]

  # Store scores:

  E_q_Zeta_1 <- Reduce(rbind, E_q_zeta_1)

  E_q_Zeta_2 <- vector("list", length = N)
  for(i in 1:N) {

    E_q_Zeta_2[[i]] <- Reduce(rbind, E_q_zeta_2[[i]])
  }
  E_q_Zeta_2 <- Reduce(rbind, E_q_Zeta_2)

  # Orthogonalization and Rotation:

  E_q_Psi_1_svd <- svd(E_q_Psi_1)
  U_Psi_1 <- E_q_Psi_1_svd$u
  D_Psi_1 <- diag(E_q_Psi_1_svd$d)
  V_Psi_1 <- E_q_Psi_1_svd$v

  E_q_Psi_2_svd <- svd(E_q_Psi_2)
  U_Psi_2 <- E_q_Psi_2_svd$u
  D_Psi_2 <- diag(E_q_Psi_2_svd$d)
  V_Psi_2 <- E_q_Psi_2_svd$v

  Zeta_1_rotn <- E_q_Zeta_1 %*% V_Psi_1 %*% D_Psi_1
  Zeta_2_rotn <- E_q_Zeta_2 %*% V_Psi_2 %*% D_Psi_2

  eigen_Zeta_1_shift <- eigen(cov(Zeta_1_rotn))
  Q_1 <- eigen_Zeta_1_shift$vectors
  Lambda_1 <- diag(eigen_Zeta_1_shift$values + 1e-10)
  Lambda_1_inv <- diag(1/(eigen_Zeta_1_shift$values + 1e-10))
  S_1 <- Q_1 %*% sqrt(Lambda_1)
  S_1_inv <- tcrossprod(sqrt(Lambda_1_inv), Q_1)

  eigen_Zeta_2_shift <- eigen(cov(Zeta_2_rotn))
  Q_2 <- eigen_Zeta_2_shift$vectors
  Lambda_2 <- diag(eigen_Zeta_2_shift$values + 1e-10)
  Lambda_2_inv <- diag(1/(eigen_Zeta_2_shift$values + 1e-10))
  S_2 <- Q_2 %*% sqrt(Lambda_2)
  S_2_inv <- tcrossprod(sqrt(Lambda_2_inv), Q_2)

  mu_hat <- E_q_mu
  Cov_mu_hat <- tcrossprod(C_g %*% Cov_q_nu[1:(K + 2), 1:(K + 2)], C_g)

  Psi_1_hat <- U_Psi_1 %*% S_1
  Zeta_1_hat <- tcrossprod(Zeta_1_rotn, S_1_inv)

  Psi_2_hat <- U_Psi_2 %*% S_2
  Zeta_2_hat <- tcrossprod(Zeta_2_rotn, S_2_inv)

  norm_const_1 <- rep(NA, L_1)
  for(l in 1:L_1) {

    norm_const_1[l] <- sqrt(trapint(time_g, (Psi_1_hat[,l])^2))
    if(norm_const_1[l] != 0) {

      Psi_1_hat[,l] <- Psi_1_hat[,l]/norm_const_1[l]
      Zeta_1_hat[,l] <- norm_const_1[l]*Zeta_1_hat[,l]

      cprod_sign <- sign(cprod(Psi_1_hat[,l], Psi_1_g[,l]))
      if(cprod_sign == -1) {

        Psi_1_hat[,l] <- -Psi_1_hat[,l]
        Zeta_1_hat[,l] <- -Zeta_1_hat[,l]
      }
    }
  }

  norm_const_2 <- rep(NA, L_2)
  for(l in 1:L_2) {

    norm_const_2[l] <- sqrt(trapint(time_g, (Psi_2_hat[,l])^2))
    if(norm_const_2[l] != 0) {

      Psi_2_hat[,l] <- Psi_2_hat[,l]/norm_const_2[l]
      Zeta_2_hat[,l] <- norm_const_2[l]*Zeta_2_hat[,l]

      cprod_sign <- sign(cprod(Psi_2_hat[,l], Psi_2_g[,l]))
      if(cprod_sign == -1) {

        Psi_2_hat[,l] <- -Psi_2_hat[,l]
        Zeta_2_hat[,l] <- -Zeta_2_hat[,l]
      }
    }
  }
  Zeta_2_hat <- lapply(split(Zeta_2_hat, rep(1:N, each = M)), function(x) matrix(x, M, L_2))

  # Update the posterior covariance matrix for each of the scores:

  scale_mat_1 <- diag(norm_const_1)
  scale_mat_2 <- diag(norm_const_2)
  transform_mat_1 <- scale_mat_1 %*% S_1_inv %*% tcrossprod(D_Psi_1, V_Psi_1)
  transform_mat_2 <- scale_mat_2 %*% S_2_inv %*% tcrossprod(D_Psi_2, V_Psi_2)
  Cov_zeta_1_hat <- vector("list", length = N)
  Cov_zeta_2_hat <- vector("list", length = N)
  Cov_zeta_12_hat <- vector("list", length = N)
  for(i in 1:N) {

    Cov_zeta_1_hat[[i]] <- tcrossprod(
      transform_mat_1 %*% Cov_q_zeta_1[[i]],
      transform_mat_1
    )

    Cov_zeta_2_hat[[i]] <- vector("list", length = M)
    Cov_zeta_12_hat[[i]] <- vector("list", length = M)
    for(j in 1:M) {

      Cov_zeta_2_hat[[i]][[j]] <- tcrossprod(
        transform_mat_2 %*% Cov_q_zeta_2[[i]][[j]],
        transform_mat_2
      )

      Cov_zeta_12_hat[[i]][[j]] <- tcrossprod(
        transform_mat_1 %*% Cov_q_zeta_12[[i]][[j]],
        transform_mat_2
      )
    }
  }

  # Summarise the results:

  Y_summary <- vector("list", length = N)
  for(i in 1:N) {

    Y_summary[[i]] <- vector("list", length = p)

    Var_11 <- tcrossprod(Psi_1_hat %*% Cov_zeta_1_hat[[i]], Psi_1_hat)

    for(j in 1:p) {

      Var_22 <- tcrossprod(Psi_2_hat %*% Cov_zeta_2_hat[[i]][[j]], Psi_2_hat)

      Cov_12 <- tcrossprod(Psi_1_hat %*% Cov_zeta_12_hat[[i]][[j]], Psi_2_hat)
      Cov_21 <- t(Cov_12)

      sd_vec <- sqrt(diag(Cov_mu_hat + Var_11 + Cov_12 + Cov_21 + Var_22))

      Psi_zeta_1 <- as.vector(Psi_1_hat %*% Zeta_1_hat[i, ])
      Psi_zeta_2 <- as.vector(Psi_2_hat %*% Zeta_2_hat[[i]][j, ])
      Y_hat <- mu_hat + Psi_zeta_1 + Psi_zeta_2

      Y_summary[[i]][[j]] <- matrix(NA, nrow = n_g, ncol = 3)
      Y_summary[[i]][[j]][, 1] <- Y_hat + qnorm(0.025)*sd_vec
      Y_summary[[i]][[j]][, 2] <- Y_hat
      Y_summary[[i]][[j]][, 3] <- Y_hat + qnorm(0.975)*sd_vec
    }
  }

  gbl_mat <- cbind(mu_hat, Psi_1_hat, Psi_2_hat)

  zeta_1_summary <- vector("list", length = N)
  zeta_2_summary <- vector("list", length = N)
  for(i in 1:N) {

    zeta_1_mean <- Zeta_1_hat[i, ][1:2]

    zeta_1_ellipse <- ellipse(
      Cov_zeta_1_hat[[i]][1:2, 1:2],
      centre = zeta_1_mean,
      level = 0.95
    )

    zeta_1_summary[[i]] <- list(zeta_1_mean, zeta_1_ellipse)
    names(zeta_1_summary[[i]]) <- c("mean", "cb")

    zeta_2_summary[[i]] <- vector("list", length = p)
    for(j in 1:p) {

      zeta_2_mean <- Zeta_2_hat[[i]][j, ]

      zeta_2_ellipse <- ellipse(
        Cov_zeta_2_hat[[i]][[j]][1:2, 1:2],
        centre = zeta_2_mean,
        level = 0.95
      )

      zeta_2_summary[[i]][[j]] <- list(zeta_2_mean, zeta_2_ellipse)
      names(zeta_2_summary[[i]][[j]]) <- c("mean", "cb")
    }
  }

  output <- list(Y_summary, gbl_mat, zeta_1_summary, zeta_2_summary)
  names(output) <- c("Y_summary", "gbl_mat", "zeta_1_summary", "zeta_2_summary")
  return(output)
}

two_level_fpc_rotation <- function(eta_in, p, time_g, C_g, Psi_1_g, Psi_2_g) {

  # order of eta_in:
  # 1. p(nu|Sigma_nu) -> nu
  # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # 3. p(zeta) -> zeta
  # 4. p(Y|nu,zeta,sigsq_eps) -> zeta

  N <- ncol(eta_in[[3]])
  L_1 <- ncol(Psi_1_g)
  L_2 <- ncol(Psi_2_g)
  L <- L_1 + L_2
  l_eta_nu <- length(eta_in[[1]])
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/(L+1) - 2

  # Compute q(nu):

  eta_nu <- list(eta_in[[1]], eta_in[[2]])
  q_nu <- gauss_q(eta_nu, use_vech = FALSE)
  E_q_nu <- q_nu[[1]]
  Cov_q_nu <- q_nu[[2]]

  # Compute q(zeta):

  E_q_zeta_1 <- vector("list", length = N)
  Cov_q_zeta_1 <- vector("list", length = N)
  E_q_zeta_2 <- vector("list", length = N)
  Cov_q_zeta_2 <- vector("list", length = N)
  Cov_q_zeta_12 <- vector("list", length = N)
  for(i in 1:N) {

    eta_in_zeta <- list(eta_in[[3]][, i], eta_in[[4]][, i])
    q_zeta <- two_level_nat_to_comm_parms(L_1, L_2, M, eta_in_zeta)
    E_q_zeta_1[[i]] <- q_zeta[[1]]
    Cov_q_zeta_1[[i]] <- q_zeta[[2]]
    E_q_zeta_2[[i]] <- q_zeta[[3]]
    Cov_q_zeta_2[[i]] <- q_zeta[[4]]
    Cov_q_zeta_12[[i]] <- q_zeta[[5]]
  }

  # Store mean function and eigenfunctions:

  E_q_V <- matrix(E_q_nu, K + 2, L + 1)
  gbl_mat <- C_g %*% E_q_V
  E_q_mu <- gbl_mat[, 1]
  E_q_Psi_1 <- gbl_mat[, 1:L_1 + 1]
  E_q_Psi_2 <- gbl_mat[, 1:L_2 + L_1 + 1]

  # Store scores:

  one_MN <- rep(1, M*N)
  one_M <- rep(1, M)

  E_q_Zeta_1 <- vector("list", length = N)
  E_q_Zeta_2 <- vector("list", length = N)
  for(i in 1:N) {

    E_q_Zeta_1[[i]] <- tcrossprod(one_M, E_q_zeta_1[[i]])
    E_q_Zeta_2[[i]] <- Reduce(rbind, E_q_zeta_2[[i]])
  }
  E_q_Zeta_1 <- Reduce(rbind, E_q_Zeta_1)
  E_q_Zeta_2 <- Reduce(rbind, E_q_Zeta_2)

  # Orthogonalization and Rotation:

  E_q_Psi_1_svd <- svd(E_q_Psi_1)
  U_Psi_1 <- E_q_Psi_1_svd$u
  D_Psi_1 <- diag(E_q_Psi_1_svd$d)
  V_Psi_1 <- E_q_Psi_1_svd$v

  E_q_Psi_2_svd <- svd(E_q_Psi_2)
  U_Psi_2 <- E_q_Psi_2_svd$u
  D_Psi_2 <- diag(E_q_Psi_2_svd$d)
  V_Psi_2 <- E_q_Psi_2_svd$v

  Zeta_1_rotn <- E_q_Zeta_1 %*% V_Psi_1 %*% D_Psi_1
  m_zeta_1 <- apply(Zeta_1_rotn, 2, mean)
  Zeta_1_shift <- Zeta_1_rotn - tcrossprod(one_MN, m_zeta_1)

  Zeta_2_rotn <- E_q_Zeta_2 %*% V_Psi_2 %*% D_Psi_2
  m_zeta_2 <- apply(Zeta_2_rotn, 2, mean)
  Zeta_2_shift <- Zeta_2_rotn - tcrossprod(one_MN, m_zeta_2)

  unique_lev_1_inds <- seq(1, M*N, by = M)
  eigen_Zeta_1_shift <- eigen(cov(Zeta_1_shift[unique_lev_1_inds, ]))
  Q_1 <- eigen_Zeta_1_shift$vectors
  Lambda_1 <- diag(eigen_Zeta_1_shift$values + 1e-10)
  Lambda_1_inv <- diag(1/(eigen_Zeta_1_shift$values + 1e-10))
  S_1 <- Q_1 %*% sqrt(Lambda_1)
  S_1_inv <- tcrossprod(sqrt(Lambda_1_inv), Q_1)

  eigen_Zeta_2_shift <- eigen(cov(Zeta_2_shift))
  Q_2 <- eigen_Zeta_2_shift$vectors
  Lambda_2 <- diag(eigen_Zeta_2_shift$values + 1e-10)
  Lambda_2_inv <- diag(1/(eigen_Zeta_2_shift$values + 1e-10))
  S_2 <- Q_2 %*% sqrt(Lambda_2)
  S_2_inv <- tcrossprod(sqrt(Lambda_2_inv), Q_2)

  mu_hat <- E_q_mu + as.vector(U_Psi_1 %*% m_zeta_1 + U_Psi_2 %*% m_zeta_2)
  Cov_mu_hat <- tcrossprod(C_g %*% Cov_q_nu[1:(K + 2), 1:(K + 2)], C_g)

  Psi_1_hat <- U_Psi_1 %*% S_1
  Zeta_1_hat <- tcrossprod(Zeta_1_shift, S_1_inv)

  Psi_2_hat <- U_Psi_2 %*% S_2
  Zeta_2_hat <- tcrossprod(Zeta_2_shift, S_2_inv)

  norm_const_1 <- rep(NA, L_1)
  for(l in 1:L_1) {

    norm_const_1[l] <- sqrt(trapint(time_g, (Psi_1_hat[,l])^2))
    if(norm_const_1[l] != 0) {

      Psi_1_hat[,l] <- Psi_1_hat[,l]/norm_const_1[l]
      Zeta_1_hat[,l] <- norm_const_1[l]*Zeta_1_hat[,l]

      cprod_sign <- sign(cprod(Psi_1_hat[,l], Psi_1_g[,l]))
      if(cprod_sign == -1) {

        Psi_1_hat[,l] <- -Psi_1_hat[,l]
        Zeta_1_hat[,l] <- -Zeta_1_hat[,l]
      }
    }
  }
  Zeta_1_hat <- Zeta_1_hat[unique_lev_1_inds, ]

  norm_const_2 <- rep(NA, L_2)
  for(l in 1:L_2) {

    norm_const_2[l] <- sqrt(trapint(time_g, (Psi_2_hat[,l])^2))
    if(norm_const_2[l] != 0) {

      Psi_2_hat[,l] <- Psi_2_hat[,l]/norm_const_2[l]
      Zeta_2_hat[,l] <- norm_const_2[l]*Zeta_2_hat[,l]

      cprod_sign <- sign(cprod(Psi_2_hat[,l], Psi_2_g[,l]))
      if(cprod_sign == -1) {

        Psi_2_hat[,l] <- -Psi_2_hat[,l]
        Zeta_2_hat[,l] <- -Zeta_2_hat[,l]
      }
    }
  }
  Zeta_2_hat <- lapply(split(Zeta_2_hat, rep(1:N, each = M)), function(x) matrix(x, M, L_2))

  # Update the posterior covariance matrix for each of the scores:

  scale_mat_1 <- diag(norm_const_1)
  scale_mat_2 <- diag(norm_const_2)
  transform_mat_1 <- scale_mat_1 %*% S_1_inv %*% tcrossprod(D_Psi_1, V_Psi_1)
  transform_mat_2 <- scale_mat_2 %*% S_2_inv %*% tcrossprod(D_Psi_2, V_Psi_2)
  Cov_zeta_1_hat <- vector("list", length = N)
  Cov_zeta_2_hat <- vector("list", length = N)
  Cov_zeta_12_hat <- vector("list", length = N)
  for(i in 1:N) {

    Cov_zeta_1_hat[[i]] <- tcrossprod(
      transform_mat_1 %*% Cov_q_zeta_1[[i]],
      transform_mat_1
    )

    Cov_zeta_2_hat[[i]] <- vector("list", length = M)
    Cov_zeta_12_hat[[i]] <- vector("list", length = M)
    for(j in 1:M) {

      Cov_zeta_2_hat[[i]][[j]] <- tcrossprod(
        transform_mat_2 %*% Cov_q_zeta_2[[i]][[j]],
        transform_mat_2
      )

      Cov_zeta_12_hat[[i]][[j]] <- tcrossprod(
        transform_mat_1 %*% Cov_q_zeta_12[[i]][[j]],
        transform_mat_2
      )
    }
  }

  # Summarise the results:

  Y_summary <- vector("list", length = N)
  for(i in 1:N) {

    Y_summary[[i]] <- vector("list", length = p)

    Var_11 <- tcrossprod(Psi_1_hat %*% Cov_zeta_1_hat[[i]], Psi_1_hat)

    for(j in 1:p) {

      Var_22 <- tcrossprod(Psi_2_hat %*% Cov_zeta_2_hat[[i]][[j]], Psi_2_hat)

      Cov_12 <- tcrossprod(Psi_1_hat %*% Cov_zeta_12_hat[[i]][[j]], Psi_2_hat)
      Cov_21 <- t(Cov_12)

      sd_vec <- sqrt(diag(Cov_mu_hat + Var_11 + Cov_12 + Cov_21 + Var_22))

      Psi_zeta_1 <- as.vector(Psi_1_hat %*% Zeta_1_hat[i, ])
      Psi_zeta_2 <- as.vector(Psi_2_hat %*% Zeta_2_hat[[i]][j, ])
      Y_hat <- mu_hat + Psi_zeta_1 + Psi_zeta_2

      Y_summary[[i]][[j]] <- matrix(NA, nrow = n_g, ncol = 3)
      Y_summary[[i]][[j]][, 1] <- Y_hat + qnorm(0.025)*sd_vec
      Y_summary[[i]][[j]][, 2] <- Y_hat
      Y_summary[[i]][[j]][, 3] <- Y_hat + qnorm(0.975)*sd_vec
    }
  }

  gbl_mat <- cbind(mu_hat, Psi_1_hat, Psi_2_hat)

  zeta_1_summary <- vector("list", length = N)
  zeta_2_summary <- vector("list", length = N)
  for(i in 1:N) {

    zeta_1_mean <- Zeta_1_hat[i, ][1:2]

    zeta_1_ellipse <- ellipse(
      Cov_zeta_1_hat[[i]][1:2, 1:2],
      centre = zeta_1_mean,
      level = 0.95
    )

    zeta_1_summary[[i]] <- list(zeta_1_mean, zeta_1_ellipse)
    names(zeta_1_summary[[i]]) <- c("mean", "cb")

    zeta_2_summary[[i]] <- vector("list", length = p)
    for(j in 1:p) {

      zeta_2_mean <- Zeta_2_hat[[i]][j, ]

      zeta_2_ellipse <- ellipse(
        Cov_zeta_2_hat[[i]][[j]][1:2, 1:2],
        centre = zeta_2_mean,
        level = 0.95
      )

      zeta_2_summary[[i]][[j]] <- list(zeta_2_mean, zeta_2_ellipse)
      names(zeta_2_summary[[i]][[j]]) <- c("mean", "cb")
    }
  }

  output <- list(Y_summary, gbl_mat, zeta_1_summary, zeta_2_summary)
  names(output) <- c("Y_summary", "gbl_mat", "zeta_1_summary", "zeta_2_summary")
  return(output)
}

##########################################
#
#  MISCELLANEOUS  FUNCTIONS
#
##########################################

is_int <- function(x, tol = .Machine$double.eps^0.5) {

  abs(x - round(x)) < tol
}

tr <- function(X) {

  if(nrow(X)!=ncol(X)) stop("X must be a square matrix.")

  ans <- sum(diag(X))
  return(ans)
}

cprod <- function(x, y) {

  if(missing(y)) {

    if(!is.vector(x)) {

      stop("Use the crossprod function for matrix inner products")
    }

    y <- x
  }

  if(!is.vector(y) & !is.vector(x)) {

    stop("Use the crossprod function for matrix inner products")
  }

  ans <- as.vector(crossprod(x, y))
  return(ans)
}




