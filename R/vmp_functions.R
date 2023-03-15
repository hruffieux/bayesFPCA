# A library that stores VMP fragment and sufficient
# statistic expectation updates:

# LIST  OF  FUNCTIONS:

# Fragment Updates:
# iter_igw_frag
# gauss_prior_frag
# gauss_pen_frag
# fpc_gauss_pen_frag
# fpc_lik_frag
# mfpc_lik_frag

# Natural Parameter to Common Parameter Updates and
# Expectatoin of Sufficient Statistic Computations:
# igw_q
# gauss_q

# Entropy Functions:
# entropy_igw
# entropy_gauss

# Cross-Entropy Functions:
# cross_entropy_igw_prior
# cross_entropy_iter_igw
# cross_entropy_gauss_prior
# cross_entropy_fpc_gauss_pen
# cross_entropy_fpc_lik_frag
# cross_entropy_mfpc_lik_frag

# FPCA functions:
# fpc_rotation
# mfpc_rotation

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

gauss_pen_frag <- function(eta_in, G_in, sigsq_beta) {

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

  eta_q_nu <- eta_in[[1]] + eta_in[[2]]
  eta_q_sigsq_nu <- eta_in[[3]] + eta_in[[4]]

  l_eta_nu <- length(eta_q_nu)
  d <- (-3 + sqrt(9 + 8* l_eta_nu))/2
  K <- d - 2
  D_nu <- duplication.matrix(K+2)

  eta_in_nu <- list(eta_in[[1]], eta_in[[2]])
  q_nu <- gauss_q(eta_in_nu, use_vech=TRUE)
  mu_q_nu <- q_nu[[1]]
  Sigma_q_nu <- q_nu[[2]]

  eta_in_sigsq_nu <- list(eta_in[[3]], eta_in[[4]])
  G_in_sigsq_nu <- c(G_in[[1]], G_in[[2]])
  q_sigsq_nu <- igw_q(eta_in_sigsq_nu, G_in_sigsq_nu)
  mu_q_recip_sigsq_nu <- q_sigsq_nu[[3]]
  E_q_inv_Sigma_nu <- adiag(1/sigsq_beta*diag(2), mu_q_recip_sigsq_nu*diag(K))

  eta_1 <- rep(0, K+2)
  eta_2 <- -0.5*cprod(D_nu, as.vector(E_q_inv_Sigma_nu))
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


mfpc_lik_frag <- function(n, eta_in, G_in, C, Y, L) {

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
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/(L+1) - 2
  D_L <- duplication.matrix(L)

  YTY <- matrix(NA, N, p)
  CTY <- vector("list", length = N)
  CTC <- vector("list", length = N)
  for(i in 1:N) {

    CTY[[i]] <- matrix(NA, K + 2, p)
    CTC[[i]] <- vector("list", length = p)
    for(j in 1:p) {

      YTY[i, j] <- cprod(Y[[i]][[j]])
      CTY[[i]][, j] <- cprod(Y[[i]][[j]], C[[i]][[j]])
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

  inds_mat <- matrix(1:d, K + 2, L + 1)
  mu_inds <- inds_mat[, 1]
  psi_inds <- as.vector(inds_mat[, 1:L + 1])

  E_q_V <- vector("list", length = p)
  E_q_nu_mu <- vector("list", length = p)
  E_q_V_psi <- vector("list", length = p)
  E_q_nu_psi <- vector("list", length = p)
  Cov_q_nu_mu <- vector("list", length = p)
  Cov_q_nu_psi <- vector("list", length = p)
  Cov_q_nu_mu_psi <- vector("list", length = p)
  for(j in 1:p) {

    E_q_V[[j]] <- matrix(E_q_nu[[j]], K + 2, L + 1)
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

      freq_scores <- cprod(E_q_V_psi[[j]], CTY[[i]][, j]) - E_q_h_mu_psi[[i]][, j]
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
      summands[2] <- - 2*cprod(E_q_nu_mu[[j]], CTY[[i]][, j])
      summands[3] <- -2*cprod(E_q_V_psi[[j]] %*% E_q_zeta[[i]], CTY[[i]][, j])
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
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/(L+1) - 2
  D_L <- duplication.matrix(L)
  n <- Reduce(rbind, lapply(Y, function(x) sapply(x, length)))

  YTY <- matrix(NA, N, p)
  CTY <- vector("list", length = N)
  CTC <- vector("list", length = N)
  for(i in 1:N) {

    CTY[[i]] <- matrix(NA, K + 2, p)
    CTC[[i]] <- vector("list", length = p)
    for(j in 1:p) {

      YTY[i, j] <- cprod(Y[[i]][[j]])
      CTY[[i]][, j] <- cprod(Y[[i]][[j]], C[[i]][[j]])
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

  inds_mat <- matrix(1:d, K + 2, L + 1)
  mu_inds <- inds_mat[, 1]
  psi_inds <- as.vector(inds_mat[, 1:L + 1])

  E_q_V <- vector("list", length = p)
  E_q_nu_mu <- vector("list", length = p)
  E_q_V_psi <- vector("list", length = p)
  E_q_nu_psi <- vector("list", length = p)
  Cov_q_nu_mu <- vector("list", length = p)
  Cov_q_nu_psi <- vector("list", length = p)
  Cov_q_nu_mu_psi <- vector("list", length = p)
  for(j in 1:p) {

    E_q_V[[j]] <- matrix(E_q_nu[[j]], K + 2, L + 1)
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
      summands[2] <- - 2*cprod(E_q_nu_mu[[j]], CTY[[i]][, j])
      summands[3] <- -2*cprod(E_q_V_psi[[j]] %*% E_q_zeta[[i]], CTY[[i]][, j])
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


##########################################
#
#  FPCA  FUNCTIONS
#
##########################################


fpc_rotation <- function(subj_names, L, eta_in, time_g, C_g, Psi_g = NULL) {

  # order of eta_in:
  # 1. p(nu|Sigma_nu) -> nu
  # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # 3. p(zeta) -> zeta
  # 4. p(Y|nu,zeta,sigsq_eps) -> zeta

  N <- ncol(eta_in[[3]])
  l_eta_nu <- length(eta_in[[1]])
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/(L+1) - 2
  n_g <- length(time_g)

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
  mu_Zeta_rotn <- apply(M_q_Zeta_rotn, 2, mean)

  mu_q_mu <- mu_q_mu + U_orth%*%mu_Zeta_rotn
  M_q_Zeta_shift <- M_q_Zeta_rotn - tcrossprod(one_N, mu_Zeta_rotn)

  eigen_M_q_Zeta_shift <- eigen(cov(M_q_Zeta_shift))
  Q <- eigen_M_q_Zeta_shift$vectors
  Lambda <- diag(eigen_M_q_Zeta_shift$values + 1e-10)
  Lambda_inv <- diag(1/(eigen_M_q_Zeta_shift$values + 1e-10))
  S <- Q%*%sqrt(Lambda)
  S_inv <- tcrossprod(sqrt(Lambda_inv), Q)

  Psi_hat <- U_orth%*%S
  Zeta_hat <- tcrossprod(M_q_Zeta_shift, S_inv)

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

  # Update the posterior covariance matrix for each zeta_i:

  scale_mat <- diag(norm_const)
  for(i in 1:N) {

    mat_transform <- S_inv%*%tcrossprod(D_diag, V_orth)
    Sigma_q_zeta[[i]] <- tcrossprod(mat_transform%*%Sigma_q_zeta[[i]], mat_transform)
    Sigma_q_zeta[[i]] <- tcrossprod(scale_mat%*%Sigma_q_zeta[[i]], scale_mat)
  }

  # Summarise the VMP results:

  Y_hat <- Y_low <- Y_upp <- vector("list", length = N)
  for(i in 1:N) {

      sd_vec_i <- sqrt(diag(tcrossprod(Psi_hat%*%Sigma_q_zeta[[i]], Psi_hat)))

      Y_hat[[i]] <- Y_mat[,i]
      Y_low[[i]] <- Y_mat[,i] + qnorm(0.025)*sd_vec_i
      Y_upp[[i]] <- Y_mat[,i] + qnorm(0.975)*sd_vec_i

  }

  names(Y_hat) <- names(Y_low) <- names(Y_upp) <- subj_names


  list_zeta_ellipse <- vector("list", length=N)

  for(i in 1:N) {

    zeta_mean <- Zeta_hat[i,][1:2]

    zeta_ellipse <- ellipse(
      Sigma_q_zeta[[i]][1:2, 1:2],
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

  outputs <- create_named_list(time_g,
                               Y_hat, Y_low, Y_upp,
                               mu_hat, list_Psi_hat,
                               Zeta_hat, list_zeta_ellipse)

}

mfpc_rotation <- function(subj_names, var_names, L, eta_in, time_g, C_g, Psi_g = NULL) {

  # order of eta_in:
  # 1. p(nu|Sigma_nu) -> nu
  # 2. p(Y|nu,zeta,sigsq_eps) -> nu
  # 3. p(zeta) -> zeta
  # 4. p(Y|nu,zeta,sigsq_eps) -> zeta

  N <- ncol(eta_in[[3]])
  p <- length(eta_in[[1]])
  l_eta_nu <- length(eta_in[[1]][[1]])
  d <- (sqrt(4*l_eta_nu + 1) - 1)/2
  K <- d/(L+1) - 2 ######################################## <-----------------
  D_L <- duplication.matrix(L)
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
    E_q_V <- matrix(E_q_nu, K + 2, L + 1) ######################################## <-----------------

    gbl_post <- C_g %*% E_q_V
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
  # zeta_summary <-
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

  outputs <- create_named_list(time_g,
                               Y_hat, Y_low, Y_upp,
                               mu_hat, list_Psi_hat,
                               Zeta_hat, list_zeta_ellipse)

  return(outputs)
}
