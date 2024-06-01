compute_mu_q_log <- function(kappa, lambda) { log(lambda/2) - digamma(kappa/2)}

e_y <- function(N, sum_T, Y, C, M_q_V, M_q_H_tilde, mu_q_zeta, Sigma_q_zeta,
                mu_q_log_sigsq_eps, mu_q_recip_sigsq_eps) {

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

  term_1 <- - sum_T/2*log(2*pi)
  term_2 <- - sum_T/2*mu_q_log_sigsq_eps
  term_3 <- - 1/2*mu_q_recip_sigsq_eps*E_q_resid

  term_1 + term_2 + term_3

}

e_nu <- function(K, L, mu_beta, # hyperparamter
                 inv_Sigma_beta, # hyperparameter
                 mu_q_nu,
                 Sigma_q_nu,
                 inv_Sigma_q_nu,
                 mu_q_recip_sigsq_mu,
                 mu_q_recip_sigsq_psi,
                 mu_q_log_sigsq_mu,
                 mu_q_log_sigsq_psi) {

  # - c_ent_nu
  #
  M_q_inv_Sigma_m <- blkdiag(inv_Sigma_beta, mu_q_recip_sigsq_mu*diag(K))
  mu_m <- c(mu_beta, rep(0, K))

  M_q_inv_Sigma_p <- vector("list", length=L)
  mu_p <- vector("list", length=L)
  for(l in 1:L) {

    M_q_inv_Sigma_p[[l]] <- blkdiag(
      inv_Sigma_beta,
      mu_q_recip_sigsq_psi[l]*diag(K)
    )

    mu_p[[l]] <- c(mu_beta, rep(0, K))
  }
  M_q_inv_Sigma_nu <- blkdiag(M_q_inv_Sigma_m, Reduce(blkdiag, M_q_inv_Sigma_p))
  mu_nu <- c(mu_m, Reduce("c", mu_p))

  d <- (L+1)*(K+2)
  M_q_nu_nuT <- Sigma_q_nu + tcrossprod(mu_q_nu)

  # term_1 <- - d/2*log(2*pi) will cancel out with term_8
  term_2 <- (L+1)/2*determinant(inv_Sigma_beta, logarithm=TRUE)$modulus[1]
  term_3 <- - K/2*mu_q_log_sigsq_mu
  term_4 <- - K/2*sum(mu_q_log_sigsq_psi)
  term_5 <- - tr(M_q_nu_nuT%*%M_q_inv_Sigma_nu)/2
  term_6 <- cprod(mu_q_nu, M_q_inv_Sigma_nu%*%mu_nu)
  term_7 <- - cprod(mu_nu, M_q_inv_Sigma_nu%*%mu_nu)/2

  E_log_p_nu_rest <- # term_1 +
    term_2 + term_3 + term_4 + term_5 + term_6 + term_7

  # - ent_nu
  #
  # term_8 <- - d/2*log(2*pi)
  term_9 <- 1/2*determinant(inv_Sigma_q_nu, logarithm=TRUE)$modulus[1]
  term_10 <- - tr(M_q_nu_nuT%*%inv_Sigma_q_nu)/2
  term_11 <- cprod(mu_q_nu, inv_Sigma_q_nu %*% mu_q_nu)/2

  E_log_q_nu <- # term_8 +
    term_9 + term_10 + term_11

  E_log_p_nu_rest -  E_log_q_nu

}


e_zeta <- function(N, inv_Sigma_zeta, # hyperparameter
                   mu_q_zeta, Sigma_q_zeta, inv_Sigma_q_zeta) {

  # - c_ent_zeta_i
  #
  E_log_p_zeta_rest <- 0
  E_log_q_zeta <- 0
  for(i in 1:N) {

    M_q_zeta_zetaT_i <- Sigma_q_zeta[[i]] + tcrossprod(mu_q_zeta[[i]])

    # - c_ent_zeta_i
    # term_1 <- - L/2*log(2*pi) # will cancel out with term_4
    term_2 <- 1/2*determinant(inv_Sigma_zeta, logarithm=TRUE)$modulus[1] # could avoid taking the determinant as diagonal matrix
    term_3 <- - tr(M_q_zeta_zetaT_i%*%inv_Sigma_zeta)/2

    E_log_p_zeta_rest <- E_log_p_zeta_rest + # term_1 +
      term_2 + term_3


    # - ent_zeta_i
    # term_4 <- - L/2*log(2*pi)
    term_5 <- 1/2*determinant(inv_Sigma_q_zeta[[i]], logarithm=TRUE)$modulus[1]
    term_6 <- - tr(M_q_zeta_zetaT_i%*%inv_Sigma_q_zeta[[i]])/2
    term_7 <- cprod(mu_q_zeta[[i]], inv_Sigma_q_zeta[[i]] %*% mu_q_zeta[[i]])/2

    E_log_q_zeta <- E_log_q_zeta + # term_4
      term_5 + term_6 + term_7

  }

  E_log_p_zeta_rest -  E_log_q_zeta

}


e_sigsq <- function(mu_q_recip_sigsq, mu_q_log_sigsq, mu_q_recip_a, mu_q_log_a,
                    kappa_vb, lambda_vb) {

  # - c_ent_sigsq
  #
  term_1 <- - (log(2) + mu_q_log_a) / 2
  term_2 <- - log(pi) / 2
  term_3 <- - mu_q_recip_sigsq * mu_q_recip_a / 2
  term_4 <- - 3 * mu_q_log_sigsq / 2

  E_log_p_sigsq_rest <- term_1 + term_2 + term_3 + term_4

  # - ent_sigsq
  #
  # q(sigma) ~ Inv-X2(alpha_vb, beta_vb)
  #
  # E_log_q_sigsq <- alpha_vb/2 * log(alpha_vb*beta_vb/2) - lgamma(alpha_vb/2) -
  #  alpha_vb*beta_vb*mu_q_recip_sigsq/2 - (1 + alpha_vb/2) * mu_q_log_sigsq
  #
  # alpha_vb <- kappa_vb; beta_vb <- lambda_vb / kappa_vb
  # E_log_q_sigsq <- kappa_vb/2 * log(lambda_vb /2) - lgamma(kappa_vb/2) -
  #  lambda_vb*mu_q_recip_sigsq/2 - (1 + kappa_vb/2) * mu_q_log_sigsq

  term_5 <- kappa_vb/2 * (log(lambda_vb) - log(2))
  term_6 <- - lgamma(kappa_vb/2)
  term_7 <- - lambda_vb*mu_q_recip_sigsq / 2
  term_8 <- - (1 + kappa_vb/2) * mu_q_log_sigsq

  E_log_q_sigsq <- term_5 + term_6 + term_7 + term_8

  E_log_p_sigsq_rest - E_log_q_sigsq

}

e_a <- function(mu_q_recip_a, mu_q_log_a, A, kappa_vb, lambda_vb) {

  # - c_ent_sigsq
  #
  term_1 <- - (log(2) + 2*log(A)) / 2
  term_2 <- - log(pi) / 2
  term_3 <- - mu_q_recip_a / (2*A^2)
  term_4 <- - 3 * mu_q_log_a / 2

  E_log_p_a_rest <- term_1 + term_2 + term_3 + term_4

  # - ent_sigsq
  #
  # q(sigma) ~ Inv-X2(alpha_vb, beta_vb)
  #
  # E_log_q_sigsq <- alpha_vb/2 * log(alpha_vb*beta_vb/2) - lgamma(alpha_vb/2) -
  #  alpha_vb*beta_vb*mu_q_recip_sigsq/2 - (1 + alpha_vb/2) * mu_q_log_sigsq
  #
  # alpha_vb <- kappa_vb; beta_vb <- lambda_vb / kappa_vb
  # E_log_q_sigsq <- kappa_vb/2 * log(lambda_vb /2) - lgamma(kappa_vb/2) -
  #  lambda_vb*mu_q_recip_sigsq/2 - (1 + kappa_vb/2) * mu_q_log_sigsq

  term_5 <- kappa_vb/2 * (log(lambda_vb) - log(2))
  term_6 <- - lgamma(kappa_vb/2)
  term_7 <- - lambda_vb*mu_q_recip_a / 2
  term_8 <- - (1 + kappa_vb/2) * mu_q_log_a

  E_log_q_a <- term_5 + term_6 + term_7 + term_8

  E_log_p_a_rest - E_log_q_a

}







