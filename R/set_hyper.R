#' @export
set_hyper <- function(sigma_zeta = 1, mu_beta = 0, sigma_beta = 1e5, A = 1e5) {

  check_structure(sigma_zeta, "vector", "double", 1)
  check_positive(sigma_zeta)

  check_structure(mu_beta, "vector", "double", c(1, 2))
  if (length(mu_beta) == 1) mu_beta <- rep(mu_beta, 2)

  check_structure(sigma_beta, "vector", "double", c(1, 2))
  check_positive(sigma_beta)

  if (length(sigma_beta) == 1) sigma_beta <- rep(sigma_beta, 2)
  Sigma_beta <- sigma_beta^2*diag(2)

  check_structure(A, "vector", "double", 1)
  check_positive(A)

  list_hyper <- create_named_list(sigma_zeta, mu_beta, Sigma_beta, A)

  class(list_hyper) <- "hyper"

  list_hyper

}
