#' Gather model hyperparameters.
#'
#' This function is used to construct an object with hyperparameter values for
#' the FPCA or mFPCA model in a format that can then be supplied to
#' \code{\link{run_mfvb_fpca}} and \code{\link{run_vmp_fpca}}.
#'
#' The \code{\link{run_mfvb_fpca}} and \code{\link{run_vmp_fpca}} functions can
#' also be used with default hyperparameter values (i.e., without using
#' \code{\link{set_hyper}}) by setting their argument \code{list_hyper} to
#' \code{NULL}.
#'
#' @param sigma_beta Vector of size 1 or 2 for the standard deviation of the
#'                   spline coefficients for the mean function.
#' @param A Positive real number for the top-level hyperparameter of the
#'          iterated inverse chi-square distributions.
#'
#' @return An object containing the hyperparameter settings to be supplied to
#'         \code{\link{run_mfvb_fpca}} or \code{\link{run_vmp_fpca}}
#'
#' @export
#'
set_hyper <- function(sigma_beta = 1e5, A = 1e5) {

  # check_structure(sigma_zeta, "vector", "double", 1)
  # check_positive(sigma_zeta)
  sigma_zeta <- 1  # fixed to 1. PVE regained through post-hoc orthonormalisation

  # Now set to zero everywhere (mfpca version MFVB assumed to be zero in the algorithm)
  # Vector of size 1 or 2 for the mean of the spline coefficients
  #                for the mean function.
  # check_structure(mu_beta, "vector", "double", c(1, 2))
  # if (length(mu_beta) == 1) mu_beta <- rep(mu_beta, 2)

  check_structure(sigma_beta, "vector", "double", c(1, 2))
  check_positive(sigma_beta)

  if (length(sigma_beta) == 1) sigma_beta <- rep(sigma_beta, 2)
  Sigma_beta <- sigma_beta^2*diag(2)

  check_structure(A, "vector", "double", 1)
  check_positive(A)

  list_hyper <- create_named_list(sigma_zeta, Sigma_beta, A)

  class(list_hyper) <- "hyper"

  list_hyper

}
