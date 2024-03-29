create_named_list <- function(...) {
  setNames(list(...), as.character(match.call()[-1]))
}

all_same <- function(x) length(unique(x)) == 1

check_natural <- function(x, eps = .Machine$double.eps^0.75){
  if (any(x < eps | abs(x - round(x)) > eps)) {
    stop(paste0(deparse(substitute(x)),
                " must be natural."))
  }
}

check_positive <- function(x, eps = .Machine$double.eps^0.75){
  if (any(x < eps)) {
    err_mess <- paste0(deparse(substitute(x)), " must be positive (larger than precision zero in R).")
    if (length(x) > 1) err_mess <- paste0("All entries of ", err_mess)
    stop(err_mess)
  }
}

check_zero_one <- function(x){
  if (any(x < 0) | any(x > 1)) {
    err_mess <- paste0(deparse(substitute(x)), " must lie between 0 and 1.")
    if (length(x) > 1) err_mess <- paste0("All entries of ", err_mess)
    stop(err_mess)
  }
}

check_structure <- function(x, struct, type, size = NULL,
                             null_ok = FALSE,  inf_ok = FALSE, na_ok = FALSE) {
  if (type == "double") {
    bool_type <-  is.double(x)
    type_mess <- "a double-precision "
  } else if (type == "integer") {
    bool_type <- is.integer(x)
    type_mess <- "an integer "
  } else if (type == "numeric") {
    bool_type <- is.numeric(x)
    type_mess <- "a numeric "
  } else if (type == "logical") {
    bool_type <- is.logical(x)
    type_mess <- "a boolean "
  } else if (type == "string") {
    bool_type <- is.character(x)
    type_mess <- "string "
  }

  bool_size <- TRUE # for case size = NULL (no assertion on the size/dimension)
  size_mess <- ""
  if (struct == "vector") {
    bool_struct <- is.vector(x) & (length(x) > 0) # not an empty vector
    if (!is.null(size)) {
      bool_size <- length(x) %in% size
      size_mess <- paste0(" of length ", paste0(size, collapse=" or "))
    }
  } else if (struct == "matrix") {
    bool_struct <- is.matrix(x) & (length(x) > 0) # not an empty matrix
    if (!is.null(size)) {
      bool_size <- all(dim(x) == size)
      size_mess <- paste0(" of dimension ", size[1], " x ", size[2])
    }
  }

  correct_obj <- bool_struct & bool_type & bool_size

  bool_null <- is.null(x)

  if (!is.list(x) & type != "string") {
    na_mess <- ""
    if (!na_ok) {
      if (!bool_null) correct_obj <- correct_obj & !any(is.na(x))
      na_mess <- " without missing value"
    }

    inf_mess <- ""
    if (!inf_ok) {
      if (!bool_null) correct_obj <- correct_obj & all(is.finite(x[!is.na(x)]))
      inf_mess <- ", finite"
    }
  } else {
    na_mess <- ""
    inf_mess <- ""
  }

  null_mess <- ""
  if (null_ok) {
    correct_obj <- correct_obj | bool_null
    null_mess <- " or must be NULL"
  }

  if(!(correct_obj)) {
    stop(paste0(deparse(substitute(x)), " must be a non-empty ", type_mess, struct,
                size_mess, inf_mess, na_mess, null_mess, "."))
  }
}


#' Create a temporal-grid objects for use in FPCA algorithms.
#'
#' This function is used to create a time grid of desired density and design
#' matrices based on an observed time grid object.
#'
#' @param time_obs Vector or list of vectors containing time of observations for
#'                 univariate or multivariate curves, respectively.
#' @param K Number of O'Sulivan spline functions to be used in the FPCA
#'          algorithms. If set to \code{NULL} will be set according to the rule
#'          of Ruppert (2002), also enforcing K >=7.
#' @param n_g Desired size for dense grid.
#' @param time_g Dense grid provided as a vector of size \code{n_g}. If provided,
#'               then \code{n_g} must be \code{NULL} as will be taken to be
#'               \code{length(time_g)}.
#' @param int_knots Position of interior knots. Default is \code{NULL} for
#'                  evenly placed.
#' @param format_univ Boolean indicating whether the univariate format is used
#'                    in case of univariate curves.
#' @return C Design matrix C(t) constructed from the set of K spline functions
#'           based on the observation time grid.
#' @return n_g Number of time points in the dense grid.
#' @return time_g Dense time grid constructed.
#' @return C_g Design matrix C(t) constructed from the set of K spline functions
#'             based on the dense time grid.
#'
#' @export
#'
get_grid_objects <- function(time_obs, K, n_g = 1000, time_g = NULL,
                             int_knots = NULL,
                             format_univ = FALSE) {


  if (is.null(int_knots)) {
  	if(format_univ) {

  	  if (is.null(K)) {  # Ruppert (2002) sets a simple default value for K as min(nobs/4,40), where nobs is the number of observations.
                    	   # here since nobs differs for each i, we take the nobs / 4 = round(median(obs_i)/4)
                    	   # and we enforce that K>=7
  	    K <- max(round(min(median(sapply(time_obs, function(time_obs_i) length(time_obs_i))/4), 40)), 7)
  	  }

  	  unique_time_obs <- sort(unique(Reduce(c, time_obs)))
  	  int_knots <- quantile(unique_time_obs, seq(0, 1, length=K)[-c(1,K)])
  	} else {

  	  p <- length(time_obs[[1]])
  	  if (is.null(K)) {   # Ruppert (2002) sets a simple default value for K as min(nobs/4,40), where nobs is the number of observations.
  	                      # here since nobs differs for each i, we take the nobs / 4 = round(median(obs_i)/4)
  	                      # and we enforce that K>=7
  	    K <- sapply(1:p, function(j) max(round(min(median(sapply(time_obs, function(time_obs_i) length(time_obs_i[[j]]))/4), 40)), 7))
  	  } else {
  	    check_structure(K, "vector", "numeric", c(1, p))
  	    if (length(K)==1) K <- rep(K, p)
  	  }

      unique_time_obs <- unname(sort(unlist(time_obs)))
  	  int_knots <- lapply(
        K,
        function(x) quantile(unique_time_obs, seq(0, 1, length = x)[-c(1, x)])
      )
  	}
  }

  N <- length(time_obs)

  if (is.null(time_g)) {
    stopifnot(!is.null(n_g))
    time_g <- seq(0, 1, length.out = n_g)
  } else {
    if(is.null(n_g)) n_g <- length(time_g)
  }

  C <- vector("list", length=N)

  if (format_univ) {
    for(i in 1:N) {
      X <- X_design(time_obs[[i]])
      Z <- ZOSull(time_obs[[i]], range.x=c(0, 1), intKnots=int_knots)
      C[[i]] <- cbind(X, Z)
    }

    X_g <- X_design(time_g)
    Z_g <- ZOSull(time_g, range.x=c(0, 1), intKnots=int_knots)
    C_g <- cbind(X_g, Z_g)
  } else {

    p <- length(time_obs[[1]])

    for(i in 1:N) {

      C[[i]] <- vector("list", length = p)
      for(j in 1:p) {

        X <- X_design(time_obs[[i]][[j]])
        Z <- ZOSull(time_obs[[i]][[j]], range.x = c(0, 1), intKnots = int_knots[[j]])
        C[[i]][[j]] <- cbind(X, Z)
      }
    }

    X_g <- X_design(time_g)
    Z_g <- lapply(
      int_knots,
      function(x) ZOSull(time_g, range.x = c(0, 1), intKnots = x)
    )
    C_g <- lapply(Z_g, function(Z) cbind(X_g, Z))
  }

  create_named_list(C, n_g, time_g, C_g)
}


X_design <- function(x) {

  if(is.list(x)) {

    x <- do.call(cbind, x)
  }

  X <- cbind(1, x)
  return(X)
}


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


normalise <- function(x) {

  ans <- x/sqrt(cprod(x))
  return(ans)
}


E_cprod <- function(mean_1, Cov_21, mean_2, A) {

  if(missing(A)) {

    A <- diag(length(mean_1))
  }

  tr_term <- tr(Cov_21 %*% A)
  cprod_term <- cprod(mean_1, A %*% mean_2)
  ans <- tr_term + cprod_term
  return(ans)
}


E_h <- function(L, mean_1, Cov_21, mean_2, A) {

  if(missing(A)) {

    A <- diag(length(mean_1))
  }

  d_1 <- length(mean_1)
  inds_1 <- matrix(1:d_1, ncol = L)

  ans <- rep(NA, L)
  for(l in 1:L) {

    mean_1_l <- mean_1[inds_1[, l]]
    Cov_21_l <- Cov_21[, inds_1[, l]]
    ans[l] <- E_cprod(mean_1_l, Cov_21_l, mean_2, A)
  }

  return(ans)
}


E_H <- function(L_1, L_2, mean_1, Cov_21, mean_2, A) {

  if(missing(A)) {

    A <- diag(length(mean_1))
  }

  d_1 <- length(mean_1)
  d_2 <- length(mean_2)
  L <- L_1 + L_2

  inds_1 <- matrix(1:d_1, ncol = L_1)
  inds_2 <- matrix(1:d_2, ncol = L_2)

  ans <- matrix(NA, L_1, L_2)
  for(l in 1:L_1) {

    mean_1_l <- mean_1[inds_1[, l]]
    for(k in 1:L_2) {

      mean_2_k <- mean_2[inds_2[, k]]

      Cov_21_kl <- Cov_21[inds_2[, k], inds_1[, l]]

      ans[l, k] <- E_cprod(mean_1_l, Cov_21_kl, mean_2_k, A)
    }
  }

  return(ans)
}


vec <- function(A) {
  return(as.vector(A))
}


vecInverse <- function(a) {
  is.wholenumber <- function(x,tol=sqrt(.Machine$double.eps))
    return(abs(x-round(x))<tol)

  a <- as.vector(a)
  if (!is.wholenumber(sqrt(length(a))))
    stop("input vector must be a perfect square in length")

  dmnVal <- round(sqrt(length(a)))

  A <- matrix(NA,dmnVal,dmnVal)
  for (j in 1:dmnVal)
    A[,j] <- a[((j-1)*dmnVal+1):(j*dmnVal)]

  return(A)
}

#' Function performing trapesoidal integration.
#'
#' @param xgrid Grid.
#' @param fgrid Function on the grid.
#' @return Integration result.
#'
#' @export
#'
trapint <- function(xgrid,fgrid) {
  ng <- length(xgrid)

  xvec <- xgrid[2:ng] - xgrid[1:(ng-1)]
  fvec <- fgrid[1:(ng-1)] + fgrid[2:ng]

  integ <- sum(xvec*fvec)/2

  return(integ)
}

wait <- function() {
  cat("Hit Enter to continue\n")
  ans <- readline()
  invisible()
}
