create_named_list <- function(...) {
  setNames(list(...), as.character(match.call()[-1]))
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


#' @export
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

#' @export
trapint <- function(xgrid,fgrid) {
  ng <- length(xgrid)

  xvec <- xgrid[2:ng] - xgrid[1:(ng-1)]
  fvec <- fgrid[1:(ng-1)] + fgrid[2:ng]

  integ <- sum(xvec*fvec)/2

  return(integ)
}

