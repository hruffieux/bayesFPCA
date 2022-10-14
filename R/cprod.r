# For constructing a vector version of the crossprod function

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
