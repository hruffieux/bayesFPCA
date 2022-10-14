# For constructing a linear regression-based design matrix

# x must be a vector or a list of vectors

X_design <- function(x) {

	if(is.list(x)) {

		x <- do.call(cbind, x)
	}

	X <- cbind(1, x)
	return(X)
}
