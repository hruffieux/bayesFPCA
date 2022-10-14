# For computing the trace of a square matrix

tr <- function(X) {

	if(nrow(X)!=ncol(X)) stop("X must be a square matrix.")

	ans <- sum(diag(X))
	return(ans)
}
