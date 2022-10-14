vecInverse <- function(a)
{
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
