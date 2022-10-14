# Obtain the vectorisation of a matrix
#
# Given an n x m matrix, returns an (nm) x 1 matrix containing the column stack
# of the matrix
# Input: A matrix to vectorize
# Output: one-column matrix; the vectorization of A
# See also vecInverse
#
vec <- function(A){
    return(as.vector(A))
}
