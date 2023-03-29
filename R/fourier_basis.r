########## R library: fourier_basis ##########

# For constructing a Fourier basis

# Created: 12 JUL 2022
# Last Updated: 12 JUL 2022

gen_sin_bf <- function(n) {
	
	sin_bf <- function(x) {
		
		ans <- sqrt(2)*sin(2*n*pi*x)
		return(ans)
	}
	
	return(sin_bf)
}

gen_cos_bf <- function(n) {
	
	cos_bf <- function(x) {
		
		ans <- sqrt(2)*cos(2*n*pi*x)
		return(ans)
	}
	
	return(cos_bf)
}

fourier_basis <- function(L) {
	
	L_even <- ((L/2) %% 1 == 0)
	
	if(L_even) {
		
		n <- 1:(L/2)
		sin_list <- lapply(n, gen_sin_bf)
		cos_list <- lapply(n, gen_cos_bf)
		fb <- do.call(c, Map(list, sin_list, cos_list))
	} else {
		
		n_sin <- 1:ceiling(L/2)
		sin_list <- lapply(n_sin, gen_sin_bf)
		
		if(floor(L/2) > 0) {
			
			n_cos <- 1:floor(L/2)
			cos_list <- lapply(n_cos, gen_cos_bf)
			
			fb <- vector("list", length = L)
			fb[c(TRUE, FALSE)] <- sin_list
			fb[c(FALSE, TRUE)] <- cos_list
		} else {
			
			fb <- sin_list
		}
	}
	
	return(fb)
}

########## End of fourier_basis ##########