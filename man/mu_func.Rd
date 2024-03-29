% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_functions.R
\name{mu_func}
\alias{mu_func}
\title{Example of mean function.}
\usage{
mu_func(t, j = 1, alpha = 3)
}
\arguments{
\item{t}{Time vector for which the function is evaluated.}

\item{j}{Phase shift and sign swap; see definition of the function above.}

\item{alpha}{Amplitude; see definition of the function above.}
}
\value{
Function evaluated at \code{t}.
}
\description{
Function to provide an example of mean function that can be used, e.g., for
simulations: \code{(-1)^j*alpha*sin((2*pi+j)*t)}
}
