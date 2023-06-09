% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cochran.R
\name{get_q_vector}
\alias{get_q_vector}
\title{get a vecotr of q values
apply q function to all of 1:n_permutations permutations}
\usage{
get_q_vector(.f, type, n_samples, n_permutations)
}
\arguments{
\item{.f}{get q function that takes one permutation number as input}

\item{type}{either "exact" or "montecarlo"}

\item{n_samples}{number of draws for montecarlo method}

\item{n_permutations}{number of possible permutations}
}
\value{
vector of q values
}
\description{
get a vecotr of q values
apply q function to all of 1:n_permutations permutations
}
\examples{
# .f should be a function created within the cochran_test function
# It should take only the permutation number as input
get_q_vector(.f = \(x){1:10[x]}, n_samples = 1000,  n_permutations = 10)


}
