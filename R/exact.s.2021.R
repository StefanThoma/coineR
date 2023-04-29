
# Example data.

# Association of workload (FTE - Full Time Equivalents in %) and
# perceived energy (1-7) of employees

fte <- c(100, 100, 90, 90, 50, 60)
energy <- c(3, 3, 2, 3, 6, 5)  # S = -7 / two-sided p = 0.15


# This function takes the vectors x and y, and computes the corresponding s value
get.s <- function(x, y){
  table.temp <- table(x,y)
  cd <- DescTools::ConDisPairs(table.temp)
  return(cd$C-cd$D)
  }

# This function takes the vectors x and y and computes the distribution of s under the assumption that the x and y vectors are independent.
# The output is a list of results, containing the counts of the s vectors, and the corresponding probabilities, and the p-value for the input vectors.
# This function does not compute the precise distribution under H0 but simulates a distribution using n.sim simulations.


get.distr.montecarlo <- function(x, y, n.sim = 10000){

  # get the length of the input vector. should be the same length for x and y
  vec.length <- length(x)


  # get the s-vector under h0
  s.vector <- replicate(n = n.sim, expr = get.s(x, y[sample(x = 1:vec.length, size  = vec.length, replace = F)]))
  # get the number of simulations (already defined)
  perm.length <- n.sim

  # get the count table
  s.count.table <- table(s.vector)
  # get the probability mass table
  ps <- s.count.table/perm.length

  # get the s value for the original input
  s.value <- get.s(x,y)

  # the probability to get such a result or more extreme under h0:
  possible.s.values <- names(s.count.table)
  which.is.it <- which(possible.s.values == paste(s.value))

  # get the p value
  p.test <- c(sum(s.count.table[1:which.is.it])/perm.length, sum(s.count.table[which.is.it:length(s.count.table)])/perm.length)

  res.list <- list("counts" = table(s.vector), "probabilities" = ps, "result" = ps[paste(s.value)], "p.values" = p.test)

  return(res.list)
}

# This function takes the vectors x and y and computes the exact distribution
# of s under the assumption that the x and y vectors are independent.
# The output is a list of results, containing the counts of the s vectors,
# the corresponding probabilities, as well as the point probability and both
# one-sided p-values for the resulting s value.
# Caution: Not recommended for n > 8 (will take much too long...)

get.distr.exact <- function(x, y){
  vec.length <- length(x)

  # get the permutation table. Its length is factorial(vec.length)
  perm.table <- gtools::permutations(n = vec.length, r = vec.length)
  perm.length <- nrow(perm.table)

  # get precise s values under h0
  s.vector <- apply(perm.table, MARGIN = 1, FUN = function(ind){get.s(x, y[ind])})


  # Alternative computation with for loop (iterate through all permutations
  # in permutation table and input new y shuffle into get.s function)
  # s.vector <- numeric(perm.length)
  # for(i in 1:perm.length){
  #  y.new <- y[perm.table[i,]]
  #  s.vector[i] <- get.s(x = x, y = y.new)
  #}



  # get the count table for s
  s.count.table <- table(s.vector)

  # get the probability mass table for s
  ps <- s.count.table/perm.length

  # get the s value for the original input
  s.value <- get.s(x,y)

  # the probability to get such a result or more extreme under h0:
  possible.s.values <- names(s.count.table)
  which.is.it <- which(possible.s.values == paste(s.value))

  # get the p values for one-sided tests (alternative = less) and  (alternative = greater)
  p.test <- c(sum(s.count.table[1:which.is.it])/perm.length, sum(s.count.table[which.is.it:length(s.count.table)])/perm.length)

  res.list <- list("vector" = s.vector, "counts" = table(s.vector), "probabilities" = ps, "result" =
                       ps[paste(s.value)], "p.values" = p.test)


  cat(c("\n", "Exact Test for Kendalls S", "\n", "\n"))
  cat(c(" S =", s.value, "\n",
        "p-values for H1 neg./pos. association:", p.test))
  invisible(res.list)
}




get.distr.exact.or.montecarlo.wo.replacement <- function(x, y, n.sim = 10000){
  vec.length <- length(x)

  # get the permutation table. Its length is factorial(vec.length)
  perm.table <- gtools::permutations(n = vec.length, r = vec.length)


  # get all combinations or a random subset of n.sim combinations from perm.table
  sample.perm.table <- perm.table[sample(nrow(perm.table),
                                         size = min(c(nrow(perm.table), n.sim)),
                                         replace = FALSE), ]
  sample.perm.length <- nrow(sample.perm.table)

  # get precise s values or a subset of n.sim s values under h0
  s.vector <- apply(sample.perm.table, MARGIN = 1, FUN = function(ind){get.s(x, y[ind])})


  # get the count table
  s.count.table <- table(s.vector)
  # get the probability mass table
  ps <- s.count.table/sample.perm.length

  # get the s value for the original input
  s.value <- get.s(x,y)

  # the probability to get such a result or more extreme under h0:
  possible.s.values <- names(s.count.table)
  which.is.it <- which(possible.s.values == paste(s.value))

  # get the p value
  p.test <- c(sum(s.count.table[1:which.is.it])/sample.perm.length, sum(s.count.table[which.is.it:length(s.count.table)])/sample.perm.length)

  res.list <- list("counts" = table(s.vector),
                   "probabilities" = ps,
                   "result" = ps[paste(s.value)],
                   "p.values" = p.test,
                   "method" = ifelse(sample.perm.length < nrow(perm.table),
                                     paste("montecarlo with", n.sim, "simulations without replacement", sep = " "),
                                     paste("exact")))

  return(res.list)
}


## some tests for get.distr.exact.or.montecarlo.wo.replacement for n = 6:8
Angst_neu6 <- c(75, 75, 77, 77, 78, 89)
IQ_neu6 <- c(100, 130, 100, 120, 130, 140)

# exact test (factorial(6) = 720 permutations)
result_neu6_exact <- get.distr.exact.or.montecarlo.wo.replacement(x = Angst_neu6, y = IQ_neu6)

# montecarlo test with 500 permutations
result_neu6_mc <- get.distr.exact.or.montecarlo.wo.replacement(x = Angst_neu6, y = IQ_neu6, n.sim = 500)

plot(result_neu6_exact$counts)
plot(result_neu6_mc$counts)

result_neu6_exact[c("p.values", "method")]
result_neu6_mc[c("p.values", "method")]



Angst_neu7 <- c(75, 75, 77, 77, 78, 89, 89)
IQ_neu7 <- c(100, 130, 100, 120, 130, 140, 140)

# exact test (factorial(7) = 5040 permutations)
result_neu7_exact <- get.distr.exact.or.montecarlo.wo.replacement(x = Angst_neu7, y = IQ_neu7)

# montecarlo test with 2500 permutations
result_neu7_mc <- get.distr.exact.or.montecarlo.wo.replacement(x = Angst_neu7, y = IQ_neu7, n.sim = 2500)

plot(result_neu7_exact$counts)
plot(result_neu7_mc$counts)

result_neu7_exact[c("p.values", "method")]
result_neu7_mc[c("p.values", "method")]



Angst_neu8 <- c(75, 75, 77, 77, 78, 89, 89, 65)
IQ_neu8 <- c(100, 130, 100, 120, 130, 140, 140, 120)

# exact test (factorial(8) = 40320 permutations)
result_neu8_exact <- get.distr.exact.or.montecarlo.wo.replacement(x = Angst_neu8, y = IQ_neu8, n.sim = 40320)

# montecarlo test with 10000 permutations
result_neu8_mc <- get.distr.exact.or.montecarlo.wo.replacement(x = Angst_neu8, y = IQ_neu8)

plot(result_neu8_exact$counts)
plot(result_neu8_mc$counts)

result_neu8_exact[c("p.values", "method")]
result_neu8_mc[c("p.values", "method")]







# compare the results for the two methods:
result <- get.distr.exact(x = Reaktionszeit, y = IQ)
result.monte <- get.distr.montecarlo(x = Reaktionszeit, n.sim = 10000)


plot(result$probabilities)
plot(result.monte$probabilities)
# Almost the same


# now lets see for larger vectors how they compare:
N <- 1:8

runtime <- data.frame(sample.size = rep(N, each = 2),
                      method = c("exact", "montecarlo"),
                      runtime = numeric(length(N)*2))

for(i in N){
  x.temp <- rnorm(i)
  y.temp <- rnorm(i)

  runtime$runtime[runtime$method == "exact" & runtime$sample.size == i] <- system.time({ get.distr.exact(x.temp, y.temp) })[3]
  runtime$runtime[runtime$method == "montecarlo" & runtime$sample.size == i] <- system.time({ get.distr.montecarlo(x.temp, y.temp) })[3]
  print(i)
}
library(ggplot2)
ggplot(data = runtime, aes(x = sample.size, y = runtime, color = method)) +
  geom_line()# +
  #scale_y_log10()

