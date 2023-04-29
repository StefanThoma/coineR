
# Example data.

# Association of workload (FTE - Full Time Equivalents in %) and
# perceived energy (1-7) of employees

fte <- c(100, 100, 90, 90, 50, 60)
energy <- c(3, 3, 2, 3, 6, 5)  # S = -7 / two-sided p = 0.1500
energy_t <- c(5, 6, 3, 2, 3, 3)  # S = 6 / two-sided p = 0.3167
energy_0 <- c(3, 6, 2, 3, 5, 3) # S = 0 /two-sided p = 1


# The function exact_kendall() takes the vectors x and y and computes a test based on
# the exact distribution of s for the H0 that the x and y vectors are independent.
# The output is a list of results, containing the counts of the s vectors,
# the corresponding probabilities, as well as the point probability and both
# one-sided and two-sided p-values for the resulting s value.
# Caution: Not recommended for n > 8 (will take much too long...

exact_kendall <- function(x, y){
    # define function to compute s
    get.s <- function(x, y){
        table.temp <- table(x,y)
        cd <- DescTools::ConDisPairs(table.temp)
        return(cd$C-cd$D)

    }
    # get n
    vec.length <- length(x)

    # get the permutation table. Its length is factorial(vec.length)
    perm.table <- gtools::permutations(n = vec.length, r = vec.length)
    perm.length <- nrow(perm.table)

    # get exact distribution of s under h0
    s.vector <- apply(perm.table, MARGIN = 1, FUN = function(ind){get.s(x, y[ind])})


    # Alternative computation with for loop (iterate through all permutations
    # in permutation table and input new y shuffle into get.s function)
    # s.vector <- numeric(perm.length)
    # for(i in 1:perm.length){
    #  y.new <- y[perm.table[i,]]
    #  s.vector[i] <- get.s(x = x, y = y.new)
    # }

    # get the count table for s
    s.count.table <- table(s.vector)

    # get the probability mass table for s
    ps <- s.count.table/perm.length

    # get the s value for the original input
    s.value <- get.s(x,y)

    # get the p values for one-sided tests (alternative = less) and (alternative = greater)
    p.1sided <- round(c(mean(s.vector <= s.value),
                        mean(s.vector >= s.value)), 4)

    # get the p value for the two-sided test
    p.2sided <- if (s.value == 0) 1L
    else {round(mean(s.vector <= if (s.value < 0) s.value else -s.value) +
                    mean(s.vector >= if (s.value > 0) s.value else -s.value), 4)}

    res.list <- list("vector" = s.vector,
                     "counts" = table(s.vector),
                     "probabilities" = ps,
                     "result" = ps[paste(s.value)],
                     "p.values" = c(p.1sided, p.2sided))


    cat(c("\n", "Exact Test for Kendalls S", "\n", "\n"))
    cat(c(" S =", s.value, "\n",
          "1-sided p-values for H1 neg./pos. association:", p.1sided, "\n",
          "2-sided p-value:", p.2sided))
    invisible(res.list)
}


# Some tests of the function
exact_kendall(fte, energy)
exact_kendall(fte, energy_0)
exact_kendall(fte, energy_t)



# The function exact_kendall_mc() takes the vectors x and y and computes a test based
# either A) on the exact distribution of s or B) on the montecarlo approximate
# distribution (with n simulations [without replacement]) of s for the H0
# that x and y are independent. The function will automatically compute the exact
# test when n_samples is larger than the number of permutations under H0, otherwise the
# montecarlo method is used.
# The output is a list of results, containing the counts of the s vectors,
# the corresponding probabilities, as well as the point probability and both
# one-sided and two-sided p-values for the resulting s value.

exact_kendall_mc <- function(x, y, n_samples = 10000){
    # define function to compute s
    get.s <- function(x, y){
        table.temp <- table(x,y)
        cd <- DescTools::ConDisPairs(table.temp)
        return(cd$C-cd$D)

    }
    # get n
    vec.length <- length(x)

    # get the permutation table. Its length is factorial(vec.length)
    perm.table <- gtools::permutations(n = vec.length, r = vec.length)


    # get all combinations or a random subset of n_samples combinations from perm.table
    sample.perm.table <- perm.table[sample(nrow(perm.table),
                                           size = min(c(nrow(perm.table), n_samples)),
                                           replace = FALSE), ]
    sample.perm.length <- nrow(sample.perm.table)

    # get precise s values or a subset of n_samples s values under h0
    s.vector <- apply(sample.perm.table, MARGIN = 1, FUN = function(ind){get.s(x, y[ind])})


    # get the count table
    s.count.table <- table(s.vector)
    # get the probability mass table
    ps <- s.count.table/sample.perm.length

    # get the s value for the original input
    s.value <- get.s(x,y)

    # get the p values for one-sided tests (alternative = less) and (alternative = greater)
    p.1sided <- round(c(mean(s.vector <= s.value),
                        mean(s.vector >= s.value)), 4)

    # get the p value for the two-sided test
    p.2sided <- if (s.value == 0) 1L
    else {round(mean(s.vector <= if (s.value < 0) s.value else -s.value) +
                    mean(s.vector >= if (s.value > 0) s.value else -s.value), 4)}

    res.list <- list("vector" = s.vector,
                     "counts" = table(s.vector),
                     "probabilities" = ps,
                     "result" = ps[paste(s.value)],
                     "p.values" = c(p.1sided, p.2sided),
                     "method" = ifelse(sample.perm.length < nrow(perm.table),
                                       paste("montecarlo with", n_samples,
                                             "simulations", sep = " "),
                                       paste("exact")))


    if(sample.perm.length < nrow(perm.table))
    {cat(c("\n", "Approximate Test for Kendalls S (Montecarlo with", n_samples, "simulations)", "\n", "\n"))}
    else {cat(c("\n", "Exact Test for Kendalls S", "\n", "\n"))}


    cat(c(" S =", s.value, "\n",
          "1-sided p-values for H1 neg./pos. association:", p.1sided, "\n",
          "2-sided p-value:", p.2sided))

    invisible(res.list)
}


# Some tests of the function
exact_kendall_mc (fte, energy, n_samples = 700)
exact_kendall_mc (fte, energy_t, n_samples = 10000)
exact_kendall_mc (fte, energy_0)
exact_kendall_mc (fte, energy_0, n_samples = 720)


