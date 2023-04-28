#' cochrans q test
#'
#' @param .data
#' @param formula
#' @param type
#'
#' @return
#' @export
#'
#' @examples
#'
autos <- matrix(data = c(0,0,0,0,0,1,0,0,0,1,
                         0,1,0,0,1,0,0,1,0,1,
                         0,1,1,0,1,1,0,1,1,1,
                         0,1,1,0,0,1,1,1,1,1),
                nrow = 10, ncol = 4, byrow = FALSE) %>% as_tibble() %>%
  mutate(family = row_number()) %>%
  pivot_longer(V1:V4, names_to = "group")
#'
formula <- value ~ group | family
.data <- autos
cochran_test <- function(.data, formula, type = c("exact", "montecarlo"), n_samples = 10000){
  # check inputs
  match.arg(type)
  stopifnot(rlang::is_scalar_integerish(n_samples))
  admiraldev::assert_atomic_vector(n_samples)
  # parse formula
  lhs <- rlang::f_lhs(formula)
  rhs <- rlang::f_rhs(formula)
  block <- rhs[[3]]
  group <- rhs[[2]]

  admiraldev::assert_data_frame(.data)


  compute_q <- create_q_f(.data = .data, formula = formula)


  group_ss <- {.data %>% group_by({{group}}) %>%
    summarize(blocksum = sum({{lhs}})) %>% pull}^2 %>% sum()
  q_value <- compute_q(group_ss)


# create a permutation list for every block
  permutation_list <- lapply(
    split(
      .data,
      .data %>% dplyr::pull({{ block }})
    ),
    function(.data) {
      do.call(rbind, combinat::permn(.data %>% dplyr::pull(lhs)))
    }
  )
  perm_names <- names(permutation_list)
# get number of blocks (n) and number of groups (p)
  np <- .data %>% dplyr::summarize(
    p = {{ group }} %>% unique() %>% length(),
    n = {{ block }} %>% unique() %>% length()
  )

  n_permutations <- factorial(np$p)^np$n
  group_vec <- .data %>% arrange({{block}}) %>% pull(group)
  get_qs <- function(n) {
    # get the permutation indices
    temp_permute <- nth_row(n = n, p = np$n, ninp = factorial(np$p))

    # get values according to permutation indices
    lhs <- vapply(1:length(perm_names), FUN = function(ind) {
      permutation_list[[ind]][temp_permute[ind], ]
    }, FUN.VALUE = numeric(np$p)) %>% as.vector()

    # compute q based on based on permuted ranks
    compute_q({aggregate(lhs, list(group_vec), FUN = sum) %>% pull(x)}^2 %>% sum)
  }
  qs <- get_q_vector(.f = get_qs, type = type, n_samples = n_samples, n_permutations)

  p_value <- mean(qs >= q_value)

  return(c(q = q_value, p = p_value))

}

get_q_vector <- function(.f, type, n_samples, n_permutations){
  if(type == "exact"){

    qs <- tryCatch(purrr::map_dbl(1:n_permutations, .f = .f, .progress = TRUE),
                   error = function(e){message(paste("the following error occured:",
                                                     e,
                                                     "Probably too many permutations for exact test.",
                                                     "Result is based on the montecarlo method with n_samples = ",
                                                     n_samples,
                                                     sep = "\n")
                   )
                     e})
    if("error" %in% class(qs) | length(qs)==0){
      if(length(qs)==0){
        message(paste("Probably too many permutations for exact test.",
                      "Result is based on the montecarlo method with n_samples = ",
                      n_samples,
                      sep = "\n"))
      }

      type <- "montecarlo"
    }

  }

  if(type == "montecarlo"){
    qs <- tryCatch(purrr::map_dbl(sample(1:n_permutations, size = 10000), .f = .f, .progress = TRUE),
                   error = function(e){warning(paste("the following error occured:",
                                                     e,
                                                     "n_samples is too large. Try smaller n_samples argument",
                                                     sep = "\n"))})
  }
  qs
}

#' create a compute q function
#'
#' @param .data dataset (may be permuted)
#' @param formula defines block, dv, and group variable
#'
#' @return q test statistic
#' @export
#'
#' @examples
create_q_f <- function(.data, formula){
  lhs <- rlang::f_lhs(formula)
  rhs <- rlang::f_rhs(formula)
  block <- rhs[[3]]
  group <- rhs[[2]]

  # compute required variables
  blocksum <- .data %>% group_by({{block}}) %>%
    summarize(blocksum = sum({{lhs}})) %>% pull

  block_ss <- sum(blocksum^2)
  total_blocksum <- sum(blocksum)
  total_blocksum_ss <- total_blocksum^2

  # number of groups
  p <- length(unique(.data %>% pull({{group}})))

  # create q function

  q_f <- function(group_ss){
    ((p-1)*(p*group_ss-total_blocksum_ss))/(p*total_blocksum - block_ss)
  }
}



