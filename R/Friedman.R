#' compute tie correction
#'
#' @param .data dataframe, may be grouped
#' @param vector name of vector in dataframe on which you want to compute the tie vector
#'
#' @return dataframe including t as tievector
#' @export
#' @importFrom magrittr "%>%"
#' @examples
#' get_tie_correction(.data = tibble(vec = c(1:6, 1:6)), vec)
get_tie_correction <- function(.data, vector){

  if (dplyr::is_grouped_df(.data)) {
    return(dplyr::do(.data, get_tie_correction(., {{vector}})))
  }

  .data %>%
    dplyr::select({{vector}}) %>%
    table() %>%
    tibble::tibble(n = .) %>%
    dplyr::mutate(
      t = 1-(sum(n^3-n)/(sum(n)^3-sum(n)))
    )
}



#' Compute q based on data and formula
#'
#' not used right now
#'
#' @param .data
#' @param formula
#' @param total_t
#'
#' @return
#' @export
#' @importFrom magrittr "%>%"
#'
#' @examples
get_q <- function(.data, formula, total_t = NULL){
  # check inputs

  # parse formula
  lhs <- rlang::f_lhs(formula)
  rhs <- rlang::f_rhs(formula)
  block <- rhs[[3]]
  group <- rhs[[2]]

  # compute rank & ranksum
  .data <- .data %>% dplyr::mutate(
    .by = {{block}},
    rank_by_block = coin::rank_trafo({{lhs}})
  ) %>% dplyr::mutate(
    .by = {{group}},
    ranksum_group = sum(rank_by_block),
  )

  # compute tie correction
  if(is.null(total_t)){
    total_t <- .data %>% dplyr::group_by({{block}}) %>%
      get_tie_correction(rank_by_block) %>%
      dplyr::summarize(t = unique(t)) %>%
      dplyr::pull(t) %>% sum()
  }
  # compute q
  # first, get number of groups and blocks
  np <- .data %>% dplyr::summarize(
    p = {{group}} %>% unique() %>% length(),
    n = {{block}} %>% unique() %>% length()
  )

  # extract ranksum vector
  rs <- .data %>% dplyr::group_by({{group}}) %>% dplyr::summarize(
    rs = max(ranksum_group)
  ) %>% dplyr::pull()

  1/total_t * ((12/(np$p*(np$p+1)) * sum(rs^2)) - 3 * np$n^2 * (np$p+1))
}


# try it out




#' compute friedman statistic and p value
#'
#'
#' @param .data dataframe including all relevant vectors
#' @param formula of the form: dv ~ group | block
#' @param type exact
#'
#' @return
#' @export
#'
#' @importFrom magrittr "%>%"
#' @examples
#' participant  <- 1:4
#' group <- 1:3
#' data <- expand.grid(
#'   "participant" = participant,
#'   "group" = group
#' )
#' data$ability <- c(7, 8, 8, 6,
#'                   10, 12, 8, 10,
#'                   11, 12, 13, 11)
#'
#'
#' friedman_test(.data = data, formula = ability ~ group | participant)
friedman_test <- function(.data, formula, type = "exact"){

  #permutations follow a pattern within block:
  # parse formula
  lhs <- rlang::f_lhs(formula)
  rhs <- rlang::f_rhs(formula)
  block <- rhs[[3]]
  group <- rhs[[2]]

  # get ranks
  .data <- .data %>% dplyr::mutate(
    .by = {{block}},
    rank = coin::rank_trafo({{lhs}})
  )
  total_t <- .data %>% dplyr::group_by({{block}}) %>%
    dplyr::mutate(
      rank_by_block = coin::rank_trafo({{lhs}})
    ) %>%
    get_tie_correction(rank_by_block) %>%
    dplyr::summarize(t = unique(t)) %>%
    dplyr::pull(t) %>% sum()

  permutation_list <- lapply(split(.data,
                                   .data %>% dplyr::pull({{block}})),
                             function(.data){do.call(rbind, combinat::permn(.data %>% pull(rank)))})

  np <- .data %>% dplyr::summarize(
    p = {{group}} %>% unique() %>% length(),
    n = {{block}} %>% unique() %>% length()
  )

  n_permutations <- factorial(np$p)^np$n

  # sorting here is essential for the get_gs fuction
  .data <- .data %>% dplyr::arrange({{block}})
  perm_names <- names(permutation_list)
  get_qs <- function(n){
    # get the permutation indices
    temp_permute <- nth_row(n = n, p = np$n, ninp = factorial(np$p))

    # get ranks according to permutation indices
    .data$rank <- vapply(1:length(perm_names), FUN = function(ind){
      permutation_list[[ind]][temp_permute[ind],]
    }, FUN.VALUE = numeric(np$p)) %>% as.vector()

    # compute q based on based on permuted ranks
    get_q_light(.data, rank ~ group | block, total_t = total_t, np = np)
  }
  qs <- sapply(1:n_permutations, FUN = get_qs)

  p <- sum(qs >= qs[1]) / n_permutations

  list(p = p, q = qs[1])
  # perms_in_block <- rbind(1:np$p, permute::allPerms(1:np$p))
  # perm_options <- nrow(perms_in_block)
  #
  # perm_matrix <- matrix(rep(1:perm_options, times = np$n), nrow = perm_options) %>% as_tibble()
  # colnames(perm_matrix) <- 1:perm_options
  #
  # # from each row of perm_grid
  # perm_grid <- expand.grid(perm_matrix)[4:1]
  #
  # total_permutations <- nrow(perm_grid)
  # # take the ranking out of the core function, otherwise it is repeated very often.
}





#' nth row
#' get nth row of permutation table of certain specification.
#' namely, the p permutation vectors are all the same and are consecutive integers
#' starting form 1 to p
#' @param n
#' @param p
#' @param ninp
#'
#' @return
#' @export
#' @importFrom magrittr "%>%"
#' @examples
#' lapply(1:1000, function(n) nth_row(n, 30, 5))
#' # or simply:
#' nth_row(10, 3, 10)
nth_row <- function(n, p = 3, ninp = 5) {
  basevec <- rep(1, p)

  x <- Rmpfr::mpfr(n-1, precBits = 64, base = 10)  ## base = 10 is default
  numb <- as.numeric(format(x, base = ninp))
  vec <- as.numeric(
    strsplit(as.character(numb),"") %>% dplyr::as_vector())

  vec <- c(rep(0, length(basevec) - length(vec) ), vec)
  basevec + vec
}





#' light version of get_q
#' does less computation
#' stuff that can be computed once is computed outside the function
#'
#' @param .data
#' @param formula
#' @param total_t
#' @importFrom magrittr "%>%"
#' @param np
#'
#' @return
#' @export
#'
#' @examples
get_q_light <- function(.data, formula, total_t, np){
  # check inputs

  # parse formula
  lhs <- rlang::f_lhs(formula)
  rhs <- rlang::f_rhs(formula)
  block <- rhs[[3]]
  group <- rhs[[2]]

  # compute ranksum vector
  rs <- .data %>% dplyr::mutate(
    .by = {{group}},
    ranksum_group = sum({{lhs}}),
  ) %>% dplyr::group_by({{group}}) %>% dplyr::summarize(
    rs = max(ranksum_group)
  ) %>% dplyr::pull()

  1/total_t * ((12/(np$p*(np$p+1)) * sum(rs^2)) - 3 * np$n^2 * (np$p+1))
}