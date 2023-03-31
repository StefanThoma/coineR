test_that("tie correction", {
  df <- dplyr::tibble(vector = c(1.5, 1.5, 3))
  testthat::expect_equal(object = get_tie_correction(.data = df, vector)$t, expected = c(.75, .75))
})
