test_that("compute cochran test", {
  autos <- matrix(data = c(0,0,0,0,0,1,0,0,0,1,
                           0,1,0,0,1,0,0,1,0,1,
                           0,1,1,0,1,1,0,1,1,1,
                           0,1,1,0,0,1,1,1,1,1),
                  nrow = 10, ncol = 4, byrow = FALSE) %>% as_tibble() %>%
    mutate(family = row_number()) %>%
    pivot_longer(V1:V4, names_to = "group")

  expect_equal(object = cochran_test(autos, formula = value ~ group | family, type = "exact"), c(q = 9, p = "notsure"))
})
