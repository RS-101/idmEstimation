# Test suite for setup_cpp_model and case-specific constructors

test_that("setup_cpp_model correctly splits data by case", {
  # Simulate mixed data with all 4 cases
  set.seed(999)
  n <- 40

  # Create status vectors
  status_dead <- c(rep(0, 10), rep(1, 10), rep(0, 10), rep(1, 10))
  status_ill <- c(rep(0, 10), rep(0, 10), rep(1, 10), rep(1, 10))

  # Create time vectors
  V_0 <- rep(0, n)
  V_healthy <- sort(runif(n, 0, 2))
  V_ill <- pmax(V_healthy, sort(runif(n, 1, 4)))
  T_obs <- pmax(V_ill, sort(runif(n, 3, 6)))

  # Setup model
  model_pointers_list <- setup_cpp_model(
    V_0, V_healthy, V_ill, T_obs,
    status_dead, status_ill
  )

    # Check that all cases are present
    expect_equal(length(model_pointers_list), 4)
    expect_equal(typeof(model_pointers_list$case_1), "externalptr")
    expect_equal(typeof(model_pointers_list$case_2), "externalptr")
    expect_equal(typeof(model_pointers_list$case_3), "externalptr")
    expect_equal(typeof(model_pointers_list$case_4), "externalptr")
})



