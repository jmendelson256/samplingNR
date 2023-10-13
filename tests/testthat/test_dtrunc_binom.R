
test_that("dtruncbinom test case works", {
  expect_equal(dtruncbinom(1:5, 5, .2),
              dbinom(1:5, size = 5, prob = .2) /
                sum(dbinom(1:5, size = 5, prob = .2)))
})


test_that("calc_moments_truncbinom moments match theoretical moments", {
  #Fn to compute theoretical moments of truncated binomial distribution
  calc_moments_truncbinom_theor <- function(size, prob) {
    theor_E_X <- size * prob / (1 - (1 - prob)^size)
    theor_E_X_sq <- size * prob * (1 - prob + size * prob) / (1 - (1 - prob)^size)
    theor_V_X <- theor_E_X_sq - theor_E_X^2

    theoretical_moments <- c(theor_E_X = theor_E_X,
                             theor_E_X_sq = theor_E_X_sq,
                             theor_V_X = theor_V_X)

    (theoretical_moments)
  }

  expect_equal(unname(calc_moments_truncbinom(30, .1)[c("E_X", "E_X_sq", "V_X")]),
               unname(calc_moments_truncbinom_theor(30, .1)),
               tolerance = testthat_tolerance())
})

