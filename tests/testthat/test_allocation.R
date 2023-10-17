
test_that("calc_zeta_discrete test cases work", {

  #Base case: expected r_h is above min
  expect_equal(calc_zeta_discrete(n_h = 20, phibar_h = .4),
               prod(calc_moments_truncbinom(size = 20, prob = .4)[c("E_X","E_one_over_X")]))

  #if phibar_h=.1, need at least 35 invitees for E(r_h)>=3.5
  #hence, make sure minimum works
  expect_equal(calc_zeta_discrete(n_h = 20, phibar_h = .1),
               prod(calc_moments_truncbinom(size = 35, prob = .1)[c("E_X","E_one_over_X")]))

  #Make sure rh_min can be turned off
  expect_equal(calc_zeta_discrete(n_h = 20, phibar_h = .1, rh_min = 0),
               prod(calc_moments_truncbinom(size = 20, prob = .1)[c("E_X","E_one_over_X")]))

  #Make sure fn works for vector of inputs
  expect_equal(calc_zeta_discrete(n_h = c(20, 30), phibar_h = c(.2, .5)),
               c(calc_zeta_discrete(n_h = 20, phibar_h = .2),
                 calc_zeta_discrete(n_h = 30, phibar_h = .5)))

  #Check rounding works
})


test_that("calc_zeta_discrete invalid args throw error", {
  expect_error(calc_zeta_discrete(n_h = 20.6, phibar_h = .4, round_flag = FALSE))
  expect_error(calc_zeta_discrete(n_h = c(10, 20), phibar_h = c(.1, .2, .3)))
  expect_error(calc_zeta_discrete(n_h = c(10, 20), phibar_h = c(.1, 0)))
  expect_error(calc_zeta_discrete(n_h = c(10, 20), phibar_h = c(1.01, .5)))
})
