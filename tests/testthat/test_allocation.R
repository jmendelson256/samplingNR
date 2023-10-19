# calc_zeta_discrete tests ==========

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
  expect_equal(calc_zeta_discrete(n_h = c(20.4999, 20.50001), c(phibar_h = c(.4, .6))),
               calc_zeta_discrete(n_h = c(20, 21), c(phibar_h = c(.4, .6))))
})

test_that("calc_zeta_discrete invalid args throw error", {
  expect_error(calc_zeta_discrete(n_h = 20.6, phibar_h = .4, round_flag = FALSE))
  expect_error(calc_zeta_discrete(n_h = c(10, 20), phibar_h = c(.1, .2, .3)))
  expect_error(calc_zeta_discrete(n_h = c(10, 20), phibar_h = c(.1, 0)))
  expect_error(calc_zeta_discrete(n_h = c(10, 20), phibar_h = c(1.01, .5)))
})

# calc_zeta tests ==========

test_that("calc_zeta test cases work", {
  main_res <- calc_zeta(n_h = c(10.2, 10.8), phibar_h = c(.5, .6))
  tmp_res1 <- calc_zeta(n_h = c(10, 10), phibar_h = c(.5, .6))
  tmp_res2 <- calc_zeta(n_h = c(11, 11), phibar_h = c(.5, .6))
  expect_equal(main_res,
               tmp_res1 * c(.8, .2) + tmp_res2 * c(.2, .8))
})

test_that("calc_zeta verbose flag works", {
  res_normal <- calc_zeta(n_h = c(10.2, 15.7), phibar_h = c(.4, .5), verbose_flag = FALSE)
  res_verbose <- calc_zeta(n_h = c(10.2, 15.7), phibar_h = c(.4, .5), verbose_flag = TRUE)
  expect_equal(res_normal, res_verbose$res)
})

# opt_nh_nonresp_oneiter tests ========


test_that("opt_nh_nonresp_oneiter correctly validates arguments", {

  #n_max or c_max need to be provided
  expect_error(opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
                                      phibar_h = c(.2, .1, .05)))

  #can't provide both n_max and c_max
  expect_error(opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
                                      phibar_h = c(.2, .1, .05),
                                      n_max = 500,
                                      c_max = 500))

  #if n_max is provided, can't also provide tau_h or c_NR_h
  expect_error(opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
                                      phibar_h = c(.2, .1, .05),
                                      n_max = 500,
                                      tau_h = 1))

  expect_error(opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
                                      phibar_h = c(.2, .1, .05),
                                      n_max = 500,
                                      c_NR_h = 1))

  #if c_max is provided, need to provide tau_h and c_NR_h
  expect_error(opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
                                      phibar_h = c(.2, .1, .05),
                                      c_max = 500,
                                      tau_h = 1))

  expect_error(opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
                                      phibar_h = c(.2, .1, .05),
                                      c_max = 500,
                                      c_NR_h = 1))


})

test_that("opt_nh_nonresp_oneiter test cases work", {
  alloc_nmax_500 <- opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
                                           phibar_h = c(.2, .1, .05),
                                           n_max = 500)

  #alloc of n=500 should be equivalent to scenario with c_NR=2 and c_max = 1000
  #  and with with const S_h and zeta thrown in
  expect_equal(alloc_nmax_500,
               opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
                                      phibar_h = c(.2, .1, .05),
                                      c_max = 1000,
                                      c_NR_h = 2,
                                      tau_h = 1,
                                      S_h = rep(2, 3),
                                      zeta_h = rep(1.01, 3)))

  #Check that roughly matches manual calculations
  expect_equal(alloc_nmax_500, c(36.1574021,
                                 102.2685769,
                                 361.574021),
               tolerance = 1e-6)

  #Check more complicated manual example
  expect_equal(opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
                                      S_h = c(1.1, 1.5, 1.8),
                                      phibar_h = c(.2, .1, .05),
                                      zeta_h = c(1.11, 1.05, 1.3),
                                      tau_h = c(1.5, 1.9, 3.1),
                                      c_NR_h = c(2.2, 2.1, 2),
                                      c_max = 1000),
               c(19.1561051, 73.8871177, 354.9839879),
               tolerance = 1e-6)
})

# opt_nh_nonresp tests =========

test_that("opt_nh_nonresp test case works", {
  #Compute 3 iterations (including base case) manually with opt_nh_nonresp_oneiter
  #See if results can be matched with opt_nh_nonresp
  myargs <- list(N_h = c(1e4, 2e4, 5e4),
                 S_h = c(1.1, 1.5, 1.8),
                 phibar_h = c(.2, .1, .05),
                 tau_h = c(1.5, 1.9, 3.1),
                 c_NR_h = c(2.2, 2.1, 2),
                 c_max = 1000)

  opt_nh1 <- do.call(opt_nh_nonresp_oneiter, myargs)
  zeta_h1 <- calc_zeta(n_h = opt_nh1, phibar_h = myargs$phibar_h)
  opt_nh2 <- do.call(opt_nh_nonresp_oneiter,
                     append(myargs, list(zeta_h = zeta_h1)))
  zeta_h2 <- calc_zeta(n_h = opt_nh2, phibar_h = myargs$phibar_h)
  opt_nh3 <- do.call(opt_nh_nonresp_oneiter,
                     append(myargs, list(zeta_h = zeta_h2)))

  opt_iter_nh3 <- do.call(opt_nh_nonresp,
                          append(myargs, list(max_iter = 3)))

  expect_equal(opt_nh3,
               c(opt_iter_nh3))
  expect_equal(attr(opt_iter_nh3, "num_iter"), 3)
  expect_equal(attr(opt_iter_nh3, "zeta_h_min_1"), zeta_h2)
})
