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

  #n_total or cost_total need to be provided
  expect_error(opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
                                      phibar_h = c(.2, .1, .05)))

  #can't provide both n_total and cost_total
  expect_error(opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
                                      phibar_h = c(.2, .1, .05),
                                      n_total = 500,
                                      cost_total = 500))

  #if n_total is provided, can't also provide tau_h or c_NR_h
  expect_error(opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
                                      phibar_h = c(.2, .1, .05),
                                      n_total = 500,
                                      tau_h = 1))

  expect_error(opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
                                      phibar_h = c(.2, .1, .05),
                                      n_total = 500,
                                      c_NR_h = 1))

  #if cost_total is provided, need to provide tau_h and c_NR_h
  expect_error(opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
                                      phibar_h = c(.2, .1, .05),
                                      cost_total = 500,
                                      tau_h = 1))

  expect_error(opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
                                      phibar_h = c(.2, .1, .05),
                                      cost_total = 500,
                                      c_NR_h = 1))


})

test_that("opt_nh_nonresp_oneiter test cases work", {
  alloc_nmax_500 <- opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
                                           phibar_h = c(.2, .1, .05),
                                           n_total = 500)

  #alloc of n=500 should be equivalent to scenario with c_NR=2 and cost_total = 1000
  #  and with with const S_h and zeta thrown in
  expect_equal(alloc_nmax_500,
               opt_nh_nonresp_oneiter(N_h = c(1e4, 2e4, 5e4),
                                      phibar_h = c(.2, .1, .05),
                                      cost_total = 1000,
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
                                      cost_total = 1000),
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
                 cost_total = 1000)

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
  expect_equal(attr(opt_iter_nh3, "zeta_h_prev"), zeta_h2)
})


# opt_nh_nonresp_oneiter argument checking/validation ==========

## Check infer_objective() ======
test_that("opt_nh_nonresp_oneiter returns correct objective", {
  expect_equal(infer_objective(cost_total = 5), "min_var")
  expect_equal(infer_objective(n_total = 5), "min_var")
  expect_equal(infer_objective(Var_target = 5), "min_n")
  expect_equal(infer_objective(CV_target = 5), "min_n")
  expect_equal(infer_objective(Var_target = 5, c_NR_h = 1, tau_h = 1), "min_cost")
  expect_equal(infer_objective(CV_target = 5, c_NR_h = 1, tau_h = 1), "min_cost")
})

test_that("opt_nh_nonresp_oneiter throws error for partial cost data", {
  expect_error(infer_objective(Var_target = 5, c_NR_h = 1),
               regexp = "Did you forget to include `tau_h")
  expect_error(infer_objective(CV_target = 5, c_NR_h = 1),
               regexp = "Did you forget to include `tau_h")
  expect_error(infer_objective(Var_target = 5, tau_h = 1),
               regexp = "Did you forget to include `c_NR_h")
  expect_error(infer_objective(CV_target = 5, tau_h = 1),
               regexp = "Did you forget to include `c_NR_h")
})

## Check validate_args_min_var() ========

test_that("validate_args_min_var throws error/warning for invalid args", {
  expect_error(validate_args_min_var(was_objective_supplied = TRUE),
               regexp = "was specified without indicating the total sample size.*or costs")

  expect_warning(validate_args_min_var(was_objective_supplied = FALSE, n_total = 5, Ybar = 5),
                 regexp = "Ybar.*was provided but is not used")

  expect_error(validate_args_min_var(was_objective_supplied = FALSE,
                                     cost_total = 100, c_NR_h = c(2,3)),
               regexp = "tau_h\` is missing")

  expect_error(validate_args_min_var(was_objective_supplied = FALSE,
                                     cost_total = 100, tau_h = c(2,3)),
               regexp = "c_NR_h\` is missing")

  expect_error(validate_args_min_var(was_objective_supplied = FALSE,
                                     n_total = 100, tau_h = c(2,3)),
               regexp = "Fixing the total sample size.*means that cost information is extraneous")

  expect_error(validate_args_min_var(was_objective_supplied = FALSE,
                                     n_total = 100, c_NR_h = c(2,3)),
               regexp = "Fixing the total sample size.*means that cost information is extraneous")

  expect_error(validate_args_min_var(was_objective_supplied = FALSE,
                                     n_total = 100, tau_h = c(2,3), c_NR_h = c(2,3)),
               regexp = "Fixing the total sample size.*means that cost information is extraneous")
})

test_that("validate_args_min_var returns correct detailed objective", {
  expect_equal(validate_args_min_var(was_objective_supplied = FALSE,
                                     cost_total = 100, c_NR_h = c(2,3), tau_h = c(1.1, 1.2)),
               "Minimize variance subject to fixed total expected costs")

  expect_equal(validate_args_min_var(was_objective_supplied = FALSE,
                                     n_total = 20),
               "Minimize variance subject to fixed total number of invitees")
})

## Check validate_args_fixed_var() ========
test_that("validate_args_fixed_var throws error/warning for invalid args", {
  expect_error(validate_args_fixed_var(was_objective_supplied = FALSE, Var_target = 5),
               "S_h\` was not specified")

  expect_error(validate_args_fixed_var(was_objective_supplied = FALSE, CV_target = 5),
               "S_h\` was not specified")

  expect_warning(validate_args_fixed_var(was_objective_supplied = FALSE, S_h = c(1,2),
                                         Var_target = 5, Ybar = 5),
                 regexp = "Ybar\` was provided but will not be used")

  expect_error(validate_args_fixed_var(was_objective_supplied = FALSE, S_h = c(1,2),
                                         CV_target = 5),
                 regexp = "CV_target\` was specified without also specifying the population mean")


})

test_that("validate_args_fixed_var() throws error for mismatch between user-provided objective and other args", {
  expect_error(validate_args_fixed_var(S_h = c(1,2), Var_target = 5,
                                       was_objective_supplied = TRUE, objective = "min_cost"),
               regexp = "min_cost.*was specified but without providing cost structure information")

  expect_error(validate_args_fixed_var(S_h = c(1,2), Var_target = 5,
                                       c_NR_h = c(1,2),
                                       was_objective_supplied = TRUE, objective = "min_cost"),
               regexp = "min_cost.*was specified but without providing cost structure information")

  expect_error(validate_args_fixed_var(S_h = c(1,2), Var_target = 5,
                                       tau_h = c(1,2),
                                       was_objective_supplied = TRUE, objective = "min_cost"),
               regexp = "min_cost.*was specified but without providing cost structure information")

  expect_error(validate_args_fixed_var(S_h = c(1,2), Var_target = 5,
                                       was_objective_supplied = TRUE, objective = "min_n", c_NR_h = c(1,2)),
               regexp = "min_n.*was specified but some cost structure information.*was included")

  expect_error(validate_args_fixed_var(S_h = c(1,2), Var_target = 5,
                                       was_objective_supplied = TRUE, objective = "min_n", tau_h = c(1,2)),
               regexp = "min_n.*was specified but some cost structure information.*was included")

  expect_error(validate_args_fixed_var(S_h = c(1,2), Var_target = 5,
                                       was_objective_supplied = TRUE, objective = "min_n", c_NR_h = c(1,2), tau_h = c(1,2)),
               regexp = "min_n.*was specified but some cost structure information.*was included")
})

test_that("validate_args_fixed_var returns correct detailed objective", {
  expect_equal(validate_args_fixed_var(was_objective_supplied = FALSE,
                                       S_h = c(1,2), Var_target = 5,
                                       c_NR_h = c(1,2), tau_h = c(1,2)),
               "Minimize total expected costs subject to fixed variance")

  expect_equal(validate_args_fixed_var(was_objective_supplied = FALSE,
                                       S_h = c(1,2), Var_target = 5),
               "Minimize total n subject to fixed variance")

  expect_equal(validate_args_fixed_var(was_objective_supplied = FALSE,
                                       S_h = c(1,2), CV_target = 5, Ybar = 2,
                                       c_NR_h = c(1,2), tau_h = c(1,2)),
               "Minimize total expected costs subject to fixed CV")

  expect_equal(validate_args_fixed_var(was_objective_supplied = FALSE,
                                       S_h = c(1,2), CV_target = 5, Ybar = 2),
               "Minimize total n subject to fixed CV")
})


## Check validate_var_target()


test_that("validate_var_target() throws error for invalid args", {
  expect_error(validate_var_target(Var_target = 0),
               regexp = "must be.* positive and finite")

  expect_error(validate_var_target(Var_target = Inf),
               regexp = "must be.* positive and finite")

  expect_error(validate_var_target(CV_target = 0),
               regexp = "must be.* positive and finite")

  expect_error(validate_var_target(CV_target = Inf),
               regexp = "must be.* positive and finite")

  expect_error(validate_var_target(CV_target = 5, Ybar = 0),
               regexp = "Ybar\` must be nonzero and finite")

  expect_error(validate_var_target(CV_target = 1e250, Ybar = 1e100),
               regexp = "Choice of.*lead to \`Var_target`=Inf")
})

test_that("validate_var_target returns correct target variance", {
  expect_equal(validate_var_target(Var_target=5),
               5)
  expect_equal(validate_var_target(CV_target=10, Ybar=2), 400)
})
