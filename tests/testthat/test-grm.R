test_that("10 iter, GP, exponential covariance", {
              ctm_fit <- grm(Y = cmaq_aqs_matched$pm25,
                             X = cmaq_aqs_matched$ctm,
                             L = cmaq_aqs_matched[, c("elevation", "forestcover",
                                                      "hwy_length", "lim_hwy_length", 
                                                      "local_rd_length", "point_emi_any")],
                             M = cmaq_aqs_matched[, c("tmp", "wind")],
                             n.iter = 10,
                             burn = 2,
                             thin = 4,
                             covariance = "exponential",
                             coords = cmaq_aqs_matched[, c("x", "y")],
                             space.id = cmaq_aqs_matched$space_id,
                             time.id = cmaq_aqs_matched$time_id,
                             spacetime.id = cmaq_aqs_matched$spacetime_id,
                             verbose.iter = 10)
              expect_equal(sum(is.na(ctm_fit$alpha.time)), 0)
              expect_equal(sum(is.na(ctm_fit$beta.time)), 0)
              expect_equal(sum(is.na(ctm_fit$alpha.space)), 0)
              expect_equal(sum(is.na(ctm_fit$beta.space)), 0)
              expect_equal(sum(is.na(ctm_fit$others)), 0)
})
