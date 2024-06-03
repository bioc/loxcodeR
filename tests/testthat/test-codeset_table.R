test_that("Data can be loaded", {
    expect_error(
        readRDS("~/Desktop/temp22/temp221/loxcodeR/LoxcodeR_app/data-2024-05-20.rds"),
        regexp = NA
    )
})

