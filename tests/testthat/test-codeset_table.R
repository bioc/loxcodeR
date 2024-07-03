test_that("Data can be loaded", {
    expect_error(
        readRDS(system.file("extdata","data-2024-05-20.rds",package="loxcodeR")),
        regexp = NA
    )
})

