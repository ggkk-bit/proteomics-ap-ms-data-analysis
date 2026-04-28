context("test-SMAD")

test_that("Test Input Data", {
    expect_is(TestDatInput, "data.frame")
    expect_equal(nrow(TestDatInput), 5000)
})