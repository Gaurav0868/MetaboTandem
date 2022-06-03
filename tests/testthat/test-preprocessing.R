test_that("data loads correctly", {

  example_files <- list.files(system.file('extdata', package = 'MetaboTandem'),
                              pattern = 'mzML', full.names = TRUE)
  metadata <- data.frame(FileName = example_files,
                         treatment = 1)

  data <- load_spectra_data(system.file('extdata', package = 'MetaboTandem'),
                    metadata = metadata,
                    format = 'mzML')

  expect_true(is(data, 'OnDiskMSnExp'))
})

test_that('data is transformed to centroided correctly', {
  data("rp_data", package = 'MetaboTandem')
  data_cent <- centroid_check(rp_data, transform = TRUE)
  is.centroided <- unique(MSnbase::fData(data_cent)$centroided)
  expect_true(is.centroided)

})
