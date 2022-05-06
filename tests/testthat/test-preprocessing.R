test_that("data loads correctly", {

  example_files <- list.files(system.file('extdata', package = 'MetaboTandem'),
                              pattern = 'mzML', full.names = TRUE)
  metadata <- data.frame(FileName = example_files,
                         treatment = c(1, 2))

  data <- load_spectra_data(system.file('extdata', package = 'MetaboTandem'),
                    metadata = metadata,
                    format = 'mzML')

  expect_true(is(data, 'OnDiskMSnExp'))
})
