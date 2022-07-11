test_that("data loads correctly", {

  example_files <- list.files(system.file('extdata', package = 'MetaboTandem'),
                              pattern = 'test.mzML', full.names = TRUE)
  metadata <- data.frame(FileName = example_files,
                         SampleID = c('RP_1_test', 'RP_2_test'),
                         treatment = 1)

  rp_data <- load_spectra_data(system.file('extdata', package = 'MetaboTandem'),
                    metadata = metadata,
                    format = 'mzML')

  expect_true(is(rp_data, 'OnDiskMSnExp'))
})

test_that('data is transformed to centroided correctly', {
  data("rp_data", package = 'MetaboTandem')
  data_cent <- centroid_check(rp_data, transform = TRUE)
  is.centroided <- unique(MSnbase::fData(data_cent)$centroided)
  expect_true(is.centroided)

})

test_that('peaks can be picked', {
  data("rp_data", package = 'MetaboTandem')
  data_cent <- centroid_check(rp_data, transform = TRUE)
  data_cent <- apply_peak_picking(data_cent,
                                  p.width = c(20, 50),
                                  snt = 0,
                                  noise = 10000)
  expect_true(xcms::hasChromPeaks(data_cent))
})

test_that('alignment and correspondence work', {
  data("rp_data", package = 'MetaboTandem')
  data_cent <- centroid_check(rp_data, transform = TRUE)
  data_cent <- apply_peak_picking(data_cent,
                                  p.width = c(20, 50),
                                  snt = 0,
                                  noise = 1000)
  data_align <- apply_alignment(data = data_cent,
                                metadata = MSnbase::pData(data_cent),
                                min_frac = 0.5,
                                min_samples = 1,
                                group_by = 'treatment',
                                plot = FALSE)

  data_grouped <- apply_correspondence(data_align,
                                       metadata = MSnbase::pData(data_cent),
                                       min_frac = 0.5,
                                       min_samples = 1,
                                       group_by = 'treatment')

  expect_true(xcms::hasAdjustedRtime(data_align))
  expect_true(xcms::hasFeatures(data_grouped))

})
