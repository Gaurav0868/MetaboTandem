#' Read mgf data
#'
#' This function is intended to read data in mgf format.
#' This function was modified from the `read_msp_mona()` function from
#' package `metID`
#'
#' @param file Database file in mgf format

read_mgf <- function(file,
                     separator = " ",
                     threads = parallelly::availableCores()) {
  print('Reading data in mgf format...')
  future::plan(strategy = future::multisession, workers = threads)
  db <- furrr::future_map(
    .x = file,
    .f = function(file) {
      data <- readr::read_lines(file, progress = FALSE)

      begin_idx <- which(data == "BEGIN IONS")
      end_idx <- which(data == "END IONS")

      data <- purrr::map2(
        .x = begin_idx,
        .y = end_idx,
        .f = function(idx1, idx2) {
          temp <- data[c(idx1:idx2)]
          temp <- temp[temp != "BEGIN IONS"]
          temp <- temp[temp != "END IONS"]
          if (length(temp) == 0) {
            return(NULL)
          }
          info <- temp[stringr::str_detect(temp, '=')] %>%
            stringr::str_split(., "=", n = 2) %>%
            do.call(rbind, .) %>%
            as.data.frame() %>%
            dplyr::rename(info = V1, value = V2)

          spec <- temp[!stringr::str_detect(temp, '=')]
          if(length(spec) != 0){
            spec <- spec %>%
              stringr::str_split(., separator, n = 2) %>%
              do.call(rbind, .) %>%
              as.data.frame() %>%
              dplyr::rename(mz = V1, intensity = V2) %>%
              dplyr::mutate(mz = as.numeric(mz),
                            intensity = as.numeric(intensity)) %>%
              as.data.frame()
          } else {
            spec <- NULL
          }

          list(info = info, spec = spec)
        }
      )
    }, .progress = TRUE
  )

  db <- Reduce(`c`, db)
  print('Finished reading db')

  return(db)
}

#' Read msp data
#'
#' This function is intended to read data in msp format.
#' This function was modified from the `read_msp_mona()` function from
#' package `metID`
#'
#' @param file Database file in msp format

read_msp <- function(file,
                     threads = parallelly::availableCores()) {
  print('Reading data in msp format..')
  future::plan(strategy = future::multisession, workers = threads)
  db <- furrr::future_map(
    .x = file,
    .f = function(file) {
      data <- readr::read_lines(file, progress = FALSE)

      if(data[1] != ""){
        data <- c("BEGIN IONS", data)
      } else {
        data[1] <- "BEGIN IONS"
      }

      if(tail(data, 1) != ""){
        data <- c(data, "END IONS")
      } else {
        data[length(data)] <- "END IONS"
      }

      data[data == ""] <- c("END IONS_split_BEGIN IONS")
      data <- unlist(stringr::str_split(data, '_split_'))

      begin_idx <- which(data == "BEGIN IONS")
      end_idx <- which(data == "END IONS")

      data <- purrr::map2(
        .x = begin_idx,
        .y = end_idx,
        .f = function(idx1, idx2) {
          temp <- data[c(idx1:idx2)]
          temp <- temp[temp != "BEGIN IONS"]
          temp <- temp[temp != "END IONS"]
          temp <- temp[temp != ""]
          if (length(temp) == 0) {
            return(NULL)
          }
          info <- temp[stringr::str_detect(temp, ':')] %>%
            stringr::str_split(., ": ", n = 2) %>%
            do.call(rbind, .) %>%
            as.data.frame() %>%
            dplyr::rename(info = V1, value = V2)

          spec <- temp[!stringr::str_detect(temp, ':')] %>%
            stringr::str_split(., " ", n = 2) %>%
            do.call(rbind, .) %>%
            as.data.frame() %>%
            dplyr::rename(mz = V1, intensity = V2) %>%
            dplyr::mutate(mz = as.numeric(mz),
                          intensity = as.numeric(intensity)) %>%
            as.data.frame()
          list(info = info, spec = spec)
        }
      )
    }, .progress = TRUE
  )

  db <- Reduce(`c`, db)
  print('Finished reading db')

  return(db)
}

parse_ms2_for_annotation <- function(sp, TITLE = NULL) {

  data <- c("BEGIN IONS",
            paste0("SCANS=", MSnbase::acquisitionNum(sp)),
            paste0("TITLE=", TITLE))

  if(MSnbase::polarity(sp) == 1){
    signal = "+"
  } else signal = "-"

  if (MSnbase::msLevel(sp) > 1) {
    data <- c(data,
              paste0("RTINSECONDS=", MSnbase::rtime(sp)),
              paste0("PEPMASS=", MSnbase::precursorMz(sp)))

    if (length(MSnbase::precursorCharge(sp)) &&
        !is.na(MSnbase::precursorCharge(sp))) {
      data <- c(data,
                paste0("CHARGE=", MSnbase::precursorCharge(sp), signal))
    }
  } else {
    data <- c(data,
              paste0("RTINSECONDS=", MSnbase::rtime(sp)))
  }

  spec <- unlist(purrr::map2(MSnbase::mz(sp),
                             MSnbase::intensity(sp),
                             ~paste(.x, .y, collapse = ' ')))

  data <- c(data, spec)
  data <- c(data, "END IONS")

  idx1 <- which(data == "BEGIN IONS")
  idx2 <- which(data == "END IONS")

  temp <- data[c(idx1:idx2)]
  temp <- temp[temp != "BEGIN IONS"]
  temp <- temp[temp != "END IONS"]
  if (length(temp) == 0) {
    return(NULL)
  }
  info <- temp[stringr::str_detect(temp, '=')] %>%
    stringr::str_split(., "=", n = 2) %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    dplyr::rename(info = V1, value = V2)

  spec <- temp[!stringr::str_detect(temp, '=')]
  if(length(spec) != 0){
    spec <- spec %>%
      stringr::str_split(., ' ', n = 2) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::rename(mz = V1, intensity = V2) %>%
      dplyr::mutate(mz = as.numeric(mz),
                    intensity = as.numeric(intensity)) %>%
      as.data.frame()
  } else {
    spec <- NULL
  }

  res <- list(info = info, spec = spec)

  return(res)
}
