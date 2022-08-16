#' Adduct annotation with CAMERA
#'
#' Function to apply peak picking parameters on a data subset
#'
#'
#' @param data An [MSnExp-class] object with peaks picked
#' @param group_by column
#'
#' @export
#'
#'
annotate_CAMERA <- function(data, group_by){
  ms1_data <- xcms::filterMsLevel(data, msLevel = 1L)
  xs <- as(ms1_data, 'xcmsSet')

  # Setting up sampclass
  sampclass(xs) <- MSnbase::pData(data)[[group_by]]
  sampnames(xs) <- MSnbase::pData(data)$SampleID

  an <- CAMERA::xsAnnotate(xs)
  anF <- CAMERA::groupFWHM(an, perfwhm = 0.6)
  anI <- CAMERA::findIsotopes(anF, mzabs = 0.01)
  anIC <- CAMERA::groupCorr(anI, cor_eic_th = 0.75)
  anFA <- CAMERA::findAdducts(anIC, polarity="positive")


  peaklist <- CAMERA::getPeaklist(anFA)

  return(peaklist)
}

#' Modified - Identify metabolites using multiple databases
#'
#' Function to identify metabolites using multiple databases at once. This
#' function was modified from the package [metid](https://github.com/tidymass/metid)
#'
#'@param ms1.data A [data.frame] with 3 columns: 'name', 'mz' and 'rt'
#'@param ms2.data A character containign the MS2 data in mgf format
#'@param parameter.list A list containing the parameters to process each database.
#'Parameters must be created using [metid::metIdentifyParam]
#'
#'@return A [mzIdentifyClass] object wth the annotation results
#'

mod_identify_all <- function (ms1.data,
                              ms2.data,
                              db_dir,
                              parameter.list){


  threads = parameter.list[[1]]$threads

  ms1.data <- ms1.data %>%
    dplyr::select(name = FeatureID,
                  mz = mzmed,
                  rt = rtmed)

  message('Parsing MS2 data')

  ms2.data <- purrr::map2(ms2.data@listData,
                          ms2.data@elementMetadata@listData[["feature_id"]],
                          ~ parse_ms2_for_annotation(.x, TITLE = .y))

  message('Preparing MS2 data for identification0')

  ms2.data <- purrr::imap(ms2.data, function(x, y){
    info <- x$info
    info <- tibble(name = y,
                   mz = info[info[1] == 'PEPMASS', 2],
                   rt = info[info[1] == 'RTINSECONDS', 2]) %>%
      dplyr::mutate(mz = as.numeric(mz),
             rt = as.numeric(rt))
    x$info <- info
    x
  })

  ms1.info <- purrr::map(ms2.data, ~.x[[1]])
  ms1.info <- do.call(rbind, ms1.info)

  ms2.info <- purrr::map(ms2.data, ~.x[[2]])

  database_names <- purrr::map(parameter.list, function(x){
    x[[1]]$database
  })

  # Load selected databases
  dbs <- purrr::map(parameter.list, function(x){
    path <- file.path(db_dir, paste0(x[[1]]$database, '.rds'))
    if(file.exists(path)){
      print(paste0('Loading ', x[[1]]$database))
      x <- readRDS(path)

    } else {
      print(paste0(x[[1]]$database, ' is not in the database_directory. Skipping'))
      x <- NULL
    }

    return(x)
  })

  # Removing parameters of databases not found in the directory

  # Create a data frame with info of the databases to be used

  identification_results <- purrr::map2(parameter.list, dbs, function(param, db){
    if(!is.null(db)){
      message(paste0('Starting to annotate with ', param[[1]]$database))

      result <-mod_metIdentify(ms1.data = ms1.data,
                               ms1.info = ms1.info,
                               ms2.data = ms2.data,
                               ms2.info = ms2.info,
                               ms1.ms2.match.mz.tol = param[[1]]$ms1.ms2.match.mz.tol,
                               ms1.ms2.match.rt.tol = param[[1]]$ms1.ms2.match.rt.tol,
                               ms1.match.ppm = param[[1]]$ms1.match.ppm,
                               ms2.match.ppm = param[[1]]$ms2.match.ppm,
                               mz.ppm.thr = param[[1]]$mz.ppm.thr,
                               ms2.match.tol = param[[1]]$ms2.match.tol,
                               fraction.weight = param[[1]]$fraction.weight,
                               dp.forward.weight = param[[1]]$dp.forward.weight,
                               dp.reverse.weight = param[[1]]$dp.reverse.weight,
                               rt.match.tol = param[[1]]$rt.match.tol,
                               polarity = param[[1]]$polarity,
                               ce = param[[1]]$ce,
                               column = param[[1]]$column,
                               ms1.match.weight = param[[1]]$ms1.match.weight,
                               rt.match.weight = param[[1]]$rt.match.weight,
                               ms2.match.weight = param[[1]]$ms2.match.weight,
                               total.score.tol =param[[1]]$total.score.tol,
                               candidate.num = param[[1]]$candidate.num,
                               database = db,
                               database_name = param[[1]]$database,
                               threads = param[[1]]$threads)


    }

  })

  print('Annotation finished')

  return(identification_results)


}

#' Modified - Identify metabolites using
#'
#' Function to identify metabolites using multiple databases at once. This
#' function was modified from the package [metid](https://github.com/tidymass/metid)
#'
#'@param ms1.data A [data.frame] with 3 columns: 'name', 'mz' and 'rt'
#'@param ms2.data A character containign the MS2 data in mgf format
#'@param parameter.list A list containing the parameters to process each database.
#'Parameters must be created using [metid::metIdentifyParam]
#'
#'@return A [mzIdentifyClass] object with the annotation results
#'

mod_metIdentify <- function(ms1.data,
                            ms1.info,
                            ms2.data = NULL,
                            ms2.info,
                            ms1.ms2.match.mz.tol = 25,
                            ms1.ms2.match.rt.tol = 10,
                            ms1.match.ppm = 25,
                            ms2.match.ppm = 30,
                            mz.ppm.thr = 400,
                            ms2.match.tol = 0.5,
                            fraction.weight = 0.3,
                            dp.forward.weight = 0.6,
                            dp.reverse.weight = 0.1,
                            rt.match.tol = 30,
                            polarity = c("positive", "negative"),
                            ce = "all",
                            column = c("rp", "hilic"),
                            ms1.match.weight = 0.25,
                            rt.match.weight = 0.25,
                            ms2.match.weight = 0.5,
                            total.score.tol = 0.5,
                            candidate.num = 3,
                            database_name,
                            database,
                            threads = 3){

  polarity <- match.arg(polarity)
  column <- match.arg(column)

  ce.list.pos <- unique(unlist(purrr::map(database@spectra.data$Spectra.positive,
                                          names)))
  ce.list.neg <- unique(unlist(purrr::map(database@spectra.data$Spectra.negative,
                                          names)))

  ce.list <- ifelse(polarity == "positive", ce.list.pos, ce.list.neg)
  if (all(ce %in% ce.list) & ce != "all") {
    stop("All ce values you set are not in database. Please check it.\n")
    ce <- ce[ce %in% ce.list]
  }
  rm(list = c("ce.list.pos", "ce.list.neg", "ce.list"))
  if (all(ce != "all")) {
    if (polarity == "positive") {
      ce.list <- unique(unlist(lapply(database@spectra.data$Spectra.positive,
                                      function(x) {
                                        names(x)
                                      })))
      if (length(grep("Unknown", ce.list)) > 0) {
        ce <- unique(c(ce, grep(pattern = "Unknown",
                                ce.list, value = TRUE)))
      }
    }
    else {
      ce.list <- unique(unlist(lapply(database@spectra.data$Spectra.negative,
                                      function(x) {
                                        names(x)
                                      })))
      if (length(grep("Unknown", ce.list)) > 0) {
        ce <- unique(c(ce, grep(pattern = "Unknown",
                                ce.list, value = TRUE)))
      }
    }
  }
  if (!database@database.info$RT) {
    message(crayon::yellow("No RT information in database.\nThe weight of RT have been set as 0."))
  }

  if (polarity == "positive" & column == "hilic") {
    data("hilic.pos", package = 'metid', envir = environment())
    adduct.table <- hilic.pos
  }
  if (polarity == "positive" & column == "rp") {
    data("rp.pos", package = 'metid', envir = environment())
    adduct.table <- rp.pos
  }
  if (polarity == "negative" & column == "hilic") {
    data("hilic.neg", package = 'metid', envir = environment())
    adduct.table <- hilic.neg
  }
  if (polarity == "negative" & column == "rp") {
    data("rp.neg", package = 'metid', envir = environment())
    adduct.table <- rp.neg
  }

  message('Matching peak table with MS2 data')
  match.result <- masstools::mz_rt_match(data1 = ms1.data[, c(2, 3)],
                                         data2 = as.data.frame(ms1.info[, c(2, 3)]),
                                         mz.tol = ms1.ms2.match.mz.tol,
                                         rt.tol = ms1.ms2.match.rt.tol,
                                         rt.error.type = "abs")

  if (is.null(match.result)) return("No peaks are matched with MS2 spectra.\n")
  if (nrow(match.result) == 0) return("No peaks are matched with MS2 spectra.\n")
  message(paste(c(length(unique(match.result[, 1])),
                  "out of", nrow(ms1.data),
                  "peaks have MS2 spectra."),
                collapse = ' '))
  message("Selecting the most intense MS2 spectrum for each peak")
  temp.idx <- unique(match.result[, 1])
  match.result <- purrr::map(temp.idx, function(idx) {
    idx2 <- match.result[which(match.result[, 1] == idx), 2]
    if (length(idx2) == 1) {
      return(c(idx, idx2))
    }
    else {
      temp.ms2.info <- ms2.info[idx2]
      return(c(idx,
               idx2[which.max(unlist(purrr::map(temp.ms2.info, function(y){
                 y <- y[order(y[, 2], decreasing = TRUE),
                        , drop = FALSE]
                 if (nrow(y) > 5) y <- y[seq_len(5), ]
                 sum(y[, 2])
               })))]))
    }
  })

  match.result <- as.data.frame(do.call(rbind, match.result)) %>%
    dplyr::rename("Index1" = V1, "Index2" = V2) %>%
    dplyr:: mutate(MS1.peak.name = ms1.data$name[Index1],
                   MS2.spectra.name = ms1.info$name[Index2])

  ms1.info <- ms1.info[unique(match.result[, 2]), , drop = FALSE]
  ms2.info <- ms2.info[unique(match.result[, 2])]

  match.result$Index.ms2.spectra <- match(match.result$MS2.spectra.name,
                                          ms1.info$name)

  ms2Matchresult <- metIdentification(ms1.info = ms1.info,
                                      ms2.info = ms2.info,
                                      polarity = polarity,
                                      ce = ce,
                                      database = database,
                                      ms1.match.ppm = ms1.match.ppm,
                                      ms2.match.ppm = ms2.match.ppm,
                                      mz.ppm.thr = mz.ppm.thr,
                                      ms2.match.tol = ms2.match.tol,
                                      rt.match.tol = rt.match.tol,
                                      column = column,
                                      ms1.match.weight = ms1.match.weight,
                                      rt.match.weight = rt.match.weight,
                                      ms2.match.weight = ms2.match.weight,
                                      total.score.tol = total.score.tol,
                                      candidate.num = candidate.num,
                                      adduct.table = adduct.table,
                                      threads = threads,
                                      fraction.weight = fraction.weight,
                                      dp.forward.weight = dp.forward.weight,
                                      dp.reverse.weight = dp.reverse.weight)
  return.result <- new(Class = "metIdentifyClass",
                       ms1.data = ms1.data,
                       ms1.info = ms1.info,
                       ms2.info = ms2.info,
                       identification.result = ms2Matchresult,
                       match.result = match.result,
                       adduct.table = adduct.table,
                       ms1.ms2.match.mz.tol = ms1.ms2.match.mz.tol,
                       ms1.ms2.match.rt.tol = ms1.ms2.match.rt.tol,
                       ms1.match.ppm = ms1.match.ppm,
                       ms2.match.ppm = ms2.match.ppm,
                       ms2.match.tol = ms2.match.tol,
                       rt.match.tol = rt.match.tol,
                       polarity = polarity,
                       ce = paste(ce, collapse = ";"),
                       column = column,
                       ms1.match.weight = ms1.match.weight,
                       rt.match.weight = rt.match.weight,
                       ms2.match.weight = ms2.match.weight,
                       total.score.tol = total.score.tol,
                       candidate.num = candidate.num,
                       database = database_name,
                       threads = threads,
                       version = "1.0.0")

}

#' Modified - Get identification table
#'
#' Function to retrieve annotation table. This
#' function was modified from the package [metid](https://github.com/tidymass/metid)
#'
#' @param ... One or multiple metIdentifyClass objects.
#' @param candidate.num The number of candidates.
#' @param type The type of identification table.
#' @return A identification table (data.frame).
#' @export

mod_get_identification_table = function(object,
                                    candidate.num = 3) {
  candidate.num <- round(candidate.num)
  if (candidate.num <= 0) {
    candidate.num <- 1
  }

  if (class(object) != "metIdentifyClass") {
    stop("Only for metIdentifyClass\n")
  }

  database <- object@database

  identification.result <- object@identification.result

  if (is.null(identification.result[[1]])) {
    return(NULL)
  }
  #
  # if (nrow(object@match.result) == 0) {
  #   message(crayon::yellow("The object is identified without MS2 spectra."))
  #   return(
  #     getIdentificationTable2(
  #       object = object,
  #       candidate.num = candidate.num,
  #       type = type,
  #       silence.deprecated = TRUE
  #     )
  #   )
  # }

  ##add database information
  identification.result <-
    purrr::map(identification.result, function(x) {
      if (nrow(x) > candidate.num) {
        x <- x[seq_len(candidate.num), , drop = FALSE]
      }
      data.frame(x,
                 "Database" = object@database)
    })

  peak.table <- object@ms1.data
  match.result <- object@match.result

  identification.table <-
    vector(mode = "list", length = nrow(peak.table))
  names(identification.table)[match.result[, 1]] <-
    object@ms1.info$name
  identification.table[match(names(identification.result),
                             names(identification.table))] <-
    identification.result

  peak.table <- purrr::pmap(peak.table, function(name, mz, rt){
    data.frame(name = name, mz = mz, rt = rt)
  })

  identification.table <- purrr::map2(
    peak.table, identification.table, function(x,y){
      if (all(is.na(y))) {
        temp <- data.frame(
          name = x$name,
          mz = x$mz,
          rt = x$rt,
          Compound.name = NA,
          CAS.ID = NA,
          HMDB.ID = NA,
          KEGG.ID = NA,
          Lab.ID = NA,
          Adduct = NA,
          mz.error = NA,
          mz.match.score = NA,
          RT.error = NA,
          RT.match.score = NA,
          CE = NA,
          SS = NA,
          Total.score = NA,
          Database = NA
        )
      } else {

        temp <- data.frame(
          name = rep(x$name, nrow(y)),
          mz = rep(x$mz, nrow(y)),
          rt = rep(x$rt, nrow(y))
        )

        temp <- cbind(temp, y)
      }

      return(temp)
    }
  )

  identification.table <-do.call(rbind, identification.table)

  ###add Candidate.number
  Candidate.number <- purrr::map2(
    .x = object@identification.result,
    .y = names(object@identification.result),
    .f = function(x, y) {
      data.frame(MS2.spectra.name = y,
                 Candidate.number = nrow(x))
    }
  ) %>%
    do.call(rbind, .)
  Candidate.number <-
    Candidate.number %>%
    dplyr::left_join(object@match.result[,c("MS1.peak.name", "MS2.spectra.name")],
                     by = "MS2.spectra.name") %>%
    dplyr::select(MS1.peak.name, dplyr::everything())

  identification.table <-
    identification.table %>%
    dplyr::left_join(Candidate.number,
                     by = c("name" = "MS1.peak.name")) %>%
    dplyr::select(colnames(object@ms1.data), MS2.spectra.name,
                  Candidate.number, dplyr::everything())

  identification.table$MS2.spectra.name[which(is.na(identification.table$Compound.name))] <- NA
  identification.table$Candidate.number[which(is.na(identification.table$Compound.name))] <- NA

  return(tibble::as_tibble(identification.table))
}
