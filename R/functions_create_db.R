#' Create massbank annotation database
#'
#' This function is intended to create a metID database using data from MassBank.
#' This function was modified from the `construct_massbank_database()` function from
#' package `metID`
#'
#' @param file Database file in mgf format
#' @param only_ms2 Keep only metabolites that have MS2 data
#' @param version Database version
#' @param link Link to the database
#' @param created_by Database creater
#' @param email Email of database creater
#' @param rt If `TRUE` database contains information of retention time
#' @param threads Number of parallel processes for reading database
#'
#' @export

create_massbank_db <- function(file,
                               only_ms2 = TRUE,
                               version = "0.0.1",
                               source = "MassBank",
                               link = "https://massbank.eu/MassBank/",
                               created_by = "MassBank",
                               email = "email",
                               rt = FALSE,
                               threads = parallelly::availableCores()){
  database <- read_msp(file, threads = threads)

  print('Startng processing database')
  all_metabolite_names <-
    purrr::map(database, function(x) {
      x$info$value[x$info$info == 'Name']
    }) %>%
    unlist() %>%
    unique()

  all_info_names <- purrr::map(
    database,
    function(x){
      x$info$info
    }
  ) %>%
    unlist() %>%
    unique()

  metabolite_info <- database %>%
    purrr::map(function(x) {
      x = as.data.frame(x$info)
      new_x = x[, 2]
      names(new_x) = x[, 1]
      new_x <- new_x[all_info_names]
      new_x
    })  %>%
    do.call(rbind, .) %>%
    as.data.frame()

  ###remove the metabolites without MS2 spectra
  if(only_ms2){
    remain_idx <- which(metabolite_info$Spectrum_type == "MS2")
    metabolite_info <- metabolite_info[remain_idx,]
    database = database[remain_idx]
  }

  metabolite_info <- metabolite_info %>%
    dplyr::select(Compound.name = Name,
                  mz = ExactMass,
                  Formula = Formula,
                  MassBank.ID = `DB#`,
                  dplyr::everything()) %>%
    dplyr::mutate(Lab.ID = paste("MassBank", 1:nrow(metabolite_info), sep = "_"),
                  RT = NA,
                  CAS.ID = NA,
                  HMDB.ID = NA,
                  KEGG.ID = NA,
                  mz.pos = NA,
                  mz.neg = NA,
                  Submitter = "MassBank",
                  Family = NA,
                  Sub.pathway = NA,
                  Note = NA) %>%
    dplyr::select(Lab.ID,
                  Compound.name,
                  mz,
                  RT,
                  CAS.ID,
                  HMDB.ID,
                  KEGG.ID,
                  Formula,
                  mz.pos,
                  mz.neg,
                  Submitter,
                  Family,
                  Sub.pathway,
                  Note,
                  dplyr::everything()
    )

  metabolite_info$Collision_energy[is.na(metabolite_info$Collision_energy)] = "not_available"
  metabolite_info$Collision_energy[metabolite_info$Collision_energy == ""] = "not_available"

  #####create metID database format
  positive_idx <- which(metabolite_info$Ion_mode == "POSITIVE")
  negative_idx <- which(metabolite_info$Ion_mode == "NEGATIVE")

  Spectra.positive <- database[positive_idx]
  Spectra.negative <- database[negative_idx]

  names(Spectra.positive) <- metabolite_info$Lab.ID[positive_idx]
  names(Spectra.negative) <- metabolite_info$Lab.ID[negative_idx]

  Spectra.positive <- purrr::map2(
    .x = Spectra.positive,
    .y = metabolite_info$Collision_energy[positive_idx],
    .f = function(x, y) {
      x = x$spec
      x = list(x)
      names(x) = y
      x
    }
  )

  Spectra.negative <- purrr::map2(
    .x = Spectra.negative,
    .y = metabolite_info$Collision_energy[negative_idx],
    .f = function(x, y) {
      x = x$spec
      x = list(x)
      names(x) = y
      x
    }
  )

  database.info <- list(
    "Version" = version,
    "Source" = source,
    "Link" = link,
    "Creater" = created_by,
    "Email" = email,
    "RT" = rt
  )

  Spectra <- list("Spectra.positive" = Spectra.positive,
                  "Spectra.negative" = Spectra.negative)

  database <- new(
    Class = "databaseClass",
    database.info = database.info,
    spectra.info = metabolite_info,
    spectra.data = Spectra
  )

  database@database.info$RT <-
    ifelse(all(is.na(database@spectra.info$RT)), FALSE, TRUE)
  print("Database done")
  return(database)
}

#' Create GNPS annotation database
#'
#' This function is intended to create a metID database using data from MassBank.
#' This function was modified from the `construct_massbank_database()` function from
#' package `metID`
#'
#' @param file Database file in mgf format
#' @param only_ms2 Keep only metabolites that have MS2 data
#' @param version Database version
#' @param link Link to the database
#' @param created_by Database creater
#' @param email Email of database creater
#' @param rt If `TRUE` database contains information of retention time
#' @param threads Number of parallel processes for reading database
#'
#' @export

create_gnps_db <- function(file,
                           only_ms2 = TRUE,
                           version = "0.0.1",
                           source = "GNPS",
                           link = "https://gnps-external.ucsd.edu/gnpslibrary/ALL_GNPS.msp",
                           created_by = "GNPS",
                           email = "email",
                           rt = FALSE,
                           threads = parallelly::availableCores()){
  database <- read_mgf(file, separator = "\t", threads = threads)

  print('Starting processing database')
  all_metabolite_names <-
    purrr::map(database, function(x) {
      x$info$value[x$info$info == 'NAME']
    }) %>%
    unlist() %>%
    unique()

  all_info_names <- purrr::map(
    database,
    function(x){
      x$info$info
    }
  ) %>%
    unlist() %>%
    unique()

  metabolite_info <- database %>%
    purrr::map(function(x) {
      x = as.data.frame(x$info)
      new_x = x[, 2]
      names(new_x) = x[, 1]
      new_x <- new_x[all_info_names]
      new_x
    })  %>%
    do.call(rbind, .) %>%
    as.data.frame()

  ###remove the metabolites without MS2 spectra
  if(only_ms2){
    remain_idx <- which(metabolite_info$MSLEVEL == "2")
    metabolite_info <- metabolite_info[remain_idx,]
    database = database[remain_idx]
  }

  metabolite_info <- metabolite_info %>%
    dplyr::select(Compound.name = NAME,
                  mz = PEPMASS,
                  GNPS.ID = SPECTRUMID,
                  dplyr::everything())

  httr::set_config(httr::config(http_version = 0))
  safe_GET <- purrr::safely(httr::GET)
  i <- 1
  formulas <- purrr::map(
    unique(metabolite_info$SMILES),
    function(x){
      if(x == "N/A" || x == ""){
        molfor <- tibble(SMILES = x, Formula = NA)
      } else {
        url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/",
                      "compound/smiles/property/MolecularFormula/TXT",
                      "?smiles=",
                      x)
        r <- safe_GET(url)

        if(is.null(r$error)){
          content <- try(httr::content(r$result, 'text', encoding = "UTF-8"))
          content <- stringr::str_remove(content, "\\\n")
        } else {
          content <- r$error
        }
        molfor <- tibble(SMILES = x, Formula = content)
      }

      print(paste(which(unique(metabolite_info$SMILES) == x), Sys.time(), content))
      i <<- i + 1

      if(i %% 5 == 0) Sys.sleep(1)
      return(molfor)
    }
  )

  formulas <- do.call(rbind, formulas)

  metabolite_info <- metabolite_info %>%
    left_join(formulas, by = 'SMILES') %>%
    dplyr::mutate(Lab.ID = paste("GNPS", 1:nrow(metabolite_info), sep = "_"),
                  RT = NA,
                  CAS.ID = NA,
                  HMDB.ID = NA,
                  KEGG.ID = NA,
                  mz.pos = NA,
                  mz.neg = NA,
                  Submitter = "MassBank",
                  Family = NA,
                  Sub.pathway = NA,
                  Note = NA) %>%
    dplyr::select(Lab.ID,
                  Compound.name,
                  mz,
                  RT,
                  CAS.ID,
                  HMDB.ID,
                  KEGG.ID,
                  Formula,
                  mz.pos,
                  mz.neg,
                  Submitter,
                  Family,
                  Sub.pathway,
                  Note,
                  dplyr::everything()
    )

  metabolite_info$Collision_energy = "not_available"

  #####create metID database format
  positive_idx <- which(metabolite_info$IONMODE == "Positive")
  negative_idx <- which(metabolite_info$IONMODE == "Negative")

  Spectra.positive <- database[positive_idx]
  Spectra.negative <- database[negative_idx]

  names(Spectra.positive) <- metabolite_info$Lab.ID[positive_idx]
  names(Spectra.negative) <- metabolite_info$Lab.ID[negative_idx]

  Spectra.positive <- purrr::map2(
    .x = Spectra.positive,
    .y = metabolite_info$Collision_energy[positive_idx],
    .f = function(x, y) {
      x = x$spec
      x = list(x)
      names(x) = y
      x
    }
  )

  Spectra.negative <- purrr::map2(
    .x = Spectra.negative,
    .y = metabolite_info$COLLISION_ENERGY[negative_idx],
    .f = function(x, y) {
      x = x$spec
      x = list(x)
      names(x) = y
      x
    }
  )

  database.info <- list(
    "Version" = version,
    "Source" = source,
    "Link" = link,
    "Creater" = created_by,
    "Email" = email,
    "RT" = rt
  )

  Spectra <- list("Spectra.positive" = Spectra.positive,
                  "Spectra.negative" = Spectra.negative)

  database <- new(
    Class = "databaseClass",
    database.info = database.info,
    spectra.info = metabolite_info,
    spectra.data = Spectra
  )

  database@database.info$RT <-
    ifelse(all(is.na(database@spectra.info$RT)), FALSE, TRUE)
  print("Database done")
  return(database)
}

#' Create MoNA annotation database
#'
#' This function is intended to create a metID database using data from MassBank of America.
#' This function was modified from the `construct_massbank_database()` function from
#' package `metID`
#'
#' @param file Database file in mgf format
#' @param only_ms2 Keep only metabolites that have MS2 data
#' @param version Database version
#' @param link Link to the database
#' @param created_by Database creater
#' @param email Email of database creater
#' @param rt If `TRUE` database contains information of retention time
#' @param threads Number of parallel processes for reading database
#'
#' @export

create_mona_db <- function(file,
                               only_ms2 = TRUE,
                               version = "0.0.1",
                               source = "MassBank of America",
                               link = "https://mona.fiehnlab.ucdavis.edu/",
                               created_by = "MoNA",
                               email = "email",
                               rt = FALSE,
                               threads = parallelly::availableCores()){
  database <- read_msp(file, threads = threads)

  print('Startng processing database')
  all_metabolite_names <-
    purrr::map(database, function(x) {
      x$info$value[x$info$info == 'Name']
    }) %>%
    unlist() %>%
    unique()

  all_info_names <- purrr::map(
    database,
    function(x){
      x$info$info
    }
  ) %>%
    unlist() %>%
    unique()

  null <- unlist(purrr::map(database, is.null))

  database <- database[!null]

  metabolite_info <- database %>%
    purrr::map(function(x) {
      x = as.data.frame(x$info) %>%
        group_by(info) %>%
        summarise(value = paste(value, collapse = ';'))
      new_x = x$value
      names(new_x) = x$info
      new_x <- new_x[all_info_names]
      new_x
    })  %>%
    do.call(rbind, .) %>%
    as.data.frame()

  ###remove the metabolites without MS2 spectra
  if(only_ms2){
    remain_idx <- which(metabolite_info$Spectrum_type == "MS2")
    metabolite_info <- metabolite_info[remain_idx,]
    database = database[remain_idx]
  }

  metabolite_info <- metabolite_info %>%
    dplyr::select(Compound.name = Name,
                  mz = ExactMass,
                  Formula = Formula,
                  MoNA.ID = `DB#`,
                  dplyr::everything()) %>%
    dplyr::mutate(Lab.ID = paste("MoNA", 1:nrow(metabolite_info), sep = "_"),
                  RT = NA,
                  CAS.ID = NA,
                  HMDB.ID = NA,
                  KEGG.ID = NA,
                  mz.pos = NA,
                  mz.neg = NA,
                  Submitter = "MoNA",
                  Family = NA,
                  Sub.pathway = NA,
                  Note = NA) %>%
    dplyr::select(Lab.ID,
                  Compound.name,
                  mz,
                  RT,
                  CAS.ID,
                  HMDB.ID,
                  KEGG.ID,
                  Formula,
                  mz.pos,
                  mz.neg,
                  Submitter,
                  Family,
                  Sub.pathway,
                  Note,
                  dplyr::everything()
    )

  metabolite_info$Collision_energy[is.na(metabolite_info$Collision_energy)] = "not_available"
  metabolite_info$Collision_energy[metabolite_info$Collision_energy == ""] = "not_available"

  #####create metID database format
  positive_idx <- which(metabolite_info$Ion_mode == "P")
  negative_idx <- which(metabolite_info$Ion_mode == "N")

  Spectra.positive <- database[positive_idx]
  Spectra.negative <- database[negative_idx]

  names(Spectra.positive) <- metabolite_info$Lab.ID[positive_idx]
  names(Spectra.negative) <- metabolite_info$Lab.ID[negative_idx]

  Spectra.positive <- purrr::map2(
    .x = Spectra.positive,
    .y = metabolite_info$Collision_energy[positive_idx],
    .f = function(x, y) {
      x = x$spec
      x = list(x)
      names(x) = y
      x
    }
  )

  Spectra.negative <- purrr::map2(
    .x = Spectra.negative,
    .y = metabolite_info$Collision_energy[negative_idx],
    .f = function(x, y) {
      x = x$spec
      x = list(x)
      names(x) = y
      x
    }
  )

  database.info <- list(
    "Version" = version,
    "Source" = source,
    "Link" = link,
    "Creater" = created_by,
    "Email" = email,
    "RT" = rt
  )

  Spectra <- list("Spectra.positive" = Spectra.positive,
                  "Spectra.negative" = Spectra.negative)

  database <- new(
    Class = "databaseClass",
    database.info = database.info,
    spectra.info = metabolite_info,
    spectra.data = Spectra
  )

  database@database.info$RT <-
    ifelse(all(is.na(database@spectra.info$RT)), FALSE, TRUE)
  print("Database done")
  return(database)
}



