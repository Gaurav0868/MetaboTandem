#' Annotate MS2 using SIRIUS
#'
#' Function to identify feature structure using SIRIUS
#' (https://bio.informatik.uni-jena.de/software/sirius/)
#'
#'@param ms1.data A [data.frame] containing the columns 'FeatureID', 'mzmed' and 'rtmed'
#'@param ms2.data An [MSpectra] object containing MS2 data

#'
#'@return A dataframe with the results of annotation using SIRIUS
#'
#'@export

annotate_SIRIUS <- function (ms2.data,
                             groupdir = 'grouped_mgf',
                             outdir = 'SIRIUS_output'){

  if(!dir.exists(outdir)) dir.create(outdir)

  separate_MS2(ms2.data, groupdir, 10)

  mgf_list <- list.files(groupdir, pattern = '*.mgf', full.names = TRUE)

  purrr::walk(mgf_list, function(file){
    name <- unlist(stringr::str_split(file, '\\/'))
    name <- stringr::str_remove(name[length(name)], '.mgf')
    outfile <- file.path(outdir, name)

    command_args <- c(paste0('-i ', file),
                      paste0('-o ', outfile),
                      'formula',
                      'fingerprint',
                      'structure',
                      'canopus',
                      'write-summaries')

    system2('sirius', args = command_args)

  })

  annotation_df <- merge_SIRIUS_files(outdir) %>%
    fix_entries() %>%
    tidyr::separate_rows('PubChem_CID', sep = ',')

  # Retrieving metabolite info from PubChem

  cid_list <- annotation_df %>%
    dplyr::distinct(PubChem_CID) %>%
    dplyr::filter(!is.na(PubChem_CID)) %>%
    dplyr::pull(PubChem_CID)

  httr::set_config(httr::config(http_version = 0))
  safe_GET <- purrr::safely(httr::GET)
  i <- 1

  info_type <- c('MolecularFormula',
                 'ExactMass',
                 'InchIKey')

  pubchem_info <- purrr::map(cid_list, function(id){

    searches <- purrr::map(info_type, function(info){

      url <- paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',
                    id,
                    '/property/',
                    info,
                    '/CSV')

      r <- safe_GET(url)

      if(is.null(r$error)){
        content <- try(httr::content(r$result, 'text', encoding = "UTF-8"))
        content <- stringr::str_remove_all(unlist(stringr::str_split(content, ','))[3],
                                           '\\\n|\\"|\\\\')
      } else {
        content <- r$error
      }

      i <<- i + 1
      if(i %% 5 == 0) Sys.sleep(1)

      return(content)
    })

    df <- tibble(CID = id,
                 MolecularFormula = searches[[1]],
                 ExactMass = searches[[2]],
                 InchIKey = searches[[3]])

    print(paste0(c('Finishing: CID- ', id)))

    return(df)
  })

  pubchem_info <- do.call(rbind, pubchem_info)

  final_df <- dplyr::left_join(annotation_df, pubchem_info,
                               by = c('PubChem_CID' = 'CID'))

  return(final_df)
}

#' Function to parse ms2 data to work with SIRIUS

parse_ms2_for_SIRIUS <- function(sp, TITLE = NULL) {

  data <- c("BEGIN IONS",
            paste0("SCANS=", MSnbase::acquisitionNum(sp)),
            paste0("TITLE=", TITLE))

  if(MSnbase::polarity(sp) == 1){
    signal = "+"
  } else signal = "-"

  if (MSnbase::msLevel(sp) > 1) {
    data <- c(data,
              paste0("RTINSECONDS=", MSnbase::rtime(sp)),
              paste0("PEPMASS=", MSnbase::precursorMz(sp)),
              paste0("MSLEVEL=", MSnbase::msLevel(sp)))

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

  return(data)

}

#' Separate MGF into groups
#'
#' This function separates MS2 data into chunks to use SIRIUS web services

separate_MS2 <- function(ms2.data, outdir, n){

  mgf <- purrr::map2(ms2.data@listData,
                            ms2.data@elementMetadata@listData[["feature_id"]],
                            ~ parse_ms2_for_SIRIUS(.x, TITLE = .y))

  mgf <- purrr::reduce(ms2.parsed, c)

  if(!dir.exists(outdir)){
    dir.create(outdir)
  }

  accum <- c()

  start <- 'FALSE'
  stop <-  'FALSE'
  group <- 1
  counter <- 0

  for(i in 1:length(mgf)){
    if(stringr::str_detect(mgf[i], 'BEGIN')){
      start <- "TRUE"
      stop <-  'FALSE'
      counter <- counter + 1
    }

    if(stringr::str_detect(mgf[i], 'END')){
      start <- "FALSE"
      stop <- "TRUE"
    }

    if(start == 'TRUE'){
      accum <- c(accum, mgf[i])
    }

    if(stop == 'TRUE'){
      accum <- c(accum, mgf[i])
      if(counter == n){
        counter <- 0
        readr::write_lines(accum, paste0(outdir, 'group_', group, '.mgf'))
        group <- group + 1
        accum <- c()
      }
    }
    if(i == length(mgf)){
      readr::write_lines(accum, paste0(outdir, 'group_', group, '.mgf'))
    }
  }

}

#' Merge SIRIUS annotation
#'
#' Function to merge the annotation files retrieved using SIRIUS


merge_SIRIUS_files <- function(annot_dir){

  id_files <- list.files(annot_dir,
                         pattern = 'compound_identifications.tsv',
                         recursive = TRUE,
                         full.names = TRUE)
  canopus_files <- list.files(annot_dir,
                              pattern = 'canopus_compound_summary.tsv',
                              recursive = TRUE,
                              full.names = TRUE)

  id_merged <- purrr::map(id_files, function(file){
    temp_df <- readr::read_tsv(file) %>%
      dplyr::mutate(FeatureID = stringr::str_extract(id, 'FT.*')) %>%
      dplyr::select(FeatureID, molecularFormula, adduct, InChI, smiles, links)
  })

  id_merged <- do.call(rbind, id_merged)

  canopus_merged <-purrr::map(canopus_files, function(file){
    temp_df <- readr::read_tsv(file) %>%
      dplyr::mutate(FeatureID = stringr::str_extract(id, 'FT.*')) %>%
      dplyr::select(FeatureID,
                    all_classifications = `ClassyFire#all classifications`,
                    superclass = `ClassyFire#superclass`,
                    class = `ClassyFire#class`,
                    subclass = `ClassyFire#subclass`)
  })

  canopus_merged <- do.call(rbind, canopus_merged)

  annotation_df <- dplyr::left_join(id_merged, canopus_merged, by = 'FeatureID')

  return(annotation_df)
}

#' Rearrange entries in annotation table
#'
#' Function to format annotation table

fix_entries <- function(annotation_df){
  fixed_df <- annotation_df %>%
    dplyr::mutate(HMDB = stringr::str_extract(links, 'HMDB:\\([0-9| ]+\\)'),
           HMDB = stringr::str_trim(stringr::str_remove_all(HMDB, 'HMDB:\\(|\\)'), side = 'both'),
           YMDB = stringr::str_extract(links, 'YMDB:\\([0-9| ]+\\)'),
           YMDB = stringr::str_remove_all(YMDB, 'YMDB:\\(|\\)'),
           KNApSAcK = stringr::str_extract(links, 'KNApSAcK:\\([0-9| ]+\\)'),
           KNApSAcK = stringr::str_remove_all(KNApSAcK, 'KNApSAcK:\\(|\\)'),
           CHEBI = stringr::str_extract(links, 'CHEBI:\\([0-9| ]+\\)'),
           CHEBI = stringr::str_remove_all(CHEBI, 'CHEBI:\\(|\\)'),
           PlantCyc = stringr::str_extract(links, 'Plantcyc:\\([A-Z|0-9| |-]+\\)'),
           PlantCyc = stringr::str_remove_all(PlantCyc, 'Plantcyc:\\(|\\)'),
           PlantCyc = stringr::str_replace_all(PlantCyc, ' ', ','),
           BioCyc = stringr::str_extract(links, 'Biocyc:\\([A-Z|0-9| |-]+\\)'),
           BioCyc = stringr::str_remove_all(BioCyc, 'Biocyc:\\(|\\)'),
           BioCyc = stringr::str_replace_all(BioCyc, ' ', ','),
           KEGG = stringr::str_extract(links, 'KEGG:\\(C[0-9| ]+\\)'),
           KEGG = stringr::str_remove_all(KEGG, 'KEGG:\\(|\\)'),
           KEGG = stringr::str_replace_all(KEGG, ' ', ','),
           COCONUT = stringr::str_extract(links, 'COCONUT:\\(CNP[0-9| ]+\\)'),
           COCONUT = stringr::str_remove_all(COCONUT, 'COCONUT:\\(|\\)'),
           COCONUT = stringr::str_replace_all(COCONUT, ' ', ','),
           PubChem_CID = stringr::str_extract(links, 'PubChem:\\([0-9| ]+\\)'),
           PubChem_CID = stringr::str_remove_all(PubChem_CID, 'PubChem:\\(|\\)'),
           PubChem_CID = stringr::str_replace_all(PubChem_CID, ' ', ','))

  for(i in 1:nrow(fixed_df)){
    if(!is.na(fixed_df$HMDB[i])){
      fixed_df$HMDB[i] <- paste(sprintf('HMDB%07d',
                                        as.numeric(unlist(stringr::str_split(fixed_df$HMDB[i], ' ')))),
                                collapse = ',')
    }
    if(!is.na(fixed_df$YMDB[i])){
      fixed_df$YMDB[i] <- paste(sprintf('YMDB%05d',
                                        as.numeric(unlist(stringr::str_split(fixed_df$YMDB[i], ' ')))),
                                collapse = ',')
    }
    if(!is.na(fixed_df$KNApSAcK[i])){
      fixed_df$KNApSAcK[i] <- paste(sprintf('C%08d',
                                            as.numeric(unlist(stringr::str_split(fixed_df$KNApSAcK[i], ' ')))),
                                    collapse = ',')
    }
    if(!is.na(fixed_df$CHEBI[i])){
      fixed_df$CHEBI[i] <- paste(paste0('CHEBI:',
                                        as.numeric(unlist(stringr::str_split(fixed_df$CHEBI[i], ' ')))),
                                 collapse = ',')
    }
  }

  return(fixed_df)
}
