# Create color vectors
create_col_vector <- function(metadata, color_by = 'AUTO'){

  if(color_by == 'AUTO'){

    color_vector <- ggpubr::get_palette(palette = 'd3',
                                        length(metadata$SampleID))
    names(color_vector) <- metadata$SampleID

  } else{

    treatments <- dplyr::pull(metadata, color_by)
    colors <- ggpubr::get_palette(palette = 'd3',
                                  length(unique(treatments)))
    names(colors) <- unique(treatments)

    color_vector <- sapply(treatments,
                           function(x)colors[x])

    names(color_vector) <- treatments
  }

  return(color_vector)
}

getVolumes <- function(){
  osSystem <- Sys.info()["sysname"]
  if (osSystem == "Darwin") {
    volumes <- list.files("/Volumes/", full.names = T)
    names(volumes) <- basename(volumes)
  }
  else if (osSystem == "Linux") {
    volumes <- c(Computer = "/")
    media <- list.files("/media/", full.names = T)
    names(media) <- basename(media)
    volumes <- c(volumes, media)
  }
  else if (osSystem == "Windows") {
    volumes <- system("wmic logicaldisk get Caption", intern = T)
    volumes <- sub(" *\\r$", "", volumes)
    keep <- !tolower(volumes) %in% c("caption", "")
    volumes <- volumes[keep]
    volNames <- system("wmic logicaldisk get VolumeName",
                       intern = T)
    volNames <- sub(" *\\r$", "", volNames)
    volNames <- volNames[keep]
    volNames <- paste0(volNames, ifelse(volNames == "", "",
                                        " "))
    volNames <- paste0(volNames, "(", volumes, ")")
    names(volumes) <- volNames
  }
  else {
    stop("unsupported OS")
  }

  return(volumes)
}
