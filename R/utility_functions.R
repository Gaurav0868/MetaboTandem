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
