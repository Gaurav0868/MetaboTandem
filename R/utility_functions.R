#' Create color vectors
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

# Create boxplot
plot_boxplot <- function(df, my_x, my_y){
  ggplot(df,
         aes(x = {{my_x}},
             y = {{my_y}})) +
    geom_boxplot() +
    theme_bw() +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))
}

# get vector of samples
get_samples <- function(metadata.df, Treatment, value){
  # Get value to filter samples
  selector <- syms({{Treatment}})

  samples <- metadata.df %>%
    dplyr::filter((!!! selector) == value)

  # Get only sample names
  samples <- samples$SampleID

  return(samples)
}

