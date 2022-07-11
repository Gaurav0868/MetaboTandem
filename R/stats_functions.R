#' Global are under the curve normalization
#'
#' Function to apply global sum normalization to the data
#'
#'
#' @param matrix Matrix or data frame with the unnormalized peak intensities
#'               with peaks as rows and samples as columns
#' @param transform_data Logical value of whether log transform data or not
#'
#' @export
#'
global.norm <- function(matrix, transform_data = TRUE){
  colsum <- colSums(matrix, na.rm = TRUE)
  colsum.median <- median(colsum)
  norm.matrix <- data.frame(matrix(nrow = nrow(matrix),
                                   ncol = ncol(matrix)))
  for(col in 1:ncol(matrix)){
    norm.matrix[,col] <- (matrix[,col] / colsum[col]) * colsum.median
  }
  colnames(norm.matrix) <- colnames(matrix)
  rownames(norm.matrix) <- rownames(matrix)
  if(transform_data == TRUE){
    norm.matrix <- log2(norm.matrix + 1)
  }

  return(norm.matrix)
}


#' Median Normalization
#'
#' Function to apply global median normalization to the data
#'
#'
#' @param matrix Matrix or data frame with the unnormalized peak intensities
#'               with peaks as rows and samples as columns
#' @param transform_data Logical value of whether log transform data or not
#'
#' @export
#'
median.norm <- function(matrix, transform_data = TRUE){

  colmedian <- apply(matrix, 2, FUN = median, na.rm = TRUE)
  colmedian.mean <- mean(colmedian)
  norm.matrix <- data.frame(matrix(nrow = nrow(matrix),
                                   ncol = ncol(matrix)))
  for(col in 1:ncol(matrix)){
    norm.matrix[,col] <- (matrix[,col] / colmedian[col]) * colmedian.mean
  }
  colnames(norm.matrix) <- colnames(matrix)
  rownames(norm.matrix) <- rownames(matrix)
  if(transform_data == TRUE){
    norm.matrix <- log2(norm.matrix + 1)
  }
  return(norm.matrix)
}

#' Mean Normalzation
#' Function to apply global mean normalization to the data
#'
#' @param matrix Matrix or data frame with the unnormalized peak intensities
#'               with peaks as rows and samples as columns
#' @param transform_data Logical value of whether log transform data or not
#'
#' @export
#'
mean.norm <- function(matrix, transform_data = TRUE){

  colmean <- colMeans(matrix, na.rm = TRUE)
  colmean.mean <- mean(colmean)
  norm.matrix <- data.frame(matrix(nrow = nrow(matrix),
                                   ncol = ncol(matrix)))
  for(col in 1:ncol(matrix)){
    norm.matrix[,col] <- (matrix[,col] / colmean[col]) * colmean.mean
  }
  colnames(norm.matrix) <- colnames(matrix)
  rownames(norm.matrix) <- rownames(matrix)
  if(transform_data == TRUE){
    norm.matrix <- log2(norm.matrix + 1)
  }
  return(norm.matrix)
}

#' Variance stabilization normalization
#' Function to apply variance stabilization normalization to the data
#'
#' @param matrix Matrix or data frame with the unnormalized peak intensities with
#'               peaks as rows and samples as columns
#' @param transform_data Logical value of whether log transform data or not
#'
#' @export
#'
vsn.norm <- function(matrix){

  norm.matrix <- suppressMessages(vsn::justvsn(as.matrix(matrix)))
  norm.matrix <- as.data.frame(norm.matrix)
  return(norm.matrix)
}

#' Cyclic LOESS normalization
#'
#' Function to apply cyclic LOESS normalization to the data
#'
#'
#' @param matrix Matrix or data frame with the unnormalized peak intensities
#'               with peaks as rows and samples as columns
#' @param transform_data Logical value of whether log transform data or not
#'
#' @export
#'
cycloess.norm <- function(matrix){

  norm.matrix <- log2(matrix)
  norm.matrix <- limma::normalizeCyclicLoess(norm.matrix, method = 'fast')
  norm.matrix <- as.data.frame(norm.matrix)
  rownames(norm.matrix) <- rownames(matrix)
  return(norm.matrix)
}

#' Max normalization
#'
#' Function to apply normalization based on most intense peak to the data
#'
#'
#' @param matrix Matrix or data frame with the unnormalized peak intensities
#'               with peaks as rows and samples as columns
#' @param transform_data Logical value of whether log transform data or not
#'
#' @export
#'
max.norm <- function(matrix, transform_data = TRUE){

  colmax <- apply(matrix, 2, FUN = max, na.rm = TRUE)
  norm.matrix <- data.frame(matrix(NA, nrow = nrow(matrix), ncol = ncol(matrix)))
  for(col in 1:ncol(matrix)){
    norm.matrix[,col] <- (matrix[,col] / colmax[col])
  }
  colnames(norm.matrix) <- colnames(matrix)
  rownames(norm.matrix) <- rownames(matrix)
  if(transform_data == TRUE){
    norm.matrix <- log2(norm.matrix + 1)
  }
  return(norm.matrix)
}


#' Test normalization methods
#'
#' Function to test all normalization methods and compare graphically
#'
#'
#' @param df Matrix or data frame with the unnormalized peak intensities
#'               with peaks as rows and samples as columns
#'
#' @export
#'
normalize_by_all <- function(df){

  no_norm.df <- tidyr::pivot_longer(df,
                                    dplyr::everything(),
                                    names_to = 'SampleID',
                                    values_to = 'AUC')
  no_norm.plot <- plot_boxplot(no_norm.df, SampleID, AUC) +
    labs(title = 'No normalization') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())

  gi_norm.df <- global.norm(df)
  gi_norm.df <- tidyr::pivot_longer(gi_norm.df,
                                    dplyr::everything(),
                                    names_to = 'SampleID',
                                    values_to = 'AUC')
  gi_norm.plot <- plot_boxplot(gi_norm.df, SampleID, AUC) +
    labs(title = 'Normalization by global AUC',
         y =  'Normalized AUC') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())


  mean_norm.df <- mean.norm(df)
  mean_norm.df <- tidyr::pivot_longer(mean_norm.df,
                                      dplyr::everything(),
                                      names_to = 'SampleID',
                                      values_to = 'AUC')
  mean_norm.plot <- plot_boxplot(mean_norm.df, SampleID, AUC) +
    labs(title = 'Normalization by mean',
         y =  'Normalized AUC') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())


  median_norm.df <- median.norm(df)
  median_norm.df <- tidyr::pivot_longer(median_norm.df,
                                        dplyr::everything(),
                                        names_to = 'SampleID',
                                        values_to = 'AUC')
  median_norm.plot <- plot_boxplot(median_norm.df, SampleID, AUC) +
    labs(title = 'Normalization by median',
         y =  'Normalized AUC') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())


  vsn_norm.df <- vsn.norm(df)
  vsn_norm.df <- tidyr::pivot_longer(vsn_norm.df,
                                     dplyr::everything(),
                                     names_to =  'SampleID',
                                     values_to = 'AUC')
  vsn_norm.plot <- plot_boxplot(vsn_norm.df, SampleID, AUC) +
    labs(title = 'VSN',
         y =  'Normalized AUC') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())


  cycloess_norm.df <- cycloess.norm(df)
  cycloess_norm.df <- tidyr::pivot_longer(cycloess_norm.df,
                                          dplyr::everything(),
                                          names_to =  'SampleID',
                                          values_to = 'AUC')
  cycloess_norm.plot <- plot_boxplot(cycloess_norm.df, SampleID, AUC) +
    labs(title = 'LOESS normalization',
         y =  'Normalized AUC') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())


  all_norm.plot <- ggpubr::ggarrange(no_norm.plot, gi_norm.plot, mean_norm.plot,
                                     median_norm.plot, vsn_norm.plot, cycloess_norm.plot,
                                     nrow = 2,
                                     ncol = 3)

  return(all_norm.plot)

}

#' Function to apply any normalization method required
#'
#' Function to test all normalization methods and compare graphically
#'
#'
#' @param df Matrix or data frame with the unnormalized peak intensities
#'               with peaks as rows and samples as columns
#'
#'

norm_methods_all <- function(df, norm_method, transform = TRUE){
  switch(norm_method,
         global.norm = global.norm(df, transform_data = transform),
         median.norm = median.norm(df, transform_data = transform),
         mean.norm = mean.norm(df, transform_data = transform),
         vsn.norm = vsn.norm(df, transform_data = transform),
         cycloes.norm = cycloes.norm(df, transform_data = transform),
         none = df)
}

#' Function to calculate NMDS
#'
#' Function to test all normalization methods and compare graphically
#'
#'
#' @param df Matrix or data frame with the unnormalized peak intensities
#'               with peaks as rows and samples as columns
#'
#'
nmds_ordination <- function(abundance_matrix, metadata, mode = 'ra', color_by){
  if(mode == 'ra'){
    nmds.matrix <- t(abundance_matrix)
    dm.method <- 'bray'
    # distance matrix by Bray because relative abundance mode was selected
    dm <- vegan::vegdist(nmds.matrix, method=dm.method)
    print('Relative abundance method selected')
  }else if(mode == 'pa'){
    nmds.matrix <- vegan::decostand(t(abundance_matrix), 'pa')
    dm.method <- 'euclidean'
    dm <- vegan::vegdist(nmds.matrix, method = dm.method)
    print('Presence/absence method selected')
  } else{
    print('Select analysis method: "pa" for presence absence or
          "ra" for relative abundance')
    stop()
  }

  set.seed(123)
  nmds <- vegan::metaMDS(dm,
                  k = 2,
                  maxit = 999,
                  trymax = 500,
                  wascores = TRUE)

  nmds.scores <- as.data.frame(vegan::scores(nmds, display = 'sites')) %>%
    tibble::rownames_to_column(var = 'SampleID') %>%
    dplyr::left_join(metadata, by = 'SampleID')

  nmds_plot <- ggplot(nmds.scores,
                      aes(x = NMDS1,
                          y = NMDS2)) +
    geom_point(aes(color = .data[[color_by]])) +
    labs(title = 'NMDS plot') +
    theme_bw()+
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))


  res <- list(nmds = nmds,
              nmds_scores = nmds.scores,
              nmds_plot = nmds_plot)

  return(nmds)

}


pca_ordination <- function(abundance_matrix, metadata, color_by){

  # Calculate PCA with prcomp
  pca <- stats::prcomp(t(abundance_matrix))

  # Get eigenvalues
  eigen <- factoextra::get_eigenvalue(pca)
  # Plot screeplot using the functions from factoextra

  scree_plot <- factoextra::fviz_eig(pca, addlabels = TRUE) +
    theme_bw() +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))

  # Extract sample coordinates for PC1 and PC2
  pca_coordinates <- tibble::as_tibble(pca$x) %>%
    dplyr::mutate(SampleID = rownames(pca$x)) %>%
    dplyr::left_join(metadata, by ='SampleID')
  # Prepare axis labels for PCA

  pc1 <- paste0('PC1 (', round(eigen$variance.percent[1], digits = 1), '%)')
  pc2 <- paste0('PC2 (', round(eigen$variance.percent[2], digits = 1), '%)')

  # Plot Individuals PCA

  pca_plot <- ggplot(pca_coordinates,
                     aes(x = PC1,
                         y = PC2)) +
    geom_point(aes(color = .data[[color_by]])) +
    labs(title = 'PCA plot',
         x = pc1,
         y = pc2) +
    theme_bw()+
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))

  res <- list(pca_coordinates = pca_coordinates,
              scree_plot = scree_plot,
              pca_plot = pca_plot)

  return(res)
}

calculate_permanova <- function(abundance_matrix, metadata, mode = 'ra', group_by){

  if(mode == 'ra'){
    nmds.matrix <- t(abundance_matrix)
    dm.method <- 'bray'
    # distance matrix by Bray because relative abundance mode was selected
    dm <- vegan::vegdist(nmds.matrix, method=dm.method)
    print('Relative abundance method selected')
  }else if(mode == 'pa'){
    nmds.matrix <- vegan::decostand(t(abundance_matrix), 'pa')
    dm.method <- 'euclidean'
    dm <- vegan::vegdist(nmds.matrix, method = dm.method)
    print('Presence/absence method selected')
  } else{
    print('Select analysis method: "pa" for presence absence or "ra" for relative abundance')
    stop()
  }

  formula <- as.formula(paste0('dm ~ ', group_by))

  permanova <- vegan::adonis2(formula,
                              data=metadata,
                              permutations=999,
                              method=dm.method)

  return(permanova)
}

get_diff_table <- function(auc_matrix,
                           control.sample_list,
                           treatment.sample_list,
                           log2_transformed = FALSE){

  # Get the AUC values per each sample and calculate the means per feature
  temp.df_control <- auc_matrix %>%
    dplyr::select(all_of(control.sample_list))
  control_means <- rowMeans(temp.df_control, na.rm = TRUE)

  temp.df_treatment <- auc_matrix %>%
    dplyr::select(all_of(treatment.sample_list))
  treatment_means <- rowMeans(temp.df_treatment, na.rm = TRUE)

  diff_table <- as.data.frame(cbind(control_means, treatment_means))

  if(log2_transformed == TRUE){
    diff_table <- diff_table %>%
      dplyr::mutate(log2FC = treatment_means - control_means)
  } else {
    diff_table <- diff_table %>%
      dplyr::mutate(ratio = treatment_means/control_means) %>% # get the control/treatment ratio
      dplyr::mutate(log2FC = log2(ratio)) # calculate log2FC
  }

  rownames(diff_table) <- rownames(auc_matrix)

  # Initialize pvalues matrix
  pvalues <- data.frame(row.names = rownames(auc_matrix),
                        pval = rep(0, length(rownames(auc_matrix))))

  #Calculate pvalue per each of the features
  for(i in 1:nrow(pvalues)){
    t.test <- t.test(as.numeric(temp.df_control[i,]),
                     as.numeric(temp.df_treatment[i,]), paired = FALSE)
    pvalues$pval[i] <- t.test$p.value
  }
  pvalues <- pvalues %>%
    tibble::rownames_to_column(var = 'FeatureID')

  diff_table <- diff_table %>%
    tibble::rownames_to_column(var = 'FeatureID') %>%
    dplyr::left_join(pvalues, by = 'FeatureID')

  diff_table$pval.adj <- p.adjust(diff_table$pval, method = 'fdr')

  return(diff_table)

}

plot_volcano <- function(df,
                         log2FC,
                         pval,
                         log2FC.threshold,
                         pval.threshold){

  #Generate label for the plot

  significant_points <- df %>%
    dplyr::select(FeatureID, {{log2FC}}, {{pval}}) %>%
    dplyr::filter(abs({{log2FC}}) > log2FC.threshold,
           -log10({{pval}}) > -log10(0.05)) %>%
    dplyr::pull(FeatureID)

  plot <- df %>%
    dplyr::mutate(color4plot = ifelse(FeatureID %in% significant_points,
                                      'significant',
                                      'non-significant')) %>%
    ggplot(aes(x = {{log2FC}},
               y = -log10({{pval}}))) +
    geom_point(aes(color = color4plot)) +
    scale_color_manual(values = c("significant" = 'red',
                                  'non-significant' = 'black')) +
    geom_vline(xintercept = c(-{{log2FC.threshold}}, {{log2FC.threshold}}),
               linetype = 'dotted',
               size = 1,
               color = 'blue') +
    geom_hline(yintercept = -log10({{pval.threshold}}),
               linetype = 'dotted',
               size = 1,
               color = 'blue') +
    theme_bw() +
    labs(title = 'Volcano plot',
         x = expression("Log"[2]*" Fold Change"),
         y = expression("-Log"[10]*" pvalue")) +
    theme(plot.title = element_text(hjust = 0.5,
                                    face = 'bold'),
          plot.subtitle = element_text(hjust = 0.5,
                                       face = 'bold'),
          legend.position = 'none')

  return(plot)
}
