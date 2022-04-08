#A series of helper functions for R code.
#TJE 2021 02 22.
library(tidyverse)
library(cowplot)
library(ggpubr)
library(gtools)

##plot grid
plot_corr_grid <- function(df, remove_lower_tri = TRUE){
  #TJE 2021 02 22
  #assume all columns are numeric
  total_datasets <- ncol(df)
  total_possibilities_mat <- t(combn(total_datasets, 2))

  pList <- NULL
  name_vec <- colnames(df)
  ndata <- length(name_vec)

  permutations_mat <- permutations(n = ndata, r = 2,v = 1:ndata, repeats.allowed = TRUE) 
  
  for(i in 1:nrow(permutations_mat)){
    a = name_vec[permutations_mat[i,1]]
    b = name_vec[permutations_mat[i,2]]
    pList[[i]] <- ggplot(df, aes(x = !!ensym(a), y = !!ensym(b))) +
      geom_point(size = 1, alpha = 0.2, shape = 16) +
      stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top",aes(label = paste(..r.label..))) +
      scale_x_continuous(name = a, expand = c(0,0)) +
      scale_y_continuous(name = b, expand = c(0,0)) +
      geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey') + 
      theme_tim_label()
  }

  #Remove the upper triangle
  if(remove_lower_tri){
  	lower_tri_no <- matrix(1:ndata^2, ncol = ndata)[lower.tri(matrix(1:ndata^2, ncol = ndata))]
  	lower_tri_no <- c(lower_tri_no, diag(matrix(1:ndata^2, ncol = ndata)))
  	pList[lower_tri_no] <- list(NULL)
  }
  pGrid <- plot_grid(
    plotlist = pList,
    labels = NULL, 
    ncol = total_datasets, 
    nrow = total_datasets,
    byrow = FALSE)
  return(pGrid)
}

##plot grid
plot_corr_grid_PRESI <- function(df, remove_lower_tri = TRUE){
  #TJE 2021 02 22
  #assume all columns are numeric
  total_datasets <- ncol(df)
  total_possibilities_mat <- t(combn(total_datasets, 2))

  pList <- NULL
  name_vec <- colnames(df)
  ndata <- length(name_vec)

  permutations_mat <- permutations(n = ndata, r = 2,v = 1:ndata, repeats.allowed = TRUE) 
  
  for(i in 1:nrow(permutations_mat)){
    a = name_vec[permutations_mat[i,1]]
    b = name_vec[permutations_mat[i,2]]
    pList[[i]] <- ggplot(df, aes(x = !!ensym(a), y = !!ensym(b))) +
      geom_point(size = 1, alpha = 0.8, shape = 16) +
      stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top",aes(label = paste(..r.label..))) +
      scale_x_continuous(name = a, expand = c(0,0)) +
      scale_y_continuous(name = b, expand = c(0,0)) +
      geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey') + 
      theme_tim_presentation()
  }

  #Remove the upper triangle
  if(remove_lower_tri){
    lower_tri_no <- matrix(1:ndata^2, ncol = ndata)[lower.tri(matrix(1:ndata^2, ncol = ndata))]
    lower_tri_no <- c(lower_tri_no, diag(matrix(1:ndata^2, ncol = ndata)))
    pList[lower_tri_no] <- list(NULL)
  }
  pGrid <- plot_grid(
    plotlist = pList,
    labels = NULL, 
    ncol = total_datasets, 
    nrow = total_datasets,
    byrow = FALSE)
  return(pGrid)
}

cat(paste("Loading tidyverse, cowplot, ggpubr, and gtools.\n"))
cat(paste("Functions: plot_corr_grid(df)\n"))
cat(paste("Assume all columns are numeric\n"))
cat(paste("Pearson or spearman?\n"))