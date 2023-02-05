library(MASS)
library(kappalab)
library(igraph)

#' Function to replace integer with 1
#' replace all alteration count with 1
#'
#' @param col column
binary_func <- function(col) {
  index.1 <-  which(col != 0)
  col[index.1] <- 1
  return(col)
}


#' Function to compute unanimity score, 
#' divides all 1's in the column
#' by the sum of all 1's in the column.
#' Unanimity coefficicent / dividend for the game.
#' Adapted from Moretti (2008).
#'
#' @param col column
compute_unanimity <- function(col) {
  
  # divide 1 by the sum of all the ones in the column
  colsum_divided <- 1 / sum(col)
  
  # replace all 1's with colsum_divided
  one_index <- which(col == 1)
  col[one_index] <- colsum_divided
  
  return(col)
}


#' Function to compute shapley value.
#'
#' @param mat boolean matrix
compute_shapley <- function(mat){
  
  # compute unanimity coefficient
  unanimity_mat <- apply(mat, MARGIN = 2, FUN = compute_unanimity)
  
  # compute the row sum of unanimity coefficients
  row_sum_unanimity <- apply(unanimity_mat, MARGIN = 1, FUN = sum)
  
  # divide through by normalizing value
  shapley_values <- row_sum_unanimity / dim(mat)[2]
  
  return(shapley_values)
}


#' Function to compute game theoretic centrality
#'
#' @param graph igraph graph
#' @param k vector of weights
gamek <- function(graph, k){
  C <- 0
  shapleyvalue <- rep(NA,length(V(graph)))
  for (i in V(graph)){
    for(u in neighborhood(graph,1,i)){
      A <- k[u] / (1 + degree(graph)[u])
      C <- sum(A)
    }
    shapleyvalue[i] <- C
  }
  return(shapleyvalue)
}


##################################################
#       Toy data (Figure 6 in the paper)         #
##################################################

# Generate graph (figure 4)
sample_nodes <- data.frame(name=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14))

sample_relations <- data.frame(from=c(1,1,1,2,3,3,4,4,4,11,12),
                               to=c(2,3,4,5,6,7,8,9,10,12,13))

g <- graph_from_data_frame(sample_relations, 
                           directed=FALSE, 
                           vertices=sample_nodes)

# Boolean matrix (figure 5)
B <- matrix(
  c(0,1,0,0,1,
    0,0,1,0,0,
    0,1,1,0,0,
    0,0,1,0,0,
    1,0,1,0,1,
    0,1,1,1,0,
    0,0,0,1,0,
    0,0,0,1,0,
    1,0,0,0,1,
    0,0,0,1,1,
    0,1,0,0,1,
    0,0,1,0,0,
    0,0,0,1,0,
    1,1,1,0,0),
  nrow=14,
  ncol=5,
  byrow = TRUE
)

# Table A: Microarray Game
shap <- compute_shapley(B)
shap_frac <- fractions(shap)
genes <- seq(1,14,1)
microgame <- data.frame(genes=genes, shap=shap)
microgame_sort <- microgame[order(microgame$shap, decreasing = TRUE),]
shap_frac_sort <- fractions(microgame_sort$shap)
print(microgame_sort)

# Table B: Game Theoretic Centrality
g_shap <- gamek(g, shap)
genes <- seq(1,14,1)
graph_game <- data.frame(genes=genes, shap=g_shap)
graph_game_sort <- graph_game[order(graph_game$shap, decreasing = TRUE),]
g_shap_frac <- fractions(graph_game_sort$shap)
print(graph_game_sort)