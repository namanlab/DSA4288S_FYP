library(tidyverse)
library(stringr)

# Function to calculate mutual information
MutInf <- function(table) {
  rowsums <- rowSums(table)
  colsums <- colSums(table)
  ent_x <- sum(-rowsums/sum(rowsums) * log2(rowsums/sum(rowsums)), na.rm = TRUE)
  ent_y <- sum(-colsums/sum(colsums) * log2(colsums/sum(colsums)), na.rm = TRUE)
  return(ent_x + ent_y - entropy_pair(table))
}

# Function to calculate entropy of a pair
entropy_pair <- function(table) {
  total <- sum(table)
  probs <- table / total
  return(sum(-probs * log2(probs), na.rm = TRUE))
}

assoc_val <- function(tab){
  tab[1, 1]*tab[2, 2] - tab[1, 2]*tab[2, 1]
}

tm <- matrix(c(1198, 557, 1493, 1278), ncol = 2, byrow = T)
assoc_val(tm) # 699443

ta <- matrix(c(512, 89, 313, 19), ncol = 2, byrow = T)
assoc_val(ta) # -18129

tb <- matrix(c(353, 17, 207, 8), ncol = 2, byrow = T)
assoc_val(tb) # - 695

