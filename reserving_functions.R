library(tidyverse)
library(matrixStats)


##################################################
# Fonction de calcul de reserving par Chain ladder
##################################################

chainLadder_reserving <- function(triangle, type){
  
  # Récupération du nombre de lignes et du nombre de colonne du triangle
  I <- nrow(triangle)
  J <- ncol(triangle)
  
  # Calcul de la matrice des paiements cumulatifs en fonction du type de la matrice d'entrée : cumuls ou incréments
  if (type == "increment") {
    triangle_cumulatif <- rowCumsums(data.matrix(triangle))
  } 
  
  if (type == "cumulative"){
    triangle_cumulatif <- data.matrix(triangle)
  }
  
  # Calcul des coeffiscients de developpements de Chain Ladder
  coefficients_de_dev_CL <- calcul_coef_de_dev_CL(triangle_cumulatif)
  
  # Matrice complete avec l'estimation de la partie inférieure avec les coeffiscients de developpement
  matrice_estimee = 0
  
  reserving_par_an = 0
  
  reserving_total_CL = 0
  
  result_list <- list(triangle, type, triangle_cumulatif, coefficients_de_dev_CL, matrice_estimee, reserving_par_an, reserving_total_CL)

}


########################################################################
# Fonction de calcul des coeffisicients de developpement de Chain Ladder
########################################################################
calcul_coef_de_dev_CL <- function(triangle_des_cumuls){
  
  I <- nrow(triangle_des_cumuls)
  J <- ncol(triangle_des_cumuls)
  
  coeffs <- c()
  for (j in 1:(J-1)) {
    coeffs <- c(coeffs, sum(triangle_des_cumuls[,j+1], na.rm = TRUE)/sum(triangle_des_cumuls[1:(I-j),j]))
  }
  
  return(coeffs)
}