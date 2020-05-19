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
  matrice_estimee = triangle_inf_estimation(triangle_cumulatif, coefficients_de_dev_CL)
  
  reserving_par_an = reserving_est_par_an(matrice_estimee)
  
  reserving_total_CL = sum(reserving_par_an)
  
  result_list <- list(triangle, type, triangle_cumulatif, coefficients_de_dev_CL, matrice_estimee, reserving_par_an, reserving_total_CL)
  
  return(result_list)
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


########################################################################
# Fonction d'estimation de la partie inferieure du triangle
########################################################################

triangle_inf_estimation <- function(triangle_cumule, coef_de_dev){
  
  I <- nrow(triangle_cumule)
  J <- ncol(triangle_cumule)
  i <- I
  
  for (j in 2:J) {
    triangle_cumule[i,j:J] <- cumprod(coef_de_dev[(j-1):length(coef_de_dev)]) * triangle_cumule[i,(j-1)]
    i = i-1
  }
  
  return(triangle_cumule)
}


########################################################################
# Fonction d'estimation des provisions à effectuer par année d'occurrence
########################################################################
reserving_est_par_an <- function(matrice_cumulative_estimee){
  I <- nrow(matrice_cumulative_estimee)
  J <- ncol(matrice_cumulative_estimee)
  
  return(matrice_cumulative_estimee[,J] - diag(matrice_cumulative_estimee[,J:1]))
}