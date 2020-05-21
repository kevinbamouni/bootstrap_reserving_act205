library(tidyverse)
library(matrixStats)


##################################################################################
# Fonction de calcul de reserving par Chain ladder
##################################################################################
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


##################################################################################
# Fonction de calcul de reserving par bootstrap (chain ladder)
##################################################################################
bootstrap_cl_reserving <- function(triangle, type, number_of_simul){
  
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
  
  ######### Algorithme de boostrap décrite dans 
  # Calcul des coeffiscients de developpements de Chain Ladder
  coefficients_de_dev_CL <- calcul_coef_de_dev_CL(triangle_cumulatif)
  
  # Estimation du triangle supérieur des cumuls à partir des coef de developpement de chain ladder
  est_cumul_past_triangle <- triangle_sup_cumul_estimation(triangle_cumulatif, coefficients_de_dev_CL)
  
  # Estimation du triangle du triangle supérieur des incrément à partir du triangle des cumuls de l'étape précédente
  est_incre_past_triangle <- triangle_sup_increment_from_cumuls(est_cumul_past_triangle)
  
  # Calcul des unscaled perason residuals
  uns_pear_res <- triangle_upr(triangle, est_incre_past_triangle)
  
  # Calcul du pearson scale parameter phi
  phi <- psp_phi(uns_pear_res)
  
  # Ajustement des pearson residuals
  adj_uns_pear_res <- triangle_upr_adjust(uns_pear_res)
  
  #Iterative loop
  
  for (l in 1:number_of_simul) {
    triangle_apr_resampled <- triangle_sup_by_resample(adj_uns_pear_res, replacement = TRUE)
    triangle_sup_incr_resampled <- triangle_apr_resampled * sqrt(est_incre_past_triangle) + est_incre_past_triangle
    bootstrapcl <- chainLadder_reserving(triangle = triangle_sup_incr_resampled, type = "increment" )
  }
  
  result_list <- list(coefficients_de_dev_CL, est_cumul_past_triangle, est_incre_past_triangle, uns_pear_res, phi, adj_uns_pear_res, bootstrapcl)
  
  return(result_list)
}


##################################################################################
# Fonction de calcul des coeffisicients de developpement de Chain Ladder
##################################################################################
calcul_coef_de_dev_CL <- function(triangle_des_cumuls){
  
  I <- nrow(triangle_des_cumuls)
  J <- ncol(triangle_des_cumuls)
  
  coeffs <- c()
  for (j in 1:(J-1)) {
    coeffs <- c(coeffs, sum(triangle_des_cumuls[,j+1], na.rm = TRUE)/sum(triangle_des_cumuls[1:(I-j),j]))
  }
  return(coeffs)
}


##################################################################################
# Fonction d'estimation de la partie inferieure du triangle
##################################################################################

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


##################################################################################
# Fonction d'estimation des provisions à effectuer par année d'occurrence
##################################################################################
reserving_est_par_an <- function(matrice_cumulative_estimee){
  I <- nrow(matrice_cumulative_estimee)
  J <- ncol(matrice_cumulative_estimee)
  
  return(matrice_cumulative_estimee[,J] - diag(matrice_cumulative_estimee[,J:1]))
}


##################################################################################
# Fonction d'estimation de la partie superieure du triangle des cumuls
##################################################################################
triangle_sup_cumul_estimation <- function(triangle_cumule, coef_de_dev){
  
  I <- nrow(triangle_cumule)
  J <- ncol(triangle_cumule)
  i <- 1
  
  for (j in 1:(J-1)) {
    triangle_cumule[i,1:(J-j)] <- triangle_cumule[i,(J-j+1)] / cumprod(coef_de_dev[(J-j):1])[length(coef_de_dev[(J-j):1]):1] 
    i = i+1
  }
  
  return(triangle_cumule)
}


##################################################################################
# Fonction de calcul du triangle des incréments à partir du triangle des cumuls.
##################################################################################
triangle_sup_increment_from_cumuls <- function(triangle_cumule){
  
  return(cbind(triangle_cumule[,1],rowDiffs(triangle_cumule)))
}


##################################################################################
# Fonction de calcul des unscaled pearson residuals upr
##################################################################################
triangle_upr <- function(triangle_inc, triangle_inc_est){
  
  return((triangle_inc - triangle_inc_est) / sqrt(triangle_inc_est))
}


##################################################################################
# Fonction de calcul pearson scale parameter phi psp_phi
##################################################################################
psp_phi <-  function(triangle_upr){
  
  I <- nrow(triangle_upr)
  
  return(sum(triangle_upr^2, na.rm = TRUE)/(0.5*I*(I+1)-2*I+1))
}


##################################################################################
# Fonction de calcul des adjust unscaled pearson residuals upr
##################################################################################
triangle_upr_adjust <- function(triangle_upr){
  
  I <- nrow(triangle_upr)
  
  return( sqrt( I / (0.5*I*(I+1)-2*I+1) ) * triangle_upr)
}

##################################################################################
# Fonction de calcul d'un triangle supérieur par rééchantillonage avec remise
##################################################################################
triangle_sup_by_resample <- function(triangle_sup_des_residus, replacement){
  
  I <- nrow(triangle_sup_des_residus)
  J <- ncol(triangle_sup_des_residus)
  
  vect <- as.numeric(na.omit(as.vector(data.matrix(triangle_sup_des_residus))))
  
  for (i in 1:I){
    triangle_sup_des_residus[i,1:(J-i+1)] <- sample(x = vect, size = (J-i+1), replace = replacement)
  }
  
  return(triangle_sup_des_residus)
}
