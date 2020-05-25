library(matrixStats)


# Fonction de calcul des provisions par la méthode de Chain ladder.
# Input : triangle, triangle des incréments, avec les années de survenance en ligne et les annees de developpement en colonnes
# Output: une liste de résultats
# 1 : triangle des incréments en input
# 2 : triangle des paiements cumulés calculée à partir dur triangle des incréments
# 3 : Coeffiscients de developpement
# 4 : Matrice, avec triangle inférieur des paiements cumulés estimé
# 5 : Vecteur des provision à consituer par année de survenance
# 6 : Provision totale à consituer, selon m
chainLadder_reserving <- function(triangle){
  
  I <- nrow(triangle)
  J <- ncol(triangle)
  
  triangle_cumulatif <- rowCumsums(data.matrix(triangle))
  
  coefficients_de_dev_CL <- calcul_coef_de_dev_CL(triangle_cumulatif)
  
  matrice_estimee = round(triangle_inf_estimation(triangle_cumulatif, coefficients_de_dev_CL))
  
  reserving_par_an = round(reserving_est_par_an(matrice_estimee))
  
  reserving_total_CL = sum(reserving_par_an)
  
  result_list <- list(triangle, triangle_cumulatif, coefficients_de_dev_CL, matrice_estimee, reserving_par_an, reserving_total_CL)
  
  return(result_list)
}



# Fonction de calcul de reserving par bootstrap (chain ladder)
# Input 1 (triangle) : triangle, triangle des incréments, avec les années de survenance en ligne et les annees de developpement en colonnes
# Input 2 (number_of_simul) : Nombre de simulation bootstrap à exécuter 
# 1 : Coeffiscientde de developpement de chain ladder
# 2 : Triangle des paiements cumulés calculée à partir dur triangle des incréments
# 3 : Triangle supérieur des incréments réestimés
# 4 : Triangle des Unscaled Pearson residuals
# 5 : Phi Pearson Scale Parameter
# 6 : Triangle des Adjust Unscaled Pearson residuals
# 7 : Matrices en 3D contenant les reserves par an à effectuer estimées pour chaque simulation bootsratp
# 8 : Reserve totale à effectuer pour chaque simulation
bootstrap_cl_reserving <- function(triangle, number_of_simul){
  
  I <- nrow(triangle)
  J <- ncol(triangle)
  
    triangle_cumulatif <- rowCumsums(data.matrix(triangle))
  
  #----------   Algorithme de boostrap décrite dans   ----------
  
  coefficients_de_dev_CL <- calcul_coef_de_dev_CL(triangle_cumulatif)
  
  est_cumul_past_triangle <- triangle_sup_cumul_estimation(triangle_cumulatif, coefficients_de_dev_CL)
  
  est_incre_past_triangle <- triangle_sup_increment_from_cumuls(est_cumul_past_triangle)
  
  uns_pear_res <- triangle_upr(triangle, est_incre_past_triangle)
  
  phi_c <- psp_phi(uns_pear_res)
  
  adj_uns_pear_res <- triangle_upr_adjust(uns_pear_res)
  
  reserve_total_cl <- c() 
  proj_matrix <- array(NA,dim = c(I,1,number_of_simul))
  
  for (l in 1:number_of_simul) {
    
    triangle_apr_resampled <- triangle_sup_by_resample(adj_uns_pear_res, replacement = TRUE)
    
    triangle_sup_incr_resampled <- triangle_apr_resampled * sqrt(est_incre_past_triangle) + est_incre_past_triangle
    
    bootstrapcl <- chainLadder_reserving(triangle = triangle_sup_incr_resampled)
    
    mat_toto <- triangle_sup_increment_from_cumuls(bootstrapcl[[4]])
    
    mat_toto <- triangle_inf_sample_norm(mat_toto, phi_c)
    
    proj_matrix[,,l] <- reserving_est_par_an(rowCumsums(mat_toto))
    
    reserve_total_cl <- c(reserve_total_cl, sum(proj_matrix[,,l]))
    
  }
  
  result_list <- list(coefficients_de_dev_CL, est_cumul_past_triangle, est_incre_past_triangle, uns_pear_res, phi_c, adj_uns_pear_res, proj_matrix, reserve_total_cl)
  
  return(result_list)
}



# Fonction de calcul des coeffisicients de developpement de Chain Ladder
# Input : triangle supérieur des paiements cumulés 
# Output : vecteur des coeffiscients de developpement selon la méthode chain ladder
calcul_coef_de_dev_CL <- function(triangle_des_cumuls){
  
  I <- nrow(triangle_des_cumuls)
  J <- ncol(triangle_des_cumuls)
  
  coeffs <- c()
  for (j in 1:(J-1)) {
    coeffs <- c(coeffs, sum(triangle_des_cumuls[,j+1], na.rm = TRUE)/sum(triangle_des_cumuls[1:(I-j),j]))
  }
  return(coeffs)
}


# Fonction d'estimation de la partie inferieure du triangle
# Input 1 (triangle_cumule) : triangle supérieur des paiements cumulés
# Input 2 (coef_de_dev) : vecteur Coeffiscients de developpement de chain ladder
# Output : Matrice avec triangle inférieur des paiements cumulés en ligne estimée selon chain ladder
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



# Fonction d'estimation des provisions à effectuer par année d'occurrence
# Input 1 (matrice_cumulative_estimee) : Matrice complete des paiements cumulés en lignes.
# Output : Vecteur des provisions à consituer par année de survenance
reserving_est_par_an <- function(matrice_cumulative_estimee){
  I <- nrow(matrice_cumulative_estimee)
  J <- ncol(matrice_cumulative_estimee)
  
  return(matrice_cumulative_estimee[,J] - diag(matrice_cumulative_estimee[,J:1]))
}



# Fonction d'estimation de la partie superieure du triangle des cumuls à l'aide des coeffiscients de developpements de chain ladder
# Input 1 (triangle_cumule) : Triangle supérieur des paiments cumulés réels en ligne
# Input 2 (coef_de_dev) : coeffiscients de developpement de chain ladder
# Output : triangle supérieur des paiements cumulés estimés
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



# Fonction de calcul par différenciation du triangle des incréments à partir du triangle des paiements cumulés.
# Input : Triangle supérieur des paiements cumulés en ligne
# Output : Triangle supérieur des incréments correspondant
triangle_sup_increment_from_cumuls <- function(triangle_cumule){
  
  return(cbind(triangle_cumule[,1],rowDiffs(triangle_cumule)))
}



# Fonction de calcul des unscaled pearson residuals upr
# Input 1 (triangle_inc) : Triangle supérieur des incréments réels
# Input 2 (triangle_inc_est) : Triangle supérieur des incréments réestimés
# Output : Triangle supérieur des Unscaled Pearson Résiduals
triangle_upr <- function(triangle_inc, triangle_inc_est){
  
  return((triangle_inc - triangle_inc_est) / sqrt(triangle_inc_est))
}



# Fonction de calcul pearson scale parameter phi psp_phi
# Input 1 (triangle) : triangle supérieur des Unscaled Pearson Résiduals
# Output : Pearson Scale Parameter 
psp_phi <-  function(triangle){
  
  I <- nrow(triangle)
  
  return(sum(triangle^2, na.rm = TRUE)/(0.5*I*(I+1)-2*I+1))
}



# Fonction de calcul des Adjust Unscaled Pearson Résiduals AUPR
# Input 1 (triangle) : triangle supérieur des Unscaled Pearson Résiduals
# Output : triangle supérieur des Adjust Unscaled Pearson Résiduals
triangle_upr_adjust <- function(triangle){
  
  I <- nrow(triangle)
  
  return( sqrt( I / (0.5*I*(I+1)-2*I+1) ) * triangle)
}



# Fonction de création d'un triangle supérieur par rééchantillonage avec remise
# Input 1 (triangle_sup_des_residus) : Triangle supérieur des résidus
# Input 2 (replacement) : Type du tirage à effectuer
# Output : Triangle Input 1 rééchantillonné avec ou sans remise
triangle_sup_by_resample <- function(triangle_sup_des_residus, replacement){
  
  I <- nrow(triangle_sup_des_residus)
  J <- ncol(triangle_sup_des_residus)
  
  vect <- as.numeric(na.omit(as.vector(data.matrix(triangle_sup_des_residus))))
  
  for (i in 1:I){
    triangle_sup_des_residus[i,1:(J-i+1)] <- sample(x = vect, size = (J-i+1), replace = replacement)
  }
  
  return(triangle_sup_des_residus)
}



# Fonction de rééchantillonnage du triangle inférieur suivant une loi normale
# Input 1 (triangle_inf) : Triangle inférieur des incréments
# Input 2 (phi) : Type du tirage à effectuer
# Output : rééchantillonnage de chaque élément C(i,j) de Input 1 suivant une loi normale de moyenne C(j,j) et variance phi*C(i,j)
triangle_inf_sample_norm <- function(triangle_inf, phi){
  
  I <- nrow(triangle_inf)
  J <- ncol(triangle_inf)
  i <- I

    for (j in 2:J) {
    triangle_inf[i,j:J] <- rnorm(length(triangle_inf[i,j:J]), mean = 0, sd=1) * sqrt(phi * abs(triangle_inf[i,j:J])) + triangle_inf[i,j:J]
   # browser()
    i = i-1
  }
  
  return(triangle_inf)
}

