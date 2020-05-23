library(matrixStats)



# Fonction de calcul de reserving par Chain ladder
chainLadder_reserving <- function(triangle){
  
  # Récupération du nombre de lignes et du nombre de colonne du triangle
  I <- nrow(triangle)
  J <- ncol(triangle)
  
  # Calcul de la matrice des paiements cumulatifs en fonction du type de la matrice d'entrée : cumuls ou incréments
    triangle_cumulatif <- rowCumsums(data.matrix(triangle))
  
  # Calcul des coeffiscients de developpements de Chain Ladder
  coefficients_de_dev_CL <- calcul_coef_de_dev_CL(triangle_cumulatif)
  
  # Matrice complete avec l'estimation de la partie inférieure avec les coeffiscients de developpement
  matrice_estimee = triangle_inf_estimation(triangle_cumulatif, coefficients_de_dev_CL)
  
  reserving_par_an = reserving_est_par_an(matrice_estimee)
  
  reserving_total_CL = sum(reserving_par_an)
  
  result_list <- list(triangle, triangle_cumulatif, coefficients_de_dev_CL, matrice_estimee, reserving_par_an, reserving_total_CL)
  
  return(result_list)
}



# Fonction de calcul de reserving par bootstrap (chain ladder)
bootstrap_cl_reserving <- function(triangle, number_of_simul){
  
  # Récupération du nombre de lignes et du nombre de colonne du triangle
  I <- nrow(triangle)
  J <- ncol(triangle)
  
  # Calcul de la matrice des paiements cumulatifs en fonction du type de la matrice d'entrée : cumuls ou incréments
    triangle_cumulatif <- rowCumsums(data.matrix(triangle))
  
  #----------   Algorithme de boostrap décrite dans   ----------
  
  # Calcul des coeffiscients de developpements de Chain Ladder
  coefficients_de_dev_CL <- calcul_coef_de_dev_CL(triangle_cumulatif)
  
  # Estimation du triangle supérieur des cumuls à partir des coef de developpement de chain ladder
  est_cumul_past_triangle <- triangle_sup_cumul_estimation(triangle_cumulatif, coefficients_de_dev_CL)
  
  # Estimation des incréments à partir du triangle des cumuls de l'étape précédente
  est_incre_past_triangle <- triangle_sup_increment_from_cumuls(est_cumul_past_triangle)
  
  # Calcul du triangle unscaled perason residuals
  uns_pear_res <- triangle_upr(triangle, est_incre_past_triangle)
  
  # Calcul du pearson scale parameter phi
  phi_c <- psp_phi(uns_pear_res)
  
  # calcul du triangle des adjust pearson residuals
  adj_uns_pear_res <- triangle_upr_adjust(uns_pear_res)
  
  #Iterative loop
  reserve_total_cl <- c() #stockage des provision totales pour chaque simulation de bootstrap
  proj_matrix <- array(NA,dim = c(I,1,number_of_simul)) #stockage des matrices incréments futures estimées par chain ladder pour chaque simulation de bootstrap
  
  for (l in 1:number_of_simul) {
    
    #Rééchantillonnage du triangle des adjust pearson residus
    triangle_apr_resampled <- triangle_sup_by_resample(adj_uns_pear_res, replacement = TRUE)
    
    #calcul du triangle supérieur d'incréments correstpondant 
    triangle_sup_incr_resampled <- triangle_apr_resampled * sqrt(est_incre_past_triangle) + est_incre_past_triangle
    
    #Application de chain ladder sur le triangle calculé à l'étape précédente
    bootstrapcl <- chainLadder_reserving(triangle = triangle_sup_incr_resampled)
    
    #stockage de la matrice des incréments obtenues à partir de la matrice complètée des cumuls estimées par chain ladder
    mat_toto <- triangle_sup_increment_from_cumuls(bootstrapcl[[4]])
    
    #reechantillonnage du triangle inferieur suivant une loi normale
    mat_toto <- triangle_inf_sample_norm(mat_toto, phi_c)
    
    #Calcul des reserves par annee de survenance
    proj_matrix[,,l] <- reserving_est_par_an(rowCumsums(mat_toto))
    
    #stockage de la provision totale estimée par chain ladder
    reserve_total_cl <- c(reserve_total_cl, sum(proj_matrix[,,l]))
    
  }
  
  result_list <- list(coefficients_de_dev_CL, est_cumul_past_triangle, est_incre_past_triangle, uns_pear_res, phi_c, adj_uns_pear_res, proj_matrix, reserve_total_cl)
  
  return(result_list)
}



# Fonction de calcul des coeffisicients de developpement de Chain Ladder
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
reserving_est_par_an <- function(matrice_cumulative_estimee){
  I <- nrow(matrice_cumulative_estimee)
  J <- ncol(matrice_cumulative_estimee)
  
  return(matrice_cumulative_estimee[,J] - diag(matrice_cumulative_estimee[,J:1]))
}



# Fonction d'estimation de la partie superieure du triangle des cumuls
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



# Fonction de calcul du triangle des incréments à partir du triangle des cumuls.
triangle_sup_increment_from_cumuls <- function(triangle_cumule){
  
  return(cbind(triangle_cumule[,1],rowDiffs(triangle_cumule)))
}



# Fonction de calcul des unscaled pearson residuals upr
triangle_upr <- function(triangle_inc, triangle_inc_est){
  
  return((triangle_inc - triangle_inc_est) / sqrt(triangle_inc_est))
}



# Fonction de calcul pearson scale parameter phi psp_phi
psp_phi <-  function(triangle){
  
  I <- nrow(triangle)
  
  return(sum(triangle^2, na.rm = TRUE)/(0.5*I*(I+1)-2*I+1))
}



# Fonction de calcul des adjust unscaled pearson residuals upr
triangle_upr_adjust <- function(triangle){
  
  I <- nrow(triangle)
  
  return( sqrt( I / (0.5*I*(I+1)-2*I+1) ) * triangle)
}



# Fonction de création d'un triangle supérieur par rééchantillonage avec remise
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



# calcul des provisions par années


