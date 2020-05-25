source("reserving_functions.R")
 

# Input
payments <- read.csv(file = "triangles_input/triangle_theorique.csv", sep = ";", header = FALSE)

# Chain Ladder
implem_theorique <-  chainLadder_reserving(triangle = payments)
implem_theorique[[3]]
implem_theorique[[5]]
implem_theorique[[6]]

# Bootstrap Chain Ladder
N = 1000
implem_theorique_bootstrap <- bootstrap_cl_reserving(triangle = payments, number_of_simul = N)

# Reporting : Chain-ladder model ; bootstrap results
report <- data.frame("annee"=c("year_1", "year_2", "year_3", "year_4", "year_5", "year_6", "year_7", "year_8", "year_9", "year_10"))
report$chain_ladder <- implem_theorique[[5]]
report$bootstrap_mean <- rowMeans(implem_theorique_bootstrap[[7]][1:10,1,1:N])
report$bootstrap_standard_dev <- rowSds(implem_theorique_bootstrap[[7]][1:10,1,1:N]) # bootstrap ecart type
report$bootstrap_rmse <- sqrt(rowMeans((implem_theorique[[5]] - implem_theorique_bootstrap[[7]][1:10,1,1:N])^2)) #bootstrap erreur quadratique moyenne
#report$bootstrap_prediction_error_pct <- (abs(report$bootstrap_prediction_error - report$bootstrap_mean)/report$bootstrap_mean)*100

report