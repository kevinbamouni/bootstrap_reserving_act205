source("reserving_functions.R")
 

############ Input #############
payments <- read.csv(file = "triangles_input/triangle_theorique.csv", sep = ";", header = FALSE)
N = 1000 # Nombre de simulations bootstrap à effectuer

# Execution de Chain Ladder
implem_theorique <-  chainLadder_reserving(triangle = payments)
# Coef de developpements
implem_theorique[[3]]
# Provisions chain ladder par année 
implem_theorique[[5]]
# Provision  chain ladder totale
implem_theorique[[6]]

# Execution de Bootstrap Chain Ladder
implem_theorique_bootstrap <- bootstrap_cl_reserving(triangle = payments, number_of_simul = N)

# Reporting par année d'occurrence : Chain-ladder model ; bootstrap results
report <- data.frame("PROVISIONS"=c("year_1", "year_2", "year_3", "year_4", "year_5", "year_6", "year_7", "year_8", "year_9", "year_10"))
report$chain_ladder <- implem_theorique[[5]]
report$bootstrap_mean <- implem_theorique_bootstrap[[9]]
report$bootstrap_prediction_error <- implem_theorique_bootstrap[[10]]
report$bootstrap_prediction_error_pct <- report$bootstrap_prediction_error/report$bootstrap_mean*100
report

report_total <- data.frame("PROVISIONS"=c("Total"))
report_total$chain_ladder_total <- implem_theorique[[6]]
report_total$bootstrap_total_mean <- mean(implem_theorique_bootstrap[[8]])
report_total$bootstrap_total_prediction_error <- sd(implem_theorique_bootstrap[[8]])
report_total$bootstrap_total_prediction_error_pct <- report_total$bootstrap_total_prediction_error/report_total$bootstrap_total_mean*100
report_total

# Export CSV des résultats
write.csv(report,"reserving_output/reporting_par_an.csv", row.names = FALSE)
write.csv(report_total,"reserving_output/reporting_total.csv", row.names = FALSE)

