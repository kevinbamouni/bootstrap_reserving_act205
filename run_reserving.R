source("reserving_functions.R")
 
############################################### Input nécéssaires ##########################################################
 
# Fichier contenant le tringle input .csv
payments <- read.csv(file = "triangles_input/triangle_theorique.csv", sep = ";", header = FALSE)

implem_theorique <-  chainLadder_reserving(triangle = payments)

implem_theorique_bootstrap <- bootstrap_cl_reserving(triangle = payments, number_of_simul = 1000)