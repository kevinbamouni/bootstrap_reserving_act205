source("reserving_functions.R")
 
############################################### Input nécéssaires ##########################################################
 
# Fichier contenant le tringle input .csv
payments <- read.csv(file = "triangles_input/triangle_theorique.csv", sep = ";", header = FALSE)

# Type du triangle input : "increment" or "cumulative"
type = "increment" #or "cumulative"

implem_theorique = chainLadder_reserving(triangle = payments, type = "increment" )
