# bootstrap_reserving_act205 (implémentation en R).

## Packages R nécéssaires

- matrixStats

## Détails des fichiers et repertoire

- Le fichier reserving_functions.R sert de bibliothèques de fonctions. Son execution ne produit aucun résultat mais charges les fonctions utiliser pour les calculs de chain ladder et de bootstrap chain ladder

- Le fichier run_reserving.R est le fichier qui contient le code qui fait appel aux fonctions du fichier reserving_functions.R pour calculer les provisions en utilisant les méthodes chain ladder et bootstrap chain ladder. Son execution, prend un input et produit les calculs de provisions pour écrire les résultats dans des fichiers CSV. Vous pouvez le modifier pour charger d'autres inputs et choisir le nombre N de simulation de bootstrap à effectuer.

- Le repertoire triangle_input/ contient les triangles en input pour le calcul de provisions au format csv. En ligne sont les années d'occurence et en colonnes les années de developpement. la partie supérieur du triangle doit être complétée. Les éléments de partie inférieure du triangle sont laissés vides.

Exemple de fichier csv input:

5012;3257;2638;898;1734;2642;1828;599;54;172
106;4179;1111;5270;3116;1817;-103;673;535;
3410;5582;4881;2268;2594;3479;649;603;;
5655;5900;4211;5500;2159;2658;984;;;
1092;8473;6271;6333;3786;225;;;;
1513;4932;5257;1233;2917;;;;;
557;3463;6926;1368;;;;;;
1351;5596;6165;;;;;;;
3133;2262;;;;;;;;
2063;;;;;;;;;

- Le repertoire reserving_output/ est destiné à contenir les résultats des calculs.

## Manuel d'exécution

- Placer le fichier csv contenant le triangle Input dans le repertoire triangle_input/

- Modifier le fichier de code run_reserving.R pour :

		- charger l'input, choisir le nombre N d'éxécutions
		- Appeler les fonctions de calcul de provisions (inspirez vous du fichier fourni)
		- Adapter le code pour le reporting.

- Exécuter le fichier de code run_reserving.R
- Exporter vers le repertoire reserving_output/ les résultats du reporting via des fichiers csv
		


## REFERENCE

STOCHASTIC CLAIMS RESERVING IN GENERAL INSURANCE
By P. D. England and R. J. Verrall
[Presented to the Institute of Actuaries, 28 January 2002]