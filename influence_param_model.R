###########################################################################################
###########################################################################################
##                         MS210 : Analyse probabiliste des structures                   ##
##                                   Projet de fin de module                             ##
##                          Optimisation d'un hélicoptère en papier                      ##
###########################################################################################
###########################################################################################
##                            Etudiant :   Mathilde Dutreuilh                            ##
##                           Encadrant :    Guillaume Perrin                             ##
##                                                                                       ##
##                            ENSTA Paris - Année 2019-2020                              ##
###########################################################################################
###########################################################################################

#             Plan OAT pour quantifier l'influence individuelle de chaque paramètre

rm(list=ls())
graphics.off()

# Définition du modele numérique utilisé
source("modele.R")
require(lhs)

# Design de base
x_model = c(0.035, 0.055, 0.035, 0.025, 0.05, 3)
# Variations autorisées des paramètres
binf = c(0.01, 0.01, 0.01, 0.01, 0.01, 3)
bsup = c(
  0.10,
  0.29 - x_model[3] - x_model[4],
  0.29 - x_model[2] - x_model[4],
  0.29 - x_model[2] - x_model[3],
  0.21 - 2 * x_model[1],
  3)
# Legendes des axes
label=c("Rw (en m)", "Rr (en m)", "Bl (en m)", "Tl (en m)", "Tw (en m)")
couleur=c("red", "blue", "green", "orange", "cyan")

for (k in 1:5) {
  #variation de Rw : k=1
  #variation de Rr : k=2
  #variation de Bl : k=3
  #variation de Tl : k=4
  #variation de Tw : k=5
  
  # Matrice des paramètres
  mat = matrix(0, 100, 5)
  
  # Remplissage des paramètres
  for (i in 1:5) {
    tmp = seq(from = binf[k],
              to = bsup[k],
              by = (bsup[k] - binf[k]) / 99)
    if (i == k) {
      mat[, i] = tmp
    }
    else {
      mat[, i] = rep(x_model[i], 100)
    }
  } #end for i
  
# Hauteur de lâcher H=3m imposée
  mat[, 6] = rep(x_model[6], 100)
  
  temps_finaux = lacher_multi(mat, beta=1.76)
  plot(tmp, temps_finaux, type = "l", xlab = label[k], ylab ="Temps de chute (s)", col = couleur[k])
} #end for k
