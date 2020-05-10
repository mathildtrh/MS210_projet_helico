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

#                    Code d'optimisation de la forme de l'hélicoptère
##                            Maximisation du temps de chute

rm(list=ls())
graphics.off()

# Définition du modele numérique utilisé
source("modele.R")
require(lhs)
require(pracma)

# Récupération des valeurs de beta_simul et de sig_beta voulues
load("resultats_calibration.RData") 


##########################################################################################
#                 Représentation graphique de notre fonction à maximiser
#                                       Par iso-niveaux

# Paramètres pré-fixés
beta = beta_simul
H = 3
N.grid=50 #taille grille isoniveaux

# Paramètres fixés artificiellement pour la première étude
Bl = 0.02
Tl = 0.02
Tw = 0.02

###########################################################################################
##                                      Fonction isoniv
## Entrées : N.grid = nombre de points d'évaluation pour chaque axe de la grille
##           Bl, Tw, Tl, H = paramètres fixés pour le tracé des iso-niveaux Rw-Rr
## Sorties : y.grid = évaluation de la fonction lacher sur les points de la grille
##           Rw.grid = axe des abscisses
##           Rr.grid = axe des ordonnées
##           tracés des iso-niveaux

isoniv <- function (N.grid, Bl, Tl, Tw, H) {
  Rw.grid = seq(X_MIN[1], X_MAX[1], len = N.grid)
  Rr.grid = seq(X_MIN[2], X_MAX[2], len = N.grid)
  par.grid <- meshgrid(x = Rw.grid,
                       y = Rr.grid)
  
  y.grid = matrix(0, N.grid, N.grid)
  i1 = 0
  for (Rw in Rw.grid) {
    i1 = i1 + 1
    i2 = 0
    #print(i1)
    for (Rr in Rr.grid) {
      i2 = i2 + 1
      y.grid[i1, i2] = lacher(x = c(Rw, Rr, Bl, Tl, Tw, H),
                              beta,
                              temps.seul = T)
    } #end for Rr
  } #end for Rw
  
  # Representation graphique
  contour(
    x = Rw.grid,
    y = Rr.grid,
    z = y.grid,
    nlevels = 20,
    xlab = "Rw",
    ylab = "Rr"
  )
  return(list(y.grid = y.grid, Rw.grid = Rw.grid, Rr.grid = Rr.grid))
} #end function isoniv

##################################################################################################
##                                  Fonction traj_optimale
## Entrées : N.grid = nombre de points d'évaluation pour chaque axe de la grille
##           y.grid = évalution des sur chaque point de la grille
##           Bl, Tw, Tl, H = paramètres fixés pour le calcul de y.grid
## Sorties : Rw.optim, Rr.optim = valeurs de Rw de Rr qui maximisent tf pour ce jeu de paramètres
##           tf = temps de chute maximal atteint
##           tracé de la trajectoire optimale

traj_optimale <- function(N.grid, y.grid, Rw.grid, Rr.grid, Bl, Tl, Tw, H) {
  n.optim = which.max(y.grid)
  n.Rw = n.optim %% N.grid + N.grid ## attention au modulo, index commence à 1
  if (n.Rw == 0) {
    n.Rw = Ngrid
  }
  n.Rr = floor(n.optim / N.grid) + 1 * (floor(n.optim / N.grid) < n.optim / N.grid)
  Rw.optim = Rw.grid[n.Rw]
  Rr.optim = Rr.grid[n.Rr]
  
  ##liste des données pour une lâcher optimal
  traj.optim = lacher(
    x = c(Rw.optim, Rr.optim,  Bl, Tl, Tw, H),
    beta,
    bool.save = T
  )
  pos = traj.optim$Zsave[1, ]
  tf = traj.optim$temps_final
  time = seq(0, tf + 0.001, by = 0.001)
  
  #tracé de la trajectoire
  plot(time, pos, type = "l")
  points(tf, 0, pch = 19, col = "red")
  
  return(list(Rw.optim = Rw.optim, Rr.optim = Rr.optim, tf=tf))
} #end function traj_optimale


############################################################################################
##                Première estimation de Rw et Rr otpimaux à Bl, Tl et Tw fixés

eval = isoniv(N.grid, Bl, Tl, Tw, H)
pales_opti_1 = traj_optimale(N.grid, eval$y.grid, eval$Rw.grid, eval$Rr.grid, Bl, Tl, Tw, H)

#############################################################################################
## Amélioration possible : itérer sur Bl, Tl, Tw (même si on se doute du résultat)
## Choisir beta avec une forme gaussienne, au lieu d'une unique valeur

############################################################################################
##              Comparaison des résultats : simulation vs experience

# EXPERIENCE

## Données expérimentales avec hélicoptère optimal
data_opti = read.csv("test_helico_opti.csv", header = TRUE, sep =";", dec=",")
y_obs=data_opti[,15] #contient la 15e colonne du fichier de données --> temps de chute pour chaque lancer
N_obs = length(yobs)

## Histogramme de la répartition des y_obs avec l'hélico optimal
nn = 11
hist(y_obs,nclass = nn)
abline(v = mean(y_obs), col = "magenta")

# SIMULATION
## Tirage gaussien de beta selon sa répartition + tracé
beta_distrib <- rnorm(1000, beta_simul, sig_beta)
## Simulations pour chaque beta tiré
y_simul = matrix(0,1,1000)
for (k in 1:1000){
  y_simul[k]=lacher(x = c(pales_opti_1$Rw.optim, pales_opti_1$Rr.optim, Bl, Tl, Tw, H),beta_distrib[k], temps.seul=T)
} #end for k

## tracer la répartition des y_simul
plot(density(y_simul),col = "red")