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

#                Code de calage du coefficient de frottement quadratique beta

rm(list=ls())
graphics.off()

# Définition du modele numérique utilisé
source("modele.R")
require(lhs)

#############################################################################################
#                     Information disponible à partir de l'expérience 

data = read.csv("experience.csv", header = TRUE, sep =";", dec=",")
yobs=data[,15] #contient la 15e colonne du fichier de données --> temps de chute pour chaque lancer
Nobs = length(yobs)
xobs = matrix(0, Nobs, 6)
xobs[,1] = data[,3]*1e-2
xobs[,2] = data[,4]*1e-2
xobs[,3] = data[,5]*1e-2
xobs[,4] = data[,6]*1e-2
xobs[,5] = data[,7]*1e-2
xobs[,6] = data[,2]

# yobs=data[,9]
# Nobs = length(yobs)
# xobs = c(0.035, 0.055, 0.035, 0.025, 0.05,3)
# xobs = matrix(rep(xobs, each=Nobs), ncol =6)

#ajout des incertitudes de mesure
sig_mes = 0.02 #en secondes
bruit_mes = matrix(rnorm(Nobs)*sig_mes, ncol = 1)
yobs = yobs+bruit_mes

#ajout des incertitudes dites constructeur
sig_dim = 0.001 #en metres
bruit_dim = matrix(rnorm(5*Nobs)*sig_dim, ncol = 5)
sig_H = 0.01 #en metres
bruit_H = matrix(rnorm(Nobs)*sig_H, ncol = 1)

xobs = xobs+cbind(bruit_dim, bruit_H)


################################################################################################
##                                 Calage bayésien de la valeur de Cx
## Via estimateur du maximum de vraisemblance

len_seq = 500
beta_test = seq(beta_MIN,beta_MAX,length = len_seq)
logvrais = matrix(0,1,len_seq)
for(i in 1:len_seq){
  logvrais[i] = LogVraisemblance(10, xobs, yobs, sig_mes, beta_test[i])
}
maxi = which(logvrais == max(logvrais))
logvrais = logvrais-max(logvrais)
PDF = exp(logvrais)

# representation graphique de la loi a posteriori
plot(beta_test, PDF, type = "l")

tmp1 = beta_test[which(PDF==1)]
beta_simul = sum(tmp1)/length(tmp1) ## moyenne des valeurs pour lesquelles on a vraisemblance = 1

abline(v = beta_simul, col = 'cyan')

# définition arbitraire d'un sigma pour beta
tmp2 = beta_test[which(0.5<=PDF)]
plus = tmp2[length(tmp2)]
moins = tmp2[1]
sig_beta = ((plus - beta_simul)  + (beta_simul - moins)) / 2 

# approximation normale pour les valeurs admissibles de beta
# c'est celle qu'on garde pour la suite
approx_normale <- rnorm(1000, beta_simul, sig_beta)
plot(function(approx_normale) 0.01*dnorm(approx_normale, beta_simul, sig_beta), xlim = c(beta_MIN, beta_MAX), col = "red", add= TRUE)

# beta_simul = 1.78466933867735
# sig_beta = 0.00529058116232473

###############################################################################################
##                           Tracé de position en fonction du temps
##                          Comparaison des beta moyenne et "extremes"

lacher_mean = lacher(xobs, beta_simul, bool.save = T)
lacher_sup = lacher(xobs, beta_simul + sig_beta, bool.save = T)
lacher_inf = lacher(xobs, beta_simul - sig_beta, bool.save = T)

z_mean = lacher_mean$Zsave[1,]
t_mean = seq(0, lacher_mean$temps_final, length.out = length(z_mean))
z_sup = lacher_sup$Zsave[1,]
t_sup = seq(0, lacher_sup$temps_final, length.out = length(z_sup))
z_inf = lacher_inf$Zsave[1,]
t_inf = seq(0, lacher_inf$temps_final, length.out = length(z_inf))

plot(t_mean, z_mean, xlab = "time", ylab = "altitude", type="l", col = "red")
lines(t_sup, z_sup, type="l", col = "blue")
lines(t_inf, z_inf, type="l", col = "cyan")
