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

#                   Modèle numérique de chute de l'hélicoptère en papier

##    - lacher
##    - lacher_multi

#                                       Aides au calcul
##    - LogVraisemblance
##    - stats_resume

###########################################################################################
##                                  Fonction lacher
##  Entrées : x = dimensions de l'hélico (Rw, Rr, Bl, Tl, Tw)
##            beta = coeffcient de frottement quadratique <-> Cx en aérodynamique
##            temps.seul = booléen if T alors ne renvoie que le temps de chute final
##            dt = pas de temps de calcul
##            bool.save = booléen if T alors on sauvegarde la trajectoire de l'hélico
##  Sorties : if temps.seul = T --> TF temps de chute
##            if temps.seul = F --> liste de taille 3 contenat Zsave, DZsave, TF

lacher <-
  function(x,
           beta,
           temps.seul = F,
           dt = 0.001,
           bool.save = F) {
    # x=(Rw, Rr, Bl, Tl, Tw) --> on en déduit la masse et le maitre couple de l'hélicoptère
    # beta est le Cx de l'hélicoptère, supposé indépendant de sa taille
    
    # description des parametres systeme
    Rw = x[1]         #largeur des rotors
    Rr = x[2]         #longueur des rotors
    Bl = x[3]         #longueur du corps
    Tl = x[4]         #longueur de la queue
    Tw = x[5]         #largeur de la queue
    
    H = x[6]          #hauteur de lancer, en metres
    
    m = (80 * (2 * Rw * (Rr + Bl) + Tw * Tl) + 25 / 35) * 1e-3    #masse de l'hélicoptère avec masse du trombone = 25/35 grammes
    #et grammage du papier = 80g/m² --> converti en kg
    S = 2 * Rw * Rr     #maitre couple
    
    # initialisation des grandeurs calculees en 1D
    # DZ[1] vitesse DZ[2] accélération
    # Z[1] position Z[2] vitesse
    # t temps écoulé pendant la chute
    
    DZ = matrix(0, 2, 1)
    Z = matrix(c(H, 0), 2, 1)
    t = 0
    
    # sauvegarde des trajectoires
    Zsave = Z
    DZsave = DZ
    Tsave = t
    
    # boucle de calcul par schema explicite d'euler en 1D
    while (Z[1] >= 0) {
      DZ = c(Z[2], 0.5 * 1.3 * S * beta / m * Z[2] * Z[2] - 9.81)
      Z = DZ * dt + Z
      t = t + dt
      # sauvegarde des trajectoires
      if (bool.save) {
        Zsave = cbind(Zsave, Z)
        DZsave = cbind(DZsave, DZ)
        Tsave = cbind(Tsave, t)
      }#end if
    }#end while
    TF = t #temps final, lorsque Z=0, sol atteint
    
    if (temps.seul) {
      return = TF
    }
    else {
      return = list(Zsave = Zsave,        # trajectoire
                    DZsave = DZsave,      # evolution temporelle de la vitesse
                    temps_final = TF)     # temps de chute
    }#end else
  }#end function

###########################################################################################
##                              Fonction lacher_multi
##  Entrées : x = dimensions de l'hélico (Rw, Rr, Bl, Tl, Tw) 1 ligne = 1 design
##            beta = coeffcient de frottement quadratique <-> Cx en aérodynamique
##            dt = pas de temps de calcul
##  Sorties : TF temps de chute


lacher_multi <- function(x, beta, dt = 0.001) {
  # x=(Rw, Rr, Bl, Tl, Tw) --> on en déduit la masse et le maitre couple de l'hélicoptère
  # beta est le Cx de l'hélicoptère, supposé indépendant de sa taille
  
  #x peut avoir plusieurs lignes!!
  param = matrix(x, ncol = 6)
  #on trie les paramètres en 1 colonne pour chaque paramètre
  TF = rep(0, nrow(param))
  
  for (i in 1:nrow(param)) {
    # description des parametres systeme
    Rw = param[i, 1]         #largeur des rotors
    Rr = param[i, 2]         #longueur des rotors
    Bl = param[i, 3]         #longueur du corps
    Tl = param[i, 4]         #longueur de la queue
    Tw = param[i, 5]         #largeur de la queue
    
    H = param[i, 6]                 #hauteur de lancer, en metres
    
    m = (80 * (2 * Rw * (Rr + Bl) + Tw * Tl) + 25 / 35) * 1e-3    #masse de l'hélicoptère
    # avec masse du trombone = 25/35 grammes
    #et grammage du papier = 80g/m² --> converti en kg
    S = 2 * Rw * Rr     #maitre couple
    
    # initialisation des grandeurs calculees en 1D
    # DZ[1] vitesse DZ[2] accélération
    # Z[1] position Z[2] vitesse
    # t temps écoulé pendant la chute
    
    DZ = matrix(0, 2, 1)
    Z = matrix(c(H, 0), 2, 1)
    t = 0
    
    # boucle de calcul par schema explicite d'euler en 1D
    while (Z[1] >= 0) {
      DZ = c(Z[2], 0.5 * 1.3 * S * beta / m * Z[2] * Z[2] - 9.81)
      Z = DZ * dt + Z
      t = t + dt
    }
    
    TF[i] = t #temps final, lorsque Z=0, sol atteint
  }#end for
  return(RES = TF)
}#end fonction


###########################################################################################
##                                    Définition des bornes

# pour les parametres systeme
# ATTENTION il faudra vérifier de la faisabilité de l'hélico avec les mesures prises
X_MIN = c(0.01, 0.01, 0.01, 0.01, 0.01, 2)
X_MAX = c(0.05, 0.18, 0.18, 0.18, 0.10, 3)
nomsx <- c("Rw", "Rr", "Bl", "Tl", "Tw", "H")

# pour le parametre de calibration
beta_MIN = 1.77
beta_MAX = 1.80
nomsbeta <- "beta"


###########################################################################################
##                              Fonction LogVraisemblance
##  Entrées : beta = coeffcient de frottement quadratique <-> Cx en aérodynamique
##            xobs = dimensions de l'hélico utilisé pour cette expérience
##            yobs = temps de chute mesuré pour cette expérience
##  Sorties : result = logvraisemblance d'un coeff beta par rapport aux données pde l'expérience

LogVraisemblance <- function(Nobs, xobs, yobs, sig_mes, beta){
  ytest = matrix(0,  nrow=Nobs, ncol=1)
  ytest = lacher_multi(xobs, beta)
  
  result = -Nobs/2*log(sig_mes^2)-1/2*sum ((ytest - yobs )*(ytest - yobs)) /(sig_mes^2)
  return(result)
}

############################################################################################
##                                Fonction stats_resume
##  Entrées : yobs = vecteur de valeurs homogènes
##  Sorties : graphiques de résumé des stats

stats_resume <- function(yobs){
  resume = c(median(yobs), mean(yobs), min(yobs), max(yobs), quantile(yobs,probs = c(0.25,0.75)))
  boxplot(yobs)
  
  nn = 15
  hist(yobs,nclass = nn)
  abline(v = resume[], col = c('blue', 'cyan', 'yellow', 'green', 'magenta', 'red'))
  ## pas encore réussi à afficher la légende
  legend(x = 0, legend = c("médiane", "moyenne", "mini", "maxi", "quantile 0.25", "quantile 0.75"))
  
  plot(density(yobs))
  abline(v = resume, col = c('blue', 'cyan', 'yellow', 'green', 'magenta', 'red'))
}







