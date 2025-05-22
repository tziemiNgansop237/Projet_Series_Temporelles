
###############################################################################
#         Fichier de test des fonctionnalités                                 #
#                PROJET DE SERIE TEMPORELLES                                  #
##        Auteurs : TZIEMI Raymond  & Ivana Kuete                             #    
##  Classe : 2A                                                               #
##  Date : Mai 2024                                                           #
###############################################################################

## Spécification du répertoire de travail

setwd(choose.dir())

## Chargement des librairies nécessaires  
library(dplyr)
library(readxl)
library(ellipse)
require(zoo)
library(ggplot2)
library(gridExtra)
library(fUnitRoots)
library(forecast)
library(tseries)
library(lubridate)
library(urca)
library(lubridate)

# Charger les bibliothèques nécessaires
library(readxl)    # pour lire les fichiers Excel
library(readr)     # pour écrire et lire les fichiers CSV

# Définir le chemin du fichier d'origine
chemin_xls <- "C:/Users/damso/Desktop/p/IPI_202501.xls"

# Lire le fichier Excel (on suppose qu'il y a une seule feuille ou que la bonne est la première)
donnees <- read_excel(chemin_xls)

# Définir le chemin de sortie pour le fichier CSV
chemin_csv <- "C:/Users/damso/Desktop/p/IPI_202501.csv"

# Sauvegarder les données en format CSV
write_csv(donnees, chemin_csv)

# Importer le fichier CSV
ipi_data <- read_csv(chemin_csv)

# Afficher un aperçu
head(ipi_data)



###### I- Importation du dataset  
base<-read_excel("IPI_202501.xls")

 ##### 1-formatage des variables d'interêts

    ## Selection de nos variable d'interêt: temps et données informatique( CI)
#data <- base[c(1,7)]
data <- base %>% select(c(`NAF rev. 2`,`CZ`))

    ## Renommer les colonnes 
colnames(data)<-c("Periode","indice")

    ## Formatage de la variable temporelle : "Année- Mois- Jour"
data$Periode <- substr(data$Periode, 2, nchar(data$Periode))
data$Periode <-  as.Date(paste0(data$Periode, "01"), format = "%Y%m%d")


 ##### 2- Définir les données d'entrainement et de test 

#    '''Etant données que ARMA est un modèle à court terme nous allons
#       Garder 12 observations pour les test '''

#       '''Test : Du 01-01-2024 au 01-01-2025
#          Entrainement: Avant le 01-01-2024
#       '''
        
      # Données d'entrainement
train_data <- data %>%
  filter(Periode < as.Date("2024-01-01"))
      
      # Données de test
test_data <- data %>%
  filter(Periode >= as.Date("2024-01-01"))


#    ''' Copie de la train_data '''

train_copie <- train_data




###### II-

windows()

ggplot(train_data, aes(x = Periode, y = indice)) +
  geom_line(color = "black") +  # Couleur de la série
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "", y = "Valeurs de l'indice", title = "") +
  theme(panel.background = element_rect(fill = "lightgrey"))  # Couleur de fond grise

#### 1- Autocorrélations simples et partielles de la série brute
windows()
par(mfrow=c(1,2))
acf(train_data$indice,40,main="");pacf(train_data$indice,40,main="")

#on suspecte l'existance d'une tendance déterministe : test de stationarité

### test de stationnarité ----
### ADF test

#    '''La fonction suivante effectue un test de ***Ljung-Box*** sur une série 
#     temporelle pour vérifier l'autocorrélation des résidus d'un modèle ajusté 
#    '''

#   '''H0 (hypothèse nulle) : Il n'y a pas d'autocorrélation dans les résidus. Les 
#   '''résidus sont indépendants (aucune relation entre eux).

#   '''H1 (hypothèse alternative) : Il existe une autocorrélation dans les 
#   '''résidus. Les résidus ne sont pas indépendants.

Qtests  <- function(var, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(var, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value # Le test de Ljung-Box teste l'hypothèse nulle selon laquelle les retards jusqu'à lag inclusivement de l'autocorrélation sont nuls
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}



#Resultat: Autocorellatoon entre les données



##### 2- Test de stationnarité de la serie
   
   ## a- Test de Dickey Fuller augmenté

#   ''' Nous selectionons tout d'abord automatiquement le nombre de lags optimal 
#       par le crictère AIC, ensuite nous verifions si les résidus de ce modèle 
#       sont indépendants , afin de renforcer la robustesse, au cas contraire  
#       les résidedus montre l'autocorellation et Nous allons réajuster le modèle
#       et augmenter les lags à k+1'''





# Fonction combinée ADF avec vérification des résidus
adfTest_valid_combined <- function(var, kmax, adftype) {
  # Sélectionner les lags via BIC
  adf_test <- ur.df(var, type = adftype, selectlags = "BIC")  # Sélection des lags selon AIC
  k <- adf_test@lags  # Nombre de lags sélectionné par AIC
  print(paste("Number of lags selected by BIC: ", k))
  adf_test_aic <- ur.df(var, type = adftype, selectlags = "BIC")
  print(paste("Number of lags selected by AIC: ", adf_test_aic@lags))
  
  # Extraire les résidus du modèle ADF
  residuals_adf <- adf_test@res
  
  # Tester les résidus du modèle ajusté pour l'autocorrélation
  pvals <- Qtests(residuals_adf, 24, fitdf = length(adf_test@teststat))[, 2]
  
  if (sum(pvals < 0.05, na.rm = TRUE) == 0) {
    # Si les résidus sont indépendants (pas d'autocorrélation)
    print("Residuals are independent (no autocorrelation detected).")
  } else {
    # Si les résidus montrent encore de l'autocorrélation, ajuster le modèle
    print("Residuals show autocorrelation, adjusting the model...")
    
    # Augmenter les lags et réajuster le modèle 
    
    pvals_new <- pvals
    while (sum(pvals_new < 0.05, na.rm = TRUE) != 0){
      k <- k + 1
      print(paste("Trying with ", k, " lags."))
    
      adf_test_new <- ur.df(var, type = adftype, lags = k)
      residuals_adf_new <- adf_test_new@res
    
      # Vérifier les résidus du modèle ajusté
      pvals_new <- Qtests(residuals_adf_new, 24, fitdf = length(adf_test_new@teststat))[, 2]
    
      if (sum(pvals_new < 0.05, na.rm = TRUE) == 0) {
        print("Adjusted model: Residuals are now independent.")
      }else {
      print("Residuals still show autocorrelation even after adjustment.")
      adf_test <- adf_test_new
      }
    }
  }
  
  return(adf_test)  # Retourner le test ADF final
}





# Exemple d'application sur une série de données
adfTest_valid_combined(train_data$indice, kmax = 24, adftype = "drift")

### Deuxieme methode

# function qui cherche le lag idéal pour la régression de Dickey Fuller
adfTest_valid <- function(var, kmax, adftype){
  k <- 0
  noautocorr <- 0 # est-ce que les résidus sont un BB faible ?
  while (noautocorr==0){
    cat(paste0("ADF with ",k," lags: residuals OK? "))
    test_adf <- adfTest(var, lags=k, type=adftype)
    pvals <- Qtests(test_adf@test$lm$residuals, 20, fitdf = length(test_adf@test$lm$coefficients))[,2] # On fait le test de potmanteau
    if (sum(pvals<0.05,na.rm=T)==0) { # si aucune p_value n'entre dans la zone de rejet
      noautocorr <- 1; cat("OK \n")
    } else cat("nope \n")
    
    
    k <- k+1
  }
  return(test_adf)
}


## Application de la fonction pour trouver le lag efficace pour le test ADF
best_adf <- adfTest_valid(var=train_data$indice,kmax=24,adftype="ct") # ADF with 4 lags: residuals OK? OK

best_adf


# On différencie la série initiale

train_data["d_indice"]=c(NA,diff(train_data$indice,1)) # différenciation d'ordre 1
train_copie["d_indice"]=c(NA,diff(train_data$indice,1)) # différenciation d'ordre 1


## Représentation des autocorrélations simples des séries indice et d_indice

window()
par(mfrow=c(1,2))
acf(train_data$indice,lag.max=24,main="",xaxt = "n")
axis(1, at = 0:24, labels = FALSE)
text(x = 0:24, y = par("usr")[3] - 0.05, labels = 0:24, srt = 90, adj = 1, xpd = TRUE)

acf(train_data$d_indice[-1],lag.max=24,main="",xaxt = "n")
axis(1, at = 0:24, labels = FALSE)
text(x = 0:24, y = par("usr")[3] - 0.05, labels = 0:24, srt = 90, adj = 1, xpd = TRUE)

## Recherche du lag idéal pour tester la stationnarité de la série d_indice

best_dadf <- adfTest_valid(var=train_data$d_indice,kmax=24,adftype="ct") # ADF with 3 lags: residuals OK? OK

best_dadf

## Utilisation des autres tests pour conforter les résultats
# KPSS test:
kpss.test(train_data$d_indice[-1]) # stationarité
pp.test(train_data$d_indice[-1])# stationnarité


### 3. Représenter graphiquement la série choisie avant et après transformation.----

windows()
# Création des deux graphiques
serie_brute <- ggplot(data = train_copie[-1,], aes(x = Periode, y = indice)) +
  geom_line() +
  labs(x = "", y = "Indice", title = "")+
  theme(panel.background = element_rect(fill = "lightgrey"))

serie_diff <- ggplot(data = train_copie[-1,], aes(x = Periode, y = d_indice)) +
  geom_line() +
  labs(x = "", y = "Indice différencié", title = "")+
  theme(panel.background = element_rect(fill = "lightgrey"))

# Affichage des deux graphiques côte à côte
grid.arrange(serie_brute, serie_diff, nrow = 1)


### 4. Choisir, en le justifiant, un modèle ARMA(p,q) pour votre série corrigée Xt. Estimer les paramètres du modèle et vérifier sa validité.

# Pour le choix des ordres p et q, nous allons nous servir des fonctions d'autocorrélation (simples et partielles)

windows()
par(mfrow = c(1, 2))

# ACF
acf(train_data$d_indice[-1], 20, main = "", xaxt = "n")
axis(1, at = 0:20, labels = FALSE)
text(x = 0:20, y = par("usr")[3] - 0.05, 
     labels = 0:20, srt = 90, adj = 1, xpd = TRUE, cex = 0.6)

# PACF
pacf(train_data$d_indice[-1], 20, main = "", xaxt = "n")
axis(1, at = 1:20, labels = FALSE)
text(x = 1:20, y = par("usr")[3] - 0.02, 
     labels = 1:20, srt = 90, adj = 1, xpd = TRUE, cex = 0.6)

pmax=2 ; qmax=2


## Estimation de tous les modèles possibles

# Initialisation de la liste des modèles valides : ceux ayant des bruits blancs comme résidus et ayant des ordres pet q bien spécifiés

modeles_valides <- list()
nb <-0
for (p in 0:pmax) {
  for (q in 0:qmax) {
    modele <- try(arima(train_data$d_indice[-1], order = c(p, 0, q), 
                        include.mean = FALSE)) # Tentative d'estimation de l'ARIMA
    if ((sum(is.na(sqrt(diag(modele$var.coef))))==0) & (all(Qtests(
              modele$residuals, 24, fitdf = p+q)[,2] > 0.05, na.rm = TRUE)
               )) { # on teste si on a de bons résidus et si les écarts types sont bien estimés
      validite <- abs(modele$coef/sqrt(diag(modele$var.coef)))>1.96 # test de validité du modèle (significativité des coefficients)
      ajout <- FALSE
      if ((p==0 | q==0) & (validite[length(validite)])) {ajout <- TRUE}
      else if ((validite[p]) & (validite[p+q])) {ajout <- TRUE}
      if (ajout){
        nb <- nb + 1
        modeles_valides <- append(modeles_valides, list(list(paste("modele_", nb, sep = ""), p, q)))
      }
    }
  }
}


modeles_valides


# Les modèles ayant réussi à ces test sont les suivants:
   #--------MA(2), ARMA(1,2), ARMA(2,1)-----------#

#==============================================================================#
# p = 0, q = 2
ma2 <- estim <- arima(train_data$d_indice,c(0,0,2))

# p=1,q=2
ar1ma2 <- estim <- arima(train_data$d_indice,c(1,0,2))

# p=2, q=1
ar2ma1 <- estim <- arima(train_data$d_indice,c(2,0,1))
#==============================================================================#

#Choix du meilleur modèles parmis 03 modèles valides par les critères AIC et BIC

models <- list(ma2 = ma2, ar1ma2 = ar1ma2, ar2ma1 = ar2ma1)
scores <- sapply(models, function(mod) {
  c(AIC = AIC(mod), BIC = BIC(mod), logLik = logLik(mod))})
print(scores)

# Trouver les colonnes avec les valeurs minimales
min_aic_model <- colnames(scores)[which.min(scores["AIC", ])]
min_bic_model <- colnames(scores)[which.min(scores["BIC", ])]

# Affichage propre
cat("meilleur modèle selon AIC :", min_aic_model, "\n")
cat("meilleur modèle selon BIC :", min_bic_model, "\n")



#-----_______Meilleur modèles:_______-------#  
#          __ AIC: ARMA(1,2)  __
#           _   BIC: MA(2)    _


# Test de vraisemblance : MA(2) VS ARMA(1,2)

#==============================================================================#

# Calcul de la statistique du test de vraisemblance (likelihood ratio test)
ll_ratio <- 2 * (scores["logLik", "ar1ma2"] - scores["logLik", "ma2"])
df_diff <- (length(coef(ar1ma2)) - length(coef(ma2))) 

# Vérification de la statistique du test et du p-value
p_value <- 1 - pchisq(ll_ratio, df_diff)
cat("P-value : ", p_value, "\n")

#  - La p-value  est significative à 5%, 
#  - Le modèle ARMA(1,2) permet de mieux expliquer les données que le modèle  MA(2).
#==============================================================================#


# ____________________________Modèle:ARMA(1,2)_________________________________#
ar1ma2

##---------5. Exprimer le modèle ARIMA(p,d,q) pour la série choisie------------#


# _____________________Modèle de base :ARIMA(1,1,2)_____________________________#
arima112 <- arima(train_data$indice,c(1,1,2),include.mean=F)


#Test de validité des résidus 
all(Qtests(arima112$residuals, 24, fitdf=3)[,2]>0.05,na.rm = TRUE)
# - Resultat: TRUE ie toute les p-values sont inférieures à 0.05
# - conclusion:les résidus ne montrent pas d'autocorrélation significative p
#   pour tous les lags testés, indiquant que le modèle ajuste bien les données. 


# Significativité des coefficients
arima112$coef/sqrt(diag(arima112$var.coef))
abs(arima112$coef/sqrt(diag(arima112$var.coef)))>1.96
# - But : Teste si les coefficients du modèle sont supérieure à 1,86
# - Resultat: TRUE, tous les coef sont supérieur à 1,96
# - conclusion:Les composantes ar1, ma1, ma2 apporte des informations 
#          significatives, donc ARIMA(1,1,2) semble bien specifié



## Représentation des résidus ainsi que des autocorrélations associés

par(mfrow=c(1,2))

plot(arima112$residuals,main="",xlab="",ylab="",col="darkblue")
acf(coredata(arima112$residuals),lag.max = 20,main="",xaxt = "n")
# Ajouter l'axe des abscisses avec des graduations de 1 en 1 sans labels
axis(1, at = 0:20, labels = FALSE)
# Ajouter les labels inclinés à 90 degrés
# Utiliser la fonction text() pour ajouter les labels inclinés
text(x = 0:20, y = par("usr")[3] - 0.05, labels = 0:20, srt = 90, adj = 1, xpd = TRUE)


##_______________Verifiecation de l'inversibilité du ARIMA(1,1,2)______________________#

#==============================================================================#
#  - Selection des coefs ma1 et ma2
#  - Calcul ds racines du polynome  c(1, -0.3903 ,-0.2266)
#  - On inverses les racines c(1.409179+0i, -3.131788-0i)
#  - On obtient les modules  c(1.409179,3.131788)
#  - toute les racines sont à l'xtérieur du cercle unité

summary(arima112)
# Extraction des coefficients 
coefs <- arima112$coef
theta <- c(1, coefs["ma1"], coefs["ma2"])

# Calcul des racines du polynôme
roots <- polyroot(theta)
print(roots)

# Calcul des racines du polynôme
roots <- polyroot(theta)
print(roots)


#___Conclusion:Modèle identifiable, correctement estimé,
#              fournit des predictions fiable et robuste aux erreurs passées____
#==============================================================================#

#******************************************************************************#
#*                        Partie 3: Prediction                                *#
#******************************************************************************#



#6. Ecrire l’équation vérifiée par la région de confiance de niveau α sur les
#     valeurs futures (XT+1,XT+2).

#**************************(voir document)*************************************#

#7. Préciser les hypothèses utilisées pour obtenir cette région.

#**************************(voir document)*************************************#

# Convertir la série en data frame pour ggplot
data_frame <- data.frame(valeurs = as.numeric(arima112$residuals))


#8. Représenter graphiquement cette région pour α = 95%. Commenter.

# ---- PRÉDICTION SUR 5 MOIS ----
# Prédiction avec intervalles de confiance
pred <- forecast(arima112, h = 5, level = 95)

# Affichage des valeurs prédites et intervalles
print(pred)

# ---- EXTRACTION DES DONNÉES TEST ----
# Assure-toi que test_data$d_indice contient au moins 5 valeurs
test_vals <- test_data$indice[1:5]

# ---- TRAÇONS LE GRAPHE ----


# Séries complètes à tracer
valeurs_observées <- c(train_data$indice, rep(NA, 5))
valeurs_test <- c(rep(NA, length(train_data$indice)), test_vals)
valeurs_prédites <- c(rep(NA, length(train_data$indice)), as.numeric(pred$mean))
IC_up <- c(rep(NA, length(train_data$indice)), pred$upper)
IC_low <- c(rep(NA, length(train_data$indice)), pred$lower)

# Créer une séquence de dates mensuelles
date_index <- seq(as.Date("1990-01-01"), by = "month", length.out = length(train_data$d_indice) + 5)

# Construire le data.frame pour le graphique
df_plot <- data.frame(
  Date =date_index ,
  Train = c(train_data$indice, rep(NA, 5)),
  Test = c(rep(NA, length(train_data$indice)), head(test_data$indice, 5)),
  Prévisions = c(rep(NA, length(train_data$indice)), as.numeric(pred$mean)),
  IC_low = c(rep(NA, length(train_data$indice)), as.numeric(pred$lower)),
  IC_up = c(rep(NA, length(train_data$indice)), as.numeric(pred$upper))
)





ggplot(df_plot, aes(x = Date)) +
  geom_line(aes(y = Train), color = "black", size = 0.4) +
  geom_line(aes(y = Test), color = "blue", size = 0.5, linetype = "dashed") +
  geom_line(aes(y = Prévisions), color = "red", size = 0.6) +
  geom_ribbon(aes(ymin = IC_low, ymax = IC_up), fill = "gray80", alpha = 0.5) +
  annotate("rect",
           xmin = df_plot$Date[length(train_data$indice) + 1],
           xmax = max(df_plot$Date),
           ymin = -Inf, ymax = Inf,
           fill = "gray90", alpha = 0.4) +
  labs(x = "Date", y = "Indice de production manufacturière") +
  theme_minimal(base_family = "serif") +
  theme(
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.position = "none"  # retirer la légende pour un style académique propre
  )


library(ellipse)
pred <- forecast(arima112, h = 2, level = 95)

# Suppose une corrélation rho entre les deux prévisions
rho <- 0.8  # corrélation supposée
sigma1 <- (pred$upper[1] - pred$lower[1]) / (2 * 1.96)
sigma2 <- (pred$upper[2] - pred$lower[2]) / (2 * 1.96)

# Calcul de la covariance à partir des écarts-types et de la corrélation
cov_val <- rho * sigma1 * sigma2

# Matrice de covariance avec corrélation
Sigma <- matrix(c(sigma1^2, cov_val,
                  cov_val, sigma2^2), nrow = 2)

# Centre des prévisions
mu <- c(as.numeric(pred$mean[1]), as.numeric(pred$mean[2]))


ellipse_pts <- ellipse(Sigma, centre = mu, level = 0.95)
ellipse_df <- as.data.frame(ellipse_pts)





################################################
# Charger les bibliothèques nécessaires
library(ggplot2)
library(patchwork)

# Ton premier graphique : ellipse
p1 <- ggplot(ellipse_df, aes(x = x, y = y)) +
  geom_path(color = "black") +
  geom_point(aes(x = mu[1], y = mu[2]), color = "black", size = 2) +
  labs(x = "Prévision Mars 2022", y = "Prévision Avril 2022",
       title = "Ellipse de confiance jointe ") +
  theme_minimal(base_family = "serif") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )



# Ton deuxième graphique : série temporelle
p2 <- ggplot(df_plot, aes(x = Date)) +
  geom_line(aes(y = Train), color = "black", size = 0.4) +
  geom_line(aes(y = Test), color = "blue", size = 0.5, linetype = "dashed") +
  geom_line(aes(y = Prévisions), color = "red", size = 0.6) +
  geom_ribbon(aes(ymin = IC_low, ymax = IC_up), fill = "gray80", alpha = 0.5) +
  annotate("rect",
           xmin = df_plot$Date[length(train_data$indice) + 1],
           xmax = max(df_plot$Date),
           ymin = -Inf, ymax = Inf,
           fill = "gray90", alpha = 0.4) +
  labs(x = "Date", y = "Indice de production manufacturière",
       title = "Prévisions  de l'indice et intervalles de confiance") +
  theme_minimal(base_family = "serif") +
  theme(
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.position = "none"
  )

# Combinaison côte à côte
p2 + p1  # ou p1 + p2


