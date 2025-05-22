
###############################################################################
#                         Projet : Séries Temporelles (2025)                 #
#                         Auteurs : TZIEMI NGANSOP & TANGOUO KUETE IVANA     #
#                         Date   : 12/05/2025                                #
#                         Objet  : Analyse et modélisation ARIMA             #
###############################################################################




#==============================================================================
#                             Chargement des packages                   
#==============================================================================

require(zoo) #format de serie temporelle pratique et facile d'utilisation
require(tseries)

library(seastests)
library(patchwork)
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

# Chargement et nettoyage de la base 
### Spécification du répertoire de travail
chemin <- "C:/Users/damso/Desktop/p/IPI_202501.csv"
base <- read_csv(chemin)

summary(base)
#==============================================================================#

# **************************************************************************#
#                                                                           #
#                                                                           #
#                         PART I : The Data                                 #
#                                                                           #
#                                                                           #
# **************************************************************************#
#==============================================================================#

# ========================================================================= #
# ============================= Question 1 ================================ #
# ========================================================================= #


###______-formatage des variables d'interêts______#

## Selection de nos variable d'interêt: temps et données informatique( CI)
data <- base %>% select(c(`NAF rev. 2`,`CZ`))

## Renommer les colonnes 
colnames(data)<-c("Periode","indice")

## Formatage de la variable temporelle : "Année- Mois- Jour"
data$Periode <- substr(data$Periode, 2, nchar(data$Periode))
data$Periode <-  as.Date(paste0(data$Periode, "01"), format = "%Y%m%d")


#####___Définition des données d'entrainement et de test__#

## Etant données que ARMA est un modèle à court terme  
## nous allons Garder 12 observations pour les test #

#     '''Test : Données Du 01-01-2024 au 01-01-2025'''
##     '''Entrainement: Données Avant le 01-01-2024'''


# Données d'entrainement
train_data <- data %>%
  filter(Periode < as.Date("2024-01-01"))

# Données de test
test_data <- data %>%
  filter(Periode >= as.Date("2024-01-01"))


#    ''' Copie de la train_data '''

train_copie <- train_data
#==============================================================================#


# ========================================================================= #
# =============================== Question 2 ============================== #
# ========================================================================= #


#==============================================================================#
### 2-a)Vérification de la saisonnalité de la série temporelle

## Puisque c'est une série mensuelle, nous allons vérifier la saisonnalité 
## sur 12 mois(1 an) en utilisant le test de QS(Q-statistic test)

ts_data <- ts(data$indice, frequency = 12)  
isSeasonal(ts_data, test = "qs")  # "qs" = test basé sur la variance intra-période

###########################################################
# Ce test nous indique que la série n'est pas saisonnière,
# nous allons donc nous concentrer sur la stationnarité .
#_____________________________#___________________________#


#==============================================================================#
### 2-a) Visualisation des données avant traitement

windows()

ggplot(train_data, aes(x = Periode, y = indice)) +
  geom_line(color = "black") +  # Couleur de la série
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "", y = "Valeurs de l'indice", title = "") +
  theme(panel.background = element_rect(fill = "lightgrey"))  # Couleur de fond grise



#### Verification des Autocorrélations simples et partielles de la série brute

windows()
par(mfrow=c(1,2))
acf(train_data$indice,40,main="");pacf(train_data$indice,40,main="")


#on suspecte l'existance d'une tendance déterministe nous allons donc 
# faire un test de stationarité de Dickey Fuller augmenté
#==============================================================================#



##### 2-b) Test de stationnarité de la serie

#Dans cette partie nous allons  écrire une fonction adftest_valid qui doit tester 
#  pour k allant à kmax la stationarité du AR a l'aide de adftest mais retient  
# uniquement celui dont les residus sont un bruit blanc faible.


### a- Fonction de verification de l'autocorellation des résidus 
Qtests  <- function(var, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(var, lag=l, type="Ljung-Box", 
    fitdf=fitdf)$p.value # Le test de Ljung-Box teste l'hypothèse nulle selon laquelle les retards jusqu'à lag inclusivement de l'autocorrélation sont nuls
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}


## b - Test de Dickey Fuller augmenté
# function qui cherche le lag idéal pour la régression de Dickey Fuller
adfTest_valid <- function(var, kmax, adftype){
  k <- 0
  noautocorr <- 0 # est-ce que les résidus sont un BB faible 
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
#==============================================================================#

## Application de la fonction pour trouver le lag efficace pour le test ADF
best_adf <- adfTest_valid(var=train_data$indice,kmax=24,adftype="ct") # ADF with 4 lags: residuals OK? OK

best_adf


  ### 2-c) Elimination de la stationnarité

# Nous allons différencier la série initiale avec un lag de 1 et nommer sa transformé dindice

train_data["d_indice"]=c(NA,diff(train_data$indice,1)) # différenciation d'ordre 1
train_copie["d_indice"]=c(NA,diff(train_data$indice,1)) # différenciation d'ordre 1



## Représentation des autocorrélations simples des séries indice et dindice(serie différenciée)

windows()

par(mfrow=c(1,2))
acf(train_data$indice,lag.max=24,main="",xaxt = "n")
axis(1, at = 0:24, labels = FALSE)
text(x = 0:24, y = par("usr")[3] - 0.05, labels = 0:24, srt = 90, adj = 1, xpd = TRUE)

acf(train_data$d_indice[-1],lag.max=24,main="",xaxt = "n")
axis(1, at = 0:24, labels = FALSE)
text(x = 0:24, y = par("usr")[3] - 0.05, labels = 0:24, srt = 90, adj = 1, xpd = TRUE)


## 2-d) Test de stationnarité de la série différenciée d_indice

best_dadf <- adfTest_valid(var=train_data$d_indice,kmax=24,adftype="ct") # ADF with 3 lags: residuals OK? OK

best_dadf

### Nous constatons que la série différencié est stationnaire


## Utilisation des autres tests pour confirmer les résultats précédents
# KPSS test:
kpss.test(train_data$d_indice[-1]) # stationarité
pp.test(train_data$d_indice[-1])# stationnarité

#la stationnarité et bel et bien vérifiée

# ========================================================================= #
# =============================== Question 3 ============================== #
# ========================================================================= #

###  Représenter graphiquement la série choisie avant et après transformation.----

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
#==============================================================================#

# **************************************************************************#
#                                                                           #
#                                                                           #
#                         PARTIE II : Modèles ARMA                          #
#                                                                           #
#                                                                           #
# **************************************************************************#

#==============================================================================#

# ========================================================================= #
# =============================== Question 4 ============================== #
# ========================================================================= #

###Choix du modèle ARMA(p,q) et vérification de la validité

# Pour le choix des ordres p et q, nous allons nous servir des fonctions d'autocorrélation (simples et partielles)

windows()
par(mfrow=c(1,2))
acf(train_data$d_indice[-1],20,main="",xaxt = "n")
axis(1, at = 0:20, labels = FALSE)
text(x = 0:20, y = par("usr")[3] - 0.05, 
     labels = 0:20, srt = 90, adj = 1, xpd = TRUE)
pacf(train_data$d_indice[-1],20,main="",xaxt = "n")#on regarde jusqu'à 20 mois de retard
axis(1, at = 1:20, labels = FALSE)
text(x = 1:20, y = par("usr")[3] - 0.02, 
     labels = 1:20, srt = 90, adj = 1, xpd = TRUE)

pmax=2 ; qmax=2


## Estimation de tous les modèles possibles

# Initialisation de la liste des modèles valides : ceux ayant des bruits blancs comme résidus et ayant des ordres pet q bien spécifiés

# Nous explores tous les couples (p, q) entre 0 et pmax/qmax, puis nous testons si le modèle ARIMA(p, 0, q) :
  
# a. Peut être estimé sans erreur.

# b. Produit des résidus proches du bruit blanc.

# c. A des coefficients significatifs.

# d. Et est raisonnablement parcimonieux (AR ou MA pur si possible).



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
   #--------MA(0,2), ARMA(1,2), ARMA(2,1)-----------#
# Nous allons projeter ces modèles sur les données d'entrainement
# afin de les comparer puis de choisir le meilleur modèle selon 
# les critères AIC et BIC.

#==============================================================================#
# p = 0, q = 2
ma2 <- estim <- arima(train_data$d_indice,c(0,0,2))

# p=1,q=2
ar1ma2 <- estim <- arima(train_data$d_indice,c(1,0,2))

# p=2, q=1
ar2ma1 <- estim <- arima(train_data$d_indice,c(2,0,1))
#==============================================================================#

#Choix du meilleur modèles parmis 03 modèles valides selon les critères AIC et BIC

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
# Afin de choisir un des deux modèles, nous allons effectuer
#un test de vraissemblance avec les hypothèses suivantes:
# Ho: le modèle ma2 suffit (le terme AR(1) n’apporte rien).
# H1: le modèle ar1ma2 est meilleur (le terme AR(1) améliore la vraisemblance).

#==============================================================================#

# Calcul de la statistique du test de vraisemblance (likelihood ratio test)
ll_ratio <- 2 * (scores["logLik", "ar1ma2"] - scores["logLik", "ma2"])
df_diff <- (length(coef(ar1ma2)) - length(coef(ma2))) 

# Vérification de la statistique du test et du p-value
p_value <- 1 - pchisq(ll_ratio, df_diff)
cat("P-value : ", p_value, "\n")

#  -- La p-value(0,04) est inférieure à 5% donc la
# différence entre les deux modèles est significatives.
#  -- Le modèle ARMA(1,2) permet de mieux
#   expliquer les données que le modèle  MA(2).
#==============================================================================#


# ____________________________Modèle:ARMA(1,2)_________________________________#


#==============================================================================
##---------5. Exprimer le modèle ARIMA(p,d,q) pour la série choisie------------#
#==============================================================================


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
# - But : Teste si les coefficients du modèle sont supérieure à 1,96
# - Resultat: TRUE, tous les coef sont supérieur à 1,96
# conclusion:Les composantes ar1, ma1, ma2 apporte des informations
# significatives, donc le modèle ARIMA(1,1,2) semble etre bien specifié


## Représentation des résidus ainsi que des autocorrélations associés

par(mfrow=c(1,2))

plot(arima112$residuals,main="",xlab="",ylab="",col="darkblue")
acf(coredata(arima112$residuals),lag.max = 20,main="",xaxt = "n")
# Ajouter l'axe des abscisses avec des graduations de 1 en 1 sans labels
axis(1, at = 0:20, labels = FALSE)
# Ajouter les labels inclinés à 90 degrés
# Utiliser la fonction text() pour ajouter les labels inclinés
text(x = 0:20, y = par("usr")[3] - 0.05, labels = 0:20, srt = 90, adj = 1, xpd = TRUE)


##_______________Verification de l'inversibilité du ARIMA(1,1,2)______________________#


#==============================================================================#
#  POur cela, nous allons vérifier si les racines du polynome ARMA(1,1,2) sont à l'extérieur du cercle unité
# - on commence par extracter les coefficients du modèle ARIMA
# -puis, on Selectionne des coefs ma1 et ma2
#  -ensuite on Calcule les racines du polynome

summary(arima112)
# Extraction des coefficients 
coefs <- arima112$coef
theta <- c(1, coefs["ma1"], coefs["ma2"])

# Calcul des racines du polynôme
roots <- polyroot(theta)
print(roots)

#  - on obtient les racines c(1.409179+0i, -3.131788-0i)
#  - leurs modules sont respectivements 1.409179 et 3.131788
#  - toute les racines sont à l'extérieur du cercle unité

#___Conclusion: Le Modèle ARIMA(1,1,2) est identifiable et correctement estimé____#
#___Il fournit des predictions fiable et robuste aux erreurs passées____#


#==============================================================================#

#******************************************************************************#
#*                        Partie 3: Prediction                                *#
#******************************************************************************#


#==============================================================================#
#6. Ecrire l’équation vérifiée par la région de confiance de niveau α sur les
#___valeurs futures (XT+1,XT+2).

#**************************(voir Rapport)*************************************#

#7. Préciser les hypothèses utilisées pour obtenir cette région.

#**************************(voir Rapport)*************************************#

# Convertir la série en data frame pour ggplot
data_frame <- data.frame(valeurs = as.numeric(arima112$residuals))


#==============================================================================#
#__8. Représenter graphiquement cette région pour α = 95%. Commenter.

#___---- PRÉDICTION SUR 5 MOIS ----___#
#__Prédiction avec intervalles de confiance
pred <- forecast(arima112, h = 5, level = 95)

#__Affichage des valeurs prédites et intervalles
print(pred)

# ---- EXTRACTION DES DONNÉES TEST ----
# Pour cela, nous vérifions d'abord si test_data$d_indice contient bien au moins 5 valeurs
test_vals <- test_data$indice[1:5]

#________---- TRAÇONS LE GRAPHE ----_______#


#__Séries complètes à tracer
valeurs_observées <- c(train_data$indice, rep(NA, 5))
valeurs_test <- c(rep(NA, length(train_data$indice)), test_vals)
valeurs_prédites <- c(rep(NA, length(train_data$indice)), as.numeric(pred$mean))
IC_up <- c(rep(NA, length(train_data$indice)), pred$upper)
IC_low <- c(rep(NA, length(train_data$indice)), pred$lower)

#__Création d'une séquence de dates mensuelles
date_index <- seq(as.Date("1990-01-01"), by = "month", length.out = length(train_data$d_indice) + 5)

#__Construction d'un DataFrame pour le graphique
df_plot <- data.frame(
  Date =date_index ,
  Train = c(train_data$indice, rep(NA, 5)),
  Test = c(rep(NA, length(train_data$indice)), head(test_data$indice, 5)),
  Prévisions = c(rep(NA, length(train_data$indice)), as.numeric(pred$mean)),
  IC_low = c(rep(NA, length(train_data$indice)), as.numeric(pred$lower)),
  IC_up = c(rep(NA, length(train_data$indice)), as.numeric(pred$upper))
)

#==============================================================================#

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

#==============================================================================#
#___Graphique de la série temporelle avec prévisions et intervalles de confiance


################################################
#==============================================================================#
#___Graphique de l'intervalle de confiance: ellipse
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

#==============================================================================#
#___Graphique de la série temporelle avec prévisions 

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

# affichage des deux graphiques dans une seule fenetre
p2 + p1 

#==============================================================================#
#__9. conditions pour améliorer la prévision

#**************************(voir Rapport)*************************************#

#------------------------------------------------------------------------------#
#==================================-FIN-=======================================#
#------------------------------------------------------------------------------#