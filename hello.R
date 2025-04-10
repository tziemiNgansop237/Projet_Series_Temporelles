###############################################################################
#                                                                             #
#                PROJET DE SERIE TEMPORELLES                                  #
##  Auteurs : TZIEMI Raymond  & Ivana Kuete                                  
##  Classe : 2A
##  Date : Mai 2024                                                                            #
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


###### I- Importation du dataset  
base<-read_excel("IPI_202501.xls")

 ##### 1-formatage des variables d'interêts

    ## Selection de nos variable d'interêt: temps et données informatique( CI)
data <- train_data[c(1,7)]

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

# Test
Qtests (train_data$indice, 20,0) 
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
  pvals <- Qtests(residuals_adf, 20, fitdf = length(adf_test@teststat))[, 2]
  
  if (sum(pvals < 0.05, na.rm = TRUE) == 0) {
    # Si les résidus sont indépendants (pas d'autocorrélation)
    print("Residuals are independent (no autocorrelation detected).")
  } else {
    # Si les résidus montrent encore de l'autocorrélation, ajuster le modèle
    print("Residuals show autocorrelation, adjusting the model...")
    
    # Augmenter les lags et réajuster le modèle 
    
    pvals_new <- pvals
    while (sum(pvals_new < 0.05, na.rm = TRUE) == 0){
      k <- k + 1
      print(paste("Trying with ", k, " lags."))
    
      adf_test_new <- ur.df(var, type = adftype, lags = k)
      residuals_adf_new <- adf_test_new@res
    
      # Vérifier les résidus du modèle ajusté
      pvals_new <- Qtests(residuals_adf_new, 20, fitdf = length(adf_test_new@teststat))[, 2]
    
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
adfTest_valid_combined(train_data$indice, kmax = 20, adftype = "drift")

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

train_data["dindice"]=c(NA,diff(train_data$indice,1)) # différenciation d'ordre 1
train_copie["dindice"]=c(NA,diff(train_data$indice,1)) # différenciation d'ordre 1


## Représentation des autocorrélations simples des séries indice et dindice

windows()
par(mfrow=c(1,2))
acf(train_data$indice,lag.max=24,main="",xaxt = "n")
axis(1, at = 0:24, labels = FALSE)
text(x = 0:24, y = par("usr")[3] - 0.05, labels = 0:24, srt = 90, adj = 1, xpd = TRUE)
acf(train_data$dindice[-1],lag.max=24,main="",xaxt = "n")
axis(1, at = 0:24, labels = FALSE)
text(x = 0:24, y = par("usr")[3] - 0.05, labels = 0:24, srt = 90, adj = 1, xpd = TRUE)

## Recherche du lag idéal pour tester la stationnarité de la série dindice

best_dadf <- adfTest_valid(var=train_data$dindice,kmax=24,adftype="ct") # ADF with 3 lags: residuals OK? OK

best_dadf

## Utilisation des autres tests pour conforter les résultats
# KPSS test:
kpss.test(train_data$dindice[-1]) # stationarité
pp.test(train_data$dindice[-1])# stationnarité


### 3. Représenter graphiquement la série choisie avant et après transformation.----

windows()
# Création des deux graphiques
serie_brute <- ggplot(data = train_copie[-1,], aes(x = Periode, y = indice)) +
  geom_line() +
  labs(x = "", y = "Indice", title = "")+
  theme(panel.background = element_rect(fill = "lightgrey"))

serie_diff <- ggplot(data = train_copie[-1,], aes(x = Periode, y = dindice)) +
  geom_line() +
  labs(x = "", y = "Indice différencié", title = "")+
  theme(panel.background = element_rect(fill = "lightgrey"))

# Affichage des deux graphiques côte à côte
grid.arrange(serie_brute, serie_diff, nrow = 1)
