###############################################################################
#                                                                             #
#                PROJET DE SERIE TEMPORELLES                                  #
#                                                                             #
###############################################################################

## Importation des libraries
require(zoo) #format de serie temporelle pratique et facile d'utilisation (mais plus volumineux)
require(tseries) #diverses fonctions sur les series temporelles

library(dplyr)
library(readxl)
library(ggplot2)

IPI <- read_excel("C:/Users/LENOVO/Desktop/PROJET_serie_temp/IPI_202501.xls")

summary(IPI)

str(IPI)

# Créer une série temporelle (exemple : ventes mensuelles de 2015 à 2024)
data_ts <- ts(IPI$BE, start = c(1990, 1), frequency = 12)

# Affichage du graphique
plot(data_ts, main = "Évolution des BE", ylab = "Ventes", xlab = "Années", col = "blue", lwd = 2)





# Création de séries temporelles à partir du dataframe IPI
data_ts1 <- ts(IPI$BE, start = c(1990, 1), frequency = 12)  # Série 1
data_ts2 <- ts(IPI$CZ, start = c(1990, 1), frequency = 12)  # Série 2
data_ts3 <- ts(IPI$DE, start = c(1990, 1), frequency = 12)  # Série 3
data_ts4 <- ts(IPI$FZ, start = c(1990, 1), frequency = 12)  # Série 4

# Tracé de la première série
plot(data_ts1, type = "l", col = "blue", lwd = 1.5, 
     ylim = range(c(data_ts1, data_ts2, data_ts3, data_ts4)), 
     main = "Comparaison de quatre séries temporelles", 
     xlab = "Années", ylab = "Valeur")

# Ajout des autres séries temporelles
lines(data_ts2, col = "red", lwd = 1.5)
lines(data_ts3, col = "yellow", lwd = 1.5)
lines(data_ts4, col = "black", lwd = 1.5)

# Ajout de la légende
legend("bottomleft", 
       legend = c("Indus ", "Indus manu", "Indus Extrac", "Construction"), 
       col = c("blue", "red", "yellow", "black"), 
       lwd = 1.5, box.lwd = 0.1, cex=0.4)


### SERIE DES PRODUITS INFORMATIQUES 

info_ts <- ts(IPI$`[CI]`,start=c(1990,1), frequency= 12)
agro_ts <- ts(IPI$`(C1)`,start=c(1990,1), frequency= 12)
plot(info_ts, type = "l", col = "blue", lwd = 1.5, 
     ylim = range(c(data_ts1, data_ts2, data_ts3, data_ts4)), 
     xlab = "Années", ylab = "Valeur")
lines(agro_ts, col = "red", lwd = 1.5)

### DAIFFERENCIATION DE LA SERIE