#=====================
# Chargement des bibliothèques nécessaires
#=====================

library(summarytools)
library(evd)
library(tseries) 
library(ismev)

# =====================
# Préparation des données
# =====================

data <- read.csv("C:/Users/hrywa/OneDrive/Bureau/M2/Statistiques_spatiales/BRINDISI.txt", sep="")
data[is.na(data)] <- 20

# =====================
# Analyse exploratoire des données
# =====================

# Histogramme des Températures
hist(data$TX, xlab='Température', main='Histogramme des Températures')

# Résumé des données
view(dfSummary(data))

# Graphique de série temporelle
plot(data$Years, data$TX, type="l",xlab='Années',ylab='Température')

# Calcul des maxima annuels
annees <- unique(data$Years)
max_annee <- numeric(length(annees))
for (i in seq_along(annees)) {
  max_annee[i] <- max(data[data$Years == annees[i], ]$TX)
}
plot(annees, max_annee, type="l",xlab='Années',ylab='Température', main="température maximale par an")

# =====================
# Test de stationnarité
# =====================

# Application du test ADF sur les données
result <- adf.test(data$TX)
print(result) 

result2 <- adf.test(max_annee)
print(result2) 

# =====================
# Analyse GEV des maximums par année
# =====================

# Ajustement et analyse GEV
gev.fit <- fgev(max_annee)
plot(profile(gev.fit, which = names(gev.fit$estimate[1])))
plot(profile(gev.fit, which = names(gev.fit$estimate[2])))
plot(profile(gev.fit, which = names(gev.fit$estimate[3])))
confint(gev.fit)

#=========================
# Calcul des niveaux de retour 
#=========================

variance=gev.fit$var.cov
scale=gev.fit$estimate[2]
shape=gev.fit$estimate[3]
loc=gev.fit$estimate[1]

# Niveaux de retour pour T=20 et T=100 ans avec méthode Delta
y20 <- -log(1 - 1/20)
y100 <- -log(1 - 1/100)

# Niveau de retour pour T=20 et T=100 ans
Z20 <- loc - (scale / shape) * (1 - (y20 ** shape))
Z100 <- loc - (scale / shape) * (1 - (y100 ** shape))

# Méthode Delta pour T=20
delta_20 <- c(1, shape ** -1 * (1 - (-log(1 - 1/20)) ** shape),
              scale * shape ** -2 * (1 - (-log(1 - 1/20)) ** shape) - scale * shape ** -1 * (-log(1 - 1/20)) ** shape * log(-log(1 - 1/20)))
var_20 <- t(delta_20) %*% variance %*% delta_20

# Méthode Delta pour T=100
delta_100 <- c(1, shape ** -1 * (1 - (-log(1 - 1/100)) ** shape),
               scale * shape ** -2 * (1 - (-log(1 - 1/100)) ** shape) - scale * shape ** -1 * (-log(1 - 1/100)) ** shape * log(-log(1 - 1/100)))
var_100 <- t(delta_100) %*% variance %*% delta_100


# Calcul des intervalles de confiance à 95%
alpha <- 0.05
z_alpha <- qnorm(1 - alpha/2)

# Intervalle de confiance pour Z20
lower_Z20 <- Z20 - z_alpha * sqrt(var_20)
upper_Z20 <- Z20 + z_alpha * sqrt(var_20)

# Intervalle de confiance pour Z100
lower_Z100 <- Z100 - z_alpha * sqrt(var_100)
upper_Z100 <- Z100 + z_alpha * sqrt(var_100)

# Affichage des intervalles de confiance
cat("Intervalle de confiance à 95% pour Z20 : [", lower_Z20, ", ", upper_Z20, "]\n")
cat("Intervalle de confiance à 95% pour Z100 : [", lower_Z100, ", ", upper_Z100, "]\n")


# =====================
# Approche GPD
# =====================

# Préparation de la série temporelle pour l'approche GPD
température_ts <- ts(max_annee, start=1951, end=2023)
mrlplot(température_ts) # Pour identifier visuellement le seuil
plot(data$TX, type='l', main="Températures avec Seuils", xlab="Temps", ylab="Température")
abline(h=37, col="blue", lty=2) 

# Modèle GPD
mod <- fpot(data$TX, 37)
mod

#===============================
#calcul a la main des niveau de retour
#===============================
# Niveau de retour pour T=20 et T=100 ans
C= sum(data$TX > 37)/length(data$TX)

R100=37+(1.96179  /(-0.05828  ))*((100*365*C)^-0.05828  -1)  # *365 car periode de retoure en années et nos données sont journaliere
R100

R20=37+(1.96179  /(-0.05828  ))*((20*365*C)^-0.05828  -1)  # *365 car periode de retoure en années et nos données sont journaliere
R20

matrice=mod$var.cov
zeta=-0.05828  #la shape de mod
V=cbind(c(zeta*(1-zeta)/length(data$TX),0,0),rbind(rep(0,2),mod$var.cov))
delta3 <- c(
  1.96179  ^(-0.05828  ) * C^(-0.05828   - 1),
  1 / -0.05828   * ((100 * 365 * C)^-0.05828   - 1),
  -1.96179   * (-0.05828  )^(-2) * ((100 * 365 * C)^-0.05828   - 1) + 1.96179   * 1 / (-0.05828  ) * (100 * 365 * C)^(-0.05828  ) * log(100 * 365 * C)
)
V_x100=t(delta3)%*% V %*%delta3
V_x100

delta4 <- c(
  1.96179  ^(-0.05828  ) * C^(-0.05828   - 1),
  1 / -0.05828   * ((20 * 365 * C)^-0.05828   - 1),
  -1.96179   * (-0.05828  )^(-2) * ((20 * 365 * C)^-0.05828   - 1) + 1.96179   * 1 / (-0.05828  ) * (20 * 365 * C)^(-0.05828  ) * log(20 * 365 * C)
)
V_x20=t(delta4)%*% V %*%delta4
V_x20

#vérification avec fpot
Mod100=fpot(data$TX,37,npp=365,mper=100)
Mod20=fpot(data$TX,37,npp=365,mper=20)


# =====================
# Affichage des niveaux de retour et variances pour T=20 et T=100 ans
# =====================
Mod100
R100
V_x100


Mod20
R20
V_x20

#============
#VAlidation du modele
#=================


#======================
#Autre modele a seuil avec cluster 
#======================

#reorésentation des données
plot(data$TX, type='l', main="Températures avec Seuils", xlab="Temps", ylab="Température")
abline(h=37, col="blue", lty=2) 


#Paramétre
r1=2
r2=4
u1=37

#représentation des clusters
clusters(data$TX,u=u1,r=r1,plot=TRUE,col="white",xlab="Années",ylab="température d'hiver minimales")
clusters(data$TX,u=u1,r=r2,plot=TRUE,col="white",xlab="Années",ylab="température d'hiver minimales")


#Estimation du modèle GPD pour les données au-dessus du seuil avec différents r
Mod1=fpot(data$TX,threshold = u1,cmax=TRUE,r=r1)
Mod2=fpot(data$TX,threshold = u1,cmax=TRUE,r=r2)


# Estimation de l'indice extremal θ pour chaque modèle
exi(data$TX,u=u1,r=r1)
exi(data$TX,u=u1,r=r2)

#les niveaux de retour
fpot(data$TX,u1,npp=365,mper=100,cmax=TRUE,r=r1)

fpot(data$TX,u1,npp=365,mper=100,cmax=TRUE,r=r2)

#===================
#Si les données sont non stationnaire
#===================


# Ajustement d'un modèle GEV avec une tendance linéaire sur les données de température
Annees <- 1951:2023  # Assurez-vous que cela correspond à la plage de vos données
mod3 <- fgev(max_annee, nsloc =1: length(max_annee))
mod3

# Tracé des données et de la tendance estimée par le modèle
plot(Annees,max_annee, type = 'o', main = "Valeurs Extrêmes avec Tendance Linéaire", xlab = "Année", ylab = "Valeurs")
Mu <- mod3$estimate[1] + mod3$estimate[2] * (Annees - min(Annees))  # Ajustement pour la tendance linéaire
lines(Annees, Mu, type = "l", col = "red")



#Calcul de vraissemblance
l1=-sum(log(mod3$estimat[3])+(1+1/mod3$estimate[4])*log(1+mod3$estimate[4]*(max_annee-Mu)/mod3$estimate[3])+(1+mod3$estimate[4]*((max_annee-Mu)/mod3$estimate[3]))^(-1/mod3$estimate[4]))

#verification avec la fonction
Mod3bis=gev.fit(max_annee,ydat=matrix(1:length(max_annee),ncol=1),mul=1)
#nllh single numeric giving the negative log-likelihood value.  117.9764 la vraissemblance   



#autre modele classique
m=length(max_annee)
Mod_iid=fgev(max_annee)
l2=-m*log(Mod_iid$estimate[2])-(1+1/Mod_iid$estimate[3])*sum(log(1+Mod_iid$estimate[3]*
                                                                   (max_annee-Mod_iid$estimate[1])/Mod_iid$estimate[2]))-sum((1+Mod_iid$estimate[3]*(max_annee-Mod_iid$estimate[1])/Mod_iid$estimate[2])^(-1/Mod_iid$estimate[3]))

# Calcul de la déviance
deviance <- 2 * (l1 - l2) #du grand moins le petit #a commparer a qchisq(0.95,1)

# Nombre de degrés de liberté (différence dans le nombre de paramètres)
df <- length(mod3$estimate) - length(Mod_iid$estimate)

# Test du Chi-carré
p_value <- pchisq(deviance, df, lower.tail = FALSE)

# Afficher les résultats
print(c("Deviance" = deviance, "P-Value" = p_value))

