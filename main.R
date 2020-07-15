############## ############## ############## ############## ############## ############## 
#### Annalyse Modele D en regression lineaire et regression de Poisson
############## ############## ############## ############## ############## ############## 
############## Data
inp <- read.csv("data_lot2.csv", header=T,sep=",")[,c(4:7,11,14,15,17,16,13,18:23,12,9)]; head(inp) # type de commune
dim(inp)
Pop <- read.csv("data_lot2.csv", header=T,sep=",")[,10]
############## Regession linear
regressor1 = lm(formula = occurrence ~ 
                        TauxBac + TauxChomage + Salaire + TauxOuvrier + TauxClub
                + TauxBacLag + TauxChomageLag + SalaireLag + TauxOuvrierLag + TauxClubLag 
                + APL + ALD30 + TauxALD30
                + APLLag + ALD30Lag + TauxALD30Lag + yLag , data = inp)
summary(regressor1)  
summary(regressor1$residuals) # Adj-R2 = 0.8411

regressor1_ = lm(formula = occurrence ~ 
                TauxBac + TauxChomage + Salaire  
                 + TauxChomageLag + SalaireLag + TauxOuvrierLag + TauxClubLag 
                 + APL + ALD30 + TauxALD30
                 + APLLag    , data = inp)
summary(regressor1_)  ## Adj-R2 = 0.8514
summary(regressor1_$residuals) 

z1 <- regressor1_$residuals

############# Regession Poisson 
exposure <- inp$occurrence/Pop
regressor2 = glm(formula = occurrence ~ 
                         TauxBac + TauxChomage + Salaire + TauxOuvrier + TauxClub
                 + TauxBacLag + TauxChomageLag + SalaireLag + TauxOuvrierLag + TauxClubLag 
                 + APL + ALD30 + TauxALD30
                 + APLLag + ALD30Lag + TauxALD30Lag + yLag 
                 + offset(log(exposure)), 
                 data = inp, family = poisson(link = log))
summary(regressor2) # AIC = 468.07
summary(residuals(regressor2))

regressor2_ = glm(formula = occurrence ~ 
                          TauxClub
                  + TauxChomage + SalaireLag + TauxClubLag 
                  + ALD30 + TauxALD30
                  + ALD30Lag   
                  + offset(log(exposure)), 
                  data = inp, family = poisson(link = log))
summary(regressor2_) # AIC = 452.77
summary(residuals(regressor2_))
z2 <- residuals(regressor2_)

#### plot each observation predicted versus actual
par(mfrow=c(1,2))

y_pred1 <- regressor1_$fitted.values
plot(inp$occurrence, type = "n",y_pred1,lwd = 3, xlab="Original", ylab = "Predicted"
     ,cex.lab=1.5, main = "Régression Linéaire ", cex.main=1.4)
text(inp$occurrence,y_pred1, col = "#5592e3",lwd = 3)
abline(a = 0, b = 1, col = "red", lwd = 2)

y_pred2 <- regressor2_$fitted.values 
plot(inp$occurrence, type = "n",y_pred2,lwd = 3, xlab="Original", ylab = "Predicted"
     ,cex.lab=1.5, main = "Régression de Poisson ", cex.main=1.4)
text(inp$occurrence,y_pred2, col = "#5592e3",lwd = 3)
abline(a = 0, b = 1, col = "red", lwd = 2)

summary(z1); summary(z2);

############## ############## ############## ############## ############## ############## 
#### La comparaison : modele C et modele D
############## ############## ############## ############## ############## ############## 
Dep <- read.csv("data_lot2.csv", header=T,sep=",")[,2]; Dep # type de commune
indD <- read.csv("indD.csv", header=T,sep=",")[,3]; indD # type de commune
###### Plot carte de Residus/Residus standardisés
###### La France 
library(cartography)
library(raster)
z <- z2; summary(z)
formes <- getData(name="GADM", country="FRA", level=2)
concordance <- z[indD]; concordance
formes$z <- concordance
mycols <- carto.pal(pal1 = "blue.pal", n1 = 9, pal2 = "orange.pal", n2 = 9)
ech <-  c(-3.45, -3.4, -2.9, -2.36, -2, -1.5, -1, 0, 1,1.5,2.3,2.5,3.7,4.5,4.77,4.9)
spplot(formes, "z", col.regions=mycols,
       at = ech,#couleurs(30), 
       main=list(label="Résidus standardisés, Durbin-spatial modèle D",cex=1.4))
spplot(formes, "z", col.regions=mycols,
       at = ech,#couleurs(30), 
       main=list(label="Résidus, Durbin-spatial modèle D",cex=1.4))

###### L'Ile-de-France
France <- getData(name="GADM", country="FRA", level=2)
FrancePartie<- subset(France, NAME_1=="Île-de-France")
#plot(FrancePartie, main="Carte de la France")
nI <- FrancePartie$NAME_2; nI
idI <- c(91, 92, 75, 77, 93, 95, 94, 78)

concordanceI <- z[match(idI, Dep)]; concordanceI #
FrancePartie$z <- concordanceI
mycols <- carto.pal(pal1 = "blue.pal", n1 = 9, pal2 = "orange.pal", n2 = 9)
spplot(FrancePartie, "z", col.regions=mycols,
       at =ech ,#couleurs(30), 
       main=list(label="Ile-de-France",cex=1.8))

###### Test Moran I
geo_dep <- read.csv("geo_dep.csv", header=T,sep=","); head(geo_dep); dim(geo_dep)
### matrix of inverse distance weights
geo.dists <- as.matrix(dist(cbind(geo_dep$lon, geo_dep$lat))); head(geo.dists)

geo.dists.inv <- 1/geo.dists
diag(geo.dists.inv) <- 0
geo.dists.inv[1:5, 1:5]

install.packages("ape")
library(ape)

z <- z2
Moran.I(z, geo.dists.inv)$p.value

############## ############## ############## ############## ############## ############## 
#### Annalyse Modele D en regression lineaire et regression de Poisson PLS
############## ############## ############## ############## ############## ############## 
### Regression PLS Poisson
install.packages("plsRglm")
library(plsRglm)
library(pls)

regressor3 <- plsRglm(inp[,18],cbind(inp[, c(1:17)], log(exposure)),7,modele="pls-glm-family",family="poisson",
                      pvals.expli=TRUE)
regressor3
summary(regressor3$InfCrit)
z3 <- regressor3$Yresidus
summary(z3)

y_pred3 <- regressor3$YChapeau

#plot each observation predicted versus actual
par(mfrow=c(1,3))
plot(inp$occurrence, type = "n",y_pred1,lwd = 3, xlab="Original", ylab = "Predicted"
     ,cex.lab=1.5, main = "Régression Linéaire ", cex.main=1.4)
text(inp$occurrence,y_pred1, col = "#5592e3",lwd = 3)
abline(a = 0, b = 1, col = "red", lwd = 2)

plot(inp$occurrence, type = "n",y_pred2,lwd = 3, xlab="Original", ylab = "Predicted"
     ,cex.lab=1.5, main = "Régression de Poisson ", cex.main=1.4)
text(inp$occurrence,y_pred2, col = "#5592e3",lwd = 3)
abline(a = 0, b = 1, col = "red", lwd = 2)

plot(inp$occurrence, type = "n",y_pred3,lwd = 3, xlab="Original", ylab = "Predicted"
     ,cex.lab=1.5, main = "Régression de Poisson PLS", cex.main=1.4)
text(inp$occurrence,y_pred3, col = "#5592e3",lwd = 3)
abline(a = 0, b = 1, col = "red", lwd = 2)

#### Boxplot des residus / residus standardises
par(mfrow=c(1,2))
summary(z1); summary(z2); summary(z3)
boxplot(cbind(z1,z2,z3), names = c("RegLineaire","RegPoisson", "RegPoissonPLS"), 
        main = "Résidus")
s1 <- scale(z1); s2 <- scale(z2); s3 <- scale(z3)
summary(s2); summary(s3)
boxplot(cbind(s1,s2,s3), names = c("RegLineaire", "RegPoisson", "RegPoissonPLS"), 
        main = "Résidus Standardisés")









