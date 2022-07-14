
###############################################
#                                             #
# Date: 14-July-2022
# Name: Analysis of the sampling sufficiency for the collection of epiphyte of the Brazilian Atlantic Forest 
# Author: ANA CAROLINA RODRIGUES DA CRUZ
# E-mail: anacruzufrj@gmail.com
#                                             #
## Description:recognize, execute and interpret the accumulation curve, richness estimators and species rarefaction curve for data from Ilha Grande, RJ
# 
##
# Dataset referring to the manuscript: Importance of the vertical gradient in the variation of epiphyte community structure in the Brazilian Atlantic Forest submitted for evaluation by a scientific journal
#
#                                           ################################################
#                                       #


# install.packages(BiodiversityR)
# install.packages(vegan)

library(vegan)
library(BiodiversityR)
library(tidyverse)

############## Acumulation curve -  plot 1 - restinga 1

# setwd("D:/Arquivos/Documents/ANA 2021/TESE")


REST1<-read.table("composição por segmento_restinga_O1.csv",header=T,sep=";", dec=",")
View(REST1)

# Acumulation curve os species

acum<-specaccum(REST1)
acum
plot(acum, col="black",xlab="Sample units", ylab="Accumulated wealth of species", main="Plot 1")
# or

plot(acum, lwd=1,ci.col="grey", xlab="Sampling effort", ylab="Accumulated species richness", col="black", main="REST1", xlim=c(0,30), ylim=c(0,70))

# Richness estimators - chao 1, Jackknife 1 e 2 e Bootstrap 

Chao22<-specpool(REST1)
Chao22 # 

Chao2<-poolaccum(REST1, permutations=100)
summary(Chao2, display = "chao")

plot(Chao2, alpha = 0.05, type = c("l","g"), lwd=1,ci.col="grey", xlab="Sampling effort", ylab="Accumulated species richness", main="")

# Plot 2 -  restinga 2


REST2<-read.table("composição por segmento_restinga_O2.csv",header=T,sep=";", dec=",")

rarecurve(REST1, col = "firebrick", cex = 0.6, xlab = "Tamanho amostral",ylab="Espécies")
# Acumulation curve os species
acum<-specaccum(REST2)
acum

plot(acum, lwd=1,ci.col="grey", xlab="Sampling effort", ylab="Accumulated species richness", col="black", main="REST2", xlim=c(0,30), ylim=c(0,70))

Chao22<-specpool(REST2)
Chao22 

Chao2<-poolaccum(REST2, permutations=100)
summary(Chao2, display = "chao")
plot(Chao2, alpha = 0.05, type = c("l","g"), lwd=1,ci.col="grey", xlab="Unidades amostrais", ylab="Especies", main="Restinga 2")

# plot 3 - parnaioca


Parna<-read.table("composição por segmento_parnaioca.csv",header=T,sep=";", dec=",")


Chao22<-specpool(Parna)
Chao22 

# Acumulation curve os species
acum<-specaccum(Parna)
acum

plot(acum, lwd=1,ci.col="grey", xlab="Sampling affort", ylab="Accumulated species richness", col="black", main="DOLF", xlim=c(0,30), ylim=c(0,70))

#

Chao2<-poolaccum(Parna, permutations=100)
summary(Chao2, display = "chao")
plot(Chao2, alpha = 0.05, type = c("l","g"), lwd=1,ci.col="grey", xlab="Unidades amostrais", ylab="Especies", main="Parnaioca")

##### plot 4 - poco do soldado 

poco<-read.table("composição por segmento_poço.csv",header=T,sep=";", dec=",")

# Acumulation curve os species
acum<-specaccum(poco)
acum

plot(acum, lwd=1,ci.col="grey", xlab="Sampling effort", ylab="Accumulated species richness", col="black", main="DOSF1", xlim=c(0,30), ylim=c(0,70))

Chao22<-specpool(poco)
Chao22 

Chao2<-poolaccum(poco, permutations=100)
summary(Chao2, display = "chao")
plot(Chao2, alpha = 0.05, type = c("l","g"), lwd=1,ci.col="grey", xlab="Unidades amostrais", ylab="Especies", main="Poco do Soldado")

### plot 5

brit<-read.table("composição por segmento_britador.csv",header=T,sep=";", dec=",")

# Acumulation curve os species
acum<-specaccum(brit)
acum

plot(acum, lwd=1,ci.col="grey", xlab="Sampling effort", ylab="Accumulated species richness", col="black", main="DOSF2", xlim=c(0,30), ylim=c(0,70))


# 
Chao22<-specpool(brit)
Chao22 

Chao2<-poolaccum(brit, permutations=100)
summary(Chao2, display = "chao")
plot(Chao2, alpha = 0.05, type = c("l","g"), lwd=1,ci.col="grey", xlab="Unidades amostrais", ylab="Especies", main="Britador")


# Diversity indices

# restinga 1 - plot 1
divers.all <- diversityresult(REST1, index = "Shannon", method="pooled")  
divers.all
divers.all <- diversityresult(REST1, index = "Simpson", method="pooled") 
divers.all

# restinga 2 - plot 2
divers.all <- diversityresult(REST2, index = "Shannon", method="pooled")  
divers.all
divers.all <- diversityresult(REST2, index = "Simpson", method="pooled") 
divers.all

# parnaioca - plot 3
divers.all <- diversityresult(Parna, index = "Shannon", method="pooled")  
divers.all
divers.all <- diversityresult(Parna, index = "Simpson", method="pooled") 
divers.all

# poco do soldado - plot 4

divers.all <- diversityresult(poco, index = "Shannon", method="pooled")  
divers.all
divers.all <- diversityresult(poco, index = "Simpson", method="pooled") 
divers.all

# britador - plot 5


divers.all <- diversityresult(brit, index = "Shannon", method="pooled")  
divers.all
divers.all <- diversityresult(brit, index = "Simpson", method="pooled") 
divers.all
