#####################################################################################################################################################
#
# Name:Analysis of grouping indicator species for the height zones of trees on Ilha Grande
# Author: ANA CAROLINA RODRIGUES DA CRUZ
# E-mail: anacruzufrj@gmail.com
# Data: 21 - September - 2021
#
#####################################################################################################################################################
#
# Dataset referring to the manuscript: Importance of the vertical gradient in the variation of epiphyte community structure in the Brazilian Atlantic Forest submitted for evaluation by a scientific journal
#
library(dplyr) #  data manipulation
library(labdsv)# to calculate the "Indicator Value"
library(indicspecies) # ditto previous
#####################################################################################################################################################

# load directory #

setwd("D:/Arquivos/Documents/ANA 2021/TESE/Tese 2021/Manuscritos/manuscrito 1/Scripts - repositorio")
dados<-read.table("composicao.geral.zonas.csv",header=T,sep=";")
#View(dados)
dados2<-dados[,-c(1:3)]
#View(dados2)
zonas<-dados[,3]
dis.comunidade <- dsvdis(dados2,'bray/curtis')
cluster<-hclust(dis.comunidade,method="complete")
plot(cluster,labels=zonas) #confused! 

# now let's see the indicator species # 

indval <- multipatt(dados2, zonas,control = how(nperm=1000))
indval
summary(indval)
summary(indval, indvalcomp=TRUE)

# no combination of groups

indvalori <- multipatt(dados2, zonas, duleg = TRUE, control = how(nperm=999)) 

summary(indvalori)
################
# Analyzing the ecological preferences of species with correlation indices # 
################

dados2 <- ifelse(dados2>0,1,0) > phi = multipatt(dados2, zonas, func = "r", control = how(nperm=999)) 

phi <- multipatt(dados2, zonas, func = "r.g", control = how(nperm=999))

summary(phi)
round(head(phi$str),3)

