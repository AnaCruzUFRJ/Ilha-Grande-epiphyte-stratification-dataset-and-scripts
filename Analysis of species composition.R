
###############################################
#
# Date: 30-August-2021
# Name: Analysis of the variation in the composition of epiphytes along the different height zones of trees on Ilha Grande
# Author: ANA CAROLINA RODRIGUES DA CRUZ
# E-mail: anacruzufrj@gmail.com
# 

# Description: Difference in epiphyte composition along the five height zones of trees in the forests of Ilha Grande, multivariate analysis (PCOA) and similarity indices

# 
#
################################################

library (vegan)# for the hellinger function
library(dplyr) # data manipulation
library(ape) # to make pcoa
library(BiodiversityR) # to make pcoa

################################################
##### load directory

setwd("D:/Arquivos/Documents/ANA 2021/TESE/Tese 2021/Manuscritos/manuscrito 1/Scripts - repositorio")

################################################

dados<-read.table("composicao_por_zona.csv",header=T,sep=";")
# View(dados)
str(dados)
dados2<-dados%>%select(-Parcela, -Forofito, -Zonas)
#View(dados2)
zonas<-dados[,3]
#View(zonas)
rowSums(dados2)
#View(dados2)
rowSums(dados2)==0 # to check if it has zero
# which(rowSums(dados2)==0) show where is the sum equal to zero 
dados3<-decostand(dados2,"hellinger") # transform to solve asymmetry problems in abundance
dados3
dist.bray <- vegdist(dados3, method = "bray", binary = F)
dist.bray
#View(as.matrix(dist.bray))
cmdscale(dist.bray)
cmdscale(dist.bray,k=2)
pcoa<-cmdscale(dist.bray,k=2,eig=F) # to see the eigenvalues of each axis 
pcoa
plot(pcoa)

######### entering the names of the zones in the pcoa
plot(pcoa,main= "PCOA para zonas de altura da arvore", type="n",xlab="PCoA1",ylab="PCoA2")
points(pcoa[1:161,1],pcoa[1:161,2],pch=20,col="dark green")# ZONA 1 
points(pcoa[162:353,1],pcoa[162:353,2],pch=20,col="dark red")# ZONA 2
points(pcoa[354:484,1],pcoa[354:484,2],pch=20,col="blue")# ZONA 3 
points(pcoa[485:590,1],pcoa[485:590,2],pch=20,col="black")# ZONA 4
points(pcoa[591:634,1],pcoa[591:634,2],pch=20,col="black")# ZONA 5

####################################################################################################################################################################################

### another way to make pcoa ###

SPP<-dados2; SP<-vegdist(SPP,"bray")
View(SPP)
coordp<-cmdscale(SP,k=2,eig=TRUE, add=F)
coordp # escores das parcelas
options(scipen=999) 
autovalores<-round((coordp$eig/sum(coordp$eig))*100,0)
ordiplot(coordp)
PCOA<-pcoa(SP)
PCOA$values[1:4,] # eigenvalues
# always use relative eigenvalues
PCOA$vectors[,1:2] # plot scores
biplot(PCOA)
PCOA $ points [ , 2 ]
SPP_SCORES <- add.spec.scores(coordp, SPP, method='cor.scores', Rscale=F, scaling=1, multi=1) # # to add species scores. See method in help
ordiplot(SPP_SCORES)
par(mfrow=c(1,1))
###########
grafico <- data.frame ( "Zonas" = as.factor(dados$Zonas), "PCoA1" = pcoa[, 1], "PCoA2" = pcoa[, 2])
#View(grafico)

colores <- c("#f4f9f3", "#9ec39a", "#8af97f", "#347c2c", "#072c04")
cores <- c("#f4f9f3", "#9ec39a", "#8af97f", "#347c2c", "#072c04")[grafico$Zonas]
# Zona 1 = green, Zona 2 = grey, Zona 3 = blue, Zona 4 = orange, Zona 5 = red
grafico$cores <- cores
par(mar = c(5, 4, 4, 2) + 0.1)

par(mar = c(5, 4, 4, 4) + 0.1, xpd = TRUE)
plot(PCoA2 ~ PCoA1, data = grafico, type = "n",
     xlab = "PCoA1 (17.7%)", ylab = "PCoA2 (13.7%)")
points(PCoA2 ~ PCoA1, data = grafico, pch = 21, bg = cores)
segments(x0 = -0.56, y0 = 0, x1 = 0.43 , y1 = 0, lty = 2, col = "grey70")
segments(x0 = 0, y0 = -.5, x1 = 0 , y1 = 0.6, lty = 2, col = "grey70")
legend(.50, .4, legend = c("1","2", "3", "4", "5"), pch = 21, pt.bg = colores, title = "Zones", bty = "n")

# ggsave("PCOA_Zonas.jpeg", units="in", dpi=600)

########################################################## 
########################################################
## testing whether the difference is significant by the Multivariate homogeneity of groups dispersions (variances) test ######################################################## 
##########################################################

ZONAS<-as.factor(zonas)
ZONAS
comparacao<- betadisper(d=dist.bray, group=ZONAS)

permutest(comparacao,pairwise =T)
boxplot(comparacao)
plot(comparacao)

######################################################################################################################################################################################################
#                                    #
#                                    #
# analyzing similarity using bray curtis distance #
#                                    #
#                                    #
######################################################################################################################################################################################################

setwd("D:/Arquivos/Documents/ANA 2021/TESE/Tese 2021/analise da composicao")
dados<-read.table("composicao_por_zona_bray.csv",header=T,sep=";")
# View(dados)
str(dados)
clust<-vegdist (dados, method = "bray")
clust # dissimilarity values
par(mfrow=c(1,1)) 
plot(hclust(clust, method="average"), hang=-2, main = "", xlab = "Zones") # we chose not to use the dendrogram in the manuscript

#########################

