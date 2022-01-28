
###############################################
#
# 
# Date: 23-August-2021
# Nome: Analysis of general data on the vertical stratification of epiphytes by tree height zones
# Author: ANA CAROLINA RODRIGUES DA CRUZ
# E-mail: anacruzufrj@gmail.com
# 
#
# Description: comparison of the richness and abundance of epiphytes along the five tree height zones in the forests of Ilha Grande
#
# Dataset referring to the manuscript: Importance of the vertical gradient in the variation of epiphyte community structure in the Brazilian Atlantic Forest submitted for evaluation by a scientific journal
#
#
################################################


library(psych) #to extract info from media, median etc
library(car) # test homoscedasticity
library(dplyr) # data manipulation
library(AID) ### need to load this package to make box cox
library(MASS) # is also used for boxcox
library(tidyverse) # data manipulation
library(ggplot2) #graphs
library(ggthemes)# graph of GGPLOT
library(RColorBrewer) # graph of GGPLOT
library(bbmle) # use AICctab function from bbmle package
library(jtools)# more detailed output of models output

##############################################
# for GLM 

library(lme4) # 
library(multcomp) # for glm tukey post hoc test (GLM)
library(emmeans) # for glm tukey post hoc test (GLMM)

##### load directory

setwd("D:/Arquivos/Documents/ANA 2021/TESE/Tese 2021/Manuscritos/manuscrito 1/Scripts - repositório")

# Analysis of richness: between line 160 and 360
# Analysis of abundance: from line 360

ZONAS<-read.table("Riqueza.por.zona.csv",header=T,sep=";")
head(ZONAS)
# View(ZONAS)
ZONAS$Zones<-as.factor(ZONAS$Zones)
str(ZONAS)
#structure(ZONAS)

# descriptive analysis # 

zonas_separadas<-unstack(ZONAS[,2:3])
#View(zonas_separadas)
describe(zonas_separadas[,1])
#View(describe(zonas_separadas[,1]))

describe(zonas_separadas[,2])
#View(describe(zonas_separadas[,2]))

describe(zonas_separadas[,3])
#View(describe(zonas_separadas[,3]))

describe(zonas_separadas[,4])
#View(describe(zonas_separadas[,4]))

describe(zonas_separadas[,5])
#View(describe(zonas_separadas[,5]))

########################################
##### Richness variation by zone #####
#########################################

plot(Richness~Zones, data=ZONAS,center=F, main = "richness of epiphytes by height zone") # to get an idea of how the data is in each zone

# testing the normality of the data # 

shapiro.test(ZONAS$Richness)
# there is no normality of the data 

leveneTest (Richness ~ as.factor(Zones), data = ZONAS)

# there is no normality of the data

### very low p-value, that is, we reject HO and accept H1: there is no normality of the data and no homogeneity of the variances

hist(ZONAS$Richness, main = "Histogram of the frequency distribution of the data", xlab = "Richness", ylab="Frequency") # note that they are shifted to the left of the graph, they do not have a gaussian or normal distribution

# making an anova to remove the residues and test the normality of the residues
ANOVA_ZONAS<-aov(Richness~Zones, data=ZONAS)
shapiro.test(resid(ANOVA_ZONAS)) 
hist(resid(ANOVA_ZONAS), main = "Histograma da distribuicao dos residuos", xlab = "Riqueza de especies", ylab="Frequencia") # note that they are shifted to the left of the graph, they do not have a gaussian or normal distribution

# there is also no normality of residues

## attempt to transform and adjust normality

BOX_COX<-ZONAS %>%
  group_by(Zones)%>%
  do(COX=boxcoxnc(1+.$Richness))

# even after the transformations are still not normal

par(mfrow=c(1,1))
# a way to present data in graphs
corboxplot<-c('1'='pink', '2'='dark red','3'='dark green', '4'='blue', '5'='black')
boxplot(Richness~Zones,data=ZONAS, main = "Riqueza de epifitas por zona de altura da arvore",
        xlab="Zonas de altura",
        ylab=parse(text=paste("Riqueza")),
        font=2,
        cex=1,
        cex.lab=1, 
        cex.axis=1, 
        col=corboxplot)

stripchart(Richness~Zones, data=ZONAS, pch=21, cex=0.5, bg="red", vertical=T, add=T, method="jitter")

#########################################################
# another way to make the graph # 
# 
pp<-ggplot(ZONAS, aes(x=Zones, y=Richness, fill=Zones)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.fill = NULL,
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = NULL) + 
 # geom_point(position=position_jitter(width=0.3, height = 0.05, seed=1),shape=21,fill="grey80",size=0.6) + 
  scale_fill_brewer(palette="greens") + 
  labs(x="Zones",y="Richness") + 
  theme_bw() 

pp
#########################################################
# graph with mean and standard deviation is better to present the information, we chose to use this one below

Richness <- c(0.78,1.25,0.74,0.52,0.20)
sd <- c(0.85,1.26,1.02,0.52,0.56)
Zones<-c(1,2,3,4,5)

dat_zones <- data.frame(Zones, Richness, sd)
dat_zones$Zones<-as.factor(dat_zones$Zones)
str(dat_zones$Zones)

p_media_desvio <- ggplot(dat_zones, aes(factor(Zones), Richness)) + 
  geom_point(stat = "identity", size = 2) + 
  geom_errorbar(aes(ymin = Richness - sd, ymax =  Richness + sd), width = 0.2) +
  xlab("Zones")+ 
  theme_apa() 

p_media_desvio

# para salvar o grafico

###########################
# ggsave("media_sd_riqueza.jpeg", units="in", dpi=600)

############################################################################################################################ 
# testing the difference in richness by zone using a generalized linear model, with the Poisson distribution, since there is no normal distribution

# testing with GLM (not mixed) #

glm_model<-glm(Richness ~ Zones, data = ZONAS, family = "poisson")

summary(glm_model)
summ(glm_model)
# there is a difference between all zones and zone 1 which is in the intercept, but we don't know about the others.
# we need to do a post hoc test
# Tukey's post hoc test for generalized linear models was chosen with the multcomp package
summary(glht (glm_model, mcp (Zones = "Tukey"))) #perform TukeyHSD
# the comparison is pair by pair, therefore it is necessary to refer in the results which differences between the pairs of height zones are significant and which are not. 

#######################################################################################################################################
#######################################################################################################################################

# now using the mixed model, considering the plots as a random effect #

head(ZONAS)
glmm.zonas<-glmer(Richness ~ Zones + (1 | Forests), data = ZONAS, family = "poisson")
summary(glmm.zonas)
summ(glmm.zonas)
# POST HOC
m_meansall <- emmeans(glmm.zonas, ~ Zones, type = "response")
RESULTADO_RIQUEZA_GLMM<-pairs(m_meansall)
RESULTADO_RIQUEZA_GLMM

cld(m_meansall,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")
str(m_meansall)
###################################################

# comparing the models in search of the most suitable #

# comparison with null model #
glm.nulo<-glm(Richness ~ 1, data = ZONAS, family = "poisson")
summary(glm.nulo)
anova(glm_model,glm.nulo)
AICctab(glm_model, glm.nulo, base = T, weights = T)

# glm x glmm

AICctab(glmm.zonas, glm_model, base = T, weights = T)
# glmm is a better model than glm #
summ(glmm.zonas)

# write.table(RESULTADO_RIQUEZA_GLMM,"RESULTADO_RIQUEZA_GLMM.csv",row.names=F,sep=";", dec = ",") # save the results table

##############################################################################################################################################

# GLM for each of the installments separately# 

##############################################################################################################################################

head(ZONAS)

# restinga 1 #

# GlM
rest1<-ZONAS %>% filter (Forests == "Restinga1")
#View(rest1)
zonas_separadas<-unstack(rest1[,2:3])
#View(zonas_separadas)
mean(zonas_separadas[,1])
describe(zonas_separadas[,1])
#View(describe(zonas_separadas[,1]))
describe(zonas_separadas[,2])
#View(describe(zonas_separadas[,2]))
describe(zonas_separadas[,3])
#View(describe(zonas_separadas[,3]))
describe(zonas_separadas[,4])
#View(describe(zonas_separadas[,4]))
describe(zonas_separadas[,5])
#View(describe(zonas_separadas[,5]))

modelo.rest1<-glm(Richness ~ Zones, data = rest1, family = "poisson")
summary(modelo.rest1)
summ(modelo.rest1, confint = TRUE, digits = 3)
 # comparison with null model
glm.nulo<-glm(Richness ~ 1, data = rest1, family = "poisson")
summary(glm.nulo)
anova(modelo.rest1,glm.nulo)
AICctab(modelo.rest1, glm.nulo, base = T, weights = T)
summary(glht (modelo.rest1, mcp (Zones = "Tukey"),digits = 3)) 

# restinga 2 #

rest2<-ZONAS %>% filter (Forests == "Restinga2")
#View(rest2)
zonas_separadas<-unstack(rest2[,2:3])
#View(zonas_separadas)
mean(zonas_separadas[,1])
describe(zonas_separadas[,1])
#View(describe(zonas_separadas[,1]))
describe(zonas_separadas[,2])
#View(describe(zonas_separadas[,2]))
describe(zonas_separadas[,3])
#View(describe(zonas_separadas[,3]))
describe(zonas_separadas[,4])
#View(describe(zonas_separadas[,4]))
describe(zonas_separadas[,5])
#View(describe(zonas_separadas[,5]))

modelo.rest2<-glm(Richness ~ Zones, data = rest2, family = "poisson")
summ(modelo.rest2, confint = TRUE, digits = 3)
# comparison with null model
glm.nulo<-glm(Richness ~ 1, data = rest2, family = "poisson")
summary(glm.nulo)
anova(modelo.rest2,glm.nulo)
AICctab(modelo.rest2, glm.nulo, base = T, weights = T)
summary(glht (modelo.rest2, mcp (Zones = "Tukey"),digits = 3))

# Lowland Forest # 

terras.baixas<-ZONAS %>% filter (Forests == "FlorestaTerrasBaixas")
zonas_separadas<-unstack(terras.baixas[,2:3])
#View(zonas_separadas)
mean(zonas_separadas[,1])
describe(zonas_separadas[,1])
#View(describe(zonas_separadas[,1]))
describe(zonas_separadas[,2])
#View(describe(zonas_separadas[,2]))
describe(zonas_separadas[,3])
#View(describe(zonas_separadas[,3]))
describe(zonas_separadas[,4])
#View(describe(zonas_separadas[,4]))
describe(zonas_separadas[,5])
modelo.terras.baixas<-glm(Richness ~ Zones, data = terras.baixas, family = "poisson")
summ(modelo.terras.baixas, confint = TRUE, digits = 3)
summary(glht (modelo.terras.baixas, mcp (Zones = "Tukey"),digits = 3))
# comparison with null model
glm.nulo<-glm(Richness ~ 1, data = terras.baixas, family = "poisson")
summary(glm.nulo)
anova(modelo.terras.baixas,glm.nulo)
AICctab(modelo.terras.baixas, glm.nulo, base = T, weights = T)

# Submontane Forest 1 #

submontana1<-ZONAS %>% filter (Forests == "FlorestaSubmontana1")
zonas_separadas<-unstack(submontana1[,2:3])
#View(zonas_separadas)
mean(zonas_separadas[,1])
describe(zonas_separadas[,1])
#View(describe(zonas_separadas[,1]))
describe(zonas_separadas[,2])
#View(describe(zonas_separadas[,2]))
describe(zonas_separadas[,3])
#View(describe(zonas_separadas[,3]))
describe(zonas_separadas[,4])
#View(describe(zonas_separadas[,4]))
describe(zonas_separadas[,5])
modelo.submontana1<-glm(Richness ~ Zones, data = submontana1, family = "poisson")
summ(modelo.submontana1, confint = TRUE, digits = 3)
summary(glht (modelo.submontana1, mcp (Zones = "Tukey"),digits = 3))
# comparison with null model
glm.nulo<-glm(Richness ~ 1, data = submontana1, family = "poisson")
summary(glm.nulo)
anova(modelo.submontana1,glm.nulo)
AICctab(modelo.submontana1, glm.nulo, base = T, weights = T)

#  Submontane Forest 2 #

submontana2<-ZONAS %>% filter (Forests == "FlorestaSubmontana2")
zonas_separadas<-unstack(submontana2[,2:3])
#View(zonas_separadas)
mean(zonas_separadas[,1])
describe(zonas_separadas[,1])
#View(describe(zonas_separadas[,1]))
describe(zonas_separadas[,2])
#View(describe(zonas_separadas[,2]))
describe(zonas_separadas[,3])
#View(describe(zonas_separadas[,3]))
describe(zonas_separadas[,4])
#View(describe(zonas_separadas[,4]))
describe(zonas_separadas[,5])
modelo.submontana2<-glm(Richness ~ Zones, data = submontana2, family = "poisson")
summ(modelo.submontana2, confint = TRUE, digits = 3)
summary(glht (modelo.submontana2, mcp (Zones = "Tukey"),digits = 3))
# comparison with null model
glm.nulo<-glm(Richness ~ 1, data = submontana2, family = "poisson")
summary(glm.nulo)
anova(modelo.submontana2,glm.nulo)
AICctab(modelo.submontana2, glm.nulo, base = T, weights = T)

#######################################################
#######################################################
# forest graphics separately # 
#######################################################
#######################################################
Parc<-c(FlorestaSubmontana2="Submontana 2", FlorestaTerrasBaixas="Terras Baixas",FlorestaSubmontana1="Submontana 1",Restinga1="Restinga 1", Restinga2="Restinga 2")

ggplot(ZONAS, aes(x=Zones, y=Richness, fill=Zones))+
  geom_boxplot(outlier.shape=NA) + # remove the outlier
  geom_point(position=position_jitter(width=0.3, height = 0.05, seed=1),shape=21,fill="grey80",size=0.6)+ scale_fill_brewer(palette="greens")+ 
  facet_wrap( ~ Forests, labeller = labeller(Forests = Parc)) + 
  labs(x="Zones",y="Richness") + 
  theme_base() +
  theme(axis.title.x=element_blank(), # removing the x-axis labels from each graph individually
              axis.text.x=element_blank(), # remove text from x axis
              axis.ticks.x=element_blank()) # remove the dashes from the axis

#######################################################################
#############################################################################################################################################################################################################################################################

# Abundance by height zone # 

#######################################################################
######################################################################
######################################################################################################################################################################################

ZONAS<-read.table("abundancia.por.zona.csv",header=T,sep=";")
head(ZONAS)
# View(ZONAS)
ZONAS$Zones<-as.factor(ZONAS$Zones)
str(ZONAS)
structure(ZONAS)
# descriptive analysis # 
zonas_separadas<-unstack(ZONAS[,2:3])
#View(zonas_separadas)
mean(zonas_separadas[,1])
describe(zonas_separadas[,1])
#View(describe(zonas_separadas[,1]))
describe(zonas_separadas[,2])
#View(describe(zonas_separadas[,2]))
describe(zonas_separadas[,3])
#View(describe(zonas_separadas[,3]))
describe(zonas_separadas[,4])
#View(describe(zonas_separadas[,4]))
describe(zonas_separadas[,5])
#View(describe(zonas_separadas[,5]))

plot(Abundance~Zones, data=ZONAS,center=F, main = "Abundance by height zone") # to get an idea of how the data is in each zone

# testing the normality of the data # 

shapiro.test(ZONAS$Abundance )
# there is no normality of the data

leveneTest (Abundance ~ as.factor(Zones), data = ZONAS)

# there is no normality of the data

### very low p-value, that is, we reject HO and accept H1: there is no normality of the data and no homogeneity of the variances ###

hist(ZONAS$Abundance , main = "Histogram of the frequency distribution of the data", xlab = "Abundance", ylab="Frequency") # note that they are shifted to the left of the graph, they do not have a gaussian or normal distribution

# making an anova to remove the residues and test the normality of the residues

ANOVA_ZONAS<-aov(Abundance~Zones, data=ZONAS)
shapiro.test(resid(ANOVA_ZONAS)) 
hist(resid(ANOVA_ZONAS), main = "Histograma da distribuicao dos residuos", xlab = "Abundancia de especies", ylab="Frequencia") 
# note that it has improved, but it is still shifted to the left, it has no gaussian or normal distribution

# also there is no normality of residues

#############################
#############################
 # analysis of abundance in tree height zones - testing
#############################
#############################

glm_model.abd<-glm(Abundance ~ Zones, data = ZONAS, family = "poisson")
summary(glm_model.abd)
summ(glm_model.abd)
# there is a difference between the zones

# we need to do a post hoc test
# Tukey's post hoc test for generalized linear models was chosen with the multcomp package

summary(glht (glm_model.abd, mcp (Zones = "Tukey"))) #perform TukeyHSD

#############################
#############################
## let's do the GLMM #####
#############################
#############################

head(ZONAS)
glmm.zonas.abd<-glmer(Abundance ~ Zones + (1 | Forests), data = ZONAS, family = "poisson")
summary(glmm.zonas.abd)
# POST HOC
m_meansall <- emmeans(glmm.zonas.abd, ~ Zones, type = "response")
RESULTADO_ABD_GLMM<-pairs(m_meansall)
RESULTADO_ABD_GLMM

cld(m_meansall,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")
###################################################

# comparing the models
# comparison with null model

glm.nulo<-glm(Abundance ~ 1, data = ZONAS, family = "poisson")
summary(glm.nulo)
anova(glm_model.abd,glm.nulo)
AICctab(glm_model.abd, glm.nulo, base = T, weights = T)

# glm x glmm

AICctab(glmm.zonas.abd, glm_model.abd, base = T, weights = T)
# glmm is a better model than glm #
summ(glmm.zonas.abd)

# write.table(RESULTADO_ABD_GLMM,"RESULTADO_ABD_GLMM.csv",row.names=F,sep=";", dec = ",") # save table with results

## GENERAL GRAPHIC # 

pp<-ggplot(ZONAS, aes(x=Zones, y=Abundance, fill=Zones)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.fill = NULL,
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = NULL, )  + 
  # geom_point(position=position_jitter(width=0.3, height = 0.05, seed=1),shape=21,fill="grey80",size=0.6) + 
  scale_fill_brewer(palette="greens") + 
  labs(x="Zones",y="Abundance") + 
  theme_bw()
# pp<-pp +theme(axis.title.x=element_blank(),
#    axis.text.x=element_blank(),
#   axis.ticks.x=element_blank())
# esse ultima tira o nome dos elementos do eixo x, optei por deixar

pp

# to save the graphic # 

ggsave("Abundancia_Zonas.jpeg", units="in", dpi=600)

###########################
###########################

# graph with mean and standard deviation # 

Abundance <- c(1.24,2.04,1.22,1.07,0.50)
sd <- c(1.81,2.57,2.00,2.31,1.84)
Zones<-c(1,2,3,4,5)

dat_zones <- data.frame(Zones, Abundance, sd)
dat_zones$Zones<-as.factor(dat_zones$Zones)
str(dat_zones)

abd_media_desvio <- ggplot(dat_zones, aes(factor(Zones), Abundance)) + 
  geom_point(stat = "identity", size = 2) + 
  geom_errorbar(aes(ymin = Abundance - sd, ymax =  Abundance + sd), width = 0.2) +
  xlab("Zones")+  
  theme_apa() 
abd_media_desvio

ggsave("media_sd_riqueza.jpeg", units="in", dpi=600) # to save the graphic # 

##############################################################################################################################################

# GLM for each of the installments separately# 

##############################################################################################################################################

head(ZONAS)
# restinga 1 #

rest1<-ZONAS %>% filter (Forests == "Restinga1")
#View(rest1)
str(rest1)
zonas_separadas<-unstack(rest1[,2:3])
#View(zonas_separadas)
mean(zonas_separadas[,1])
describe(zonas_separadas[,1])
#View(describe(zonas_separadas[,1]))
describe(zonas_separadas[,2])
#View(describe(zonas_separadas[,2]))
describe(zonas_separadas[,3])
#View(describe(zonas_separadas[,3]))
describe(zonas_separadas[,4])
#View(describe(zonas_separadas[,4]))
describe(zonas_separadas[,5])
modelo.rest1<-glm(Abundance ~ Zones, data = rest1, family = "poisson")
summ(modelo.rest1, confint = TRUE, digits = 3)
summary(modelo.rest1)
summary(glht (modelo.rest1, mcp (Zones = "Tukey"))) 
# comparison with null model
glm.nulo<-glm(Abundance ~ 1, data = rest1, family = "poisson")
summary(glm.nulo)
anova(modelo.rest1,glm.nulo)
AICctab(modelo.rest1, glm.nulo, base = T, weights = T)

# restinga 2 #

rest2<-ZONAS %>% filter (Forests == "Restinga2")
#View(rest2)
str(rest2)
zonas_separadas<-unstack(rest2[,2:3])
#View(zonas_separadas)
mean(zonas_separadas[,1])
describe(zonas_separadas[,1])
#View(describe(zonas_separadas[,1]))
describe(zonas_separadas[,2])
#View(describe(zonas_separadas[,2]))
describe(zonas_separadas[,3])
#View(describe(zonas_separadas[,3]))
describe(zonas_separadas[,4])
#View(describe(zonas_separadas[,4]))
describe(zonas_separadas[,5])
modelo.rest2<-glm(Abundance ~ Zones, data = rest2, family = "poisson")
summ(modelo.rest2, confint = TRUE, digits = 3)
summary(glht (modelo.rest2, mcp (Zones = "Tukey"),digits = 3))
# comparison with null model
glm.nulo<-glm(Abundance ~ 1, data = rest2, family = "poisson")
summary(glm.nulo)
anova(modelo.rest2,glm.nulo)
AICctab(modelo.rest2, glm.nulo, base = T, weights = T)

# lowland forest # 

terras.baixas<-ZONAS %>% filter (Forests == "FlorestaTerrasBaixas")
str(terras.baixas)
zonas_separadas<-unstack(terras.baixas[,2:3])
#View(zonas_separadas)
mean(zonas_separadas[,1])
describe(zonas_separadas[,1])
#View(describe(zonas_separadas[,1]))
describe(zonas_separadas[,2])
#View(describe(zonas_separadas[,2]))
describe(zonas_separadas[,3])
#View(describe(zonas_separadas[,3]))
describe(zonas_separadas[,4])
#View(describe(zonas_separadas[,4]))
describe(zonas_separadas[,5])
modelo.terras.baixas<-glm(Abundance ~ Zones, data = terras.baixas, family = "poisson")
summ(modelo.terras.baixas, confint = TRUE, digits = 3)
summary(glht (modelo.terras.baixas, mcp (Zones = "Tukey"),digits = 3))
# comparison with null model
glm.nulo<-glm(Abundance ~ 1, data = terras.baixas, family = "poisson")
summary(glm.nulo)
anova(modelo.terras.baixas,glm.nulo)
AICctab(modelo.terras.baixas, glm.nulo, base = T, weights = T)

# submontane forest 1 #

submontana1<-ZONAS %>% filter (Forests == "FlorestaSubmontana1")
str(submontana1)
zonas_separadas<-unstack(submontana1[,2:3])
#View(zonas_separadas)
mean(zonas_separadas[,1])
describe(zonas_separadas[,1])
#View(describe(zonas_separadas[,1]))
describe(zonas_separadas[,2])
#View(describe(zonas_separadas[,2]))
describe(zonas_separadas[,3])
#View(describe(zonas_separadas[,3]))
describe(zonas_separadas[,4])
#View(describe(zonas_separadas[,4]))
describe(zonas_separadas[,5])
modelo.submontana1<-glm(Abundance ~ Zones, data = submontana1, family = "poisson")
summ(modelo.submontana1, confint = TRUE, digits = 3)
summary(glht (modelo.submontana1, mcp (Zones = "Tukey"),digits = 3))
# comparison with null model
glm.nulo<-glm(Abundance ~ 1, data = submontana1, family = "poisson")
summary(glm.nulo)
anova(modelo.submontana1,glm.nulo)
AICctab(modelo.submontana1, glm.nulo, base = T, weights = T)

# submontane forest 2 #
submontana2<-ZONAS %>% filter (Forests == "Florestasubmontana2")
#View(submontana2)
zonas_separadas<-unstack(submontana2[,2:3])
#View(zonas_separadas)
mean(zonas_separadas[,1])
describe(zonas_separadas[,1])
#View(describe(zonas_separadas[,1]))
describe(zonas_separadas[,2])
#View(describe(zonas_separadas[,2]))
describe(zonas_separadas[,3])
#View(describe(zonas_separadas[,3]))
describe(zonas_separadas[,4])
#View(describe(zonas_separadas[,4]))
describe(zonas_separadas[,5])
modelo.submontana2<-glm(Abundance ~ Zones, data = submontana2, family = "poisson")
summ(modelo.submontana2, confint = TRUE, digits = 3)
summary(glht (modelo.submontana2, mcp (Zones = "Tukey"),digits = 3))
# comparison with null model
glm.nulo<-glm(Abundance ~ 1, data = submontana2, family = "poisson")
summary(glm.nulo)
anova(modelo.submontana2,glm.nulo)
AICctab(modelo.submontana2, glm.nulo, base = T, weights = T)

#######################################################
#######################################################
# GRAPHICS of the forests separately # 
#######################################################
#######################################################

Parc<-c(FlorestaSubmontana2="Submontana 2", FlorestaTerrasBaixas="Terras Baixas",FlorestaSubmontana1="Submontana 1",Restinga1="Restinga 1", Restinga2="Restinga 2")

ggplot(ZONAS, aes(x=Zones, y=Abundance, fill=Zones))+
  geom_boxplot(outlier.shape=NA) + # remove the outlier
  geom_point(position=position_jitter(width=0.3, height = 0.05, seed=1),shape=21,fill="grey80",size=0.6)+ scale_fill_brewer(palette="greens")+ 
  facet_wrap( ~ Forests, labeller = labeller(Forests = Parc)) + 
  labs(x="Zones",y="Abundance") + 
  theme_base() +
  theme(axis.title.x=element_blank(), # removing the x-axis labels from each graph individually
        axis.text.x=element_blank(), # remove text from x axis
        axis.ticks.x=element_blank()) # remove the dashes from the axis


###########################################################################################################################################################

