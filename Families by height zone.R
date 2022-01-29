
###############################################
#
# Date: 23-August-2021
# Name: Analysis of the distribution of botanical families of epiphytes along the arboreal height zones
# Author: ANA CAROLINA RODRIGUES DA CRUZ
# E-mail: anacruzufrj@gmail.com
# 

# Description: comparison of the richness and abundance of the main epiphyte families along the tree height zones of the forests of Ilha Grande
#
################################################
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
library(tidyr) # data manipulation
library(ggplot2) #graphs
library(ggthemes)# graph of GGPLOT
library(RColorBrewer) # graph of GGPLOT
library(bbmle) # use the AICctab function from the bbmle package
library(jtools)# more detailed output of models output

##############################################

# for GLM 

library(lme4) # to do GLM
library(multcomp) # for glm tukey post hoc test 
library(emmeans) # ditto previous

##########################################

##### load directory
setwd("D:/Arquivos/Documents/ANA 2021/TESE/Tese 2021/Manuscritos/manuscrito 1/Scripts - repositorio")

########## ORICHIDACEAE LINE 53 #####
########## BROMELIACEAE LINE 203 ######
########## POLYPODIACEAE LINE 335 #####
########## ARACEAE LINE 463 #####
########## CACTACEAE LINE 589 #####

######################################################################

ZONAS<-read.table("orchidaceae_por_zona.csv",header=T,sep=";")
ZONAS$Zonas<-as.factor(ZONAS$Zonas)
#View(ZONAS)
str(ZONAS)
#structure(ZONAS)

# isolating the zones to see the richness of each one of them # #

zona1<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "1")
describe(zona1$Riqueza)
zona2<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "2")
describe(zona2$Riqueza)
zona3<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "3")
describe(zona3$Riqueza)
zona4<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "4")
describe(zona4$Riqueza)
zona5<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "5")
describe(zona5$Riqueza)
plot(Riqueza~Zonas, data=ZONAS,center=F, main = "Riqueza por zona de altura") # just to get an idea of how the data is in each zone

########################################
##### richness variation by zone  #####

# testing the normality of the data # 

shapiro.test(ZONAS$Riqueza) #there is no normality of the data
leveneTest (Riqueza ~ as.factor(Zonas), data = ZONAS) #there is no homoscedasticity of the data
hist(ZONAS$Riqueza, main = "Histogram of the frequency distribution of the data", xlab = "Richness", ylab="Frequency") #n ote that they are shifted to the left of the graph, they do not have a gaussian or normal distribution
# making an anova to remove the residuals and test the normality of the residuals
ANOVA_ZONAS<-aov(Riqueza~Zonas, data=ZONAS)
shapiro.test(resid(ANOVA_ZONAS)) 
hist(resid(ANOVA_ZONAS),  main = "Histogram of the frequency distribution of the data", xlab = "Richness", ylab="Frequency") # note that it has improved, but they are still shifted to the left, no gaussian or normal distribution
# there is also no normality of the residues
## GRAPHIC

pp<-ggplot(ZONAS, aes(x=Zonas, y=Riqueza, fill=Zonas)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitter(width=0.3, height = 0.05, seed=1),shape=21,fill="grey80",size=0.6) + 
  scale_fill_brewer(palette="greens") + 
  labs(x="Zones",y="Richness") + 
  theme_bw()
pp

############ testing with GLM ###########

glm_model<-glm(Riqueza ~ Zonas, data = ZONAS, family = "poisson")

summary(glm_model)
summary(glht (glm_model, mcp (Zonas = "Tukey"))) #perform TukeyHSD

# now mixed model, considering the plots as a random effect #
head(ZONAS)
glmm.zonas<-glmer(Riqueza ~ Zonas + (1 | Parcela), data = ZONAS, family = "poisson")
summary(glmm.zonas)
# POST HOC
m_meansall <- emmeans(glmm.zonas, ~ Zonas, type = "response")
pairs(m_meansall)
cld(m_meansall,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")

###################################################

# comparing the models

# comparison with null model
glm.nulo<-glm(Riqueza ~ 1, data = ZONAS, family = "poisson")
summary(glm.nulo)
anova(glm_model,glm.nulo)
anova(glm_model, glm.nulo, refit = FALSE, test = "Chisq")

AICctab(glm_model, glm.nulo, base = T, weights = T)
# glm x glmm

AICctab(glmm.zonas, glm_model, base = T, weights = F)
summ(glmm.zonas)
anova (glmm.zonas, glm_model)
# Glmm is a better model than glm 
# comparison with null model
glmm.nulo<-glmer(Riqueza ~ 1 + (1 | Parcela), data = ZONAS, family = "poisson")
anova(glmm.zonas,glmm.nulo)

# there is a difference in the richness of orchids along the tree height zones #
# zone 5 in relation to others #

#################################################################################################################################################################################################                 
# now let's test the abundance #
################################################################################################################################################################################################# 

head(ZONAS)

# isolating the zones to see the abundance of each one #
zona1<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "1")
describe(zona1$Abundancia)
zona2<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "2")
describe(zona2$Abundancia)
zona3<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "3")
describe(zona3$Abundancia)
zona4<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "4")
describe(zona4$Abundancia)
zona5<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "5")
describe(zona5$Abundancia)

#####################################

shapiro.test(ZONAS$Abundancia )
# #there is no normality of the data
leveneTest (Abundancia  ~ as.factor(Zonas), data = ZONAS)
# THERE IS NO HOMOCEDASTICITY OF THE DATA #
# GRAPHIC# 

pp<-ggplot(ZONAS, aes(x=Zonas, y=Abundancia, fill=Zonas)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitter(width=0.3, height = 0.05, seed=1),shape=21,fill="grey80",size=0.6) + 
  scale_fill_brewer(palette="greens") + 
  labs(x="Zones",y="Abundance") + 
  theme_bw()
pp
# glm # 
glm_model_abd<-glm(Abundancia ~ Zonas, data = ZONAS, family = "poisson")
summary(glm_model_abd)
summary(glht (glm_model_abd, mcp (Zonas = "Tukey"))) #perform TukeyHSD

# now mixed model, considering the plots as a random effect #
head(ZONAS)
glmm.zonas.abd<-glmer(Abundancia ~ Zonas + (1 | Parcela), data = ZONAS, family = "poisson")
summary(glmm.zonas.abd)
summ(glmm.zonas.abd)
# POST HOC
m_meansall <- emmeans(glmm.zonas.abd, ~ Zonas, type = "response")
pairs(m_meansall)
cld(m_meansall,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")
###################################################

# comparing the models

# comparison with null model
glm.nulo<-glm(Abundancia ~ 1, data = ZONAS, family = "poisson")
summary(glm.nulo)
anova(glm_model_abd, glm.nulo, refit = FALSE, test = "Chisq")
AICctab(glm_model_abd, glm.nulo, base = T, weights = T)
# glm x glmm
AICctab(glmm.zonas.abd, glm_model_abd, base = T, weights = T)
summ(glmm.zonas.abd)
# glmm is a better model than glm

############################################################################################################################################
# BROMELIACEAE #
###########################################################################################################################################

ZONAS<-read.table("Bromeliaceae_por_zona.csv",header=T,sep=";")
ZONAS$Zonas<-as.factor(ZONAS$Zonas)
str(ZONAS)
# isolating the zones to see the richness of each one of them #
zona1<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "1")
describe(zona1$Riqueza)
zona2<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "2")
describe(zona2$Riqueza)
zona3<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "3")
describe(zona3$Riqueza)
zona4<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "4")
describe(zona4$Riqueza)
zona5<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "5")
describe(zona5$Riqueza)
#####################################33
plot(Riqueza~Zonas, data=ZONAS,center=F, main = "Riqueza por zona de altura") 
########################################
##### Richness variation by zone #####

# testing data for normality #

shapiro.test(ZONAS$Riqueza)
# there is no normality of the data
leveneTest (Riqueza ~ as.factor(Zonas), data = ZONAS)
#there is no homoscedasticity
hist(ZONAS$Riqueza, main = "Histogram of the frequency distribution of the data", xlab = "Richness", ylab="Frequency") # note that it has improved, but they are still shifted to the left, no gaussian or normal distribution
ANOVA_ZONAS<-aov(Riqueza~Zonas, data=ZONAS)
shapiro.test(resid(ANOVA_ZONAS)) 
hist(resid(ANOVA_ZONAS), main = "Histograma da distribuicao dos residuos", xlab = "Riqueza de especies", ylab="Frequencia")
# there is no normality of the data
## Graphics # 
pp<-ggplot(ZONAS, aes(x=Zonas, y=Riqueza, fill=Zonas)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitter(width=0.3, height = 0.05, seed=1),shape=21,fill="grey80",size=0.6) + 
  scale_fill_brewer(palette="greens") + 
  labs(x="Zones",y="Richness") + 
  theme_bw()
pp

 # GLM
glm_model<-glm(Riqueza ~ Zonas, data = ZONAS, family = "poisson")
summary(glm_model)
summary(glht (glm_model, mcp (Zonas = "Tukey"))) #perform TukeyHSD
# now mixed model, considering the plots as a random effecthead(ZONAS)
glmm.zonas<-glmer(Riqueza ~ Zonas + (1 | Parcela), data = ZONAS, family = "poisson")
summary(glmm.zonas)
# POST HOC
m_meansall <- emmeans(glmm.zonas, ~ Zonas, type = "response")
pairs(m_meansall)
cld(m_meansall,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")
###################################################

# comparing the models
# glm x glmm
AICctab(glmm.zonas, glm_model, base = T, weights = T)
# glmm is a better model than glm #
summary(glmm.zonas)
summ(glmm.zonas)

# comparison with null model

glm.nulo<-glm(Riqueza ~ 1, data = ZONAS, family = "poisson")
summary(glm.nulo)
anova(glm_model, glm.nulo, refit = FALSE, test = "Chisq")
AICctab(glm_model, glm.nulo, base = T, weights = T)
# null model comparison to extract p-value
glmm.nulo<-glmer(Riqueza ~ 1 + (1 | Parcela), data = ZONAS, family = "poisson")
anova(glmm.zonas,glm_model)
#################################################################################################################################################################################################                 
# now let's test the abundance #
################################################################################################################################################################################################# 

head(ZONAS)
zona1<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "1")
describe(zona1$Abundancia)
zona2<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "2")
describe(zona2$Abundancia)
zona3<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "3")
describe(zona3$Abundancia)
zona4<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "4")
describe(zona4$Abundancia)
zona5<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "5")
describe(zona5$Abundancia)
#####################################
shapiro.test(ZONAS$Abundancia )
# there is no data normality
leveneTest (Abundancia  ~ as.factor(Zonas), data = ZONAS)
# there is no homoscedasticity
## Graphic 
pp<-ggplot(ZONAS, aes(x=Zonas, y=Abundancia, fill=Zonas)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitter(width=0.3, height = 0.05, seed=1),shape=21,fill="grey80",size=0.6) + 
  scale_fill_brewer(palette="greens") + 
  labs(x="Zones",y="Abundance") + 
  theme_bw()
pp
# glm # 
glm_model_abd<-glm(Abundancia ~ Zonas, data = ZONAS, family = "poisson")
summary(glm_model_abd)
summary(glht (glm_model_abd, mcp (Zonas = "Tukey"))) #perform TukeyHSD

# now mixed model, considering the plots as a random effect#
head(ZONAS)
glmm.zonas.abd<-glmer(Abundancia ~ Zonas + (1 | Parcela), data = ZONAS, family = "poisson")
summary(glmm.zonas.abd)
summ(glmm.zonas.abd)
# POST HOC
m_meansall <- emmeans(glmm.zonas.abd, ~ Zonas, type = "response")
pairs(m_meansall)
cld(m_meansall,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")
###################################################
# comparing the models

# comparison with null model
glm.nulo<-glm(Abundancia ~ 1, data = ZONAS, family = "poisson")
summary(glm.nulo)
anova(glm_model_abd, glm.nulo, refit = FALSE, test = "Chisq")
# glm x glmm
AICctab(glmm.zonas.abd, glm_model_abd, base = T, weights = T)
summ(glmm.zonas.abd)
# glmm is a better model than glm #

############################################################################################################################################
# POLYPODIACEAE #
############################################################################################################################################

ZONAS<-read.table("polypodiaceae_por_zona.csv",header=T,sep=";")
#View(ZONAS)
ZONAS$Zonas<-as.factor(ZONAS$Zonas)
str(ZONAS)
#structure(ZONAS)
zona1<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "1")
describe(zona1$Riqueza)
zona2<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "2")
describe(zona2$Riqueza)
zona3<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "3")
describe(zona3$Riqueza)
zona4<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "4")
describe(zona4$Riqueza)
zona5<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "5")
describe(zona5$Riqueza)
#####################################33
plot(Riqueza~Zonas, data=ZONAS,center=F, main = "Riqueza por zona de altura")

##### Richness variation by zone #####

# testing data for normality #
shapiro.test(ZONAS$Riqueza)
# there is no data normality
leveneTest (Riqueza ~ as.factor(Zonas), data = ZONAS)
# there is no homoscedasticity
hist(ZONAS$Riqueza, main = "", xlab = "Richness", ylab="Frequency") 
ANOVA_ZONAS<-aov(Riqueza~Zonas, data=ZONAS)
shapiro.test(resid(ANOVA_ZONAS)) 
hist(resid(ANOVA_ZONAS), main = "", xlab = "", ylab="Frequency") # there is no data normality

### Graphic # 

pp<-ggplot(ZONAS, aes(x=Zonas, y=Riqueza, fill=Zonas)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitter(width=0.3, height = 0.05, seed=1),shape=21,fill="grey80",size=0.6) + 
  scale_fill_brewer(palette="greens") + 
  labs(x="Zones",y="Richness") + 
  theme_bw()
pp

glm_model<-glm(Riqueza ~ Zonas, data = ZONAS, family = "poisson")
summary(glm_model)
summary(glht (glm_model, mcp (Zonas = "Tukey"))) #perform TukeyHSD
# now mixed model, considering the plots as a random effect #
head(ZONAS)
glmm.zonas<-glmer(Riqueza ~ Zonas + (1 | Parcela), data = ZONAS, family = "poisson")
summary(glmm.zonas)
# POST HOC
m_meansall <- emmeans(glmm.zonas, ~ Zonas, type = "response")
pairs(m_meansall)
cld(m_meansall,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")

summ(glmm.zonas)
###################################################
# comparing the models
# comparison with null model
glm.nulo<-glm(Riqueza ~ 1, data = ZONAS, family = "poisson")
summary(glm.nulo)
anova(glm_model, glm.nulo, refit = FALSE, test = "Chisq")
# glm x glmm
AICctab(glmm.zonas, glm_model, base = T, weights = T)
# glmm is a better model than glm #
# there is a difference in the richness of polypodiaceae along the tree height zones

#################################################################################################################################################################################################                 
# agora vamos testar a abundancia #
#################################################################################################################################################################################################                 
head(ZONAS)
zona1<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "1")
describe(zona1$Abundancia)
zona2<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "2")
describe(zona2$Abundancia)
zona3<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "3")
describe(zona3$Abundancia)
zona4<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "4")
describe(zona4$Abundancia)
zona5<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "5")
describe(zona5$Abundancia)
#####################################
shapiro.test(ZONAS$Abundancia )
# there is no data normality
leveneTest (Abundancia  ~ as.factor(Zonas), data = ZONAS)
# there is no homogeneity#
# Graphic # 

pp<-ggplot(ZONAS, aes(x=Zonas, y=Abundancia, fill=Zonas)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitter(width=0.3, height = 0.05, seed=1),shape=21,fill="grey80",size=0.6) + 
  scale_fill_brewer(palette="greens") + 
  labs(x="Zones",y="Abundance") + 
  theme_bw()
pp

# glm # 
glm_model_abd<-glm(Abundancia ~ Zonas, data = ZONAS, family = "poisson")
summary(glm_model_abd)
summary(glht (glm_model_abd, mcp (Zonas = "Tukey"))) #perform TukeyHSD
# now mixed model, considering the plots as a random effect #
head(ZONAS)
glmm.zonas.abd<-glmer(Abundancia ~ Zonas + (1 | Parcela), data = ZONAS, family = "poisson")
summary(glmm.zonas.abd)
summ(glmm.zonas.abd)
# POST HOC
m_meansall <- emmeans(glmm.zonas.abd, ~ Zonas, type = "response")
pairs(m_meansall)
cld(m_meansall,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")
###################################################
# comparing the models
# comparison with null model
glm.nulo<-glm(Abundancia ~ 1, data = ZONAS, family = "poisson")
summary(glm.nulo)
anova(glm_model_abd, glm.nulo, refit = FALSE, test = "Chisq")
# glm x glmm
AICctab(glmm.zonas.abd, glm_model_abd, base = T, weights = T)
summ(glmm.zonas.abd)
# glmm is a better model than glm #

############################################################################################################################################
# Araceae #
############################################################################################################################################
ZONAS<-read.table("araceae_por_zona.csv",header=T,sep=";")
ZONAS$Zonas<-as.factor(ZONAS$Zonas)
str(ZONAS)
#structure(ZONAS)
zona1<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "1")
describe(zona1$Riqueza)
zona2<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "2")
describe(zona2$Riqueza)
zona3<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "3")
describe(zona3$Riqueza)
zona4<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "4")
describe(zona4$Riqueza)
zona5<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "5")
describe(zona5$Riqueza)
####################################33
plot(Riqueza~Zonas, data=ZONAS,center=F, main = "Richness") 
## Richness variation by zone #####
# testing data for normality #
shapiro.test(ZONAS$Riqueza)
# there is no data normality
leveneTest (Riqueza ~ as.factor(Zonas), data = ZONAS)
# there is no homoscedasticity
hist(ZONAS$Riqueza, main = "", xlab = "Richness", ylab="Frequency")
ANOVA_ZONAS<-aov(Riqueza~Zonas, data=ZONAS)
shapiro.test(resid(ANOVA_ZONAS)) 
hist(resid(ANOVA_ZONAS), main = "", xlab = "", ylab="Frequency") 

# also there is no normality#
# Graphics # 

pp<-ggplot(ZONAS, aes(x=Zonas, y=Riqueza, fill=Zonas)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitter(width=0.3, height = 0.05, seed=1),shape=21,fill="grey80",size=0.6) + 
  scale_fill_brewer(palette="greens") + 
  labs(x="Zones",y="Richness") + 
  theme_bw()
pp

glm_model<-glm(Riqueza ~ Zonas, data = ZONAS, family = "poisson")
summary(glm_model)
summary(glht (glm_model, mcp (Zonas = "Tukey"))) #perform TukeyHSD
# now mixed model, considering the plots as a random effect #
head(ZONAS)
glmm.zonas<-glmer(Riqueza ~ Zonas + (1 | Parcela), data = ZONAS, family = "poisson")
summary(glmm.zonas)
summ(glmm.zonas)
# POST HOC
m_meansall <- emmeans(glmm.zonas, ~ Zonas, type = "response")
pairs(m_meansall)
cld(m_meansall,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")
###################################################
# comparing the models
# comparison with null model
glm.nulo<-glm(Riqueza ~ 1, data = ZONAS, family = "poisson")
summary(glm.nulo)
anova(glm_model, glm.nulo, refit = FALSE, test = "Chisq")
# glm x glmm
AICctab(glmm.zonas, glm_model, base = T, weights = T)
summ(glmm.zonas)
# glmm is a better model than glm #

# there is a difference in the richness of polypodiaceae along the tree height zones #

############################################################################################################################################
# now let's test the abundance ############################################################################################################################################
head(ZONAS)
zona1<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "1")
describe(zona1$Abundancia)
zona2<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "2")
describe(zona2$Abundancia)
zona3<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "3")
describe(zona3$Abundancia)
zona4<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "4")
describe(zona4$Abundancia)
zona5<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "5")
describe(zona5$Abundancia)
#####################################
shapiro.test(ZONAS$Abundancia )
# there is no data normality
leveneTest (Abundancia  ~ as.factor(Zonas), data = ZONAS)
# there is no homogeneity
## Graphics ##
pp<-ggplot(ZONAS, aes(x=Zonas, y=Abundancia, fill=Zonas)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitter(width=0.3, height = 0.05, seed=1),shape=21,fill="grey80",size=0.6) + 
  scale_fill_brewer(palette="greens") + 
  labs(x="Zones",y="Abundance") + 
  theme_bw()
pp

# glm # 
glm_model_abd<-glm(Abundancia ~ Zonas, data = ZONAS, family = "poisson")
summary(glm_model_abd)
summary(glht (glm_model_abd, mcp (Zonas = "Tukey"))) #perform TukeyHSD
# now mixed model, considering the plots as a random effect #
head(ZONAS)
glmm.zonas.abd<-glmer(Abundancia ~ Zonas + (1 | Parcela), data = ZONAS, family = "poisson")
summary(glmm.zonas.abd)
summ(glmm.zonas.abd)
# POST HOC
m_meansall <- emmeans(glmm.zonas.abd, ~ Zonas, type = "response")
pairs(m_meansall)
cld(m_meansall,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")
###################################################
# comparing the models

# comparison with null model
glm.nulo<-glm(Abundancia ~ 1, data = ZONAS, family = "poisson")
summary(glm.nulo)
anova(glm_model_abd, glm.nulo, refit = FALSE, test = "Chisq")
# glm x glmm

AICctab(glmm.zonas.abd, glm_model_abd, base = T, weights = T)
# glmm is a better model than glm #
summ(glmm.zonas.abd)

############################################################################################################################################
# cactaceae #
############################################################################################################################################

ZONAS<-read.table("cactaceae_por_zona.csv",header=T,sep=";")
ZONAS$Zonas<-as.factor(ZONAS$Zonas)
str(ZONAS)
#structure(ZONAS)
View(ZONAS)
zona1<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "1")
describe(zona1$Riqueza)
zona2<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "2")
describe(zona2$Riqueza)
zona3<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "3")
describe(zona3$Riqueza)
zona4<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "4")
describe(zona4$Riqueza)
zona5<-ZONAS%>%select(Riqueza, Zonas) %>% filter(Zonas == "5")
describe(zona5$Riqueza)
########################################
##### Richness variation by zone #####
# testing data for normality #
shapiro.test(ZONAS$Riqueza)
# there is no normality of the data
leveneTest (Riqueza ~ as.factor(Zonas), data = ZONAS)
# there is no homoscedasticity
### very low p-value, that is, we reject HO and accept H1: there is no normality of the data and no homogeneity of the variances
hist(ZONAS$Riqueza, main = "", xlab = "Richness", ylab="Frequency")
ANOVA_ZONAS<-aov(Riqueza~Zonas, data=ZONAS)
shapiro.test(resid(ANOVA_ZONAS)) 
hist(resid(ANOVA_ZONAS), main = "", xlab = "Riqueza de especies", ylab="Frequencia") # note that it has improved, but it is still shifted to the left, it has no gaussian or normal distribution
# also there is no normality of residuals
## Graphics# 

pp<-ggplot(ZONAS, aes(x=Zonas, y=Riqueza, fill=Zonas)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitter(width=0.3, height = 0.05, seed=1),shape=21,fill="grey80",size=0.6) + 
  scale_fill_brewer(palette="greens") + 
  labs(x="Zones",y="Richness") + 
  theme_bw()
pp

glm_model<-glm(Riqueza ~ Zonas, data = ZONAS, family = "poisson")
summary(glm_model)
summary(glht (glm_model, mcp (Zonas = "Tukey"))) #perform TukeyHSD
# now mixed model, considering the plots as a random effect #
head(ZONAS)
glmm.zonas<-glmer(Riqueza ~ Zonas + (1 | Parcela), data = ZONAS, family = "poisson")
summary(glmm.zonas)
summ(glmm.zonas)
# POST HOC
m_meansall <- emmeans(glmm.zonas, ~ Zonas, type = "response")
pairs(m_meansall)
cld(m_meansall,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")
###################################################

# comparing the models
# comparison with null model
glm.nulo<-glm(Riqueza ~ 1, data = ZONAS, family = "poisson")
summary(glm.nulo)
anova(glm_model, glm.nulo, refit = FALSE, test = "Chisq")
# glm x glmm
AICctab(glmm.zonas, glm_model, base = T, weights = T)
# glmm is a better model than glm #
summ(glmm.zonas)
# there is a difference in cactaceae richness along the tree height zones #

#################################################################################################################################################################################################        
# now let's test the abundance # ##################################################################################################################################################################################################        
head(ZONAS)
zona1<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "1")
describe(zona1$Abundancia)
zona2<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "2")
describe(zona2$Abundancia)
zona3<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "3")
describe(zona3$Abundancia)
zona4<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "4")
describe(zona4$Abundancia)
zona5<-ZONAS%>%select(Abundancia, Zonas) %>% filter(Zonas == "5")
describe(zona5$Abundancia)
#####################################
shapiro.test(ZONAS$Abundancia )
# there is no data normality
leveneTest (Abundancia  ~ as.factor(Zonas), data = ZONAS)
# there is no homogeneity
## Graphics # 

pp<-ggplot(ZONAS, aes(x=Zonas, y=Abundancia, fill=Zonas)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitter(width=0.3, height = 0.05, seed=1),shape=21,fill="grey80",size=0.6) + 
  scale_fill_brewer(palette="greens") + 
  labs(x="Zones",y="Abundance") + 
  theme_bw()
pp
# glm # 
glm_model_abd<-glm(Abundancia ~ Zonas, data = ZONAS, family = "poisson")
summary(glm_model_abd)
summary(glht (glm_model_abd, mcp (Zonas = "Tukey"))) #perform TukeyHSD
# now mixed model, considering the plots as a random effect #
head(ZONAS)
glmm.zonas.abd<-glmer(Abundancia ~ Zonas + (1 | Parcela), data = ZONAS, family = "poisson")
summary(glmm.zonas.abd)
summ(glmm.zonas.abd)
# POST HOC
m_meansall <- emmeans(glmm.zonas.abd, ~ Zonas, type = "response")
pairs(m_meansall)
cld(m_meansall,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")
###################################################
# comparing the models

# comparison with null model
glm.nulo<-glm(Riqueza ~ 1, data = ZONAS, family = "poisson")
summary(glm.nulo)
anova(glm_model_abd, glm.nulo, refit = FALSE, test = "Chisq")
# glm x glmm
AICctab(glmm.zonas.abd, glm_model_abd, base = T, weights = T)
# glmm is a better model than glm #
summ(glmm.zonas.abd)
######################################################