
########### 1.  Calculation proportion of indirect links per each pollinator and seed disperser species

########### 2.  Comparison of proportion of indirect links per pollinator and seed disperser species between treatments

########### 3. Estimation of dissimilarity of indirect links between treatments

#-----------------------------------------------------------------------------------------------------------------
########### 1.  Calculation proportion of indirect links per each pollinator and seed disperser species

library(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      NI TREATMENT                            
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Proportion of indirect links per pollinator sp--

#upload matrix
NI<- read.csv("Paths_NI.csv", sep=",", row.names = 1)
NI_disp <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/NI_disp.csv", sep=",",row.names=1)
tot_disp = nrow(NI_disp) #denominator of the proportion

NI_pol_ind<-NI %>% #calculate the number and proportion of indirect links
  mutate(Num_int = rowSums(NI, na.rm=TRUE),
         Prop_int = Num_int / tot_disp) %>% 
  select(Num_int,Prop_int)

NI_pol_ind <- cbind(Treat = rep("NI",nrow(NI_pol_ind)), Species = rownames(NI_pol_ind),NI_pol_ind)
rownames(NI_pol_ind) = NULL


#Proportion of indirect links per seed disperser sp--

NI <- t(NI) # transpose the matrix to consider the seed disperser as rows
NI<-as.data.frame(NI)

NI_pol <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/NI_pol.csv", sep=",",row.names=1)
tot_pol = nrow(NI_pol) #denominator of the proportion 

NI_seeddisp_ind<-NI %>% ##calculate the number and proportion of indirect links
  mutate(Num_int = rowSums(NI, na.rm=TRUE),
         Prop_int = Num_int / tot_pol) %>% 
  select(Num_int,Prop_int)

NI_seeddisp_ind <- cbind(Treat = rep("NI",nrow(NI_seeddisp_ind)), Species = rownames(NI_seeddisp_ind),NI_seeddisp_ind)
rownames(NI_seeddisp_ind) = NULL

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      I TREATMENT                            
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Proportion of indirect links per pollinator sp--

#upload matrix
I<- read.csv("Paths_I.csv", sep=",", row.names = 1)
I_disp <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/I_disp.csv", sep=",",row.names=1)
tot_disp = nrow(I_disp)#denominator of the proportion

I_pol_ind<-I %>% #calculate the number and proportion of indirect links
  mutate(Num_int = rowSums(I, na.rm=TRUE),
         Prop_int = Num_int / tot_disp) %>% 
  select(Num_int,Prop_int)

I_pol_ind <- cbind(Treat = rep("I",nrow(I_pol_ind)),Species = rownames(I_pol_ind),I_pol_ind)
rownames(I_pol_ind) = NULL

#Calculate the proportion of indirect links for seed disperser species--

I <- t(I) # transpose the matrix to consider the seed disperser as rows
I<-as.data.frame(I)

I_pol <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/I_pol.csv", sep=",",row.names=1)
tot_pol = nrow(I_pol) #denominator of the proportion

I_seeddisp_ind<-I %>% #calculate the number and proportion of indirect links
  mutate(Num_int = rowSums(I, na.rm=TRUE),
         Prop_int = Num_int / tot_pol) %>% 
  select(Num_int,Prop_int)

I_seeddisp_ind <- cbind(Treat = rep("I",nrow(I_seeddisp_ind)), Species = rownames(I_seeddisp_ind),I_seeddisp_ind)
rownames(I_seeddisp_ind) = NULL

#######################################################
####### Create data frame and summarize results #######
#######################################################

Pol_indirect <-rbind(NI_pol_ind,I_pol_ind) #indirect links per pollinator species for both treatments

Disp_indirect <-rbind(NI_seeddisp_ind,I_seeddisp_ind)#indirect links per seed disperser species for both treatments

#-----------------------------------------------------------------------------------------------------------------
########### 2.  Comparison of proportion of indirect links per pollinator and seed disperser species between treatments

library (vegan)
library(ggplot2)
library(lme4)
library(glmmTMB)
library("stats4")
library("bbmle")


##Proportion of indirect links per pollinator species--

#Prepare dataframe
Pol_indirect_new <- Pol_indirect %>% 
  mutate(Prop_int=replace(Prop_int, Prop_int==1.00, 0.9999999)) %>%
  mutate(Prop_int=replace(Prop_int, Prop_int==0.00, 0.0000001))

m1<-glmmTMB(Prop_int~ Treat + ( 1| Species),beta_family(link = "logit"), data = Pol_indirect_new)
summary(m1) 

#Homogeneity 
EM<-resid(m1, type= "pearson") 
FM<-fitted(m1)
plot(x=FM, y=EM, xlab = "Ajustados", ylab = "Residuales normalizados")
abline(0,0, col="red", lwd= 3) 

#Independence
boxplot(EM~Pol_indirect_new$Treat, main="Tratamiento")


##Proportion of indirect links per seed disperser species--
Disp_indirect_new <- Disp_indirect %>% 
  mutate(Prop_int=replace(Prop_int, Prop_int==1.00, 0.9999999)) %>%
  mutate(Prop_int=replace(Prop_int, Prop_int==0.00, 0.0000001)) 

m2<-glmmTMB(Prop_int~ Treat + ( 1| Species),beta_family(link = "logit"), data = Disp_indirect_new)
summary(m4) 

#Homogeneity 
EM<-resid(m2, type= "pearson") 
FM<-fitted(m2)
plot(x=FM, y=EM, xlab = "Ajustados", ylab = "Residuales normalizados")
abline(0,0, col="red", lwd= 3) 

#Independence 
boxplot(EM~Disp_indirect_new$Treat, main="Tratamiento")

#-----------------------------------------------------------------------------------------------------------------
########### 3. Estimation of dissimilarity of indirect links between treatments

#Here we calculated the dissimilarity of indirect links identities between treatments using the Jaccard distance. Jaccard Distance values range from 0 to 1. Higher values indicate less similarity between samples.

library(vegan)

#upload the matrix
NI<-read.csv("Paths_diss_NI.csv", row.names = 1)#matrix containing presence/absence of indirect links between pollinator and seed disp after preparing the dataframe to calculate dissimilarity with the other treatment (same rows and columns)
NI<-as.matrix(NI)
I<-read.csv("Paths_diss_I.csv", row.names = 1) #matrix containing presence/absence of indirect links between pollinator and seed disp after preparing the dataframe to calculate dissimilarity with the other treatment (same rows and columns)
I<-as.matrix(I)

#array
vNI<-array(as.matrix(NI))
vI<-array(as.matrix(I))

Intpol_disp<-rbind(vNI,vI)
Intpol_dispst<-sqrt(Intpol_disp) #standarize

DisIntpol_disp<- vegdist(Intpol_dispst, binary = TRUE, method = "jaccard") #calculate the distance matrix using jaccard index

DisIntpol_disp 

