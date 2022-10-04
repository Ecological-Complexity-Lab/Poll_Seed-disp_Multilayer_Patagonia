########### 1. Estimation of the propagation after a particular disturbance event for each treatment

########### 2. Comparison between treatments
-----------------------------------------------------------------------------------------------------------------
  
###### 1. Estimation of the propagation after a particular disturbance event for each treatment --
library(tidyverse)
library(dplyr)
source('/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Extra_Functions/Extrafunctions.R')

#Upload matrix
NI_disp <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/NI_disp.csv", sep=",",row.names=1)
NI_pol  <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/NI_pol.csv", sep=",",row.names=1)
I_disp<- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/I_disp.csv", sep=",",row.names=1)
I_pol <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/I_pol.csv", sep=";",row.names=1)

# Load Ri values 
Rdisp_NI <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/Rdisp_NI.csv")
Rdisp_I <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/Rdisp_I.csv")
Rpol_NI <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/Rpol_NI.csv")
Rpol_I <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/Rpol_I.csv")
Rplant_pol_NI <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/Rplant_pol_NI.csv")
Rplant_pol_I <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/Rplant_pol_I.csv")
Rplant_disp_NI <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/Rplant_disp_NI.csv")
Rplant_disp_I <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/Rplant_disp_I.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      NI TREATMENT                            
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Ndisp=NI_disp
  Npol=NI_pol
  Rdisp=Rdisp_NI
  Rpol=Rpol_NI
  rplants_pol=Rplant_pol_NI
  rplants_disp=Rplant_disp_NI

Npol  <- Npol [ order(row.names(Npol)), ] #order pollinators

#Ri values for each species and trophic group
rdisp=Rdisp$x
rpol=Rpol$x
rplants_pol=rplants_pol$x
rplants_disp=rplants_disp$x

ndisp=dim(Ndisp)[1]
npol=dim(Npol)[1]
nplants=dim(Ndisp)[2]
nanim=npol+ndisp
ntot=nplants+npol+ndisp #Total number of nodes

# create final matrix
pol_indx=1:npol
disp_indx=(npol+1):(ndisp+npol)

imatrix=rbind(Npol,Ndisp)
IM= as.matrix(imatrix)

####################
####  Scenarios ####
####################

# Removing species (we remove the species one hundred times and estimate the average)

# a. random plants 
nsims_p = nplants*100

PL.degree=numeric(nsims_p)
PL.dead.plants=numeric(nsims_p)
PL.dead.pols=numeric(nsims_p)
PL.dead.disp=numeric(nsims_p)

set.seed(999)
bye.plants=rep(1:nplants,100) #

#Call function Stochastic Coextinction model
for (i in 1:nsims_p){
  NC=netcascade_multi(as.matrix(IM),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=NULL, deadPols=NULL, deadDisps=NULL,targetGuild='plant',target=bye.plants[i],return.matrix=F)
  PL.degree[i]=max(NC$cascade_data$degree)
  if(is.data.frame(NC$plant_species_data)) {
  PL.dead.plants[i]=length(NC$plant_species_data$lost_plants)
  }
  if (is.data.frame(NC$pol_species_data)) {
    PL.dead.pols[i]=length(NC$pol_species_data$lost_pols)
  }
  if (is.data.frame(NC$disp_species_data)){
  PL.dead.disp[i]=length(NC$disp_species_data$lost_disps)
  } 
}

#b. seed dispersers
nsims_d = length(disp_indx)*100

DP.degree=numeric(nsims_d)
DP.dead.plants=numeric(nsims_d)
DP.dead.pols=numeric(nsims_d)
DP.dead.disp=numeric(nsims_d)

set.seed(9999)
bye.disp=rep(disp_indx,100)


for (i in 1:nsims_d)
{
  NC=netcascade_multi(as.matrix(IM),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=NULL, deadPols=NULL, deadDisps=NULL,targetGuild='disperser',target=bye.disp[i],return.matrix=F)
  DP.degree[i]=max(NC$cascade_data$degree)
  if(is.data.frame(NC$plant_species_data)) {
    DP.dead.plants[i]=length(NC$plant_species_data$lost_plants)
  }
  if (is.data.frame(NC$pol_species_data)) {
    DP.dead.pols[i]=length(NC$pol_species_data$lost_pols)
  }
  if (is.data.frame(NC$disp_species_data)){
    DP.dead.disp[i]=length(NC$disp_species_data$lost_disps)
  } 
}

#c. pollinators
nsims_pol = npol*100

PO.degree=numeric(nsims_pol)
PO.dead.plants=numeric(nsims_pol)
PO.dead.pols=numeric(nsims_pol)
PO.dead.disp=numeric(nsims_pol)

set.seed(99999)
bye.pol= rep(pol_indx,100) 

for (i in 1:nsims_pol)
{
  NC=netcascade_multi(as.matrix(IM),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=NULL, deadPols=NULL, deadDisps=NULL,targetGuild='pollinator',target=bye.pol[i],return.matrix=F)
  PO.degree[i]=max(NC$cascade_data$degree)
  if(is.data.frame(NC$plant_species_data)) {
    PO.dead.plants[i]=length(NC$plant_species_data$lost_plants)
  }
  if (is.data.frame(NC$pol_species_data)) {
    PO.dead.pols[i]=length(NC$pol_species_data$lost_pols)
  }
  if (is.data.frame(NC$disp_species_data)){
    PO.dead.disp[i]=length(NC$disp_species_data$lost_disps)
  } 
}

#######################################################
####### Create data frame and summarize results #######
#######################################################

dead.plants=c(PL.dead.plants,DP.dead.plants,PO.dead.plants)
dead.disps=c(PL.dead.disp,DP.dead.disp,PO.dead.disp)
dead.pols=c(PL.dead.pols,DP.dead.pols,PO.dead.pols)
degrees=c(PL.degree,DP.degree,PO.degree)
initial.random=c(rep('Plant',nsims_p),rep('Disp',nsims_d),rep('Pol',nsims_pol))
Species_ID = c(rep(1:nplants,100), rep(disp_indx,100), rep(pol_indx,100))
Results=data.frame(Species_ID,initial.random,dead.plants,dead.disps,dead.pols,degrees)
Results=Results %>% mutate(porcDeadTot= ((dead.plants+dead.disps+dead.pols)/ntot)*100,porcDeadPlants=(dead.plants/nplants)*100,porcDeadDisp=(dead.disps/ndisp)*100,porcDeadPols=(dead.pols/npol)*100)


Results_perspeciesID<- Results %>% 
  group_by(Species_ID,initial.random) %>% 
  summarise_at(vars(c("dead.plants","dead.disps","dead.pols","degrees","porcDeadPlants","porcDeadTot",
                      "porcDeadDisp","porcDeadPols")),funs(mean(.))) %>% 
  arrange(initial.random)
Treat = rep("NI",ntot)

Results_final<- cbind(Treat = Treat,Results_perspeciesID)
Results_final2<- Results_final %>% 
  arrange(factor(initial.random, levels = c("Plant","Pol","Disp"))) %>% 
  mutate(Species_ID = 1:ntot) #we replace the Species ID of pollinator and seed disperser number to match their name in the node list

#change species ID for species name
Nodes_NI<- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/NI_nodelist.csv", sep=",",row.names=1)

Nodes=Nodes_NI

Results_final2$Species_ID<-as.factor(Results_final2$Species_ID)
levels(Results_final2$Species_ID) <- Nodes$name_node
factor(Results_final2$Species_ID)

#add the structural role to each species ID (previously calculated in modularity folder)
NI_sprole <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Modularity/NI_super_list.csv", sep=",",row.names=1)

NI_roles<- NI_sprole %>%
  group_by(node_id)%>% 
  dplyr::slice(1) %>% 
  dplyr::select(role)

Propagation_NI<-cbind(Results_final2,roles = NI_roles$role)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      I TREATMENT                            
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ndisp=I_disp
Npol=I_pol
Rdisp=Rdisp_I
Rpol=Rpol_I
rplants_pol=Rplant_pol_I
rplants_disp=Rplant_disp_I

Npol  <- Npol [ order(row.names(Npol)), ] #order pollinators

#Ri values for each species and trophic group
rdisp=Rdisp$x
rpol=Rpol$x
rplants_pol=rplants_pol$x
rplants_disp=rplants_disp$x

ndisp=dim(Ndisp)[1]
npol=dim(Npol)[1]
nplants=dim(Ndisp)[2]
nanim=npol+ndisp
ntot=nplants+npol+ndisp #Total number of nodes

# create final matrix
pol_indx=1:npol
disp_indx=(npol+1):(ndisp+npol)

imatrix=rbind(Npol,Ndisp) 
IM= as.matrix(imatrix)

####################
####  Scenarios ####
####################

# Removing species (we remove the species one hundred times and estimate the average)

# a. random plants 
nsims_p = nplants*100

PL.degree=numeric(nsims_p)
PL.dead.plants=numeric(nsims_p)
PL.dead.pols=numeric(nsims_p)
PL.dead.disp=numeric(nsims_p)

set.seed(999)
bye.plants=rep(1:nplants,100) #

#Call the function Stochastic Coextinction model
for (i in 1:nsims_p){
  NC=netcascade_multi(as.matrix(IM),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=NULL, deadPols=NULL, deadDisps=NULL,targetGuild='plant',target=bye.plants[i],return.matrix=F)
  PL.degree[i]=max(NC$cascade_data$degree)
  if(is.data.frame(NC$plant_species_data)) {
    PL.dead.plants[i]=length(NC$plant_species_data$lost_plants)
  }
  if (is.data.frame(NC$pol_species_data)) {
    PL.dead.pols[i]=length(NC$pol_species_data$lost_pols)
  }
  if (is.data.frame(NC$disp_species_data)){
    PL.dead.disp[i]=length(NC$disp_species_data$lost_disps)
  } 
}

#b. seed dispersers
nsims_d = length(disp_indx)*100

DP.degree=numeric(nsims_d)
DP.dead.plants=numeric(nsims_d)
DP.dead.pols=numeric(nsims_d)
DP.dead.disp=numeric(nsims_d)

set.seed(9999)
bye.disp=rep(disp_indx,100)


for (i in 1:nsims_d)
{
  NC=netcascade_multi(as.matrix(IM),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=NULL, deadPols=NULL, deadDisps=NULL,targetGuild='disperser',target=bye.disp[i],return.matrix=F)
  DP.degree[i]=max(NC$cascade_data$degree)
  if(is.data.frame(NC$plant_species_data)) {
    DP.dead.plants[i]=length(NC$plant_species_data$lost_plants)
  }
  if (is.data.frame(NC$pol_species_data)) {
    DP.dead.pols[i]=length(NC$pol_species_data$lost_pols)
  }
  if (is.data.frame(NC$disp_species_data)){
    DP.dead.disp[i]=length(NC$disp_species_data$lost_disps)
  } 
}

#c. pollinators
nsims_pol = npol*100

PO.degree=numeric(nsims_pol)
PO.dead.plants=numeric(nsims_pol)
PO.dead.pols=numeric(nsims_pol)
PO.dead.disp=numeric(nsims_pol)

set.seed(99999)
bye.pol= rep(pol_indx,100) 

for (i in 1:nsims_pol)
{
  NC=netcascade_multi(as.matrix(IM),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=NULL, deadPols=NULL, deadDisps=NULL,targetGuild='pollinator',target=bye.pol[i],return.matrix=F)
  PO.degree[i]=max(NC$cascade_data$degree)
  if(is.data.frame(NC$plant_species_data)) {
    PO.dead.plants[i]=length(NC$plant_species_data$lost_plants)
  }
  if (is.data.frame(NC$pol_species_data)) {
    PO.dead.pols[i]=length(NC$pol_species_data$lost_pols)
  }
  if (is.data.frame(NC$disp_species_data)){
    PO.dead.disp[i]=length(NC$disp_species_data$lost_disps)
  } 
}

#######################################################
####### Create data frame and summarize results #######
#######################################################

dead.plants=c(PL.dead.plants,DP.dead.plants,PO.dead.plants)
dead.disps=c(PL.dead.disp,DP.dead.disp,PO.dead.disp)
dead.pols=c(PL.dead.pols,DP.dead.pols,PO.dead.pols)
degrees=c(PL.degree,DP.degree,PO.degree)
initial.random=c(rep('Plant',nsims_p),rep('Disp',nsims_d),rep('Pol',nsims_pol))
Species_ID = c(rep(1:nplants,100), rep(disp_indx,100), rep(pol_indx,100))
Results=data.frame(Species_ID,initial.random,dead.plants,dead.disps,dead.pols,degrees)
Results=Results %>% mutate(porcDeadTot= ((dead.plants+dead.disps+dead.pols)/ntot)*100,porcDeadPlants=(dead.plants/nplants)*100,porcDeadDisp=(dead.disps/ndisp)*100,porcDeadPols=(dead.pols/npol)*100)


Results_perspeciesID<- Results %>% 
  group_by(Species_ID,initial.random) %>% 
  summarise_at(vars(c("dead.plants","dead.disps","dead.pols","degrees","porcDeadPlants","porcDeadTot",
                      "porcDeadDisp","porcDeadPols")),funs(mean(.))) %>% 
  arrange(initial.random)
Treat = rep("I",ntot)

Results_final<- cbind(Treat = Treat,Results_perspeciesID)
Results_final2<- Results_final %>% 
  arrange(factor(initial.random, levels = c("Plant","Pol","Disp"))) %>% 
  mutate(Species_ID = 1:ntot) #we replace the Species ID of pollinator and seed disperser number to match their name in the node list

#change species ID for species name
Nodes_I<- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/I_nodelist.csv", sep=",",row.names=1)

Nodes=Nodes_I

Results_final2$Species_ID<-as.factor(Results_final2$Species_ID)
levels(Results_final2$Species_ID) <- Nodes$name_node
factor(Results_final2$Species_ID)

#add the structural role to each species ID (previously calculated in modularity folder)
I_sprole <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Modularity/I_super_list.csv", sep=",",row.names=1)

I_roles<- I_sprole %>%
  group_by(node_id)%>% 
  dplyr::slice(1) %>% 
  dplyr::select(role)

Propagation_I<-cbind(Results_final2,roles = I_roles$role)

#######################################################
####### Create FINAL data frame  #######
#######################################################

Propagation<- rbind(Propagation_NI,Propagation_I)

-------------------------------------------------------------------------------
########### 2. Comparison between treatments
library(lme4)
library("stats4")
library("bbmle")
library(ggplot2)
library(car)
library(lsmeans)
library(multcomp)

#Prepare the dataframe
Propagation[Propagation == 0] <- 0.0000001 #reeplace for a positive value to run gamma
Propagation$Species_ID = as.factor(Propagation$Species_ID)
Propagation$Treat<- factor(Propagation$Treat, levels=c("NI","I"))
Propagation$initial.random <- factor(Propagation$initial.random , levels=c("Pol","Plant", "Disp"))
Propagation$roles <- factor(Propagation$roles , levels=c("network hub","module hub","connector", "peripheral"))


#Compare extinct percentage of species according to Treat and structural role of species----
m1<- glmer (porcDeadTot~ Treat* roles +( 1| Species_ID)-1, family = Gamma, data = Propagation)
summary(m1)
Anova(m1)#Interaction between fixed factors was not significant, so we incorporated them in a separated way

#Homogeneity
E_lme<-resid(m1, type= "deviance") 
F_lme<-fitted(m1) 
plot(x=F_lme, y=E_lme, xlab = "Ajustados", ylab = "Residuales normalizados")
abline(0,0, col="red", lwd= 3)

##Independence 
boxplot(E_lme~Propagation$Treat, main="Tratamiento")
boxplot(E_lme~Propagation$roles, main="roles") 

#posteriori test
library(lsmeans)
lsm <- lsmeans(m1, ~ Treat*roles)
summary(pairs(lsm), type = "response") 
cld(lsm, 
    alpha=.05,
    Letters=letters)

#Plot--
ext_order=ggplot(Propagation)+geom_boxplot(aes(y=porcDeadTot,x=roles,fill=Treat),  position = position_dodge2(preserve = "single"))+
  scale_fill_manual(values=c("cornflowerblue","orange2"))+
  scale_x_discrete(name="Topological role")+
  scale_y_continuous(name = "Percentage of extinct species")+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line (colour = "black",linetype = "solid", size = 1.4 ),
        axis.title = element_text(size= 15),
        axis.ticks = element_line(colour = "black", size = 1.4),
        axis.text = element_text (colour = "black", size = 14 ),
        legend.background = element_blank(),
        legend.key = element_rect(colour = "white", fill = NA))
ext_order


####Extinction percentage for each trophic group-- (Extra analysis appendix S3)

#Compare if pollinators that goes extinct after the removal of a species change according to the 
#trophic group and treatment----

m2<- glmer (porcDeadPols~ Treat+initial.random + ( 1| Species_ID), family = Gamma, data = Propagation)
summary(m2)
Anova(m2)

#Homogeneity 
E_lme<-resid(m2, type= "deviance") 
F_lme<-fitted(m2) 
plot(x=F_lme, y=E_lme, xlab = "Ajustados", ylab = "Residuales normalizados")
abline(0,0, col="red", lwd= 3)

##Independence 
boxplot(E_lme~Propagation$Treat, main="Tratamiento")
boxplot(E_lme~Propagation$initial.random, main="initial.random") 

#posteriori test
t2<-glht(m2, linfct=mcp(initial.random="Tukey"), adjust.method="fdr")
summary(t2)

#Plot-
#function to calculate the min, mean-1SEM, mean, mean+1SEM, and Max.Then put them into a boxplot
parameters <- function(x) {
  v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}

g2<-ggplot(Propagation, aes(x= Treat ,y=porcDeadPols, fill = Treat)) +
  stat_summary(fun.data=parameters, geom="boxplot", position = position_dodge(1))+
  facet_wrap(~initial.random) +
  scale_fill_manual(values=c("cornflowerblue","orange2"))+
  scale_x_discrete(name="Initial extinction")+
  scale_y_continuous(name = "Extinct pollinators (%)")+
  theme(panel.grid.minor =  element_line(color = "gray"),
        panel.background = element_blank(),
        axis.line = element_line (colour = "black",linetype = "solid", size = 1.4 ),
        axis.title = element_text(size= 15),
        axis.text = element_text (colour = "black", size = 14 ),
        axis.ticks = element_line(colour = "black", size = 1.4),
        legend.background = element_blank(),
        legend.key = element_rect(colour = "white", fill = NA))
g2 


#Compare if plants that goes extinct after the removal of a species change according to the
#trophic group and treatment----

m3<- glmer (porcDeadPlants~ Treat+initial.random + ( 1| Species_ID), family = Gamma, data = Propagation)
summary(m3)
Anova(m3)

#Homogeneity 
E_lme<-resid(m3, type= "deviance") 
F_lme<-fitted(m3) 
plot(x=F_lme, y=E_lme, xlab = "Ajustados", ylab = "Residuales normalizados")
abline(0,0, col="red", lwd= 3)

##Independence 
boxplot(E_lme~Propagation$Treat, main="Tratamiento")
boxplot(E_lme~Propagation$initial.random, main="initial.random") 

#posteriori test
t3<-glht(m3, linfct=mcp(initial.random="Tukey"), adjust.method="fdr")
summary(t3)

#Plot-
g3<-ggplot(Propagation, aes(x= Treat ,y=porcDeadPlants, fill = Treat)) +
  stat_summary(fun.data=parameters, geom="boxplot", position = position_dodge(1))+
  facet_wrap(~initial.random) +
  scale_fill_manual(values=c("cornflowerblue","orange2"))+
  scale_x_discrete(name="Initial extinction")+
  scale_y_continuous(name = "Extinct plants (%)")+
  theme(panel.grid.minor =  element_line(color = "gray"),
        panel.background = element_blank(),
        axis.line = element_line (colour = "black",linetype = "solid", size = 1.4 ),
        axis.title = element_text(size= 15),
        axis.text = element_text (colour = "black", size = 14 ),
        axis.ticks = element_line(colour = "black", size = 1.4),
        legend.background = element_blank(),
        legend.key = element_rect(colour = "white", fill = NA))
g3


#Compare if seed disperser that goes extinct after the removal of a species change according to the
#trophic group and treatment----
m4 <-glmer (porcDeadDisp~ Treat+initial.random + ( 1| Species_ID), family = Gamma, data = Propagation)
summary(m4)
Anova(m4)

#Homogeneity 
E_lme<-resid(m4, type= "deviance") 
F_lme<-fitted(m4) 
plot(x=F_lme, y=E_lme, xlab = "Ajustados", ylab = "Residuales normalizados")
abline(0,0, col="red", lwd= 3)

##Independence 
boxplot(E_lme~Propagation$Treat, main="Tratamiento")
boxplot(E_lme~Propagation$initial.random, main="initial.random") 

#posteriori test
t4<-glht(m4, linfct=mcp(initial.random="Tukey"), adjust.method="fdr")
summary(t4)

#Plot--
g4<-ggplot(Propagation, aes(x= Treat ,y=porcDeadDisp, fill = Treat)) +
  stat_summary(fun.data=parameters, geom="boxplot", position = position_dodge(1))+
  facet_wrap(~initial.random) +
  scale_fill_manual(values=c("cornflowerblue","orange2"))+
  scale_x_discrete(name="Initial extinction")+
  scale_y_continuous(name = "Extinct seed dispersers (%)")+
  theme(panel.grid.minor =  element_line(color = "gray"),
        panel.background = element_blank(),
        axis.line = element_line (colour = "black",linetype = "solid", size = 1.4 ),
        axis.title = element_text(size= 15),
        axis.text = element_text (colour = "black", size = 14 ),
        axis.ticks = element_line(colour = "black", size = 1.4),
        legend.background = element_blank(),
        legend.key = element_rect(colour = "white", fill = NA))
g4


