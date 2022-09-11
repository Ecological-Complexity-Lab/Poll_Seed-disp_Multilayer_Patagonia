#This file has the codes to estimate the robustness of each network  and compare it between treatments

########### 1. Estimation of the robustness for each treatment

########### 2. Comparison between treatments

#-----------------------------------------------------------------------------------------------------------------
###### 1. Estimation of the robustness for each treatment --

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

#add the structural role to each species 

#upload species role 
NI_sprole <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Modularity/NI_super_list.csv", sep=",",row.names=1)
roles<- NI_sprole %>%
  group_by(node_id)%>% 
  dplyr::slice(1) %>% 
  dplyr::select(role)

#assign a value (1:4) to each role classification
n.roles <- roles %>% 
  mutate(role, role = case_when (role == "network hub"~ 4,
                                 role == "module hub"~ 3,
                                 role == "connector"~ 2,
                                 role == "peripheral"~ 1)) 
n.roles <- n.roles$role

####################
####  Scenarios ####
####################
nsims=1000

# a) Random Species 
# b) Most to least 
# c) least to most 

scenario=c('random','MtoL','LtoM') #scenarios a,b and c
Scenario=c()
reps=c()
AUC=matrix(NA,ncol=3,nrow=nsims)
sp.rm=c()
sp.present=c()

for (j in 1:3)
{
  if (scenario[j]=='random')
  {
    bye.s=sample(1:length(n.roles),size=nsims, replace = TRUE)#select species at random to remove
  }
  
  if(scenario[j]=='MtoL')
  {
    ind.max = which(n.roles==max(n.roles))#select species with the maximum role classification to remove
    bye.s=sample(ind.max,nsims, replace =T) 
  }
  
  if(scenario[j]=='LtoM')
  {
    ind.min = which(n.roles==min(n.roles))#species with the minimum role classification
    bye.s=sample(ind.min,nsims, replace =T)   
  }
  
  
  for (i in 1:nsims)
  {
   
    all.dead=1
    deadPlants=NULL
    deadPols=NULL
    deadDisps=NULL
    
    n.remove=c()
    n.dead=c()
    deads=NULL
    
    while(all.dead<(ntot-1))#select the next species to remove after finishing the cascade effects of the previous one
    {
    
      if (scenario[j]=='random')
      {
        quienes=(1:ntot)[!(1:ntot)%in% deads]
        next.rem=sample(quienes,1)#randomly select an alive species to remove
      }
      
      if(scenario[j]=='MtoL')
      {
        n.rolesNA[deads]=-999 #replace the classification of the dead species for "-999", so the algorithm doesn't chose dead species to remove.
        ind.max = which(n.rolesNA==max(n.rolesNA))#select an alive species with the maximum classification to remove
        next.rem=ind.max[sample(length(ind.max),size=1)]
       
      }
      
      if(scenario[j]=='LtoM')
      {
        n.rolesNA[deads]=9999 #replace the classification of the dead species for "9999", so the , so the algorithm doesn't chose dead species to remove.
        ind.min = which(n.rolesNA==min(n.rolesNA))#select an alive species with the minimum classification to remove
        next.rem=ind.min[sample(length(ind.min),size=1)]
      }
     #identify the trophic group of the species to remove according to the node ID
      if(next.rem>(length(rpol)+nplants))
      {tg="disperser"
      next.rem.fc = next.rem-nplants
      }
      if(nplants < next.rem & next.rem <= (length(rpol)+nplants))
      {tg="pollinator"
      next.rem.fc = next.rem-nplants 
      }
      if(nplants >= next.rem)
      { tg="plant"
        next.rem.fc = next.rem
      }
      n.rolesNA=n.roles 
      
      #Call the function Stochastic Coextinction model
      
      NC=netcascade_multi(as.matrix(IM),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=deadPlants, deadPols=deadPols, deadDisps=deadDisps, targetGuild=tg ,target=next.rem.fc,return.matrix=T)
      
      deadPlants=c(deadPlants,NC$lost_plants)#deadPlants represent the accumulated dead plants, while lost plants refers to the plant species lost in the current step
      deadPols=c(deadPols,NC$lost_pols)#same for pollinators
      deadDisps=c(deadDisps,NC$lost_disps)#same for seed dispersers
      
      all.dead=sum(length(deadPlants),length(deadPols),length(deadDisps))#total number of dead species
      deads=c(deadPols+nplants,deadDisps+nplants,deadPlants)
      n.remove=c(n.remove,1)#number of species removed
      n.dead=c(n.dead,sum(length(NC$lost_plants),length(NC$lost_pols),length(NC$lost_disps)))
      
    }

    sp.rms=c(0,cumsum(n.remove))/max(cumsum(n.remove))#Prop of removed species per simulation
    sp.presents=c(ntot,sum(n.dead)-cumsum(n.dead))/ntot #Prop of alive species after each removal
    
    Scenario=c(Scenario,rep(scenario[j],length(sp.rms)))
    reps=c(reps,rep(i,length(sp.rms)))
    
    #Call the function AUC
    AUC[i,j]=auc(sp.rms,sp.presents)#calculate the area under the curve (proxy of robustness)
    sp.rm=c(sp.rm,sp.rms)
    sp.present=c(sp.present,sp.presents)
    
  }
  
}


#######################################################
####### Create data frame and summarize results #######
#######################################################

reps2=rep(1:nsims,3)
Scenario2=rep(scenario,each=nsims)
Auc=as.vector(AUC)
Results.AUC=data.frame(reps=reps2,Scenario=Scenario2,Auc)  
Results.AUC_NI = Results.AUC


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

#add the structural role to each species 

#upload species role 
I_sprole <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Modularity/I_super_list.csv", sep=",",row.names=1)
roles<- I_sprole %>%
  group_by(node_id)%>% 
  dplyr::slice(1) %>% 
  dplyr::select(role)

#assign a value (1:4) to each role classification
n.roles <- roles %>% 
  mutate(role, role = case_when (role == "network hub"~ 4,
                                 role == "module hub"~ 3,
                                 role == "connector"~ 2,
                                 role == "peripheral"~ 1)) 
n.roles <- n.roles$role

nsims=1000

####################
####  Scenarios ####
####################

# a) Random Species 
# b) Most to least 
# c) least to most 

scenario=c('random','MtoL','LtoM') #scenarios a,b and c
Scenario=c()
reps=c()
AUC=matrix(NA,ncol=3,nrow=nsims)
sp.rm=c()
sp.present=c()

for (j in 1:3)
{
  if (scenario[j]=='random')
  {
    bye.s=sample(1:length(n.roles),size=nsims, replace = TRUE)#select species at random to remove
  }
  
  if(scenario[j]=='MtoL')
  {
    ind.max = which(n.roles==max(n.roles))#select species with the maximum role classification to remove
    bye.s=sample(ind.max,nsims, replace =T) 
  }
  
  if(scenario[j]=='LtoM')
  {
    ind.min = which(n.roles==min(n.roles))#species with the minimum role classification
    bye.s=sample(ind.min,nsims, replace =T)   
  }
  
  
  for (i in 1:nsims)
  {
    
    all.dead=1
    deadPlants=NULL
    deadPols=NULL
    deadDisps=NULL
    
    n.remove=c()
    n.dead=c()
    deads=NULL
    
    while(all.dead<(ntot-1))#select the next species to remove after finishing the cascade effects of the previous one
    {
      
      if (scenario[j]=='random')
      {
        quienes=(1:ntot)[!(1:ntot)%in% deads]
        next.rem=sample(quienes,1)#randomly select an alive species to remove
      }
      
      if(scenario[j]=='MtoL')
      {
        n.rolesNA[deads]=-999 #replace the classification of the dead species for "-999", so the algorithm doesn't chose dead species to remove.
        ind.max = which(n.rolesNA==max(n.rolesNA))#select an alive species with the maximum classification to remove
        next.rem=ind.max[sample(length(ind.max),size=1)]
        
      }
      
      if(scenario[j]=='LtoM')
      {
        n.rolesNA[deads]=9999 #replace the classification of the dead species for "9999", so the , so the algorithm doesn't chose dead species to remove.
        ind.min = which(n.rolesNA==min(n.rolesNA))##select an alive species with the minimum classification to remove
        next.rem=ind.min[sample(length(ind.min),size=1)]
      }
      #identify the trophic group of the species to remove according to the node ID
      if(next.rem>(length(rpol)+nplants))
      {tg="disperser"
      next.rem.fc = next.rem-nplants
      }
      if(nplants < next.rem & next.rem <= (length(rpol)+nplants))
      {tg="pollinator"
      next.rem.fc = next.rem-nplants 
      }
      if(nplants >= next.rem)
      { tg="plant"
      next.rem.fc = next.rem
      }
      n.rolesNA=n.roles 
      
      #Call the function Stochastic Coextinction model
      
      NC=netcascade_multi(as.matrix(IM),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=deadPlants, deadPols=deadPols, deadDisps=deadDisps, targetGuild=tg ,target=next.rem.fc,return.matrix=T)
      
      deadPlants=c(deadPlants,NC$lost_plants)#deadPlants represens the accumulated dead plants, while lost plants refers to the plant species lost in the current step
      deadPols=c(deadPols,NC$lost_pols)#same for pollinators
      deadDisps=c(deadDisps,NC$lost_disps)#same for seed dispersers
      
      all.dead=sum(length(deadPlants),length(deadPols),length(deadDisps))#total number of dead species
      deads=c(deadPols+nplants,deadDisps+nplants,deadPlants)
      n.remove=c(n.remove,1)#number of species removed
      n.dead=c(n.dead,sum(length(NC$lost_plants),length(NC$lost_pols),length(NC$lost_disps)))
      
    }
    
    sp.rms=c(0,cumsum(n.remove))/max(cumsum(n.remove))#Prop of removed species per simulation
    sp.presents=c(ntot,sum(n.dead)-cumsum(n.dead))/ntot #Prop of alive species after each removal
    
    Scenario=c(Scenario,rep(scenario[j],length(sp.rms)))
    reps=c(reps,rep(i,length(sp.rms)))
    
    #Call the function AUC
    AUC[i,j]=auc(sp.rms,sp.presents)#calculate the area under the curve (proxy of robustness)
    sp.rm=c(sp.rm,sp.rms)
    sp.present=c(sp.present,sp.presents)
    
  }
  
}

#######################################################
####### Create data frame and summarize results #######
#######################################################

reps2=rep(1:nsims,3)
Scenario2=rep(scenario,each=nsims)
Auc=as.vector(AUC)
Results.AUC=data.frame(reps=reps2,Scenario=Scenario2,Auc)  
Results.AUC_I = Results.AUC

#######################################################
####### Create FINAL data frame  #######
#######################################################

Results_AUC<-rbind(Results.AUC_NI,Results.AUC_I)
Treat<-rep(c("NI","I"), each = 3000)
Robustness<-cbind(Treat,Results_AUC)


#-------------------------------------------------------------------------------
########### 2. Comparison between treatments

library(lme4)
library(lsmeans)
library("stats4")
library("bbmle")
library(car)
library(lsmeans)
library(multcomp)

#Prepare data
Robustness$Treat <- factor(Robustness$Treat , levels=c("NI","I"))
Robustness$Scenario <- factor(Robustness$Scenario , levels=c("MtoL","LtoM", "random"))

#Model
m1<- glm (Auc~ Treat*Scenario, family = Gamma, data = Robustness)
summary(m1)
Anova(m1)

#Homogeneity (residuals vs fitted values)
E_lme<-resid(m1, type= "deviance") 
F_lme<-fitted(m1) 
plot(x=F_lme, y=E_lme, xlab = "Ajustados", ylab = "Residuales normalizados")
abline(0,0, col="red", lwd= 3)

#Independence 
boxplot(E_lme~Robustness$Treat, main="Tratamiento")
boxplot(E_lme~Robustness$Scenario, main="Scenario") 

#posterior with interaction
lsm <- lsmeans(m1, ~ Treat*Scenario)
summary(pairs(lsm), type = "response") #comparision
cld(lsm, 
    alpha=.05,
    Letters=letters)

#Plot- 
aucplot=ggplot(Robustness)+geom_boxplot(aes(y=Auc,x=Scenario,fill=Treat))+
  scale_fill_manual(values=c("cornflowerblue","orange2"))+
  scale_x_discrete(name="order of species removal")+
  scale_y_continuous(name = "Area under the curve (AUC)")+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line (colour = "black",linetype = "solid", size = 1.4 ),
        axis.title = element_text(size= 15),
        axis.ticks = element_line(colour = "black", size = 1.4),
        axis.text = element_text (colour = "black", size = 14 ),
        legend.background = element_blank(),
        legend.key = element_rect(colour = "white", fill = NA))
aucplot






