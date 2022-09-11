########### 1. Estimation of robustness in networks controlled by size (NI treatment) 

########### 2. Comparison with matrix of NI without controlling by size
#-----------------------------------------------------------------------------------------------------------------

########### 1. Estimation of robustness in networks controlled by size (NI treatment) 

library(tidyverse)
library(dplyr)

source('/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Extra_Functions/Extrafunctions.R')
source("Simulation_control_size.R") #simulated networks controlled by size


# Load Ri values to assign of each species in each network
tdisp <- read.csv("traits_frug.csv", sep = ";")
tpol<- read.csv("traits_pol.csv", sep = ",")
tplants_pol<- read.csv("traits_plants_pol_2.csv")
tplants_disp<- read.csv("traits_plants_disp_2.csv")
tpol %>% arrange(Species_Ref) #order the pollinators in the R list according alphabetic order

# List of roles to assign of each species
sprole <- read.csv("/Users/agustin/Desktop/Papers/Paper_multicapa/Matrices_consecos/ISRAEL/Analisis_pooled/Modularity/Modularity_sproles/Change_sprole/super_list_treat.csv", sep=",",row.names=1)

#loops
Final=data.frame()
for(m in res){
  
#we assign the R for each species in each matrix --
  
  #plants_pols
  plant.common=colnames(m)[colnames(m)%in%tplants_pol[,1]]
  Rplants_pol<- tplants_pol %>% 
    filter(Plant_name%in%plant.common)
  rplants_pol<- Rplants_pol[,5]
  
  #plants_disps
  plant.common=colnames(m)[colnames(m)%in%tplants_disp[,1]]
  Rplants_disp<- tplants_disp %>% 
    filter(Plant_name%in%plant.common)
  rplants_disp<- Rplants_disp[,4]
  
  #pollinators
  pol.common = rownames(m)[rownames(m)%in%tpol[,2]]
  Rpol<-tpol %>% 
    filter(Species_Ref%in%pol.common) %>% 
    arrange(Species_Ref)
  rpol<-Rpol[,3]
  
  #seed dispersers
  disp.common=row.names(m)[row.names(m)%in%tdisp[,1]]
  Rdisp<- tdisp %>% 
    filter(Seeddisp_name%in%disp.common)
  rdisp<-Rdisp[,7]
  
  ndisp= length(rdisp)
  npol= length(rpol)
  nplants=ncol(m)
  nanim= nrow(m)
  ntot=nplants+npol+ndisp 
  pol_indx=1:npol
  disp_indx=(npol+1):(ndisp+npol)
  plants_names<-colnames(m)
  pol_names<-rownames(m[1:npol,]) 
  disp_names <- rownames(m[disp_indx,])
  
#assign the species role to nodes
  roles_NI<- sprole %>%
    filter(Treat == "NI") %>% 
    group_by(node_id)%>% 
    slice(1) %>% 
    select(role)
  
  roles_plants<- roles_NI %>% 
    filter(node_id%in%plants_names) %>% 
    select(role) 
  
  roles_pol<- roles_NI %>% 
    filter(node_id%in%pol_names) %>% 
    select(role) 
  
  roles_disp<- roles_NI%>% 
    filter(node_id%in%disp_names) %>% 
    select(role) 
  
  roles_m<-rbind(roles_plants,roles_pol,roles_disp)
  
  #assign a value (1:4) to each role classification
  n.roles <- roles_m %>% 
    mutate(role, role = case_when (role == "network hub"~ 4,
                                   role == "module hub"~ 3,
                                   role == "connector"~ 2,
                                   role == "peripheral"~ 1)) 
  n.roles <- n.roles$role

  ####################
  ####  Scenarios ####
  ####################
  
  # a) Random Species 
  # b) Most to least 
  # c) least to most 
  
  scenario=c('random','MtoL','LtoM') #scenarios a,b and c
  Scenario=c()
  AUC=matrix(NA,ncol=3,nrow=1)
  sp.rm=c()
  sp.present=c()  
  
  for (j in 1:3)
  {
    if (scenario[j]=='random')
    {
      bye.s=sample(1:length(n.roles),size=nsims, replace = TRUE) #select species at random to remove
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
        n.rolesNA[deads]=-999 #replace the classification of the dead species for "-999", so the algorithm doesn't chose dead species to remove
        ind.max = which(n.rolesNA==max(n.rolesNA))#select an alive species with the maximum classification to remove
        next.rem=ind.max[sample(length(ind.max),size=1)]
        
      }
      
      if(scenario[j]=='LtoM')
      {
        n.rolesNA[deads]=9999#replace the classification of the dead species for "9999", so the , so the algorithm doesn't chose dead species to remove.
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
      NC=netcascade_multi(as.matrix(m),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=deadPlants, deadPols=deadPols, deadDisps=deadDisps, targetGuild=tg ,target=next.rem.fc,return.matrix=T)
      
      deadPlants=c(deadPlants,NC$lost_plants)#deadPlants represent the accumulated dead plants, while lost plants refers to the plant species lost in the current step
      deadPols=c(deadPols,NC$lost_pols)#same for pollinators
      deadDisps=c(deadDisps,NC$lost_disps)#same for seed dispersers
      
      all.dead=sum(length(deadPlants),length(deadPols),length(deadDisps))#total number of dead species
      deads=c(deadPols+nplants,deadDisps+nplants,deadPlants)
      n.remove=c(n.remove,1)#number of species removed
      n.dead=c(n.dead,sum(length(NC$lost_plants),length(NC$lost_pols),length(NC$lost_disps)))
    
    
    sp.rms=c(0,cumsum(n.remove))/max(cumsum(n.remove))#Prop of removed species per simulation
    sp.presents=c(ntot,sum(n.dead)-cumsum(n.dead))/ntot #Prop of alive species after each removal
    
    #Call the function AUC
    Scenario=c(Scenario,rep(scenario[j],length(sp.rms)))
    AUC[j]=auc(sp.rms,sp.presents)#calculate the area under the curve (proxy of robustness)
    sp.rm=c(sp.rm,sp.rms)
    sp.present=c(sp.present,sp.presents)
    
  }
  Final = rbind(Final,AUC)
  }
  
  Rob_NI_size = Final #Robustness of 300 networks controlled by size
  
  #######################################################
  ####### Create data frame and summarize results #######
  #######################################################
  
colnames(Rob_NI_size)<-c("random","MtoL","LtoM")
  
NI_c2<- Rob_NI_size %>%
  gather(key = "Scenario", value = "Auc")

#-------------------------------------------------------------------------------
########### 2. Comparison with matrix of NI without controlling by size

#upload matrix of NI without controlling for size
NI_nc<-read.csv() #should upload the file containing the robustness of NI networks generated in the code titled "Extinction_robustness"
NI_nc2<-NI_nc %>% 
  filter(Treat=="NI") %>% 
  select(-reps) 
NI_nc3<-NI_nc2 %>% 
  group_by(Scenario) %>% 
  slice_sample(n =300) %>% #select randomly 300 simulation for each scenario
select(-(Treat))

#create final dataframe with controlled and not controlled data
Final<-rbind(NI_c2,NI_nc3)
Final2<-cbind(Treat = rep(c("Controlled","Uncontrolled"), each = 900), Final)
Final2$Treat <- factor(Final2$Treat , levels=c("Controlled","Uncontrolled"))
Final2$Scenario <- factor(Final2$Scenario , levels=c("MtoL","LtoM", "random"))

##Model to compare network robustness according to uncontrolled vs controlled data -
library(lme4)
library("stats4")
library("bbmle")
library(car)
library(lsmeans)
library(multcomp)

m1<- glm (Auc~ Treat, family = Gamma, data = Final2)
summary(m1)
Anova(m1, test.statistic = "F")

#Homogeneity 
E_lme<-resid(m1, type= "deviance") 
F_lme<-fitted(m1) 
plot(x=F_lme, y=E_lme, xlab = "Ajustados", ylab = "Residuales normalizados")
abline(0,0, col="red", lwd= 3)

##Independence
boxplot(E_lme~Final2$Treat, main="Tratamiento")
boxplot(E_lme~Final2$Scenario, main="Scenario") 


###Plot (Appendix S3) -
aucplot=ggplot(Final2)+geom_boxplot(aes(y=Auc,x=Scenario,fill=Treat))+
  scale_fill_manual(values=c("#a69880", "cornflowerblue"))+
  scale_x_discrete(name="Order of species removal")+
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




